#!/usr/bin/env python3
"""
Step 7b: MetaNeighbor Pass 2 — resolve unresolved clusters from Pass 1.

Loads the lists of unresolved mouse and human clusters produced by
07_cluster_bridge.py (Pass 1), subsamples cells from the original atlases,
runs MetaNeighbor AUROC, and merges Pass 1 + Pass 2 results into a single
all_cluster_correspondences.csv.

Inputs (from --bridge-dir, defaults to results/cluster_bridge/)
------
  unresolved_mouse.csv   : list of mouse cluster IDs not resolved in Pass 1
  unresolved_human.csv   : list of human cluster IDs not resolved in Pass 1
  rbh_resolved_pairs.csv : Pass 1 RBH resolved pairs

Outputs (written to --outdir, defaults to results/cluster_bridge/)
-------
  metaneighbor_auroc_matrix.csv   : AUROC matrix (n_unresolved_mouse x n_unresolved_human)
  pass2_resolved_pairs.csv        : best MetaNeighbor match per unresolved mouse cluster
  all_cluster_correspondences.csv : Pass 1 + Pass 2 combined
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from cge_subtype.src.cluster_correspondence import compute_metaneighbor_auroc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------
_MOUSE_H5AD_FILES = [
    Path(
        "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices"
        "/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-raw.h5ad"
    ),
    Path(
        "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices"
        "/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-2-raw.h5ad"
    ),
]
_MOUSE_META_CSV = Path(
    "/mnt/data0/AllenMouseSC/abc_download_root/metadata"
    "/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv"
)
_HUMAN_H5AD_DIR = Path("/mnt/data0/HumanBrainCellType/SuperTypeRawDat")
_DEFAULT_BRIDGE_DIR = PROJECT_ROOT / "cge_subtype" / "results" / "cluster_bridge"
_DEFAULT_OUTDIR = PROJECT_ROOT / "cge_subtype" / "results" / "cluster_bridge"
_DEFAULT_ORTHOLOG = Path(
    "/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results"
    "/cross_species/orthologs/ortholog_mapping.csv"
)

CHUNK_SIZE = 1000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_ortholog_map(ortholog_path: Path) -> dict[str, str]:
    """Return {mouse_symbol -> human_symbol} for one-to-one orthologs."""
    df = pd.read_csv(ortholog_path)
    one2one = df[df["is_one_to_one"].astype(bool)]
    return dict(zip(one2one["mouse_symbol"], one2one["human_symbol"]))


def _get_gene_names(adata: ad.AnnData) -> list[str]:
    for col in ("Gene", "gene_symbol", "gene_name"):
        if col in adata.var.columns:
            return adata.var[col].tolist()
    return adata.var_names.tolist()


def _get_cluster_col(adata: ad.AnnData) -> str | None:
    for col in ("cluster_id", "cluster", "supertype", "cell_type"):
        if col in adata.obs.columns:
            return col
    return None


def subsample_cells(
    cell_to_cluster: pd.Series,
    target_clusters: list[str],
    max_per_cluster: int,
    seed: int = 42,
) -> list[str]:
    """Return a list of cell IDs subsampled to ``max_per_cluster`` per cluster."""
    rng = np.random.default_rng(seed)
    target_set = set(target_clusters)
    sampled: list[str] = []

    for cluster in target_clusters:
        cells = cell_to_cluster[cell_to_cluster == cluster].index.tolist()
        if len(cells) > max_per_cluster:
            chosen = rng.choice(cells, size=max_per_cluster, replace=False).tolist()
        else:
            chosen = cells
        sampled.extend(chosen)

    log.info(
        "  Subsampled %d cells across %d clusters (max %d/cluster).",
        len(sampled),
        len(target_clusters),
        max_per_cluster,
    )
    return sampled


# ---------------------------------------------------------------------------
# Load expression for a set of cells from backed h5ad files
# ---------------------------------------------------------------------------
def _load_cells_from_h5ads(
    h5ad_files: list[Path],
    target_cells: set[str],
    gene_names_ref: list[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Read expression rows for ``target_cells`` from a list of h5ad files.

    Returns (expr_df, gene_names) where expr_df is (cells x genes).
    """
    rows: dict[str, np.ndarray] = {}
    gene_names: list[str] | None = gene_names_ref

    for h5ad_path in h5ad_files:
        if not h5ad_path.exists():
            log.warning("h5ad not found: %s — skipping.", h5ad_path)
            continue

        adata = ad.read_h5ad(h5ad_path, backed="r")
        if gene_names is None:
            gene_names = _get_gene_names(adata)

        overlap = list(set(adata.obs_names.tolist()) & target_cells)
        if not overlap:
            adata.file.close()
            continue

        log.info("  Loading %d cells from %s ...", len(overlap), h5ad_path.name)
        positions = [adata.obs_names.get_loc(c) for c in overlap]
        chunk = adata[positions].X
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        else:
            chunk = np.asarray(chunk)

        for i, cell_id in enumerate(overlap):
            rows[cell_id] = chunk[i].astype(np.float32)

        adata.file.close()

    if not rows:
        return pd.DataFrame(), gene_names or []

    expr_df = pd.DataFrame(rows, index=gene_names).T  # cells x genes
    expr_df.index.name = "cell_id"
    return expr_df, gene_names


def _load_human_cells_from_h5ads(
    h5ad_dir: Path,
    target_cells: set[str],
    target_clusters: set[str],
    max_per_cluster: int = 200,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.Series]:
    """Load expression and cluster labels for target human cells (subsampled)."""
    h5ad_files = sorted(h5ad_dir.glob("Supercluster_*.h5ad"))
    rng = np.random.default_rng(seed)

    # Track how many cells we've collected per cluster
    cluster_cell_counts: dict[str, int] = {}
    rows: dict[str, np.ndarray] = {}
    labels: dict[str, str] = {}
    gene_names: list[str] | None = None

    for h5ad_path in h5ad_files:
        adata = ad.read_h5ad(h5ad_path, backed="r")
        cluster_col = _get_cluster_col(adata)
        if cluster_col is None:
            adata.file.close()
            continue

        if gene_names is None:
            gene_names = _get_gene_names(adata)

        # Find cells in target clusters
        obs_clusters = adata.obs[cluster_col].astype(str)
        mask = obs_clusters.isin(target_clusters)
        if target_cells:
            mask = mask | adata.obs_names.isin(target_cells)

        matching_idx = np.where(mask.values)[0]
        if len(matching_idx) == 0:
            adata.file.close()
            continue

        # Subsample per cluster within this file
        selected_idx = []
        matching_clusters = obs_clusters.iloc[matching_idx].values
        for clust in np.unique(matching_clusters):
            already = cluster_cell_counts.get(clust, 0)
            remaining = max_per_cluster - already
            if remaining <= 0:
                continue
            clust_positions = matching_idx[matching_clusters == clust]
            if len(clust_positions) > remaining:
                clust_positions = rng.choice(clust_positions, remaining, replace=False)
            selected_idx.extend(clust_positions.tolist())
            cluster_cell_counts[clust] = already + len(clust_positions)

        if not selected_idx:
            adata.file.close()
            continue

        selected_idx = sorted(selected_idx)
        log.info("  Loading %d human cells from %s ...", len(selected_idx), h5ad_path.name)
        chunk = adata[selected_idx].X
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        else:
            chunk = np.asarray(chunk)

        selected_cell_ids = adata.obs_names[selected_idx].tolist()
        selected_labels = obs_clusters.iloc[selected_idx].values
        for i, cell_id in enumerate(selected_cell_ids):
            rows[cell_id] = chunk[i].astype(np.float32)
            labels[cell_id] = selected_labels[i]

        adata.file.close()

    if not rows:
        return pd.DataFrame(), pd.Series(dtype=str)

    log.info("  Total human cells loaded: %d across %d clusters", len(rows), len(cluster_cell_counts))
    expr_df = pd.DataFrame(rows, index=gene_names).T
    labels_series = pd.Series(labels)
    return expr_df, labels_series


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pass2(
    bridge_dir: Path,
    outdir: Path,
    max_per_cluster: int,
    n_hvgs: int,
    ortholog_path: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------- Load unresolved lists
    unresolved_mouse_path = bridge_dir / "unresolved_mouse.csv"
    unresolved_human_path = bridge_dir / "unresolved_human.csv"

    if not unresolved_mouse_path.exists() or not unresolved_human_path.exists():
        raise FileNotFoundError(
            f"Pass 1 outputs not found in {bridge_dir}. Run 07_cluster_bridge.py first."
        )

    unresolved_mouse = pd.read_csv(unresolved_mouse_path)["mouse_cluster"].tolist()
    unresolved_human = pd.read_csv(unresolved_human_path)["human_cluster"].tolist()

    log.info(
        "Unresolved clusters — mouse: %d, human: %d",
        len(unresolved_mouse),
        len(unresolved_human),
    )

    if not unresolved_mouse or not unresolved_human:
        log.info("No unresolved clusters — nothing to do for Pass 2.")
        _merge_and_save(bridge_dir, outdir, pd.DataFrame())
        return

    # -------------------------------------------------- Ortholog map
    orth_map = load_ortholog_map(ortholog_path)
    human_symbols_set = set(orth_map.values())

    # -------------------------------------------------- Subsample mouse cells
    log.info("Loading mouse metadata ...")
    meta = pd.read_csv(_MOUSE_META_CSV, usecols=["cell_label", "cluster_alias"])
    meta = meta.set_index("cell_label")["cluster_alias"]
    sampled_mouse_cells = subsample_cells(meta, unresolved_mouse, max_per_cluster)

    log.info("Loading mouse expression for %d cells ...", len(sampled_mouse_cells))
    mouse_expr, mouse_genes = _load_cells_from_h5ads(
        _MOUSE_H5AD_FILES, set(sampled_mouse_cells)
    )
    if mouse_expr.empty:
        log.error("Failed to load mouse expression. Aborting Pass 2.")
        return

    # Rename mouse genes to human orthologs
    rename_map = {g: orth_map[g] for g in mouse_expr.columns if g in orth_map}
    mouse_expr = mouse_expr.rename(columns=rename_map)
    mouse_expr = mouse_expr[[c for c in mouse_expr.columns if c in human_symbols_set]]

    mouse_labels = meta.loc[mouse_expr.index]

    # -------------------------------------------------- Subsample human cells
    log.info("Loading human expression for unresolved clusters ...")
    human_expr, human_labels = _load_human_cells_from_h5ads(
        _HUMAN_H5AD_DIR,
        target_cells=set(),
        target_clusters=set(unresolved_human),
    )
    if human_expr.empty:
        log.error("Failed to load human expression. Aborting Pass 2.")
        return

    # Restrict human to ortholog genes
    human_expr = human_expr[[c for c in human_expr.columns if c in human_symbols_set]]

    # Subsample human cells per cluster
    sampled_human_cells = subsample_cells(
        human_labels[human_labels.isin(unresolved_human)],
        unresolved_human,
        max_per_cluster,
    )
    human_expr = human_expr.loc[human_expr.index.isin(sampled_human_cells)]
    human_labels = human_labels.loc[human_expr.index]

    # -------------------------------------------------- Shared genes + HVGs
    shared_genes = mouse_expr.columns.intersection(human_expr.columns).tolist()
    log.info("Shared ortholog genes: %d", len(shared_genes))

    mouse_expr = mouse_expr[shared_genes]
    human_expr = human_expr[shared_genes]

    # Select HVGs by variance
    combined = pd.concat([mouse_expr, human_expr], axis=0)
    gene_var = combined.var(axis=0)
    n_hvgs_actual = min(n_hvgs, len(shared_genes))
    top_genes = gene_var.nlargest(n_hvgs_actual).index.tolist()
    log.info("Using %d HVGs for MetaNeighbor.", len(top_genes))

    mouse_expr = mouse_expr[top_genes]
    human_expr = human_expr[top_genes]

    # -------------------------------------------------- Combine + species labels
    all_cells = list(mouse_expr.index) + list(human_expr.index)
    expr_combined = pd.concat([mouse_expr, human_expr], axis=0)

    species = pd.Series(
        ["mouse"] * len(mouse_expr) + ["human"] * len(human_expr),
        index=all_cells,
    )

    # -------------------------------------------------- MetaNeighbor
    log.info(
        "Running MetaNeighbor AUROC (%d mouse cells x %d human cells) ...",
        len(mouse_expr),
        len(human_expr),
    )
    auroc_matrix = compute_metaneighbor_auroc(
        expr_df=expr_combined,
        labels_mouse=mouse_labels,
        labels_human=human_labels,
        species=species,
    )

    auroc_out = outdir / "metaneighbor_auroc_matrix.csv"
    auroc_matrix.to_csv(auroc_out)
    log.info("  AUROC matrix saved: %s  %s", auroc_out, auroc_matrix.shape)

    # -------------------------------------------------- Pass 2 resolved pairs
    # Best human cluster per unresolved mouse cluster (argmax AUROC)
    best_human_idx = auroc_matrix.values.argmax(axis=1)
    pass2_records = []
    for mi, mouse_cl in enumerate(auroc_matrix.index):
        hi = best_human_idx[mi]
        human_cl = auroc_matrix.columns[hi]
        auroc_val = float(auroc_matrix.values[mi, hi])
        pass2_records.append(
            {
                "mouse_cluster": mouse_cl,
                "human_cluster": human_cl,
                "auroc": auroc_val,
                "is_rbh": False,  # MetaNeighbor pass — no reciprocal check here
                "method": "metaneighbor_pass2",
            }
        )
    pass2_df = pd.DataFrame(pass2_records)

    pass2_out = outdir / "pass2_resolved_pairs.csv"
    pass2_df.to_csv(pass2_out, index=False)
    log.info("  Pass 2 pairs saved: %s", pass2_out)

    _merge_and_save(bridge_dir, outdir, pass2_df)


def _merge_and_save(bridge_dir: Path, outdir: Path, pass2_df: pd.DataFrame) -> None:
    """Merge Pass 1 RBH pairs with Pass 2 MetaNeighbor pairs."""
    pass1_path = bridge_dir / "rbh_resolved_pairs.csv"
    if not pass1_path.exists():
        log.warning("Pass 1 resolved pairs not found at %s — saving Pass 2 only.", pass1_path)
        all_df = pass2_df
    else:
        pass1_df = pd.read_csv(pass1_path)
        # Align columns so concat works cleanly
        all_df = pd.concat(
            [pass1_df, pass2_df],
            axis=0,
            ignore_index=True,
            sort=False,
        )

    out_path = outdir / "all_cluster_correspondences.csv"
    all_df.to_csv(out_path, index=False)
    log.info("All cluster correspondences saved: %s  (%d rows)", out_path, len(all_df))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="MetaNeighbor Pass 2: resolve unresolved clusters from Pass 1."
    )
    parser.add_argument(
        "--bridge-dir",
        type=Path,
        default=_DEFAULT_BRIDGE_DIR,
        help="Directory containing Pass 1 outputs (default: results/cluster_bridge/).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help="Output directory (default: same as --bridge-dir).",
    )
    parser.add_argument(
        "--max-per-cluster",
        type=int,
        default=200,
        help="Max cells to subsample per cluster (default: 200).",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=3000,
        help="Number of HVGs for MetaNeighbor (default: 3000).",
    )
    parser.add_argument(
        "--ortholog-map",
        type=Path,
        default=_DEFAULT_ORTHOLOG,
        help="Ortholog mapping CSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_pass2(
        bridge_dir=args.bridge_dir,
        outdir=args.outdir,
        max_per_cluster=args.max_per_cluster,
        n_hvgs=args.n_hvgs,
        ortholog_path=args.ortholog_map,
    )


if __name__ == "__main__":
    main()
