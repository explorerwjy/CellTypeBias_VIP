#!/usr/bin/env python3
"""
Recompute MetaNeighbor AUROC matrix for the cross-species CGE comparison
using the canonical Bakken 2021 implementation.

Inputs
------
- Mouse: Allen WMB-10Xv3 cortex h5ad files (Isocortex-1, Isocortex-2). Filter
  to GABAergic CGE/MGE classes; subsample to ~max_per_cluster cells per
  Allen cluster_alias.
- Human: Siletti et al. CGE supercluster h5ad. Filter to clusters 276-296.
  Subsample to ~max_per_cluster cells per cluster.
- Ortholog map: mouse symbol → human symbol (one-to-one).

Output
------
cge_subtype/results/cluster_bridge/metaneighbor_auroc_cge.csv
  AUROC matrix indexed by mouse cluster_alias (rows) × human cluster_id (cols).

Algorithm: see cge_subtype/src/cluster_correspondence.py
::compute_best_hits_metaneighbor (canonical Bakken 2021 port).
"""

import argparse
import gc
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from cge_subtype.src.cluster_correspondence import compute_best_hits_metaneighbor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
MOUSE_H5AD_FILES = [
    Path(
        "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices"
        "/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-raw.h5ad"
    ),
    Path(
        "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices"
        "/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-2-raw.h5ad"
    ),
]
MOUSE_META_CSV = Path(
    "/mnt/data0/AllenMouseSC/abc_download_root/metadata"
    "/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv"
)
HUMAN_H5AD_DIR = Path("/mnt/data0/HumanBrainCellType/SuperTypeRawDat")
ORTHOLOG_PATH = Path(
    "/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results"
    "/cross_species/orthologs/ortholog_mapping.csv"
)

DEFAULT_OUT = PROJECT_ROOT / "cge_subtype" / "results" / "cluster_bridge" / "metaneighbor_auroc_cge.csv"

# Mouse classes to include (CGE + MGE GABAergic from cortex)
MOUSE_GABA_CLASSES = ["06 CTX-CGE GABA", "07 CTX-MGE GABA"]

# Human CGE clusters
HUMAN_CGE_CLUSTERS = [str(c) for c in range(276, 297)]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_ortholog_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    one2one = df[df["is_one_to_one"].astype(bool)]
    return dict(zip(one2one["mouse_symbol"], one2one["human_symbol"]))


def _get_gene_names(adata):
    for col in ("Gene", "gene_symbol", "gene_name"):
        if col in adata.var.columns:
            return adata.var[col].tolist()
    return adata.var_names.tolist()


def load_mouse_cells(
    max_per_cluster: int,
    seed: int,
    ortholog_map: dict[str, str],
) -> tuple[pd.DataFrame, pd.Series]:
    """Load mouse expression for GABAergic cortical cells, subsampled per cluster_alias.

    Returns (cells × human-symbol-genes DataFrame, cluster_alias Series).
    """
    log.info("Loading mouse metadata ...")
    meta = pd.read_csv(
        MOUSE_META_CSV,
        usecols=["cell_label", "cluster_alias", "subclass", "class", "dataset_label"],
        low_memory=False,
    )
    meta = meta[meta["dataset_label"] == "WMB-10Xv3"]
    meta = meta[meta["class"].isin(MOUSE_GABA_CLASSES)]
    meta = meta.set_index("cell_label")
    log.info(f"  GABAergic mouse cells available: {len(meta):,}")

    # Subsample per cluster_alias
    rng = np.random.default_rng(seed)
    sampled_cells = []
    for clust, group in meta.groupby("cluster_alias"):
        if len(group) > max_per_cluster:
            chosen = rng.choice(group.index.values, max_per_cluster, replace=False)
        else:
            chosen = group.index.values
        sampled_cells.extend(chosen.tolist())
    sampled_cells = set(sampled_cells)
    log.info(
        f"  After subsampling (max {max_per_cluster}/cluster): "
        f"{len(sampled_cells):,} cells across {meta['cluster_alias'].nunique()} clusters"
    )

    # Load expression from h5ad files
    rows: dict[str, np.ndarray] = {}
    gene_names: list[str] | None = None

    for h5ad_path in MOUSE_H5AD_FILES:
        if not h5ad_path.exists():
            log.warning(f"  {h5ad_path} not found, skipping.")
            continue
        log.info(f"  Loading {h5ad_path.name} ...")
        adata = ad.read_h5ad(h5ad_path, backed="r")
        if gene_names is None:
            gene_names = _get_gene_names(adata)
            log.info(f"    {len(gene_names)} genes")
        # Find matching cells
        obs_names = adata.obs_names.values
        keep_mask = np.array([c in sampled_cells for c in obs_names])
        keep_idx = np.where(keep_mask)[0]
        if len(keep_idx) == 0:
            adata.file.close()
            continue
        log.info(f"    {len(keep_idx):,} matching cells")
        chunk = adata[keep_idx].X
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        else:
            chunk = np.asarray(chunk)
        for i, idx in enumerate(keep_idx):
            cell_id = obs_names[idx]
            rows[cell_id] = chunk[i].astype(np.float32)
        adata.file.close()
        del adata
        gc.collect()

    if not rows:
        raise RuntimeError("No mouse expression loaded")

    expr_df = pd.DataFrame(rows, index=gene_names).T  # cells × genes
    log.info(f"  Mouse expression: {expr_df.shape}")

    # Map mouse gene symbols → human symbols (one-to-one orthologs only)
    rename = {g: ortholog_map[g] for g in expr_df.columns if g in ortholog_map}
    expr_df = expr_df.rename(columns=rename)
    # Drop columns that didn't map (still mouse symbols)
    human_symbols_set = set(ortholog_map.values())
    expr_df = expr_df[[c for c in expr_df.columns if c in human_symbols_set]]
    # Average duplicate columns (multiple mouse genes can map to same human symbol
    # if the ortholog map has many-to-one entries)
    if expr_df.columns.duplicated().any():
        n_dup = expr_df.columns.duplicated().sum()
        log.info(f"  Averaging {n_dup} duplicate human symbol columns ...")
        expr_df = expr_df.T.groupby(level=0).mean().T
    log.info(f"  After ortholog renaming: {expr_df.shape}")

    cluster_labels = meta.loc[expr_df.index, "cluster_alias"]
    return expr_df, cluster_labels


def load_human_cge_cells(
    max_per_cluster: int,
    seed: int,
    ortholog_human_genes: set[str],
) -> tuple[pd.DataFrame, pd.Series]:
    """Load human CGE supercluster cells, subsampled to max_per_cluster per cluster_id."""
    cge_h5ad = HUMAN_H5AD_DIR / "Supercluster_CGE_interneuron.h5ad"
    if not cge_h5ad.exists():
        raise FileNotFoundError(f"Human CGE h5ad not found: {cge_h5ad}")

    log.info(f"Loading human CGE cells from {cge_h5ad.name} ...")
    adata = ad.read_h5ad(cge_h5ad, backed="r")
    log.info(f"  CGE supercluster total: {adata.n_obs:,}")

    # Find cluster column
    cluster_col = None
    for col in ("cluster_id", "cluster"):
        if col in adata.obs.columns:
            cluster_col = col
            break
    if cluster_col is None:
        adata.file.close()
        raise RuntimeError("No cluster_id column in human CGE h5ad")

    cluster_str = adata.obs[cluster_col].astype(str)
    target_clusters = set(HUMAN_CGE_CLUSTERS)
    keep_mask = cluster_str.isin(target_clusters)
    log.info(f"  Cells in target CGE clusters (276-296): {keep_mask.sum():,}")

    # Subsample per cluster
    rng = np.random.default_rng(seed)
    sampled_idx_list = []
    for clust in HUMAN_CGE_CLUSTERS:
        c_mask = (cluster_str == clust).values & keep_mask.values
        c_pos = np.where(c_mask)[0]
        if len(c_pos) == 0:
            log.warning(f"  Cluster {clust} has 0 cells")
            continue
        if len(c_pos) > max_per_cluster:
            chosen = rng.choice(c_pos, max_per_cluster, replace=False)
        else:
            chosen = c_pos
        sampled_idx_list.append(chosen)
    sampled_idx = np.sort(np.concatenate(sampled_idx_list))
    log.info(f"  After subsampling: {len(sampled_idx):,} cells")

    # Load expression
    chunk = adata[sampled_idx].X
    if hasattr(chunk, "toarray"):
        chunk = chunk.toarray()
    else:
        chunk = np.asarray(chunk)

    cell_ids = adata.obs_names[sampled_idx].tolist()
    cluster_labels = cluster_str.iloc[sampled_idx].values

    # Gene names — Siletti uses Ensembl IDs as index, gene symbols in 'Gene'
    gene_col = None
    for col in ("Gene", "gene_symbol", "gene_name"):
        if col in adata.var.columns:
            gene_col = col
            break
    if gene_col is not None:
        gene_names = adata.var[gene_col].tolist()
    else:
        gene_names = adata.var_names.tolist()
    log.info(f"  Gene column: {gene_col}, {len(gene_names)} genes")

    adata.file.close()

    expr_df = pd.DataFrame(chunk, index=cell_ids, columns=gene_names)

    # Filter to one-to-one ortholog human symbols
    expr_df = expr_df.loc[:, [c for c in expr_df.columns if c in ortholog_human_genes]]
    # Average duplicate columns (in case multiple Ensembl IDs map to the same symbol)
    if expr_df.columns.duplicated().any():
        expr_df = expr_df.T.groupby(level=0).mean().T

    log.info(f"  Human expression: {expr_df.shape}")
    cluster_series = pd.Series(cluster_labels, index=cell_ids)
    return expr_df, cluster_series


def select_hvgs(combined: pd.DataFrame, n_hvgs: int) -> list[str]:
    """Top-N highest-variance genes (computed across all combined cells)."""
    var = combined.var(axis=0)
    return var.nlargest(n_hvgs).index.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-per-cluster", type=int, default=50)
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    log.info("=" * 60)
    log.info("Recompute MetaNeighbor AUROC for cross-species CGE")
    log.info(f"  max_per_cluster: {args.max_per_cluster}")
    log.info(f"  n_hvgs:          {args.n_hvgs}")
    log.info(f"  out:             {args.out}")
    log.info("=" * 60)

    log.info("Loading ortholog map ...")
    ortholog_map = load_ortholog_map(ORTHOLOG_PATH)
    human_symbols_set = set(ortholog_map.values())
    log.info(f"  {len(ortholog_map)} one-to-one orthologs")

    mouse_expr, mouse_labels = load_mouse_cells(
        args.max_per_cluster, args.seed, ortholog_map
    )
    human_expr, human_labels = load_human_cge_cells(
        args.max_per_cluster, args.seed, human_symbols_set
    )

    # Restrict to shared genes
    shared = mouse_expr.columns.intersection(human_expr.columns)
    log.info(f"Shared ortholog genes: {len(shared)}")
    mouse_expr = mouse_expr[shared]
    human_expr = human_expr[shared]

    # Concatenate and select HVGs
    log.info("Selecting HVGs ...")
    combined = pd.concat([mouse_expr, human_expr], axis=0)
    hvgs = select_hvgs(combined, args.n_hvgs)
    log.info(f"  HVGs: {len(hvgs)}")
    combined = combined[hvgs]

    # Build inputs for compute_best_hits_metaneighbor
    # Function expects (n_genes, n_cells); transpose
    cluster_labels = pd.concat(
        [mouse_labels.astype(str), human_labels.astype(str)]
    )
    cluster_labels = cluster_labels.loc[combined.index]
    study_ids = pd.Series(
        ["mouse"] * len(mouse_expr) + ["human"] * len(human_expr),
        index=mouse_expr.index.tolist() + human_expr.index.tolist(),
    )
    study_ids = study_ids.loc[combined.index]

    log.info(
        f"Running compute_best_hits_metaneighbor on "
        f"{len(combined)} cells × {combined.shape[1]} genes ..."
    )
    expr_genes_x_cells = combined.values.T.astype(float)
    full_auroc = compute_best_hits_metaneighbor(
        expr_genes_x_cells,
        cluster_labels.values,
        study_ids.values,
    )
    log.info(f"  Full AUROC: {full_auroc.shape}")

    # Save the FULL square AUROC matrix (both directions) for downstream analysis
    full_out = args.out.parent / "metaneighbor_auroc_cge_full.csv"
    full_auroc.to_csv(full_out)
    log.info(f"Saved full square matrix: {full_out}")

    # ----- Build the "best mouse per human cluster" matrix -----
    #
    # Bakken's compute_best_hits returns an asymmetric matrix:
    #   result[A, B] = "in the test study where A lives, do A-cluster cells get
    #                   higher votes for voter group B than other test cells"
    #
    # Two directions are possible for cross-species correspondence:
    #   Direction A: result[mouse|X, human|Y]
    #     mouse-as-test: do mouse-X cells get higher votes for human-Y centroid
    #     than other mouse cells? Biased by mouse clusters whose expression is
    #     broadly similar to all human cells (they win against any human Y).
    #   Direction B: result[human|Y, mouse|X]
    #     human-as-test: do human-Y cells get higher votes for mouse-X centroid
    #     than other human cells? Specific to "which mouse centroid most
    #     uniquely identifies this human cluster".
    #
    # We use Direction B as the primary view because it directly answers the
    # question "for each human CGE cluster, which mouse cluster best matches?".
    # The matrix is then transposed to (rows=mouse, cols=human) so the file
    # format matches the existing downstream pipeline (idxmax over rows per col).
    mouse_groups = [g for g in full_auroc.index if g.startswith("mouse|")]
    human_groups = [g for g in full_auroc.index if g.startswith("human|")]

    dir_b = full_auroc.loc[human_groups, mouse_groups].copy()  # rows=human, cols=mouse
    dir_b_T = dir_b.T.copy()  # rows=mouse, cols=human (matches downstream format)

    # Strip prefixes
    dir_b_T.index = [r.split("|", 1)[1] for r in dir_b_T.index]
    dir_b_T.columns = [c.split("|", 1)[1] for c in dir_b_T.columns]

    # Sort: rows by mouse cluster_alias (numeric), columns by human cluster_id (numeric)
    dir_b_T.index = pd.Index([float(x) for x in dir_b_T.index])
    dir_b_T = dir_b_T.sort_index()
    dir_b_T.columns = [str(c) for c in dir_b_T.columns]
    dir_b_T = dir_b_T[sorted(dir_b_T.columns, key=int)]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    dir_b_T.to_csv(args.out)
    log.info(f"Saved Direction B (transposed) matrix: {args.out}  shape={dir_b_T.shape}")

    # Quick summary
    print()
    print("=" * 60)
    print("Top mouse cluster matches per human CGE cluster (canonical Bakken Direction B)")
    print("=" * 60)
    for c in sorted(dir_b_T.columns, key=int):
        best = dir_b_T[c].idxmax()
        print(f"  {c}: best mouse cluster_alias={best}, AUROC={dir_b_T[c].max():.3f}")


if __name__ == "__main__":
    main()
