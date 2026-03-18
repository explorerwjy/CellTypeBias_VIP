#!/usr/bin/env python3
"""
Step 6: Compute pseudobulk centroids for mouse and/or human atlases.

Loads each atlas in backed mode and processes cells in chunks of 1000 to keep
peak RAM manageable.  Outputs one CSV per species: clusters x genes.

Mouse data
----------
  Files   : /mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/
              WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-{1,2}-raw.h5ad
  Metadata: /mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/
              20230830/views/cell_metadata_with_cluster_annotation.csv
  Cluster column : cluster_alias
  Cell ID column : cell_label

Human data
----------
  Files   : /mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_*.h5ad
  Cluster column (tried in order): cluster_id, cluster, supertype, cell_type
  Gene name (tried in order):      Gene, gene_symbol, gene_name, var_names

Outputs
-------
  <outdir>/mouse_pseudobulk.csv
  <outdir>/human_pseudobulk.csv
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
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
HUMAN_H5AD_GLOB = (
    "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_*.h5ad"
)

CHUNK_SIZE = 1000  # cells per chunk to cap RAM


# ---------------------------------------------------------------------------
# Helper: resolve gene names from adata.var
# ---------------------------------------------------------------------------
def _get_gene_names(adata: ad.AnnData) -> list[str]:
    """Return gene symbol list from adata.var, trying common column names."""
    for col in ("Gene", "gene_symbol", "gene_name"):
        if col in adata.var.columns:
            return adata.var[col].tolist()
    # Fall back to var_names
    return adata.var_names.tolist()


# ---------------------------------------------------------------------------
# Pseudobulk: accumulate sum + count in chunks, then divide
# ---------------------------------------------------------------------------
def _accumulate_pseudobulk_chunked(
    adata: ad.AnnData,
    labels: pd.Series,
    chunk_size: int = CHUNK_SIZE,
) -> dict[str, np.ndarray]:
    """Return {cluster: sum_vector} dict, processing adata in row-chunks.

    ``labels`` must be indexed by the same cell IDs as ``adata.obs_names``.
    Only cells present in both ``adata.obs_names`` and ``labels.index`` are used.
    """
    # Align labels to adata obs order; NaN-labelled cells are skipped
    common_cells = adata.obs_names[adata.obs_names.isin(labels.index)]
    if len(common_cells) == 0:
        log.warning("No overlap between adata.obs_names and labels index — skipping.")
        return {}

    n_genes = adata.n_vars
    cluster_sums: dict[str, np.ndarray] = {}
    cluster_counts: dict[str, int] = {}

    n_cells = len(common_cells)
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        chunk_cells = common_cells[start:end]

        # Positional indices into adata for backed mode
        cell_positions = [adata.obs_names.get_loc(c) for c in chunk_cells]
        # Read chunk — backed mode returns a view so we convert to dense
        chunk_data = adata[cell_positions].X
        if hasattr(chunk_data, "toarray"):
            chunk_data = chunk_data.toarray()
        else:
            chunk_data = np.asarray(chunk_data)

        chunk_labels = labels.loc[chunk_cells].values

        for i, cell_label in enumerate(chunk_labels):
            if pd.isna(cell_label):
                continue
            lbl = str(cell_label)
            if lbl not in cluster_sums:
                cluster_sums[lbl] = np.zeros(n_genes, dtype=np.float64)
                cluster_counts[lbl] = 0
            cluster_sums[lbl] += chunk_data[i].astype(np.float64)
            cluster_counts[lbl] += 1

        if (chunk_idx + 1) % 20 == 0 or chunk_idx == n_chunks - 1:
            log.info(
                "  Processed chunk %d/%d (%d/%d cells)",
                chunk_idx + 1, n_chunks, end, n_cells,
            )

    # Divide sums by counts
    mean_dict: dict[str, np.ndarray] = {
        lbl: cluster_sums[lbl] / cluster_counts[lbl]
        for lbl in cluster_sums
    }
    return mean_dict


# ---------------------------------------------------------------------------
# Mouse
# ---------------------------------------------------------------------------
def compute_mouse_pseudobulk(outdir: Path) -> None:
    """Load Allen WMB Isocortex h5ads, compute pseudobulk, save CSV."""
    log.info("=== Mouse pseudobulk ===")

    # Load metadata
    log.info("Loading mouse metadata from %s", MOUSE_META_CSV)
    meta = pd.read_csv(MOUSE_META_CSV, usecols=["cell_label", "cluster_alias"])
    meta = meta.set_index("cell_label")["cluster_alias"].dropna()
    log.info("  %d cells with cluster_alias in metadata", len(meta))

    all_mean_dicts: list[dict[str, np.ndarray]] = []
    gene_names: list[str] | None = None

    for h5ad_path in MOUSE_H5AD_FILES:
        log.info("Loading (backed) %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path, backed="r")
        log.info("  Shape: %s", adata.shape)

        if gene_names is None:
            gene_names = _get_gene_names(adata)

        mean_dict = _accumulate_pseudobulk_chunked(adata, meta)
        all_mean_dicts.append(mean_dict)
        adata.file.close()
        log.info("  Done. %d clusters from this file.", len(mean_dict))

    # Merge across files: clusters may appear in both files
    merged: dict[str, list[np.ndarray]] = {}
    for d in all_mean_dicts:
        for lbl, vec in d.items():
            merged.setdefault(lbl, []).append(vec)

    # Average the per-file means (weighted equally; exact merge would need counts)
    final_means = {lbl: np.mean(vecs, axis=0) for lbl, vecs in merged.items()}

    pb = pd.DataFrame(final_means, index=gene_names).T  # clusters x genes
    pb.index.name = "cluster"

    out_path = outdir / "mouse_pseudobulk.csv"
    pb.to_csv(out_path)
    log.info("Mouse pseudobulk saved: %s  (%d clusters x %d genes)", out_path, *pb.shape)


# ---------------------------------------------------------------------------
# Human
# ---------------------------------------------------------------------------
def compute_human_pseudobulk(outdir: Path) -> None:
    """Load Supercluster h5ads, compute pseudobulk, save CSV."""
    log.info("=== Human pseudobulk ===")

    h5ad_files = sorted(Path("/mnt/data0/HumanBrainCellType/SuperTypeRawDat").glob("Supercluster_*.h5ad"))
    if not h5ad_files:
        log.error("No Supercluster_*.h5ad files found. Check HUMAN_H5AD_GLOB path.")
        return

    log.info("Found %d Supercluster h5ad files", len(h5ad_files))

    all_rows: list[pd.Series] = []
    gene_names: list[str] | None = None
    cluster_col_used: str | None = None

    for h5ad_path in h5ad_files:
        log.info("Loading (backed) %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path, backed="r")
        log.info("  Shape: %s", adata.shape)

        # Determine gene names (deduplicate to avoid DataFrame alignment errors)
        if gene_names is None:
            raw_names = _get_gene_names(adata)
            # Keep first occurrence of each gene name
            seen = set()
            unique_mask = []
            for g in raw_names:
                if g not in seen:
                    seen.add(g)
                    unique_mask.append(True)
                else:
                    unique_mask.append(False)
            gene_names = [g for g, keep in zip(raw_names, unique_mask) if keep]
            _gene_mask = unique_mask  # save mask for subsetting expression chunks
            log.info("  Gene names: %d unique / %d total", len(gene_names), len(raw_names))

        # Determine cluster column
        cluster_col: str | None = None
        for col in ("cluster_id", "cluster", "supertype", "cell_type"):
            if col in adata.obs.columns:
                cluster_col = col
                break
        if cluster_col is None:
            log.warning("  No cluster column found in %s — skipping.", h5ad_path.name)
            adata.file.close()
            continue

        if cluster_col_used is None:
            cluster_col_used = cluster_col
            log.info("  Using cluster column: %s", cluster_col)

        labels = adata.obs[cluster_col].dropna()
        mean_dict = _accumulate_pseudobulk_chunked(adata, labels)
        adata.file.close()

        for lbl, vec in mean_dict.items():
            # Subset to unique genes if there were duplicates
            if len(vec) > len(gene_names):
                vec = vec[np.array(_gene_mask)]
            row = pd.Series(vec, index=gene_names, name=lbl)
            all_rows.append(row)

        log.info("  Done. %d clusters from this file.", len(mean_dict))

    if not all_rows:
        log.error("No human pseudobulk data computed.")
        return

    pb = pd.DataFrame(all_rows)  # clusters x genes
    # Collapse duplicate cluster labels by averaging
    pb = pb.groupby(pb.index).mean()
    pb.index.name = "cluster"

    out_path = outdir / "human_pseudobulk.csv"
    pb.to_csv(out_path)
    log.info("Human pseudobulk saved: %s  (%d clusters x %d genes)", out_path, *pb.shape)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute pseudobulk centroids for mouse and/or human atlases."
    )
    parser.add_argument(
        "--species",
        choices=["mouse", "human", "both"],
        default="both",
        help="Which species to process (default: both).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=PROJECT_ROOT / "cge_subtype" / "results" / "pseudobulk",
        help="Output directory for pseudobulk CSVs (default: cge_subtype/results/pseudobulk/).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", args.outdir)

    if args.species in ("mouse", "both"):
        compute_mouse_pseudobulk(args.outdir)

    if args.species in ("human", "both"):
        compute_human_pseudobulk(args.outdir)

    log.info("Pseudobulk computation complete.")


if __name__ == "__main__":
    main()
