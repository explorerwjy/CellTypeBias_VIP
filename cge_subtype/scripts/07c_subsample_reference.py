#!/usr/bin/env python3
"""
Step 7c: Create subsampled reference h5ad for Module B's Harmony mapping.

Loads the Allen WMB metadata, subsamples up to ``--max-per-cluster`` cells per
cluster, then reads those cells from the two Isocortex h5ad files in backed
mode.  Adds ``cluster_alias``, ``subclass``, and ``class`` columns to
adata.obs, then writes the result to:

  <outdir>/reference_subsampled.h5ad

This compact reference is the input to the Harmony-based cross-species mapping
in Module B.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

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
_META_COLS = ["cell_label", "cluster_alias", "subclass_label", "class_label"]

_DEFAULT_OUTDIR = PROJECT_ROOT / "cge_subtype" / "results" / "preprocessed"
CHUNK_SIZE = 1000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def subsample_cell_ids(
    meta: pd.DataFrame,
    max_per_cluster: int,
    seed: int = 42,
) -> list[str]:
    """Return a subsampled list of cell IDs (max_per_cluster per cluster_alias)."""
    rng = np.random.default_rng(seed)
    sampled: list[str] = []
    for cluster, grp in meta.groupby("cluster_alias"):
        cells = grp["cell_label"].tolist()
        if len(cells) > max_per_cluster:
            cells = rng.choice(cells, size=max_per_cluster, replace=False).tolist()
        sampled.extend(cells)
    return sampled


def load_cells_from_h5ads(
    h5ad_files: list[Path],
    target_cells: set[str],
) -> tuple[ad.AnnData | None, list[str]]:
    """Load expression rows for ``target_cells`` from backed h5ad files.

    Returns (concatenated AnnData, ordered cell IDs actually loaded).
    """
    parts: list[ad.AnnData] = []

    for h5ad_path in h5ad_files:
        if not h5ad_path.exists():
            log.warning("h5ad not found: %s — skipping.", h5ad_path)
            continue

        adata = ad.read_h5ad(h5ad_path, backed="r")
        overlap = [c for c in adata.obs_names if c in target_cells]
        if not overlap:
            log.info("  No target cells in %s — skipping.", h5ad_path.name)
            adata.file.close()
            continue

        log.info("  Reading %d cells from %s ...", len(overlap), h5ad_path.name)
        positions = [adata.obs_names.get_loc(c) for c in overlap]

        # Read in chunks to cap RAM
        n_chunks = (len(positions) + CHUNK_SIZE - 1) // CHUNK_SIZE
        chunk_arrays: list[np.ndarray] = []
        for ci in range(n_chunks):
            pos_chunk = positions[ci * CHUNK_SIZE: (ci + 1) * CHUNK_SIZE]
            chunk_x = adata[pos_chunk].X
            if hasattr(chunk_x, "toarray"):
                chunk_x = chunk_x.toarray()
            else:
                chunk_x = np.asarray(chunk_x)
            chunk_arrays.append(chunk_x.astype(np.float32))

        x_full = np.vstack(chunk_arrays)  # (n_cells_this_file, n_genes)

        part = ad.AnnData(
            X=sparse.csr_matrix(x_full),
            obs=pd.DataFrame(index=overlap),
            var=adata.var.copy(),
        )
        parts.append(part)
        adata.file.close()
        log.info("    Done. Part shape: %s", part.shape)

    if not parts:
        return None, []

    if len(parts) == 1:
        combined = parts[0]
    else:
        combined = ad.concat(parts, axis=0, join="inner", merge="same")

    return combined, combined.obs_names.tolist()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def create_subsampled_reference(outdir: Path, max_per_cluster: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ 1. Metadata
    log.info("Loading mouse metadata: %s", _MOUSE_META_CSV)
    # Read only the columns we need; gracefully handle missing columns
    avail_cols = pd.read_csv(_MOUSE_META_CSV, nrows=0).columns.tolist()
    read_cols = [c for c in _META_COLS if c in avail_cols]
    missing = [c for c in _META_COLS if c not in avail_cols]
    if missing:
        log.warning("Metadata columns not found (will be NaN): %s", missing)

    meta = pd.read_csv(_MOUSE_META_CSV, usecols=read_cols)
    log.info("  Total cells in metadata: %d", len(meta))

    # Drop rows without a cluster_alias
    meta = meta.dropna(subset=["cluster_alias"])
    n_clusters = meta["cluster_alias"].nunique()
    log.info("  Clusters: %d", n_clusters)

    # ------------------------------------------------------------------ 2. Subsample
    log.info("Subsampling up to %d cells per cluster ...", max_per_cluster)
    sampled_ids = subsample_cell_ids(meta, max_per_cluster)
    log.info("  Subsampled cell count: %d", len(sampled_ids))

    # ------------------------------------------------------------------ 3. Load expression
    log.info("Loading expression for subsampled cells from h5ad files ...")
    adata, loaded_ids = load_cells_from_h5ads(_MOUSE_H5AD_FILES, set(sampled_ids))

    if adata is None or adata.n_obs == 0:
        log.error("No cells loaded. Check that h5ad files exist at expected paths.")
        return

    log.info("  Loaded adata shape: %s", adata.shape)

    # ------------------------------------------------------------------ 4. Add metadata to obs
    # Index meta by cell_label
    meta_idx = meta.set_index("cell_label")

    # cluster_alias
    if "cluster_alias" in meta_idx.columns:
        adata.obs["cluster_alias"] = meta_idx.reindex(adata.obs_names)["cluster_alias"].values

    # subclass (may be stored as subclass_label or subclass)
    for src_col, dst_col in [("subclass_label", "subclass"), ("subclass", "subclass")]:
        if src_col in meta_idx.columns:
            adata.obs["subclass"] = meta_idx.reindex(adata.obs_names)[src_col].values
            break

    # class (may be class_label or class)
    for src_col, dst_col in [("class_label", "class"), ("class", "class")]:
        if src_col in meta_idx.columns:
            adata.obs["class"] = meta_idx.reindex(adata.obs_names)[src_col].values
            break

    log.info(
        "  obs columns added: %s",
        [c for c in ("cluster_alias", "subclass", "class") if c in adata.obs.columns],
    )

    # ------------------------------------------------------------------ 5. Save
    out_path = outdir / "reference_subsampled.h5ad"
    log.info("Saving subsampled reference to %s ...", out_path)
    adata.write_h5ad(out_path)
    log.info(
        "  Done. Final shape: %d cells x %d genes.",
        adata.n_obs,
        adata.n_vars,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create subsampled Allen WMB reference h5ad for Harmony mapping."
    )
    parser.add_argument(
        "--max-per-cluster",
        type=int,
        default=200,
        help="Maximum cells to sample per cluster_alias (default: 200).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help="Output directory (default: results/preprocessed/).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    create_subsampled_reference(outdir=args.outdir, max_per_cluster=args.max_per_cluster)


if __name__ == "__main__":
    main()
