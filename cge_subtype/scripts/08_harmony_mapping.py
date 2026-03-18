#!/usr/bin/env python3
"""
Step 08: Single-run Harmony mapping of a patch-seq dataset to the Allen WMB reference.

Loads the subsampled reference (results/preprocessed/reference_subsampled.h5ad),
loads a patch-seq query dataset (M1 or V1), runs Harmony-based label transfer via
run_harmony_mapping(), evaluates with evaluate_mapping(), and saves predictions,
metrics, and the confusion matrix as CSV files tagged with the hyperparameter values.

Usage
-----
    python 08_harmony_mapping.py --dataset M1 --n-pcs 50 --n-hvgs 3000 \
        --theta 2.0 --n-neighbors 30 --outdir results/harmony/

    python 08_harmony_mapping.py --dataset V1 --neuronal-only \
        --n-pcs 50 --n-hvgs 3000 --theta 2.0 --n-neighbors 30 \
        --outdir results/harmony/
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Project root on sys.path so cge_subtype.src is importable
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from cge_subtype.src.harmony_mapping import run_harmony_mapping, evaluate_mapping  # noqa: E402

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Dataset configurations
# ---------------------------------------------------------------------------
DATASET_CONFIGS: dict[str, dict] = {
    "M1": {
        "counts_path": Path(
            "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_counts.csv.gz"
        ),
        "meta_path": Path(
            "/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv"
        ),
        "meta_sep": "\t",
        "ground_truth_col": "RNA family",
    },
    "V1": {
        "counts_path": Path(
            "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_counts.csv"
        ),
        "meta_path": Path(
            "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_metadata.csv"
        ),
        "meta_sep": ",",
        "ground_truth_col": "corresponding_AIT2.3.1_alias",
    },
}

# ---------------------------------------------------------------------------
# Default reference path
# ---------------------------------------------------------------------------
_DEFAULT_REF = PROJECT_ROOT / "cge_subtype" / "results" / "preprocessed" / "reference_subsampled.h5ad"

# Column in reference obs that holds the cluster labels to transfer
_REF_CLUSTER_COL = "cluster_alias"

# Column in reference obs that holds the broad class annotation
_REF_CLASS_COL = "class"

# Patterns identifying neuronal cells in the class column (case-insensitive)
_NEURONAL_PATTERNS = re.compile(r"glut|gaba", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_query_adata(dataset: str) -> tuple[ad.AnnData, pd.DataFrame, str]:
    """Load query counts + metadata for *dataset*.

    Returns
    -------
    query_adata : AnnData  (cells x genes, raw counts)
    metadata    : DataFrame indexed by cell ID
    gt_col      : ground-truth column name in metadata
    """
    cfg = DATASET_CONFIGS[dataset]
    counts_path: Path = cfg["counts_path"]
    meta_path: Path = cfg["meta_path"]
    sep: str = cfg["meta_sep"]
    gt_col: str = cfg["ground_truth_col"]

    log.info("Loading %s counts from %s …", dataset, counts_path)
    counts_df = pd.read_csv(counts_path, index_col=0)
    # counts_df: genes x cells  OR  cells x genes — detect by shape heuristic
    # Convention: if #rows >> #cols we assume genes×cells and transpose
    if counts_df.shape[0] > counts_df.shape[1]:
        log.info(
            "  Counts shape %s looks like genes×cells — transposing.", counts_df.shape
        )
        counts_df = counts_df.T

    log.info("  Counts shape (cells × genes): %s", counts_df.shape)

    query_adata = ad.AnnData(
        X=counts_df.values.astype(np.float32),
        obs=pd.DataFrame(index=counts_df.index),
        var=pd.DataFrame(index=counts_df.columns),
    )

    log.info("Loading %s metadata from %s …", dataset, meta_path)
    meta = pd.read_csv(meta_path, sep=sep, index_col=0)
    log.info("  Metadata shape: %s", meta.shape)

    return query_adata, meta, gt_col


def filter_neuronal(ref_adata: ad.AnnData, class_col: str) -> ad.AnnData:
    """Return a copy of *ref_adata* retaining only neuronal cells."""
    if class_col not in ref_adata.obs.columns:
        log.warning(
            "Class column '%s' not found in reference obs; skipping neuronal filter.",
            class_col,
        )
        return ref_adata

    mask = ref_adata.obs[class_col].astype(str).str.contains(
        _NEURONAL_PATTERNS, regex=True
    )
    n_before = ref_adata.n_obs
    filtered = ref_adata[mask].copy()
    log.info(
        "Neuronal filter: %d → %d cells (removed %d non-neuronal).",
        n_before,
        filtered.n_obs,
        n_before - filtered.n_obs,
    )
    return filtered


def build_output_tag(
    dataset: str,
    n_pcs: int,
    n_hvgs: int,
    theta: float,
    n_neighbors: int,
    neuronal_only: bool,
) -> str:
    """Build a filename tag encoding the hyperparameter values."""
    tag = f"{dataset}_pcs{n_pcs}_hvg{n_hvgs}_theta{theta}_k{n_neighbors}"
    if neuronal_only:
        tag += "_neuronly"
    return tag


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    tag = build_output_tag(
        dataset=args.dataset,
        n_pcs=args.n_pcs,
        n_hvgs=args.n_hvgs,
        theta=args.theta,
        n_neighbors=args.n_neighbors,
        neuronal_only=args.neuronal_only,
    )
    log.info("Run tag: %s", tag)

    # ------------------------------------------------------------------
    # 1. Load reference
    # ------------------------------------------------------------------
    ref_path: Path = args.ref_path
    log.info("Loading reference from %s …", ref_path)
    ref_adata = ad.read_h5ad(ref_path)
    log.info("  Reference shape: %s", ref_adata.shape)

    # If var_names are Ensembl IDs but gene_symbol column exists, switch to gene symbols
    if "gene_symbol" in ref_adata.var.columns:
        symbols = ref_adata.var["gene_symbol"].astype(str)
        # Drop genes with missing/duplicate symbols
        valid = (symbols != "") & (symbols != "nan") & (~symbols.duplicated(keep="first"))
        ref_adata = ref_adata[:, valid].copy()
        ref_adata.var_names = ref_adata.var["gene_symbol"][valid].values
        ref_adata.var_names_make_unique()
        log.info("  Converted var_names to gene symbols: %d genes", ref_adata.n_vars)

    # Subsample per cluster if requested
    if args.max_per_cluster > 0:
        from cge_subtype.src.harmony_mapping import subsample_reference  # noqa: E402
        log.info("Subsampling reference to max %d cells per cluster …", args.max_per_cluster)
        ref_adata = subsample_reference(
            ref_adata,
            cluster_col=_REF_CLUSTER_COL,
            max_per_cluster=args.max_per_cluster,
        )
        log.info("  After subsampling: %s", ref_adata.shape)

    # Optionally restrict to neuronal cells
    if args.neuronal_only:
        ref_adata = filter_neuronal(ref_adata, class_col=_REF_CLASS_COL)

    if ref_adata.n_obs == 0:
        log.error("Reference is empty after filtering. Aborting.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 2. Load query
    # ------------------------------------------------------------------
    query_adata, metadata, gt_col = load_query_adata(args.dataset)

    # ------------------------------------------------------------------
    # 3. Run Harmony mapping
    # ------------------------------------------------------------------
    log.info(
        "Running Harmony mapping (n_pcs=%d, n_hvgs=%d, theta=%.1f, k=%d) …",
        args.n_pcs,
        args.n_hvgs,
        args.theta,
        args.n_neighbors,
    )
    predictions = run_harmony_mapping(
        ref_adata=ref_adata,
        query_adata=query_adata,
        cluster_col=_REF_CLUSTER_COL,
        n_pcs=args.n_pcs,
        n_hvgs=args.n_hvgs,
        theta=args.theta,
        n_neighbors=args.n_neighbors,
    )
    log.info("  Predictions: %d cells.", len(predictions))

    # ------------------------------------------------------------------
    # 4. Evaluate
    # ------------------------------------------------------------------
    log.info("Evaluating predictions against ground truth column '%s' …", gt_col)
    eval_results = evaluate_mapping(
        predictions=predictions,
        ground_truth_col=gt_col,
        metadata=metadata,
    )
    log.info(
        "  Overall accuracy: %.4f | Non-neuronal rate: %.4f | N cells: %d",
        eval_results["overall_accuracy"],
        eval_results["non_neuronal_rate"],
        eval_results["n_cells"],
    )

    # ------------------------------------------------------------------
    # 5. Save outputs
    # ------------------------------------------------------------------
    # 5a. Predictions CSV
    pred_path = outdir / f"{tag}_predictions.csv"
    predictions.to_csv(pred_path, index=False)
    log.info("Saved predictions → %s", pred_path)

    # 5b. Metrics CSV
    metrics_df = pd.DataFrame(
        [
            {
                "dataset": args.dataset,
                "n_pcs": args.n_pcs,
                "n_hvgs": args.n_hvgs,
                "theta": args.theta,
                "n_neighbors": args.n_neighbors,
                "neuronal_only": args.neuronal_only,
                "overall_accuracy": eval_results["overall_accuracy"],
                "non_neuronal_rate": eval_results["non_neuronal_rate"],
                "n_cells": eval_results["n_cells"],
                "tag": tag,
            }
        ]
    )
    metrics_path = outdir / f"{tag}_metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)
    log.info("Saved metrics → %s", metrics_path)

    # 5c. Confusion matrix CSV
    confusion_path = outdir / f"{tag}_confusion.csv"
    eval_results["confusion_matrix"].to_csv(confusion_path)
    log.info("Saved confusion matrix → %s", confusion_path)

    log.info("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Single-run Harmony mapping of a patch-seq dataset to the WMB reference."
    )
    parser.add_argument(
        "--dataset",
        choices=["M1", "V1"],
        required=True,
        help="Patch-seq dataset to map (M1 or V1).",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of PCA components (default: 50).",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=3000,
        help="Number of highly variable genes (default: 3000).",
    )
    parser.add_argument(
        "--theta",
        type=float,
        default=2.0,
        help="Harmony diversity penalty theta (default: 2.0).",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of kNN neighbors for label transfer (default: 30).",
    )
    parser.add_argument(
        "--max-per-cluster",
        type=int,
        default=0,
        help=(
            "If > 0, subsample reference to this many cells per cluster before mapping "
            "(default: 0 = no additional subsampling, use the file as-is)."
        ),
    )
    parser.add_argument(
        "--neuronal-only",
        action="store_true",
        default=False,
        help=(
            "If set, filter reference to neuronal classes only "
            "(match 'Glut' or 'GABA' in the class column, case-insensitive)."
        ),
    )
    parser.add_argument(
        "--ref-path",
        type=Path,
        default=_DEFAULT_REF,
        help=f"Path to subsampled reference h5ad (default: {_DEFAULT_REF}).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=PROJECT_ROOT / "cge_subtype" / "results" / "harmony",
        help="Output directory for predictions / metrics / confusion matrix CSVs.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
