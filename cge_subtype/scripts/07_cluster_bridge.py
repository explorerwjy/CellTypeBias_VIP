#!/usr/bin/env python3
"""
Step 7: Spearman RBH Pass 1 — cross-species cluster bridge.

Loads mouse and human pseudobulk centroids (clusters x genes), restricts to
one-to-one ortholog pairs, selects top HVGs, computes a Spearman correlation
matrix, determines a permutation-based threshold, and identifies reciprocal
best hit (RBH) pairs.

Inputs
------
  --mouse-pb   : mouse_pseudobulk.csv   (clusters x genes, mouse gene symbols)
  --human-pb   : human_pseudobulk.csv   (clusters x genes, human gene symbols)
  --ortholog-map : ortholog_mapping.csv
                   Columns: mouse_symbol, human_symbol, is_one_to_one

Outputs (written to --outdir)
-------------------------------
  spearman_corr_matrix.csv    full (n_mouse x n_human) correlation matrix
  rbh_resolved_pairs.csv      RBH pairs that passed the threshold
  unresolved_mouse.csv        mouse clusters without a confirmed RBH
  unresolved_human.csv        human clusters without a confirmed RBH
  pass1_summary.csv           aggregate statistics
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from cge_subtype.src.cluster_correspondence import (
    compute_spearman_corr_matrix,
    determine_rbh_threshold,
    find_reciprocal_best_hits,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Default paths
_DEFAULT_ORTHOLOG = Path(
    "/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results"
    "/cross_species/orthologs/ortholog_mapping.csv"
)
_DEFAULT_MOUSE_PB = PROJECT_ROOT / "results" / "pseudobulk" / "mouse_pseudobulk.csv"
_DEFAULT_HUMAN_PB = PROJECT_ROOT / "results" / "pseudobulk" / "human_pseudobulk.csv"
_DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "cluster_bridge"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_ortholog_map(ortholog_path: Path) -> dict[str, str]:
    """Load one-to-one ortholog mapping; return {mouse_symbol -> human_symbol}."""
    log.info("Loading ortholog map: %s", ortholog_path)
    df = pd.read_csv(ortholog_path)
    required = {"mouse_symbol", "human_symbol", "is_one_to_one"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ortholog CSV missing columns: {missing}")

    one2one = df[df["is_one_to_one"].astype(bool)]
    mapping = dict(zip(one2one["mouse_symbol"], one2one["human_symbol"]))
    log.info("  %d one-to-one ortholog pairs loaded.", len(mapping))
    return mapping


def select_hvgs(
    mouse_pb: pd.DataFrame,
    human_pb: pd.DataFrame,
    n_hvgs: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select top ``n_hvgs`` highly variable genes by combined variance.

    Genes must already be in a shared namespace (human symbols on both sides).
    """
    shared = mouse_pb.columns.intersection(human_pb.columns)
    log.info("Shared genes (post-ortholog conversion): %d", len(shared))

    mouse_sub = mouse_pb[shared]
    human_sub = human_pb[shared]

    # Variance across clusters in the combined set
    combined = pd.concat([mouse_sub, human_sub], axis=0)
    gene_var = combined.var(axis=0)

    n_hvgs = min(n_hvgs, len(shared))
    top_genes = gene_var.nlargest(n_hvgs).index
    log.info("Selected %d HVGs (requested %d).", len(top_genes), n_hvgs)

    return mouse_sub[top_genes], human_sub[top_genes]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pass1(
    mouse_pb_path: Path,
    human_pb_path: Path,
    ortholog_path: Path,
    n_hvgs: int,
    n_perm: int,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ 1. Load
    log.info("Loading mouse pseudobulk: %s", mouse_pb_path)
    mouse_pb = pd.read_csv(mouse_pb_path, index_col=0)
    log.info("  Mouse: %d clusters x %d genes", *mouse_pb.shape)

    log.info("Loading human pseudobulk: %s", human_pb_path)
    human_pb = pd.read_csv(human_pb_path, index_col=0)
    log.info("  Human: %d clusters x %d genes", *human_pb.shape)

    # ------------------------------------------------- 2. Ortholog conversion
    orth_map = load_ortholog_map(ortholog_path)
    # Rename mouse columns to human symbols (only for mapped genes)
    mouse_mapped_cols = {
        col: orth_map[col] for col in mouse_pb.columns if col in orth_map
    }
    mouse_pb_human = mouse_pb.rename(columns=mouse_mapped_cols)
    # Drop any mouse columns not in the ortholog map
    mouse_pb_human = mouse_pb_human[[c for c in mouse_pb_human.columns if c in orth_map.values()]]
    log.info(
        "Mouse genes after ortholog conversion: %d -> %d",
        mouse_pb.shape[1],
        mouse_pb_human.shape[1],
    )

    # -------------------------------------------------- 3. HVG selection
    mouse_hvg, human_hvg = select_hvgs(mouse_pb_human, human_pb, n_hvgs)

    # -------------------------------------------------- 4. Spearman correlation
    log.info(
        "Computing Spearman correlation matrix (%d mouse x %d human clusters) "
        "using %d genes ...",
        len(mouse_hvg),
        len(human_hvg),
        mouse_hvg.shape[1],
    )
    corr_matrix = compute_spearman_corr_matrix(mouse_hvg, human_hvg)
    log.info("  Correlation matrix: %s", corr_matrix.shape)

    corr_out = outdir / "spearman_corr_matrix.csv"
    corr_matrix.to_csv(corr_out)
    log.info("  Saved: %s", corr_out)

    # -------------------------------------------------- 5. RBH threshold
    log.info("Determining RBH threshold via permutation (n_perm=%d) ...", n_perm)
    threshold = determine_rbh_threshold(corr_matrix, method="permutation", n_perm=n_perm)
    log.info("  RBH threshold (95th pct null): %.4f", threshold)

    # -------------------------------------------------- 6. RBH pairs
    log.info("Finding reciprocal best hits ...")
    rbh_df = find_reciprocal_best_hits(corr_matrix, threshold=threshold)

    resolved = rbh_df[rbh_df["is_rbh"]].copy()
    unresolved_mouse_clusters = rbh_df[~rbh_df["is_rbh"]]["mouse_cluster"].tolist()

    # Human clusters not appearing in resolved pairs
    resolved_human = set(resolved["human_cluster"])
    unresolved_human_clusters = [
        c for c in corr_matrix.columns if c not in resolved_human
    ]

    log.info(
        "  RBH resolved: %d pairs | unresolved mouse: %d | unresolved human: %d",
        len(resolved),
        len(unresolved_mouse_clusters),
        len(unresolved_human_clusters),
    )

    # -------------------------------------------------- 7. Save outputs
    resolved_out = outdir / "rbh_resolved_pairs.csv"
    resolved.to_csv(resolved_out, index=False)
    log.info("  Saved: %s", resolved_out)

    unresolved_mouse_out = outdir / "unresolved_mouse.csv"
    pd.DataFrame({"mouse_cluster": unresolved_mouse_clusters}).to_csv(
        unresolved_mouse_out, index=False
    )
    log.info("  Saved: %s", unresolved_mouse_out)

    unresolved_human_out = outdir / "unresolved_human.csv"
    pd.DataFrame({"human_cluster": unresolved_human_clusters}).to_csv(
        unresolved_human_out, index=False
    )
    log.info("  Saved: %s", unresolved_human_out)

    # Summary
    summary = pd.DataFrame(
        [
            {
                "n_mouse_clusters": corr_matrix.shape[0],
                "n_human_clusters": corr_matrix.shape[1],
                "n_shared_genes_hvg": mouse_hvg.shape[1],
                "rbh_threshold": threshold,
                "n_rbh_resolved": len(resolved),
                "n_unresolved_mouse": len(unresolved_mouse_clusters),
                "n_unresolved_human": len(unresolved_human_clusters),
            }
        ]
    )
    summary_out = outdir / "pass1_summary.csv"
    summary.to_csv(summary_out, index=False)
    log.info("  Saved: %s", summary_out)

    log.info("Pass 1 complete.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Spearman RBH Pass 1: cross-species cluster bridge."
    )
    parser.add_argument(
        "--mouse-pb",
        type=Path,
        default=_DEFAULT_MOUSE_PB,
        help="Mouse pseudobulk CSV (clusters x genes).",
    )
    parser.add_argument(
        "--human-pb",
        type=Path,
        default=_DEFAULT_HUMAN_PB,
        help="Human pseudobulk CSV (clusters x genes).",
    )
    parser.add_argument(
        "--ortholog-map",
        type=Path,
        default=_DEFAULT_ORTHOLOG,
        help="Ortholog mapping CSV with columns mouse_symbol, human_symbol, is_one_to_one.",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=3000,
        help="Number of highly variable genes to use (default: 3000).",
    )
    parser.add_argument(
        "--n-perm",
        type=int,
        default=1000,
        help="Permutations for threshold estimation (default: 1000).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help="Output directory (default: results/cluster_bridge/).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_pass1(
        mouse_pb_path=args.mouse_pb,
        human_pb_path=args.human_pb,
        ortholog_path=args.ortholog_map,
        n_hvgs=args.n_hvgs,
        n_perm=args.n_perm,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
