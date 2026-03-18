"""Aggregate per-cell ephys features from processed DANDI datasets.

Usage
-----
    python 10_aggregate_ephys.py [--outdir results/ephys_convergence/]

Reads:
    /mnt/data0/DANDI/Processed/000008/  (mouse, DANDI 000008)
    /mnt/data0/DANDI/Processed/000636/  (human, DANDI 000636)

Writes:
    <outdir>/mouse_ephys_features.csv
    <outdir>/human_ephys_features.csv
    <outdir>/combined_ephys_features.csv

NOTE: Do not run this script interactively — it can take several hours to
walk all cell directories.  Use a cluster job or screen session.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure the project root is on the import path when running as a script
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from cge_subtype.src.ephys_harmonization import (
    aggregate_cell_features,
    LOG_FEATURES,
    ALL_FEATURES,
)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

MOUSE_PROCESSED_DIR = "/mnt/data0/DANDI/Processed/000008"
HUMAN_PROCESSED_DIR = "/mnt/data0/DANDI/Processed/000636"
DEFAULT_OUTDIR = "results/ephys_convergence"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def log_transform_skewed(df: pd.DataFrame) -> pd.DataFrame:
    """Apply log1p to LOG_FEATURES columns that are present in *df*.

    Negative values (e.g. resting Vm) are skipped with a warning.
    Returns a copy of *df* with transformed columns.
    """
    df = df.copy()
    for feat in LOG_FEATURES:
        if feat not in df.columns:
            continue
        col = df[feat]
        if (col.dropna() < 0).any():
            logger.warning(
                "Feature '%s' has negative values; skipping log-transform.", feat
            )
            continue
        df[feat] = np.log1p(col)
        logger.debug("Log-transformed '%s'.", feat)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Aggregate per-cell ephys features from processed DANDI datasets."
    )
    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=(
            f"Output directory for CSV files (default: {DEFAULT_OUTDIR}). "
            "Will be created if it does not exist."
        ),
    )
    parser.add_argument(
        "--mouse-dir",
        default=MOUSE_PROCESSED_DIR,
        help=f"Path to processed mouse DANDI directory (default: {MOUSE_PROCESSED_DIR}).",
    )
    parser.add_argument(
        "--human-dir",
        default=HUMAN_PROCESSED_DIR,
        help=f"Path to processed human DANDI directory (default: {HUMAN_PROCESSED_DIR}).",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("Output directory: %s", outdir.resolve())

    # ------------------------------------------------------------------
    # 1. Aggregate mouse features (DANDI 000008)
    # ------------------------------------------------------------------
    logger.info("Aggregating mouse ephys features from %s …", args.mouse_dir)
    mouse_df = aggregate_cell_features(args.mouse_dir, species_name="mouse")
    logger.info("Mouse: %d cells, %d features", len(mouse_df), mouse_df.shape[1] - 2)

    # ------------------------------------------------------------------
    # 2. Aggregate human features (DANDI 000636)
    # ------------------------------------------------------------------
    logger.info("Aggregating human ephys features from %s …", args.human_dir)
    human_df = aggregate_cell_features(args.human_dir, species_name="human")
    logger.info("Human: %d cells, %d features", len(human_df), human_df.shape[1] - 2)

    # ------------------------------------------------------------------
    # 3. Log-transform skewed features
    # ------------------------------------------------------------------
    logger.info("Log-transforming skewed features: %s", LOG_FEATURES)
    mouse_df = log_transform_skewed(mouse_df)
    human_df = log_transform_skewed(human_df)

    # ------------------------------------------------------------------
    # 4. Combine and save
    # ------------------------------------------------------------------
    combined_df = pd.concat([mouse_df, human_df], ignore_index=True)
    logger.info(
        "Combined: %d cells total (%d mouse + %d human)",
        len(combined_df),
        (combined_df["species"] == "mouse").sum(),
        (combined_df["species"] == "human").sum(),
    )

    mouse_out = outdir / "mouse_ephys_features.csv"
    human_out = outdir / "human_ephys_features.csv"
    combined_out = outdir / "combined_ephys_features.csv"

    mouse_df.to_csv(mouse_out, index=False)
    human_df.to_csv(human_out, index=False)
    combined_df.to_csv(combined_out, index=False)

    logger.info("Saved mouse features   → %s", mouse_out)
    logger.info("Saved human features   → %s", human_out)
    logger.info("Saved combined features → %s", combined_out)

    # ------------------------------------------------------------------
    # 5. Summary statistics
    # ------------------------------------------------------------------
    logger.info("Summary of combined data:")
    logger.info("  Columns: %s", list(combined_df.columns))
    logger.info("  NaN rates (top features):")
    feat_cols = [c for c in ALL_FEATURES if c in combined_df.columns]
    nan_rates = combined_df[feat_cols].isna().mean().sort_values(ascending=False)
    for feat, rate in nan_rates.head(10).items():
        if rate > 0:
            logger.info("    %s: %.1f%%", feat, rate * 100)

    return combined_df


if __name__ == "__main__":
    main()
