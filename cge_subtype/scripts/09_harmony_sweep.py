#!/usr/bin/env python3
"""
Step 09: Hyperparameter sweep over Harmony mapping configurations.

Runs 08_harmony_mapping.py in parallel for every combination of:
    n_pcs      in [30, 50, 75]
    n_hvgs     in [2000, 3000, 5000]
    theta      in [1.0, 2.0]
    n_neighbors in [15, 30, 50]

That gives 54 Harmony combos per dataset.  With --include-neuronal-only the
same 54 combos are also run with --neuronal-only, giving 108 total per dataset.

After all runs complete, *_metrics.csv files are collected into a single
sweep_summary.csv, sorted by overall_accuracy (descending).

Usage
-----
    # Sweep both datasets, default 10 workers
    python 09_harmony_sweep.py --dataset both

    # Sweep M1 only, 20 workers, include neuronal-only variants
    python 09_harmony_sweep.py --dataset M1 --n-processes 20 --include-neuronal-only
"""

from __future__ import annotations

import argparse
import itertools
import logging
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_08 = Path(__file__).resolve().parent / "08_harmony_mapping.py"

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
# Hyperparameter grid
# ---------------------------------------------------------------------------
GRID: dict[str, list] = {
    "n_pcs": [30, 50, 75],
    "n_hvgs": [2000, 3000, 5000],
    "theta": [1.0, 2.0],
    "n_neighbors": [15, 30, 50],
}


# ---------------------------------------------------------------------------
# Build full list of configs
# ---------------------------------------------------------------------------
def build_configs(
    datasets: list[str],
    include_neuronal_only: bool,
    outdir: Path,
) -> list[dict]:
    """Return a list of config dicts, one per (dataset, hyperparams, neuronal_only) combo."""
    configs: list[dict] = []

    keys = list(GRID.keys())
    values = list(GRID.values())

    for dataset in datasets:
        for combo in itertools.product(*values):
            params = dict(zip(keys, combo))
            for neuronal_only in ([False, True] if include_neuronal_only else [False]):
                configs.append(
                    {
                        "dataset": dataset,
                        "neuronal_only": neuronal_only,
                        "outdir": outdir,
                        **params,
                    }
                )

    return configs


# ---------------------------------------------------------------------------
# Worker function: runs 08_harmony_mapping.py as a subprocess
# ---------------------------------------------------------------------------
def _run_single(cfg: dict) -> dict:
    """Run one Harmony mapping configuration.

    Returns a dict with keys ``tag``, ``returncode``, ``stdout``, ``stderr``.
    """
    tag_parts = [
        cfg["dataset"],
        f"pcs{cfg['n_pcs']}",
        f"hvg{cfg['n_hvgs']}",
        f"theta{cfg['theta']}",
        f"k{cfg['n_neighbors']}",
    ]
    if cfg["neuronal_only"]:
        tag_parts.append("neuronly")
    tag = "_".join(tag_parts)

    cmd = [
        sys.executable,
        str(SCRIPT_08),
        "--dataset", cfg["dataset"],
        "--n-pcs", str(cfg["n_pcs"]),
        "--n-hvgs", str(cfg["n_hvgs"]),
        "--theta", str(cfg["theta"]),
        "--n-neighbors", str(cfg["n_neighbors"]),
        "--outdir", str(cfg["outdir"]),
    ]
    if cfg["neuronal_only"]:
        cmd.append("--neuronal-only")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    status = "OK" if result.returncode == 0 else "FAILED"
    log.info("[%s] %s  (rc=%d)", status, tag, result.returncode)
    if result.returncode != 0:
        log.warning("  stderr: %s", result.stderr[-500:] if result.stderr else "(none)")

    return {
        "tag": tag,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


# ---------------------------------------------------------------------------
# Collect metrics
# ---------------------------------------------------------------------------
def collect_metrics(outdir: Path) -> pd.DataFrame:
    """Glob all *_metrics.csv files in outdir and concatenate them."""
    metric_files = sorted(outdir.glob("*_metrics.csv"))
    if not metric_files:
        log.warning("No *_metrics.csv files found in %s.", outdir)
        return pd.DataFrame()

    dfs = []
    for fp in metric_files:
        try:
            dfs.append(pd.read_csv(fp))
        except Exception as exc:
            log.warning("Could not read %s: %s", fp, exc)

    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values("overall_accuracy", ascending=False).reset_index(
        drop=True
    )
    return combined


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    if args.dataset == "both":
        datasets = ["M1", "V1"]
    else:
        datasets = [args.dataset]

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    configs = build_configs(
        datasets=datasets,
        include_neuronal_only=args.include_neuronal_only,
        outdir=outdir,
    )

    n_total = len(configs)
    log.info(
        "Sweep: %d configurations | datasets=%s | neuronal_only=%s | workers=%d",
        n_total,
        datasets,
        args.include_neuronal_only,
        args.n_processes,
    )
    log.info("Output directory: %s", outdir)

    # ------------------------------------------------------------------
    # Run in parallel
    # ------------------------------------------------------------------
    if args.n_processes > 1:
        with Pool(processes=args.n_processes) as pool:
            results = pool.map(_run_single, configs)
    else:
        results = [_run_single(cfg) for cfg in configs]

    n_ok = sum(1 for r in results if r["returncode"] == 0)
    n_failed = n_total - n_ok
    log.info("Sweep complete: %d/%d succeeded, %d failed.", n_ok, n_total, n_failed)

    if n_failed > 0:
        failed_tags = [r["tag"] for r in results if r["returncode"] != 0]
        log.warning("Failed configurations: %s", failed_tags)

    # ------------------------------------------------------------------
    # Collect and summarise metrics
    # ------------------------------------------------------------------
    log.info("Collecting metrics from %s …", outdir)
    summary = collect_metrics(outdir)

    if summary.empty:
        log.warning("No metrics collected — sweep_summary.csv not written.")
        return

    summary_path = outdir / "sweep_summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("Saved sweep summary → %s  (%d rows)", summary_path, len(summary))

    # Log best config per dataset
    for dataset in datasets:
        ds_rows = summary[summary["dataset"] == dataset] if "dataset" in summary.columns else summary
        if ds_rows.empty:
            continue
        best = ds_rows.iloc[0]
        log.info(
            "Best config for %s: n_pcs=%s  n_hvgs=%s  theta=%s  k=%s  neuronal_only=%s"
            "  → accuracy=%.4f",
            dataset,
            best.get("n_pcs", "?"),
            best.get("n_hvgs", "?"),
            best.get("theta", "?"),
            best.get("n_neighbors", "?"),
            best.get("neuronal_only", "?"),
            best.get("overall_accuracy", float("nan")),
        )

    log.info("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    default_outdir = PROJECT_ROOT / "cge_subtype" / "results" / "harmony"

    parser = argparse.ArgumentParser(
        description=(
            "Hyperparameter sweep over Harmony mapping configurations "
            "(runs 08_harmony_mapping.py in parallel)."
        )
    )
    parser.add_argument(
        "--dataset",
        choices=["M1", "V1", "both"],
        default="both",
        help="Dataset(s) to sweep (default: both).",
    )
    parser.add_argument(
        "--n-processes",
        type=int,
        default=10,
        help="Number of parallel worker processes (default: 10).",
    )
    parser.add_argument(
        "--include-neuronal-only",
        action="store_true",
        default=False,
        help=(
            "If set, also run each config with --neuronal-only, doubling the total "
            "number of configurations."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=default_outdir,
        help=f"Output directory for all run results (default: {default_outdir}).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
