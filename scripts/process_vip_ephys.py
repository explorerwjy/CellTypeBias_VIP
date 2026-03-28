#!/usr/bin/env python3
"""
Process VIP interneuron ABF recordings through the EphysSumStats pipeline.

For each cell directory in IV_step_recordings/:
  1. Find the Intrinsic_IV curve_VIP ABF file
  2. Convert to parquet bundles (mV + pA)
  3. Run sweep classification (sweep_config.json)
  4. Run full analysis (spike detection, input resistance, etc.)

Usage:
    python scripts/process_vip_ephys.py [--skip-excluded] [--dry-run]
"""
import sys
from pathlib import Path

# Add EphysSumStats to path
EPHYS_ROOT = Path("/home/jw3514/Work/NeurSim/EphysSumStats")
sys.path.insert(0, str(EPHYS_ROOT))

import json
import os
import re
import traceback
import pyabf
import pandas as pd
import numpy as np

# Override artifact thresholds for mouse VIP data before importing pipeline modules.
# At 20 kHz, a normal AP upstroke of 300 V/s creates a 15 mV per-sample jump,
# which exceeds the default 10 mV threshold designed for human NWB data.
import config.analysis_config as cfg
cfg.VOLTAGE_JUMP_THRESHOLD = 50.0   # mV — well above any real AP upstroke
cfg.SECOND_DERIV_THRESHOLD = 50e9   # increase to avoid false artifact flags

from src.sweep_classifier import process_bundle
from src.run_analysis import run_for_bundle

# ----- paths -----
DATA_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/VIP_Ephys/IV_step_recordings")
OUTPUT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/VIP_Ephys/analysis_bundles")


def find_iv_step_abf(cell_dir: Path) -> Path | None:
    """Find the Intrinsic_IV curve_VIP ABF file in a cell directory."""
    abfs = sorted(cell_dir.glob("*.abf"))
    for abf_path in abfs:
        abf = pyabf.ABF(str(abf_path))
        if abf.protocol == "Intrinsic_IV curve_VIP":
            return abf_path
    return None


def abf_to_bundle(abf_path: Path, cell_id: str, out_dir: Path):
    """
    Convert ABF file to parquet bundle compatible with EphysSumStats.

    Creates:
      out_dir/{cell_id}/
        mv_{cell_id}.parquet   (voltage traces)
        pa_{cell_id}.parquet   (current traces)
        manifest.json
    """
    abf = pyabf.ABF(str(abf_path))

    rows = []
    for sweep in abf.sweepList:
        for ch in range(abf.channelCount):
            abf.setSweep(sweepNumber=sweep, channel=ch)
            name = abf.adcNames[ch].strip() if hasattr(abf, "adcNames") else f"ADC{ch}"
            unit = abf.adcUnits[ch] if hasattr(abf, "adcUnits") else ""
            rows.append(pd.DataFrame({
                "sweep": sweep,
                "kind": "ADC",
                "channel_index": ch,
                "channel_name": name,
                "unit": unit,
                "t_s": abf.sweepX,
                "value": abf.sweepY,
            }))

    long_df = pd.concat(rows, ignore_index=True)
    df_mV = long_df[long_df["unit"] == "mV"].reset_index(drop=True)
    df_pA = long_df[long_df["unit"] == "pA"].reset_index(drop=True)

    # Subtract amplifier current offset so baseline reads ~0 pA.
    # The VIP protocol has ~1.5s pre-stimulus baseline (t < 1.5s).
    # The amplifier secondary output has a constant offset (~-1.5 pA)
    # that prevents the sweep classifier from detecting baseline periods.
    if not df_pA.empty:
        pre_stim = df_pA[df_pA["t_s"] < 1.0]
        if not pre_stim.empty:
            offset = pre_stim["value"].median()
            df_pA["value"] = df_pA["value"] - offset

    bundle_dir = out_dir / cell_id
    bundle_dir.mkdir(parents=True, exist_ok=True)

    mv_name = f"mv_{cell_id}.parquet"
    pa_name = f"pa_{cell_id}.parquet"
    df_mV.to_parquet(bundle_dir / mv_name, index=False)
    df_pA.to_parquet(bundle_dir / pa_name, index=False)

    # Parse genotype / mouse / cell from cell_id
    parts = cell_id.split("_")
    genotype = parts[0]  # WT or Df16

    manifest = {
        "file_id": cell_id,
        "abf_path": str(abf_path.resolve()),
        "tables": {"mv": mv_name, "pa": pa_name},
        "meta": {
            "genotype": genotype,
            "cell_id": cell_id,
            "protocol": abf.protocol,
            "sampleRate_Hz": int(abf.sampleRate),
            "sweepCount": int(abf.sweepCount),
            "sweepLength_sec": float(abf.sweepLengthSec),
            "source_file": abf_path.name,
        },
    }
    with open(bundle_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    return bundle_dir


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-excluded", action="store_true",
                        help="Skip cells marked as excluded or lostcell")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be processed without running")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cell_dirs = sorted(
        d for d in DATA_DIR.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )

    # Filter out excluded / lost cells
    if args.skip_excluded:
        cell_dirs = [d for d in cell_dirs if "excluded" not in d.name and "lostcell" not in d.name]
    else:
        cell_dirs = [d for d in cell_dirs if "lostcell" not in d.name]  # always skip empty dirs

    print(f"Found {len(cell_dirs)} cell directories to process")
    print(f"Output: {OUTPUT_DIR}\n")

    results = {"success": [], "failed": [], "skipped": []}

    for i, cell_dir in enumerate(cell_dirs, 1):
        cell_id = cell_dir.name
        print(f"\n{'='*70}")
        print(f"[{i}/{len(cell_dirs)}] {cell_id}")
        print(f"{'='*70}")

        if args.dry_run:
            abf_path = find_iv_step_abf(cell_dir)
            print(f"  Would process: {abf_path.name if abf_path else 'NO IV FILE FOUND'}")
            continue

        # Step 1: Find ABF file
        abf_path = find_iv_step_abf(cell_dir)
        if abf_path is None:
            print(f"  SKIP: No Intrinsic_IV curve_VIP file found")
            results["skipped"].append(cell_id)
            continue

        print(f"  ABF: {abf_path.name}")

        try:
            # Step 2: Create bundle
            print(f"  Creating bundle...")
            bundle_dir = abf_to_bundle(abf_path, cell_id, OUTPUT_DIR)
            print(f"  Bundle: {bundle_dir}")

            # Step 3: Sweep classification
            print(f"  Running sweep classification...")
            sweep_config = process_bundle(str(bundle_dir))
            n_valid = sweep_config.get("valid_sweeps", 0)
            n_total = len(sweep_config.get("sweeps", {}))
            print(f"  Sweeps: {n_valid}/{n_total} valid")

            # Step 4: Full analysis
            print(f"  Running full analysis...")
            run_for_bundle(str(bundle_dir))
            print(f"  DONE: {cell_id}")
            results["success"].append(cell_id)

        except Exception as e:
            print(f"  ERROR: {e}")
            traceback.print_exc()
            results["failed"].append((cell_id, str(e)))

    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    print(f"  Success: {len(results['success'])}")
    print(f"  Failed:  {len(results['failed'])}")
    print(f"  Skipped: {len(results['skipped'])}")
    if results["failed"]:
        print(f"\n  Failed cells:")
        for cell_id, err in results["failed"]:
            print(f"    {cell_id}: {err}")


if __name__ == "__main__":
    main()
