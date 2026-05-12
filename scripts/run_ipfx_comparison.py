"""Re-extract VIP interneuron ephys features using Allen Institute's IPFX pipeline.

Validates whether the rise slope / AP width / peak voltage differences between
WT and Df16 VIP interneurons are robust across feature extraction pipelines.
"""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyabf
from scipy.stats import ttest_ind
from ipfx.feature_extractor import SpikeFeatureExtractor

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
REC_DIR = PROJ / "dat" / "VIP_Ephys" / "IV_step_recordings"
OUT_CSV = PROJ / "dat" / "VIP_Ephys" / "ipfx_features.csv"


def load_iv_protocol_abf(cell_dir):
    for f in sorted(cell_dir.glob("*.abf")):
        a = pyabf.ABF(str(f))
        if a.protocol == "Intrinsic_IV curve_VIP":
            return a
    return None


def extract_rheobase_features(cell_dir, filter_hz=4.0, stim_start=0.5, stim_end=2.0):
    abf = load_iv_protocol_abf(cell_dir)
    if abf is None:
        return None

    ext = SpikeFeatureExtractor(start=stim_start, end=stim_end, filter=filter_hz)

    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0)
        t = abf.sweepX.astype(float)
        v = abf.sweepY.astype(float)
        abf.setSweep(sw, channel=1)
        i = abf.sweepY.astype(float)

        try:
            spikes = ext.process(t, v, i)
        except Exception:
            continue

        if len(spikes) > 0:
            first = spikes.iloc[0]
            stim_mask = (t >= stim_start) & (t <= stim_end)
            pre_mask = (t >= 0.1) & (t < stim_start - 0.05)
            return {
                "rheobase_sweep": sw,
                "rheobase_pA": float(i[stim_mask].mean()),
                "RMP_mV": float(v[pre_mask].mean()),
                "n_spikes": len(spikes),
                "threshold_v": float(first["threshold_v"]),
                "peak_v": float(first["peak_v"]),
                "upstroke": float(first["upstroke"]),
                "downstroke": float(first["downstroke"]),
                "width_ms": float(first["width"]) * 1000,
                "ud_ratio": float(first["upstroke_downstroke_ratio"]),
                "trough_v": float(first["trough_v"]),
                "fast_trough_v": float(first["fast_trough_v"]),
            }
    return None


def main():
    rows = []
    cell_dirs = sorted(REC_DIR.iterdir())
    print(f"Processing {len(cell_dirs)} directories...")

    for cd in cell_dirs:
        if not cd.is_dir():
            continue
        cid = cd.name
        if cid.startswith("WT"):
            g = "WT"
        elif cid.startswith("Df16"):
            g = "Df16"
        else:
            continue

        try:
            feats = extract_rheobase_features(cd)
        except Exception as e:
            print(f"  {cid}: ERROR {e}")
            continue
        if feats is None:
            print(f"  {cid}: no spike detected")
            continue

        feats["cell_id"] = cid
        feats["genotype"] = g
        rows.append(feats)
        print(
            f"  {cid}: rheo={feats['rheobase_pA']:5.0f} pA, "
            f"upstroke={feats['upstroke']:6.1f} V/s, "
            f"width={feats['width_ms']:5.2f} ms, "
            f"peak={feats['peak_v']:6.1f} mV"
        )

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")
    print(f"Cells: {len(df)} ({(df.genotype=='WT').sum()} WT, {(df.genotype=='Df16').sum()} Df16)")

    print("\n" + "=" * 100)
    print("IPFX FEATURE COMPARISON: WT vs Df16")
    print("=" * 100)

    features_to_test = [
        ("upstroke", "Rise slope / peak dV/dt (V/s)"),
        ("peak_v", "Peak voltage (mV)"),
        ("width_ms", "AP width (ms)"),
        ("RMP_mV", "RMP (mV)"),
        ("downstroke", "Decay slope (V/s)"),
        ("threshold_v", "Threshold (mV)"),
        ("ud_ratio", "Up/down ratio"),
        ("rheobase_pA", "Rheobase (pA)"),
    ]

    print(
        f"\n{'Feature':35s} {'WT (mean±SD)':>22s} {'Df16 (mean±SD)':>22s} "
        f"{'t':>6s} {'p':>8s}  sig"
    )
    print("-" * 105)
    for col, label in features_to_test:
        wt = df[df.genotype == "WT"][col].dropna()
        df16 = df[df.genotype == "Df16"][col].dropna()
        if len(wt) < 3 or len(df16) < 3:
            continue
        t, p = ttest_ind(wt, df16)
        sig = (
            "***" if p < 0.001
            else "**" if p < 0.01
            else "*" if p < 0.05
            else "." if p < 0.1
            else ""
        )
        print(
            f"{label:35s} {wt.mean():9.2f}±{wt.std():7.2f} (n={len(wt):2d}) "
            f"{df16.mean():9.2f}±{df16.std():7.2f} (n={len(df16):2d}) "
            f"{t:6.2f} {p:8.4f}  {sig}"
        )


if __name__ == "__main__":
    main()
