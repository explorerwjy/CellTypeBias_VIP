"""Extract AP features at 3 reference points to test robustness:
  1. Rheobase first spike (single AP at threshold)
  2. Rheobase + 20 pA first spike (defined suprathreshold)
  3. Hero sweep mean of first 5 spikes (robust averaging, most standardized)

Then run the same WT vs Df16 comparison for each.
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

STIM_START = 0.5
STIM_END = 2.0
FILT = 4.0  # kHz

AP_FEATURES = [
    "threshold_v", "peak_v", "upstroke", "downstroke",
    "width", "upstroke_downstroke_ratio",
    "trough_v", "fast_trough_v",
]


def load_iv_protocol_abf(cell_dir):
    for f in sorted(cell_dir.glob("*.abf")):
        a = pyabf.ABF(str(f))
        if a.protocol == "Intrinsic_IV curve_VIP":
            return a
    return None


def get_sweep(abf, sw):
    abf.setSweep(sw, channel=0)
    t = abf.sweepX.astype(float)
    v = abf.sweepY.astype(float)
    abf.setSweep(sw, channel=1)
    i = abf.sweepY.astype(float)
    return t, v, i


def stim_amplitude(t, i):
    bm = (t >= 0.0) & (t < 0.4)
    sm = (t >= STIM_START) & (t <= STIM_END)
    return float(i[sm].mean() - i[bm].mean())


def spike_row_to_dict(spike_row, prefix=""):
    """Extract AP features from an IPFX spike row."""
    out = {}
    for feat in AP_FEATURES:
        key = f"{prefix}_{feat}" if prefix else feat
        if feat == "width":
            out[key] = float(spike_row[feat]) * 1000  # ms
        else:
            out[key] = float(spike_row[feat])
    # Derived
    ahp_key = f"{prefix}_ahp_mV" if prefix else "ahp_mV"
    amp_key = f"{prefix}_amplitude_mV" if prefix else "amplitude_mV"
    out[ahp_key] = float(spike_row["threshold_v"]) - float(spike_row["trough_v"])
    out[amp_key] = float(spike_row["peak_v"]) - float(spike_row["threshold_v"])
    return out


def mean_spike_features(spikes_df, n_spikes=5, prefix=""):
    """Average AP features across the first N spikes of a sweep."""
    take = spikes_df.head(n_spikes)
    out = {}
    for feat in AP_FEATURES:
        key = f"{prefix}_{feat}" if prefix else feat
        vals = take[feat].astype(float)
        out[key] = float(vals.mean()) * 1000 if feat == "width" else float(vals.mean())
    ahp_key = f"{prefix}_ahp_mV" if prefix else "ahp_mV"
    amp_key = f"{prefix}_amplitude_mV" if prefix else "amplitude_mV"
    ahps = (take["threshold_v"] - take["trough_v"]).astype(float)
    amps = (take["peak_v"] - take["threshold_v"]).astype(float)
    out[ahp_key] = float(ahps.mean())
    out[amp_key] = float(amps.mean())
    out[f"{prefix}_n_spikes_avg" if prefix else "n_spikes_avg"] = len(take)
    return out


def extract_all(cell_dir):
    abf = load_iv_protocol_abf(cell_dir)
    if abf is None:
        return None

    ext = SpikeFeatureExtractor(start=STIM_START, end=STIM_END, filter=FILT)

    sweeps = []
    for sw in range(abf.sweepCount):
        t, v, i = get_sweep(abf, sw)
        amp = stim_amplitude(t, i)
        try:
            spikes = ext.process(t, v, i)
        except Exception:
            continue
        sweeps.append({"sw": sw, "amp": amp, "n_spikes": len(spikes), "spikes": spikes})

    if not sweeps:
        return None

    spiking = [s for s in sweeps if s["n_spikes"] > 0]
    if not spiking:
        return None

    # --- Reference 1: Rheobase first spike ---
    rheo = spiking[0]
    out = {"rheobase_pA": rheo["amp"]}
    out.update(spike_row_to_dict(rheo["spikes"].iloc[0], prefix="R1"))

    # --- Reference 2: Rheobase + ~20 pA first spike ---
    # Find the sweep whose amp is closest to rheobase + 20 pA and has at least 1 spike
    target_amp = rheo["amp"] + 20
    # prefer sweeps with amp >= target
    candidates = [s for s in spiking if s["amp"] >= target_amp - 2]
    if candidates:
        chosen = min(candidates, key=lambda s: abs(s["amp"] - target_amp))
        out["R2_stim_pA"] = chosen["amp"]
        out.update(spike_row_to_dict(chosen["spikes"].iloc[0], prefix="R2"))
    else:
        out["R2_stim_pA"] = np.nan
        for f in AP_FEATURES + ["ahp_mV", "amplitude_mV"]:
            out[f"R2_{f}"] = np.nan

    # --- Reference 3: Hero sweep (max spikes), mean of first 5 spikes ---
    hero = max(spiking, key=lambda s: s["n_spikes"])
    out["R3_stim_pA"] = hero["amp"]
    out["R3_n_spikes_total"] = hero["n_spikes"]
    out.update(mean_spike_features(hero["spikes"], n_spikes=5, prefix="R3"))

    return out


def compare(df, label, feature_cols):
    print(f"\n{'=' * 100}")
    print(f"{label}")
    print(f"{'=' * 100}")
    print(f"{'Feature':28s} {'WT':>20s} {'Df16':>20s} {'d':>7s} {'p':>9s}  sig")
    print("-" * 100)

    rows = []
    for col in feature_cols:
        if col not in df.columns:
            continue
        wt = df[df.genotype == "WT"][col].dropna()
        df16 = df[df.genotype == "Df16"][col].dropna()
        if len(wt) < 3 or len(df16) < 3:
            continue
        t, p = ttest_ind(wt, df16)
        pooled = np.sqrt(((len(wt) - 1) * wt.std() ** 2 + (len(df16) - 1) * df16.std() ** 2) / (len(wt) + len(df16) - 2))
        d = (wt.mean() - df16.mean()) / pooled if pooled > 0 else 0
        rows.append((col, wt, df16, t, p, d))

    rows.sort(key=lambda r: r[4])  # by p
    for col, wt, df16, t, p, d in rows:
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "." if p < 0.1 else ""
        name = col.split("_", 1)[1] if col.startswith(("R1_", "R2_", "R3_")) else col
        print(f"{name:28s} {wt.mean():8.2f}±{wt.std():6.2f}({len(wt):2d}) "
              f"{df16.mean():8.2f}±{df16.std():6.2f}({len(df16):2d}) "
              f"{d:+7.2f} {p:9.4f}  {sig}")
    return rows


def main():
    rows = []
    for cd in sorted(REC_DIR.iterdir()):
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
            f = extract_all(cd)
        except Exception as e:
            print(f"  {cid}: ERROR {e}")
            continue
        if f is None:
            continue
        f["cell_id"] = cid
        f["genotype"] = g
        rows.append(f)
        r1u = f.get("R1_upstroke", np.nan)
        r2u = f.get("R2_upstroke", np.nan)
        r3u = f.get("R3_upstroke", np.nan)
        r3n = f.get("R3_n_spikes_total", 0)
        print(f"  {cid:30s}: rheo={f['rheobase_pA']:5.1f} pA  "
              f"R1={r1u:6.1f}  R2={r2u:6.1f} (@{f.get('R2_stim_pA', np.nan):5.1f} pA)  "
              f"R3={r3u:6.1f} (hero={r3n}sp @{f.get('R3_stim_pA', np.nan):5.1f})")

    df = pd.DataFrame(rows)
    df.to_csv(PROJ / "dat/VIP_Ephys/ipfx_3way.csv", index=False)
    print(f"\nSaved. N = {len(df)}")

    # Feature lists
    ap_cols = ["upstroke", "downstroke", "width", "peak_v", "threshold_v",
               "upstroke_downstroke_ratio", "trough_v", "fast_trough_v",
               "ahp_mV", "amplitude_mV"]

    r1_cols = [f"R1_{c}" for c in ap_cols]
    r2_cols = [f"R2_{c}" for c in ap_cols]
    r3_cols = [f"R3_{c}" for c in ap_cols]

    rheobase_pa = df.groupby("genotype")["rheobase_pA"].agg(["mean", "std", "count"])
    print(f"\nRheobase: WT {rheobase_pa.loc['WT','mean']:.1f}±{rheobase_pa.loc['WT','std']:.1f}, "
          f"Df16 {rheobase_pa.loc['Df16','mean']:.1f}±{rheobase_pa.loc['Df16','std']:.1f}")

    r2_valid = df["R2_upstroke"].notna()
    print(f"R2 (rheobase+20 pA) available for {r2_valid.sum()}/{len(df)} cells")

    compare(df, "REFERENCE 1: Rheobase first spike", r1_cols)
    compare(df, "REFERENCE 2: Rheobase + 20 pA first spike", r2_cols)
    compare(df, "REFERENCE 3: Hero sweep, mean of first 5 spikes", r3_cols)

    # Summary: which features are robust across all 3 references?
    print("\n" + "=" * 100)
    print("ROBUSTNESS SUMMARY: features significant (p < 0.05) across reference points")
    print("=" * 100)
    print(f"{'Feature':28s} {'R1 (rheo)':>12s} {'R2 (rheo+20)':>14s} {'R3 (hero avg5)':>16s}")
    print("-" * 75)
    for base_col in ap_cols:
        ps = []
        for prefix in ["R1", "R2", "R3"]:
            col = f"{prefix}_{base_col}"
            if col not in df.columns:
                ps.append(np.nan); continue
            wt = df[df.genotype == "WT"][col].dropna()
            df16 = df[df.genotype == "Df16"][col].dropna()
            if len(wt) < 3 or len(df16) < 3:
                ps.append(np.nan); continue
            _, p = ttest_ind(wt, df16)
            ps.append(p)

        def fmt(p):
            if np.isnan(p):
                return "n/a"
            s = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "." if p < 0.1 else ""
            return f"{p:.4f} {s}"

        print(f"{base_col:28s} {fmt(ps[0]):>12s} {fmt(ps[1]):>14s} {fmt(ps[2]):>16s}")


if __name__ == "__main__":
    main()
