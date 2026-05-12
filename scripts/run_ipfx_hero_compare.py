"""Compare AP features using two alternative 'hero' sweep definitions:

  Option A (Allen convention): rheobase + 40 pA, average of first 5 spikes
  Option B: lowest-current sweep with >=5 spikes, average of first 5 spikes

Against the reference rheobase first-spike analysis (R1).
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
FILT = 4.0
N_AVG = 5  # average over first N spikes

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


def stim_amp(t, i):
    bm = (t >= 0.0) & (t < 0.4)
    sm = (t >= STIM_START) & (t <= STIM_END)
    return float(i[sm].mean() - i[bm].mean())


def spike_features_from_row(row):
    out = {}
    for feat in AP_FEATURES:
        out[feat] = float(row[feat]) * 1000 if feat == "width" else float(row[feat])
    out["ahp_mV"] = float(row["threshold_v"]) - float(row["trough_v"])
    out["amplitude_mV"] = float(row["peak_v"]) - float(row["threshold_v"])
    return out


def mean_spike_features(spikes_df, n=N_AVG):
    take = spikes_df.head(n)
    out = {}
    for feat in AP_FEATURES:
        vals = take[feat].astype(float)
        out[feat] = float(vals.mean()) * 1000 if feat == "width" else float(vals.mean())
    ahps = (take["threshold_v"] - take["trough_v"]).astype(float)
    amps = (take["peak_v"] - take["threshold_v"]).astype(float)
    out["ahp_mV"] = float(ahps.mean())
    out["amplitude_mV"] = float(amps.mean())
    out["n_spikes_avg"] = len(take)
    return out


def extract(cell_dir):
    abf = load_iv_protocol_abf(cell_dir)
    if abf is None:
        return None

    ext = SpikeFeatureExtractor(start=STIM_START, end=STIM_END, filter=FILT)
    sweeps = []
    for sw in range(abf.sweepCount):
        t, v, i = get_sweep(abf, sw)
        amp = stim_amp(t, i)
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

    rheo = spiking[0]
    out = {"rheobase_pA": rheo["amp"]}

    # --- R1: rheobase first spike ---
    rf = spike_features_from_row(rheo["spikes"].iloc[0])
    out.update({f"R1_{k}": v for k, v in rf.items()})

    # --- Hero A: rheobase + 40 pA, sweep closest to that current, average first 5 spikes ---
    target = rheo["amp"] + 40
    candidates_A = [s for s in spiking if s["n_spikes"] >= 2]
    if candidates_A:
        heroA = min(candidates_A, key=lambda s: abs(s["amp"] - target))
        out["A_stim_pA"] = heroA["amp"]
        out["A_n_spikes_total"] = heroA["n_spikes"]
        mf = mean_spike_features(heroA["spikes"])
        out.update({f"A_{k}": v for k, v in mf.items()})
    else:
        out["A_stim_pA"] = np.nan
        out["A_n_spikes_total"] = np.nan
        for k in AP_FEATURES + ["ahp_mV", "amplitude_mV", "n_spikes_avg"]:
            out[f"A_{k}"] = np.nan

    # --- Hero B: lowest-current sweep with >=5 spikes, average first 5 spikes ---
    candidates_B = [s for s in spiking if s["n_spikes"] >= N_AVG]
    if candidates_B:
        heroB = min(candidates_B, key=lambda s: s["amp"])
        out["B_stim_pA"] = heroB["amp"]
        out["B_n_spikes_total"] = heroB["n_spikes"]
        mf = mean_spike_features(heroB["spikes"])
        out.update({f"B_{k}": v for k, v in mf.items()})
    else:
        out["B_stim_pA"] = np.nan
        out["B_n_spikes_total"] = np.nan
        for k in AP_FEATURES + ["ahp_mV", "amplitude_mV", "n_spikes_avg"]:
            out[f"B_{k}"] = np.nan

    return out


def compare(df, label, prefix, feats):
    print(f"\n{'=' * 100}")
    print(f"{label}")
    print(f"{'=' * 100}")
    print(f"{'Feature':25s} {'WT':>20s} {'Df16':>20s} {'d':>7s} {'p':>9s}  sig")
    print("-" * 100)
    rows = []
    for f in feats:
        col = f"{prefix}_{f}"
        if col not in df.columns:
            continue
        wt = df[df.genotype == "WT"][col].dropna()
        df16 = df[df.genotype == "Df16"][col].dropna()
        if len(wt) < 3 or len(df16) < 3:
            continue
        t, p = ttest_ind(wt, df16)
        pooled = np.sqrt(((len(wt) - 1) * wt.std() ** 2 + (len(df16) - 1) * df16.std() ** 2)
                         / (len(wt) + len(df16) - 2))
        d = (wt.mean() - df16.mean()) / pooled if pooled > 0 else 0
        rows.append({"feat": f, "wt": wt, "df16": df16, "t": t, "p": p, "d": d})
    rows.sort(key=lambda r: r["p"])
    for r in rows:
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else "." if r["p"] < 0.1 else ""
        print(f"{r['feat']:25s} {r['wt'].mean():8.2f}±{r['wt'].std():6.2f}({len(r['wt']):2d}) "
              f"{r['df16'].mean():8.2f}±{r['df16'].std():6.2f}({len(r['df16']):2d}) "
              f"{r['d']:+7.2f} {r['p']:9.4f}  {sig}")
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
            f = extract(cd)
        except Exception as e:
            print(f"  {cid}: ERROR {e}")
            continue
        if f is None:
            continue
        f["cell_id"] = cid
        f["genotype"] = g
        rows.append(f)

    df = pd.DataFrame(rows)
    df.to_csv(PROJ / "dat/VIP_Ephys/ipfx_hero_compare.csv", index=False)

    # Report target current choices
    print(f"\nCells: {len(df)} ({(df.genotype=='WT').sum()} WT, {(df.genotype=='Df16').sum()} Df16)")
    print("\nHero A (rheobase+40 pA) usable:")
    a_ok = df["A_upstroke"].notna()
    print(f"  {a_ok.sum()}/{len(df)} cells")
    print(f"  WT target stim range: {df[a_ok & (df.genotype=='WT')]['A_stim_pA'].min():.1f}–{df[a_ok & (df.genotype=='WT')]['A_stim_pA'].max():.1f} pA, "
          f"mean n_spikes: {df[a_ok & (df.genotype=='WT')]['A_n_spikes_total'].mean():.1f}")
    print(f"  Df16 target stim range: {df[a_ok & (df.genotype=='Df16')]['A_stim_pA'].min():.1f}–{df[a_ok & (df.genotype=='Df16')]['A_stim_pA'].max():.1f} pA, "
          f"mean n_spikes: {df[a_ok & (df.genotype=='Df16')]['A_n_spikes_total'].mean():.1f}")

    print("\nHero B (lowest sweep with ≥5 spikes) usable:")
    b_ok = df["B_upstroke"].notna()
    print(f"  {b_ok.sum()}/{len(df)} cells")
    if b_ok.sum() > 0:
        print(f"  WT stim range: {df[b_ok & (df.genotype=='WT')]['B_stim_pA'].min():.1f}–{df[b_ok & (df.genotype=='WT')]['B_stim_pA'].max():.1f} pA, "
              f"mean n_spikes: {df[b_ok & (df.genotype=='WT')]['B_n_spikes_total'].mean():.1f}")
        print(f"  Df16 stim range: {df[b_ok & (df.genotype=='Df16')]['B_stim_pA'].min():.1f}–{df[b_ok & (df.genotype=='Df16')]['B_stim_pA'].max():.1f} pA, "
              f"mean n_spikes: {df[b_ok & (df.genotype=='Df16')]['B_n_spikes_total'].mean():.1f}")

    feats = ["upstroke", "downstroke", "width", "peak_v", "threshold_v",
             "amplitude_mV", "ahp_mV", "upstroke_downstroke_ratio", "trough_v", "fast_trough_v"]

    compare(df, "R1: Rheobase first spike (reference)", "R1", feats)
    compare(df, "Hero A: Rheobase + 40 pA, mean of first 5 spikes", "A", feats)
    compare(df, "Hero B: Lowest ≥5-spike sweep, mean of first 5 spikes", "B", feats)

    # Robustness summary
    print("\n" + "=" * 95)
    print("ROBUSTNESS: p-values across definitions")
    print("=" * 95)
    print(f"{'Feature':25s} {'R1 (rheo 1 spike)':>20s} {'A (rheo+40, avg5)':>22s} {'B (≥5sp sweep, avg5)':>24s}")
    print("-" * 95)
    for f in feats:
        ps = []
        for prefix in ["R1", "A", "B"]:
            col = f"{prefix}_{f}"
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
        print(f"{f:25s} {fmt(ps[0]):>20s} {fmt(ps[1]):>22s} {fmt(ps[2]):>24s}")


if __name__ == "__main__":
    main()
