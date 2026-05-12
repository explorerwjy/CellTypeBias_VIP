"""Rigorous hero sweep analysis with three alternative definitions:

  A1 (strict): rheobase + 40 pA, mean of first 5 spikes, EXCLUDE cells with <5 spikes
  A2 (walk-up): lowest sweep >= rheobase + 40 pA with >= 5 spikes, mean of first 5
  A3 (all spikes): rheobase + 40 pA, mean of ALL spikes in sweep (Allen style)

Compares WT vs Df16 in each case.
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
N_AVG = 5
TARGET_OFFSET = 40  # pA above rheobase

AP_FEATURES = ["threshold_v", "peak_v", "upstroke", "downstroke",
               "width", "upstroke_downstroke_ratio", "trough_v", "fast_trough_v"]


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


def mean_spike_features(spikes_df, n=None):
    """Average AP features across first n spikes (or all if n is None)."""
    take = spikes_df.head(n) if n is not None else spikes_df
    out = {}
    for feat in AP_FEATURES:
        vals = take[feat].astype(float)
        out[feat] = float(vals.mean()) * 1000 if feat == "width" else float(vals.mean())
    ahps = (take["threshold_v"] - take["trough_v"]).astype(float)
    amps = (take["peak_v"] - take["threshold_v"]).astype(float)
    out["ahp_mV"] = float(ahps.mean())
    out["amplitude_mV"] = float(amps.mean())
    out["n_spikes_used"] = len(take)
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

    target = rheo["amp"] + TARGET_OFFSET

    # --- A1: strict rheo+40, first 5 spikes, EXCLUDE if <5 available ---
    cands = [s for s in spiking if abs(s["amp"] - target) <= 10]  # within ±10 pA of target
    cands_5 = [s for s in cands if s["n_spikes"] >= N_AVG]
    if cands_5:
        # pick closest to target
        chosen = min(cands_5, key=lambda s: abs(s["amp"] - target))
        out["A1_stim_pA"] = chosen["amp"]
        out["A1_n_total"] = chosen["n_spikes"]
        mf = mean_spike_features(chosen["spikes"], n=N_AVG)
        out.update({f"A1_{k}": v for k, v in mf.items()})
    else:
        out["A1_stim_pA"] = np.nan
        out["A1_n_total"] = np.nan
        for k in AP_FEATURES + ["ahp_mV", "amplitude_mV", "n_spikes_used"]:
            out[f"A1_{k}"] = np.nan

    # --- A2: walk up from rheo+40 until find a sweep with >=5 spikes ---
    walk = [s for s in spiking if s["amp"] >= target - 2 and s["n_spikes"] >= N_AVG]
    if not walk:
        # Fall back to any sweep above target
        walk = [s for s in spiking if s["amp"] >= target - 2]
    if walk:
        chosen = min(walk, key=lambda s: s["amp"])
        out["A2_stim_pA"] = chosen["amp"]
        out["A2_n_total"] = chosen["n_spikes"]
        mf = mean_spike_features(chosen["spikes"], n=min(N_AVG, chosen["n_spikes"]))
        out.update({f"A2_{k}": v for k, v in mf.items()})
    else:
        out["A2_stim_pA"] = np.nan
        out["A2_n_total"] = np.nan
        for k in AP_FEATURES + ["ahp_mV", "amplitude_mV", "n_spikes_used"]:
            out[f"A2_{k}"] = np.nan

    # --- A3: sweep closest to rheo+40, use ALL spikes ---
    cands_any = [s for s in spiking if s["n_spikes"] >= 1]
    chosen = min(cands_any, key=lambda s: abs(s["amp"] - target))
    out["A3_stim_pA"] = chosen["amp"]
    out["A3_n_total"] = chosen["n_spikes"]
    mf = mean_spike_features(chosen["spikes"], n=None)
    out.update({f"A3_{k}": v for k, v in mf.items()})

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
    df.to_csv(PROJ / "dat/VIP_Ephys/ipfx_hero_fixed.csv", index=False)

    # Summaries
    print(f"\nTotal cells: {len(df)}")
    print("\nA1 (strict rheo+40 ±10 pA, ≥5 spikes required):")
    a1 = df["A1_upstroke"].notna()
    print(f"  {a1.sum()}/{len(df)} usable | WT {(a1 & (df.genotype=='WT')).sum()}/10, Df16 {(a1 & (df.genotype=='Df16')).sum()}/19")
    if a1.sum() > 0:
        print(f"  Stim range: WT {df[a1 & (df.genotype=='WT')]['A1_stim_pA'].agg(['min','max','mean']).to_dict()}")
        print(f"  Stim range: Df16 {df[a1 & (df.genotype=='Df16')]['A1_stim_pA'].agg(['min','max','mean']).to_dict()}")

    print("\nA2 (walk up from rheo+40 to first ≥5-spike sweep):")
    a2 = df["A2_upstroke"].notna()
    print(f"  {a2.sum()}/{len(df)} usable | WT {(a2 & (df.genotype=='WT')).sum()}/10, Df16 {(a2 & (df.genotype=='Df16')).sum()}/19")
    if a2.sum() > 0:
        wta2 = df[a2 & (df.genotype=='WT')]
        da2 = df[a2 & (df.genotype=='Df16')]
        print(f"  Stim range: WT {wta2['A2_stim_pA'].min():.1f}-{wta2['A2_stim_pA'].max():.1f} (mean {wta2['A2_stim_pA'].mean():.1f})")
        print(f"  Stim range: Df16 {da2['A2_stim_pA'].min():.1f}-{da2['A2_stim_pA'].max():.1f} (mean {da2['A2_stim_pA'].mean():.1f})")
        print(f"  n_spikes used: WT mean {wta2['A2_n_spikes_used'].mean():.1f}, Df16 mean {da2['A2_n_spikes_used'].mean():.1f}")

    print("\nA3 (rheo+40, use ALL spikes in sweep):")
    a3 = df["A3_upstroke"].notna()
    print(f"  {a3.sum()}/{len(df)} usable")
    if a3.sum() > 0:
        wt3 = df[a3 & (df.genotype=='WT')]
        d3 = df[a3 & (df.genotype=='Df16')]
        print(f"  n_spikes used: WT mean {wt3['A3_n_spikes_used'].mean():.1f}, Df16 mean {d3['A3_n_spikes_used'].mean():.1f}")

    feats = ["upstroke", "downstroke", "width", "peak_v", "threshold_v",
             "amplitude_mV", "ahp_mV", "upstroke_downstroke_ratio", "trough_v", "fast_trough_v"]

    compare(df, "A1: Strict rheo+40 pA, first 5 spikes (exclude <5)", "A1", feats)
    compare(df, "A2: Walk up from rheo+40 to ≥5 spikes, first 5", "A2", feats)
    compare(df, "A3: Rheo+40 pA, mean of ALL spikes in sweep", "A3", feats)

    # Robustness table
    print("\n" + "=" * 105)
    print("ROBUSTNESS across hero definitions (p-values)")
    print("=" * 105)
    print(f"{'Feature':25s} {'A1 (strict, first5)':>22s} {'A2 (walkup, first5)':>24s} {'A3 (rheo+40, all)':>22s}")
    print("-" * 105)
    for f in feats:
        ps = []
        for prefix in ["A1", "A2", "A3"]:
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
        print(f"{f:25s} {fmt(ps[0]):>22s} {fmt(ps[1]):>24s} {fmt(ps[2]):>22s}")


if __name__ == "__main__":
    main()
