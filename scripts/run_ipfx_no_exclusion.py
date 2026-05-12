"""Inclusive hero-sweep analyses that keep all 29 cells:

  H1: rheobase + 40 pA ± 10 pA, first AP only
  H2: rheobase + 40 pA ± 10 pA, mean of all APs in sweep (2 or more)
  H3: rheobase + 40 pA ± 10 pA, mean of up to first 5 APs

All three keep all cells. H1 uses 1 spike per cell (most standardized),
H2 uses all available spikes per cell (maximum averaging),
H3 averages up to 5 spikes but takes fewer if fewer exist.
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
TARGET_OFFSET = 40  # pA above rheobase
WINDOW = 10  # ±10 pA from target

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


def first_spike_features(spikes_df):
    row = spikes_df.iloc[0]
    out = {}
    for f in AP_FEATURES:
        out[f] = float(row[f]) * 1000 if f == "width" else float(row[f])
    out["ahp_mV"] = float(row["threshold_v"]) - float(row["trough_v"])
    out["amplitude_mV"] = float(row["peak_v"]) - float(row["threshold_v"])
    out["n_spikes_used"] = 1
    return out


def mean_spike_features(spikes_df, n=None):
    take = spikes_df.head(n) if n is not None else spikes_df
    out = {}
    for f in AP_FEATURES:
        vals = take[f].astype(float)
        out[f] = float(vals.mean()) * 1000 if f == "width" else float(vals.mean())
    ahps = (take["threshold_v"] - take["trough_v"]).astype(float)
    amps = (take["peak_v"] - take["threshold_v"]).astype(float)
    out["ahp_mV"] = float(ahps.mean())
    out["amplitude_mV"] = float(amps.mean())
    out["n_spikes_used"] = len(take)
    return out


def pick_hero_sweep(spiking, rheobase_amp):
    """Pick the sweep closest to rheobase + TARGET_OFFSET within ±WINDOW."""
    target = rheobase_amp + TARGET_OFFSET
    cands = [s for s in spiking if abs(s["amp"] - target) <= WINDOW]
    if not cands:
        # Fall back: closest spiking sweep at or above rheobase
        cands = spiking
    return min(cands, key=lambda s: abs(s["amp"] - target))


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

    hero = pick_hero_sweep(spiking, rheo["amp"])
    out["hero_stim_pA"] = hero["amp"]
    out["hero_n_spikes"] = hero["n_spikes"]
    out["hero_offset_pA"] = hero["amp"] - rheo["amp"]

    # H1: first spike only
    h1 = first_spike_features(hero["spikes"])
    out.update({f"H1_{k}": v for k, v in h1.items()})

    # H2: mean of all spikes
    h2 = mean_spike_features(hero["spikes"], n=None)
    out.update({f"H2_{k}": v for k, v in h2.items()})

    # H3: mean of up to first 5 spikes
    h3 = mean_spike_features(hero["spikes"], n=5)
    out.update({f"H3_{k}": v for k, v in h3.items()})

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
    df.to_csv(PROJ / "dat/VIP_Ephys/ipfx_no_exclusion.csv", index=False)

    # Report the distribution of hero sweep stim levels and spike counts
    print(f"\nTotal cells: {len(df)} ({(df.genotype=='WT').sum()} WT, {(df.genotype=='Df16').sum()} Df16)")
    print(f"\nHero sweep: target = rheobase + {TARGET_OFFSET} pA, window ± {WINDOW} pA")
    print(f"\nPer-group summary:")
    for g in ["WT", "Df16"]:
        sub = df[df.genotype == g]
        print(f"  {g} (n={len(sub)}):")
        print(f"    rheobase: mean {sub['rheobase_pA'].mean():.1f} ± {sub['rheobase_pA'].std():.1f} pA")
        print(f"    hero stim: mean {sub['hero_stim_pA'].mean():.1f} ± {sub['hero_stim_pA'].std():.1f} pA")
        print(f"    hero offset from rheobase: mean {sub['hero_offset_pA'].mean():.1f} ± {sub['hero_offset_pA'].std():.1f} pA")
        print(f"    hero spike count: mean {sub['hero_n_spikes'].mean():.1f} ± {sub['hero_n_spikes'].std():.1f} (min {sub['hero_n_spikes'].min():.0f}, max {sub['hero_n_spikes'].max():.0f})")

    # Confirm no cells are dropped
    n_drop_h1 = df["H1_upstroke"].isna().sum()
    n_drop_h2 = df["H2_upstroke"].isna().sum()
    n_drop_h3 = df["H3_upstroke"].isna().sum()
    print(f"\nCells dropped: H1 {n_drop_h1}, H2 {n_drop_h2}, H3 {n_drop_h3} (should all be 0)")

    feats = ["upstroke", "downstroke", "width", "peak_v", "threshold_v",
             "amplitude_mV", "ahp_mV", "upstroke_downstroke_ratio", "trough_v", "fast_trough_v"]

    compare(df, f"H1: First AP at rheobase + {TARGET_OFFSET} pA (1 spike/cell, all 29 cells)", "H1", feats)
    compare(df, f"H2: Mean of ALL APs in sweep (variable n per cell, all 29 cells)", "H2", feats)
    compare(df, f"H3: Mean of up to first 5 APs (≤5 per cell, all 29 cells)", "H3", feats)

    # Robustness table
    print("\n" + "=" * 105)
    print("ROBUSTNESS across inclusive definitions (p-values)")
    print("=" * 105)
    print(f"{'Feature':25s} {'H1 (1st AP)':>15s} {'H2 (all APs)':>17s} {'H3 (≤5 APs)':>17s}")
    print("-" * 75)
    for f in feats:
        ps = []
        ds = []
        for prefix in ["H1", "H2", "H3"]:
            col = f"{prefix}_{f}"
            if col not in df.columns:
                ps.append(np.nan); ds.append(np.nan); continue
            wt = df[df.genotype == "WT"][col].dropna()
            df16 = df[df.genotype == "Df16"][col].dropna()
            if len(wt) < 3 or len(df16) < 3:
                ps.append(np.nan); ds.append(np.nan); continue
            _, p = ttest_ind(wt, df16)
            pooled = np.sqrt(((len(wt) - 1) * wt.std() ** 2 + (len(df16) - 1) * df16.std() ** 2)
                             / (len(wt) + len(df16) - 2))
            d = (wt.mean() - df16.mean()) / pooled if pooled > 0 else 0
            ps.append(p); ds.append(d)

        def fmt(p, d):
            if np.isnan(p):
                return "n/a"
            s = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "." if p < 0.1 else ""
            return f"p={p:.4f}{s}"
        print(f"{f:25s} {fmt(ps[0],ds[0]):>15s} {fmt(ps[1],ds[1]):>17s} {fmt(ps[2],ds[2]):>17s}")


if __name__ == "__main__":
    main()
