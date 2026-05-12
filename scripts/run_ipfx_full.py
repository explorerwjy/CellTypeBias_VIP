"""Extract comprehensive ephys features using IPFX, mirroring the Allen Cell Types pipeline.

Per-cell features:
- AP waveform (rheobase first spike): upstroke, downstroke, width, peak, threshold, trough
- Spike train (suprathreshold sweep): firing rate, ISI, CV, adaptation index, latency, F-I slope
- Subthreshold: input resistance, sag ratio, membrane tau, baseline (RMP)
"""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyabf
from scipy.stats import ttest_ind

from ipfx.feature_extractor import SpikeFeatureExtractor, SpikeTrainFeatureExtractor
import ipfx.spike_train_features as stf
import ipfx.subthresh_features as subf

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
REC_DIR = PROJ / "dat" / "VIP_Ephys" / "IV_step_recordings"
OUT_CSV = PROJ / "dat" / "VIP_Ephys" / "ipfx_features_full.csv"

STIM_START = 0.5
STIM_END = 2.0
FILT = 4.0  # kHz Bessel for 20 kHz sampling


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


def stim_amplitude(t, i, baseline_window=(0.0, 0.4)):
    """Mean injected current during stim minus baseline current."""
    bm = (t >= baseline_window[0]) & (t < baseline_window[1])
    sm = (t >= STIM_START) & (t <= STIM_END)
    return float(i[sm].mean() - i[bm].mean())


def extract_cell_features(cell_dir):
    abf = load_iv_protocol_abf(cell_dir)
    if abf is None:
        return None

    spike_ext = SpikeFeatureExtractor(start=STIM_START, end=STIM_END, filter=FILT)
    train_ext = SpikeTrainFeatureExtractor(start=STIM_START, end=STIM_END)

    out = {"n_sweeps": abf.sweepCount}

    # ---- Pass 1: catalog all sweeps ----
    sweep_info = []  # list of dicts
    for sw in range(abf.sweepCount):
        t, v, i = get_sweep(abf, sw)
        amp = stim_amplitude(t, i)
        try:
            spikes = spike_ext.process(t, v, i)
        except Exception:
            continue
        n = len(spikes)
        sweep_info.append({"sw": sw, "amp": amp, "n_spikes": n, "spikes": spikes, "t": t, "v": v, "i": i})

    if not sweep_info:
        return None

    swdf = pd.DataFrame([{"sw": s["sw"], "amp": s["amp"], "n_spikes": s["n_spikes"]} for s in sweep_info])

    # ---- RMP from baseline of sweep 0 ----
    s0 = sweep_info[0]
    pre_mask = (s0["t"] >= 0.1) & (s0["t"] < STIM_START - 0.05)
    out["RMP_mV"] = float(s0["v"][pre_mask].mean())

    # ---- Subthreshold features: input resistance from negative-current sweeps ----
    # Gather all sub-threshold sweeps with negative current injection
    neg_sweeps = [s for s in sweep_info if s["amp"] < -2 and s["n_spikes"] == 0]
    if len(neg_sweeps) >= 2:
        amps = []
        deflections = []
        for s in neg_sweeps:
            try:
                defl, _ = subf.voltage_deflection(s["t"], s["v"], s["i"], STIM_START, STIM_END, deflect_type="min")
                amps.append(s["amp"])
                deflections.append(defl - out["RMP_mV"])
            except Exception:
                pass
        if len(amps) >= 2:
            # Linear fit V = R * I → slope is Rin in MΩ if I is in pA and V in mV
            slope, _ = np.polyfit(amps, deflections, 1)
            out["Rin_MOhm"] = float(slope * 1000.0)  # mV/pA → MΩ

        # Sag from the most hyperpolarizing sweep
        try:
            most_neg = min(neg_sweeps, key=lambda s: s["amp"])
            sag, _ = subf.sag(most_neg["t"], most_neg["v"], most_neg["i"], STIM_START, STIM_END,
                              peak_width=0.005, baseline_interval=0.1)
            out["sag"] = float(sag)
        except Exception:
            out["sag"] = np.nan

        # Membrane tau
        try:
            tau = subf.time_constant(most_neg["t"], most_neg["v"], most_neg["i"], STIM_START, STIM_END,
                                      frac=0.1, baseline_interval=0.1)
            out["tau_s"] = float(tau)
        except Exception:
            out["tau_s"] = np.nan

    # ---- Rheobase: first sweep with a spike ----
    spiking_sweeps = [s for s in sweep_info if s["n_spikes"] > 0]
    if not spiking_sweeps:
        return out

    rheo = spiking_sweeps[0]
    out["rheobase_pA"] = rheo["amp"]

    first = rheo["spikes"].iloc[0]
    out["rheo_threshold_v"] = float(first["threshold_v"])
    out["rheo_peak_v"] = float(first["peak_v"])
    out["rheo_upstroke"] = float(first["upstroke"])
    out["rheo_downstroke"] = float(first["downstroke"])
    out["rheo_width_ms"] = float(first["width"]) * 1000
    out["rheo_ud_ratio"] = float(first["upstroke_downstroke_ratio"])
    out["rheo_trough_v"] = float(first["trough_v"])
    out["rheo_fast_trough_v"] = float(first["fast_trough_v"])
    out["rheo_amplitude_mV"] = out["rheo_peak_v"] - out["rheo_threshold_v"]
    # AHP relative to threshold
    out["rheo_ahp_mV"] = out["rheo_threshold_v"] - out["rheo_trough_v"]

    # ---- Hero sweep: highest firing rate ----
    hero = max(spiking_sweeps, key=lambda s: s["n_spikes"])
    out["hero_amp_pA"] = hero["amp"]
    out["hero_n_spikes"] = hero["n_spikes"]
    out["hero_freq_Hz"] = hero["n_spikes"] / (STIM_END - STIM_START)

    # Spike-train features for the hero sweep
    try:
        train_features = train_ext.process(hero["t"], hero["v"], hero["i"], hero["spikes"], extra_features=["adapt"])
    except Exception:
        train_features = {}

    if isinstance(train_features, dict):
        for k in ["adapt", "latency", "isi_cv", "mean_isi", "first_isi"]:
            if k in train_features:
                val = train_features[k]
                out[f"hero_{k}"] = float(val) if val is not None and not np.isnan(val) else np.nan

    # ---- F-I curve slope ----
    fi_amps = [s["amp"] for s in spiking_sweeps]
    fi_rates = [s["n_spikes"] / (STIM_END - STIM_START) for s in spiking_sweeps]
    if len(fi_amps) >= 3:
        try:
            fi_slope = stf.fit_fi_slope(np.array(fi_amps), np.array(fi_rates))
            out["fi_slope_HzpA"] = float(fi_slope)
        except Exception:
            out["fi_slope_HzpA"] = np.nan
    out["max_freq_Hz"] = max(fi_rates) if fi_rates else np.nan

    # ---- Latency to first spike (rheobase) ----
    out["rheo_latency_ms"] = (float(first["threshold_t"]) - STIM_START) * 1000

    return out


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
            f = extract_cell_features(cd)
        except Exception as e:
            print(f"  {cid}: ERROR {e}")
            continue
        if f is None:
            continue
        f["cell_id"] = cid
        f["genotype"] = g
        rows.append(f)
        print(f"  {cid}: rheo={f.get('rheobase_pA', 'NA')!s:>8}, hero_freq={f.get('hero_freq_Hz', 'NA')!s:>6}, Rin={f.get('Rin_MOhm', 'NA')!s:>8}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")
    print(f"Cells: {len(df)} ({(df.genotype=='WT').sum()} WT, {(df.genotype=='Df16').sum()} Df16)")

    # Test ALL numeric features
    print("\n" + "=" * 110)
    print("WT vs Df16: ALL features (sorted by p-value)")
    print("=" * 110)

    feat_cols = [c for c in df.columns if c not in ("cell_id", "genotype") and df[c].dtype != object]

    results = []
    for col in feat_cols:
        wt = df[df.genotype == "WT"][col].dropna()
        df16 = df[df.genotype == "Df16"][col].dropna()
        if len(wt) < 3 or len(df16) < 3:
            continue
        t, p = ttest_ind(wt, df16)
        # Cohen's d
        pooled = np.sqrt(((len(wt) - 1) * wt.std() ** 2 + (len(df16) - 1) * df16.std() ** 2) / (len(wt) + len(df16) - 2))
        d = (wt.mean() - df16.mean()) / pooled if pooled > 0 else 0
        results.append({
            "feature": col, "wt_mean": wt.mean(), "wt_sd": wt.std(), "n_wt": len(wt),
            "df16_mean": df16.mean(), "df16_sd": df16.std(), "n_df16": len(df16),
            "t": t, "p": p, "cohens_d": d,
        })

    res_df = pd.DataFrame(results).sort_values("p")
    print(f"\n{'Feature':25s} {'WT':>20s} {'Df16':>20s} {'d':>7s} {'t':>6s} {'p':>9s}")
    print("-" * 95)
    for _, r in res_df.iterrows():
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else "." if r["p"] < 0.1 else ""
        print(f"{r['feature']:25s} {r['wt_mean']:8.2f}±{r['wt_sd']:6.2f}({r['n_wt']:2d}) "
              f"{r['df16_mean']:8.2f}±{r['df16_sd']:6.2f}({r['n_df16']:2d}) "
              f"{r['cohens_d']:+7.2f} {r['t']:+6.2f} {r['p']:9.4f} {sig}")


if __name__ == "__main__":
    main()
