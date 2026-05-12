"""
Extract passive + intrinsic features (non-AP-shape) per cell using Allen IPFX:

  Subthreshold I-V protocol (Intrinsic_IV curve_VIP):
    - RMP (mV)         baseline voltage at 0 pA holding
    - IR (MΩ)          slope of V-I in subthreshold linear range
    - tau (ms)         single-exp fit of V decay during most-negative step
    - sag (mV)         min V - steady-state V during most-negative step
    - sag_ratio        sag / (RMP - V_steady)

  F-I protocol (Intrinsic_IV curve_VIP_Erica):
    - rheobase (pA)    smallest stim with at least one spike (Allen detection)
    - threshold (mV)   threshold_v of rheobase first AP
    - latency (ms)     threshold_t - stim_start at rheobase first AP
    - ahp (mV)         RMP - trough_v of rheobase first AP

Output: dat/VIP_Ephys/allen_ipfx_vip/passive_intrinsic.csv
"""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyabf
from ipfx.feature_extractor import SpikeFeatureExtractor

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
REC_DIR = PROJ / "dat" / "VIP_Ephys" / "IV_step_recordings"
OUT_CSV = PROJ / "dat" / "VIP_Ephys" / "allen_ipfx_vip" / "passive_intrinsic.csv"


def detect_stim_window(t, i, threshold_pA=2.0):
    baseline = np.percentile(i, 25)
    above = np.abs(i - baseline) > threshold_pA
    if not above.any():
        return None
    edges = np.diff(above.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1
    if above[0]:
        starts = np.concatenate([[0], starts])
    if above[-1]:
        ends = np.concatenate([ends, [len(above)]])
    runs = list(zip(starts, ends))
    if not runs:
        return None
    s, e = max(runs, key=lambda r: r[1] - r[0])
    return s, min(e, len(t) - 1), float(np.median(i[s:e] - baseline))


def fit_exp_tau(t_arr, v_arr):
    """Single exponential to V(t) = V_inf + (V0 - V_inf) * exp(-t/tau)."""
    from scipy.optimize import curve_fit
    v_inf = v_arr[-int(len(v_arr) * 0.1):].mean()
    v0 = v_arr[0]
    if abs(v_inf - v0) < 0.5:
        return np.nan
    try:
        def f(t, tau, v0_, vinf_):
            return vinf_ + (v0_ - vinf_) * np.exp(-t / tau)
        p0 = [0.05, v0, v_inf]
        popt, _ = curve_fit(f, t_arr - t_arr[0], v_arr, p0=p0, maxfev=2000,
                            bounds=([0.001, -120, -120], [1.0, 0, 0]))
        return float(popt[0]) * 1000  # ms
    except Exception:
        return np.nan


def subthreshold_features(abf):
    """RMP, IR, tau, sag from sub-threshold I-V protocol."""
    out = {"RMP_mV": np.nan, "IR_MOhm": np.nan, "tau_ms": np.nan,
           "sag_mV": np.nan, "sag_ratio": np.nan}
    # Per-sweep V_baseline, V_steady, V_min, stim_pA
    rows = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0)
        v = abf.sweepY.astype(float); t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1)
        i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None:
            continue
        s, e, amp = win
        # Use 100 ms before stim onset for baseline (auto-detected)
        pre_mask = (t >= max(0.0, t[s] - 0.4)) & (t < t[s] - 0.05)
        v_base = float(v[pre_mask].mean()) if pre_mask.any() else float(v[:int(0.05 * abf.dataRate)].mean())
        # last 30% of stim window for steady-state
        stim_dur = t[e] - t[s]
        steady_mask = (t >= t[s] + stim_dur * 0.7) & (t <= t[e])
        v_steady = float(v[steady_mask].mean()) if steady_mask.any() else np.nan
        v_min = float(v[s:e].min())
        # Spike detection over the actual stim window
        try:
            sfe = SpikeFeatureExtractor(start=t[s] - 0.005,
                                        end=min(t[-1], t[e] + 0.005),
                                        filter=4.0)
            spikes = sfe.process(t, v, i)
            in_win = (spikes["threshold_t"] >= t[s]) & (spikes["threshold_t"] <= t[e])
            n_sp = int(in_win.sum()) if len(spikes) > 0 else 0
        except Exception:
            n_sp = 0
        rows.append({"sweep": sw, "stim_pA": amp, "v_base": v_base,
                     "v_steady": v_steady, "v_min": v_min, "n_spikes": n_sp,
                     "v_trace": v[s:e].copy(), "t_trace": t[s:e].copy()})
    if not rows:
        return out
    sw_df = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith(("v_t", "t_t"))}
                          for r in rows])
    # RMP: average baseline across all sweeps
    out["RMP_mV"] = float(sw_df["v_base"].mean())
    # IR: hyperpolarizing-only sweeps (cleanest, no near-threshold rectification),
    # at least 3 sweeps, expand range until we have enough
    for cap in (8, 12, 20, 30):
        sub = sw_df[(sw_df["n_spikes"] == 0) & (sw_df["stim_pA"] <= 1) & (sw_df["stim_pA"] >= -cap)]
        if len(sub) >= 3:
            slope, _ = np.polyfit(sub["stim_pA"], sub["v_steady"], 1)
            out["IR_MOhm"] = float(slope * 1000)  # mV/pA → MΩ
            out["IR_n_sweeps"] = int(len(sub))
            break
    # tau: most-negative subthreshold step (must be hyperpolarizing)
    sw_df["delta"] = sw_df["v_steady"] - sw_df["v_base"]
    hyper = sw_df[(sw_df["stim_pA"] < -3) & (sw_df["n_spikes"] == 0)]
    if len(hyper) >= 1:
        sw_idx = hyper.iloc[hyper["stim_pA"].abs().argmax()]["sweep"]
        row = next(r for r in rows if r["sweep"] == sw_idx)
        # Use first 200 ms after step onset for fit
        rel_t = row["t_trace"] - row["t_trace"][0]
        keep = rel_t <= 0.2
        if keep.sum() > 20:
            out["tau_ms"] = fit_exp_tau(rel_t[keep], row["v_trace"][keep])
        # sag = peak (most negative) V minus steady-state V
        out["sag_mV"] = row["v_min"] - row["v_steady"]
        denom = (out["RMP_mV"] - row["v_steady"])
        out["sag_ratio"] = out["sag_mV"] / denom if abs(denom) > 0.5 else np.nan
    return out


def fi_protocol_features(abf):
    """Rheobase, threshold, latency, AHP from F-I protocol."""
    out = {"rheobase_pA": np.nan, "threshold_mV": np.nan,
           "latency_ms": np.nan, "ahp_mV": np.nan, "RMP_FI_mV": np.nan}
    rows = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0)
        v = abf.sweepY.astype(float); t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1)
        i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None:
            continue
        s, e, amp = win
        v_base = float(v[:int(0.05 * abf.dataRate)].mean())
        ext = SpikeFeatureExtractor(start=t[s] - 0.005, end=t[e] + 0.005, filter=4.0)
        try:
            spk = ext.process(t, v, i)
            in_win = (spk["threshold_t"] >= t[s]) & (spk["threshold_t"] <= t[e])
            spk = spk[in_win].reset_index(drop=True)
        except Exception:
            spk = pd.DataFrame()
        rows.append({"sweep": sw, "stim_pA": amp, "v_base": v_base,
                     "stim_start": t[s], "n_spikes": len(spk),
                     "first_spike": spk.iloc[0] if len(spk) > 0 else None})
    df = pd.DataFrame([{"sweep": r["sweep"], "stim_pA": r["stim_pA"],
                        "v_base": r["v_base"], "stim_start": r["stim_start"],
                        "n_spikes": r["n_spikes"]} for r in rows])
    if df.empty:
        return out
    out["RMP_FI_mV"] = float(df["v_base"].mean())
    spiking = df[df["n_spikes"] > 0].sort_values("stim_pA")
    if spiking.empty:
        return out
    rheo_row = spiking.iloc[0]
    rheo_data = next(r for r in rows if r["sweep"] == rheo_row["sweep"])
    out["rheobase_pA"] = float(rheo_row["stim_pA"])
    fa = rheo_data["first_spike"]
    if fa is not None:
        out["threshold_mV"] = float(fa["threshold_v"])
        out["latency_ms"] = (float(fa["threshold_t"]) - float(rheo_row["stim_start"])) * 1000
        out["ahp_mV"] = float(rheo_row["v_base"]) - float(fa["trough_v"])  # RMP - trough
    return out


def process_cell(cell_dir):
    cell_id = cell_dir.name
    if cell_id.startswith("WT"):
        geno = "WT"
    elif cell_id.startswith("Df16"):
        geno = "Df16"
    else:
        return None

    rec = {"cell_id": cell_id, "genotype": geno}

    # Subthreshold protocol
    sub_abf = None; fi_abf = None
    for f in sorted(cell_dir.glob("*.abf")):
        try:
            a = pyabf.ABF(str(f))
        except Exception:
            continue
        if a.protocol == "Intrinsic_IV curve_VIP" and sub_abf is None:
            sub_abf = a
        elif a.protocol == "Intrinsic_IV curve_VIP_Erica" and fi_abf is None:
            fi_abf = a
    if sub_abf is not None:
        rec.update(subthreshold_features(sub_abf))
    else:
        rec.update({"RMP_mV": np.nan, "IR_MOhm": np.nan, "tau_ms": np.nan,
                    "sag_mV": np.nan, "sag_ratio": np.nan})
    if fi_abf is not None:
        rec.update(fi_protocol_features(fi_abf))
    else:
        rec.update({"rheobase_pA": np.nan, "threshold_mV": np.nan,
                    "latency_ms": np.nan, "ahp_mV": np.nan, "RMP_FI_mV": np.nan})

    return rec


def main():
    rows = []
    for cd in sorted(REC_DIR.iterdir()):
        if not cd.is_dir():
            continue
        rec = process_cell(cd)
        if rec is None:
            continue
        rows.append(rec)
        print(f"  {rec['cell_id']}: RMP={rec.get('RMP_mV', np.nan):.1f}, "
              f"IR={rec.get('IR_MOhm', np.nan):.0f}, tau={rec.get('tau_ms', np.nan):.1f}, "
              f"sag={rec.get('sag_mV', np.nan):.1f}, "
              f"rheo={rec.get('rheobase_pA', np.nan):.0f}, "
              f"thresh={rec.get('threshold_mV', np.nan):.1f}, "
              f"latency={rec.get('latency_ms', np.nan):.1f}, "
              f"AHP={rec.get('ahp_mV', np.nan):.1f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {len(df)} cells to {OUT_CSV}")


if __name__ == "__main__":
    main()
