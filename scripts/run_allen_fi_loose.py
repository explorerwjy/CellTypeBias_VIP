"""
Re-extract F-I curve using a permissive Allen IPFX configuration:
  - dv_cutoff = 5 V/s (default 20)        — allows attenuated upstrokes
  - min_peak = -25 mV (default -30)        — keeps the basic peak floor
  - min_height = 2 mV (default 2)          — unchanged
  - thresh_frac = 0.05 (default 0.05)      — unchanged

Tests the hypothesis that Erica's higher F-I rates come from counting
attenuated spikes that fail Allen's default dV/dt cutoff.
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
OUT_DIR = PROJ / "dat" / "VIP_Ephys" / "allen_ipfx_vip" / "features_fi_loose"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PROTOCOL = "Intrinsic_IV curve_VIP_Erica"


def detect_stim_window(t, i, threshold_pA=2.0):
    bm = np.percentile(i, 25); above = np.abs(i - bm) > threshold_pA
    if not above.any(): return None
    edges = np.diff(above.astype(int))
    starts = np.where(edges == 1)[0] + 1; ends = np.where(edges == -1)[0] + 1
    if above[0]: starts = np.concatenate([[0], starts])
    if above[-1]: ends = np.concatenate([ends, [len(above)]])
    runs = list(zip(starts, ends))
    if not runs: return None
    s, e = max(runs, key=lambda r: r[1] - r[0]); e = min(e, len(t) - 1)
    return s, e, float(np.median(i[s:e] - bm))


def process_cell(cd):
    cell_id = cd.name
    if cell_id.startswith("WT"): geno = "WT"
    elif cell_id.startswith("Df16"): geno = "Df16"
    else: return None
    abf = None
    for f in sorted(cd.glob("*.abf")):
        try: a = pyabf.ABF(str(f))
        except: continue
        if a.protocol == PROTOCOL: abf = a; break
    if abf is None: return None

    rows = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0); v = abf.sweepY.astype(float); t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1); i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None:
            rows.append({"cell_id": cell_id, "genotype": geno, "sweep_idx": sw,
                         "stim_pA": np.nan, "stim_dur_s": np.nan, "n_spikes": 0, "freq_Hz": 0.0})
            continue
        s, e, amp = win
        ext = SpikeFeatureExtractor(start=t[s] - 0.005, end=min(t[-1], t[e] + 0.005),
                                    filter=4.0, dv_cutoff=5.0, min_peak=-25.0,
                                    min_height=2.0, thresh_frac=0.05)
        try:
            spk = ext.process(t, v, i)
            in_w = (spk["threshold_t"] >= t[s]) & (spk["threshold_t"] <= t[e])
            n = int(in_w.sum()) if len(spk) > 0 else 0
        except Exception:
            n = 0
        dur = t[e] - t[s]
        rows.append({"cell_id": cell_id, "genotype": geno, "sweep_idx": sw,
                     "stim_pA": round(amp, 1), "stim_dur_s": round(dur, 4),
                     "n_spikes": n, "freq_Hz": n / dur if dur > 0 else 0.0})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / f"{cell_id}_allen_fi_loose.csv", index=False)
    return df


def main():
    n = 0
    for cd in sorted(REC_DIR.iterdir()):
        if not cd.is_dir(): continue
        out = process_cell(cd)
        if out is None: continue
        max_n = out["n_spikes"].max()
        max_amp = out["stim_pA"].max()
        print(f"  {cd.name}: max stim={max_amp:.0f} pA, max n_spikes={max_n}")
        n += 1
    print(f"\nProcessed {n} cells.")


if __name__ == "__main__":
    main()
