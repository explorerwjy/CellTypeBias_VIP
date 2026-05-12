"""
Run Allen IPFX on the F-I protocol (Intrinsic_IV curve_VIP_Erica) and save
ALL per-spike features. Used for downstream ISI / adaptation / train-level
feature extraction.

Output: dat/VIP_Ephys/allen_ipfx_vip/features_fi_perspike/<cell_id>_fi_spikes.csv
        Columns: cell_id, genotype, sweep_idx, stim_pA, stim_start_s, stim_end_s,
                 spike_idx_in_sweep, threshold_t, threshold_v, peak_t, peak_v,
                 trough_t, trough_v, fast_trough_t, fast_trough_v,
                 upstroke, downstroke, width, upstroke_downstroke_ratio
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
OUT_DIR = PROJ / "dat" / "VIP_Ephys" / "allen_ipfx_vip" / "features_fi_perspike"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PROTOCOL = "Intrinsic_IV curve_VIP_Erica"


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
    e = min(e, len(t) - 1)
    return s, e, float(np.median(i[s:e] - baseline))


def process_cell(cell_dir):
    cell_id = cell_dir.name
    if cell_id.startswith("WT"):
        geno = "WT"
    elif cell_id.startswith("Df16"):
        geno = "Df16"
    else:
        return None

    abf = None
    for f in sorted(cell_dir.glob("*.abf")):
        try:
            a = pyabf.ABF(str(f))
        except Exception:
            continue
        if a.protocol == PROTOCOL:
            abf = a
            break
    if abf is None:
        return None

    rows = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0)
        v = abf.sweepY.astype(float)
        t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1)
        i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None:
            continue
        s, e, amp = win
        ext = SpikeFeatureExtractor(start=t[s] - 0.005, end=min(t[-1], t[e] + 0.005),
                                    filter=4.0)
        try:
            spikes = ext.process(t, v, i)
            spikes = spikes[(spikes["threshold_t"] >= t[s])
                            & (spikes["threshold_t"] <= t[e])].reset_index(drop=True)
        except Exception:
            spikes = pd.DataFrame()
        for ix, row in spikes.iterrows():
            d = row.to_dict()
            d["cell_id"] = cell_id
            d["genotype"] = geno
            d["sweep_idx"] = sw
            d["stim_pA"] = round(amp, 1)
            d["stim_start_s"] = round(t[s], 4)
            d["stim_end_s"] = round(t[e], 4)
            d["spike_idx_in_sweep"] = ix
            rows.append(d)
    if not rows:
        return None
    df = pd.DataFrame(rows)
    out_cols = ["cell_id", "genotype", "sweep_idx", "stim_pA",
                "stim_start_s", "stim_end_s", "spike_idx_in_sweep"]
    out_cols += [c for c in df.columns if c not in out_cols]
    df = df[out_cols]
    df.to_csv(OUT_DIR / f"{cell_id}_fi_spikes.csv", index=False)
    return df


def main():
    n_done = 0
    for cd in sorted(REC_DIR.iterdir()):
        if not cd.is_dir():
            continue
        df = process_cell(cd)
        if df is None:
            continue
        n_sw = df["sweep_idx"].nunique()
        n_sp = len(df)
        print(f"  {cd.name}: {n_sp} spikes across {n_sw} sweeps")
        n_done += 1
    print(f"\nProcessed {n_done} cells.")


if __name__ == "__main__":
    main()
