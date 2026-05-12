"""
Run Allen IPFX SpikeFeatureExtractor on the F-I protocol used to generate
manuscript Figure 4D (Intrinsic_IV curve_VIP_Erica).

Each sweep in this protocol delivers a single ~200 ms current pulse at one
amplitude. The stim window varies in time across cells, so we detect it from
the current channel of each sweep before invoking the spike extractor.

Output: dat/VIP_Ephys/allen_ipfx_vip/features_fi/<cell_id>_allen_fi.csv
        Columns: cell_id, genotype, sweep_idx, stim_pA, stim_start_s,
                 stim_end_s, stim_dur_s, n_spikes, freq_Hz
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
OUT_DIR = PROJ / "dat" / "VIP_Ephys" / "allen_ipfx_vip" / "features_fi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROTOCOL = "Intrinsic_IV curve_VIP_Erica"
FILTER_KHZ = 4.0


def detect_stim_window(t, i, threshold_pA=2.0, min_dur_s=0.05):
    """Find the stim window: largest contiguous region where |I - baseline| > threshold."""
    baseline = np.percentile(i, 25)
    above = np.abs(i - baseline) > threshold_pA
    if not above.any():
        return None
    # find runs of True
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
    # largest run
    start, end = max(runs, key=lambda r: r[1] - r[0])
    dur = (end - start) / len(t) * (t[-1] - t[0])
    if dur < min_dur_s:
        return None
    return float(t[start]), float(t[min(end, len(t) - 1)]), float(np.median(i[start:end] - baseline))


def process_cell(cell_dir):
    cell_id = cell_dir.name
    if cell_id.startswith("WT"):
        geno = "WT"
    elif cell_id.startswith("Df16"):
        geno = "Df16"
    else:
        return None

    abf_file = None
    for f in sorted(cell_dir.glob("*.abf")):
        try:
            a = pyabf.ABF(str(f))
        except Exception:
            continue
        if a.protocol == PROTOCOL:
            abf_file = a
            break
    if abf_file is None:
        return None

    rows = []
    for sw in range(abf_file.sweepCount):
        abf_file.setSweep(sw, channel=0)
        v = abf_file.sweepY.astype(float)
        t = abf_file.sweepX.astype(float)
        abf_file.setSweep(sw, channel=1)
        i = abf_file.sweepY.astype(float)

        win = detect_stim_window(t, i)
        if win is None:
            rows.append({"cell_id": cell_id, "genotype": geno, "sweep_idx": sw,
                         "stim_pA": np.nan, "stim_start_s": np.nan, "stim_end_s": np.nan,
                         "stim_dur_s": np.nan, "n_spikes": 0, "freq_Hz": 0.0})
            continue
        s0, s1, amp = win

        # Pad spike window slightly so threshold/peak detection has room
        ext = SpikeFeatureExtractor(start=max(0.0, s0 - 0.005),
                                    end=min(t[-1], s1 + 0.005),
                                    filter=FILTER_KHZ)
        try:
            spk = ext.process(t, v, i)
            # Only count spikes whose threshold falls inside the stim window
            if len(spk) > 0 and "threshold_t" in spk.columns:
                in_win = (spk["threshold_t"] >= s0) & (spk["threshold_t"] <= s1)
                n = int(in_win.sum())
            else:
                n = int(len(spk))
        except Exception:
            n = 0

        dur = s1 - s0
        rows.append({"cell_id": cell_id, "genotype": geno, "sweep_idx": sw,
                     "stim_pA": round(amp, 1), "stim_start_s": round(s0, 4),
                     "stim_end_s": round(s1, 4), "stim_dur_s": round(dur, 4),
                     "n_spikes": n, "freq_Hz": n / dur if dur > 0 else 0.0})

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / f"{cell_id}_allen_fi.csv", index=False)
    return out


def main():
    n_cells = 0
    n_skipped = 0
    for cd in sorted(REC_DIR.iterdir()):
        if not cd.is_dir():
            continue
        out = process_cell(cd)
        if out is None:
            n_skipped += 1
            continue
        max_amp = out["stim_pA"].max()
        max_freq = out["freq_Hz"].max()
        print(f"  {cd.name}: {len(out)} sweeps, max stim={max_amp:.0f} pA, max freq={max_freq:.1f} Hz")
        n_cells += 1
    print(f"\nDone. Processed {n_cells} cells, skipped {n_skipped}.")


if __name__ == "__main__":
    main()
