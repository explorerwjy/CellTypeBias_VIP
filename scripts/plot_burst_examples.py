"""
Plot example traces at the max-firing sweep for cells with the highest #bursts
in WT and Df16, so the burstiness comparison can be judged by eye.

Output: dat/VIP_Ephys/trace_pdfs/burst_examples_max_sweep.pdf
"""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pyabf
from ipfx.feature_extractor import SpikeFeatureExtractor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
REC_DIR = PROJ / "dat" / "VIP_Ephys" / "IV_step_recordings"
OUT_DIR = PROJ / "dat" / "VIP_Ephys" / "trace_pdfs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SPIKE_DIR = PROJ / "dat" / "VIP_Ephys" / "allen_ipfx_vip" / "features_fi_perspike"

PROTOCOLS = ("Intrinsic_IV curve_VIP_Erica",
             "Intrinsic_Thrsehold_VIP_Erica",
             "Intrinsic_IV curve_VIP")


def detect_stim_window(t, i, threshold_pA=2.0):
    baseline = np.percentile(i, 25)
    above = np.abs(i - baseline) > threshold_pA
    if not above.any(): return None
    edges = np.diff(above.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1
    if above[0]: starts = np.concatenate([[0], starts])
    if above[-1]: ends = np.concatenate([ends, [len(above)]])
    runs = list(zip(starts, ends))
    if not runs: return None
    s, e = max(runs, key=lambda r: r[1] - r[0])
    e = min(e, len(t) - 1)
    return s, e, float(np.median(i[s:e] - baseline))


def cell_burstiness_summary(cell_csv):
    df = pd.read_csv(cell_csv)
    if df.empty: return None
    cell_id = df.cell_id.iloc[0]
    geno = df.genotype.iloc[0]
    sw_counts = (df.groupby("sweep_idx")
                   .agg(stim=("stim_pA", "first"), n=("spike_idx_in_sweep", "count"))
                   .reset_index())
    sw_counts = sw_counts[sw_counts.n >= 5]
    if sw_counts.empty: return None
    rheo = float(sw_counts.sort_values("stim").iloc[0].stim)
    sw_idx = int(sw_counts.loc[sw_counts.n.idxmax()].sweep_idx)
    spike_t = df[df.sweep_idx == sw_idx].sort_values("threshold_t").threshold_t.values
    isi = np.diff(spike_t)
    if len(isi) < 2: return None
    burst_thr = min(0.025, np.median(isi) / 2)
    in_burst = isi < burst_thr
    n_bursts = 0; run = 0
    for v in in_burst:
        if v: run += 1
        else:
            if run >= 2: n_bursts += 1
            run = 0
    if run >= 2: n_bursts += 1
    return {"cell_id": cell_id, "genotype": geno, "rheo": rheo, "sweep_idx": sw_idx,
            "stim_pA": float(sw_counts.loc[sw_counts.sweep_idx == sw_idx, "stim"].iloc[0]),
            "n_spikes": int(len(spike_t)),
            "n_bursts": int(n_bursts),
            "lv": float(np.mean(3 * (isi[:-1] - isi[1:])**2 / (isi[:-1] + isi[1:])**2)),
            "isi_cv": float(isi.std() / isi.mean()),
            "spike_times": spike_t.tolist(),
            "isi_first_two": isi[:2].tolist() if len(isi) >= 2 else None}


def load_sweep_trace(cell_id, sweep_idx):
    cd = REC_DIR / cell_id
    if not cd.exists():
        for d in REC_DIR.iterdir():
            if d.name.startswith(cell_id):
                cd = d; break
    candidates = []
    for f in sorted(cd.glob("*.abf")):
        try: a = pyabf.ABF(str(f))
        except Exception: continue
        if a.protocol in PROTOCOLS:
            candidates.append((PROTOCOLS.index(a.protocol), a))
    if not candidates: return None
    candidates.sort(key=lambda x: x[0])
    abf = candidates[0][1]
    abf.setSweep(sweep_idx, channel=0)
    t = abf.sweepX.astype(float); v = abf.sweepY.astype(float)
    abf.setSweep(sweep_idx, channel=1)
    i = abf.sweepY.astype(float)
    win = detect_stim_window(t, i)
    if win is None: return t, v, None, None
    s, e, _ = win
    return t, v, t[s], t[e]


def plot_examples():
    summaries = []
    for f in sorted(SPIKE_DIR.glob("*_fi_spikes.csv")):
        if "lostcell" in f.name: continue
        s = cell_burstiness_summary(f)
        if s is None: continue
        summaries.append(s)
    df = pd.DataFrame(summaries)
    print(df[["cell_id", "genotype", "n_spikes", "n_bursts", "lv", "isi_cv", "stim_pA"]]
          .sort_values(["genotype", "n_bursts", "lv"], ascending=[True, False, False])
          .to_string(index=False))

    # Pick: top-3 burstiest WT and top-3 most-irregular Df16 (highest LV)
    wt = df[df.genotype == "WT"].sort_values(["n_bursts", "lv", "n_spikes"],
                                             ascending=[False, False, False]).head(3)
    d16 = df[df.genotype == "Df16"].sort_values(["lv", "n_bursts", "n_spikes"],
                                                ascending=[False, False, False]).head(3)

    fig, axes = plt.subplots(3, 2, figsize=(13, 9))
    color = {"WT": "#2166ac", "Df16": "#b2182b"}
    for col_idx, (label, sub) in enumerate([("WT (top-3 burstiest)", wt),
                                             ("Df(16) (top-3 most irregular)", d16)]):
        for row_idx, (_, r) in enumerate(sub.reset_index().iterrows()):
            ax = axes[row_idx, col_idx]
            res = load_sweep_trace(r["cell_id"], int(r["sweep_idx"]))
            if res is None or res[0] is None:
                ax.set_title(f"{r['cell_id']} (no trace)"); continue
            t, v, s_t, e_t = res
            if s_t is None:
                ax.plot(t, v, color=color[r["genotype"]], linewidth=0.7)
            else:
                m = (t >= s_t - 0.05) & (t <= e_t + 0.1)
                ax.plot((t[m] - s_t) * 1000, v[m], color=color[r["genotype"]], linewidth=0.7)
                ax.axvspan(0, (e_t - s_t) * 1000, color="gray", alpha=0.10, linewidth=0)
            ax.set_title(
                f"{r['cell_id']}\n{r['stim_pA']:.0f} pA, {r['n_spikes']} spikes, "
                f"#bursts={r['n_bursts']}, LV={r['lv']:.2f}", fontsize=9)
            ax.set_xlabel("Time from stim onset (ms)", fontsize=8)
            ax.set_ylabel("Vm (mV)", fontsize=8)
            ax.tick_params(labelsize=7)
            for s in ("top", "right"): ax.spines[s].set_visible(False)
        axes[0, col_idx].text(0.5, 1.40, label, ha="center", fontsize=12, fontweight="bold",
                              transform=axes[0, col_idx].transAxes)

    fig.suptitle("Burstiness examples at max-firing sweep (Allen IPFX detection)",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = OUT_DIR / "burst_examples_max_sweep.pdf"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    plot_examples()
