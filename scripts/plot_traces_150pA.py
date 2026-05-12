"""
Plot every Erica-included cell's trace at the sweep closest to +150 pA, the
same current shown in manuscript Figure 4C. Two PDFs (WT and Df16), one
panel per cell, annotated with cell ID, actual stim_pA, spike count, and
burstiness metrics (n_bursts, LV).

Output:
  dat/VIP_Ephys/trace_pdfs/WT_150pA.pdf
  dat/VIP_Ephys/trace_pdfs/Df16_150pA.pdf
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
# Prefer LONG-PULSE protocols for burst observation (1.5 s steps)
# over the short-pulse F-I protocol (200 ms). Manuscript Figure 4C
# representative traces come from the long-pulse IV protocol.
PROTOCOLS = ("Intrinsic_IV curve_VIP",            # 1.5 s pulses (preferred)
             "Intrinsic_IV curve_VIP_Erica",       # 200 ms F-I (fallback)
             "Intrinsic_Thrsehold_VIP_Erica")      # 50 ms threshold (last)
TARGET_CURRENT_PA = 150.0


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


def map_erica_to_cell_dir(group, mouse, cell):
    geno = "WT" if str(group).strip().upper() == "A" else "Df16"
    mouse_str = str(mouse).strip()
    if mouse_str == "?": mouse_id = "A10916"
    elif mouse_str == "??": mouse_id = "B10916"
    else:
        try: mouse_id = f"{int(float(mouse_str)):02d}"
        except Exception: return None, geno
    cand = f"{geno}_mouse{mouse_id}_cell{int(cell):02d}"
    cd = REC_DIR / cand
    if cd.exists(): return cd, geno
    for d in REC_DIR.iterdir():
        if d.name.startswith(cand): return d, geno
    return None, geno


def load_abf(cell_dir, target_pA=TARGET_CURRENT_PA, tol_pA=15):
    """Pick the protocol that ACTUALLY contains a sweep near target_pA, by priority."""
    cands = []
    for f in sorted(cell_dir.glob("*.abf")):
        try: a = pyabf.ABF(str(f))
        except Exception: continue
        if a.protocol in PROTOCOLS:
            cands.append((PROTOCOLS.index(a.protocol), a))
    cands.sort(key=lambda x: x[0])
    # First try in priority order, picking the protocol that actually has a sweep at target
    for _, a in cands:
        for sw in range(a.sweepCount):
            a.setSweep(sw, channel=1)
            i = a.sweepY.astype(float); t = a.sweepX.astype(float)
            win = detect_stim_window(t, i)
            if win is None: continue
            _, _, amp = win
            if abs(amp - target_pA) <= tol_pA:
                return a
    # If no protocol has the target, return the highest-priority one anyway
    return cands[0][1] if cands else None


def extract_target_sweep(abf, target_pA, tol_pA=15):
    """Find sweep with stim closest to target_pA; return None if no sweep within tol."""
    sw_info = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0)
        v = abf.sweepY.astype(float); t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1)
        i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None: continue
        s, e, amp = win
        sw_info.append({"sweep": sw, "stim": amp, "s": s, "e": e,
                        "t": t.copy(), "v": v.copy()})
    if not sw_info: return None
    df = pd.DataFrame([{k: r[k] for k in ("sweep", "stim", "s", "e")} for r in sw_info])
    closest = df.iloc[(df["stim"] - target_pA).abs().argmin()]
    if abs(closest["stim"] - target_pA) > tol_pA:
        return None
    return next(r for r in sw_info if r["sweep"] == int(closest["sweep"]))


def burst_metrics(t, v, i, s_idx, e_idx):
    """Run Allen IPFX detection on the sweep, compute burst metrics."""
    abf_t = np.asarray(t); abf_v = np.asarray(v); abf_i = np.asarray(i)
    ext = SpikeFeatureExtractor(start=abf_t[s_idx] - 0.005,
                                end=min(abf_t[-1], abf_t[e_idx] + 0.005),
                                filter=4.0)
    try:
        spk = ext.process(abf_t, abf_v, abf_i)
        in_win = (spk["threshold_t"] >= abf_t[s_idx]) & (spk["threshold_t"] <= abf_t[e_idx])
        spk = spk[in_win].reset_index(drop=True)
    except Exception:
        spk = pd.DataFrame()
    n = len(spk)
    if n < 2:
        return {"n_spikes": n, "lv": np.nan, "isi_cv": np.nan, "n_bursts": 0,
                "spike_t": spk["threshold_t"].values if n > 0 else np.array([])}
    isi = np.diff(spk["threshold_t"].values)
    lv = float(np.mean(3 * (isi[:-1] - isi[1:])**2 / (isi[:-1] + isi[1:])**2)) if len(isi) >= 2 else np.nan
    cv = float(isi.std() / isi.mean()) if isi.mean() > 0 else np.nan
    burst_thr = min(0.025, np.median(isi) / 2)
    in_burst = isi < burst_thr
    n_bursts = 0; run = 0
    for vv in in_burst:
        if vv: run += 1
        else:
            if run >= 2: n_bursts += 1
            run = 0
    if run >= 2: n_bursts += 1
    return {"n_spikes": n, "lv": lv, "isi_cv": cv, "n_bursts": n_bursts,
            "spike_t": spk["threshold_t"].values}


def main():
    excel = PROJ / "dat" / "VIP_Ephys" / "erica_prism" / "raw_features" / "Stephanie_Data_PrePrism.xlsx"
    df = pd.read_excel(excel, sheet_name="Intrinsic Properties")
    erica = df[df["Half-width (ms)"].notna() & df["Rise (V/s)"].notna()]

    wt_data, df16_data = [], []
    for _, row in erica.iterrows():
        cd, geno = map_erica_to_cell_dir(row["Group"], row["Mouse"], row["Cell"])
        if cd is None: continue
        abf = load_abf(cd)
        if abf is None: continue
        sw = extract_target_sweep(abf, TARGET_CURRENT_PA, tol_pA=15)
        if sw is None:
            print(f"  [{geno:4s}] {cd.name}: NO sweep within ±15 pA of 150 pA")
            (wt_data if geno == "WT" else df16_data).append((cd.name, None))
            continue
        bm = burst_metrics(sw["t"], sw["v"], None, sw["s"], sw["e"])
        # actually we need current to invoke IPFX correctly — re-load i channel
        abf.setSweep(sw["sweep"], channel=1)
        i_ch = abf.sweepY.astype(float)
        bm = burst_metrics(sw["t"], sw["v"], i_ch, sw["s"], sw["e"])
        sw.update(bm)
        sw["stim"] = float(sw["stim"])
        print(f"  [{geno:4s}] {cd.name}: {sw['stim']:.1f} pA, "
              f"{sw['n_spikes']} spikes, #bursts={sw['n_bursts']}, LV={sw['lv']:.2f}")
        (wt_data if geno == "WT" else df16_data).append((cd.name, sw))

    plot_pdf(wt_data,   "WT",   OUT_DIR / "WT_150pA.pdf")
    plot_pdf(df16_data, "Df16", OUT_DIR / "Df16_150pA.pdf")


def plot_pdf(cell_data, genotype, out_path):
    valid = [(c, d) for c, d in cell_data if d is not None]
    if not valid:
        print(f"  ! No cells for {genotype}; skipping {out_path}"); return
    color = "#2166ac" if genotype == "WT" else "#b2182b"
    ncols = 3
    nrows = int(np.ceil(len(valid) / ncols))
    with PdfPages(out_path) as pdf:
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.7 * nrows), squeeze=False)
        for ax_idx, (cid, d) in enumerate(valid):
            ax = axes[ax_idx // ncols][ax_idx % ncols]
            t = d["t"]; v = d["v"]
            stim_t0, stim_t1 = t[d["s"]], t[d["e"]]
            mask = (t >= stim_t0 - 0.05) & (t <= stim_t1 + 0.1)
            ax.plot((t[mask] - stim_t0) * 1000, v[mask], color=color, linewidth=0.7)
            ax.axvspan(0, (stim_t1 - stim_t0) * 1000, color="gray", alpha=0.10, linewidth=0)
            ax.set_title(
                f"{cid}\n{d['stim']:.0f} pA  |  {d['n_spikes']} spk, "
                f"#bursts={d['n_bursts']}, LV={d['lv']:.2f}", fontsize=8)
            ax.set_xlabel("Time from stim onset (ms)", fontsize=8)
            ax.set_ylabel("Vm (mV)", fontsize=8)
            ax.tick_params(labelsize=7)
            for s in ("top", "right"): ax.spines[s].set_visible(False)
        for k in range(len(valid), nrows * ncols):
            axes[k // ncols][k % ncols].axis("off")
        title = f"{genotype} — sweep closest to +150 pA  (n = {len(valid)} cells)"
        fig.suptitle(title, fontsize=12, fontweight="bold", y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.985])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
    print(f"  saved {out_path}")


if __name__ == "__main__":
    main()
