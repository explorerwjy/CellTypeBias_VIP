"""
Generate 12 trace PDFs for both protocols × 3 conditions × 2 genotypes:

  Protocols:
    - Intrinsic_IV curve_VIP_Erica  (200 ms pulses, F-I curve)
    - Intrinsic_IV curve_VIP        (1500 ms pulses, long-pulse / IV)

  Conditions per protocol:
    - rheobase        — smallest stim that fires ≥ 1 spike
    - hero            — closest sweep to rheobase + 40 pA (within ±10 pA, then ±20)
    - 150pA           — closest sweep to +150 pA (within ±15 pA)

Cells: those Erica included in her analysis (30 cells from
Stephanie_Data_PrePrism.xlsx with non-NaN Half-width and Rise Slope).

Output: dat/VIP_Ephys/trace_pdfs/{WT|Df16}_{rheobase|hero|150pA}_{200ms|1500ms}.pdf
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

PROTOCOL_TAGS = {
    "200ms":  ["Intrinsic_IV curve_VIP_Erica", "Intrinsic_Thrsehold_VIP_Erica"],
    "1500ms": ["Intrinsic_IV curve_VIP"],
}


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


def load_abf(cell_dir, allowed_protocols):
    cands = []
    for f in sorted(cell_dir.glob("*.abf")):
        try: a = pyabf.ABF(str(f))
        except Exception: continue
        if a.protocol in allowed_protocols:
            cands.append((allowed_protocols.index(a.protocol), a))
    if not cands: return None
    cands.sort(key=lambda x: x[0])
    return cands[0][1]


def collect_sweep_info(abf):
    info = []
    for sw in range(abf.sweepCount):
        abf.setSweep(sw, channel=0); v = abf.sweepY.astype(float); t = abf.sweepX.astype(float)
        abf.setSweep(sw, channel=1); i = abf.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None: continue
        s, e, amp = win
        ext = SpikeFeatureExtractor(start=t[s] - 0.005, end=min(t[-1], t[e] + 0.005),
                                    filter=4.0)
        try:
            spk = ext.process(t, v, i)
            in_w = (spk["threshold_t"] >= t[s]) & (spk["threshold_t"] <= t[e])
            n_sp = int(in_w.sum()) if len(spk) > 0 else 0
        except Exception:
            n_sp = 0
        info.append({"sweep": sw, "stim": amp, "s": s, "e": e,
                     "n_sp": n_sp, "t": t, "v": v})
    return info


def select_condition(sweep_info, condition):
    """Return chosen sweep dict for the requested condition, or None."""
    if not sweep_info: return None
    df = pd.DataFrame([{k: r[k] for k in ("sweep", "stim", "n_sp")} for r in sweep_info])

    if condition == "rheobase":
        spiking = df[df["n_sp"] > 0].sort_values("stim")
        if spiking.empty: return None
        sw_idx = int(spiking.iloc[0]["sweep"])
        return next(r for r in sweep_info if r["sweep"] == sw_idx)

    if condition == "hero":
        spiking = df[df["n_sp"] > 0].sort_values("stim")
        if spiking.empty: return None
        rheo = float(spiking.iloc[0]["stim"])
        target = rheo + 40
        cand = df[(df["stim"] >= rheo + 30) & (df["stim"] <= rheo + 50)]
        if cand.empty:
            cand = df[(df["stim"] >= rheo + 20) & (df["stim"] <= rheo + 60)]
        if cand.empty: return None
        sw_idx = int(cand.iloc[(cand["stim"] - target).abs().argmin()]["sweep"])
        return next(r for r in sweep_info if r["sweep"] == sw_idx)

    if condition == "150pA":
        cand = df.iloc[(df["stim"] - 150).abs().argmin()]
        if abs(float(cand["stim"]) - 150) > 15: return None
        sw_idx = int(cand["sweep"])
        return next(r for r in sweep_info if r["sweep"] == sw_idx)

    return None


def plot_grid(cells, genotype, condition, protocol_tag, out_path):
    valid = [(c, d) for c, d in cells if d is not None]
    if not valid:
        print(f"  ! No valid cells for {genotype}/{condition}/{protocol_tag}; skipping")
        return
    color = "#2166ac" if genotype == "WT" else "#b2182b"
    ncols = 3
    nrows = int(np.ceil(len(valid) / ncols))
    with PdfPages(out_path) as pdf:
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.7 * nrows), squeeze=False)
        for k, (cid, d) in enumerate(valid):
            ax = axes[k // ncols][k % ncols]
            t = d["t"]; v = d["v"]
            stim_t0, stim_t1 = t[d["s"]], t[d["e"]]
            mask = (t >= stim_t0 - 0.05) & (t <= stim_t1 + 0.15)
            ax.plot((t[mask] - stim_t0) * 1000, v[mask], color=color, linewidth=0.65)
            ax.axvspan(0, (stim_t1 - stim_t0) * 1000, color="gray", alpha=0.10, linewidth=0)
            ax.set_title(f"{cid}\n{d['stim']:.0f} pA, {d['n_sp']} spikes", fontsize=8)
            ax.set_xlabel("Time from stim onset (ms)", fontsize=8)
            ax.set_ylabel("Vm (mV)", fontsize=8)
            ax.tick_params(labelsize=7)
            for s in ("top", "right"): ax.spines[s].set_visible(False)
        for k in range(len(valid), nrows * ncols):
            axes[k // ncols][k % ncols].axis("off")
        title = f"{genotype} — {condition} sweep, {protocol_tag} protocol  (n = {len(valid)} cells)"
        fig.suptitle(title, fontsize=12, fontweight="bold", y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.985])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
    print(f"  saved {out_path.name}  ({len(valid)} cells)")


def main():
    excel = PROJ / "dat" / "VIP_Ephys" / "erica_prism" / "raw_features" / "Stephanie_Data_PrePrism.xlsx"
    df = pd.read_excel(excel, sheet_name="Intrinsic Properties")
    erica = df[df["Half-width (ms)"].notna() & df["Rise (V/s)"].notna()]

    for proto_tag, proto_list in PROTOCOL_TAGS.items():
        print(f"\n=== Protocol: {proto_tag} ({proto_list[0]}) ===")
        # Cache sweep_info per cell to avoid re-reading
        cell_info = {}
        for _, row in erica.iterrows():
            cd, geno = map_erica_to_cell_dir(row["Group"], row["Mouse"], row["Cell"])
            if cd is None: continue
            abf = load_abf(cd, proto_list)
            if abf is None:
                cell_info[cd.name] = (geno, None); continue
            cell_info[cd.name] = (geno, collect_sweep_info(abf))

        for condition in ("rheobase", "hero", "150pA"):
            wt_data, df16_data = [], []
            for cid, (geno, info) in cell_info.items():
                if info is None:
                    chosen = None
                else:
                    chosen = select_condition(info, condition)
                if geno == "WT": wt_data.append((cid, chosen))
                else: df16_data.append((cid, chosen))
            plot_grid(wt_data,   "WT",   condition, proto_tag,
                      OUT_DIR / f"WT_{condition}_{proto_tag}.pdf")
            plot_grid(df16_data, "Df16", condition, proto_tag,
                      OUT_DIR / f"Df16_{condition}_{proto_tag}.pdf")


if __name__ == "__main__":
    main()
