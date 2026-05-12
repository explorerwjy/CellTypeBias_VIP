"""
Generate 4 PDFs of voltage traces, one per (genotype, sweep type):
  - WT_rheobase.pdf
  - WT_hero.pdf       (rheobase + 40 ± 10 pA)
  - Df16_rheobase.pdf
  - Df16_hero.pdf

Cells: those Erica included in her analysis (30 cells, from the
Stephanie_Data_PrePrism.xlsx 'Intrinsic Properties' sheet, with non-NaN
Half-width and Rise Slope).

Sweep selection: F-I protocol "Intrinsic_IV curve_VIP_Erica".
  - Rheobase = smallest stim amplitude that fires ≥1 spike (Allen detection).
  - Hero = closest sweep with stim ∈ [rheobase + 30, rheobase + 50] pA;
           fallback widens to ±20 pA from rheobase + 40.
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
PROTOCOLS = (
    "Intrinsic_IV curve_VIP_Erica",
    "Intrinsic_Thrsehold_VIP_Erica",
    "Intrinsic_IV curve_VIP",  # fallback — sometimes the same name delivers F-I-range
)


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


def map_erica_to_cell_dir(group, mouse, cell):
    """Map Erica's (Group, Mouse, Cell) to our directory name."""
    geno = "WT" if str(group).strip().upper() == "A" else "Df16"
    mouse_str = str(mouse).strip()
    if mouse_str == "?":
        mouse_id = "A10916"
    elif mouse_str == "??":
        mouse_id = "B10916"
    else:
        try:
            mouse_id = f"{int(float(mouse_str)):02d}"
        except Exception:
            return None
    cell_id = f"{int(cell):02d}"
    cand = f"{geno}_mouse{mouse_id}_cell{cell_id}"
    cd = REC_DIR / cand
    if cd.exists():
        return cd
    # try with suffixes (e.g., _excluded, _lostcell)
    for d in REC_DIR.iterdir():
        if d.name.startswith(cand):
            return d
    return None


def get_cell_sweeps(cell_dir):
    """Return: rheo_data, hero_data dicts; falls back through PROTOCOLS in order.

    Picks the protocol that has the widest stim range with at least one spiking
    sweep, so we maximize the chance of finding hero (rheobase + 40 pA).
    """
    candidates = []
    for f in sorted(cell_dir.glob("*.abf")):
        try:
            a = pyabf.ABF(str(f))
        except Exception:
            continue
        if a.protocol not in PROTOCOLS:
            continue
        candidates.append((PROTOCOLS.index(a.protocol), a))
    if not candidates:
        return None, None
    # Prefer protocols by priority (earlier in PROTOCOLS list = preferred)
    candidates.sort(key=lambda x: x[0])
    abf_file = candidates[0][1]
    used_protocol = abf_file.protocol

    sweep_info = []
    for sw in range(abf_file.sweepCount):
        abf_file.setSweep(sw, channel=0)
        v = abf_file.sweepY.astype(float)
        t = abf_file.sweepX.astype(float)
        abf_file.setSweep(sw, channel=1)
        i = abf_file.sweepY.astype(float)
        win = detect_stim_window(t, i)
        if win is None:
            continue
        s, e, amp = win
        ext = SpikeFeatureExtractor(start=t[s] - 0.005, end=min(t[-1], t[e] + 0.005),
                                    filter=4.0)
        try:
            spk = ext.process(t, v, i)
            in_win = (spk["threshold_t"] >= t[s]) & (spk["threshold_t"] <= t[e])
            n_sp = int(in_win.sum()) if len(spk) > 0 else 0
        except Exception:
            n_sp = 0
        sweep_info.append({"sweep": sw, "stim": amp, "s": s, "e": e,
                           "n_sp": n_sp, "t": t.copy(), "v": v.copy()})

    if not sweep_info:
        return None, None
    sw_df = pd.DataFrame([{k: r[k] for k in ("sweep", "stim", "s", "e", "n_sp")}
                          for r in sweep_info])
    spiking = sw_df[sw_df["n_sp"] > 0].sort_values("stim")
    if spiking.empty:
        return None, None
    rheo_amp = float(spiking.iloc[0]["stim"])
    rheo_sw = int(spiking.iloc[0]["sweep"])
    rheo_data = next(r for r in sweep_info if r["sweep"] == rheo_sw)

    # Hero: closest sweep with stim ∈ [rheo + 30, rheo + 50]
    target = rheo_amp + 40
    cand = sw_df[(sw_df["stim"] >= rheo_amp + 30) & (sw_df["stim"] <= rheo_amp + 50)]
    if cand.empty:
        cand = sw_df[(sw_df["stim"] >= rheo_amp + 20) & (sw_df["stim"] <= rheo_amp + 60)]
    if cand.empty:
        return rheo_data, None
    hero_sw = int(cand.iloc[(cand["stim"] - target).abs().argmin()]["sweep"])
    hero_data = next(r for r in sweep_info if r["sweep"] == hero_sw)

    return rheo_data, hero_data


def plot_pdf(cell_data, sweep_type, genotype, out_path):
    """Grid plot of voltage traces, one panel per cell."""
    valid = [(cid, d) for cid, d in cell_data if d is not None]
    n = len(valid)
    if n == 0:
        print(f"  ! No cells for {genotype} / {sweep_type}; skipping")
        return
    # 3 columns, ceil(n/3) rows
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    color = "#2166ac" if genotype == "WT" else "#b2182b"

    with PdfPages(out_path) as pdf:
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.6 * nrows),
                                 squeeze=False)
        for ax_idx, (cell_id, data) in enumerate(valid):
            ax = axes[ax_idx // ncols][ax_idx % ncols]
            t = data["t"]
            v = data["v"]
            stim_t0 = t[data["s"]]
            stim_t1 = t[data["e"]]
            # show 100 ms before and 200 ms after stim
            xlim_lo = max(t[0], stim_t0 - 0.1)
            xlim_hi = min(t[-1], stim_t1 + 0.2)
            mask = (t >= xlim_lo) & (t <= xlim_hi)
            ax.plot((t[mask] - stim_t0) * 1000, v[mask],
                    color=color, linewidth=0.7)
            # stim window shading
            ax.axvspan(0, (stim_t1 - stim_t0) * 1000,
                       color="gray", alpha=0.10, linewidth=0)
            ax.set_title(f"{cell_id}\n{data['stim']:.0f} pA, {data['n_sp']} spk",
                         fontsize=9)
            ax.set_xlabel("Time from stim onset (ms)", fontsize=8)
            ax.set_ylabel("Vm (mV)", fontsize=8)
            ax.tick_params(labelsize=7)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        # blank unused axes
        for k in range(len(valid), nrows * ncols):
            axes[k // ncols][k % ncols].axis("off")

        title = f"{genotype} — {'rheobase' if sweep_type == 'rheo' else 'hero (rheobase + 40 ± 10 pA)'} sweep"
        title += f"  (n = {len(valid)} cells)"
        fig.suptitle(title, fontsize=12, fontweight="bold", y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.985])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
    print(f"  saved {out_path}  ({len(valid)} cells)")


def main():
    # 1) Erica's included cells from the Excel
    excel = PROJ / "dat" / "VIP_Ephys" / "erica_prism" / "raw_features" / "Stephanie_Data_PrePrism.xlsx"
    df = pd.read_excel(excel, sheet_name="Intrinsic Properties")
    erica = df[df["Half-width (ms)"].notna() & df["Rise (V/s)"].notna()].copy()
    erica["geno"] = erica["Group"].apply(lambda g: "WT" if str(g).strip().upper() == "A" else "Df16")

    # 2) Map and load
    wt_rheo, wt_hero = [], []
    df16_rheo, df16_hero = [], []
    for _, row in erica.iterrows():
        cd = map_erica_to_cell_dir(row["Group"], row["Mouse"], row["Cell"])
        if cd is None:
            print(f"  ! Could not map: {row['Group']} mouse {row['Mouse']} cell {row['Cell']}")
            continue
        rheo, hero = get_cell_sweeps(cd)
        bucket_rheo = wt_rheo if row["geno"] == "WT" else df16_rheo
        bucket_hero = wt_hero if row["geno"] == "WT" else df16_hero
        bucket_rheo.append((cd.name, rheo))
        bucket_hero.append((cd.name, hero))
        rh_str = f"rheo={rheo['stim']:.0f}pA" if rheo else "rheo=N/A"
        he_str = f"hero={hero['stim']:.0f}pA ({hero['n_sp']} spk)" if hero else "hero=N/A"
        print(f"  [{row['geno']:4s}] {cd.name}  {rh_str}, {he_str}")

    # 3) Plot 4 PDFs
    plot_pdf(wt_rheo,   "rheo", "WT",   OUT_DIR / "WT_rheobase.pdf")
    plot_pdf(wt_hero,   "hero", "WT",   OUT_DIR / "WT_hero.pdf")
    plot_pdf(df16_rheo, "rheo", "Df16", OUT_DIR / "Df16_rheobase.pdf")
    plot_pdf(df16_hero, "hero", "Df16", OUT_DIR / "Df16_hero.pdf")


if __name__ == "__main__":
    main()
