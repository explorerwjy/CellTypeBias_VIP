"""
Replicate manuscript Figure 4D (firing frequency vs current step) using
Allen IPFX spike detection on Protocol 2 (200 ms pulses, Intrinsic_IV curve_VIP_Erica).

Compare to Erica's reported numbers (Current x Frequency.csv) on the same axes.

Two outputs:
  - notebooks_rebuttal/Figure4D_replicate_Allen.png     — Allen detection only
  - notebooks_rebuttal/Figure4D_replicate_compare.png   — Allen + Erica overlaid
"""
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
DATA = PROJ / "dat" / "VIP_Ephys"
OUT = PROJ / "notebooks_rebuttal"

# Manuscript's color palette (matched from Figure 4D image)
COL = {"WT": "#4a8a82", "Df16": "#b2557e"}


def load_allen_fi():
    files = sorted((DATA / "allen_ipfx_vip" / "features_fi").glob("*_allen_fi.csv"))
    rows = []
    for f in files:
        if "lostcell" in f.name:
            continue
        df = pd.read_csv(f)
        rows.append(df)
    sweep = pd.concat(rows, ignore_index=True)
    # Bin to 10 pA grid (matching Erica)
    sweep["cur_bin"] = (np.round(sweep["stim_pA"] / 10) * 10).astype("Int64")
    sweep = sweep.dropna(subset=["cur_bin"])
    sweep["cur_bin"] = sweep["cur_bin"].astype(int)
    fi = sweep.groupby(["cell_id", "genotype", "cur_bin"]).agg(
        freq_Hz=("freq_Hz", "mean")).reset_index()
    return fi[(fi.cur_bin >= 0) & (fi.cur_bin <= 200)]


def load_erica_fi():
    fn = DATA / "erica_prism" / "raw_features" / "Current x Frequency.csv"
    df = pd.read_csv(fn, encoding="utf-8-sig", header=None)
    headers = df.iloc[0].tolist()
    body = df.iloc[1:].reset_index(drop=True)
    current = pd.to_numeric(body.iloc[:, 0], errors="coerce").values
    rec = []
    for ci, h in enumerate(headers[1:], 1):
        h_str = str(h).strip()
        if h_str.upper().startswith("WT"):
            geno = "WT"
        elif "DF" in h_str.upper() and "16" in h_str:
            geno = "Df16"
        else:
            continue
        col = pd.to_numeric(body.iloc[:, ci], errors="coerce").values
        for cur, v in zip(current, col):
            if not np.isnan(v) and not np.isnan(cur):
                rec.append({"cell": f"c{ci}", "genotype": geno,
                            "cur_bin": int(cur), "freq_Hz": float(v)})
    fi = pd.DataFrame(rec)
    return fi[(fi.cur_bin >= 0) & (fi.cur_bin <= 200)]


def per_bin_stats(fi, cur_col="cur_bin", freq_col="freq_Hz"):
    """Compute mean, SEM, n, MWU p per current bin per genotype."""
    rows = []
    for cur in sorted(fi[cur_col].unique()):
        wt = fi[(fi[cur_col] == cur) & (fi.genotype == "WT")][freq_col].values
        df16 = fi[(fi[cur_col] == cur) & (fi.genotype == "Df16")][freq_col].values
        p = np.nan
        if len(wt) >= 3 and len(df16) >= 3:
            try:
                _, p = mannwhitneyu(wt, df16, alternative="two-sided")
            except Exception:
                p = np.nan
        rows.append({
            "cur": cur,
            "wt_mean": wt.mean() if len(wt) else np.nan,
            "wt_sem": wt.std(ddof=1) / np.sqrt(len(wt)) if len(wt) > 1 else 0,
            "n_wt": len(wt),
            "df16_mean": df16.mean() if len(df16) else np.nan,
            "df16_sem": df16.std(ddof=1) / np.sqrt(len(df16)) if len(df16) > 1 else 0,
            "n_df16": len(df16),
            "mwu_p": p,
        })
    return pd.DataFrame(rows)


def plot_fi(stats_df, ax, title, max_n_wt, max_n_df16, show_n_label=True):
    cur = stats_df["cur"].values
    # Shaded SEM band + line
    for geno, mean_col, sem_col in [("WT", "wt_mean", "wt_sem"),
                                     ("Df16", "df16_mean", "df16_sem")]:
        m = stats_df[mean_col].values
        s = stats_df[sem_col].values
        ax.fill_between(cur, m - s, m + s, color=COL[geno], alpha=0.18,
                        linewidth=0)
        label = f"{'Wildtype' if geno == 'WT' else 'Df(16)A+/-'} (n={max_n_wt if geno == 'WT' else max_n_df16})" if show_n_label else None
        ax.plot(cur, m, color=COL[geno], linewidth=2.0, marker="o",
                markersize=4, label=label)

    # Significance markers — asterisks at p < 0.05 bins
    ymax = max(np.nanmax(stats_df["wt_mean"] + stats_df["wt_sem"]),
               np.nanmax(stats_df["df16_mean"] + stats_df["df16_sem"]))
    for _, r in stats_df.iterrows():
        if pd.notna(r["mwu_p"]) and r["mwu_p"] < 0.05:
            ax.text(r["cur"], ymax * 1.05, "*", ha="center",
                    fontsize=12, fontweight="bold", color="black")

    ax.set_xlabel("Current (pA)")
    ax.set_ylabel("AP frequency (Hz)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlim(-5, 205)
    ax.set_ylim(-2, ymax * 1.15)


def main():
    allen = load_allen_fi()
    erica = load_erica_fi()
    a_stats = per_bin_stats(allen)
    e_stats = per_bin_stats(erica)

    print("Allen IPFX detection cohort (200 ms protocol):")
    print(allen.groupby("genotype")["cell_id"].nunique())
    print()
    print("Erica reported cohort:")
    print(erica.groupby("genotype")["cell"].nunique())
    print()

    # Print per-bin MWU comparison
    print(f'{"current(pA)":>11s} {"Allen WT m":>10s} {"Df16 m":>9s} {"MWU p":>9s}  {"Erica WT m":>10s} {"Df16 m":>9s} {"MWU p":>9s}')
    for c in range(0, 210, 10):
        a = a_stats[a_stats.cur == c]
        e = e_stats[e_stats.cur == c]
        ap = f"{a.iloc[0].mwu_p:.3f}" if len(a) and pd.notna(a.iloc[0].mwu_p) else "--"
        ep = f"{e.iloc[0].mwu_p:.3f}" if len(e) and pd.notna(e.iloc[0].mwu_p) else "--"
        a_wt = f"{a.iloc[0].wt_mean:6.2f}" if len(a) else "--"
        a_d = f"{a.iloc[0].df16_mean:6.2f}" if len(a) else "--"
        e_wt = f"{e.iloc[0].wt_mean:6.2f}" if len(e) else "--"
        e_d = f"{e.iloc[0].df16_mean:6.2f}" if len(e) else "--"
        print(f"{c:11d} {a_wt:>10s} {a_d:>9s} {ap:>9s}  {e_wt:>10s} {e_d:>9s} {ep:>9s}")

    n_wt_a = allen[allen.genotype == "WT"]["cell_id"].nunique()
    n_d_a = allen[allen.genotype == "Df16"]["cell_id"].nunique()
    n_wt_e = erica[erica.genotype == "WT"]["cell"].nunique()
    n_d_e = erica[erica.genotype == "Df16"]["cell"].nunique()

    # Plot 1: Allen IPFX only (clean replicate)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    plot_fi(a_stats, ax,
            "Figure 4D replicate — Allen IPFX (Protocol 2, 200 ms pulses)",
            n_wt_a, n_d_a)
    fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    plt.tight_layout()
    plt.savefig(OUT / "Figure4D_replicate_Allen.png",
                dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nSaved {OUT / 'Figure4D_replicate_Allen.png'}")

    # Plot 2: Allen vs Erica side by side
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.7), sharey=False)
    plot_fi(a_stats, axes[0],
            "Allen IPFX detection (Protocol 2)", n_wt_a, n_d_a)
    plot_fi(e_stats, axes[1],
            "Erica's CSV (manuscript Figure 4D source)", n_wt_e, n_d_e)
    fig.patch.set_alpha(0)
    for ax in axes:
        ax.patch.set_alpha(0)
    plt.tight_layout()
    plt.savefig(OUT / "Figure4D_replicate_compare.png",
                dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"Saved {OUT / 'Figure4D_replicate_compare.png'}")


if __name__ == "__main__":
    main()
