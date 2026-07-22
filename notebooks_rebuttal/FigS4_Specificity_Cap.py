# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: gencic
#     language: python
#     name: gencic
# ---

# %% [markdown]
# # Figure S4 — Specificity cap validation (canonical 2-panel generator).
# Panel A: specificity inflation vs library size. Panel B: mutation-bias
# robustness to cap level. Merged 2026-07-01 from Specificity_Cap_Analysis +
# Cap_Sensitivity_Figure (both now archived to dev_notebooks/).

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from matplotlib.ticker import FixedLocator, NullLocator, LogFormatterSciNotation
from scipy import stats
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import yaml
cfg = yaml.safe_load(open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml"))
PROJ_DIR = Path(cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import *

FIG_DIR = PROJ_DIR / "results" / "figures" / "specificity_cap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Font setup
font_path = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
if Path(font_path).exists():
    fm.fontManager.addfont(font_path)
    fm._load_fontmanager(try_read_cache=False)

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "font.size": 18,
    "font.family": "Arial",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 19.5,
    "xtick.labelsize": 16.5,
    "ytick.labelsize": 16.5,
})

# %% [markdown]
# ---
# ## Load shared data
#
# Unclipped specificity matrix (`clip100.0` is effectively uncapped; max ~97×).
# The cap threshold used throughout the study is `2 × mean(all values)`.

# %%
spec_unclip = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv",
    index_col=0,
)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]

mean_spec = np.mean(spec_unclip.values.flatten())
clip_threshold = mean_spec * 2
n_genes, n_cts = spec_unclip.shape
print(f"Matrix: {n_genes} genes × {n_cts} cell types")
print(f"Global mean specificity: {mean_spec:.4f}")
print(f"Clip threshold (2× mean): {clip_threshold:.4f}")
print(f"Max unclipped specificity: {spec_unclip.values.max():.1f}×")

# %% [markdown]
# ---
# ## Panel A — Fraction of genes exceeding cap vs library size
#
# Per cell type, the fraction of genes whose specificity exceeds the 2×-mean
# cap. Low-UMI (shallow) cell types carry more extreme-specificity genes, i.e.
# extreme specificity is largely a technical artifact of library depth.
# `Total_UMI`, `Supercluster`, and neuronal-vs-non-neuronal class are taken
# from the study's `Anno` table (`Total UMI` / `Supercluster` columns; neuronal
# membership from `Neur_idx`).

# %%
ct_stats = pd.DataFrame(index=Anno.index)
for ct in ct_stats.index:
    vals = spec_unclip[ct].values
    ct_stats.loc[ct, "frac_clipped"] = np.mean(vals > clip_threshold)
    ct_stats.loc[ct, "max_spec"] = np.max(vals)

ct_stats["Total_UMI"] = Anno["Total UMI"]
ct_stats["N_cells"] = Anno["Number of cells"]
ct_stats["Supercluster"] = Anno["Supercluster"]
ct_stats["is_neuronal"] = ct_stats.index.isin(Neur_idx)

for col in ["frac_clipped", "max_spec", "Total_UMI", "N_cells"]:
    ct_stats[col] = pd.to_numeric(ct_stats[col])

neur_mask = ct_stats["is_neuronal"]
rho_umi, p_umi = stats.spearmanr(ct_stats["Total_UMI"], ct_stats["frac_clipped"])
print(f"Spearman R (Total UMI, frac_clipped): rho = {rho_umi:.3f}, p = {p_umi:.2e}")
print(f"Neuronal: {neur_mask.sum()}, Non-neuronal: {(~neur_mask).sum()}")

# %% [markdown]
# ---
# ## Panel B — Cap-sweep robustness
#
# Recompute mutation bias for SCZ / ASD (HIQ) / DDD at a sweep of cap levels
# (clip the specificity matrix at `cap × mean`, then compute weighted bias with
# `HumanCT_AvgZ_Weighted`), and correlate each cap's cell-type bias vector
# against the reference cap (2×) via Spearman R.

# %%
# Gene weights for three disorders
gw_files = {
    "SCZ": PROJ_DIR / "dat" / "GeneWeights" / "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "ASD (HIQ)": PROJ_DIR / "dat" / "GeneWeights" / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "DDD": PROJ_DIR / "dat" / "GeneWeights" / "DDD.top61.gw.bgmr.csv",
}
gw_dicts = {}
for label, fpath in gw_files.items():
    gw_dicts[label] = Fil2Dict(str(fpath))
    print(f"{label}: {len(gw_dicts[label])} genes")

# %%
CAP_LEVELS = [1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 100.0]
ref_cap = 2.0

# bias_all[disorder][cap] = annotated bias DataFrame
bias_all = {}
for disorder, gw_dict in gw_dicts.items():
    bias_all[disorder] = {}
    for cap in CAP_LEVELS:
        threshold = mean_spec * cap
        # Clip specificity (no mean-centering — raw clipped values)
        spec_clipped = spec_unclip.clip(lower=0, upper=threshold)
        bias_df = HumanCT_AvgZ_Weighted(spec_clipped, gw_dict)
        bias_df = AnnotateCTDat(bias_df, Anno)
        bias_all[disorder][cap] = bias_df
    print(f"{disorder}: computed bias at {len(CAP_LEVELS)} cap levels")

# %%
# Spearman R vs cap = 2x. AnnotateCTDat sorts by EFFECT, so align by cell-type
# ID (the DataFrame index) before correlating.
corr_results = {}
for disorder in gw_dicts:
    ref_bias = bias_all[disorder][ref_cap]["EFFECT"]
    ref_bias.index = bias_all[disorder][ref_cap].index
    corr_results[disorder] = {}
    for cap in CAP_LEVELS:
        other_bias = bias_all[disorder][cap]["EFFECT"]
        other_bias.index = bias_all[disorder][cap].index
        common = ref_bias.index.intersection(other_bias.index)
        rho, _ = stats.spearmanr(ref_bias.loc[common], other_bias.loc[common])
        corr_results[disorder][cap] = rho

corr_df = pd.DataFrame(corr_results)
corr_df.index.name = "Cap level"
print("Spearman R vs cap = 2x:")
print(corr_df.round(4))

# %% [markdown]
# ---
# ## Assemble 2-panel Figure S4

# %%
disorder_colors = {"SCZ": "#2563eb", "ASD (HIQ)": "#dc2626", "DDD": "#16a34a"}
disorder_markers = {"SCZ": "o", "ASD (HIQ)": "s", "DDD": "D"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=150)

# ===== Panel A: fraction exceeding cap vs library size =====
axA = axes[0]
axA.scatter(ct_stats.loc[neur_mask, "Total_UMI"],
            ct_stats.loc[neur_mask, "frac_clipped"],
            color="red", alpha=0.5, s=22, edgecolors="white", lw=0.3,
            label="Neuronal", zorder=3)
axA.scatter(ct_stats.loc[~neur_mask, "Total_UMI"],
            ct_stats.loc[~neur_mask, "frac_clipped"],
            color="blue", alpha=0.6, s=22, edgecolors="white", lw=0.3,
            label="Non-neuronal", zorder=4)

axA.set_xscale("log")
axA.set_xlabel("Library size (total UMI)")
axA.set_ylabel("Fraction of genes exceeding cap")
# Show 10^3, 10^4, 10^5 on the log x-axis (extend limits if needed to include them)
_xlo, _xhi = axA.get_xlim()
axA.set_xlim(min(_xlo, 9e2), max(_xhi, 1.1e5))
axA.xaxis.set_major_locator(FixedLocator([1e3, 1e4, 1e5]))
axA.xaxis.set_minor_locator(NullLocator())
axA.xaxis.set_major_formatter(LogFormatterSciNotation())
axA.legend(fontsize=15, framealpha=0.8, loc="upper right")
axA.text(0.03, 0.03, f"$\\rho$ = {rho_umi:.3f}\np = {p_umi:.1e}",
         transform=axA.transAxes, ha="left", va="bottom", fontsize=16.5,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
axA.set_title("A", fontweight="bold", loc="left", fontsize=24, pad=10)

# ===== Panel B: Spearman R vs cap level =====
axB = axes[1]
# Shade the 1.5-3x range
axB.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)

for disorder in corr_df.columns:
    axB.plot(CAP_LEVELS, corr_df[disorder].values,
             color=disorder_colors[disorder],
             marker=disorder_markers[disorder],
             markersize=7, lw=2.2, label=disorder, zorder=3)

axB.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
axB.set_xscale("log")
axB.set_xticks(CAP_LEVELS)
axB.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=13.5,
                    rotation=45, ha="right", rotation_mode="anchor")
axB.minorticks_off()
axB.set_xlabel("Specificity cap level (× mean)")
axB.set_ylabel("Spearmans' R vs. cap = 2×")
axB.set_ylim(min(corr_df.values.min() - 0.03, 0.35), 1.01)
axB.legend(fontsize=15, frameon=True, framealpha=0.9, edgecolor="none",
           loc="lower left", handlelength=1.5)
axB.set_title("B", fontweight="bold", loc="left", fontsize=24, pad=10)

fig.tight_layout(w_pad=3.0)

fig.savefig(FIG_DIR / "FigS4_combined.pdf", dpi=300, transparent=True, bbox_inches="tight")
fig.savefig(FIG_DIR / "FigS4_combined.png", dpi=300, transparent=True, bbox_inches="tight")
print(f"Saved combined figure to {FIG_DIR / 'FigS4_combined.pdf'}")
plt.show()

# %% [markdown]
# ---
# ## Save each panel individually

# %%
# Panel A standalone
figA, axA = plt.subplots(figsize=(6.0, 5.0), dpi=150)
axA.scatter(ct_stats.loc[neur_mask, "Total_UMI"],
            ct_stats.loc[neur_mask, "frac_clipped"],
            color="red", alpha=0.5, s=22, edgecolors="white", lw=0.3,
            label="Neuronal", zorder=3)
axA.scatter(ct_stats.loc[~neur_mask, "Total_UMI"],
            ct_stats.loc[~neur_mask, "frac_clipped"],
            color="blue", alpha=0.6, s=22, edgecolors="white", lw=0.3,
            label="Non-neuronal", zorder=4)
axA.set_xscale("log")
axA.set_xlabel("Library size (total UMI)")
axA.set_ylabel("Fraction of genes exceeding cap")
# Show 10^3, 10^4, 10^5 on the log x-axis (extend limits if needed to include them)
_xlo, _xhi = axA.get_xlim()
axA.set_xlim(min(_xlo, 9e2), max(_xhi, 1.1e5))
axA.xaxis.set_major_locator(FixedLocator([1e3, 1e4, 1e5]))
axA.xaxis.set_minor_locator(NullLocator())
axA.xaxis.set_major_formatter(LogFormatterSciNotation())
axA.legend(fontsize=15, framealpha=0.8, loc="upper right")
axA.text(0.03, 0.03, f"$\\rho$ = {rho_umi:.3f}\np = {p_umi:.1e}",
         transform=axA.transAxes, ha="left", va="bottom", fontsize=16.5,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
figA.tight_layout()
figA.savefig(FIG_DIR / "FigS4_panelA.png", dpi=300, transparent=True, bbox_inches="tight")
print(f"Saved {FIG_DIR / 'FigS4_panelA.png'}")
plt.close(figA)

# %%
# Panel B standalone
figB, axB = plt.subplots(figsize=(6.0, 5.0), dpi=150)
axB.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)
for disorder in corr_df.columns:
    axB.plot(CAP_LEVELS, corr_df[disorder].values,
             color=disorder_colors[disorder],
             marker=disorder_markers[disorder],
             markersize=7, lw=2.2, label=disorder, zorder=3)
axB.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
axB.set_xscale("log")
axB.set_xticks(CAP_LEVELS)
axB.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=13.5,
                    rotation=45, ha="right", rotation_mode="anchor")
axB.minorticks_off()
axB.set_xlabel("Specificity cap level (× mean)")
axB.set_ylabel("Spearmans' R vs. cap = 2×")
axB.set_ylim(min(corr_df.values.min() - 0.03, 0.35), 1.01)
axB.legend(fontsize=15, frameon=True, framealpha=0.9, edgecolor="none",
           loc="lower left", handlelength=1.5)
figB.tight_layout()
figB.savefig(FIG_DIR / "FigS4_panelB.png", dpi=300, transparent=True, bbox_inches="tight")
print(f"Saved {FIG_DIR / 'FigS4_panelB.png'}")
plt.close(figB)

# %% [markdown]
# ---
# ## Summary statistics

# %%
print("=" * 60)
print("PANEL A — specificity inflation vs library size")
print("=" * 60)
total_vals = spec_unclip.values.size
total_clipped = int(np.sum(spec_unclip.values > clip_threshold))
print(f"Max unclipped specificity: {spec_unclip.values.max():.1f}×")
print(f"Clip threshold (2× mean): {clip_threshold:.4f}")
print(f"Values exceeding cap: {total_clipped:,} / {total_vals:,} "
      f"({100 * total_clipped / total_vals:.2f}%)")
print(f"Spearman R (library size, frac exceeding cap): "
      f"rho = {rho_umi:.3f}, p = {p_umi:.2e}")
print(f"Mean frac exceeding — Neuronal: {ct_stats.loc[neur_mask, 'frac_clipped'].mean():.2%}, "
      f"Non-neuronal: {ct_stats.loc[~neur_mask, 'frac_clipped'].mean():.2%}")

print("\n" + "=" * 60)
print("PANEL B — Spearman R vs cap = 2x")
print("=" * 60)
for disorder in corr_df.columns:
    rhos = [corr_results[disorder][c] for c in [1.5, 2.0, 2.5, 3.0]]
    print(f"{disorder}: rho at 1.5x-3x = {min(rhos):.4f}-{max(rhos):.4f} | "
          f"1.5x = {corr_results[disorder][1.5]:.4f} | "
          f"100x = {corr_results[disorder][100.0]:.4f}")
