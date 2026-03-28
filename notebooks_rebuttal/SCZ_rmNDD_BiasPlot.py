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
# ---

# %% [markdown]
# # SCZ Mutation Bias After Removing NDD Genes
#
# **R3 Q7 Reviewer Figure**: Cell-type mutation bias for schizophrenia genes
# after removing genes overlapping with NDD gene lists.
#
# - SCZ full: 53 genes (main pipeline, Mis2-excluded)
# - SCZ ∩ NDD top 61: **1 gene overlap (GRIN2A)** → 52 genes remain
# - SCZ ∩ NDD top 285 (bgmr): **8 genes overlap** → 45 genes remain
# - Spearman correlation (full vs rm NDD61): rho = 0.967
# - Spearman correlation (full vs rm NDD285): rho = 0.923

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import SuperClusterBias_BoxPlot

# %% [markdown]
# ## Configuration

# %%
BIAS_DIR = PROJ_DIR / "dat" / "Spec_Bias_Jul24_Centered"
FIGURES_DIR = PROJ_DIR / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Load Data

# %%
# SCZ full bias (recompute from main pipeline gene weights for consistency)
from CellType_PSY import Fil2Dict, HumanCT_AvgZ_Weighted, AnnotateCTDat

spec = pd.read_csv(PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv", index_col=0)
anno = pd.read_excel(PROJ_DIR / "dat" / "annotation.xlsx", index_col=0)
anno = anno[anno.index.notna()]
anno.index = anno.index.astype(int)

scz_gw = Fil2Dict(str(PROJ_DIR / "dat" / "GeneWeights" / "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw"))
scz_full_bias = HumanCT_AvgZ_Weighted(spec, scz_gw)
scz_full_bias = AnnotateCTDat(scz_full_bias, anno)
print(f"SCZ full: {len(scz_gw)} genes, {scz_full_bias.shape[0]} cell types")

# SCZ rm NDD61 (1 gene removed: GRIN2A)
scz_rm61_bias = pd.read_csv(BIAS_DIR / "HCT.SCZ61.ExNDD61.bgmr.csv", index_col=0)
print(f"SCZ rm NDD61: {scz_rm61_bias.shape[0]} cell types")

# SCZ rm NDD285 (8 genes removed)
scz_rm285_bias = pd.read_csv(BIAS_DIR / "HCT.SCZ61.ExNDD285.csv", index_col=0)
print(f"SCZ rm NDD285: {scz_rm285_bias.shape[0]} cell types")

# Spearman correlations
rho61, p61 = spearmanr(scz_full_bias["EFFECT"], scz_rm61_bias.loc[scz_full_bias.index, "EFFECT"])
rho285, p285 = spearmanr(scz_full_bias["EFFECT"], scz_rm285_bias.loc[scz_full_bias.index, "EFFECT"])
print(f"\nSpearman correlation (full vs rm NDD61):  rho={rho61:.4f}, p={p61:.2e}")
print(f"Spearman correlation (full vs rm NDD285): rho={rho285:.4f}, p={p285:.2e}")

# %% [markdown]
# ## Figure 1: Supercluster Bias — SCZ Full vs SCZ rm NDD61 vs SCZ rm NDD285

# %%
fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor="none")

for ax, (bias_df, title) in zip(axes, [
    (scz_full_bias, "SCZ (53 genes)"),
    (scz_rm61_bias, "SCZ rm NDD61 (52 genes)"),
    (scz_rm285_bias, "SCZ rm NDD285 (45 genes)"),
]):
    ax.patch.set_alpha(0)
    SuperClusterBias_BoxPlot(bias_df, title, ax=ax)

fig.patch.set_alpha(0)
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "FigSX_SCZ_rmNDD_bias.pdf",
    dpi=300, bbox_inches="tight", transparent=True,
)
plt.show()
print(f"Saved: {FIGURES_DIR / 'FigSX_SCZ_rmNDD_bias.pdf'}")

# %% [markdown]
# ## Figure 2: Combined — Bias Correlation + CGE Comparison

# %%
# Load ASD w/o ID bias from main results
MAIN_RESULTS = PROJ_DIR / "results" / "main_results" / "random" / "Centering"
asd_hiq_bias = pd.read_csv(MAIN_RESULTS / "ASD_HIQ_bias_addP.csv", index_col=0)

from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300, facecolor="none")

# --- Panel A: Bias correlation scatter ---
ax1.patch.set_alpha(0)

x = scz_full_bias["EFFECT"]
y = scz_rm285_bias.loc[scz_full_bias.index, "EFFECT"]

colors = []
for idx in scz_full_bias.index:
    sc = scz_full_bias.loc[idx, "Supercluster"]
    if sc == "CGE interneuron":
        colors.append("#E74C3C")
    elif sc == "MGE interneuron":
        colors.append("#3498DB")
    else:
        colors.append("#AAAAAA")

ax1.scatter(x, y, c=colors, s=15, alpha=0.7, edgecolor="none")

lims = [min(x.min(), y.min()), max(x.max(), y.max())]
ax1.plot(lims, lims, "k--", alpha=0.4, lw=1)

ax1.set_xlabel("SCZ Full (53 genes) — Mutation Bias", fontsize=12)
ax1.set_ylabel("SCZ rm NDD285 (45 genes) — Mutation Bias", fontsize=12)
ax1.set_title(f"Spearman rho = {rho285:.3f}", fontsize=13)

legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#E74C3C", markersize=8, label="CGE interneuron"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#3498DB", markersize=8, label="MGE interneuron"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#AAAAAA", markersize=8, label="Other"),
]
ax1.legend(handles=legend_elements, loc="upper left", framealpha=0)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.text(-0.12, 1.05, "A", transform=ax1.transAxes, fontsize=18, fontweight="bold", va="top")

# --- Panel B: CGE interneuron bias comparison ---
ax2.patch.set_alpha(0)

datasets = {
    "ASD w/o ID": asd_hiq_bias,
    "SCZ": scz_full_bias,
    "SCZ rm\nNDD61": scz_rm61_bias,
    "SCZ rm\nNDD285": scz_rm285_bias,
}

CT = "CGE interneuron"
CT_idx = anno[anno["Supercluster"] == CT].index.values

palette = ["#7FB3D8", "#E74C3C", "#F39C12", "#F5D76E"]
positions = range(1, len(datasets) + 1)

for i, (pos, (label, bias_df)) in enumerate(zip(positions, datasets.items())):
    vals = bias_df.loc[CT_idx, "EFFECT"]
    jitter_x = np.random.default_rng(42 + i).normal(pos, 0.04, size=len(vals))
    ax2.scatter(jitter_x, vals, color=palette[i], edgecolor="black", s=20, alpha=0.9, linewidth=0.3)

ax2.boxplot(
    [df.loc[CT_idx, "EFFECT"] for df in datasets.values()],
    positions=list(positions),
    showfliers=False,
    patch_artist=True,
    widths=0.4,
    boxprops=dict(facecolor="white", alpha=0.5, edgecolor="black", linewidth=1),
    medianprops=dict(color="grey", linewidth=1),
    whiskerprops=dict(color="grey", linewidth=1),
    capprops=dict(color="grey", linewidth=1),
)

test_pairs = [
    (1, 2, "ASD w/o ID", "SCZ"),
    (1, 3, "ASD w/o ID", "SCZ rm\nNDD61"),
    (1, 4, "ASD w/o ID", "SCZ rm\nNDD285"),
    (2, 3, "SCZ", "SCZ rm\nNDD61"),
    (2, 4, "SCZ", "SCZ rm\nNDD285"),
]

y_max = max(df.loc[CT_idx, "EFFECT"].max() for df in datasets.values())
y_step = y_max * 0.08

for idx, (x1, x2, lab1, lab2) in enumerate(test_pairs):
    vals1 = datasets[lab1].loc[CT_idx, "EFFECT"]
    vals2 = datasets[lab2].loc[CT_idx, "EFFECT"]
    _, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
    y = y_max + y_step * (idx + 1)
    ax2.plot([x1, x1, x2, x2], [y, y + y_step * 0.3, y + y_step * 0.3, y],
             lw=0.8, c="k", ls="--", alpha=0.7)
    pstr = f"P = {pval:.1e}" if pval < 0.01 else f"P = {pval:.2f}"
    ax2.text((x1 + x2) / 2, y + y_step * 0.35, pstr,
             ha="center", va="bottom", fontsize=9)

ax2.set_xticks(list(positions))
ax2.set_xticklabels(datasets.keys(), fontsize=11)
ax2.set_ylabel("Mutation Bias", fontsize=12)
ax2.set_title("CGE interneuron", fontsize=13)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.text(-0.12, 1.05, "B", transform=ax2.transAxes, fontsize=18, fontweight="bold", va="top")

fig.patch.set_alpha(0)
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "FigSX_SCZ_rmNDD_combined.pdf",
    dpi=300, bbox_inches="tight", transparent=True,
)
plt.show()
print(f"Saved: {FIGURES_DIR / 'FigSX_SCZ_rmNDD_combined.pdf'}")

# %% [markdown]
# ## Summary Statistics

# %%
print("=" * 60)
print("SUMMARY: SCZ rm NDD Bias")
print("=" * 60)

for name, df in [
    ("SCZ full (53)", scz_full_bias),
    ("SCZ rm NDD61 (52)", scz_rm61_bias),
    ("SCZ rm NDD285 (45)", scz_rm285_bias),
]:
    cge = df[df["Supercluster"] == "CGE interneuron"]
    mge = df[df["Supercluster"] == "MGE interneuron"]
    print(f"\n{name}:")
    print(f"  CGE IN mean bias: {cge['EFFECT'].mean():.4f}")
    print(f"  MGE IN mean bias: {mge['EFFECT'].mean():.4f}")
    sc_means = df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
    cge_rank = list(sc_means.index).index("CGE interneuron") + 1
    print(f"  CGE rank among {len(sc_means)} superclusters: {cge_rank}")

print(f"\nSpearman rho (full vs rm NDD61, 461 cell types):  {rho61:.4f}")
print(f"Spearman rho (full vs rm NDD285, 461 cell types): {rho285:.4f}")
print(f"Gene overlap: SCZ ∩ NDD61 = 1 (GRIN2A), SCZ ∩ NDD285 = 8")
print("\n" + "=" * 60)
