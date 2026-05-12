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
#     name: python3
# ---

# %% [markdown]
# # Cap Sensitivity Figure
#
# Demonstrates that SCZ/ASD/DDD mutation bias is robust to the specificity
# capping threshold across a range of levels (1x, 1.5x, 2x, 2.5x, 3x, 5x).
#
# **Key conclusions:**
# - At 1x, signal is over-compressed (Spearmans' R drops).
# - At 1.5x--3x, results are highly stable (rho > 0.95).
# - The default 2x cap preserves biological signal while preventing
#   technical inflation from low-UMI cell types.

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
import seaborn as sns
from scipy import stats

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import HumanCT_AvgZ_Weighted, AnnotateCTDat, Fil2Dict, Anno

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
    "font.size": 11,
    "font.family": "Arial",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# %% [markdown]
# ## Load data

# %%
# Load unclipped specificity matrix (clip100 is effectively unclipped)
spec_unclip = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv",
    index_col=0,
)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]

mean_spec = np.mean(spec_unclip.values.flatten())
n_genes, n_cts = spec_unclip.shape
print(f"Matrix: {n_genes} genes x {n_cts} cell types")
print(f"Global mean specificity: {mean_spec:.4f}")

# %%
# Load gene weights for three disorders
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
# Load annotation
anno = pd.read_excel(PROJ_DIR / "dat" / "annotation.xlsx", index_col=0)
anno = anno[anno.index.notna()]
anno.index = anno.index.astype(int)

# %% [markdown]
# ## Compute bias at each cap level

# %%
CAP_LEVELS = [1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0]

# Store all results: bias_all[disorder][cap_multiplier] = annotated DataFrame
bias_all = {}

for disorder, gw_dict in gw_dicts.items():
    bias_all[disorder] = {}
    for cap in CAP_LEVELS:
        threshold = mean_spec * cap
        # Clip specificity (no mean-centering — raw clipped values)
        spec_clipped = spec_unclip.clip(lower=0, upper=threshold)
        # Compute weighted bias
        bias_df = HumanCT_AvgZ_Weighted(spec_clipped, gw_dict)
        bias_df = AnnotateCTDat(bias_df, Anno)
        bias_all[disorder][cap] = bias_df

    print(f"{disorder}: computed bias at {len(CAP_LEVELS)} cap levels")

# %% [markdown]
# ## Panel A: Supercluster-level heatmap (SCZ)

# %%
# Compute mean EFFECT per supercluster at each cap level for SCZ
disorder_for_heatmap = "SCZ"
sc_bias = {}

for cap in CAP_LEVELS:
    df = bias_all[disorder_for_heatmap][cap]
    sc_means = df.groupby("Supercluster")["EFFECT"].mean()
    sc_bias[cap] = sc_means

sc_bias_df = pd.DataFrame(sc_bias)
sc_bias_df.columns = [f"{c}x" for c in sc_bias_df.columns]

# Sort rows by bias at cap=2x (descending)
sort_col = "2.0x"
sc_bias_df = sc_bias_df.sort_values(sort_col, ascending=True)

print(f"Superclusters: {len(sc_bias_df)}")
print("\nTop 5 at cap=2x:")
print(sc_bias_df.sort_values(sort_col, ascending=False).head())

# %% [markdown]
# ## Panel B: Spearmans' R vs cap=2x

# %%
# For each disorder and cap level, compute Spearmans' R (cell-type-level)
# between EFFECT at that cap and EFFECT at cap=2x
ref_cap = 2.0
corr_results = {}

for disorder in gw_dicts.keys():
    ref_bias = bias_all[disorder][ref_cap].set_index(
        bias_all[disorder][ref_cap].index
    )["EFFECT"]
    corr_results[disorder] = {}
    for cap in CAP_LEVELS:
        other_bias = bias_all[disorder][cap]["EFFECT"]
        # Align by index (cell type IDs)
        other_bias.index = bias_all[disorder][cap].index
        common = ref_bias.index.intersection(other_bias.index)
        rho, pval = stats.spearmanr(ref_bias.loc[common], other_bias.loc[common])
        corr_results[disorder][cap] = rho

corr_df = pd.DataFrame(corr_results)
corr_df.index.name = "Cap level"
print(corr_df.round(4))

# %% [markdown]
# ## Panel C: Key supercluster bias across cap levels (SCZ)

# %%
# Extract mean EFFECT per supercluster across cap levels for key superclusters
key_superclusters = [
    "CGE interneuron",
    "MGE interneuron",
    "Deep-layer intratelencephalic",
    "Upper-layer intratelencephalic",
    "Astrocyte",
    "LAMP5-LHX6 and Chandelier",
]

sc_lines = {}
for sc in key_superclusters:
    sc_lines[sc] = []
    for cap in CAP_LEVELS:
        df = bias_all[disorder_for_heatmap][cap]
        vals = df.loc[df["Supercluster"] == sc, "EFFECT"]
        sc_lines[sc].append(vals.mean())

sc_lines_df = pd.DataFrame(sc_lines, index=CAP_LEVELS)
print(sc_lines_df.round(4))

# %% [markdown]
# ## Panel B: Disorder-specific supercluster bias across cap levels
#
# Show how the cap level affects discriminability between disorders
# for key cell types: MSN (ASD > SCZ) and CGE (SCZ > ASD).

# %%
# Compute mean supercluster bias for key superclusters across all cap levels
# Include neuronal and non-neuronal types for contrast
panel_scs = ["Medium spiny neuron", "CGE interneuron", "MGE interneuron",
             "Vascular", "Microglia", "Astrocyte", "Ependymal"]
comparison_disorders = ["SCZ", "ASD (HIQ)"]

sc_disorder_bias = {}
for sc in panel_scs:
    for disorder in comparison_disorders:
        key = f"{sc} | {disorder}"
        sc_disorder_bias[key] = []
        for cap in CAP_LEVELS:
            df = bias_all[disorder][cap]
            vals = df.loc[df["Supercluster"] == sc, "EFFECT"]
            sc_disorder_bias[key].append(vals.mean())

sc_disorder_df = pd.DataFrame(sc_disorder_bias, index=CAP_LEVELS)

# Compute mean supercluster RANK (across all superclusters) at each cap level
# Lower rank = higher bias; rank 1 = most biased supercluster
sc_rank_data = {}
for disorder in comparison_disorders:
    sc_rank_data[disorder] = {}
    for cap in CAP_LEVELS:
        df = bias_all[disorder][cap]
        sc_means = df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
        n_sc = len(sc_means)
        for rank_i, (sc_name, _) in enumerate(sc_means.items(), 1):
            if sc_name not in sc_rank_data[disorder]:
                sc_rank_data[disorder][sc_name] = {}
            sc_rank_data[disorder][sc_name][cap] = rank_i

print("Supercluster ranks at cap=2x vs 10x:")
for disorder in comparison_disorders:
    print(f"\n--- {disorder} ---")
    for sc in panel_scs:
        r2 = sc_rank_data[disorder][sc][2.0]
        r10 = sc_rank_data[disorder][sc][10.0]
        arrow = "↑" if r10 < r2 else ("↓" if r10 > r2 else "=")
        print(f"  {sc:30s}: rank {r2:2d} (2x) → {r10:2d} (10x)  {arrow}")

# %% [markdown]
# ## Assemble 3-panel figure

# %%
fig, axes = plt.subplots(1, 3, figsize=(21, 6), dpi=150)

# --- Shared style definitions ---
disorder_colors = {"SCZ": "#2563eb", "ASD (HIQ)": "#dc2626", "DDD": "#16a34a"}
disorder_markers = {"SCZ": "o", "ASD (HIQ)": "s", "DDD": "D"}

# Neuronal = solid, non-neuronal = dashed
sc_style = {
    "CGE interneuron":      {"color": "#e11d48", "marker": "o", "ls": "-",  "lw": 2.5},
    "MGE interneuron":      {"color": "#7c3aed", "marker": "s", "ls": "-",  "lw": 2.5},
    "Medium spiny neuron":  {"color": "#0891b2", "marker": "D", "ls": "-",  "lw": 2.5},
    "Vascular":             {"color": "#64748b", "marker": "v", "ls": "--", "lw": 1.8},
    "Microglia":            {"color": "#92400e", "marker": "^", "ls": "--", "lw": 1.8},
    "Astrocyte":            {"color": "#d97706", "marker": "P", "ls": "--", "lw": 1.8},
    "Ependymal":            {"color": "#6b7280", "marker": "X", "ls": "--", "lw": 1.8},
}
sc_short = {
    "CGE interneuron": "CGE IN", "MGE interneuron": "MGE IN",
    "Medium spiny neuron": "MSN", "Vascular": "Vascular",
    "Microglia": "Microglia", "Astrocyte": "Astrocyte", "Ependymal": "Ependymal",
}
show_disorder = "SCZ"

# Shared x-axis formatting helper
def fmt_xaxis(ax):
    ax.set_xticks(CAP_LEVELS)
    ax.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=9)
    ax.set_xlabel("Specificity cap level (× mean)", fontsize=11)
    ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)

# Optimal range shading (1.5–3×) on all panels
def shade_optimal(ax):
    ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)

# ===== Panel A: Spearmans' R =====
ax = axes[0]
shade_optimal(ax)

for disorder in corr_df.columns:
    rho_vals = corr_df[disorder].values
    ax.plot(CAP_LEVELS, rho_vals,
            color=disorder_colors[disorder],
            marker=disorder_markers[disorder],
            markersize=7, lw=2.2, label=disorder, zorder=3)

fmt_xaxis(ax)
ax.set_ylabel("Spearmans' R vs. cap = 2×", fontsize=11)
ax.set_ylim(min(corr_df.values.min() - 0.03, 0.35), 1.01)
ax.legend(fontsize=9, frameon=True, framealpha=0.9, edgecolor="none",
          loc="lower left", handlelength=1.5)
ax.set_title("A", fontweight="bold", loc="left", fontsize=16, pad=10)

# ===== Panel B: Raw bias with neuronal + non-neuronal lines (SCZ) =====
ax = axes[1]
shade_optimal(ax)

# Plot neuronal first (solid), then non-neuronal (dashed) for legend grouping
neuronal_scs = [sc for sc in panel_scs if sc_style[sc]["ls"] == "-"]
nonneuronal_scs = [sc for sc in panel_scs if sc_style[sc]["ls"] == "--"]

for sc in neuronal_scs + nonneuronal_scs:
    key = f"{sc} | {show_disorder}"
    sty = sc_style[sc]
    ax.plot(CAP_LEVELS, sc_disorder_df[key].values,
            color=sty["color"], marker=sty["marker"], markersize=5,
            lw=sty["lw"], ls=sty["ls"], label=sc_short[sc], zorder=3)

fmt_xaxis(ax)
ax.set_ylabel("Mean supercluster bias", fontsize=11)
ax.set_title("B", fontweight="bold", loc="left", fontsize=16, pad=10)

# Two-column legend: neuronal | non-neuronal
leg = ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, edgecolor="none",
                loc="lower right", ncol=2, columnspacing=1.0, handlelength=1.8,
                title="Neuronal          Non-neuronal", title_fontsize=7.5)
leg._legend_box.align = "left"

# ===== Panel C: Supercluster rank across cap levels (SCZ) =====
ax = axes[2]
shade_optimal(ax)

for sc in neuronal_scs + nonneuronal_scs:
    sty = sc_style[sc]
    ranks = [sc_rank_data[show_disorder][sc][cap] for cap in CAP_LEVELS]
    ax.plot(CAP_LEVELS, ranks,
            color=sty["color"], marker=sty["marker"], markersize=5,
            lw=sty["lw"], ls=sty["ls"], label=sc_short[sc], zorder=3)

fmt_xaxis(ax)
ax.invert_yaxis()
ax.set_ylabel("Supercluster rank (1 = highest bias)", fontsize=11)
ax.set_title("C", fontweight="bold", loc="left", fontsize=16, pad=10)

leg = ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, edgecolor="none",
                loc="lower right", ncol=2, columnspacing=1.0, handlelength=1.8,
                title="Neuronal          Non-neuronal", title_fontsize=7.5)
leg._legend_box.align = "left"

fig.tight_layout(w_pad=3.5)

# ===== Save =====
fig.savefig(
    FIG_DIR / "cap_sensitivity_figure.pdf",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)
fig.savefig(
    FIG_DIR / "cap_sensitivity_figure.png",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)
plt.show()
print(f"\nSaved to {FIG_DIR / 'cap_sensitivity_figure.pdf'}")

# %% [markdown]
# ## Summary statistics for manuscript text

# %%
print("=" * 60)
print("SPEARMAN CORRELATION TABLE (vs. cap = 2x)")
print("=" * 60)
for disorder in corr_df.columns:
    print(f"\n{disorder}:")
    for cap in CAP_LEVELS:
        rho = corr_results[disorder][cap]
        marker = " <--" if cap == ref_cap else ""
        print(f"  Cap {cap}x: rho = {rho:.4f}{marker}")

print("\n" + "=" * 60)
print("KEY FINDING: Stability range")
print("=" * 60)
for disorder in corr_df.columns:
    rhos = [corr_results[disorder][c] for c in [1.5, 2.0, 2.5, 3.0]]
    print(f"{disorder} at 1.5x-3x: rho = {min(rhos):.4f} - {max(rhos):.4f}")
    print(f"  at 1x: rho = {corr_results[disorder][1.0]:.4f}")
    print(f"  at 5x: rho = {corr_results[disorder][5.0]:.4f}")
