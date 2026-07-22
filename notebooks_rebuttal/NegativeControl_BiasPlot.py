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
# # Negative Control Traits: Cell-Type Mutation Bias
#
# **R3 Q3 Reviewer Figure**: Cell-type mutation bias profiles for non-brain negative control
# traits (HDL cholesterol, IBD, Alanine aminotransferase) derived from UKBB rare-variant
# burden analyses.
#
# These traits show glia-dominated bias profiles with no CGE interneuron enrichment,
# demonstrating that the CGE signal is specific to psychiatric/cognitive disorders.

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

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import SuperClusterBias_BoxPlot, plot_mutation_bias_comparison_V2, Anno

# %% [markdown]
# ## Configuration

# %%
BIAS_DIR = PROJ_DIR / "dat" / "Spec_Bias_Jul24_Centered"
RESULTS_DIR = PROJ_DIR / "results" / "systematic_best_of_n_WB_mean_phastCons_n_CDS_bases_Best1000" / "Centering"
FIGURES_DIR = PROJ_DIR / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Three negative controls: glia cells on top
TRAITS = {
    "HDL Cholesterol": {
        "bias": BIAS_DIR / "HCT.HDL_C.top61.csv",
        "bias_addP": RESULTS_DIR / "NegCtrl_HDL_bias_addP.csv",
    },
    "IBD": {
        "bias": BIAS_DIR / "HCT.IBD.top61.csv",
        "bias_addP": RESULTS_DIR / "NegCtrl_IBD_bias_addP.csv",
    },
    "Alanine AT": {
        "bias": BIAS_DIR / "HCT.Alanine.top61.csv",
        "bias_addP": RESULTS_DIR / "NegCtrl_Alanine_bias_addP.csv",
    },
}

trait_names = list(TRAITS.keys())

# %% [markdown]
# ## Load Data

# %%
# Load bias (EFFECT) results
bias_data = {}
for name, paths in TRAITS.items():
    if paths["bias"].exists():
        bias_data[name] = pd.read_csv(paths["bias"], index_col=0)
        print(f"Loaded bias for {name}: {bias_data[name].shape}")
    else:
        print(f"MISSING bias: {name}")

# Load bias + P-value results
pval_data = {}
for name, paths in TRAITS.items():
    if paths["bias_addP"].exists():
        pval_data[name] = pd.read_csv(paths["bias_addP"], index_col=0)
        print(f"Loaded bias+P for {name}: {pval_data[name].shape}, cols: {[c for c in pval_data[name].columns if c in ['P-value','q-value','-logP','Z-score']]}")
    else:
        print(f"MISSING bias+P: {name} — will skip p-value panel")

# %% [markdown]
# ## Figure 1: Mutation Bias (EFFECT) — 1×3 Panels

# %%
fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor="none")

for ax, name in zip(axes, trait_names):
    ax.patch.set_alpha(0)
    if name in bias_data:
        SuperClusterBias_BoxPlot(bias_data[name], name, ax=ax)
    else:
        ax.text(0.5, 0.5, f"No data for {name}", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(name, fontsize=14, fontweight="bold")

fig.suptitle("Negative Control Traits — Cell-Type Mutation Bias", fontsize=16, fontweight="bold", y=1.02)
fig.patch.set_alpha(0)
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "FigSX_negative_controls_bias.pdf",
    dpi=300, bbox_inches="tight", transparent=True,
)
plt.show()
print(f"Saved: {FIGURES_DIR / 'FigSX_negative_controls_bias.pdf'}")

# %% [markdown]
# ## Figure 2: Significance (−log₁₀P) — 1×3 Panels

# %%
if len(pval_data) == 3:
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor="none")

    for ax, name in zip(axes, trait_names):
        ax.patch.set_alpha(0)
        if name in pval_data:
            SuperClusterBias_BoxPlot(pval_data[name], name, EffectCol="-logP", fdr_cut=0.1, ax=ax)
        else:
            ax.text(0.5, 0.5, f"No data for {name}", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(name, fontsize=14, fontweight="bold")

    fig.suptitle("Negative Control Traits — Significance (−log₁₀ P-value)", fontsize=16, fontweight="bold", y=1.02)
    fig.patch.set_alpha(0)
    plt.tight_layout()
    plt.savefig(
        FIGURES_DIR / "FigSX_negative_controls_pval.pdf",
        dpi=300, bbox_inches="tight", transparent=True,
    )
    plt.show()
    print(f"Saved: {FIGURES_DIR / 'FigSX_negative_controls_pval.pdf'}")
else:
    print(f"Only {len(pval_data)}/3 traits have p-values. Skipping p-value figure.")
    print("Missing:", [n for n in trait_names if n not in pval_data])

# %% [markdown]
# ## Summary Statistics

# %%
print("=" * 60)
print("SUMMARY: Negative Control Bias Statistics")
print("=" * 60)

for name in trait_names:
    if name not in bias_data:
        continue
    df = bias_data[name]
    cge = df[df["Supercluster"] == "CGE interneuron"]
    mge = df[df["Supercluster"] == "MGE interneuron"]
    glia_scs = ["Microglia", "Astrocyte", "Oligodendrocyte", "OPC"]
    glia = df[df["Supercluster"].isin(glia_scs)]

    print(f"\n{name}:")
    print(f"  Top 3 superclusters: {df.groupby('Supercluster')['EFFECT'].mean().sort_values(ascending=False).head(3).to_dict()}")
    print(f"  CGE IN mean bias: {cge['EFFECT'].mean():.4f}")
    print(f"  MGE IN mean bias: {mge['EFFECT'].mean():.4f}")
    print(f"  Glia mean bias:   {glia['EFFECT'].mean():.4f}")

    if name in pval_data:
        df_p = pval_data[name]
        n_sig = (df_p["q-value"] <= 0.1).sum() if "q-value" in df_p.columns else "N/A"
        print(f"  Significant cell types (q ≤ 0.1): {n_sig}")
        cge_p = df_p[df_p["Supercluster"] == "CGE interneuron"]
        if len(cge_p) > 0 and "q-value" in cge_p.columns:
            cge_sig = (cge_p["q-value"] <= 0.1).sum()
            print(f"  CGE IN significant (q ≤ 0.1): {cge_sig}/{len(cge_p)}")

print("\n" + "=" * 60)

# %% [markdown]
# ## Figure 3 (Reviewer Response): CGE Interneuron — VNR vs Non-Brain Negative Controls
#
# Fig 3A-style strip plot showing per-cell-type CGE-interneuron mutation bias for
# VNR−, VNR+, and the four non-brain negative-control traits shown in the QQ-plot
# panel (HDL, Alanine, RBC, IBD).
#
# **Pairwise tests:** Mann–Whitney U comparing VNR− to each other cohort across the
# 31 CGE-interneuron cell types, BH-FDR corrected within the 5 pairs.

# %%
import scipy.stats as _sstats
from statsmodels.stats.multitest import multipletests as _multipletests

BIAS_ADDP_DIR = PROJ_DIR / "results" / "main_results" / "random" / "Centering"

cge_datasets = {
    "VNR-":     pd.read_csv(BIAS_ADDP_DIR / "UKBB_VNR_Neg_bias_addP.csv", index_col=0),
    "VNR+":     pd.read_csv(BIAS_ADDP_DIR / "UKBB_VNR_Pos_bias_addP.csv", index_col=0),
    "HDL-C":    pd.read_csv(BIAS_ADDP_DIR / "NegCtrl_HDL_bias_addP.csv", index_col=0),
    "Alanine":  pd.read_csv(BIAS_ADDP_DIR / "NegCtrl_Alanine_bias_addP.csv", index_col=0),
    "IBD":      pd.read_csv(BIAS_ADDP_DIR / "NegCtrl_IBD_bias_addP.csv", index_col=0),
}

CT = "CGE interneuron"
# RBC dropped for consistency with the Fig S7 negative-control set.
test_pairs = [("VNR-", "VNR+"), ("VNR-", "HDL-C"), ("VNR-", "Alanine"),
              ("VNR-", "IBD")]

# Compute on-the-fly Mann-Whitney + BH-FDR over the 5 specified pairs
_CT_idx = Anno[Anno["Supercluster"] == CT].index.values
_pvals_raw = []
for a, b in test_pairs:
    _, p = _sstats.mannwhitneyu(cge_datasets[a].loc[_CT_idx, "EFFECT"],
                                cge_datasets[b].loc[_CT_idx, "EFFECT"])
    _pvals_raw.append(p)
_pvals_fdr = _multipletests(_pvals_raw, method="fdr_bh")[1]

cge_PvalDF = pd.DataFrame([
    {"Pair": f"{a} - {b}", "SuperCluster": CT, "MWU_FDR": p}
    for (a, b), p in zip(test_pairs, _pvals_fdr)
])

print("Pairwise BH-FDR p-values (n=5 tests):")
for (a, b), p_raw, p_fdr in zip(test_pairs, _pvals_raw, _pvals_fdr):
    print(f"  {a} vs {b}: raw p={p_raw:.3e}, BH-FDR={p_fdr:.3e}")

# %%
# Inline Fig 3A-style plot (custom bracket spacing for 5-pair comparison)
def _format_pval(p):
    if p < 1e-3:
        return f"p={p:.1e}"
    return f"p={p:.3f}"

# Build per-cell-type values and sort cohorts by mean
data = {name: df.loc[_CT_idx, "EFFECT"] for name, df in cge_datasets.items()}
sorted_keys = sorted(data.keys(), key=lambda k: data[k].mean())
sorted_data = {k: data[k] for k in sorted_keys}

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 13

# Taller figure to accommodate 5 stacked brackets
fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=300, facecolor="none")
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

positions = range(1, len(sorted_data) + 1)
colors = plt.cm.Set2(np.linspace(0, 1, len(sorted_data)))

for i, (pos, (label, values)) in enumerate(zip(positions, sorted_data.items())):
    x = np.random.normal(pos, 0.04, size=len(values))
    ax.scatter(x, values, color=colors[i], edgecolor="black", s=22, alpha=1.0, label=label)

bp = ax.boxplot(
    [v.values for v in sorted_data.values()],
    positions=list(positions), showfliers=False, patch_artist=True, widths=0.4,
    boxprops=dict(facecolor="white", alpha=0, edgecolor="black", linewidth=1),
    medianprops=dict(color="grey", linewidth=1),
    whiskerprops=dict(color="grey", linewidth=1),
    capprops=dict(color="grey", linewidth=1),
)
for box in bp["boxes"]:
    box.set(facecolor="white", alpha=0.6, edgecolor="black", linewidth=1)

# Bracket placement — generous spacing for 5 stacked brackets
global_y_max = max(v.max() for v in sorted_data.values())
global_y_min = min(v.min() for v in sorted_data.values())
y_range = global_y_max - global_y_min
y_offset = y_range * 0.06       # bracket arm height
min_sep   = y_range * 0.11      # vertical gap between consecutive brackets

annotation_heights = []
fdr_lookup = dict(zip([(a, b) for a, b in test_pairs], _pvals_fdr))

# Sort pairs by x-distance (shortest first → bottom, longest → top) to reduce visual crossing
pair_with_dist = []
labels = list(sorted_data.keys())
for a, b in test_pairs:
    x1 = labels.index(a) + 1
    x2 = labels.index(b) + 1
    pair_with_dist.append(((a, b), x1, x2, abs(x2 - x1)))
pair_with_dist.sort(key=lambda t: t[3])

for (a, b), x1, x2, _ in pair_with_dist:
    p_fdr = fdr_lookup[(a, b)]
    x_center = (x1 + x2) / 2
    y_base = max(data[a].max(), data[b].max())
    if annotation_heights:
        y = max(max(annotation_heights) + min_sep, y_base + y_offset)
    else:
        y = y_base + y_offset
    annotation_heights.append(y + y_offset)
    ax.plot([x1, x1, x2, x2], [y, y + y_offset / 2, y + y_offset / 2, y],
            lw=0.9, c="k", ls="--", alpha=0.75)
    ax.text(x_center, y + y_offset / 2 + 0.01 * y_range, _format_pval(p_fdr),
            ha="center", va="bottom", fontsize=10)

# Y-limits — leave headroom above the topmost bracket
ax.set_ylim(top=max(annotation_heights) + y_offset * 2)

ax.set_xticks(list(positions))
ax.set_xticklabels(list(sorted_data.keys()), rotation=45, ha="right",
                   rotation_mode="anchor", fontsize=13)
ax.set_ylabel("Mutation Bias", fontsize=14)
ax.set_title(f"{CT}", fontsize=14, pad=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("black")
    ax.spines[s].set_linewidth(1)
ax.grid(linestyle="--", alpha=0.3)
plt.tight_layout()

# Save BEFORE show (avoid inline backend clearing the figure)
fig.savefig(FIGURES_DIR / "FigR_CGE_VNR_vs_NegCtrls.pdf",
            dpi=300, bbox_inches="tight", transparent=True)
fig.savefig(FIGURES_DIR / "FigR_CGE_VNR_vs_NegCtrls.png",
            dpi=300, bbox_inches="tight", transparent=True)
print(f"Saved: {FIGURES_DIR / 'FigR_CGE_VNR_vs_NegCtrls.pdf'}")
plt.show()

# %%
