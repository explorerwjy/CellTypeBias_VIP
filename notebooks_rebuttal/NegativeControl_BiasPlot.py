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
from CellType_PSY import SuperClusterBias_BoxPlot

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
