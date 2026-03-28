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
# # SCZ Protective Genes (OR < 1): Cell-Type Mutation Bias
#
# **R3 Q4 Reviewer Figure**: Cell-type mutation bias profile for schizophrenia genes
# with protective direction of effect (more mutations in controls than cases, OR < 1)
# from the SCHEMA exome sequencing study.
#
# These genes show a depleted CGE interneuron bias and positive bias toward non-neuronal
# cell types (Ependymal, Astrocyte, Vascular), the inverse of the SCZ risk gene pattern.

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

BIAS_FILE = BIAS_DIR / "HCT.SCZ500.protect.Z2.Spec.csv"
PVAL_FILE = RESULTS_DIR / "SCZ_protect_bias_addP.csv"

# %% [markdown]
# ## Load Data

# %%
bias_df = pd.read_csv(BIAS_FILE, index_col=0)
print(f"Bias data: {bias_df.shape}")
print(f"Columns: {list(bias_df.columns)}")

pval_df = pd.read_csv(PVAL_FILE, index_col=0)
print(f"\nBias+P data: {pval_df.shape}")
print(f"P-value columns: {[c for c in pval_df.columns if c in ['P-value','q-value','-logP','Z-score']]}")

# %% [markdown]
# ## Figure: Mutation Bias (EFFECT) and Significance (−log₁₀P) — 1×2 Panels

# %%
fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="none")

# Panel A: Bias (EFFECT)
axes[0].patch.set_alpha(0)
SuperClusterBias_BoxPlot(bias_df, "SCZ Protective (OR < 1)", ax=axes[0])

# Panel B: Significance (-logP)
axes[1].patch.set_alpha(0)
SuperClusterBias_BoxPlot(pval_df, "SCZ Protective (OR < 1)", EffectCol="-logP", fdr_cut=0.1, ax=axes[1])

fig.patch.set_alpha(0)
plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "FigSX_SCZ_protective_bias.pdf",
    dpi=300, bbox_inches="tight", transparent=True,
)
plt.show()
print(f"Saved: {FIGURES_DIR / 'FigSX_SCZ_protective_bias.pdf'}")

# %% [markdown]
# ## Summary Statistics

# %%
print("=" * 60)
print("SUMMARY: SCZ Protective Gene Bias")
print("=" * 60)

# Supercluster-level summary
sc_summary = bias_df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
print("\nTop 5 superclusters (highest bias):")
for sc, val in sc_summary.head(5).items():
    print(f"  {sc}: {val:.4f}")

print("\nBottom 5 superclusters (most depleted):")
for sc, val in sc_summary.tail(5).items():
    print(f"  {sc}: {val:.4f}")

cge = bias_df[bias_df["Supercluster"] == "CGE interneuron"]
mge = bias_df[bias_df["Supercluster"] == "MGE interneuron"]
print(f"\nCGE IN mean bias: {cge['EFFECT'].mean():.4f}")
print(f"MGE IN mean bias: {mge['EFFECT'].mean():.4f}")

# P-value summary
if "q-value" in pval_df.columns:
    n_sig = (pval_df["q-value"] <= 0.1).sum()
    print(f"\nSignificant cell types (q <= 0.1): {n_sig}/{len(pval_df)}")

    cge_p = pval_df[pval_df["Supercluster"] == "CGE interneuron"]
    if len(cge_p) > 0:
        print(f"CGE IN Z-scores: min={cge_p['Z-score'].min():.2f}, max={cge_p['Z-score'].max():.2f}, mean={cge_p['Z-score'].mean():.2f}")
        cge_sig = (cge_p["q-value"] <= 0.1).sum()
        print(f"CGE IN significant (q <= 0.1): {cge_sig}/{len(cge_p)}")

print("\n" + "=" * 60)
