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
# # Expression Matrix Preprocessing
#
# This notebook preprocesses raw cell type expression data into specificity matrices
# used by the CellTypeBias pipeline.
#
# ## Workflow Overview
#
# 1. **Load raw expression** - Gene × Cell Type matrix (TPM or UMI counts)
# 2. **TPM normalization** - Normalize within each cell type
# 3. **Specificity calculation** - Convert to relative specificity scores
# 4. **Filtering** - Remove lowly-expressed genes
# 5. **Clipping** - Cap extreme values
# 6. **Centering** - Mean-center per cell type (optional)
# 7. **Quantile normalization** - Alternative normalization (optional)
#
# ## Output Files
#
# All outputs go to `dat/ExpMats/`:
# - `HumanCT.TPM.{tpm_cut}.Filt.Spec.clip.lowexp.cut{exp_cut}.csv` - Raw specificity
# - `HumanCT.TPM.{tpm_cut}.Filt.Spec.clip.lowexp.cut{exp_cut}.mean_centered.csv` - Mean-centered (default)
# - `HumanCT.TPM.{tpm_cut}.Filt.Spec.clip.lowexp.cut{exp_cut}.qn.csv` - Quantile normalized

# %% [markdown]
# ## Setup

# %%
import sys
import os
from pathlib import Path
import yaml

# Project paths
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% [markdown]
# ## Configuration

# %%
# === INPUT FILES ===
# Raw expression matrix: genes (rows) x cell types (columns)
# Values should be raw counts or TPM
RAW_EXPRESSION_FILE = "/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv"

# Gene-level expression summary (for filtering low-expression genes)
GENE_EXPRESSION_FILE = "/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/HumanCT.MatchDF.csv"

# === OUTPUT DIRECTORY ===
OUTPUT_DIR = PROJ_DIR / "dat" / "ExpMats"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === PARAMETERS ===
TPM_CUT = 0.1          # Minimum TPM to keep a gene in a cell type
LOW_EXP_CUT = 10000    # Remove genes with total expression below this
CLIP_MULTIPLIER = 2    # Clip specificity at N × mean (override via CLI: python preprocessing.py --clip 3.0)

# Allow CLI override of CLIP_MULTIPLIER for sensitivity analysis
import argparse
if "__file__" in dir() or not hasattr(sys, "ps1"):  # running as script, not interactive
    _parser = argparse.ArgumentParser(add_help=False)
    _parser.add_argument("--clip", type=float, default=None)
    _args, _ = _parser.parse_known_args()
    if _args.clip is not None:
        CLIP_MULTIPLIER = _args.clip
        print(f"[CLI override] CLIP_MULTIPLIER = {CLIP_MULTIPLIER}")

print(f"Project directory: {PROJ_DIR}")
print(f"Output directory: {OUTPUT_DIR}")

# %% [markdown]
# ## Load Data

# %%
# Load raw expression matrix
print("Loading raw expression matrix...")
expr_raw = pd.read_csv(RAW_EXPRESSION_FILE, index_col=0)
print(f"  Shape: {expr_raw.shape[0]:,} genes × {expr_raw.shape[1]} cell types")

# Load gene expression levels (for filtering)
print("Loading gene expression summary...")
gene_exp = pd.read_csv(GENE_EXPRESSION_FILE, index_col=0)
print(f"  {len(gene_exp):,} genes with expression data")

# %% [markdown]
# ## Step 1: TPM Normalization and Specificity Calculation

# %%
def calculate_specificity_with_filtering(expression_df, min_tpm=0.1):
    """
    Calculate cell type specificity as fold-enrichment with TPM filtering.

    For each gene, specificity = TPM / mean(TPM across all cell types).
    This produces fold-enrichment values centered around 1.0 (not proportions).
    A gene exclusively expressed in 1 of N cell types has specificity = N in that type.

    Parameters
    ----------
    expression_df : DataFrame
        Expression matrix (genes × cell types)
    min_tpm : float
        Minimum TPM threshold for filtering

    Returns
    -------
    specificity_df : DataFrame
        Fold-enrichment specificity scores (genes × cell types)
    tpm_df : DataFrame
        TPM-normalized expression
    """
    # TPM normalization: scale each cell type to sum to 1M
    cell_type_sums = expression_df.sum(axis=0)
    tpm_df = expression_df.div(cell_type_sums, axis=1) * 1_000_000

    # Filter genes below TPM threshold (set to 0)
    tpm_filtered = tpm_df.where(tpm_df >= min_tpm, 0)

    # Specificity as fold-enrichment: divide by row mean (gene's mean TPM across cell types)
    gene_means = tpm_filtered.mean(axis=1)
    specificity_df = tpm_filtered.div(gene_means, axis=0)

    # Handle genes with zero mean (would cause NaN)
    specificity_df = specificity_df.fillna(0)

    return specificity_df, tpm_df

# %%
# Calculate specificity
print(f"Calculating specificity (TPM threshold = {TPM_CUT})...")
specificity_df, tpm_df = calculate_specificity_with_filtering(expr_raw, min_tpm=TPM_CUT)
print(f"  Specificity matrix shape: {specificity_df.shape}")

# %% [markdown]
# ## Step 2: Clip Extreme Values

# %%
# Clip at 2× mean specificity to reduce outlier influence
mean_spec = np.mean(specificity_df.values.flatten())
upper_clip = mean_spec * CLIP_MULTIPLIER
print(f"Mean specificity: {mean_spec:.6f}")
print(f"Clipping at: {upper_clip:.6f} ({CLIP_MULTIPLIER}× mean)")

specificity_clipped = specificity_df.clip(lower=0, upper=upper_clip)

# %% [markdown]
# ## Step 3: Filter Low-Expression Genes

# %%
# Identify genes below expression threshold
low_exp_genes = gene_exp[gene_exp["Exp"] < LOW_EXP_CUT].index.tolist()
print(f"Genes below expression threshold ({LOW_EXP_CUT:,}): {len(low_exp_genes):,}")

# Filter them out
specificity_filtered = specificity_clipped.loc[~specificity_clipped.index.isin(low_exp_genes)]
print(f"Final gene count: {len(specificity_filtered):,}")

# Ensure column names are integers (cell type IDs)
specificity_filtered.columns = specificity_filtered.columns.astype(int)

# %% [markdown]
# ## Step 4: Create Output Variants

# %%
def mean_center_by_celltype(spec_mat):
    """Mean-center specificity within each cell type."""
    return spec_mat.subtract(spec_mat.mean(axis=0), axis=1)

def quantile_normalize(df):
    """Quantile normalize across cell types."""
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

# %%
# Raw specificity (filtered + clipped)
spec_raw = specificity_filtered.copy()

# Mean-centered
spec_centered = mean_center_by_celltype(specificity_filtered)

# Quantile normalized
spec_qn = quantile_normalize(specificity_filtered)

print("Created three variants:")
print(f"  Raw:           {spec_raw.shape}")
print(f"  Mean-centered: {spec_centered.shape}")
print(f"  Quantile norm: {spec_qn.shape}")

# %% [markdown]
# ## Step 5: Save Outputs

# %%
# Build output filenames (include clip multiplier when non-default)
clip_tag = "clip" if CLIP_MULTIPLIER == 2 else f"clip{CLIP_MULTIPLIER}"
_exp = len(str(int(LOW_EXP_CUT))) - 1
_coeff = int(LOW_EXP_CUT) // (10 ** _exp)
base_name = f"HumanCT.TPM.{TPM_CUT}.Filt.Spec.{clip_tag}.lowexp.cut{_coeff}e{_exp}"

out_raw = OUTPUT_DIR / f"{base_name}.csv"
out_centered = OUTPUT_DIR / f"{base_name}.mean_centered.csv"
out_qn = OUTPUT_DIR / f"{base_name}.qn.csv"

print("Saving output files...")
spec_raw.to_csv(out_raw)
print(f"  {out_raw.name}")

spec_centered.to_csv(out_centered)
print(f"  {out_centered.name}")

spec_qn.to_csv(out_qn)
print(f"  {out_qn.name}")

print("\nDone!")

# %% [markdown]
# ## Diagnostic Plots

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Distribution of raw specificity
axes[0].hist(spec_raw.values.flatten(), bins=100, alpha=0.7, edgecolor='black')
axes[0].set_xlabel('Specificity')
axes[0].set_ylabel('Count')
axes[0].set_title('Raw Specificity Distribution')
axes[0].axvline(x=mean_spec, color='red', linestyle='--', label=f'Mean: {mean_spec:.4f}')
axes[0].legend()

# Distribution of mean-centered
axes[1].hist(spec_centered.values.flatten(), bins=100, alpha=0.7, edgecolor='black', color='green')
axes[1].set_xlabel('Specificity (centered)')
axes[1].set_ylabel('Count')
axes[1].set_title('Mean-Centered Distribution')
axes[1].axvline(x=0, color='red', linestyle='--')

# Expression vs fraction of zeros
common_genes = list(set(gene_exp.index) & set(spec_raw.index))
exp_vals = np.log10(gene_exp.loc[common_genes, "Exp"] + 1)
frac_zeros = (spec_raw.loc[common_genes] == 0).sum(axis=1) / spec_raw.shape[1]

axes[2].scatter(exp_vals, frac_zeros, s=1, alpha=0.3)
axes[2].axvline(x=np.log10(LOW_EXP_CUT + 1), color='red', linestyle='--', label=f'Cutoff: {LOW_EXP_CUT:,}')
axes[2].set_xlabel('Log10(Expression + 1)')
axes[2].set_ylabel('Fraction of Zeros')
axes[2].set_title('Expression vs Zero Fraction')
axes[2].legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "preprocessing_diagnostics.png", dpi=150)
plt.show()

# %% [markdown]
# ## Summary
#
# The preprocessing pipeline produces three expression specificity matrices:
#
# | File | Description | Use Case |
# |------|-------------|----------|
# | `*.csv` | Raw specificity (filtered, clipped) | Baseline |
# | `*.mean_centered.csv` | Mean-centered per cell type | **Default for bias analysis** |
# | `*.qn.csv` | Quantile normalized | Alternative normalization |
#
# The **mean-centered** version is recommended for bias analysis as it removes
# cell-type-specific baseline differences while preserving relative gene rankings.
