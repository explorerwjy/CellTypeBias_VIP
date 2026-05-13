# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (gencic)
#     language: python
#     name: gencic
# ---

# %% [markdown]
# # Module E: Multi-Modal CCKBC Classification and Updated 22q Bias
#
# This notebook implements the Module E pipeline for the CGE subtype analysis.
#
# ## Overview
#
# Previous 22q bias analyses used marker-based classification alone (CCK+, SNCG+, M2R-, CR-).
# Module E integrates evidence from four independent modules:
#
# - **Module A** (mouse homology): Allen Mouse CCF Sncg subclass correspondence
# - **Module C** (transcriptomic mapping): Direct/indirect scVI integration paths
# - **Module D** (electrophysiology): Fast-spiking ephys profile via scANVI mapping
# - **Existing markers**: CCK and SNCG expression thresholds
#
# Each CGE cluster receives a `cckbc_confidence` score (0–6) and a tier:
#
# | Score | Tier |
# |-------|------|
# | ≥ 3   | high-confidence CCKBC |
# | 1–2   | tentative CCKBC |
# | 0     | ISI VIP |
#
# The updated 22q bias comparison tests whether CCKBC clusters show systematically
# lower 22q deletion bias than ISI VIP clusters across four 22q gene sets.

# %% [markdown]
# ## 1. Setup

# %%
# %load_ext autoreload
# %autoreload 2

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Project paths (relative to notebook location: cge_subtype/notebooks/)
_NOTEBOOK_DIR = Path(os.path.abspath(""))
_PROJ_DIR = _NOTEBOOK_DIR.parent.parent  # CellTypeBias_VIP/

# Add project src to path for CellType_PSY
sys.path.insert(0, str(_PROJ_DIR / "src"))
from CellType_PSY import HumanCT_AvgZ_Weighted, AnnotateCTDat

# Add project root for cge_subtype package
sys.path.insert(0, str(_PROJ_DIR))
from cge_subtype.src.multimodal_classification import classify_clusters, EVIDENCE_COLUMNS

# Plotting style
sns.set_style("ticks")
plt.rcParams["figure.dpi"] = 120

print(f"Project root: {_PROJ_DIR}")

# %% [markdown]
# ## 2. Load Expression Specificity Matrix

# %%
SPEC_MAT_PATH = _PROJ_DIR / "dat" / "ExpMats" / "HumanCT.spec.csv"

# TODO: update path if a mean-centered or quantile-normalised variant is preferred
# e.g. HumanCT.spec.mean_centered.csv
print(f"Loading specificity matrix from:\n  {SPEC_MAT_PATH}")

SpecMat = pd.read_csv(SPEC_MAT_PATH, index_col=0)
print(f"Specificity matrix: {SpecMat.shape[0]} genes x {SpecMat.shape[1]} cell types")

# %% [markdown]
# ## 3. Load 22q Gene Sets

# %%
GW_DIR = _PROJ_DIR / "dat" / "GeneWeights"

# Four 22q gene sets
# TODO: confirm that X22q_HighExp.gw.csv and X22q_DEG_d75.gw.csv exist, or update paths
GENE_SET_FILES = {
    "22q_del":      GW_DIR / "X22q.gw.csv",
    "22q_mouse":    GW_DIR / "X22q.mousemodel.gw.csv",
    "22q_HighExp":  GW_DIR / "X22q_HighExp.gw.csv",    # TODO: verify filename
    "22q_DEG_d75":  GW_DIR / "X22q_DEG_d75.gw.csv",    # TODO: verify filename
}

gene_sets = {}
for name, path in GENE_SET_FILES.items():
    if path.exists():
        gw = pd.read_csv(path, header=None, names=["Entrez_ID", "Weight"])
        gw["Entrez_ID"] = gw["Entrez_ID"].astype(str)
        gene_sets[name] = gw
        print(f"  {name}: {len(gw)} genes loaded from {path.name}")
    else:
        print(f"  WARNING: {name} not found at {path}")

print(f"\nLoaded {len(gene_sets)} of {len(GENE_SET_FILES)} gene sets.")

# %% [markdown]
# ## 4. Assemble Multi-Modal Evidence Table
#
# Each row is one CGE cluster (index = cluster integer ID).
# Boolean columns represent evidence from each module.
#
# **TODO**: populate each evidence column from the relevant module outputs:
#
# - `module_a_is_sncg`        ← Module A output CSV (e.g. `cge_subtype/results/module_a_sncg_flag.csv`)
# - `module_c_direct_cckbc`   ← Module C output (direct path, `cge_subtype/results/module_c_direct.csv`)
# - `module_c_indirect_cckbc` ← Module C output (indirect path)
# - `module_d_fast_spiking`   ← Module D output (e.g. `cge_subtype/results/module_d_ephys_flag.csv`)
# - `marker_cck_positive`     ← From existing CGE_VIP_annotation*.csv EFFECT columns or expression
# - `marker_sncg_positive`    ← Same source

# %%
# TODO: replace this placeholder with actual module outputs

# Expected column structure shown here; all values set to False until modules provide data
_N_CLUSTERS = 21  # total human CGE clusters

# Placeholder: replace with real cluster IDs from annotation
# TODO: load from results/main_results/contrasts/CGE_VIP_annotation.csv
CGE_CLUSTER_IDS = list(range(1, _N_CLUSTERS + 1))  # placeholder IDs

evidence_df = pd.DataFrame(
    {col: False for col in EVIDENCE_COLUMNS},
    index=CGE_CLUSTER_IDS,
)
evidence_df.index.name = "cluster_id"
evidence_df = evidence_df.reset_index()

print("Evidence table (placeholder — all False):")
print(evidence_df.head())
print(f"\nExpected evidence columns: {EVIDENCE_COLUMNS}")

# %% [markdown]
# ### 4a. Load Module A Evidence (Mouse Sncg Flag)

# %%
# TODO: load Module A Sncg correspondence flag
# MODULE_A_PATH = _PROJ_DIR / "cge_subtype" / "results" / "module_a_sncg_flag.csv"
# if MODULE_A_PATH.exists():
#     module_a = pd.read_csv(MODULE_A_PATH)
#     # Merge into evidence_df on cluster_id
#     evidence_df = evidence_df.merge(
#         module_a[["cluster_id", "module_a_is_sncg"]], on="cluster_id", how="left", suffixes=("_old", "")
#     )
#     evidence_df["module_a_is_sncg"] = evidence_df["module_a_is_sncg"].fillna(False).astype(bool)
#     print(f"Module A: {evidence_df['module_a_is_sncg'].sum()} clusters flagged as Sncg")
# else:
#     print(f"Module A output not found: {MODULE_A_PATH}")

# %% [markdown]
# ### 4b. Load Module C Evidence (scVI Mapping Paths)

# %%
# TODO: load Module C direct and indirect CCKBC mapping flags
# MODULE_C_PATH = _PROJ_DIR / "cge_subtype" / "results" / "module_c_mapping.csv"
# if MODULE_C_PATH.exists():
#     module_c = pd.read_csv(MODULE_C_PATH)
#     for col in ["module_c_direct_cckbc", "module_c_indirect_cckbc"]:
#         evidence_df[col] = evidence_df["cluster_id"].map(
#             module_c.set_index("cluster_id")[col]
#         ).fillna(False).astype(bool)
#     print(f"Module C direct:   {evidence_df['module_c_direct_cckbc'].sum()} clusters")
#     print(f"Module C indirect: {evidence_df['module_c_indirect_cckbc'].sum()} clusters")

# %% [markdown]
# ### 4c. Load Module D Evidence (Ephys Fast-Spiking Flag)

# %%
# TODO: load Module D fast-spiking classification flag
# MODULE_D_PATH = _PROJ_DIR / "cge_subtype" / "results" / "module_d_ephys_flag.csv"
# if MODULE_D_PATH.exists():
#     module_d = pd.read_csv(MODULE_D_PATH)
#     evidence_df["module_d_fast_spiking"] = evidence_df["cluster_id"].map(
#         module_d.set_index("cluster_id")["module_d_fast_spiking"]
#     ).fillna(False).astype(bool)
#     print(f"Module D: {evidence_df['module_d_fast_spiking'].sum()} clusters flagged as fast-spiking")

# %% [markdown]
# ### 4d. Load Existing Marker Evidence (CCK+, SNCG+)

# %%
# TODO: load CCK and SNCG marker expression from CGE_VIP_annotation.csv
# ANNOTATION_PATH = _PROJ_DIR / "results" / "main_results" / "contrasts" / "CGE_VIP_annotation.csv"
# if ANNOTATION_PATH.exists():
#     anno = pd.read_csv(ANNOTATION_PATH, index_col="Index")
#     # Apply cutoffs: CCK > 2.0, SNCG > 0.5 (log2 UMI+1 scale)
#     CCK_CUTOFF = 2.0
#     SNCG_CUTOFF = 0.5
#     evidence_df["marker_cck_positive"] = evidence_df["cluster_id"].map(
#         (anno["CCK_ExpL"] > CCK_CUTOFF).astype(bool)
#     ).fillna(False)
#     evidence_df["marker_sncg_positive"] = evidence_df["cluster_id"].map(
#         (anno["SNCG_ExpL"] > SNCG_CUTOFF).astype(bool)
#     ).fillna(False)
#     print(f"CCK+:  {evidence_df['marker_cck_positive'].sum()} clusters")
#     print(f"SNCG+: {evidence_df['marker_sncg_positive'].sum()} clusters")

# %% [markdown]
# ## 5. Classify Clusters

# %%
classified = classify_clusters(evidence_df)

print("Classification results:")
print(classified[["cluster_id", "cckbc_confidence", "cckbc_tier"]].to_string(index=False))

tier_counts = classified["cckbc_tier"].value_counts()
print(f"\nTier distribution:\n{tier_counts.to_string()}")

# Convenience subsets
high_conf = classified[classified["cckbc_tier"] == "high-confidence CCKBC"]["cluster_id"].tolist()
tentative = classified[classified["cckbc_tier"] == "tentative CCKBC"]["cluster_id"].tolist()
isi_vip   = classified[classified["cckbc_tier"] == "ISI VIP"]["cluster_id"].tolist()

print(f"\nhigh-confidence CCKBC clusters: {high_conf}")
print(f"tentative CCKBC clusters:       {tentative}")
print(f"ISI VIP clusters:               {isi_vip}")

# %% [markdown]
# ## 6. Calculate 22q Bias Per Cluster
#
# For each gene set, compute the weighted mean specificity (`HumanCT_AvgZ_Weighted`)
# for every cell type, then extract CGE cluster rows.

# %%
# TODO: define which cell-type columns in SpecMat correspond to CGE clusters
# Replace this with the actual CGE cluster column indices from the annotation file
# CGE_COL_INDICES = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values
CGE_COL_INDICES = None  # placeholder

bias_results = {}

for gs_name, gw_df in gene_sets.items():
    print(f"\nComputing bias for {gs_name}...")
    bias = HumanCT_AvgZ_Weighted(gw_df, SpecMat)
    bias_results[gs_name] = bias
    if CGE_COL_INDICES is not None:
        cge_bias = bias.loc[CGE_COL_INDICES]
        print(f"  CGE bias range: [{cge_bias['EFFECT'].min():.4f}, {cge_bias['EFFECT'].max():.4f}]")
    else:
        print(f"  Bias computed for {len(bias)} cell types (TODO: set CGE_COL_INDICES)")

# %%
# Build per-cluster bias table
# TODO: once CGE_COL_INDICES is defined, pivot to cluster x gene_set table
# bias_cge = pd.DataFrame({
#     gs_name: bias_results[gs_name].loc[CGE_COL_INDICES, "EFFECT"]
#     for gs_name in bias_results
# })
# bias_cge.index.name = "cluster_id"
#
# # Merge with classification
# cge_full = classified.set_index("cluster_id").join(bias_cge, how="inner")
# print(f"Full CGE table: {cge_full.shape}")
# cge_full.head()

# Placeholder for downstream cells
cge_full = classified.copy()
print("NOTE: cge_full is a placeholder. Populate CGE_COL_INDICES to enable bias merging.")

# %% [markdown]
# ## 7. Group Comparisons

# %% [markdown]
# ### 7a. Primary: Mann-Whitney CCKBC vs ISI VIP (per gene set)

# %%
# TODO: run after bias_cge is populated
# primary_results = {}
# for gs_name in gene_sets:
#     cckbc_bias = cge_full.loc[
#         cge_full["cckbc_tier"].isin(["high-confidence CCKBC", "tentative CCKBC"]), gs_name
#     ].dropna().values
#     isi_bias = cge_full.loc[cge_full["cckbc_tier"] == "ISI VIP", gs_name].dropna().values
#
#     if len(cckbc_bias) >= 2 and len(isi_bias) >= 2:
#         stat, pval = stats.mannwhitneyu(isi_bias, cckbc_bias, alternative="greater")
#         primary_results[gs_name] = {
#             "U": stat,
#             "p_value": pval,
#             "n_cckbc": len(cckbc_bias),
#             "n_isi": len(isi_bias),
#             "median_cckbc": np.median(cckbc_bias),
#             "median_isi": np.median(isi_bias),
#         }
#         print(f"  {gs_name}: ISI > CCKBC p = {pval:.4e} "
#               f"(n_CCKBC={len(cckbc_bias)}, n_ISI={len(isi_bias)})")
#     else:
#         print(f"  {gs_name}: insufficient data (n_CCKBC={len(cckbc_bias)}, n_ISI={len(isi_bias)})")
#
# primary_df = pd.DataFrame(primary_results).T
# primary_df

# %% [markdown]
# ### 7b. Secondary: 3-Way Comparison (Kruskal-Wallis)

# %%
# TODO: run after bias_cge is populated
# for gs_name in gene_sets:
#     high_vals = cge_full.loc[cge_full["cckbc_tier"] == "high-confidence CCKBC", gs_name].dropna().values
#     tent_vals = cge_full.loc[cge_full["cckbc_tier"] == "tentative CCKBC", gs_name].dropna().values
#     isi_vals  = cge_full.loc[cge_full["cckbc_tier"] == "ISI VIP", gs_name].dropna().values
#
#     groups = [g for g in [high_vals, tent_vals, isi_vals] if len(g) >= 2]
#     if len(groups) >= 2:
#         h_stat, p_kw = stats.kruskal(*groups)
#         print(f"  {gs_name}: Kruskal-Wallis H={h_stat:.3f}, p={p_kw:.4e}")

# %% [markdown]
# ### 7c. Correlation: Confidence Score vs 22q Bias (Spearman)

# %%
# TODO: run after bias_cge is populated
# for gs_name in gene_sets:
#     if gs_name in cge_full.columns:
#         rho, pval = stats.spearmanr(
#             cge_full["cckbc_confidence"],
#             cge_full[gs_name],
#             nan_policy="omit",
#         )
#         print(f"  {gs_name}: Spearmans' R={rho:.3f}, p={pval:.4e}")

# %% [markdown]
# ## 8. Sensitivity Analysis
#
# Repeat primary comparisons using only high-confidence CCKBC clusters (score >= 3)
# vs. only ISI VIP clusters (score == 0), excluding tentative assignments.

# %%
# TODO: run after bias_cge is populated
# sensitivity_results = {}
# for gs_name in gene_sets:
#     hc_bias  = cge_full.loc[cge_full["cckbc_tier"] == "high-confidence CCKBC", gs_name].dropna().values
#     isi_bias = cge_full.loc[cge_full["cckbc_tier"] == "ISI VIP", gs_name].dropna().values
#
#     if len(hc_bias) >= 2 and len(isi_bias) >= 2:
#         stat, pval = stats.mannwhitneyu(isi_bias, hc_bias, alternative="greater")
#         sensitivity_results[gs_name] = {
#             "p_value_sensitivity": pval,
#             "n_hc_cckbc": len(hc_bias),
#             "n_isi": len(isi_bias),
#         }
#
# if sensitivity_results:
#     sens_df = pd.DataFrame(sensitivity_results).T
#     print("Sensitivity analysis (high-confidence CCKBC only):")
#     print(sens_df)

# %% [markdown]
# ## 9. Visualization: 3-Panel Figure

# %%
# TODO: run after bias_cge is populated and gene sets are available
# This section generates a 3-panel composite figure:
#   Panel A: Boxplot — 22q bias by CCKBC tier (for first available gene set)
#   Panel B: Scatter — CCKBC confidence score vs 22q bias
#   Panel C: Evidence heatmap — clusters x evidence columns

# Panel C: evidence heatmap (works even with placeholder data)
fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

ev_cols = [c for c in EVIDENCE_COLUMNS if c in classified.columns]
if ev_cols:
    heatmap_data = classified.set_index("cluster_id")[ev_cols].astype(int)
    sns.heatmap(
        heatmap_data,
        ax=ax,
        cmap="Blues",
        linewidths=0.5,
        linecolor="lightgray",
        cbar_kws={"shrink": 0.6, "label": "Evidence present"},
        vmin=0, vmax=1,
    )
    ax.set_title("Module E Evidence Heatmap (CGE clusters)", fontsize=11)
    ax.set_xlabel("Evidence Column", fontsize=9)
    ax.set_ylabel("Cluster ID", fontsize=9)
    ax.tick_params(axis="x", rotation=30, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=7)
else:
    ax.text(0.5, 0.5, "No evidence columns present\n(all placeholder data)",
            ha="center", va="center", transform=ax.transAxes, fontsize=12)
    ax.set_title("Evidence Heatmap (placeholder)")

plt.tight_layout()
plt.show()

# TODO: full 3-panel figure code
# fig, axes = plt.subplots(1, 3, figsize=(15, 5))
# fig.patch.set_alpha(0)
# for ax in axes:
#     ax.patch.set_alpha(0)
#
# # Panel A: Boxplot — 22q bias by tier
# gs_name = list(gene_sets.keys())[0]  # first available gene set
# tier_order = ["high-confidence CCKBC", "tentative CCKBC", "ISI VIP"]
# tier_data = [
#     cge_full.loc[cge_full["cckbc_tier"] == t, gs_name].dropna().values
#     for t in tier_order
# ]
# bp = axes[0].boxplot(tier_data, patch_artist=True, showfliers=False,
#                      medianprops=dict(color="black", linewidth=1.4))
# for patch, color in zip(bp["boxes"], ["#d62728", "#ff7f0e", "#1f77b4"]):
#     patch.set_facecolor(color)
#     patch.set_alpha(0.6)
# axes[0].set_xticks([1, 2, 3])
# axes[0].set_xticklabels(["High-conf\nCCKBC", "Tentative\nCCKBC", "ISI VIP"], fontsize=8)
# axes[0].set_ylabel(f"{gs_name} Bias", fontsize=9)
# axes[0].set_title(f"Panel A: {gs_name} Bias by CCKBC Tier", fontsize=10)
# axes[0].spines[["top", "right"]].set_visible(False)
#
# # Panel B: Scatter — confidence vs bias
# axes[1].scatter(
#     cge_full["cckbc_confidence"],
#     cge_full[gs_name],
#     c=cge_full["cckbc_confidence"],
#     cmap="RdYlBu_r",
#     s=60, edgecolors="black", linewidths=0.5,
# )
# axes[1].set_xlabel("CCKBC Confidence Score", fontsize=9)
# axes[1].set_ylabel(f"{gs_name} Bias", fontsize=9)
# axes[1].set_title("Panel B: Confidence vs Bias", fontsize=10)
# axes[1].spines[["top", "right"]].set_visible(False)
#
# # Panel C: Evidence heatmap (in ax[2])
# sns.heatmap(heatmap_data, ax=axes[2], cmap="Blues", linewidths=0.5, cbar=False, vmin=0, vmax=1)
# axes[2].set_title("Panel C: Evidence Heatmap", fontsize=10)
#
# fig.suptitle("Module E: Multi-Modal CCKBC Classification and 22q Bias", fontsize=13, y=1.02)
# plt.tight_layout()
# # TODO: save figure
# # fig.savefig(_PROJ_DIR / "cge_subtype" / "results" / "figures" / "module_e_3panel.pdf",
# #             bbox_inches="tight", transparent=True)
# plt.show()

# %% [markdown]
# ## 10. Save Results

# %%
# TODO: save classification table
# OUT_DIR = _PROJ_DIR / "cge_subtype" / "results"
# OUT_DIR.mkdir(parents=True, exist_ok=True)
#
# # Classification table
# classified.to_csv(OUT_DIR / "module_e_cckbc_classification.csv", index=False)
# print(f"Saved classification to: {OUT_DIR / 'module_e_cckbc_classification.csv'}")
#
# # Bias table
# if "cge_full" in dir() and set(gene_sets.keys()).issubset(cge_full.columns):
#     cge_full.to_csv(OUT_DIR / "module_e_22q_bias_by_cluster.csv")
#     print(f"Saved bias table to: {OUT_DIR / 'module_e_22q_bias_by_cluster.csv'}")
#
# # Primary comparison statistics
# if "primary_df" in dir():
#     primary_df.to_csv(OUT_DIR / "module_e_primary_stats.csv")
#     print(f"Saved primary stats to: {OUT_DIR / 'module_e_primary_stats.csv'}")

print("NOTE: Uncomment save section once cge_full and bias data are populated.")
print(f"Expected output paths:")
print(f"  cge_subtype/results/module_e_cckbc_classification.csv")
print(f"  cge_subtype/results/module_e_22q_bias_by_cluster.csv")
print(f"  cge_subtype/results/module_e_primary_stats.csv")
print(f"  cge_subtype/results/figures/module_e_3panel.pdf")
