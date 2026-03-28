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

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# %% [markdown]
# # GWAS vs Rare Mutation Cell-Type Bias Comparison
#
# This notebook compares cell-type bias profiles derived from **rare mutations**
# (ASD, SCZ, DDD, 22q11.2) with **GWAS common-variant** results from a
# Nature Communications paper (NC). The NC data provides MAGMA-style cell-type
# association betas for 55 GWAS traits across the same 461 cell types.
#
# **Goals:**
# 1. Correlate rare-mutation bias with GWAS bias for matching disorders (e.g., SCZ rare vs SCZ GWAS)
# 2. Systematically compute correlations of each rare-mutation bias with all 55 GWAS traits
# 3. Examine whether neuronal vs non-neuronal cell types drive the correlations
# 4. Visualize cross-trait correlation structure via heatmap and clustering

# %% [markdown]
# ## 1. Load Data

# %%
# Expression specificity matrix
HumanCT_Z2_HCT = pd.read_csv(
    "/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/"
    "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    index_col=0
)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
print("Expression matrix:", HumanCT_Z2_HCT.shape)

# %%
# Gene weight files
GeneWeightDIR = str(PROJ_DIR / "dat/GeneWeights/") + "/"

# ASD (HIQ and LIQ)
HIQ_GW = Fil2Dict("{}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))
LIQ_GW = Fil2Dict("{}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))

# SCZ
SCZ_GW = Fil2Dict("{}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(GeneWeightDIR))

# DDD
DDD_top61_GW = Fil2Dict("{}/DDD.hc.gw".format(GeneWeightDIR))

# 22q11.2
X22q_GW = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/X22q.gw.csv"))

# Negative controls (non-brain rare-variant associations)
NegCtrl_HDL_GW = Fil2Dict("{}/NegCtrl_HDL.gw.csv".format(GeneWeightDIR))
NegCtrl_Alanine_GW = Fil2Dict("{}/NegCtrl_Alanine.gw.csv".format(GeneWeightDIR))
NegCtrl_HbA1c_GW = Fil2Dict("{}/NegCtrl_HbA1c.gw.csv".format(GeneWeightDIR))
NegCtrl_T2D_GW = Fil2Dict("{}/NegCtrl_T2D.gw.csv".format(GeneWeightDIR))
NegCtrl_IBD_GW = Fil2Dict("{}/NegCtrl_IBD.gw.csv".format(GeneWeightDIR))
NegCtrl_RBC_GW = Fil2Dict("{}/NegCtrl_RBC.gw.csv".format(GeneWeightDIR))
NegCtrl_Parkinson_GW = Fil2Dict("{}/NegCtrl_Parkinson.gw.csv".format(GeneWeightDIR))
NegCtrl_Alzheimer_GW = Fil2Dict("{}/NegCtrl_Alzheimer.gw.csv".format(GeneWeightDIR))

print(f"HIQ: {len(HIQ_GW)} genes, LIQ: {len(LIQ_GW)} genes, "
      f"SCZ: {len(SCZ_GW)} genes, DDD: {len(DDD_top61_GW)} genes, "
      f"22q: {len(X22q_GW)} genes")
print(f"Negative controls — HDL: {len(NegCtrl_HDL_GW)}, Alanine: {len(NegCtrl_Alanine_GW)}, "
      f"HbA1c: {len(NegCtrl_HbA1c_GW)}, T2D: {len(NegCtrl_T2D_GW)}, "
      f"IBD: {len(NegCtrl_IBD_GW)}, RBC: {len(NegCtrl_RBC_GW)}, "
      f"Parkinson: {len(NegCtrl_Parkinson_GW)}, Alzheimer: {len(NegCtrl_Alzheimer_GW)}")

# %%
# Compute cell-type bias for each gene set
# Psychiatric disorders (rare mutations)
disorder_sets = [
    ("ASD_HIQ", HIQ_GW), ("ASD_LIQ", LIQ_GW),
    ("SCZ", SCZ_GW), ("DDD", DDD_top61_GW), ("22q", X22q_GW),
]
# Negative controls (rare-variant associations for non-brain traits)
negctrl_sets = [
    ("NC_HDL", NegCtrl_HDL_GW), ("NC_Alanine", NegCtrl_Alanine_GW),
    ("NC_HbA1c", NegCtrl_HbA1c_GW), ("NC_T2D", NegCtrl_T2D_GW),
    ("NC_IBD", NegCtrl_IBD_GW), ("NC_RBC", NegCtrl_RBC_GW),
    ("NC_Parkinson", NegCtrl_Parkinson_GW), ("NC_Alzheimer", NegCtrl_Alzheimer_GW),
]

rare_biases = {}
for name, gw in disorder_sets + negctrl_sets:
    bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
    bias = AnnotateCTDat(bias, Anno)
    rare_biases[name] = bias
    print(f"{name}: {bias.shape[0]} cell types, {len(gw)} genes")

disorder_names = [n for n, _ in disorder_sets]
negctrl_names = [n for n, _ in negctrl_sets]

# %%
# Load NC GWAS cell-type bias (Nature Communications)
NC_Bias = pd.read_csv(
    "/home/jw3514/Work/NeuralP/notebooks/data/NC_Bias_Beta.csv", index_col=0
)
print(f"NC GWAS bias: {NC_Bias.shape[0]} cell types × {NC_Bias.shape[1]} traits")
print(f"\nTraits: {list(NC_Bias.columns)}")

# %% [markdown]
# ## 2. Direct Comparisons: Rare vs GWAS for Matching Disorders

# %%
# SCZ rare vs SCZ GWAS
values1 = rare_biases["SCZ"].sort_index()["EFFECT"]
values2 = NC_Bias.sort_index()["scz2022"]
plot_correlation(values1, values2, "SCZ Rare Mutation", "SCZ GWAS", None, dpi=80)
plt.tight_layout()
plt.show()

# %%
# EDU GWAS vs SCZ GWAS (cognitive–psychiatric overlap)
values1 = NC_Bias.sort_index()["educational_attainment"]
values2 = NC_Bias.sort_index()["scz2022"]
plot_correlation(values1, values2, "EDU GWAS", "SCZ GWAS", None, dpi=80)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Systematic Trait Correlations
#
# For each rare-mutation disorder, compute Spearman correlations with all 55 GWAS traits,
# both across all cell types and restricted to neuronal cell types only.

# %%
def compute_trait_correlations(values1, nc_bias_df, neur_idx=Neur_idx):
    """
    Compute Spearman correlations between a rare-mutation bias vector
    and each GWAS trait column in nc_bias_df.

    Returns DataFrame with columns: Trait, r_all, p_all, r_neur, p_neur
    """
    results = []
    for trait in nc_bias_df.columns:
        values2 = nc_bias_df.sort_index()[trait].values
        r_all, p_all = spearmanr(values1, values2)
        r_neur, p_neur = spearmanr(values1[neur_idx], values2[neur_idx])
        results.append({
            "Trait": trait,
            "r_all": r_all, "p_all": p_all,
            "r_neur": r_neur, "p_neur": p_neur
        })
    df = pd.DataFrame(results)
    # FDR correction
    df["q_all"] = multipletests(df["p_all"], method="fdr_bh")[1]
    df["q_neur"] = multipletests(df["p_neur"], method="fdr_bh")[1]
    return df.sort_values("r_all", ascending=False).reset_index(drop=True)


# %%
# Compute for each disorder
corr_results = {}
for name, bias_df in rare_biases.items():
    values1 = bias_df.sort_index()["EFFECT"].values
    corr_results[name] = compute_trait_correlations(values1, NC_Bias)
    print(f"\n=== {name} — Top 10 correlated GWAS traits (all cell types) ===")
    display(corr_results[name].head(10))

# %% [markdown]
# ## 4. Correlation Heatmap: Disorders vs Negative Controls × GWAS Traits

# %%
# Build r_all matrix: rows=GWAS traits, cols=all gene sets (disorders + neg controls)
r_matrix = pd.DataFrame({
    name: corr_results[name].set_index("Trait")["r_all"]
    for name in rare_biases
})

# Separate columns for ordered display
r_matrix_ordered = r_matrix[disorder_names + negctrl_names]

# Select GWAS traits significant (q < 0.1) for at least one gene set
sig_traits = set()
for name in rare_biases:
    sig = corr_results[name].loc[corr_results[name]["q_all"] < 0.1, "Trait"]
    sig_traits.update(sig)
print(f"{len(sig_traits)} GWAS traits significant (q<0.1) for ≥1 gene set")

# If none pass threshold, use top 25 by max absolute r across all gene sets
if len(sig_traits) < 5:
    r_abs_max = r_matrix.abs().max(axis=1)
    sig_traits = set(r_abs_max.nlargest(25).index)
    print(f"Using top 25 traits by max |r| instead ({len(sig_traits)} traits)")

r_matrix_sig = r_matrix_ordered.loc[sorted(sig_traits)]

# %%
fig, ax = plt.subplots(figsize=(12, max(6, len(r_matrix_sig) * 0.35)), dpi=150)
sns.heatmap(
    r_matrix_sig,
    cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
    annot=True, fmt=".2f", linewidths=0.5,
    ax=ax
)
# Add vertical line to separate disorders from negative controls
ax.axvline(x=len(disorder_names), color="black", linewidth=2)
ax.set_title("Spearman r: Rare-Variant Bias × GWAS Trait Bias\n"
             "(left: psychiatric disorders | right: negative controls)")
ax.set_xlabel("Gene Set")
ax.set_ylabel("GWAS Trait")
ax.patch.set_alpha(0)
fig.patch.set_alpha(0)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Disorder vs Negative Control: Distribution of |r| with Brain vs Non-Brain GWAS Traits
#
# Do psychiatric rare-mutation gene sets have stronger correlations with
# brain/psychiatric GWAS traits than negative controls?

# %%
# Categorize GWAS traits as brain-related or non-brain
brain_gwas = [
    "adhd2023", "bip2021", "bip2021_type1", "bip2021_type2",
    "educational_attainment", "iq", "mdd2019", "mdd_ect", "mdd_ppd",
    "neuroticism", "scz2022", "insomnia", "epilepsy",
    "alzheimers_disease", "parkinson2019", "amyotrophic_lateral_sclerosis",
    "suicide_attempt", "suicide_death", "with_suicidal_thought", "without_suicidal_thought",
    "cortical_surface_area", "cortical_thickness",
    "volume_accumbens", "volume_amygdala", "volume_brainstem",
    "volume_caudate", "volume_pallidum", "volume_putamen", "volume_thalamus",
    "multiple_sclerosis", "migraine", "hearing_loss",
]
nonbrain_gwas = [t for t in NC_Bias.columns if t not in brain_gwas]
print(f"Brain/psych GWAS: {len(brain_gwas)}, Non-brain GWAS: {len(nonbrain_gwas)}")
print(f"Non-brain: {nonbrain_gwas}")

# %%
# For each gene set, compute mean |r| with brain vs non-brain GWAS traits
summary_rows = []
for name in disorder_names + negctrl_names:
    r_vals = corr_results[name].set_index("Trait")["r_all"]
    mean_r_brain = r_vals.loc[[t for t in brain_gwas if t in r_vals.index]].abs().mean()
    mean_r_nonbrain = r_vals.loc[[t for t in nonbrain_gwas if t in r_vals.index]].abs().mean()
    group = "Disorder" if name in disorder_names else "Neg Control"
    summary_rows.append({
        "Gene Set": name, "Group": group,
        "mean_|r|_brain_GWAS": mean_r_brain,
        "mean_|r|_nonbrain_GWAS": mean_r_nonbrain,
    })
summary_df = pd.DataFrame(summary_rows)
display(summary_df)

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=150)

# Panel A: mean |r| with brain GWAS
colors = ["#d62728" if g == "Disorder" else "#7f7f7f" for g in summary_df["Group"]]
axes[0].barh(summary_df["Gene Set"], summary_df["mean_|r|_brain_GWAS"], color=colors)
axes[0].set_xlabel("Mean |Spearman r|")
axes[0].set_title("Brain/Psychiatric GWAS Traits")
axes[0].invert_yaxis()
axes[0].patch.set_alpha(0)

# Panel B: mean |r| with non-brain GWAS
axes[1].barh(summary_df["Gene Set"], summary_df["mean_|r|_nonbrain_GWAS"], color=colors)
axes[1].set_xlabel("Mean |Spearman r|")
axes[1].set_title("Non-Brain GWAS Traits")
axes[1].invert_yaxis()
axes[1].patch.set_alpha(0)

fig.suptitle("Disorders vs Negative Controls: Correlation with GWAS Traits", y=1.02)
fig.patch.set_alpha(0)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Neuronal vs Non-Neuronal Correlations

# %%
# Compare r_all vs r_neur — disorders only
fig, axes = plt.subplots(1, len(disorder_names), figsize=(4 * len(disorder_names), 4), dpi=150)
for ax, name in zip(axes, disorder_names):
    df = corr_results[name]
    ax.scatter(df["r_all"], df["r_neur"], alpha=0.6, s=30, c="steelblue")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    lim = max(df["r_all"].abs().max(), df["r_neur"].abs().max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5, alpha=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("r (all cell types)")
    ax.set_ylabel("r (neuronal only)")
    ax.set_title(name)
    ax.set_aspect("equal")
    ax.patch.set_alpha(0)

    # Label top 3 outliers (largest |r_all - r_neur|)
    df_tmp = df.copy()
    df_tmp["diff"] = (df_tmp["r_all"] - df_tmp["r_neur"]).abs()
    for _, row in df_tmp.nlargest(3, "diff").iterrows():
        ax.annotate(row["Trait"], (row["r_all"], row["r_neur"]),
                    fontsize=6, alpha=0.8)

fig.patch.set_alpha(0)
fig.suptitle("All vs Neuronal-Only Correlations with GWAS Traits (Disorders)", y=1.02)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 7. Clustered Heatmap: All GWAS Traits × All Gene Sets

# %%
# Full heatmap with hierarchical clustering
g = sns.clustermap(
    r_matrix_ordered.T,
    cmap="RdBu_r", center=0, vmin=-0.5, vmax=0.5,
    figsize=(max(14, len(r_matrix) * 0.25), 6),
    dendrogram_ratio=(0.1, 0.15),
    linewidths=0.3,
    yticklabels=True,
    xticklabels=True,
    cbar_kws={"label": "Spearman r"}
)
g.ax_heatmap.set_xlabel("GWAS Trait")
g.ax_heatmap.set_ylabel("Gene Set (Disorder / Neg Control)")
g.ax_heatmap.tick_params(axis='x', labelsize=7, rotation=90)
plt.suptitle("Clustered: Rare-Variant × GWAS Trait Correlations", y=1.02)
plt.show()

# %% [markdown]
# ## 8. Per-Disorder Tables (Disorders)

# %%
for name in disorder_names:
    print(f"\n{'='*60}")
    print(f"  {name}: Top 15 correlated GWAS traits")
    print(f"{'='*60}")
    display(corr_results[name].head(15).style.format({
        "r_all": "{:.3f}", "p_all": "{:.2e}", "q_all": "{:.3f}",
        "r_neur": "{:.3f}", "p_neur": "{:.2e}", "q_neur": "{:.3f}"
    }).background_gradient(subset=["r_all", "r_neur"], cmap="RdBu_r", vmin=-0.5, vmax=0.5))

# %% [markdown]
# ## 9. Per-Disorder Tables (Negative Controls)

# %%
for name in negctrl_names:
    print(f"\n{'='*60}")
    print(f"  {name}: Top 10 correlated GWAS traits")
    print(f"{'='*60}")
    display(corr_results[name].head(10).style.format({
        "r_all": "{:.3f}", "p_all": "{:.2e}", "q_all": "{:.3f}",
        "r_neur": "{:.3f}", "p_neur": "{:.2e}", "q_neur": "{:.3f}"
    }).background_gradient(subset=["r_all", "r_neur"], cmap="RdBu_r", vmin=-0.5, vmax=0.5))

# %% [markdown]
# ## 10. Scatter Plots: Top Correlated GWAS Traits per Disorder

# %%
for name in disorder_names:
    bias_df = rare_biases[name]
    top_traits = corr_results[name].head(3)["Trait"].values
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=100)
    for ax, trait in zip(axes, top_traits):
        values1 = bias_df.sort_index()["EFFECT"]
        values2 = NC_Bias.sort_index()[trait]
        plt.sca(ax)
        plot_correlation(values1, values2, f"{name} Rare", f"{trait} GWAS", None, dpi=100)
    fig.suptitle(f"{name}: Top 3 correlated GWAS traits", fontsize=14, y=1.02)
    fig.patch.set_alpha(0)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# ## 11. Cross-Gene-Set Correlation Matrix (Disorders + Neg Controls)

# %%
all_names = disorder_names + negctrl_names
cross_corr = pd.DataFrame(np.ones((len(all_names), len(all_names))),
                           index=all_names, columns=all_names)
for i, n1 in enumerate(all_names):
    for j, n2 in enumerate(all_names):
        if i < j:
            v1 = rare_biases[n1].sort_index()["EFFECT"].values
            v2 = rare_biases[n2].sort_index()["EFFECT"].values
            r, _ = spearmanr(v1, v2)
            cross_corr.loc[n1, n2] = r
            cross_corr.loc[n2, n1] = r

fig, ax = plt.subplots(figsize=(10, 9), dpi=150)
sns.heatmap(cross_corr, annot=True, fmt=".2f", cmap="RdBu_r",
            center=0, vmin=-0.5, vmax=1, square=True, ax=ax,
            linewidths=0.5)
# Divider between disorders and neg controls
ax.axhline(y=len(disorder_names), color="black", linewidth=2)
ax.axvline(x=len(disorder_names), color="black", linewidth=2)
ax.set_title("Cell-Type Bias Correlation: Disorders vs Negative Controls")
ax.patch.set_alpha(0)
fig.patch.set_alpha(0)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 12. Summary: Do Negative Controls Preferentially Correlate with Non-Brain GWAS?
#
# Key question: psychiatric disorder gene sets should correlate with brain/psychiatric
# GWAS traits, while negative controls should either show no pattern or correlate
# with their matching non-brain GWAS trait.

# %%
# For each gene set, find the top 5 most correlated GWAS traits
print("=" * 80)
print("TOP 5 CORRELATED GWAS TRAITS PER GENE SET")
print("=" * 80)
for name in disorder_names + negctrl_names:
    top5 = corr_results[name].head(5)
    group = "DISORDER" if name in disorder_names else "NEG CTRL"
    print(f"\n[{group}] {name}:")
    for _, row in top5.iterrows():
        sig = "*" if row["q_all"] < 0.1 else ""
        brain_flag = "brain" if row["Trait"] in brain_gwas else "non-brain"
        print(f"  {row['Trait']:35s}  r={row['r_all']:+.3f}  q={row['q_all']:.3f}{sig}  ({brain_flag})")
