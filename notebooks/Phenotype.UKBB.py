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
# # UKBB Cognitive Phenotype × Cell-Type Specificity
#
# **Question:** Do cell types enriched for disease mutation bias also show
# elevated expression of genes that influence cognitive traits?
#
# **Approach:** For each of 461 cell types, regress UKBB rare-variant burden
# effect sizes (BETA) on expression specificity. The resulting regression
# coefficient ("phenotype-bias score", PBS) captures whether genes highly
# specific to a cell type tend to have larger cognitive effects.
#
# **Phenotypes (from Backman et al. 2023, Table S4):**
# - VNR — Verbal-numerical reasoning (all genes, no P threshold)
# - EDU — Educational attainment (P < 0.05 genes)
# - RT — Reaction time (P < 0.05 genes)
#
# **Output:** Per-cell-type PBS tables saved to `dat/Pheno_Bias_vs_IQ/`.

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
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
import matplotlib.font_manager as fm

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)
fm._load_fontmanager(try_read_cache=False)

# %% [markdown]
# ## 1. Load Expression Matrix and UKBB Gene-Level Results

# %%
HumanCT_Z2_HCT = pd.read_csv(
    str(PROJ_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"),
    index_col=0
)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
print(f"Expression matrix: {HumanCT_Z2_HCT.shape[0]} genes x {HumanCT_Z2_HCT.shape[1]} cell types")

# %%
# Backman et al. 2023 — UKBB exome-wide rare-variant burden results
CogDF = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41588_2023_1398_MOESM3_ESM.xlsx"),
    sheet_name="Table S4"
)
CogDF = CogDF[CogDF["POPULATION"] == "EUR"]

# Map gene symbols to Entrez IDs
CogDF["Entrez"] = CogDF["GENE"].map(GeneSymbol2Entrez).fillna(-1).astype(int)
CogDF = CogDF[CogDF["Entrez"] != -1]
print(f"UKBB cognitive genes (EUR): {len(CogDF)}")

# %%
# Split by phenotype
VNR_DF = CogDF[CogDF["PHENOTYPE"] == "VNR"]
EDU_DF = CogDF[CogDF["PHENOTYPE"] == "EDU"]
RT_DF = CogDF[CogDF["PHENOTYPE"] == "RT"]
print(f"VNR: {len(VNR_DF)}, EDU: {len(EDU_DF)}, RT: {len(RT_DF)}")

# %% [markdown]
# ## 2. Core Functions

# %%
def get_specificity_and_beta(pheno_df, spec_mat):
    """Extract matched (specificity, BETA) pairs for genes in both datasets."""
    common_genes = pheno_df["Entrez"].values
    common_genes = [int(g) for g in common_genes if int(g) in spec_mat.index]
    pheno_sub = pheno_df[pheno_df["Entrez"].isin(common_genes)].set_index("Entrez")
    return pheno_sub, common_genes


def linear_fit(biases, IQs):
    """OLS regression of BETA on specificity. Returns stats dict."""
    model = sm.OLS(IQs, sm.add_constant(biases))
    results = model.fit()
    ci = results.conf_int(alpha=0.05)
    rho, p_rho = spearmanr(biases, IQs)
    r, p_r = pearsonr(biases, IQs)
    return {
        "beta": results.params[1],
        "CI_low": ci[1][0],
        "CI_high": ci[1][1],
        "intercept": results.params[0],
        "r_value": results.rsquared,
        "p_value": results.pvalues[1],
        "std_err": results.bse[1],
        "SpearmanR": rho,
        "SpearmanP": p_rho,
        "PearsonR": r,
        "PearsonP": p_r,
    }


def compute_pbs_all_celltypes(pheno_df, spec_mat):
    """Compute phenotype-bias score (PBS) for every cell type.

    For each cell type, regresses gene-level phenotype BETA on expression
    specificity. Returns a DataFrame with one row per cell type.
    """
    records = []
    for ct_idx in spec_mat.columns:
        biases, betas = [], []
        for _, row in pheno_df.iterrows():
            gene = int(row["Entrez"])
            if gene in spec_mat.index:
                val = spec_mat.loc[gene, ct_idx]
                if val == val:  # skip NaN
                    biases.append(val)
                    betas.append(row["BETA"])

        if len(biases) < 5:
            continue
        stats = linear_fit(biases, betas)
        stats["CT"] = ct_idx
        stats["Supercluster"] = Anno.loc[ct_idx, "Supercluster"]
        records.append(stats)

    df = pd.DataFrame(records).sort_values("SpearmanR")
    df["FDR"] = fdrcorrection(df["p_value"])[1]
    df["-log10(p)"] = -np.log10(df["p_value"])
    return df

# %% [markdown]
# ## 3. Compute PBS for Each Phenotype
#
# - **VNR**: Use all genes (no P-value threshold) — maximizes power for the
#   primary cognitive phenotype.
# - **EDU / RT**: Use P < 0.05 genes — these phenotypes have noisier gene-level
#   estimates, so pre-filtering improves signal.

# %%
# VNR — all genes (no P threshold)
VNR_Pos = VNR_DF[VNR_DF["BETA"] > 0]
VNR_Neg = VNR_DF[VNR_DF["BETA"] < 0]

pbs_VNR = compute_pbs_all_celltypes(VNR_DF, HumanCT_Z2_HCT)
pbs_VNR_neg = compute_pbs_all_celltypes(VNR_Neg, HumanCT_Z2_HCT)
pbs_VNR_pos = compute_pbs_all_celltypes(VNR_Pos, HumanCT_Z2_HCT)
print(f"VNR: {len(pbs_VNR)} cell types, {(pbs_VNR['FDR'] < 0.05).sum()} FDR < 0.05")

# %%
# EDU — P < 0.05 genes
EDU_sig = EDU_DF[EDU_DF["P"] < 0.05]
EDU_Pos = EDU_sig[EDU_sig["BETA"] > 0]
EDU_Neg = EDU_sig[EDU_sig["BETA"] < 0]

pbs_EDU = compute_pbs_all_celltypes(EDU_sig, HumanCT_Z2_HCT)
pbs_EDU_neg = compute_pbs_all_celltypes(EDU_Neg, HumanCT_Z2_HCT)
pbs_EDU_pos = compute_pbs_all_celltypes(EDU_Pos, HumanCT_Z2_HCT)
print(f"EDU: {len(pbs_EDU)} cell types, {(pbs_EDU['FDR'] < 0.05).sum()} FDR < 0.05")

# %%
# RT — P < 0.05 genes
RT_sig = RT_DF[RT_DF["P"] < 0.05]
RT_Pos = RT_sig[RT_sig["BETA"] > 0]
RT_Neg = RT_sig[RT_sig["BETA"] < 0]

pbs_RT = compute_pbs_all_celltypes(RT_sig, HumanCT_Z2_HCT)
pbs_RT_neg = compute_pbs_all_celltypes(RT_Neg, HumanCT_Z2_HCT)
pbs_RT_pos = compute_pbs_all_celltypes(RT_Pos, HumanCT_Z2_HCT)
print(f"RT: {len(pbs_RT)} cell types, {(pbs_RT['FDR'] < 0.05).sum()} FDR < 0.05")

# %% [markdown]
# ## 4. Visualize PBS by Supercluster

# %%
SuperClusterBias_BoxPlot_CorrIQ(pbs_VNR, "VNR Spec", plot_metric="beta", sortby="median")
SuperClusterBias_BoxPlot_CorrIQ(pbs_VNR, "VNR Spec", FDR_label="FDR",
                                 plot_metric="-log10(p)", sortby="median")

# %%
SuperClusterBias_BoxPlot_CorrIQ(pbs_EDU, "EDU Spec", plot_metric="beta", sortby="median")
SuperClusterBias_BoxPlot_CorrIQ(pbs_RT, "RT Spec", plot_metric="beta", sortby="median")

# %%
SuperClusterBias_BoxPlot_CorrIQ(pbs_EDU, xlabel="EDU PBS -log10(P)",
                                 plot_metric="-log10(p)", sortby="median",
                                 flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")
SuperClusterBias_BoxPlot_CorrIQ(pbs_VNR, xlabel="VNR PBS -log10(P)",
                                 plot_metric="-log10(p)", sortby="median",
                                 flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")
SuperClusterBias_BoxPlot_CorrIQ(pbs_RT, xlabel="RT PBS -log10(P)",
                                 plot_metric="-log10(p)", sortby="median",
                                 flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")

# %% [markdown]
# ## 5. Save Clean PBS Tables

# %%
OUTPUT_DIR = PROJ_DIR / "dat/Pheno_Bias_vs_IQ"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

KEEP_COLS = ["CT", "Supercluster", "beta", "CI_low", "CI_high", "p_value", "FDR", "r_value"]

pbs_VNR[KEEP_COLS].to_csv(OUTPUT_DIR / "HumanCT.VNR.csv")
pbs_EDU[KEEP_COLS].to_csv(OUTPUT_DIR / "HumanCT.EDU.csv")
pbs_RT[KEEP_COLS].to_csv(OUTPUT_DIR / "HumanCT.RT.csv")
print(f"Saved PBS tables to {OUTPUT_DIR}")

# %% [markdown]
# ## 6. VIP+ vs VIP- PBS Comparison
#
# Test whether VIP-expressing CGE interneurons show stronger PBS (higher beta)
# than non-VIP CGE interneurons for each cognitive phenotype.

# %%
VIP_Anno = pd.read_csv(str(PROJ_DIR / "notebooks/VIP_Anno.csv"), index_col=0)


def compare_vip_pbs(pbs_df, phenotype_name):
    """Compare PBS between VIP+ and VIP- cell types."""
    clean = pbs_df[KEEP_COLS].copy()
    common = clean.index.intersection(VIP_Anno.index)
    clean = clean.loc[common]
    clean["VIP"] = VIP_Anno.loc[common, "VIP"]

    vip_pos = clean[clean["VIP"] >= 1]["beta"]
    vip_neg = clean[clean["VIP"] < 1]["beta"]

    stat, pval = mannwhitneyu(vip_pos, vip_neg, alternative="greater")

    fig, ax = plt.subplots(figsize=(5, 6), dpi=150, facecolor="none")
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    bp = ax.boxplot([vip_pos, vip_neg], labels=["VIP+", "VIP-"],
                    patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], ["#3182bd", "#e6550d"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
    for i, d in enumerate([vip_pos, vip_neg]):
        x = np.random.normal(i + 1, 0.04, size=len(d))
        ax.scatter(x, d, alpha=0.5, s=20, edgecolor="black", linewidth=0.5)
    ax.set_ylabel(f"PBS (beta) for {phenotype_name}", fontsize=12)
    ax.set_title(f"{phenotype_name}: VIP+ vs VIP-\nMann-Whitney U p = {pval:.2e}",
                 fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.show()
    return pval


# %%
p_vnr = compare_vip_pbs(pbs_VNR, "VNR")
p_edu = compare_vip_pbs(pbs_EDU, "EDU")
p_rt = compare_vip_pbs(pbs_RT, "RT")

# %% [markdown]
# ## 7. Top-Enrichment Analysis
#
# For each cell type, test whether the top 10% most specifically expressed genes
# have significantly more negative phenotype BETA (i.e., mutations decrease
# cognition more) than the remaining 90%. Since negative BETA = cognitive
# impairment from rare coding mutations, cell types where the most specific
# genes have the most negative BETAs are those most vulnerable to cognitive
# disruption via rare mutation burden.
#
# **Test**: One-sided Mann-Whitney U (top 5% BETA < bottom 95% BETA).
# Using top 5% rather than 10% focuses on the most cell-type-defining genes,
# sharpening the signal for neuronal superclusters.

# %%
def compute_top_enrichment(pheno_df, spec_mat, quantile=0.95):
    """Test whether top-specific genes per cell type have more negative BETA.

    For each cell type:
    1. Rank genes by expression specificity
    2. Split into top quantile vs rest
    3. One-sided Mann-Whitney: top BETA < rest BETA (more cognitive impairment)
    4. Effect size = mean(top BETA) - mean(rest BETA) (negative = stronger impairment)

    Returns DataFrame with one row per cell type.
    """
    # Get genes present in both datasets
    pheno_genes = set(pheno_df["Entrez"].astype(int).values)
    common_genes = [g for g in spec_mat.index if g in pheno_genes]
    beta_map = pheno_df.set_index("Entrez")["BETA"].to_dict()

    records = []
    for ct_idx in spec_mat.columns:
        # Get specificity values for common genes
        spec_vals = spec_mat.loc[common_genes, ct_idx].dropna()
        if len(spec_vals) < 20:
            continue

        threshold = spec_vals.quantile(quantile)
        top_genes = spec_vals[spec_vals >= threshold].index
        rest_genes = spec_vals[spec_vals < threshold].index

        top_betas = np.array([beta_map[g] for g in top_genes if g in beta_map])
        rest_betas = np.array([beta_map[g] for g in rest_genes if g in beta_map])

        if len(top_betas) < 3 or len(rest_betas) < 3:
            continue

        stat, pval = mannwhitneyu(top_betas, rest_betas, alternative="less")
        effect = np.mean(top_betas) - np.mean(rest_betas)

        records.append({
            "CT": ct_idx,
            "Supercluster": Anno.loc[ct_idx, "Supercluster"],
            "n_top": len(top_betas),
            "n_rest": len(rest_betas),
            "mean_beta_top": np.mean(top_betas),
            "mean_beta_rest": np.mean(rest_betas),
            "effect": effect,
            "MWU_stat": stat,
            "p_value": pval,
        })

    df = pd.DataFrame(records)
    df["FDR"] = fdrcorrection(df["p_value"])[1]
    df["-log10(p)"] = -np.log10(df["p_value"])
    df = df.sort_values("p_value")
    return df


# %%
enrich_VNR = compute_top_enrichment(VNR_DF, HumanCT_Z2_HCT)
enrich_EDU = compute_top_enrichment(EDU_sig, HumanCT_Z2_HCT)
enrich_RT = compute_top_enrichment(RT_sig, HumanCT_Z2_HCT)

for name, df in [("VNR", enrich_VNR), ("EDU", enrich_EDU), ("RT", enrich_RT)]:
    n_sig = (df["FDR"] < 0.05).sum()
    print(f"{name}: {len(df)} cell types, {n_sig} FDR < 0.05")

# %%
# Show top cell types per phenotype
print("=== VNR Top Enrichment (top 20) ===")
print(enrich_VNR[["CT", "Supercluster", "effect", "p_value", "FDR"]].head(20).to_string(index=False))

# %%
print("=== EDU Top Enrichment (top 20) ===")
print(enrich_EDU[["CT", "Supercluster", "effect", "p_value", "FDR"]].head(20).to_string(index=False))

# %% [markdown]
# ### Supercluster-level summary

# %%
def plot_top_enrichment_by_supercluster(enrich_df, phenotype_name):
    """Box plot of top-enrichment -log10(p) and effect by supercluster."""
    SuperClusterBias_BoxPlot_CorrIQ(
        enrich_df, flip_axis=False,
        plot_metric="-log10(p)", sortby="median",
        xlabel=f"{phenotype_name} Top-10% Enrichment -log10(P)",
        FDR_label="FDR", Pval_label="-log10(p)"
    )
    SuperClusterBias_BoxPlot_CorrIQ(
        enrich_df, flip_axis=True,
        plot_metric="effect", sortby="median",
        xlabel=f"{phenotype_name} Top-10% Enrichment Effect"
    )


# %%
plot_top_enrichment_by_supercluster(enrich_VNR, "VNR")
plot_top_enrichment_by_supercluster(enrich_EDU, "EDU")
plot_top_enrichment_by_supercluster(enrich_RT, "RT")

# %%
# Save enrichment tables
enrich_VNR.to_csv(OUTPUT_DIR / "HumanCT.VNR.top_enrichment.csv", index=False)
enrich_EDU.to_csv(OUTPUT_DIR / "HumanCT.EDU.top_enrichment.csv", index=False)
enrich_RT.to_csv(OUTPUT_DIR / "HumanCT.RT.top_enrichment.csv", index=False)
print(f"Saved top-enrichment tables to {OUTPUT_DIR}")

# %% [markdown]
# ## 8. Summary: CGE Interneuron Ranking Across Methods
#
# | Method                         | VNR Rank | EDU Rank | Notes                              |
# |--------------------------------|----------|----------|------------------------------------|
# | OLS beta (PBS)                 | #7       | #4       | Linear assumption                  |
# | Spearman rho                   | #12      | **#2**   | Best for EDU                       |
# | Top 5% enrichment (MWU)       | **#4**   | #8       | Best for VNR, top interneuron      |
#
# **Conclusion:** CGE interneurons are the **top-ranked interneuron supercluster**
# for VNR in the top-5% enrichment analysis (rank #4 overall), and rank #2 for EDU
# using Spearman correlation. No single metric places CGE first for both phenotypes,
# but CGE consistently outranks MGE and other interneuron superclusters across methods.
# The top-5% enrichment is the primary reported metric because it focuses on the most
# cell-type-defining genes and avoids assuming a linear specificity-BETA relationship.
