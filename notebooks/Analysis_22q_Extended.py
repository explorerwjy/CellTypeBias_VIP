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
# # Extended 22q11.2 Computational Analysis
#
# This notebook extends the 22q11.2 analysis with four sections:
# 1. **Per-gene contribution** — which 22q genes drive the CGE/VIP bias?
# 2. **CNV comparison** — is VIP bias unique to 22q or shared with other psychiatric CNVs?
# 3. **22q DEG overlap** — do downstream transcriptomic consequences also show VIP/CGE bias?
# 4. **Jackknife stability** — how robust is the VIP bias to removing individual genes?

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
import yaml

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import (
    HumanCT_AvgZ_Weighted, Fil2Dict, AnnotateCTDat, LoadGeneINFO
)

plt.style.use('seaborn-v0_8-whitegrid')

# %%
# --- Load core data ---
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()
Anno = pd.read_excel(str(PROJ_DIR / "dat/annotation.xlsx"), index_col=0)

# Expression specificity matrix (mean-centered, clipped)
SpecMat = pd.read_csv(
    str(PROJ_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv"),
    index_col=0
)
SpecMat.columns = [int(c) for c in SpecMat.columns]

# 22q gene weights (all uniform = 1)
X22q_GW = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/X22q.gw.csv"))
X22q_genes = list(X22q_GW.keys())
X22q_small_GW = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/X22q.mousemodel.gw.csv"))

# CGE cell type indices and VIP grouping
CGE_idx = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values.astype(int)

# Raw expression matrix for VIP marker expression
ExpL = pd.read_csv(
    str(PROJ_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"),
    index_col=0
)
ExpL.columns = [int(c) for c in ExpL.columns]

# VIP expression per CGE cluster → define VIP+/VIP- groups
VIP_CUTOFF = 1.0
vip_expr = ExpL.loc[GeneSymbol2Entrez["VIP"], CGE_idx]
vip_pos_idx = vip_expr[vip_expr >= VIP_CUTOFF].index.values.astype(int)
vip_neg_idx = vip_expr[vip_expr < VIP_CUTOFF].index.values.astype(int)
print(f"CGE clusters: {len(CGE_idx)} total, {len(vip_pos_idx)} VIP+, {len(vip_neg_idx)} VIP-")

# Pre-computed 22q bias results
Bias_Save_Dir = str(PROJ_DIR / "results/main_results/random/Centering/") + "/"
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)

# %% [markdown]
# ---
# ## 1. Per-gene contribution to CGE/VIP bias
#
# For each 22q gene, compute its expression specificity in VIP+ vs VIP- CGE clusters.
# The difference reveals which genes *drive* the overall VIP+ > VIP- bias signal.

# %%
# Compute per-gene mean specificity in VIP+ vs VIP- CGE clusters
valid_22q = [g for g in X22q_genes if g in SpecMat.index]
invalid_22q = [g for g in X22q_genes if g not in SpecMat.index]
if invalid_22q:
    print(f"  {len(invalid_22q)} genes not in expression matrix (dropped)")
print(f"  {len(valid_22q)} valid 22q genes for analysis")

gene_contrib = pd.DataFrame(index=valid_22q)
gene_contrib["symbol"] = [Entrez2Symbol.get(g, str(g)) for g in valid_22q]
gene_contrib["spec_VIP_pos"] = SpecMat.loc[valid_22q, vip_pos_idx].mean(axis=1)
gene_contrib["spec_VIP_neg"] = SpecMat.loc[valid_22q, vip_neg_idx].mean(axis=1)
gene_contrib["delta_VIP"] = gene_contrib["spec_VIP_pos"] - gene_contrib["spec_VIP_neg"]
gene_contrib["spec_CGE"] = SpecMat.loc[valid_22q, CGE_idx].mean(axis=1)

# Overall CGE bias (mean across all cell types vs CGE)
all_idx = SpecMat.columns.values.astype(int)
gene_contrib["spec_all"] = SpecMat.loc[valid_22q, all_idx].mean(axis=1)
gene_contrib["delta_CGE"] = gene_contrib["spec_CGE"] - gene_contrib["spec_all"]

# In the 1.5 Mb deletion?
small_del_genes = set(X22q_small_GW.keys())
gene_contrib["in_1.5Mb"] = [g in small_del_genes for g in valid_22q]

gene_contrib = gene_contrib.sort_values("delta_VIP", ascending=False)

# %%
# --- Figure 1A: Per-gene VIP+ vs VIP- specificity difference (bar chart) ---
fig, ax = plt.subplots(figsize=(14, 5), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

colors = ['#D7191C' if d > 0 else '#2C7BB6' for d in gene_contrib["delta_VIP"]]
bars = ax.bar(range(len(gene_contrib)), gene_contrib["delta_VIP"], color=colors, edgecolor='none', alpha=0.85)

# Mark 1.5Mb deletion genes with edge
for i, (_, row) in enumerate(gene_contrib.iterrows()):
    if row["in_1.5Mb"]:
        bars[i].set_edgecolor('black')
        bars[i].set_linewidth(1.5)

ax.set_xticks(range(len(gene_contrib)))
ax.set_xticklabels(gene_contrib["symbol"], rotation=90, fontsize=8)
ax.set_ylabel("Specificity difference\n(VIP+ − VIP−)", fontsize=12)
ax.set_title("Per-gene contribution to VIP+ bias within CGE interneurons (22q11.2 genes)", fontsize=13)
ax.axhline(0, color='black', lw=0.8)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#D7191C', alpha=0.85, label='VIP+ enriched'),
    Patch(facecolor='#2C7BB6', alpha=0.85, label='VIP− enriched'),
    Patch(facecolor='grey', edgecolor='black', linewidth=1.5, label='In 1.5 Mb deletion'),
]
ax.legend(handles=legend_elements, frameon=True, fontsize=9, loc='upper right')

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_pergene_VIP_contribution.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# Print top/bottom contributors
print("=== Top 10 genes driving VIP+ bias ===")
for _, row in gene_contrib.head(10).iterrows():
    flag = " *1.5Mb" if row["in_1.5Mb"] else ""
    print(f"  {row['symbol']:>10s}  delta_VIP={row['delta_VIP']:+.4f}  CGE_spec={row['spec_CGE']:.4f}{flag}")

print("\n=== Bottom 5 genes (VIP- enriched) ===")
for _, row in gene_contrib.tail(5).iterrows():
    flag = " *1.5Mb" if row["in_1.5Mb"] else ""
    print(f"  {row['symbol']:>10s}  delta_VIP={row['delta_VIP']:+.4f}  CGE_spec={row['spec_CGE']:.4f}{flag}")

# %%
# --- Figure 1B: Scatter plot CGE specificity vs VIP delta ---
fig, ax = plt.subplots(figsize=(7, 6), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

colors_scatter = ['#D7191C' if v else '#4575B4' for v in gene_contrib["in_1.5Mb"]]
ax.scatter(gene_contrib["spec_CGE"], gene_contrib["delta_VIP"],
           c=colors_scatter, s=50, alpha=0.7, edgecolor='white', linewidth=0.5)

# Label top contributors
top_genes = gene_contrib.head(5)
bottom_genes = gene_contrib.tail(3)
for _, row in pd.concat([top_genes, bottom_genes]).iterrows():
    ax.annotate(row["symbol"], (row["spec_CGE"], row["delta_VIP"]),
                fontsize=8, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points')

ax.axhline(0, color='grey', lw=0.8, ls='--')
ax.set_xlabel("Mean CGE specificity", fontsize=11)
ax.set_ylabel("VIP+ − VIP− specificity", fontsize=11)
ax.set_title("22q gene CGE specificity vs VIP selectivity", fontsize=12)

rho, pval = spearmanr(gene_contrib["spec_CGE"], gene_contrib["delta_VIP"])
ax.text(0.05, 0.95, f"Spearman ρ = {rho:.2f}, P = {pval:.2g}",
        transform=ax.transAxes, fontsize=9, va='top')

legend_elements = [
    plt.scatter([], [], c='#D7191C', s=50, label='1.5 Mb deletion'),
    plt.scatter([], [], c='#4575B4', s=50, label='3 Mb only'),
]
ax.legend(frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_CGEspec_vs_VIPdelta.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %% [markdown]
# ---
# ## 2. Comparison with other neuropsychiatric CNVs
#
# Test whether CGE/VIP bias is unique to 22q or shared across psychiatric CNVs.
# We define gene lists for major CNVs from the literature and compute bias profiles.

# %%
# --- Define CNV gene lists (Entrez IDs) ---
# Sources: UCSC Genome Browser, OMIM, Marshall et al. 2017, Kendall et al. 2017
# Gene lists curated from genomic coordinates of recurrent CNV intervals

# 16p11.2 (~600 kb, BP4-BP5): 29 genes
# chr16:29.6-30.2 Mb (hg38)
CNV_16p11_symbols = [
    "SPN", "QPRT", "C16orf54", "ZG16", "KIF22", "MAZ", "PRRT2", "PAGR1",
    "MVP", "CDIPT", "SEZ6L2", "ASPHD1", "KCTD13", "TMEM219", "TAOK2",
    "HIRIP3", "INO80E", "DOC2A", "C16orf92", "FAM57B", "ALDOA", "PPP4C",
    "TBX6", "YPEL3", "GDPD3", "MAPK3", "CORO1A", "BOLA2", "SLX1A"
]

# 15q13.3 (~1.5 Mb, BP4-BP5): ~7 core genes
# chr15:30.9-32.5 Mb (hg38)
CNV_15q13_symbols = [
    "FAN1", "MTMR10", "TRPM1", "KLF13", "OTUD7A", "CHRNA7"
]

# 1q21.1 (~1.35 Mb): recurrent CNV
# chr1:146.5-147.9 Mb (hg38)
CNV_1q21_symbols = [
    "GJA5", "GJA8", "GPR89A", "BCL9", "ACP6", "PDZK1", "CHD1L",
    "FMO5", "PRKAB2", "NBPF11"
]

# 3q29 (~1.6 Mb): high SCZ risk CNV
# chr3:195.7-197.3 Mb (hg38)
CNV_3q29_symbols = [
    "TFRC", "ZDHHC19", "SLC51A", "TM4SF19", "UBXN7", "RNF168",
    "WDR53", "FBXO45", "PIGX", "PAK2", "SENP5", "NCBP2",
    "DLG1", "BDH1", "PCYT1A"
]

# 15q11.2 (BP1-BP2): ~4 genes, mild risk
CNV_15q11_symbols = [
    "NIPA1", "NIPA2", "CYFIP1", "TUBGCP5"
]


def symbols_to_gw(symbols, gene_map):
    """Convert gene symbols to uniform gene weight dict {entrez: 1}."""
    gw = {}
    missing = []
    for s in symbols:
        eid = gene_map.get(s)
        if eid is not None and eid in SpecMat.index:
            gw[eid] = 1
        else:
            missing.append(s)
    if missing:
        print(f"  Missing/not in expr matrix: {missing}")
    return gw


cnv_definitions = {
    "22q11.2\n(46 genes)": X22q_GW,
    "16p11.2\n(29 genes)": symbols_to_gw(CNV_16p11_symbols, GeneSymbol2Entrez),
    "15q13.3\n(6 genes)": symbols_to_gw(CNV_15q13_symbols, GeneSymbol2Entrez),
    "1q21.1\n(10 genes)": symbols_to_gw(CNV_1q21_symbols, GeneSymbol2Entrez),
    "3q29\n(15 genes)": symbols_to_gw(CNV_3q29_symbols, GeneSymbol2Entrez),
    "15q11.2\n(4 genes)": symbols_to_gw(CNV_15q11_symbols, GeneSymbol2Entrez),
}

print("\nCNV gene counts (valid in expression matrix):")
for name, gw in cnv_definitions.items():
    print(f"  {name.replace(chr(10), ' ')}: {len(gw)} genes")

# %%
# --- Compute bias profiles for each CNV ---
cnv_biases = {}
for name, gw in cnv_definitions.items():
    if len(gw) < 3:
        print(f"Skipping {name} — too few genes ({len(gw)})")
        continue
    bias = HumanCT_AvgZ_Weighted(SpecMat, gw)
    bias = AnnotateCTDat(bias, Anno)
    cnv_biases[name] = bias

# %%
# --- Figure 2A: CGE interneuron bias comparison across CNVs ---
fig, ax = plt.subplots(figsize=(10, 5), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

cnv_names = list(cnv_biases.keys())
cge_data = []
for name in cnv_names:
    bias = cnv_biases[name]
    cge_bias = bias.loc[bias.index.isin(CGE_idx), "EFFECT"]
    cge_data.append(cge_bias.values)

bp = ax.boxplot(cge_data, patch_artist=True, widths=0.5)
colors_cnv = ['#D7191C', '#FDAE61', '#ABD9E9', '#2C7BB6', '#A6761D', '#999999']
for patch, color in zip(bp['boxes'], colors_cnv[:len(cnv_names)]):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_xticklabels(cnv_names, fontsize=9)
ax.set_ylabel("CGE interneuron mutation bias (EFFECT)", fontsize=11)
ax.set_title("CGE interneuron bias across neuropsychiatric CNVs", fontsize=13)
ax.axhline(0, color='grey', lw=0.8, ls='--')

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_CNV_comparison_CGE.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# --- Figure 2B: VIP+ vs VIP- bias for each CNV ---
fig, axes = plt.subplots(1, len(cnv_biases), figsize=(3.2 * len(cnv_biases), 4.5),
                         dpi=150, facecolor='none', sharey=True)
fig.patch.set_alpha(0)

for i, (name, bias) in enumerate(cnv_biases.items()):
    ax = axes[i]
    ax.patch.set_alpha(0)

    vip_pos_bias = bias.loc[bias.index.isin(vip_pos_idx), "EFFECT"]
    vip_neg_bias = bias.loc[bias.index.isin(vip_neg_idx), "EFFECT"]

    bp = ax.boxplot([vip_pos_bias, vip_neg_bias], patch_artist=True, widths=0.55,
                    positions=[0, 1])
    bp['boxes'][0].set_facecolor('#D7191C')
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor('#2C7BB6')
    bp['boxes'][1].set_alpha(0.5)

    # Jitter
    for vals, pos, c in [(vip_pos_bias, 0, '#D7191C'), (vip_neg_bias, 1, '#2C7BB6')]:
        jitter = np.random.normal(pos, 0.06, size=len(vals))
        ax.scatter(jitter, vals, c=c, s=12, alpha=0.5, edgecolor='none', zorder=3)

    stat, pval = mannwhitneyu(vip_pos_bias, vip_neg_bias, alternative='two-sided')
    stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'n.s.'
    ax.set_title(f"{name}\nP={pval:.2g} {stars}", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["VIP+", "VIP−"], fontsize=9)
    if i == 0:
        ax.set_ylabel("Mutation bias (EFFECT)", fontsize=10)

fig.suptitle("VIP+ vs VIP− CGE interneuron bias across psychiatric CNVs", fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_CNV_VIP_comparison.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# --- Figure 2C: Supercluster-level bias heatmap for all CNVs ---
superclusters_of_interest = [
    "CGE interneuron", "MGE interneuron", "LAMP5-LHX6 and Chandelier",
    "Medium spiny neuron", "Upper-layer intratelencephalic",
    "Deep-layer intratelencephalic", "Hippocampal CA1-3",
    "Amygdala excitatory"
]

sc_mean_bias = pd.DataFrame(index=superclusters_of_interest, columns=list(cnv_biases.keys()))
for name, bias in cnv_biases.items():
    for sc in superclusters_of_interest:
        sc_idx = Anno[Anno["Supercluster"] == sc].index.values.astype(int)
        sc_vals = bias.loc[bias.index.isin(sc_idx), "EFFECT"]
        sc_mean_bias.loc[sc, name] = sc_vals.mean()

sc_mean_bias = sc_mean_bias.astype(float)

fig, ax = plt.subplots(figsize=(10, 5), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)
sns.heatmap(sc_mean_bias, annot=True, fmt=".3f", cmap="RdBu_r", center=0,
            linewidths=0.5, ax=ax, cbar_kws={"label": "Mean mutation bias"})
ax.set_title("Supercluster-level mutation bias across neuropsychiatric CNVs", fontsize=12)
ax.set_ylabel("")
plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_CNV_heatmap.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %% [markdown]
# ---
# ## 3. 22q11.2 DEG overlap with VIP bias
#
# Do the genes *dysregulated* by 22q deletion also show VIP/CGE bias?
#
# **Key test**: after **excluding 22q cis-genes** from the DEG lists, do the
# *trans*-regulated DEGs still show VIP+ > VIP- bias? If yes, the downstream
# transcriptomic cascade of the deletion independently converges on VIP interneurons.
#
# Two independent datasets:
# - Lin et al. 2021 (Nat Comms): iPSC-derived neurons, days 0/4/28
# - Khan et al. 2020 (Nat Med): cerebral spheroids (hCS), days 25/50/75/100

# %%
# --- Load Lin et al. 2021 DEG data (Nature Comms) ---
deg_files = {
    "Day 0 (iPSC)": str(PROJ_DIR / "dat/suppl.data/41467_2022_31436_MOESM4_ESM.xlsx"),
    "Day 4 (NPC)": str(PROJ_DIR / "dat/suppl.data/41467_2022_31436_MOESM5_ESM.xlsx"),
    "Day 28 (neuron)": str(PROJ_DIR / "dat/suppl.data/41467_2022_31436_MOESM6_ESM.xlsx"),
}

FDR_CUT = 0.05
deg_data = {}
for label, path in deg_files.items():
    df = pd.read_excel(path, skiprows=1)
    df = df[df["gene_biotype"] == "protein_coding"]
    df["EntrezID"] = df["external_gene_name"].map(GeneSymbol2Entrez)
    df = df[df["EntrezID"].notnull()].copy()
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")
    df = df.sort_values("padj")
    deg_data[label] = df
    sig = df[df["padj"] < FDR_CUT]
    n_cis = sig.index.isin(X22q_genes).sum()
    print(f"{label}: {len(sig)} DEGs (FDR<{FDR_CUT}) — {len(sig[sig['log2FoldChange']>0])} up, "
          f"{len(sig[sig['log2FoldChange']<0])} down — {n_cis} are 22q cis-genes")

# %%
# --- Load Khan et al. 2020 DEG data (Nature Medicine) ---
natmed_file = str(PROJ_DIR / "dat/suppl.data/41591_2020_1043_MOESM3_ESM.xlsx")
natmed_sheets = {
    "Day 25 (hCS)": "Day25",
    "Day 50 (hCS)": "Day50",
    "Day 75 (hCS)": "Day75",
    "Day 100 (hCS)": "Day100",
}

FDR_CUT_NM = 0.1  # This dataset uses FDR<0.1 as in original notebook
natmed_data = {}
for label, sheet in natmed_sheets.items():
    df = pd.read_excel(natmed_file, sheet_name=sheet)
    df["EntrezID"] = df["Gene"].map(GeneSymbol2Entrez)
    df = df[df["EntrezID"].notnull()].copy()
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")
    df = df.sort_values("FDR")
    natmed_data[label] = df
    sig = df[df["FDR"] < FDR_CUT_NM]
    n_cis = sig.index.isin(X22q_genes).sum()
    print(f"{label}: {len(sig)} DEGs (FDR<{FDR_CUT_NM}) — "
          f"{len(sig[sig['beta']>0])} up, {len(sig[sig['beta']<0])} down — {n_cis} are 22q cis-genes")

# %%
# --- Overlap between VIP-biased 22q cis-genes and DEGs ---
# Which of the top VIP-driving 22q genes are also DEGs?
print("=== Overlap: top VIP-driving 22q cis-genes × DEGs ===\n")
top_vip_genes = gene_contrib[gene_contrib["delta_VIP"] > 0]

print("--- Lin et al. 2021 (Nat Comms) ---")
for label, deg_df in deg_data.items():
    sig = deg_df[deg_df["padj"] < FDR_CUT]
    overlap = set(top_vip_genes.index) & set(sig.index)
    if overlap:
        print(f"{label}:")
        for g in overlap:
            sym = Entrez2Symbol.get(g, str(g))
            lfc = sig.loc[g, "log2FoldChange"]
            padj = sig.loc[g, "padj"]
            dvip = gene_contrib.loc[g, "delta_VIP"]
            print(f"  {sym:>10s}  delta_VIP={dvip:+.4f}  log2FC={lfc:+.3f}  padj={padj:.2g}")
    else:
        print(f"{label}: no overlap with VIP-driving genes")
    print()

print("--- Khan et al. 2020 (Nat Med) ---")
for label, deg_df in natmed_data.items():
    sig = deg_df[deg_df["FDR"] < FDR_CUT_NM]
    overlap = set(top_vip_genes.index) & set(sig.index)
    if overlap:
        print(f"{label}:")
        for g in overlap:
            sym = Entrez2Symbol.get(g, str(g))
            lfc = sig.loc[g, "beta"]
            fdr = sig.loc[g, "FDR"]
            dvip = gene_contrib.loc[g, "delta_VIP"]
            print(f"  {sym:>10s}  delta_VIP={dvip:+.4f}  beta={lfc:+.3f}  FDR={fdr:.2g}")
    else:
        print(f"{label}: no overlap with VIP-driving genes")
    print()

# %% [markdown]
# ### 3B. Trans-effect analysis: DEGs EXCLUDING 22q cis-genes
#
# **This is the key test.** If we remove all 22q cis-genes from the DEG lists,
# do the remaining *trans*-regulated genes still show VIP+ > VIP- bias?
# A positive result means the deletion's downstream transcriptomic footprint
# independently converges on VIP interneurons — not just a circular effect
# of the deleted genes themselves being VIP-biased.

# %%
# Set of 22q cis-gene Entrez IDs to exclude
X22q_cis_set = set(X22q_genes)
print(f"Excluding {len(X22q_cis_set)} 22q cis-genes from DEG lists\n")


def plot_deg_vip_bias(all_deg_data, fdr_col, lfc_col, fdr_cut, dataset_label,
                      exclude_cis=False, save_suffix=""):
    """Plot VIP+ vs VIP- bias for DEGs, optionally excluding 22q cis-genes.

    Parameters
    ----------
    all_deg_data : dict of label -> DataFrame
    fdr_col : str, column name for FDR/padj
    lfc_col : str, column name for log2FC/beta
    fdr_cut : float
    dataset_label : str, for figure title
    exclude_cis : bool, if True remove 22q cis-genes
    save_suffix : str, appended to output filename
    """
    ncols = len(all_deg_data)
    fig, axes = plt.subplots(2, ncols, figsize=(3.5 * ncols, 8), dpi=150, facecolor='none')
    fig.patch.set_alpha(0)
    if ncols == 1:
        axes = axes.reshape(2, 1)

    results_summary = []

    for col, (label, deg_df) in enumerate(all_deg_data.items()):
        sig = deg_df[deg_df[fdr_col] < fdr_cut]
        if exclude_cis:
            n_before = len(sig)
            sig = sig[~sig.index.isin(X22q_cis_set)]
            n_removed = n_before - len(sig)
        sig_up = sig[sig[lfc_col] > 0]
        sig_down = sig[sig[lfc_col] < 0]

        for row, (sub_label, sub_df) in enumerate([
            ("Down (trans)", sig_down) if exclude_cis else ("Downregulated", sig_down),
            ("Up (trans)", sig_up) if exclude_cis else ("Upregulated", sig_up),
        ]):
            ax = axes[row, col]
            ax.patch.set_alpha(0)

            deg_gw = {eid: 1 for eid in sub_df.index if eid in SpecMat.index}
            if len(deg_gw) < 5:
                ax.text(0.5, 0.5, f"n={len(deg_gw)} genes\n(too few)", transform=ax.transAxes,
                        ha='center', va='center', fontsize=10)
                ax.set_title(f"{label}\n{sub_label} (n={len(deg_gw)})", fontsize=9)
                continue

            bias = HumanCT_AvgZ_Weighted(SpecMat, deg_gw)
            vp = bias.loc[bias.index.isin(vip_pos_idx), "EFFECT"]
            vn = bias.loc[bias.index.isin(vip_neg_idx), "EFFECT"]

            bp = ax.boxplot([vp, vn], patch_artist=True, widths=0.55, positions=[0, 1])
            bp['boxes'][0].set_facecolor('#D7191C')
            bp['boxes'][0].set_alpha(0.5)
            bp['boxes'][1].set_facecolor('#2C7BB6')
            bp['boxes'][1].set_alpha(0.5)

            for vals, pos, c in [(vp, 0, '#D7191C'), (vn, 1, '#2C7BB6')]:
                jitter = np.random.normal(pos, 0.06, size=len(vals))
                ax.scatter(jitter, vals, c=c, s=10, alpha=0.4, edgecolor='none', zorder=3)

            stat, pval = mannwhitneyu(vp, vn, alternative='two-sided')
            stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'n.s.'
            ax.set_title(f"{label}\n{sub_label} (n={len(deg_gw)})\nP={pval:.2g} {stars}", fontsize=9)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["VIP+", "VIP-"], fontsize=9)
            if col == 0:
                ax.set_ylabel("Mutation bias (EFFECT)", fontsize=10)

            results_summary.append({
                "dataset": dataset_label, "timepoint": label,
                "direction": sub_label, "n_genes": len(deg_gw),
                "P": pval, "VIP_pos_mean": vp.mean(), "VIP_neg_mean": vn.mean()
            })

    cis_tag = " (22q cis-genes EXCLUDED)" if exclude_cis else ""
    fig.suptitle(f"VIP+ vs VIP- bias: {dataset_label} DEGs{cis_tag}",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    fname = f"fig_22q_DEG_{save_suffix}.png"
    plt.savefig(str(PROJ_DIR / "results" / fname), dpi=150, transparent=True, bbox_inches='tight')
    plt.show()
    return pd.DataFrame(results_summary)


# %% [markdown]
# #### Lin et al. 2021 — all DEGs (including cis-genes)

# %%
res_lin_all = plot_deg_vip_bias(
    deg_data, fdr_col="padj", lfc_col="log2FoldChange", fdr_cut=FDR_CUT,
    dataset_label="Lin et al. 2021", exclude_cis=False,
    save_suffix="Lin_all")

# %% [markdown]
# #### Lin et al. 2021 — trans-DEGs only (22q cis-genes EXCLUDED)

# %%
res_lin_trans = plot_deg_vip_bias(
    deg_data, fdr_col="padj", lfc_col="log2FoldChange", fdr_cut=FDR_CUT,
    dataset_label="Lin et al. 2021", exclude_cis=True,
    save_suffix="Lin_trans")

# %% [markdown]
# #### Khan et al. 2020 — all DEGs (including cis-genes)

# %%
res_khan_all = plot_deg_vip_bias(
    natmed_data, fdr_col="FDR", lfc_col="beta", fdr_cut=FDR_CUT_NM,
    dataset_label="Khan et al. 2020", exclude_cis=False,
    save_suffix="Khan_all")

# %% [markdown]
# #### Khan et al. 2020 — trans-DEGs only (22q cis-genes EXCLUDED)

# %%
res_khan_trans = plot_deg_vip_bias(
    natmed_data, fdr_col="FDR", lfc_col="beta", fdr_cut=FDR_CUT_NM,
    dataset_label="Khan et al. 2020", exclude_cis=True,
    save_suffix="Khan_trans")

# %%
# --- Summary table: all vs trans comparison ---
print("=== Summary: VIP+ vs VIP- P-values (all DEGs vs trans-only) ===\n")
all_results = pd.concat([res_lin_all, res_lin_trans, res_khan_all, res_khan_trans],
                        ignore_index=True)
for _, row in all_results.iterrows():
    tag = "TRANS" if "trans" in row["direction"] else "all  "
    sig = "***" if row["P"] < 0.001 else "** " if row["P"] < 0.01 else "*  " if row["P"] < 0.05 else "n.s"
    print(f"  {row['dataset']:20s} {row['timepoint']:18s} {row['direction']:18s} "
          f"n={row['n_genes']:4d}  P={row['P']:.4g}  {sig}")

# %% [markdown]
# ### 3C. Cross-dataset comparison: CGE ranking, overlap, and interpretation
#
# The two DEG datasets come from very different model systems:
# - **Lin et al. 2021**: 2D monolayer iPSC-derived neurons (mostly excitatory)
# - **Khan et al. 2020**: 3D cerebral spheroids (mixed cell types incl. interneurons)
#
# After excluding 22q cis-genes, do the trans-DEGs overlap?
# And do they show the same CGE bias profile?

# %%
# --- Trans-DEG overlap between datasets ---
print("=" * 65)
print("Trans-DEG overlap (22q cis-genes excluded)")
print("=" * 65)

# Build trans-DEG sets for each timepoint
lin_trans_sets = {}
for label, deg_df in deg_data.items():
    sig = deg_df[deg_df["padj"] < FDR_CUT]
    sig_trans = sig[~sig.index.isin(X22q_cis_set)]
    lin_trans_sets[label] = set(sig_trans.index)

khan_trans_sets = {}
for label, deg_df in natmed_data.items():
    sig = deg_df[deg_df["FDR"] < FDR_CUT_NM]
    sig_trans = sig[~sig.index.isin(X22q_cis_set)]
    khan_trans_sets[label] = set(sig_trans.index)

for ll, lg in lin_trans_sets.items():
    for kl, kg in khan_trans_sets.items():
        ov = lg & kg
        un = lg | kg
        j = len(ov) / len(un) if un else 0
        print(f"  Lin {ll:18s} ({len(lg):3d}) vs Khan {kl:15s} ({len(kg):3d}): "
              f"overlap={len(ov):3d}  Jaccard={j:.3f}")

print(f"\n  --> Trans-DEGs are almost entirely non-overlapping between the two systems.")

# %%
# --- CGE supercluster ranking for trans-DEGs ---
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

print("=" * 65)
print("CGE interneuron ranking among 31 superclusters (trans-DEGs only)")
print("=" * 65)


def compute_sc_ranking(gene_set, label):
    """Compute bias, return CGE rank and mean bias."""
    gw = {g: 1 for g in gene_set if g in SpecMat.index}
    if len(gw) < 5:
        return len(gw), None, None
    bias = HumanCT_AvgZ_Weighted(SpecMat, gw)
    bias = AnnotateCTDat(bias, Anno)
    sc = bias.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
    rank = list(sc.index).index("CGE interneuron") + 1
    return len(gw), rank, sc["CGE interneuron"]


print(f"\n{'Dataset':<32s} {'Dir':<6s} {'n':>5s} {'CGE rank':>10s} {'CGE bias':>10s}")
print("-" * 68)

cge_ranking_rows = []
for label, deg_df in deg_data.items():
    sig = deg_df[deg_df["padj"] < FDR_CUT]
    sig_trans = sig[~sig.index.isin(X22q_cis_set)]
    for d, filt in [("Down", sig_trans["log2FoldChange"] < 0), ("Up", sig_trans["log2FoldChange"] > 0)]:
        sub = sig_trans[filt]
        n, rank, val = compute_sc_ranking(sub.index, f"Lin {label} {d}")
        if rank is not None:
            print(f"Lin {label:<22s}      {d:<6s} {n:5d} {rank:8d}/31    {val:+.4f}")
            cge_ranking_rows.append({"dataset": f"Lin {label}", "dir": d, "n": n,
                                     "CGE_rank": rank, "CGE_bias": val, "source": "Lin et al."})
        else:
            print(f"Lin {label:<22s}      {d:<6s} {n:5d}   too few")

for label, deg_df in natmed_data.items():
    sig = deg_df[deg_df["FDR"] < FDR_CUT_NM]
    sig_trans = sig[~sig.index.isin(X22q_cis_set)]
    for d, filt in [("Down", sig_trans["beta"] < 0), ("Up", sig_trans["beta"] > 0)]:
        sub = sig_trans[filt]
        n, rank, val = compute_sc_ranking(sub.index, f"Khan {label} {d}")
        if rank is not None:
            print(f"Khan {label:<21s}      {d:<6s} {n:5d} {rank:8d}/31    {val:+.4f}")
            cge_ranking_rows.append({"dataset": f"Khan {label}", "dir": d, "n": n,
                                     "CGE_rank": rank, "CGE_bias": val, "source": "Khan et al."})
        else:
            print(f"Khan {label:<21s}      {d:<6s} {n:5d}   too few")

cge_ranking_df = pd.DataFrame(cge_ranking_rows)

# %%
# --- Figure 3C: CGE ranking comparison ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150, facecolor='none')
fig.patch.set_alpha(0)

for i, (src, color, marker) in enumerate([("Lin et al.", "#2C7BB6", "o"), ("Khan et al.", "#D7191C", "s")]):
    ax = axes[i]
    ax.patch.set_alpha(0)
    sub = cge_ranking_df[cge_ranking_df["source"] == src]
    for _, row in sub.iterrows():
        fc = color if row["dir"] == "Down" else "none"
        ax.scatter(row["n"], row["CGE_rank"], c=color, s=80, marker=marker,
                   edgecolor=color, facecolor=fc, linewidth=1.5, zorder=3)
        ax.annotate(f"{row['dataset'].split('(')[0].strip()}\n{row['dir']}",
                    (row["n"], row["CGE_rank"]), fontsize=7,
                    xytext=(6, 0), textcoords='offset points', va='center')

    ax.axhline(5, color='red', ls='--', lw=0.8, alpha=0.5, label='Top 5')
    ax.axhline(15.5, color='grey', ls=':', lw=0.8, alpha=0.5, label='Median')
    ax.set_xlabel("Number of trans-DEGs", fontsize=11)
    ax.set_ylabel("CGE interneuron rank (1=highest bias)", fontsize=11)
    ax.set_title(f"{src}\n({'2D monolayer neurons' if 'Lin' in src else '3D cerebral spheroids'})",
                 fontsize=11)
    ax.set_ylim(0.5, 31.5)
    ax.invert_yaxis()
    ax.legend(fontsize=8, loc='lower right')

fig.suptitle("CGE supercluster ranking of trans-DEGs (22q cis-genes excluded)", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_DEG_CGE_ranking.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
print("""
INTERPRETATION:
  Lin et al. (2D monolayer): Directed differentiation → mostly excitatory neurons.
    Trans-DEGs are NOT CGE-biased (ranks 16-26/31, negative bias).
    Yet some subsets still show VIP+ > VIP- *within* CGE → subtle specificity.

  Khan et al. (3D spheroids): Self-organized tissue → mixed cell types incl. interneurons.
    Trans-DEGs ARE CGE-biased (ranks 3-8/31, positive bias).
    The 3D model captures interneuron-relevant transcriptomic programs.

  The two datasets have near-zero gene overlap (Jaccard < 0.01) yet BOTH
  show VIP+ > VIP- bias after excluding cis-genes. This convergence from
  completely independent gene sets and model systems strongly supports
  a genuine VIP+ interneuron vulnerability downstream of 22q deletion.
""")

# %% [markdown]
# ### 3D. Gordon et al. 2019 (Mol Psychiatry) — Df(h22q11)/+ mouse brain tissue
#
# DEGs from **real mouse brain** (cortex and hippocampus) of a 22q11.2 deletion model.
# Unlike iPSC cultures, this captures in vivo transcriptomic consequences in intact tissue.
# Mouse gene symbols are mapped to human orthologs via uppercasing (standard for 1:1 orthologs).

# %%
# --- Load Gordon et al. 22q11 DEGs (mouse brain) ---
gordon_file = str(PROJ_DIR / "dat/22q_DEGs/Gordon_et_al/41380_2019_576_MOESM3_ESM.xlsx")

gordon_data = {}
for sheet, label in [("cortex 22q11", "Cortex (mouse)"), ("hippocampus 22q11", "Hippocampus (mouse)")]:
    df = pd.read_excel(gordon_file, sheet_name=sheet)
    # Map mouse symbols to human Entrez IDs (uppercase first letter → human symbol)
    df["human_symbol"] = df["mgi_symbol"].str.upper()
    # Special cases where mouse/human symbols differ
    symbol_fixes = {"SEPT5": "SEPTIN5", "UFD1L": "UFD1"}
    df["human_symbol"] = df["human_symbol"].replace(symbol_fixes)
    df["EntrezID"] = df["human_symbol"].map(GeneSymbol2Entrez)
    df = df[df["EntrezID"].notnull()].copy()
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")

    gordon_data[label] = df
    # Use P<0.005 as threshold (as in the original paper)
    sig = df[df["P.Value"] < 0.005]
    n_cis = sig.index.isin(X22q_cis_set).sum()
    print(f"{label}: {len(df)} genes mapped, {len(sig)} DEGs (P<0.005) — "
          f"{len(sig[sig['logFC']>0])} up, {len(sig[sig['logFC']<0])} down — "
          f"{n_cis} are 22q cis-genes")

# %%
# --- VIP+ vs VIP- bias: Gordon et al., all DEGs and trans-only ---
PCUT_GORDON = 0.005  # Paper's threshold

fig, axes = plt.subplots(2, 2, figsize=(10, 9), dpi=150, facecolor='none')
fig.patch.set_alpha(0)

gordon_results = []
for col, (label, deg_df) in enumerate(gordon_data.items()):
    sig = deg_df[deg_df["P.Value"] < PCUT_GORDON]

    for row, (exclude, row_label) in enumerate([(False, "All DEGs"), (True, "Trans-only (excl 22q cis)")]):
        ax = axes[row, col]
        ax.patch.set_alpha(0)

        sub = sig[~sig.index.isin(X22q_cis_set)] if exclude else sig
        # Use all DEGs (not split by direction) since sample sizes are small
        deg_gw = {eid: 1 for eid in sub.index if eid in SpecMat.index}

        if len(deg_gw) < 5:
            ax.text(0.5, 0.5, f"n={len(deg_gw)} genes\n(too few)", transform=ax.transAxes,
                    ha='center', va='center', fontsize=10)
            ax.set_title(f"{label}\n{row_label} (n={len(deg_gw)})", fontsize=10)
            continue

        bias = HumanCT_AvgZ_Weighted(SpecMat, deg_gw)
        vp = bias.loc[bias.index.isin(vip_pos_idx), "EFFECT"]
        vn = bias.loc[bias.index.isin(vip_neg_idx), "EFFECT"]

        bp = ax.boxplot([vp, vn], patch_artist=True, widths=0.55, positions=[0, 1])
        bp['boxes'][0].set_facecolor('#D7191C')
        bp['boxes'][0].set_alpha(0.5)
        bp['boxes'][1].set_facecolor('#2C7BB6')
        bp['boxes'][1].set_alpha(0.5)
        for vals, pos, c in [(vp, 0, '#D7191C'), (vn, 1, '#2C7BB6')]:
            jitter = np.random.normal(pos, 0.06, size=len(vals))
            ax.scatter(jitter, vals, c=c, s=12, alpha=0.5, edgecolor='none', zorder=3)

        stat, pval = mannwhitneyu(vp, vn, alternative='two-sided')
        stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'n.s.'
        ax.set_title(f"{label}\n{row_label} (n={len(deg_gw)})\nP={pval:.2g} {stars}", fontsize=10)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["VIP+", "VIP-"], fontsize=10)
        if col == 0:
            ax.set_ylabel("Mutation bias (EFFECT)", fontsize=10)

        # Also compute CGE ranking
        bias_ann = AnnotateCTDat(bias, Anno)
        sc = bias_ann.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
        cge_rank = list(sc.index).index("CGE interneuron") + 1

        gordon_results.append({
            "tissue": label, "filter": row_label, "n_genes": len(deg_gw),
            "P_VIP": pval, "CGE_rank": cge_rank, "CGE_bias": sc["CGE interneuron"],
        })

fig.suptitle("VIP+ vs VIP- bias: Gordon et al. 2019 — Df(h22q11)/+ mouse brain\n(P < 0.005)",
             fontsize=12, y=1.03)
plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_DEG_Gordon_VIP.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# --- Summary ---
print("=== Gordon et al. 2019 — mouse brain 22q11 DEGs ===\n")
print(f"{'Tissue':<22s} {'Filter':<30s} {'n':>5s} {'VIP P':>10s} {'CGE rank':>10s} {'CGE bias':>10s}")
print("-" * 90)
for r in gordon_results:
    sig = '***' if r['P_VIP'] < 0.001 else '** ' if r['P_VIP'] < 0.01 else '*  ' if r['P_VIP'] < 0.05 else 'n.s'
    print(f"{r['tissue']:<22s} {r['filter']:<30s} {r['n_genes']:5d} {r['P_VIP']:10.4g} {sig}"
          f" {r['CGE_rank']:7d}/31  {r['CGE_bias']:+10.4f}")

# %% [markdown]
# ---
# ## 4. Jackknife stability of VIP+ > VIP- bias
#
# Remove each 22q gene one at a time and recompute the VIP+ vs VIP- bias difference.
# This reveals which genes are critical for maintaining the signal (and which, like TBX1,
# when removed cause the signal to weaken).

# %%
# --- Jackknife: leave-one-out ---
# Full-set VIP bias difference
full_bias = HumanCT_AvgZ_Weighted(SpecMat, X22q_GW)
full_vip_pos_mean = full_bias.loc[full_bias.index.isin(vip_pos_idx), "EFFECT"].mean()
full_vip_neg_mean = full_bias.loc[full_bias.index.isin(vip_neg_idx), "EFFECT"].mean()
full_delta = full_vip_pos_mean - full_vip_neg_mean

# Also compute MWU p-value for full set
full_vp = full_bias.loc[full_bias.index.isin(vip_pos_idx), "EFFECT"]
full_vn = full_bias.loc[full_bias.index.isin(vip_neg_idx), "EFFECT"]
_, full_pval = mannwhitneyu(full_vp, full_vn, alternative='two-sided')

print(f"Full 22q gene set: VIP+ mean = {full_vip_pos_mean:.4f}, VIP- mean = {full_vip_neg_mean:.4f}")
print(f"  Delta (VIP+ - VIP-) = {full_delta:.4f}, MWU P = {full_pval:.4g}")

# %%
# Leave-one-out
jackknife_results = []
for gene in valid_22q:
    loo_gw = {g: w for g, w in X22q_GW.items() if g != gene}
    loo_bias = HumanCT_AvgZ_Weighted(SpecMat, loo_gw)

    vp_mean = loo_bias.loc[loo_bias.index.isin(vip_pos_idx), "EFFECT"].mean()
    vn_mean = loo_bias.loc[loo_bias.index.isin(vip_neg_idx), "EFFECT"].mean()
    delta = vp_mean - vn_mean

    vp = loo_bias.loc[loo_bias.index.isin(vip_pos_idx), "EFFECT"]
    vn = loo_bias.loc[loo_bias.index.isin(vip_neg_idx), "EFFECT"]
    _, pval = mannwhitneyu(vp, vn, alternative='two-sided')

    jackknife_results.append({
        "entrez": gene,
        "symbol": Entrez2Symbol.get(gene, str(gene)),
        "delta_VIP": delta,
        "pval": pval,
        "delta_change": delta - full_delta,
        "in_1.5Mb": gene in small_del_genes,
    })

jk_df = pd.DataFrame(jackknife_results)
jk_df = jk_df.sort_values("delta_change", ascending=True)  # most destabilizing first

# %%
# --- Figure 4A: Jackknife bar chart ---
fig, ax = plt.subplots(figsize=(14, 5), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

# Color by whether removal hurts (red) or helps (blue) the VIP signal
colors = ['#D7191C' if d < 0 else '#2C7BB6' for d in jk_df["delta_change"]]
bars = ax.bar(range(len(jk_df)), jk_df["delta_change"], color=colors, alpha=0.8)

# Mark 1.5Mb genes
for i, (_, row) in enumerate(jk_df.iterrows()):
    if row["in_1.5Mb"]:
        bars[i].set_edgecolor('black')
        bars[i].set_linewidth(1.5)

ax.set_xticks(range(len(jk_df)))
ax.set_xticklabels(jk_df["symbol"], rotation=90, fontsize=8)
ax.set_ylabel("Δ(VIP+ − VIP−) vs full set", fontsize=11)
ax.set_title("Jackknife: effect of removing each 22q gene on VIP+ > VIP− bias", fontsize=12)
ax.axhline(0, color='black', lw=0.8)

legend_elements = [
    Patch(facecolor='#D7191C', alpha=0.8, label='Removal weakens VIP signal'),
    Patch(facecolor='#2C7BB6', alpha=0.8, label='Removal strengthens VIP signal'),
    Patch(facecolor='grey', edgecolor='black', linewidth=1.5, label='In 1.5 Mb deletion'),
]
ax.legend(handles=legend_elements, frameon=True, fontsize=9, loc='lower left')

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_jackknife_VIP.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# --- Figure 4B: Jackknife p-value stability ---
fig, ax = plt.subplots(figsize=(14, 4), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

jk_sorted = jk_df.sort_values("pval", ascending=False)
colors_p = ['#D7191C' if p < 0.05 else '#FDAE61' if p < 0.1 else '#999999'
            for p in jk_sorted["pval"]]
ax.bar(range(len(jk_sorted)), -np.log10(jk_sorted["pval"]), color=colors_p, alpha=0.8)
ax.axhline(-np.log10(0.05), color='red', ls='--', lw=1, label='P = 0.05')
ax.axhline(-np.log10(full_pval), color='black', ls=':', lw=1.5, label=f'Full set P = {full_pval:.3g}')

ax.set_xticks(range(len(jk_sorted)))
ax.set_xticklabels(jk_sorted["symbol"], rotation=90, fontsize=8)
ax.set_ylabel(r"$-\log_{10}$(P) VIP+ vs VIP$-$", fontsize=11)
ax.set_title("Jackknife: VIP+ vs VIP- significance after removing each gene", fontsize=12)
ax.legend(frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_jackknife_pval.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# Print genes whose removal most destabilizes the signal
print("=== Genes whose removal WEAKENS VIP+ > VIP- signal (most destabilizing) ===")
for _, row in jk_df.head(10).iterrows():
    flag = " *1.5Mb" if row["in_1.5Mb"] else ""
    print(f"  {row['symbol']:>10s}  Δ = {row['delta_change']:+.5f}  "
          f"LOO delta = {row['delta_VIP']:.5f}  P = {row['pval']:.3g}{flag}")

print(f"\n  Full set:    delta = {full_delta:.5f}  P = {full_pval:.3g}")

print("\n=== Genes whose removal STRENGTHENS VIP signal ===")
for _, row in jk_df.tail(5).iterrows():
    flag = " *1.5Mb" if row["in_1.5Mb"] else ""
    print(f"  {row['symbol']:>10s}  Δ = {row['delta_change']:+.5f}  "
          f"LOO delta = {row['delta_VIP']:.5f}  P = {row['pval']:.3g}{flag}")

# %% [markdown]
# ---
# ## Summary
#
# ### Key findings:
# 1. **Per-gene contribution**: A small number of 22q genes drive the VIP+ bias.
#    Identifying these "driver" genes connects the CNV-level signal to specific
#    molecular mechanisms.
#
# 2. **CNV comparison**: Tests whether VIP/CGE bias is specific to 22q or a
#    general property of neuropsychiatric CNVs.
#
# 3. **DEG analysis**: Whether the *downstream* transcriptomic footprint of the
#    22q deletion also shows VIP-preferential expression.
#
# 4. **Jackknife**: The VIP signal is robust to most single-gene removals,
#    but certain genes (e.g., TBX1) are critical — removing them weakens the signal.
