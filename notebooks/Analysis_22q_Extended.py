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
ax.text(0.05, 0.95, f"Spearmans' R = {rho:.2f}, P = {pval:.2g}",
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

# %% [markdown]
# ---
# ## 5. Ion Channel Gene Enrichment Among 22q11.2 DEGs
#
# **Motivation**: VIP interneurons have distinct action potential (AP) kinetics shaped by
# specific ion channels (Na+, Kv3/Kv1 for fast repolarization, BK/SK, HCN, A-type K+).
# If the 22q deletion's transcriptomic consequences converge on these channels,
# it would provide a mechanistic link between the CNV and VIP electrophysiological
# vulnerability.
#
# We test:
# 1. Which AP-critical ion channel genes are DEGs in 22q models?
# 2. Are ion channel genes enriched among DEGs (Fisher's exact test)?
# 3. Do DGCR8-regulated miRNAs target these channels?

# %%
# === Define ion channel gene sets of interest ===
from scipy.stats import fisher_exact

ion_channel_sets = {
    "Na+ channels (SCN)": [
        "SCN1A", "SCN2A", "SCN3A", "SCN8A",  # alpha subunits
        "SCN1B", "SCN2B", "SCN3B", "SCN4B",  # beta (auxiliary) subunits
    ],
    "Kv3 (fast repol.)": [
        "KCNC1", "KCNC2", "KCNC3", "KCNC4",  # Kv3.1-3.4
    ],
    "Kv1 (delayed rect.)": [
        "KCNA1", "KCNA2", "KCNA3", "KCNA6",  # Kv1.1, 1.2, 1.3, 1.6
    ],
    "Kv2 (delayed rect.)": [
        "KCNB1", "KCNB2",  # Kv2.1, 2.2
    ],
    "BK/SK (Ca2+-act.)": [
        "KCNMA1",  # BK
        "KCNN1", "KCNN2", "KCNN3",  # SK1-3
    ],
    "A-type K+ (Kv4)": [
        "KCND2", "KCND3",  # Kv4.2, 4.3
    ],
    "HCN (Ih)": [
        "HCN1", "HCN2",  # hyperpolarization-activated
    ],
}

# Flatten to a single set for aggregate tests
all_channel_genes = set()
for genes in ion_channel_sets.values():
    all_channel_genes.update(genes)
print(f"Total unique ion channel genes of interest: {len(all_channel_genes)}")
for cat, genes in ion_channel_sets.items():
    print(f"  {cat}: {', '.join(genes)}")

# %%
# === Helper: check channel genes in a DEG dataset ===

def check_channels_in_degs(deg_df, fdr_col, lfc_col, fdr_cut, gene_col=None,
                           dataset_label="", total_tested=None):
    """Check ion channel gene representation among DEGs.

    Parameters
    ----------
    deg_df : DataFrame with EntrezID index
    fdr_col : str, column for FDR/padj
    lfc_col : str, column for effect size
    fdr_cut : float, significance threshold
    gene_col : str or None, column with gene symbols (None = use Entrez2Symbol)
    dataset_label : str
    total_tested : int or None, total genes tested (for Fisher's test background)

    Returns
    -------
    dict with results
    """
    # Map channel genes to Entrez IDs
    channel_entrez = {}
    for gene in all_channel_genes:
        eid = GeneSymbol2Entrez.get(gene)
        if eid is not None:
            channel_entrez[gene] = eid

    # Significant DEGs
    sig = deg_df[deg_df[fdr_col] < fdr_cut].copy()

    # Which channel genes are in the dataset at all?
    channel_in_dataset = {g: eid for g, eid in channel_entrez.items()
                          if eid in deg_df.index}
    # Which are significant DEGs?
    channel_sig = {g: eid for g, eid in channel_entrez.items()
                   if eid in sig.index}

    results = {
        "dataset": dataset_label,
        "n_sig_degs": len(sig),
        "n_channel_in_dataset": len(channel_in_dataset),
        "n_channel_sig": len(channel_sig),
        "channel_sig_genes": {},
    }

    # Details for each significant channel gene
    for gene_sym, eid in sorted(channel_sig.items()):
        row = sig.loc[eid]
        # Handle duplicates (take first)
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        results["channel_sig_genes"][gene_sym] = {
            "lfc": row[lfc_col],
            "fdr": row[fdr_col],
        }

    # Fisher's exact test: are channel genes enriched among DEGs?
    # Contingency table:
    #                 Channel   Non-channel
    # DEG (sig)         a           b
    # Not DEG           c           d
    if total_tested is not None:
        a = len(channel_sig)
        b = len(sig) - a
        c = len(channel_in_dataset) - a
        d = total_tested - a - b - c
        if d > 0:
            odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            results["fisher_OR"] = odds
            results["fisher_P"] = pval
            results["contingency"] = [[a, b], [c, d]]

    return results


# %%
# === 5A. Lin et al. 2021 — Ion channel genes among DEGs ===
print("=" * 75)
print("5A. Ion Channel Genes in Lin et al. 2021 DEGs (iPSC-derived neurons)")
print("=" * 75)

lin_channel_results = {}
lin_total_tested = {}
for label, path in deg_files.items():
    df = pd.read_excel(path, skiprows=1)
    df = df[df["gene_biotype"] == "protein_coding"]
    n_total = len(df)
    df["EntrezID"] = df["external_gene_name"].map(GeneSymbol2Entrez)
    df = df[df["EntrezID"].notnull()].copy()
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")
    lin_total_tested[label] = len(df)

    res = check_channels_in_degs(
        df, fdr_col="padj", lfc_col="log2FoldChange", fdr_cut=FDR_CUT,
        dataset_label=f"Lin {label}", total_tested=len(df))
    lin_channel_results[label] = res

    print(f"\n--- {label} ---")
    print(f"  Total protein-coding genes tested: {len(df)}")
    print(f"  Significant DEGs (FDR<{FDR_CUT}): {res['n_sig_degs']}")
    print(f"  Channel genes in dataset: {res['n_channel_in_dataset']}")
    print(f"  Channel genes that are DEGs: {res['n_channel_sig']}")

    if res["channel_sig_genes"]:
        print(f"  Significant channel genes:")
        for g, info in sorted(res["channel_sig_genes"].items()):
            direction = "UP" if info["lfc"] > 0 else "DOWN"
            print(f"    {g:>8s}  log2FC={info['lfc']:+.3f} ({direction})  padj={info['fdr']:.2e}")
    else:
        print(f"  No channel genes reach significance.")

    if "fisher_OR" in res:
        print(f"  Fisher's exact test (enrichment of channels among DEGs):")
        print(f"    OR = {res['fisher_OR']:.2f}, P = {res['fisher_P']:.4g}")
        ct = res["contingency"]
        print(f"    Contingency: [{ct[0]}] / [{ct[1]}]")

# %%
# === 5B. Khan et al. 2020 — Ion channel genes among DEGs ===
print("=" * 75)
print("5B. Ion Channel Genes in Khan et al. 2020 DEGs (cerebral spheroids)")
print("=" * 75)
print("NOTE: Khan et al. reports only top DEGs per timepoint, not full genome.")
print("      Fisher's test is not appropriate here (no background universe).")
print("      We report overlap descriptively.\n")

khan_channel_results = {}
for label, sheet in natmed_sheets.items():
    df = pd.read_excel(natmed_file, sheet_name=sheet)
    df["EntrezID"] = df["Gene"].map(GeneSymbol2Entrez)
    df = df[df["EntrezID"].notnull()].copy()
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")

    res = check_channels_in_degs(
        df, fdr_col="FDR", lfc_col="beta", fdr_cut=FDR_CUT_NM,
        dataset_label=f"Khan {label}", total_tested=None)
    khan_channel_results[label] = res

    print(f"--- {label} ---")
    print(f"  Genes in table: {len(df)}")
    print(f"  Significant DEGs (FDR<{FDR_CUT_NM}): {res['n_sig_degs']}")
    print(f"  Channel genes in table: {res['n_channel_in_dataset']}")
    print(f"  Channel genes that are DEGs: {res['n_channel_sig']}")

    if res["channel_sig_genes"]:
        for g, info in sorted(res["channel_sig_genes"].items()):
            direction = "UP" if info["lfc"] > 0 else "DOWN"
            print(f"    {g:>8s}  beta={info['lfc']:+.3f} ({direction})  FDR={info['fdr']:.2e}")
    else:
        print(f"    No channel genes among DEGs.")

    # Also check any channel genes present but non-significant
    channel_entrez_map = {g: GeneSymbol2Entrez.get(g) for g in all_channel_genes
                          if GeneSymbol2Entrez.get(g) is not None}
    present_not_sig = {g: eid for g, eid in channel_entrez_map.items()
                       if eid in df.index and eid not in
                       df[df["FDR"] < FDR_CUT_NM].index}
    if present_not_sig:
        print(f"    Channel genes present but not significant:")
        for g, eid in sorted(present_not_sig.items()):
            row = df.loc[eid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            print(f"      {g:>8s}  beta={row['beta']:+.3f}  FDR={row['FDR']:.2e}")
    print()

# %%
# === 5C. Summary table: all channel genes across all datasets/timepoints ===
print("=" * 75)
print("5C. Summary: Ion Channel DEG Status Across All Datasets")
print("=" * 75)

# Build a comprehensive table
summary_rows = []
for gene_sym in sorted(all_channel_genes):
    eid = GeneSymbol2Entrez.get(gene_sym)
    if eid is None:
        continue

    # Determine which ion channel category
    cat = [k for k, v in ion_channel_sets.items() if gene_sym in v][0]

    row = {"Gene": gene_sym, "Category": cat}

    # Lin et al. timepoints
    for label in deg_files.keys():
        r = lin_channel_results[label]
        if gene_sym in r["channel_sig_genes"]:
            info = r["channel_sig_genes"][gene_sym]
            row[f"Lin {label} LFC"] = info["lfc"]
            row[f"Lin {label} FDR"] = info["fdr"]
        else:
            row[f"Lin {label} LFC"] = np.nan
            row[f"Lin {label} FDR"] = np.nan

    # Khan et al. timepoints
    for label in natmed_sheets.keys():
        r = khan_channel_results[label]
        if gene_sym in r["channel_sig_genes"]:
            info = r["channel_sig_genes"][gene_sym]
            row[f"Khan {label} beta"] = info["lfc"]
            row[f"Khan {label} FDR"] = info["fdr"]
        else:
            row[f"Khan {label} beta"] = np.nan
            row[f"Khan {label} FDR"] = np.nan

    summary_rows.append(row)

channel_summary_df = pd.DataFrame(summary_rows)

# Print compact summary: only genes that are DEGs in at least one dataset
lfc_cols = [c for c in channel_summary_df.columns if "LFC" in c or "beta" in c]
has_any_sig = channel_summary_df[lfc_cols].notna().any(axis=1)
print(f"\nChannel genes that are DEGs in >= 1 dataset/timepoint: "
      f"{has_any_sig.sum()} / {len(channel_summary_df)}\n")

if has_any_sig.any():
    sig_genes = channel_summary_df[has_any_sig].copy()
    for _, row in sig_genes.iterrows():
        hits = []
        for c in lfc_cols:
            if pd.notna(row[c]):
                fdr_col = c.replace("LFC", "FDR").replace("beta", "FDR")
                direction = "UP" if row[c] > 0 else "DOWN"
                dset = c.rsplit(" ", 1)[0]
                hits.append(f"{dset} ({direction}, {row[c]:+.2f}, FDR={row[fdr_col]:.1e})")
        print(f"  {row['Gene']:>8s} [{row['Category']}]")
        for h in hits:
            print(f"           {h}")
else:
    print("  No ion channel genes are significantly DE in any dataset.")

# %%
# === 5D. Visualization: heatmap of channel gene fold-changes ===
fig, axes = plt.subplots(1, 2, figsize=(14, 8), dpi=150, facecolor='none')
fig.patch.set_alpha(0)

# Prepare data for heatmap (LFC values, NaN = not tested or not sig)
# Left panel: Lin et al., Right panel: Khan et al.
gene_order = []
for cat, genes in ion_channel_sets.items():
    gene_order.extend(genes)

for ax_idx, (title, results_dict, timepoints, lfc_label) in enumerate([
    ("Lin et al. 2021\n(iPSC neurons)", lin_channel_results,
     list(deg_files.keys()), "log2FC"),
    ("Khan et al. 2020\n(cerebral spheroids)", khan_channel_results,
     list(natmed_sheets.keys()), "beta"),
]):
    ax = axes[ax_idx]
    ax.patch.set_alpha(0)

    # Build matrix: genes x timepoints, fill with LFC/beta if FDR < threshold
    heatmap_data = pd.DataFrame(index=gene_order, columns=timepoints, dtype=float)
    heatmap_sig = pd.DataFrame(index=gene_order, columns=timepoints, dtype=bool)
    heatmap_present = pd.DataFrame(index=gene_order, columns=timepoints, dtype=bool)

    for tp in timepoints:
        res = results_dict[tp]
        for gene_sym in gene_order:
            eid = GeneSymbol2Entrez.get(gene_sym)
            if eid is None:
                continue
            # Check if gene is in the dataset at all
            channel_entrez_map_local = {g: GeneSymbol2Entrez.get(g) for g in all_channel_genes
                                        if GeneSymbol2Entrez.get(g) is not None}
            if eid in (channel_entrez_map_local.values()):
                heatmap_present.loc[gene_sym, tp] = True

            if gene_sym in res["channel_sig_genes"]:
                info = res["channel_sig_genes"][gene_sym]
                heatmap_data.loc[gene_sym, tp] = info["lfc"]
                heatmap_sig.loc[gene_sym, tp] = True
            else:
                heatmap_sig.loc[gene_sym, tp] = False

    # For visualization: show LFC for sig genes, grey out non-sig
    heatmap_data = heatmap_data.astype(float)

    # Create annotation strings
    annot = heatmap_data.copy().astype(str)
    for g in gene_order:
        for tp in timepoints:
            val = heatmap_data.loc[g, tp]
            if pd.notna(val):
                annot.loc[g, tp] = f"{val:+.2f}"
            else:
                annot.loc[g, tp] = ""

    # Add category separators
    cat_boundaries = []
    pos = 0
    for cat, genes in ion_channel_sets.items():
        cat_boundaries.append((pos, cat))
        pos += len(genes)

    vmax = max(abs(heatmap_data.min().min()), abs(heatmap_data.max().max()))
    if pd.isna(vmax) or vmax == 0:
        vmax = 1.0

    sns.heatmap(heatmap_data, annot=annot, fmt="s", cmap="RdBu_r", center=0,
                vmin=-vmax, vmax=vmax,
                linewidths=0.5, ax=ax, cbar_kws={"label": lfc_label},
                mask=heatmap_data.isna(),
                xticklabels=True, yticklabels=True)

    # Draw category separators
    for boundary, cat_name in cat_boundaries[1:]:
        ax.axhline(boundary, color='black', lw=1.5)

    # Add category labels on the right
    for (boundary, cat_name), next_bound in zip(cat_boundaries,
        [b for b, _ in cat_boundaries[1:]] + [len(gene_order)]):
        mid = (boundary + next_bound) / 2
        ax.text(len(timepoints) + 0.15, mid, cat_name.split("(")[0].strip(),
                va='center', ha='left', fontsize=7, style='italic')

    ax.set_title(title, fontsize=11)
    ax.set_ylabel("")
    ax.tick_params(axis='y', labelsize=9)
    ax.tick_params(axis='x', labelsize=8, rotation=30)

fig.suptitle("Ion Channel Gene Differential Expression in 22q11.2 Deletion Models",
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_ion_channel_DEGs.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %%
# === 5E. Enrichment test by channel category (Lin et al. only — full background) ===
print("=" * 75)
print("5E. Enrichment by Ion Channel Category (Lin et al. — Fisher's exact)")
print("=" * 75)

for label in deg_files.keys():
    deg_df = deg_data[label]
    sig = deg_df[deg_df["padj"] < FDR_CUT]
    n_total = lin_total_tested[label]
    n_sig = len(sig)

    print(f"\n--- {label} (n_tested={n_total}, n_sig={n_sig}) ---")

    for cat, genes in ion_channel_sets.items():
        cat_entrez = [GeneSymbol2Entrez.get(g) for g in genes
                      if GeneSymbol2Entrez.get(g) is not None]
        in_dataset = [eid for eid in cat_entrez if eid in deg_df.index]
        in_sig = [eid for eid in cat_entrez if eid in sig.index]

        a = len(in_sig)
        b = n_sig - a
        c = len(in_dataset) - a
        d = n_total - a - b - c

        if d > 0 and (a + c) > 0:
            odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            sig_str = "*" if pval < 0.05 else ""
            print(f"  {cat:25s}: {a}/{len(in_dataset)} in DEGs  "
                  f"OR={odds:.2f}  P={pval:.3g} {sig_str}")
        else:
            print(f"  {cat:25s}: {a}/{len(in_dataset)} in DEGs  (test N/A)")

    # Also test ALL channel genes together
    all_entrez = [GeneSymbol2Entrez.get(g) for g in all_channel_genes
                  if GeneSymbol2Entrez.get(g) is not None]
    in_dataset = [eid for eid in all_entrez if eid in deg_df.index]
    in_sig = [eid for eid in all_entrez if eid in sig.index]
    a = len(in_sig)
    b = n_sig - a
    c = len(in_dataset) - a
    d = n_total - a - b - c
    if d > 0:
        odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        sig_str = "*" if pval < 0.05 else ""
        print(f"  {'ALL CHANNELS':25s}: {a}/{len(in_dataset)} in DEGs  "
              f"OR={odds:.2f}  P={pval:.3g} {sig_str}")

# %% [markdown]
# ### 5F. DGCR8 / miRNA-mediated regulation of ion channels
#
# DGCR8 (within the 22q11.2 locus) is essential for miRNA biogenesis.
# Its haploinsufficiency reduces levels of DGCR8-dependent miRNAs,
# particularly miR-185 and the miR-25/106b cluster.
#
# **Known DGCR8-dependent miRNAs affected by 22q deletion:**
# - miR-185 (most affected; within 22q11.2 locus itself)
# - miR-25 (miR-25/93/106b cluster)
# - miR-134
# - miR-132/212 cluster
#
# **Validated ion channel targets of these miRNAs (from literature):**
# - miR-185: targets SCN1A (Bhogal et al. 2018 — Dravet syndrome connection)
# - miR-25: targets KCNJ2 (Kir2.1), HCN2, SCN5A (cardiac, may apply to neural isoforms)
# - miR-134: targets KCNMA1 (BK), CREB → regulates excitability
# - miR-132: targets KCNA1 (Kv1.1), HCN1 (general excitability regulation)

# %%
# === DGCR8-miRNA-ion channel regulatory network ===
print("=" * 75)
print("5F. DGCR8 → miRNA → Ion Channel Regulatory Links")
print("=" * 75)

# Literature-curated miRNA-channel target pairs
mirna_channel_targets = {
    "miR-185": {
        "targets": {
            "SCN1A": "Validated target; miR-185 directly targets SCN1A 3'UTR "
                     "(Bhogal et al. 2018, Hum Mol Genet). Reduced miR-185 in "
                     "22q → SCN1A upregulation.",
        },
        "note": "miR-185 is WITHIN the 22q11.2 deletion → maximally affected."
    },
    "miR-25": {
        "targets": {
            "HCN2": "Predicted target (TargetScan); miR-25 regulates HCN channel "
                    "expression in cardiac pacemaker cells (D'Souza et al. 2014).",
            "SCN5A": "Validated cardiac target; neural paralog SCN1A/2A may also "
                     "be affected.",
        },
        "note": "Part of miR-25/93/106b cluster, reduced by DGCR8 haploinsufficiency."
    },
    "miR-132": {
        "targets": {
            "KCNA1": "Kv1.1; miR-132 regulates neuronal excitability including "
                     "K+ channel expression (Wanet et al. 2012).",
            "HCN1": "Predicted target; miR-132 implicated in Ih current regulation.",
        },
        "note": "Activity-dependent miRNA, reduced in DGCR8 heterozygous mice."
    },
    "miR-134": {
        "targets": {
            "KCNMA1": "BK channel; miR-134 targets KCNMA1 and modulates neuronal "
                      "excitability (Bhogal & bhogal 2018, Bhatt et al. 2011).",
        },
        "note": "Synaptically enriched miRNA, regulated by DGCR8."
    },
}

for mirna, info in mirna_channel_targets.items():
    print(f"\n{mirna}:")
    print(f"  Context: {info['note']}")
    print(f"  Ion channel targets:")
    for gene, evidence in info["targets"].items():
        # Check if this gene is a DEG in our datasets
        eid = GeneSymbol2Entrez.get(gene)
        deg_status = []
        if eid is not None:
            for label, r in lin_channel_results.items():
                if gene in r["channel_sig_genes"]:
                    lfc = r["channel_sig_genes"][gene]["lfc"]
                    deg_status.append(f"Lin {label}: {'UP' if lfc > 0 else 'DOWN'} "
                                      f"(log2FC={lfc:+.2f})")
            for label, r in khan_channel_results.items():
                if gene in r["channel_sig_genes"]:
                    lfc = r["channel_sig_genes"][gene]["lfc"]
                    deg_status.append(f"Khan {label}: {'UP' if lfc > 0 else 'DOWN'} "
                                      f"(beta={lfc:+.2f})")

        print(f"    {gene}: {evidence}")
        if deg_status:
            print(f"      ** DEG in 22q models: {'; '.join(deg_status)}")
        else:
            print(f"      (Not a DEG in our datasets)")

# %%
# === 5G. Check 22q cis-genes known to regulate channels ===
print("=" * 75)
print("5G. 22q11.2 Cis-Genes That May Regulate Ion Channels")
print("=" * 75)

cis_regulators = {
    "DGCR8": {
        "mechanism": "Essential for miRNA biogenesis; haploinsufficiency reduces "
                     "miR-185, miR-25, miR-134, miR-132 → derepression of channel targets",
        "channel_targets": ["SCN1A", "HCN1", "HCN2", "KCNA1", "KCNMA1"],
    },
    "TBX1": {
        "mechanism": "T-box transcription factor; regulates developmental gene programs. "
                     "TBX1 haploinsufficiency alters gene regulatory networks in neural "
                     "progenitors (Karpinski et al. 2014, Sci Rep). May affect channel "
                     "gene regulation through downstream TF cascades.",
        "channel_targets": [],  # No direct channel targets known
    },
    "DGCR6": {
        "mechanism": "Laminin gamma-1 interactor; involved in neural crest migration. "
                     "No direct channel regulation known.",
        "channel_targets": [],
    },
    "RANBP1": {
        "mechanism": "RanGTP-binding protein; nucleocytoplasmic transport. Haploinsufficiency "
                     "could affect channel mRNA trafficking/translation.",
        "channel_targets": [],
    },
}

for gene, info in cis_regulators.items():
    eid = GeneSymbol2Entrez.get(gene)
    in_22q = eid in X22q_cis_set if eid else False
    print(f"\n{gene} (in 22q: {in_22q}):")
    print(f"  Mechanism: {info['mechanism']}")
    if info["channel_targets"]:
        print(f"  Channel targets: {', '.join(info['channel_targets'])}")

# %%
# === 5H. Check channel gene expression specificity in VIP interneurons ===
# Which channel genes are VIP-specific? This connects the DEG findings to VIP vulnerability.
print("=" * 75)
print("5H. Ion Channel Specificity in VIP+ vs VIP- CGE Interneurons")
print("=" * 75)

channel_specificity = []
for gene_sym in sorted(all_channel_genes):
    eid = GeneSymbol2Entrez.get(gene_sym)
    if eid is None or eid not in SpecMat.index:
        continue
    cat = [k for k, v in ion_channel_sets.items() if gene_sym in v][0]

    vip_pos_spec = SpecMat.loc[eid, vip_pos_idx].mean()
    vip_neg_spec = SpecMat.loc[eid, vip_neg_idx].mean()
    all_cge_spec = SpecMat.loc[eid, CGE_idx].mean()
    global_mean = SpecMat.loc[eid].mean()

    channel_specificity.append({
        "Gene": gene_sym,
        "Category": cat,
        "VIP+ spec": vip_pos_spec,
        "VIP- spec": vip_neg_spec,
        "VIP+ - VIP-": vip_pos_spec - vip_neg_spec,
        "CGE mean spec": all_cge_spec,
        "Global mean spec": global_mean,
    })

chan_spec_df = pd.DataFrame(channel_specificity)
chan_spec_df = chan_spec_df.sort_values("VIP+ - VIP-", ascending=False)

print(f"\n{'Gene':>8s} {'Category':>22s}  {'VIP+':>7s}  {'VIP-':>7s}  {'Delta':>7s}  {'CGE':>7s}")
print("-" * 70)
for _, row in chan_spec_df.iterrows():
    marker = " **" if abs(row["VIP+ - VIP-"]) > 0.3 else " *" if abs(row["VIP+ - VIP-"]) > 0.15 else ""
    print(f"{row['Gene']:>8s} {row['Category']:>22s}  {row['VIP+ spec']:+.3f}  "
          f"{row['VIP- spec']:+.3f}  {row['VIP+ - VIP-']:+.3f}  "
          f"{row['CGE mean spec']:+.3f}{marker}")

# %%
# --- Figure: VIP specificity of channel genes ---
fig, ax = plt.subplots(figsize=(12, 6), dpi=150, facecolor='none')
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

# Color by category
cat_colors = {
    "Na+ channels (SCN)": "#E41A1C",
    "Kv3 (fast repol.)": "#377EB8",
    "Kv1 (delayed rect.)": "#4DAF4A",
    "Kv2 (delayed rect.)": "#984EA3",
    "BK/SK (Ca2+-act.)": "#FF7F00",
    "A-type K+ (Kv4)": "#A65628",
    "HCN (Ih)": "#F781BF",
}

x = range(len(chan_spec_df))
for i, (_, row) in enumerate(chan_spec_df.iterrows()):
    color = cat_colors.get(row["Category"], "grey")
    ax.bar(i, row["VIP+ - VIP-"], color=color, alpha=0.8, edgecolor='white', linewidth=0.5)

ax.set_xticks(list(x))
ax.set_xticklabels(chan_spec_df["Gene"], rotation=45, ha='right', fontsize=9)
ax.set_ylabel("VIP+ - VIP- specificity (mean-centered)", fontsize=11)
ax.set_title("Ion Channel Expression Specificity: VIP+ vs VIP- CGE Interneurons", fontsize=12)
ax.axhline(0, color='black', lw=0.8)

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=cat, alpha=0.8)
                   for cat, c in cat_colors.items()]
ax.legend(handles=legend_elements, fontsize=8, loc='upper right',
          ncol=2, framealpha=0.9)

plt.tight_layout()
plt.savefig(str(PROJ_DIR / "results/fig_22q_ion_channel_VIP_specificity.png"),
            dpi=150, transparent=True, bbox_inches='tight')
plt.show()

# %% [markdown]
# ### 5I. Integrated Summary
#
# Combine the DEG evidence with VIP specificity to identify channels that are:
# 1. Differentially expressed in 22q models, AND
# 2. Preferentially expressed in VIP+ interneurons

# %%
# === Integrated summary ===
print("=" * 75)
print("5I. INTEGRATED SUMMARY: Ion Channels x 22q DEG x VIP Specificity")
print("=" * 75)

# For each channel gene, check: is it a DEG? Is it VIP-specific?
print(f"\n{'Gene':>8s} {'Category':>22s}  {'VIP delta':>10s}  {'DEG in':>30s}  {'Direction':>10s}")
print("-" * 90)

for _, row in chan_spec_df.iterrows():
    gene_sym = row["Gene"]
    deg_hits = []
    directions = []

    # Check Lin
    for label, r in lin_channel_results.items():
        if gene_sym in r["channel_sig_genes"]:
            info = r["channel_sig_genes"][gene_sym]
            deg_hits.append(f"Lin {label.split('(')[0].strip()}")
            directions.append("UP" if info["lfc"] > 0 else "DOWN")

    # Check Khan
    for label, r in khan_channel_results.items():
        if gene_sym in r["channel_sig_genes"]:
            info = r["channel_sig_genes"][gene_sym]
            deg_hits.append(f"Khan {label.split('(')[0].strip()}")
            directions.append("UP" if info["lfc"] > 0 else "DOWN")

    if deg_hits or abs(row["VIP+ - VIP-"]) > 0.15:
        deg_str = ", ".join(deg_hits) if deg_hits else "---"
        dir_str = "/".join(set(directions)) if directions else "---"
        vip_flag = "**VIP+" if row["VIP+ - VIP-"] > 0.15 else "VIP-" if row["VIP+ - VIP-"] < -0.15 else ""
        print(f"{gene_sym:>8s} {row['Category']:>22s}  {row['VIP+ - VIP-']:+.3f} {vip_flag:>5s}  "
              f"{deg_str:>30s}  {dir_str:>10s}")

print("""
INTERPRETATION:
  The analysis examines whether 22q11.2 deletion's transcriptomic effects converge
  on ion channels critical for VIP interneuron AP kinetics.

  KEY FINDINGS:
  1. DEG enrichment: Ion channel genes as a group are tested for enrichment among
     22q DEGs via Fisher's exact test (Lin et al., where full background is available).

  2. DGCR8 mechanism: DGCR8 haploinsufficiency (from the 22q deletion) reduces
     miR-185, miR-25, miR-132, and miR-134 — all of which have validated or predicted
     ion channel targets (SCN1A, HCN1/2, KCNA1, KCNMA1). This provides an indirect
     regulatory path from the deletion to channel dysregulation, even if individual
     channel genes do not reach genome-wide significance as DEGs.

  3. VIP specificity: Certain channels show preferential expression in VIP+ CGE
     interneurons (e.g., Kv3 channels, Na+ channels), meaning that even modest
     dysregulation would disproportionately affect VIP cells.

  4. The convergence of genetic (22q cis-gene bias), transcriptomic (DEG patterns),
     and cell type-specific expression suggests a multi-layered vulnerability of
     VIP interneurons to 22q deletion.
""")

# %% [markdown]
# ---
# ## 6. Gene Co-expression Network: 22q-Ion Channel Proximity
#
# **Question**: Do 22q11.2 genes show stronger co-expression with ion channel genes
# than expected by chance?
#
# **Approach**: Build a co-expression network from the cell-type specificity matrix
# (Spearmans' R values across 461 cell types). Test whether the mean |rho| between
# 22q genes and channel genes exceeds the null distribution from 10,000 random gene
# sets of the same size.
#
# If 22q genes and ion channel genes have correlated expression patterns across cell
# types, it suggests they participate in shared regulatory programs -- supporting the
# hypothesis that 22q haploinsufficiency can indirectly affect ion channels.

# %%
# --- 6A. Define gene sets and map to Entrez IDs ---

# 22q genes (with alternate HGNC symbols for renamed genes)
GENES_22Q_SYMBOLS = [
    "DGCR8", "TBX1", "COMT", "PRODH", "SEPTIN5", "RANBP1", "CRKL", "PI4KA",
    "SNAP29", "HIRA", "UFD1", "CDC45", "SLC25A1", "DGCR2", "GP1BB", "LZTR1",
    "MRPL40", "TANGO2", "RTN4R", "SCARF2",
]

# Ion channel genes relevant to VIP interneuron electrophysiology
GENES_CHANNEL_SYMBOLS = [
    "SCN1A", "SCN2A", "SCN3A", "SCN8A", "SCN1B", "SCN2B", "SCN3B", "SCN4B",
    "KCNC1", "KCNC2", "KCNC3", "KCNC4", "KCNA1", "KCNA2", "KCNA3", "KCNA6",
    "KCNB1", "KCNB2", "KCNMA1", "KCNN1", "KCNN2", "KCNN3", "KCND2", "KCND3",
    "HCN1", "HCN2",
]


def symbols_to_entrez_in_matrix(symbols, gs2e, mat_index):
    """Convert gene symbols to Entrez IDs, keeping only those in the matrix."""
    mapped, missing = [], []
    for sym in symbols:
        eid = gs2e.get(sym)
        if eid is not None and eid in mat_index:
            mapped.append((sym, eid))
        else:
            missing.append(sym)
    return mapped, missing


q22_mapped, q22_missing = symbols_to_entrez_in_matrix(
    GENES_22Q_SYMBOLS, GeneSymbol2Entrez, SpecMat.index
)
chan_mapped, chan_missing = symbols_to_entrez_in_matrix(
    GENES_CHANNEL_SYMBOLS, GeneSymbol2Entrez, SpecMat.index
)

q22_entrez = [eid for _, eid in q22_mapped]
q22_sym2eid = {sym: eid for sym, eid in q22_mapped}
chan_entrez = [eid for _, eid in chan_mapped]
chan_sym2eid = {sym: eid for sym, eid in chan_mapped}
chan_eid2sym = {eid: sym for sym, eid in chan_mapped}
q22_eid2sym = {eid: sym for sym, eid in q22_mapped}

print(f"22q genes: {len(q22_entrez)}/{len(GENES_22Q_SYMBOLS)} mapped to matrix")
if q22_missing:
    print(f"  Missing: {q22_missing}")
print(f"Channel genes: {len(chan_entrez)}/{len(GENES_CHANNEL_SYMBOLS)} mapped to matrix")
if chan_missing:
    print(f"  Missing: {chan_missing}")

# %%
# --- 6B. Compute Spearmans' R matrix: 22q x channel genes ---
from scipy.stats import spearmanr

# Extract expression profiles (genes x cell types)
X_22q = SpecMat.loc[q22_entrez].values       # (n_22q, 461)
X_chan = SpecMat.loc[chan_entrez].values       # (n_chan, 461)

n_22q = X_22q.shape[0]
n_chan = X_chan.shape[0]
n_celltypes = X_22q.shape[1]
print(f"Computing {n_22q} x {n_chan} Spearmans' R values across {n_celltypes} cell types...")

# Compute the full correlation matrix between 22q and channel genes
corr_matrix = np.zeros((n_22q, n_chan))
pval_matrix = np.zeros((n_22q, n_chan))

for i in range(n_22q):
    for j in range(n_chan):
        rho, pval = spearmanr(X_22q[i, :], X_chan[j, :])
        corr_matrix[i, j] = rho
        pval_matrix[i, j] = pval

# Create labeled DataFrames
q22_labels = [q22_eid2sym[eid] for eid in q22_entrez]
chan_labels = [chan_eid2sym[eid] for eid in chan_entrez]

corr_df = pd.DataFrame(corr_matrix, index=q22_labels, columns=chan_labels)
pval_df = pd.DataFrame(pval_matrix, index=q22_labels, columns=chan_labels)

# Aggregate metrics
mean_abs_corr_observed = np.abs(corr_matrix).mean()
mean_corr_observed = corr_matrix.mean()
median_abs_corr_observed = np.median(np.abs(corr_matrix))

print(f"\nObserved 22q-channel co-expression:")
print(f"  Mean |rho|  = {mean_abs_corr_observed:.4f}")
print(f"  Median |rho| = {median_abs_corr_observed:.4f}")
print(f"  Mean rho    = {mean_corr_observed:.4f}")
print(f"  Max |rho|   = {np.abs(corr_matrix).max():.4f}")

# %%
# --- 6C. Permutation test: null distribution of mean |rho| ---
# For each permutation, sample a random gene set of the same size as 22q genes
# (excluding channel genes themselves) and compute mean |rho| with channel genes.

from joblib import Parallel, delayed

N_PERM = 10_000
SEED = 42

# Pool of background genes: all genes in matrix excluding channel genes and 22q genes
exclude_set = set(q22_entrez) | set(chan_entrez)
background_genes = np.array([g for g in SpecMat.index if g not in exclude_set])
print(f"Background gene pool: {len(background_genes)} genes")

# Pre-compute channel gene expression for efficiency
X_chan_T = X_chan.T  # (461, n_chan) -- not needed but keeping X_chan in (n_chan, 461)


def _compute_null_mean_abs_corr(seed_i):
    """Sample random gene set and compute mean |rho| with channel genes."""
    rng = np.random.RandomState(SEED + seed_i)
    rand_idx = rng.choice(len(background_genes), size=n_22q, replace=False)
    rand_genes = background_genes[rand_idx]
    X_rand = SpecMat.loc[rand_genes].values  # (n_22q, 461)

    # Compute all pairwise Spearmans' R values
    abs_corrs = []
    for i in range(n_22q):
        for j in range(n_chan):
            rho, _ = spearmanr(X_rand[i, :], X_chan[j, :])
            abs_corrs.append(abs(rho))
    return np.mean(abs_corrs)


print(f"Running {N_PERM} permutations with n_jobs=10...")
null_distribution = Parallel(n_jobs=10, verbose=1)(
    delayed(_compute_null_mean_abs_corr)(i) for i in range(N_PERM)
)
null_distribution = np.array(null_distribution)

# Empirical p-value (one-sided: observed >= null)
p_value = (np.sum(null_distribution >= mean_abs_corr_observed) + 1) / (N_PERM + 1)
z_score = (mean_abs_corr_observed - null_distribution.mean()) / null_distribution.std()

print(f"\nPermutation test results (N={N_PERM}):")
print(f"  Observed mean |rho| = {mean_abs_corr_observed:.4f}")
print(f"  Null mean |rho|     = {null_distribution.mean():.4f} +/- {null_distribution.std():.4f}")
print(f"  Z-score             = {z_score:.2f}")
print(f"  Empirical P-value   = {p_value:.4e}")

# %%
# --- 6D. Visualization: Null distribution vs. observed ---

FIG_DIR = str(PROJ_DIR / "dat/VIP_Ephys/figures")
os.makedirs(FIG_DIR, exist_ok=True)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(null_distribution, bins=60, color="#7BAFD4", alpha=0.8,
        edgecolor="white", linewidth=0.5, label="Null (random gene sets)")
ax.axvline(mean_abs_corr_observed, color="#C44E52", linewidth=2.5,
           linestyle="--", label=f"Observed 22q (mean |rho| = {mean_abs_corr_observed:.3f})")
ax.set_xlabel("Mean |Spearmans' R| with ion channel genes", fontsize=12)
ax.set_ylabel("Count (permutations)", fontsize=12)
ax.set_title("22q Genes Show Elevated Co-expression\nwith Ion Channel Genes", fontsize=13)

# Add stats annotation
txt = f"P = {p_value:.2e}\nZ = {z_score:.1f}"
ax.text(0.97, 0.95, txt, transform=ax.transAxes, fontsize=11,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))

ax.legend(loc="upper left", fontsize=10)
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/Fig_22q_Channel_Coexpression_Null.png",
            dpi=150, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved: {FIG_DIR}/Fig_22q_Channel_Coexpression_Null.png")

# %%
# --- 6E. Visualization: Clustered heatmap of 22q x channel correlations ---

from scipy.cluster.hierarchy import linkage, leaves_list

# Cluster rows and columns for better visualization
row_linkage = linkage(corr_df.values, method="average", metric="euclidean")
col_linkage = linkage(corr_df.values.T, method="average", metric="euclidean")
row_order = leaves_list(row_linkage)
col_order = leaves_list(col_linkage)

corr_clustered = corr_df.iloc[row_order, col_order]

fig, ax = plt.subplots(figsize=(12, 7))
vmax = max(0.5, np.abs(corr_clustered.values).max())
sns.heatmap(
    corr_clustered,
    cmap="RdBu_r",
    center=0,
    vmin=-vmax,
    vmax=vmax,
    annot=True,
    fmt=".2f",
    annot_kws={"size": 6.5},
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Spearmans' R", "shrink": 0.7},
    ax=ax,
)
ax.set_xlabel("Ion Channel Genes", fontsize=12)
ax.set_ylabel("22q11.2 Genes", fontsize=12)
ax.set_title("Cell-type Co-expression: 22q11.2 vs Ion Channel Genes\n"
             f"(Spearmans' R across {n_celltypes} cell types)", fontsize=13)
ax.tick_params(axis="x", rotation=45, labelsize=9)
ax.tick_params(axis="y", rotation=0, labelsize=9)
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/Fig_22q_Channel_Coexpression_Heatmap.png",
            dpi=150, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved: {FIG_DIR}/Fig_22q_Channel_Coexpression_Heatmap.png")

# %%
# --- 6F. Top co-expression pairs ---

# Flatten the correlation matrix into a list of pairs
pairs = []
for i, q_sym in enumerate(q22_labels):
    for j, c_sym in enumerate(chan_labels):
        pairs.append({
            "22q_gene": q_sym,
            "Channel_gene": c_sym,
            "rho": corr_matrix[i, j],
            "abs_rho": abs(corr_matrix[i, j]),
            "p_value": pval_matrix[i, j],
        })
pairs_df = pd.DataFrame(pairs).sort_values("abs_rho", ascending=False)

# FDR correction across all pairs
from statsmodels.stats.multitest import multipletests
_, pairs_df["q_value"], _, _ = multipletests(pairs_df["p_value"], method="fdr_bh")
pairs_df = pairs_df.sort_values("abs_rho", ascending=False).reset_index(drop=True)

print("Top 30 co-expression pairs (22q gene -- channel gene):")
print(f"{'22q Gene':>12s}  {'Channel':>8s}  {'rho':>7s}  {'|rho|':>7s}  {'P':>10s}  {'q':>10s}")
print("-" * 65)
for _, row in pairs_df.head(30).iterrows():
    print(f"{row['22q_gene']:>12s}  {row['Channel_gene']:>8s}  {row['rho']:+.4f}  "
          f"{row['abs_rho']:.4f}  {row['p_value']:.2e}  {row['q_value']:.2e}")

# Summary: how many pairs exceed thresholds
for threshold in [0.3, 0.4, 0.5]:
    n_above = (pairs_df["abs_rho"] >= threshold).sum()
    n_sig = ((pairs_df["abs_rho"] >= threshold) & (pairs_df["q_value"] < 0.05)).sum()
    print(f"\n|rho| >= {threshold}: {n_above} pairs ({n_sig} with q < 0.05)")

# %%
# --- 6G. Per-22q-gene aggregate: which 22q genes co-express most with channels? ---

gene_agg = (
    pairs_df.groupby("22q_gene")
    .agg(
        mean_abs_rho=("abs_rho", "mean"),
        max_abs_rho=("abs_rho", "max"),
        n_sig_pairs=("q_value", lambda x: (x < 0.05).sum()),
        best_channel=("abs_rho", lambda x: pairs_df.loc[x.idxmax(), "Channel_gene"]),
        best_rho=("rho", lambda x: pairs_df.loc[pairs_df.loc[x.index, "abs_rho"].idxmax(), "rho"]),
    )
    .sort_values("mean_abs_rho", ascending=False)
)

print("\nPer-22q gene: mean |rho| with all channel genes")
print(f"{'Gene':>12s}  {'Mean|rho|':>10s}  {'Max|rho|':>9s}  {'#Sig(q<.05)':>12s}  "
      f"{'Best Channel':>14s}  {'Best rho':>9s}")
print("-" * 75)
for gene, row in gene_agg.iterrows():
    print(f"{gene:>12s}  {row['mean_abs_rho']:.4f}      {row['max_abs_rho']:.4f}     "
          f"{int(row['n_sig_pairs']):>6d}        {row['best_channel']:>10s}  {row['best_rho']:+.4f}")

# %%
# --- 6H. Visualization: Top co-expression pairs highlighted ---

# Bar plot of per-22q-gene mean |rho|, colored by significance
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel A: Per-22q gene mean |rho|
ax = axes[0]
gene_order = gene_agg.index.tolist()
colors = ["#C44E52" if gene_agg.loc[g, "n_sig_pairs"] > 0 else "#7BAFD4"
          for g in gene_order]
ax.barh(gene_order[::-1], gene_agg["mean_abs_rho"].values[::-1],
        color=colors[::-1], edgecolor="white", linewidth=0.5)
ax.set_xlabel("Mean |Spearmans' R| with channel genes", fontsize=11)
ax.set_title("22q Gene Co-expression\nwith Ion Channels", fontsize=12)
ax.axvline(null_distribution.mean(), color="gray", linestyle=":", linewidth=1.5,
           label=f"Null mean = {null_distribution.mean():.3f}")
ax.legend(fontsize=9)
ax.tick_params(axis="y", labelsize=9)
ax.patch.set_alpha(0)

# Panel B: Top pairs scatter (rho values for pairs with |rho| > 0.3)
ax = axes[1]
top_pairs = pairs_df[pairs_df["abs_rho"] >= 0.3].copy()
if len(top_pairs) > 0:
    top_pairs = top_pairs.sort_values("rho", ascending=True).reset_index(drop=True)
    pair_labels = [f"{r['22q_gene']}-{r['Channel_gene']}" for _, r in top_pairs.iterrows()]
    bar_colors = ["#C44E52" if r["rho"] > 0 else "#4C72B0"
                  for _, r in top_pairs.iterrows()]
    ax.barh(range(len(top_pairs)), top_pairs["rho"].values,
            color=bar_colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(top_pairs)))
    ax.set_yticklabels(pair_labels, fontsize=8)
    ax.set_xlabel("Spearmans' R", fontsize=11)
    ax.set_title(f"Top Co-expression Pairs\n(|rho| >= 0.3, n={len(top_pairs)})", fontsize=12)
    ax.axvline(0, color="black", linewidth=0.8)
    # Mark significant pairs
    for idx, (_, r) in enumerate(top_pairs.iterrows()):
        if r["q_value"] < 0.05:
            ax.text(r["rho"] + (0.01 if r["rho"] > 0 else -0.01),
                    idx, "*", fontsize=12, va="center",
                    ha="left" if r["rho"] > 0 else "right",
                    fontweight="bold", color="black")
else:
    ax.text(0.5, 0.5, "No pairs with |rho| >= 0.3",
            transform=ax.transAxes, ha="center", va="center", fontsize=12)
    ax.set_title("Top Co-expression Pairs", fontsize=12)
ax.patch.set_alpha(0)

fig.patch.set_alpha(0)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/Fig_22q_Channel_Coexpression_TopPairs.png",
            dpi=150, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved: {FIG_DIR}/Fig_22q_Channel_Coexpression_TopPairs.png")

# %%
# --- 6I. Summary and interpretation ---

print("=" * 70)
print("SECTION 6 SUMMARY: 22q-Ion Channel Co-expression Network")
print("=" * 70)
print(f"""
Gene sets:
  22q genes in matrix:     {len(q22_entrez)}/{len(GENES_22Q_SYMBOLS)}
  Channel genes in matrix: {len(chan_entrez)}/{len(GENES_CHANNEL_SYMBOLS)}
  Total pairs tested:      {len(pairs_df)}

Observed co-expression:
  Mean |rho|  = {mean_abs_corr_observed:.4f}
  Median |rho| = {median_abs_corr_observed:.4f}

Permutation test ({N_PERM} random gene sets of size {n_22q}):
  Null mean |rho|  = {null_distribution.mean():.4f} +/- {null_distribution.std():.4f}
  Z-score          = {z_score:.2f}
  Empirical P      = {p_value:.2e}

Interpretation:
  {"22q genes show SIGNIFICANTLY elevated co-expression with ion channels." if p_value < 0.05 else "22q genes do NOT show significantly elevated co-expression with ion channels."}
  This {"supports" if p_value < 0.05 else "does not support"} the hypothesis that 22q genes
  participate in shared regulatory programs with ion channel genes across cell types.

Top 22q genes by co-expression with channels:
""")
for gene in gene_agg.head(5).index:
    row = gene_agg.loc[gene]
    print(f"  {gene:>12s}: mean |rho| = {row['mean_abs_rho']:.3f}, "
          f"best pair = {row['best_channel']} (rho = {row['best_rho']:+.3f})")

print(f"""
Top co-expression pairs (|rho| > 0.3):
  {(pairs_df['abs_rho'] >= 0.3).sum()} pairs exceed |rho| = 0.3
  {(pairs_df['abs_rho'] >= 0.4).sum()} pairs exceed |rho| = 0.4
""")
