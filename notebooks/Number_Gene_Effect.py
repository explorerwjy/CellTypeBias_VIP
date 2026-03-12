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
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
sys.path.insert(1, '/home/jw3514/Work/UNIMED/src')
from CellType_PSY import *
#from UNIMED import *

# Load config.yaml  
import yaml
with open(ProjDIR + '/config/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

#import scanpy as sc
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

try:
    os.chdir(f"{ProjDIR}/notebooks/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:
    print(f"Unexpected error: {e}")

# %%
# Load expression matrix
expression_matrix = config['analysis_types']['Centering']
print(expression_matrix)

# %%
HumanCT_Z2_HCT = pd.read_csv(ProjDIR + expression_matrix, index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
HumanCT_Z2_HCT.shape

# %%
GeneWeightDIR = "../dat/GeneWeights/"
Bias_Save_Dir = ProjDIR + "results/main_results/random/Centering/"
if not os.path.exists(Bias_Save_Dir): # make dir if not exists
    os.makedirs(Bias_Save_Dir)

# Load BGMR for consistent gene weight correction across all disorders
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv", low_memory=False)
BGMR["entrez_id"] = BGMR["entrez_id"].astype(int)
BGMR = BGMR.set_index("entrez_id")

# %%
# Compute SCZ reference inline from top-61 of the sweep's gene weight file
# (ensures R=1 at N=61; pipeline file has only 53 genes after nopLI/exclude_Mis2 filtering)

# %%
# The file has no header row, so treat all rows as data and provide a custom header
custom_header = ["Weight"]
GeneWeights = pd.read_csv(
    f"{GeneWeightDIR}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    index_col=0,
    header=None,
    names=custom_header,
)

GeneWeights.head()

# %%
GeneDF = pd.read_csv("../dat/SCZ.ALLGENE.MutCountModified.csv", index_col=0)

# %%
for i, row in GeneWeights.iterrows():
    if i in GeneDF.index.values:
        GeneWeights.loc[i, "Pval"] = GeneDF.loc[i, "P meta"]
        GeneWeights.loc[i, "FDR"] = GeneDF.loc[i, "Q meta"]
        GeneWeights.loc[i, "Gene Symbol"] = GeneDF.loc[i, "Gene Symbol"]

# %%
# Reference: top-61 from the same gene weight file used in the sweep
N_REF = 61
ref_SCZ_gw = dict(zip(GeneWeights.index[:N_REF], GeneWeights["Weight"].values[:N_REF]))
SCZ_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_SCZ_gw)
SCZ_Bias = AnnotateCTDat(SCZ_Bias, Anno)

# %%
# Calculate correlation as genes are removed, sorted by LofZ (descending)
GeneIdx, Corr_SCZ, Pval_SCZ = [], [], []
Corr_SCZ_unweighted = []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp_SCZ_GW = dict(zip(GeneWeights.index.values[:i], GeneWeights["Weight"].values[:i]))
    tmp_SCZ_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_SCZ_GW)
    tmp_SCZ_Bias = AnnotateCTDat(tmp_SCZ_Bias, Anno)
    r, p = GetSingeCellBiasCorr(SCZ_Bias, tmp_SCZ_Bias, efflabel="EFFECT")
        
    Corr_SCZ.append(r)
    Pval_SCZ.append(p)
    
    # Unweighted version (each gene has weight 1)
    tmp_SCZ_GW_unweighted = dict(zip(GeneWeights.index.values[:i], np.ones(i)))
    tmp_SCZ_Bias_unweighted = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_SCZ_GW_unweighted)
    tmp_SCZ_Bias_unweighted = AnnotateCTDat(tmp_SCZ_Bias_unweighted, Anno)
    r_unweighted, _ = GetSingeCellBiasCorr(SCZ_Bias, tmp_SCZ_Bias_unweighted, efflabel="EFFECT")
    Corr_SCZ_unweighted.append(r_unweighted)

# %%
import matplotlib.ticker as mticker

def plot_gene_set_correlation(
    GeneIdx_top,
    Corr_top,
    GeneSig_top,
    corr_label="Bias Correlation",
    pval_label=r"$-\log_{10}$(min p-value of Included Genes)",
    xlabel="Number of Genes",
    ylabel_corr="Bias Correlation with Main ASD with ID Set",
    ylabel_pval=None,
    figsize=(6, 4),
    dpi=120,
    color_corr='#1f77b4',
    color_pval='#d62728',
    legend_loc='lower left',
    legend_fontsize=10,
    ylim_corr=(0, 1.01),
    nbins_x=8,
    title=None,
    ax1=None,
    ax2=None,
    show=True,
    Corr_unweighted=None,
    corr_unweighted_label="Bias Correlation (Unweighted)",
    color_unweighted='#2ca02c',
):
    """
    Plot gene set correlation and -log10(p-value) as a function of number of genes.

    Parameters
    ----------
    GeneIdx_top : array-like
        X-axis values (number of genes included).
    Corr_top : array-like
        Correlation values to plot on left y-axis.
    GeneSig_top : array-like
        -log10(p-value) values to plot on right y-axis.
    corr_label : str
        Label for the correlation line.
    pval_label : str
        Label for the p-value line.
    xlabel : str
        X-axis label.
    ylabel_corr : str
        Y-axis label for correlation.
    ylabel_pval : str or None
        Y-axis label for p-value. If None, uses pval_label.
    figsize : tuple
        Figure size.
    dpi : int
        Figure DPI.
    color_corr : str
        Color for correlation line.
    color_pval : str
        Color for p-value line.
    legend_loc : str
        Location for legend.
    legend_fontsize : int
        Font size for legend.
    ylim_corr : tuple
        Y-axis limits for correlation.
    nbins_x : int
        Number of bins for x-axis major ticks.
    title : str or None
        Optional plot title.
    ax1, ax2 : matplotlib axes or None
        Optionally provide axes to plot on.
    show : bool
        Whether to call plt.show().
    """
    import matplotlib.pyplot as plt

    if ax1 is None or ax2 is None:
        fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax1.figure

    # Plot correlation
    ax1.plot(GeneIdx_top, Corr_top, color=color_corr, linewidth=2, marker='o', markersize=4, label=corr_label)
    
    # Plot unweighted correlation if provided
    if Corr_unweighted is not None:
        ax1.plot(GeneIdx_top, Corr_unweighted, color=color_unweighted, linewidth=2, marker='^', markersize=4, label=corr_unweighted_label, linestyle='--')
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel_corr, color=color_corr, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_corr)
    ax1.set_ylim(*ylim_corr)
    ax1.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=nbins_x))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax1.xaxis.set_minor_formatter(mticker.NullFormatter())

    # Only draw horizontal grid lines, and only behind the lines, not over text
    ax1.grid(True, which='major', axis='y', linestyle='--', alpha=0.5, zorder=1)
    # Remove x grid to avoid covering x-labels/ticks
    ax1.grid(False, which='major', axis='x', zorder=1, alpha=0.5)

    # Plot -log10(P-value) on secondary axis
    if ax2 is None:
        ax2 = ax1.twinx()
    ax2.plot(GeneIdx_top, GeneSig_top, color=color_pval, linewidth=2, marker='s', markersize=4, label=pval_label, zorder=3)
    ax2.set_ylabel(ylabel_pval if ylabel_pval is not None else pval_label, color=color_pval, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_pval)
    ax2.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    # Add legend in the bottom left
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc=legend_loc, fontsize=legend_fontsize, frameon=False)

    # Publication style tweaks
    fig.tight_layout(pad=2)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.2)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.2)
    plt.subplots_adjust(top=0.95, right=0.88)
    if title is not None:
        ax1.set_title(title, fontsize=14, pad=10)
    if show:
        plt.show()
    return fig, ax1, ax2


# %%
plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_SCZ,
    GeneSig_top=-np.log10(np.array(GeneWeights["Pval"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main SCZ Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_SCZ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)"
)

# %% [markdown]
# # ASD

# %%
# ASD references will be computed inline after weights are built (below)

# %%
Spark_Denovo = pd.read_excel("../dat/41588_2022_1148_MOESM4_ESM.xlsx",
                           skiprows=2, sheet_name="Table S7")
Spark_Denovo = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]!="."]
Spark_Denovo.shape

# %%
Mut_n_IQ = pd.read_csv("../dat/ASD_IQ_Mut.csv")
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"]>70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"]<=70]

# %%
Spark_Denovo_sub = Spark_Denovo[["EntrezID", "HGNC", "pDenovoWEST_Meta", "AutismMerged_LoF", "AutismMerged_Dmis_REVEL0.5"]]
Spark_Denovo_sub = Spark_Denovo_sub.set_index("EntrezID")
for i, row in Spark_Denovo_sub.iterrows():
    LIQ_counts = LowIQMuts_ALL[LowIQMuts_ALL["Entrez"]==i].shape[0]
    HIQ_counts = HighIQMuts_ALL[HighIQMuts_ALL["Entrez"]==i].shape[0]
    Spark_Denovo_sub.loc[i, "LIQ_counts"] = LIQ_counts
    Spark_Denovo_sub.loc[i, "HIQ_counts"] = HIQ_counts
    if row["pDenovoWEST_Meta"] == 0:
        Spark_Denovo_sub.loc[i, "pDenovoWEST_Meta"] = 1e-6


# %%
# BGMR-correct ASD IQ mutation counts to match pipeline gene weights
# Pipeline uses wLGD=1, wMis=1, REVEL>0.5 filtering, BGMR subtraction
# Here HIQ_counts/LIQ_counts are raw mutation counts (all variant types combined).
# To match the pipeline, we need BGMR-corrected weights per gene.
# Load the BGMR-corrected gene weight files directly (same files the pipeline uses).
HIQ_GW_bgmr = pd.read_csv(GeneWeightDIR + "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
                           index_col=0, header=None, names=["Weight"])
LIQ_GW_bgmr = pd.read_csv(GeneWeightDIR + "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
                           index_col=0, header=None, names=["Weight"])

# For the sweep we need weights for genes beyond 61 too.
# Use raw counts for genes outside the top-61 (they weren't BGMR-corrected),
# but ensure top-61 matches pipeline exactly.
# Build extended weight tables: BGMR-corrected for genes in the bgmr file,
# raw counts for the rest.
for label_col, bgmr_gw in [("HIQ_weight_bgmr", HIQ_GW_bgmr), ("LIQ_weight_bgmr", LIQ_GW_bgmr)]:
    Spark_Denovo_sub[label_col] = np.nan
    for g in Spark_Denovo_sub.index:
        if g in bgmr_gw.index:
            Spark_Denovo_sub.loc[g, label_col] = bgmr_gw.loc[g, "Weight"]

# For genes not in the bgmr file, fall back to raw counts
Spark_Denovo_sub["HIQ_weight_bgmr"] = Spark_Denovo_sub["HIQ_weight_bgmr"].fillna(
    Spark_Denovo_sub["HIQ_counts"])
Spark_Denovo_sub["LIQ_weight_bgmr"] = Spark_Denovo_sub["LIQ_weight_bgmr"].fillna(
    Spark_Denovo_sub["LIQ_counts"])

# %%
# Compute ASD references inline from top-61 of the same gene weights used in the sweep
ref_HIQ_gw = dict(zip(Spark_Denovo_sub.index[:N_REF], Spark_Denovo_sub["HIQ_weight_bgmr"].values[:N_REF]))
HIQ_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_HIQ_gw)
HIQ_Bias = AnnotateCTDat(HIQ_Bias, Anno)

ref_LIQ_gw = dict(zip(Spark_Denovo_sub.index[:N_REF], Spark_Denovo_sub["LIQ_weight_bgmr"].values[:N_REF]))
LIQ_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_LIQ_gw)
LIQ_Bias = AnnotateCTDat(LIQ_Bias, Anno)

# %%
# Calculate correlation as genes are added — HIQ (BGMR-corrected)
GeneIdx_HIQ, Corr_HIQ, Pval_HIQ = [], [], []
Corr_HIQ_unweighted = []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp_HIQ_ASD_GW = dict(zip(Spark_Denovo_sub.index.values[:i], Spark_Denovo_sub["HIQ_weight_bgmr"].values[:i]))
    tmp_HIQ_ASD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_HIQ_ASD_GW)
    tmp_HIQ_ASD_Bias = AnnotateCTDat(tmp_HIQ_ASD_Bias, Anno)
    r, p = GetSingeCellBiasCorr(HIQ_Bias, tmp_HIQ_ASD_Bias, efflabel="EFFECT")

    Corr_HIQ.append(r)
    Pval_HIQ.append(p)

    # Unweighted version (each gene has weight 1)
    tmp_HIQ_ASD_GW_unweighted = dict(zip(Spark_Denovo_sub.index.values[:i], np.ones(i)))
    tmp_HIQ_ASD_Bias_unweighted = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_HIQ_ASD_GW_unweighted)
    tmp_HIQ_ASD_Bias_unweighted = AnnotateCTDat(tmp_HIQ_ASD_Bias_unweighted, Anno)
    r_unweighted, _ = GetSingeCellBiasCorr(HIQ_Bias, tmp_HIQ_ASD_Bias_unweighted, efflabel="EFFECT")
    Corr_HIQ_unweighted.append(r_unweighted)

# %%
# Calculate correlation as genes are added — LIQ (BGMR-corrected)
GeneIdx_LIQ, Corr_LIQ, Pval_LIQ = [], [], []
Corr_LIQ_unweighted = []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp_LIQ_ASD_GW = dict(zip(Spark_Denovo_sub.index.values[:i], Spark_Denovo_sub["LIQ_weight_bgmr"].values[:i]))
    tmp_LIQ_ASD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_LIQ_ASD_GW)
    tmp_LIQ_ASD_Bias = AnnotateCTDat(tmp_LIQ_ASD_Bias, Anno)
    r, p = GetSingeCellBiasCorr(LIQ_Bias, tmp_LIQ_ASD_Bias, efflabel="EFFECT")

    Corr_LIQ.append(r)
    Pval_LIQ.append(p)

    # Unweighted version (each gene has weight 1)
    tmp_LIQ_ASD_GW_unweighted = dict(zip(Spark_Denovo_sub.index.values[:i], np.ones(i)))
    tmp_LIQ_ASD_Bias_unweighted = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_LIQ_ASD_GW_unweighted)
    tmp_LIQ_ASD_Bias_unweighted = AnnotateCTDat(tmp_LIQ_ASD_Bias_unweighted, Anno)
    r_unweighted, _ = GetSingeCellBiasCorr(LIQ_Bias, tmp_LIQ_ASD_Bias_unweighted, efflabel="EFFECT")
    Corr_LIQ_unweighted.append(r_unweighted)

# %%
plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_HIQ,
    GeneSig_top=-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main ASD Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_HIQ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)"
)

plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_LIQ,
    GeneSig_top=-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main ASD ID Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_LIQ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)"
)

# %% [markdown]
# # DDD

# %%
# DDD reference will be computed inline after DDD_Genes_sub is built (below)

# %%
DDD_Genes = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
DDD_Genes = DDD_Genes.sort_values("denovoWEST_p_full")
entrez_ids = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_Genes["symbol"].values]
DDD_Genes["EntrezID"] = entrez_ids
DDD_Genes.shape

# %%
DDD_Genes_sub = DDD_Genes[["symbol", "EntrezID", "denovoWEST_p_full", "missense_variant", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost"]].reset_index(drop=True)
# Keep all genes — don't filter by BGMR membership. Genes without valid EntrezID
# or BGMR data get weight=0 (padded), preserving the true significance rank count.

# %%
# Helper: compute BGMR-corrected weights, padding invalid genes with weight=0
def _ddd_gene_weights(df, Nproband, BGMR):
    """Compute DDD gene weights; genes not in BGMR get weight=0 (preserves rank count)."""
    valid = df[df["EntrezID"].isin(BGMR.index)]
    if len(valid) > 0:
        return Aggregate_Gene_Weights_NDD(valid, Nproband=Nproband, BGMR=BGMR, wLGD=1, wMis=1)
    return {}

# %%
# Compute DDD reference inline from top-61 (preserving original rank, padding filtered genes)
N_DDD = 31058
ref_DDD_gw = _ddd_gene_weights(DDD_Genes_sub.head(N_REF), Nproband=N_DDD, BGMR=BGMR)
DDD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_DDD_gw)
DDD_Bias = AnnotateCTDat(DDD_Bias, Anno)

# %%
# Calculate correlation as genes are added, sorted by significance (descending)
Corr_DDD, Pval_DDD = [], []
Corr_DDD_unweighted = []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp = DDD_Genes_sub.head(i)
    tmp_DDD_GW = _ddd_gene_weights(tmp, Nproband=N_DDD, BGMR=BGMR)
    tmp_DDD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_DDD_GW)
    tmp_DDD_Bias = AnnotateCTDat(tmp_DDD_Bias, Anno)
    r, p = GetSingeCellBiasCorr(DDD_Bias, tmp_DDD_Bias, efflabel="EFFECT")

    Corr_DDD.append(r)
    Pval_DDD.append(p)

    # Unweighted version: valid genes get weight=1, invalid get weight=0
    valid_ids = tmp[tmp["EntrezID"].isin(BGMR.index)]["EntrezID"].values
    tmp_DDD_GW_unweighted = dict(zip(valid_ids, np.ones(len(valid_ids))))
    tmp_DDD_Bias_unweighted = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_DDD_GW_unweighted)
    tmp_DDD_Bias_unweighted = AnnotateCTDat(tmp_DDD_Bias_unweighted, Anno)
    r_unweighted, _ = GetSingeCellBiasCorr(DDD_Bias, tmp_DDD_Bias_unweighted, efflabel="EFFECT")
    Corr_DDD_unweighted.append(r_unweighted)

# %%
plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_DDD,
    GeneSig_top=-np.log10(np.array(DDD_Genes["denovoWEST_p_full"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main DDD Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_DDD_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)"
)

# %% [markdown]
# # Reviewer Figure R3.2a: SCZ Split-Half Bias Comparison
#
# Do added SCZ genes (62-200) carry concordant cell-type signal,
# or are they noise? Compare bias profiles of top-61 vs added vs random.
#
# Supports reply to Reviewer #3, Point 2, Paragraph 1 (concordant signal).

# %%
from scipy.stats import spearmanr

# Compute bias for three gene subsets
top61_gw = dict(zip(GeneWeights.index[:61], GeneWeights["Weight"].values[:61]))
added_gw = dict(zip(GeneWeights.index[61:200], GeneWeights["Weight"].values[61:200]))

top61_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, top61_gw)
top61_bias = AnnotateCTDat(top61_bias, Anno)

added_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, added_gw)
added_bias = AnnotateCTDat(added_bias, Anno)

# Random control: 139 genes not in SCZ top-200, with average weight of added genes
np.random.seed(42)
exclude_set = set(GeneWeights.index[:200].astype(int))
gene_pool = [g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set]
random_genes = np.random.choice(gene_pool, size=139, replace=False)
avg_added_weight = GeneWeights["Weight"].values[61:200].mean()
random_gw = dict(zip(random_genes, np.ones(139) * avg_added_weight))
random_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, random_gw)
random_bias = AnnotateCTDat(random_bias, Anno)

# Correlations
common_idx = top61_bias.index
neur_mask = top61_bias.index.isin(Neur_idx)

r_added_all, p_added_all = spearmanr(top61_bias["EFFECT"], added_bias.loc[common_idx, "EFFECT"])
r_added_neur, p_added_neur = spearmanr(top61_bias.loc[neur_mask, "EFFECT"], added_bias.loc[top61_bias.loc[neur_mask].index, "EFFECT"])
r_rand_all, p_rand_all = spearmanr(top61_bias["EFFECT"], random_bias.loc[common_idx, "EFFECT"])

print(f"Top-61 vs Added (62-200):  all CTs r={r_added_all:.3f} (p={p_added_all:.2e}), neurons r={r_added_neur:.3f} (p={p_added_neur:.2e})")
print(f"Top-61 vs Random (n=139):  all CTs r={r_rand_all:.3f} (p={p_rand_all:.2e})")

print("\nTop 5 superclusters — Top-61 genes:")
for sc, val in top61_bias.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False).head(5).items():
    print(f"  {sc}: {val:.4f}")

print("\nTop 5 superclusters — Added genes (62-200):")
for sc, val in added_bias.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False).head(5).items():
    print(f"  {sc}: {val:.4f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=120)

for ax_idx, (comp_bias, comp_label, r_val, p_val) in enumerate([
    (added_bias, "Added genes (62-200)", r_added_all, p_added_all),
    (random_bias, "Random genes (n=139)", r_rand_all, p_rand_all),
]):
    ax = axes[ax_idx]
    x = top61_bias["EFFECT"].values
    y = comp_bias.loc[common_idx, "EFFECT"].values

    ax.scatter(x[neur_mask], y[neur_mask], color="red", alpha=0.4, s=20,
               edgecolors="white", lw=0.3, label="Neuronal", zorder=3)
    ax.scatter(x[~neur_mask], y[~neur_mask], color="blue", alpha=0.5, s=20,
               edgecolors="white", lw=0.3, label="Non-neuronal", zorder=4)

    # Reference line
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, "k--", lw=0.8, alpha=0.4)
    ax.axhline(0, color="gray", lw=0.5, alpha=0.3)
    ax.axvline(0, color="gray", lw=0.5, alpha=0.3)

    ax.set_xlabel("Top-61 SCZ genes — EFFECT", fontsize=11)
    ax.set_ylabel(f"{comp_label} — EFFECT", fontsize=11)
    ax.set_title(chr(65 + ax_idx), fontweight="bold", loc="left", fontsize=14)
    ax.legend(fontsize=9, framealpha=0.8)
    ax.text(0.97, 0.03, f"ρ = {r_val:.3f}\np = {p_val:.1e}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

fig.suptitle("R3.2a — SCZ: Added genes carry concordant cell-type signal", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()

# %% [markdown]
# # Reviewer Figure R3.2b: Sliding Window Correlation Decay (SCZ)
#
# Use a sliding window of 61 genes to show how far down the ranked list
# the cell-type bias signal persists.
#
# Supports reply to Reviewer #3, Point 2, Paragraph 1 (signal persistence).

# %%
GeneDF = pd.read_csv("../dat/SCZ.ALLGENE.MutCountModified.csv", index_col=0)

# Add P-values to full gene weight list (if not already present)
for i in GeneWeights.index:
    if i in GeneDF.index.values:
        GeneWeights.loc[i, "Pval"] = GeneDF.loc[i, "P meta"]

window_size = 61
start_ranks = np.arange(0, 200 - window_size + 1, 10)
window_corrs = []
window_mean_pvals = []

for start in start_ranks:
    end = start + window_size
    w_gw = dict(zip(GeneWeights.index[start:end], GeneWeights["Weight"].values[start:end]))
    w_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, w_gw)
    w_bias = AnnotateCTDat(w_bias, Anno)
    r, _ = spearmanr(top61_bias["EFFECT"], w_bias.loc[common_idx, "EFFECT"])
    window_corrs.append(r)

    pvals_in_window = GeneWeights.iloc[start:end]["Pval"].dropna().astype(float)
    window_mean_pvals.append(pvals_in_window.mean() if len(pvals_in_window) > 0 else np.nan)

# %%
fig, ax1 = plt.subplots(figsize=(7, 4.5), dpi=120)

color_corr = "#1f77b4"
color_pval = "#d62728"

ax1.plot(start_ranks + 1, window_corrs, color=color_corr, marker="o", markersize=5, lw=2, label="Spearman ρ with top-61")
ax1.set_xlabel("Starting gene rank", fontsize=12)
ax1.set_ylabel("Spearman ρ (window bias vs top-61 bias)", color=color_corr, fontsize=11)
ax1.tick_params(axis="y", labelcolor=color_corr)
ax1.set_ylim(-0.1, 1.05)
ax1.axhline(0, color="gray", lw=0.5, alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(start_ranks + 1, -np.log10(np.array(window_mean_pvals)),
         color=color_pval, marker="s", markersize=4, lw=1.5, ls="--",
         label=r"$-\log_{10}$(mean P-value)")
ax2.set_ylabel(r"$-\log_{10}$(mean P-value in window)", color=color_pval, fontsize=11)
ax2.tick_params(axis="y", labelcolor=color_pval)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=9, framealpha=0.8)

ax1.set_title("R3.2b — SCZ: Sliding window (61 genes) correlation decay", fontweight="bold", fontsize=12)
for spine in ["top"]:
    ax1.spines[spine].set_visible(False)
    ax2.spines[spine].set_visible(False)

fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()

idx_below = np.searchsorted(-np.array(window_corrs), -0.3)
if idx_below < len(start_ranks):
    print(f"Signal persistence: ρ > 0.3 until rank ~{start_ranks[idx_below] + 1}")
else:
    print(f"Signal persistence: ρ > 0.3 across ALL windows (min ρ = {min(window_corrs):.3f})")

# %% [markdown]
# # Reviewer Figure R3.2c: Real vs Random Gene Additions (SCZ, ASD, DDD)
#
# Does expanding each disorder's gene set beyond top-61 add signal or noise?
# For each disorder we compare:
# - **Real**: append the next ranked disease genes (62–N)
# - **Random**: append random genes to the fixed top-61 core (100 reps, 95% CI)
#
# Metric: Spearman ρ between expanded bias profile and top-61 reference.

# %%
from joblib import Parallel, delayed

# Load expanded gene weight files
LIQ_GW_full = pd.read_csv(GeneWeightDIR + "LIQ.top500.gw", index_col=0, header=None, names=["Weight"])
DDD_GW_full = pd.read_csv(GeneWeightDIR + "DDD.hc.gw", index_col=0, header=None, names=["Weight"])

print(f"SCZ: {len(GeneWeights)} genes, ASD(ID): {len(LIQ_GW_full)} genes, DDD: {len(DDD_GW_full)} genes")

# %%
n_add_values = np.arange(0, 150, 10)  # 0, 10, ..., 140 → total 61–201
n_reps = 100

disorder_configs = [
    ("SCZ",         GeneWeights,  "#ff7f0e"),
    ("ASD with ID", LIQ_GW_full,  "#1f77b4"),
    ("DDD",         DDD_GW_full,  "#2ca02c"),
]

all_results = {}

for disorder_name, gw_df, color in disorder_configs:
    print(f"\n{'='*50}")
    print(f"Processing {disorder_name} ({len(gw_df)} genes available)")

    # Reference: top-61 bias profile for this disorder
    n_ref = min(61, len(gw_df))
    ref_gw = dict(zip(gw_df.index[:n_ref], gw_df["Weight"].values[:n_ref]))
    ref_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_gw)
    ref_bias = AnnotateCTDat(ref_bias, Anno)
    ref_effect = ref_bias["EFFECT"].values
    common_idx = ref_bias.index

    # Gene pool for random: exclude top-200 of this disorder
    n_exclude = min(200, len(gw_df))
    exclude_set = set(gw_df.index[:n_exclude].astype(int))
    gene_pool = np.array([g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set])

    # Full weight vector (for transferring rank-position weights to random genes)
    all_weights = gw_df["Weight"].values

    # ---- Real ranked genes ----
    real_corrs = []
    for n_add in n_add_values:
        n_total = min(61 + n_add, len(gw_df))
        gw = dict(zip(gw_df.index[:n_total], all_weights[:n_total]))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        real_corrs.append(r)

    # ---- Random gene additions (weight-transferred) ----
    # Random genes inherit the weight of the real gene at the same rank position,
    # so the weight structure is matched and the only difference is gene identity.
    def _one_random_trial(n_add, seed, top61_idx, top61_wt, gene_pool, rank_weights):
        rng = np.random.default_rng(seed)
        rand_genes = rng.choice(gene_pool, size=n_add, replace=False)
        gw = dict(zip(top61_idx, top61_wt))
        gw.update(dict(zip(rand_genes, rank_weights)))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        return r

    top61_idx = gw_df.index[:n_ref]
    top61_wt = gw_df["Weight"].values[:n_ref]

    rand_mean = [1.0]
    rand_lo = [1.0]
    rand_hi = [1.0]

    for n_add in n_add_values[1:]:
        # Weights for ranks 62 to 61+n_add (transfer real rank-position weights)
        n_total = min(61 + n_add, len(gw_df))
        rank_weights = all_weights[61:n_total]
        # If n_add exceeds available genes, pad with the last available weight
        if len(rank_weights) < n_add:
            rank_weights = np.concatenate([
                rank_weights,
                np.full(n_add - len(rank_weights), rank_weights[-1] if len(rank_weights) > 0 else 1.0)
            ])
        rs = Parallel(n_jobs=10)(
            delayed(_one_random_trial)(
                n_add, seed=hash((disorder_name, n_add, rep)) % (2**31),
                top61_idx=top61_idx, top61_wt=top61_wt,
                gene_pool=gene_pool, rank_weights=rank_weights,
            )
            for rep in range(n_reps)
        )
        rs = np.array(rs)
        rand_mean.append(rs.mean())
        rand_lo.append(np.percentile(rs, 2.5))
        rand_hi.append(np.percentile(rs, 97.5))

    all_results[disorder_name] = {
        "real": np.array(real_corrs),
        "rand_mean": np.array(rand_mean),
        "rand_lo": np.array(rand_lo),
        "rand_hi": np.array(rand_hi),
        "color": color,
    }
    print(f"  Done. Real ρ at N=201: {real_corrs[-1]:.3f}, Random mean: {rand_mean[-1]:.3f}")

# %%
total_genes = 61 + n_add_values

fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=120, sharey=True)

for ax, (disorder_name, gw_df, color) in zip(axes, disorder_configs):
    res = all_results[disorder_name]

    # Random: mean + 95% CI band
    ax.fill_between(total_genes, res["rand_lo"], res["rand_hi"],
                    color="#999999", alpha=0.25, label="Random genes (95% CI)")
    ax.plot(total_genes, res["rand_mean"], color="#999999", lw=2, ls="--",
            marker="s", markersize=3, label="Random genes (mean)")

    # Real ranked genes
    ax.plot(total_genes, res["real"], color=color, lw=2.5,
            marker="o", markersize=5, label=f"Ranked {disorder_name} genes", zorder=5)

    ax.axvline(61, color="gray", ls=":", lw=1.5, alpha=0.5)
    ax.set_xlabel("Total number of genes", fontsize=12)
    ax.set_title(disorder_name, fontweight="bold", fontsize=14)
    ax.legend(fontsize=8, framealpha=0.8, loc="lower left")
    ax.set_ylim(0.3, 1.02)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

axes[0].set_ylabel("Spearman ρ with top-61 bias profile", fontsize=12)

fig.suptitle("R3.2c — Expanding gene sets: real ranked genes vs random additions",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()

# Print summary for all disorders
for disorder_name, _, _ in disorder_configs:
    res = all_results[disorder_name]
    print(f"\n=== {disorder_name}: Correlation with top-61 reference ===")
    print(f"{'N total':>8s}  {'Real':>7s}  {'Rand mean':>10s}  {'Rand 95%CI':>14s}")
    for i, n_add in enumerate(n_add_values):
        nt = 61 + n_add
        if n_add == 0:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {'—':>10s}  {'—':>14s}")
        else:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {res['rand_mean'][i]:10.3f}  "
                  f"[{res['rand_lo'][i]:.3f}, {res['rand_hi'][i]:.3f}]")

# %%
