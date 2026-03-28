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
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
sys.path.insert(1, '/home/jw3514/Work/UNIMED/src')
from CellType_PSY import *
#from UNIMED import *

config = _cfg

#import scanpy as sc
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# %%
# Load expression matrix
expression_matrix = config['analysis_types']['Centering']
print(expression_matrix)

# %%
HumanCT_Z2_HCT = pd.read_csv(str(PROJ_DIR / expression_matrix), index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
HumanCT_Z2_HCT.shape

# %%
GeneWeightDIR = str(PROJ_DIR / "dat" / "GeneWeights") + "/"
Bias_Save_Dir = str(PROJ_DIR / "results/main_results/random/Centering/")
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
GeneDF = pd.read_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"), index_col=0)

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
_fig_scz, _, _ = plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_SCZ,
    GeneSig_top=-np.log10(np.array(GeneWeights["Pval"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main SCZ Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_SCZ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)",
    show=False,
)
_figdir = str(PROJ_DIR / "results" / "figures") + "/"
os.makedirs(_figdir, exist_ok=True)
_fig_scz.savefig(_figdir + "gene_sweep_SCZ.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# # ASD

# %%
# ASD references will be computed inline after weights are built (below)

# %%
Spark_Denovo = pd.read_excel(str(PROJ_DIR / "dat/suppl.data/41588_2022_1148_MOESM4_ESM.xlsx"),
                           skiprows=2, sheet_name="Table S7")
Spark_Denovo = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]!="."]
Spark_Denovo.shape

# %%
Mut_n_IQ = pd.read_csv(str(PROJ_DIR / "dat/ASD_IQ_Mut.csv"))
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
HIQ_GW_bgmr = pd.read_csv(f"{GeneWeightDIR}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
                           index_col=0, header=None, names=["Weight"])
LIQ_GW_bgmr = pd.read_csv(f"{GeneWeightDIR}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
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
_fig_hiq, _, _ = plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_HIQ,
    GeneSig_top=-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main ASD Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_HIQ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)",
    show=False,
)
_fig_hiq.savefig(_figdir + "gene_sweep_ASD_woID.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

_fig_liq, _, _ = plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_LIQ,
    GeneSig_top=-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main ASD ID Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_LIQ_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)",
    show=False,
)
_fig_liq.savefig(_figdir + "gene_sweep_ASD_wID.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

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
_fig_ddd, _, _ = plot_gene_set_correlation(
    GeneIdx_top=GeneIdx,
    Corr_top=Corr_DDD,
    GeneSig_top=-np.log10(np.array(DDD_Genes["denovoWEST_p_full"].values[GeneIdx], dtype=np.float64)),
    corr_label="Bias Correlation",
    ylabel_corr = "Bias Correlation with Main DDD Set",
    pval_label=r"$-\log_{10}$(max p-value of Included Genes)",
    Corr_unweighted=Corr_DDD_unweighted,
    corr_unweighted_label="Bias Correlation (Unweighted)",
    show=False,
)
_fig_ddd.savefig(_figdir + "gene_sweep_DDD.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

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
GeneDF = pd.read_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"), index_col=0)

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

# Build gene weight DataFrames for R3.2c using the SAME gene ranking as sweeps (Panels A-D).
# Panel C/D sweeps use Spark_Denovo_sub (ASD) and DDD_Genes_sub (DDD) ranked by significance,
# with BGMR-corrected weights for top-61 and raw counts beyond.
# We build identical DataFrames here for consistency.

# ASD with ID: from Spark_Denovo_sub, same ranking as Panel C sweep
LIQ_GW_full = pd.DataFrame({
    "Weight": Spark_Denovo_sub["LIQ_weight_bgmr"].values
}, index=Spark_Denovo_sub.index)
LIQ_GW_full["Weight"] = LIQ_GW_full["Weight"].clip(lower=0)

# DDD: from DDD_Genes_sub, same ranking as Panel D sweep
# Compute BGMR-corrected weights for ALL genes (not just top-61)
_ddd_all_gw = _ddd_gene_weights(DDD_Genes_sub, Nproband=N_DDD, BGMR=BGMR)
DDD_GW_full = pd.DataFrame([
    {"EntrezID": row["EntrezID"], "Weight": _ddd_all_gw.get(row["EntrezID"], 0)}
    for _, row in DDD_Genes_sub.iterrows()
    if row["EntrezID"] > 0
]).set_index("EntrezID")
DDD_GW_full["Weight"] = DDD_GW_full["Weight"].clip(lower=0)

print(f"SCZ: {len(GeneWeights)} genes, ASD(ID): {len(LIQ_GW_full)} genes, DDD: {len(DDD_GW_full)} genes")
print("(Same gene ranking and weights as sweep Panels A-D)")

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

    # Full weight vector — only use genes with nonzero weight for expansion.
    # Zero-weight genes (no mutations) contribute nothing to bias but inflate n_add,
    # causing random noise (which averages to ~0 on mean-centered data) to appear
    # to "preserve" correlation better than real signal.
    all_weights = gw_df["Weight"].values

    # Build list of (index, weight) for genes beyond top-61 with nonzero weight
    added_pool = [(gw_df.index[i], all_weights[i])
                  for i in range(n_ref, len(gw_df)) if all_weights[i] > 0]
    added_indices = [x[0] for x in added_pool]
    added_weights = np.array([x[1] for x in added_pool])
    print(f"  Genes beyond top-{n_ref} with nonzero weight: {len(added_pool)}")

    # ---- Real ranked genes (only nonzero-weight) ----
    real_corrs = []
    for n_add in n_add_values:
        n_use = min(n_add, len(added_pool))
        gw = dict(zip(gw_df.index[:n_ref], all_weights[:n_ref]))
        gw.update(dict(zip(added_indices[:n_use], added_weights[:n_use])))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        real_corrs.append(r)

    # ---- Random gene additions (matched nonzero weights) ----
    # Random genes get the same weights as the real nonzero-weight genes,
    # so the only difference is gene identity.
    def _one_random_trial(n_use, seed, top61_idx, top61_wt, gene_pool, rank_weights):
        rng = np.random.default_rng(seed)
        rand_genes = rng.choice(gene_pool, size=n_use, replace=False)
        gw = dict(zip(top61_idx, top61_wt))
        gw.update(dict(zip(rand_genes, rank_weights[:n_use])))
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
        n_use = min(n_add, len(added_pool))
        if n_use == 0:
            rand_mean.append(1.0)
            rand_lo.append(1.0)
            rand_hi.append(1.0)
            continue
        rs = Parallel(n_jobs=10)(
            delayed(_one_random_trial)(
                n_use, seed=hash((disorder_name, n_add, rep)) % (2**31),
                top61_idx=top61_idx, top61_wt=top61_wt,
                gene_pool=gene_pool, rank_weights=added_weights,
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
_figdir = str(PROJ_DIR / "results" / "figures") + "/"
os.makedirs(_figdir, exist_ok=True)
fig.savefig(_figdir + "FigSX_gene_expansion_R3.2c.pdf",
            dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
fig.savefig(_figdir + "FigSX_gene_expansion_R3.2c.png",
            dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
print(f"Saved: {_figdir}FigSX_gene_expansion_R3.2c.pdf/png")
plt.show()

# %% [markdown]
# ## R3.2c-slice: Incremental Gene Correlation (Slice Only, No Top-61 Core)
#
# The combined-set figure above is confounded by a dilution-vs-distortion effect
# on mean-centered data: random genes contribute ~zero signal, so the combined
# bias is just the top-61 scaled down (preserving rank → high correlation).
#
# This panel computes bias from ONLY the incremental slice (no top-61 core)
# and correlates with the top-61 reference. This directly tests whether the
# added genes carry concordant cell-type signal.

# %%
# Compute slice-only correlations for all disorders
slice_results = {}

for disorder_name, gw_df, color in disorder_configs:
    print(f"\n{'='*50}")
    print(f"Slice analysis: {disorder_name}")

    n_ref = min(61, len(gw_df))
    ref_gw = dict(zip(gw_df.index[:n_ref], gw_df["Weight"].values[:n_ref]))
    ref_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_gw)
    ref_bias = AnnotateCTDat(ref_bias, Anno)
    ref_effect = ref_bias["EFFECT"].values
    common_idx = ref_bias.index

    all_weights = gw_df["Weight"].values
    added_pool = [(gw_df.index[i], all_weights[i])
                  for i in range(n_ref, len(gw_df)) if all_weights[i] > 0]
    added_indices = [x[0] for x in added_pool]
    added_weights = np.array([x[1] for x in added_pool])

    n_exclude = min(200, len(gw_df))
    exclude_set = set(gw_df.index[:n_exclude].astype(int))
    gene_pool = np.array([g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set])

    # ---- Real slice only (no top-61 core) ----
    real_slice_corrs = [np.nan]  # n_add=0 → no genes → undefined
    for n_add in n_add_values[1:]:
        n_use = min(n_add, len(added_pool))
        if n_use == 0:
            real_slice_corrs.append(np.nan)
            continue
        gw = dict(zip(added_indices[:n_use], added_weights[:n_use]))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        real_slice_corrs.append(r)

    # ---- Random slice only (no top-61 core) ----
    def _one_random_slice_trial(n_use, seed, gene_pool, rank_weights, ref_effect, common_idx):
        rng = np.random.default_rng(seed)
        rand_genes = rng.choice(gene_pool, size=n_use, replace=False)
        gw = dict(zip(rand_genes, rank_weights[:n_use]))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        return r

    rand_slice_mean = [np.nan]
    rand_slice_lo = [np.nan]
    rand_slice_hi = [np.nan]

    for n_add in n_add_values[1:]:
        n_use = min(n_add, len(added_pool))
        if n_use == 0:
            rand_slice_mean.append(np.nan)
            rand_slice_lo.append(np.nan)
            rand_slice_hi.append(np.nan)
            continue
        rs = Parallel(n_jobs=10)(
            delayed(_one_random_slice_trial)(
                n_use, seed=hash((disorder_name, "slice", n_add, rep)) % (2**31),
                gene_pool=gene_pool, rank_weights=added_weights,
                ref_effect=ref_effect, common_idx=common_idx,
            )
            for rep in range(n_reps)
        )
        rs = np.array(rs)
        rand_slice_mean.append(rs.mean())
        rand_slice_lo.append(np.percentile(rs, 2.5))
        rand_slice_hi.append(np.percentile(rs, 97.5))

    slice_results[disorder_name] = {
        "real": np.array(real_slice_corrs, dtype=float),
        "rand_mean": np.array(rand_slice_mean, dtype=float),
        "rand_lo": np.array(rand_slice_lo, dtype=float),
        "rand_hi": np.array(rand_slice_hi, dtype=float),
        "color": color,
    }
    print(f"  Real slice at N=140: rho={real_slice_corrs[-1]:.3f}, "
          f"Random mean: {rand_slice_mean[-1]:.3f}")

# %%
# Plot slice-only panels
fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=120, sharey=True)

for ax, (disorder_name, gw_df, color) in zip(axes, disorder_configs):
    res = slice_results[disorder_name]
    valid = ~np.isnan(res["real"])
    x = n_add_values[valid]

    # Random: mean + 95% CI band
    ax.fill_between(x, res["rand_lo"][valid], res["rand_hi"][valid],
                    color="#999999", alpha=0.25, label="Random genes (95% CI)")
    ax.plot(x, res["rand_mean"][valid], color="#999999", lw=2, ls="--",
            marker="s", markersize=3, label="Random genes (mean)")

    # Real ranked genes
    ax.plot(x, res["real"][valid], color=color, lw=2.5,
            marker="o", markersize=5, label=f"Ranked {disorder_name} genes", zorder=5)

    ax.axhline(0, color="gray", ls=":", lw=1, alpha=0.5)
    ax.set_xlabel("Number of added genes (slice only)", fontsize=12)
    ax.set_title(disorder_name, fontweight="bold", fontsize=14)
    ax.legend(fontsize=8, framealpha=0.8, loc="lower right")
    ax.set_ylim(-0.6, 1.02)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

axes[0].set_ylabel("Spearman ρ with top-61 bias (slice only, no core)", fontsize=12)

fig.suptitle("Incremental gene correlation with main set (slice only, no top-61 core)",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(_figdir + "FigSX_gene_expansion_slice.pdf",
            dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
fig.savefig(_figdir + "FigSX_gene_expansion_slice.png",
            dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
print(f"Saved: {_figdir}FigSX_gene_expansion_slice.pdf/png")
plt.show()

# Save plot data for Figures_Supp inline plotting
import pickle
PLOT_DATA_DIR = str(PROJ_DIR / "results" / "figures" / "plot_data") + "/"
os.makedirs(PLOT_DATA_DIR, exist_ok=True)

gene_sweep_data = {
    'GeneIdx': GeneIdx.tolist(),
    'SCZ': {
        'Corr': Corr_SCZ,
        'Corr_unweighted': Corr_SCZ_unweighted,
        'GeneSig': (-np.log10(np.array(GeneWeights["Pval"].values[GeneIdx], dtype=np.float64))).tolist(),
        'ylabel_corr': "Bias Correlation with Main SCZ Set",
    },
    'ASD_woID': {
        'Corr': Corr_HIQ,
        'Corr_unweighted': Corr_HIQ_unweighted,
        'GeneSig': (-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64))).tolist(),
        'ylabel_corr': "Bias Correlation with Main ASD Set",
    },
    'ASD_wID': {
        'Corr': Corr_LIQ,
        'Corr_unweighted': Corr_LIQ_unweighted,
        'GeneSig': (-np.log10(np.array(Spark_Denovo_sub["pDenovoWEST_Meta"].values[GeneIdx], dtype=np.float64))).tolist(),
        'ylabel_corr': "Bias Correlation with Main ASD ID Set",
    },
    'DDD': {
        'Corr': Corr_DDD,
        'Corr_unweighted': Corr_DDD_unweighted,
        'GeneSig': (-np.log10(np.array(DDD_Genes["denovoWEST_p_full"].values[GeneIdx], dtype=np.float64))).tolist(),
        'ylabel_corr': "Bias Correlation with Main DDD Set",
    },
}

gene_expansion_data = {
    'n_add_values': n_add_values.tolist(),
    'total_genes': (61 + n_add_values).tolist(),
    'disorders': {},
}
for disorder_name, gw_df, color in disorder_configs:
    res = all_results[disorder_name]
    gene_expansion_data['disorders'][disorder_name] = {
        'real': res['real'].tolist(),
        'rand_mean': res['rand_mean'].tolist(),
        'rand_lo': res['rand_lo'].tolist(),
        'rand_hi': res['rand_hi'].tolist(),
        'color': res['color'],
    }

with open(PLOT_DATA_DIR + "gene_sweep_data.pkl", "wb") as f:
    pickle.dump(gene_sweep_data, f)
with open(PLOT_DATA_DIR + "gene_expansion_data.pkl", "wb") as f:
    pickle.dump(gene_expansion_data, f)
print(f"Saved plot data to {PLOT_DATA_DIR}")

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

# %% [markdown]
# # DDD Verification: Real Top-200 vs Random-140 Permutation Test
#
# The R3.2c figure shows DDD real top-200 genes don't clearly exceed random additions.
# Here we run a direct permutation test:
# 1. Keep the top-61 DDD genes with their real BGMR-corrected weights
# 2. For the real set: add genes ranked 62-201 with their real weights → correlation with top-61 ref
# 3. For random draws: replace genes 62-201 with 140 random genes, using the same weight vector
# 4. Repeat 100 times → empirical p-value = fraction of random draws ≥ real correlation

# %%
# --- DDD real top-201 correlation ---
# Reuse DDD_GW_full built in R3.2c above
n_ref = 61
n_add = 140
n_reps_perm = 100

# Reference: top-61 DDD bias
ref_gw_ddd = dict(zip(DDD_GW_full.index[:n_ref], DDD_GW_full["Weight"].values[:n_ref]))
ref_bias_ddd = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_gw_ddd)
ref_bias_ddd = AnnotateCTDat(ref_bias_ddd, Anno)
ref_effect_ddd = ref_bias_ddd["EFFECT"].values
common_idx_ddd = ref_bias_ddd.index

# Real expansion: top-61 + ranked 62-201 (nonzero-weight only, same as R3.2c)
all_weights_ddd = DDD_GW_full["Weight"].values
added_pool_ddd = [(DDD_GW_full.index[i], all_weights_ddd[i])
                  for i in range(n_ref, len(DDD_GW_full)) if all_weights_ddd[i] > 0]
added_indices_ddd = [x[0] for x in added_pool_ddd]
added_weights_ddd = np.array([x[1] for x in added_pool_ddd])
n_real_added = min(n_add, len(added_pool_ddd))

print(f"DDD genes beyond top-61 with nonzero weight: {len(added_pool_ddd)}")
print(f"Using {n_real_added} real genes for expansion")

# Compute real top-(61+n_real_added) bias
real_gw = dict(zip(DDD_GW_full.index[:n_ref], all_weights_ddd[:n_ref]))
real_gw.update(dict(zip(added_indices_ddd[:n_real_added], added_weights_ddd[:n_real_added])))
real_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, real_gw)
real_bias = AnnotateCTDat(real_bias, Anno)
real_r, _ = spearmanr(ref_effect_ddd, real_bias.loc[common_idx_ddd, "EFFECT"].values)
print(f"\nReal top-{61 + n_real_added} DDD correlation with top-61: rho = {real_r:.4f}")

# --- Random permutation draws ---
# Gene pool: all genes in expression matrix except top-200 of DDD
n_exclude = min(200, len(DDD_GW_full))
exclude_set_ddd = set(DDD_GW_full.index[:n_exclude].astype(int))
gene_pool_ddd = np.array([g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set_ddd])

top61_idx_ddd = DDD_GW_full.index[:n_ref]
top61_wt_ddd = DDD_GW_full["Weight"].values[:n_ref]

def _one_ddd_random_trial(rep, n_use, seed):
    rng = np.random.default_rng(seed)
    rand_genes = rng.choice(gene_pool_ddd, size=n_use, replace=False)
    gw = dict(zip(top61_idx_ddd, top61_wt_ddd))
    gw.update(dict(zip(rand_genes, added_weights_ddd[:n_use])))
    bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
    bias = AnnotateCTDat(bias, Anno)
    r, _ = spearmanr(ref_effect_ddd, bias.loc[common_idx_ddd, "EFFECT"].values)
    return r

rand_rs = Parallel(n_jobs=10)(
    delayed(_one_ddd_random_trial)(rep, n_real_added, seed=42 + rep)
    for rep in range(n_reps_perm)
)
rand_rs = np.array(rand_rs)

# Empirical p-value: fraction of random ≥ real
p_empirical = np.mean(rand_rs >= real_r)
print(f"\n{'='*60}")
print(f"DDD Permutation Test: Top-{61+n_real_added} Real vs 61+{n_real_added} Random")
print(f"{'='*60}")
print(f"Real correlation:      rho = {real_r:.4f}")
print(f"Random mean:           rho = {rand_rs.mean():.4f} +/- {rand_rs.std():.4f}")
print(f"Random [2.5%, 97.5%]:  [{np.percentile(rand_rs, 2.5):.4f}, {np.percentile(rand_rs, 97.5):.4f}]")
print(f"Fraction random >= real: {p_empirical:.2f} ({int(np.sum(rand_rs >= real_r))}/{n_reps_perm})")
if p_empirical > 0.05:
    print("-> Real top-200 DDD genes do NOT significantly exceed random additions")
else:
    print("-> Real top-200 DDD genes significantly exceed random additions")

# %%
# Histogram of random correlations vs real
fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
ax.hist(rand_rs, bins=20, color="#2ca02c", alpha=0.6, edgecolor="white", label="Random draws (n=100)")
ax.axvline(real_r, color="red", lw=2.5, ls="--", label=f"Real top-{61+n_real_added} (rho={real_r:.3f})")
ax.axvline(rand_rs.mean(), color="gray", lw=1.5, ls=":", label=f"Random mean (rho={rand_rs.mean():.3f})")
ax.set_xlabel("Spearman rho with top-61 DDD reference", fontsize=12)
ax.set_ylabel("Count", fontsize=12)
ax.set_title(f"DDD: Real top-{61+n_real_added} vs 61 + {n_real_added} random genes\n"
             f"p(random >= real) = {p_empirical:.2f}", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, framealpha=0.8)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()

# %% [markdown]
# # DDD: Isolated Slice 62-200 Bias vs Top-61
#
# Key question: does the ADDED slice (genes 62-200) itself carry concordant signal?
# Compare bias from ONLY genes 62-200 vs bias from ONLY 140 random genes.

# %%
# Load DDD gene weights (now p-value ordered)
DDD_GW_top285 = pd.read_csv(f"{GeneWeightDIR}/DDD.top285.gw.bgmr.csv",
                           index_col=0, header=None, names=["Weight"])

# Top-61 reference bias
DDD_MainBias = pd.read_csv(
    str(PROJ_DIR / "results/main_results/random/Centering/DDD_61_bias_addP.csv"),
    index_col=0)





# %%
# --- Real slice 62-200: bias from ONLY those genes ---
real_slice_gw = DDD_GW_top285.iloc[60:200]["Weight"].to_dict()
real_slice_gw_nz = {k: v for k, v in real_slice_gw.items() if v > 0}
real_slice_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, real_slice_gw_nz)
real_slice_bias = AnnotateCTDat(real_slice_bias, Anno)
real_slice_r, real_slice_p = GetSingeCellBiasCorr(real_slice_bias, DDD_MainBias)

print(f"Real genes 62-200 ONLY (n={len(real_slice_gw_nz)} nonzero-weight):")
print(f"  Bias correlation with top-61: rho = {real_slice_r:.4f}, p = {real_slice_p:.2e}")

# --- Random 140 genes: bias from ONLY random genes (no top-61 core) ---
n_rand = len(real_slice_gw_nz)
exclude_set_ddd2 = set(DDD_GW_top285.index[:200].astype(int))
gene_pool_ddd2 = np.array([g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set_ddd2])
rank_weights = np.array(list(real_slice_gw_nz.values()))

def _rand_slice_trial(rep, n_use, rank_weights, gene_pool):
    rng = np.random.default_rng(42 + rep)
    rand_genes = rng.choice(gene_pool, size=n_use, replace=False)
    gw = dict(zip(rand_genes, rank_weights[:n_use]))
    bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
    bias = AnnotateCTDat(bias, Anno)
    r, _ = GetSingeCellBiasCorr(bias, DDD_MainBias)
    return r

rand_slice_rs = Parallel(n_jobs=10)(
    delayed(_rand_slice_trial)(rep, n_rand, rank_weights, gene_pool_ddd2)
    for rep in range(100)
)
rand_slice_rs = np.array(rand_slice_rs)

p_slice = np.mean(rand_slice_rs >= real_slice_r)
print(f"\nRandom {n_rand} genes ONLY (no top-61 core):")
print(f"  Mean rho = {rand_slice_rs.mean():.4f} +/- {rand_slice_rs.std():.4f}")
print(f"  95% CI: [{np.percentile(rand_slice_rs, 2.5):.4f}, {np.percentile(rand_slice_rs, 97.5):.4f}]")
print(f"\nFraction random >= real: {p_slice:.2f} ({int(np.sum(rand_slice_rs >= real_slice_r))}/100)")
if real_slice_r > np.percentile(rand_slice_rs, 97.5):
    print("-> Real slice 62-200 carries SIGNIFICANT concordant signal with top-61")
elif real_slice_r > rand_slice_rs.mean():
    print("-> Real slice 62-200 trends above random but NOT significant")
else:
    print("-> Real slice 62-200 does NOT exceed random genes")

# %%
# Histogram: isolated slice comparison
fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
ax.hist(rand_slice_rs, bins=20, color="#2ca02c", alpha=0.6, edgecolor="white",
        label=f"Random {n_rand} genes (n=100 draws)")
ax.axvline(real_slice_r, color="red", lw=2.5, ls="--",
           label=f"Real genes 62-200 (rho={real_slice_r:.3f})")
ax.axvline(rand_slice_rs.mean(), color="gray", lw=1.5, ls=":",
           label=f"Random mean (rho={rand_slice_rs.mean():.3f})")
ax.set_xlabel("Spearman rho with top-61 DDD bias (slice ONLY, no core)", fontsize=12)
ax.set_ylabel("Count", fontsize=12)
ax.set_title(f"DDD: Isolated slice 62-200 vs {n_rand} random genes\n"
             f"p(random >= real) = {p_slice:.2f}", fontsize=12, fontweight="bold")
ax.legend(fontsize=9, framealpha=0.8)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()
