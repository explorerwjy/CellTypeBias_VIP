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
Bias_Save_Dir = ProjDIR + "results/random/Centering/"
if not os.path.exists(Bias_Save_Dir): # make dir if not exists
    os.makedirs(Bias_Save_Dir)

# %%
SCZ_Bias = pd.read_csv("{}/SCZ_bias_addP.csv".format(Bias_Save_Dir), index_col=0)

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
HIQ_Bias = pd.read_csv("{}/ASD_HIQ_bias_addP.csv".format(Bias_Save_Dir), index_col=0)
LIQ_Bias = pd.read_csv("{}/ASD_LIQ_bias_addP.csv".format(Bias_Save_Dir), index_col=0)

HIQ_Bias.head()
LIQ_Bias.head()

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
# Calculate correlation as genes are removed, sorted by LofZ (descending)
GeneIdx_HIQ, Corr_HIQ, Pval_HIQ = [], [], []
Corr_HIQ_unweighted = []
#GeneIdx_LIQ, Corr_LIQ, Pval_LIQ = [], [], []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp_HIQ_ASD_GW = dict(zip(Spark_Denovo_sub.index.values[:i], Spark_Denovo_sub["HIQ_counts"].values[:i]))
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
# Calculate correlation as genes are removed, sorted by LofZ (descending)
GeneIdx_LIQ, Corr_LIQ, Pval_LIQ = [], [], []
Corr_LIQ_unweighted = []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp_LIQ_ASD_GW = dict(zip(Spark_Denovo_sub.index.values[:i], Spark_Denovo_sub["LIQ_counts"].values[:i]))
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
DDD_Bias = pd.read_csv("{}/DDD_61_bias_addP.csv".format(Bias_Save_Dir), index_col=0)

# %%
DDD_Genes = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
DDD_Genes = DDD_Genes.sort_values("denovoWEST_p_full")
entrez_ids = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_Genes["symbol"].values]
DDD_Genes["EntrezID"] = entrez_ids
DDD_Genes.shape

# %%
DDD_Genes_sub = DDD_Genes[["symbol", "EntrezID", "denovoWEST_p_full", "missense_variant", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost"]]
#DDD_Genes_sub = 

# %%
# Calculate correlation as genes are removed, sorted by LofZ (descending)
Corr_DDD, Pval_DDD = [], []
Corr_DDD_unweighted = []
#GeneIdx_LIQ, Corr_LIQ, Pval_LIQ = [], [], []
GeneIdx = np.arange(10, 200)
for i in GeneIdx:
    tmp = DDD_Genes_sub.head(i)
    tmp_DDD_GW = Aggregate_Gene_Weights_NDD(tmp)
    tmp_DDD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, tmp_DDD_GW)
    tmp_DDD_Bias = AnnotateCTDat(tmp_DDD_Bias, Anno)
    r, p = GetSingeCellBiasCorr(DDD_Bias, tmp_DDD_Bias, efflabel="EFFECT")
        
    Corr_DDD.append(r)
    Pval_DDD.append(p)
    
    # Unweighted version (each gene has weight 1)
    tmp_DDD_GW_unweighted = dict(zip(tmp["EntrezID"].values, np.ones(len(tmp))))
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

# %%

# %%

# %%
