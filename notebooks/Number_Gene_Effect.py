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
# # Reviewer Figure R3.2c: Cross-Disorder Effect Size vs Gene Set Size
#
# At N=61, SCZ/ASD(ID)/DDD have comparable CGE effect sizes.
# Beyond 61, SCZ dilutes faster — justifying N=61 as the regime
# where cross-disorder comparisons are fair.
#
# Supports reply to Reviewer #3, Point 2, Paragraph 3 (effect size divergence).

# %%
# Load expanded gene weight files for ASD(ID) and DDD
LIQ_GW_full = pd.read_csv(GeneWeightDIR + "LIQ.top500.gw", index_col=0, header=None, names=["Weight"])
DDD_GW_full = pd.read_csv(GeneWeightDIR + "DDD.hc.gw", index_col=0, header=None, names=["Weight"])

print(f"SCZ: {len(GeneWeights)} genes, ASD(ID): {len(LIQ_GW_full)} genes, DDD: {len(DDD_GW_full)} genes")

# %%
gene_counts_sweep = [20, 30, 40, 50, 61, 80, 100, 120, 150, 175, 200]

def compute_bias_stats(gw_df, n_genes):
    n = min(n_genes, len(gw_df))
    gw = dict(zip(gw_df.index[:n], gw_df["Weight"].values[:n]))
    bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
    bias = AnnotateCTDat(bias, Anno)
    neur = bias.loc[bias.index.isin(Neur_idx)]
    cge = bias.loc[bias["Supercluster"] == "CGE interneuron", "EFFECT"]
    return {
        "n": n,
        "mean_cge": cge.mean(),
        "std_neur": neur["EFFECT"].std(),
    }

results_scz = [compute_bias_stats(GeneWeights, n) for n in gene_counts_sweep]
results_asd = [compute_bias_stats(LIQ_GW_full, n) for n in gene_counts_sweep]
results_ddd = [compute_bias_stats(DDD_GW_full, n) for n in gene_counts_sweep]

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=120)

colors = {"SCZ": "#ff7f0e", "ASD with ID": "#1f77b4", "DDD": "#2ca02c"}

# Panel A: CGE interneuron mean effect
for res_list, label in [(results_scz, "SCZ"), (results_asd, "ASD with ID"), (results_ddd, "DDD")]:
    ns = [r["n"] for r in res_list]
    vals = [r["mean_cge"] for r in res_list]
    ax1.plot(ns, vals, marker="o", markersize=5, lw=2, color=colors[label], label=label)

ax1.axvline(61, color="gray", ls="--", lw=1.5, alpha=0.6, label="N = 61")
ax1.set_xlabel("Number of genes", fontsize=12)
ax1.set_ylabel("Mean CGE interneuron effect", fontsize=12)
ax1.set_title("A", fontweight="bold", loc="left", fontsize=14)
ax1.legend(fontsize=9, framealpha=0.8)
for spine in ["top", "right"]:
    ax1.spines[spine].set_visible(False)

# Panel B: Neuronal effect std (discrimination)
for res_list, label in [(results_scz, "SCZ"), (results_asd, "ASD with ID"), (results_ddd, "DDD")]:
    ns = [r["n"] for r in res_list]
    vals = [r["std_neur"] for r in res_list]
    ax2.plot(ns, vals, marker="o", markersize=5, lw=2, color=colors[label], label=label)

ax2.axvline(61, color="gray", ls="--", lw=1.5, alpha=0.6, label="N = 61")
ax2.set_xlabel("Number of genes", fontsize=12)
ax2.set_ylabel("Std of neuronal effects (discrimination)", fontsize=12)
ax2.set_title("B", fontweight="bold", loc="left", fontsize=14)
ax2.legend(fontsize=9, framealpha=0.8)
for spine in ["top", "right"]:
    ax2.spines[spine].set_visible(False)

fig.suptitle("R3.2c — Effect size comparable at N=61, SCZ dilutes faster beyond", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.patch.set_alpha(0)
plt.show()

# Print key numbers
print("\n=== CGE effect at key gene set sizes ===")
print(f"{'N':>5s}  {'SCZ':>8s}  {'ASD(ID)':>8s}  {'DDD':>8s}")
for i, n in enumerate(gene_counts_sweep):
    if n in [61, 100, 200]:
        print(f"{n:5d}  {results_scz[i]['mean_cge']:8.4f}  {results_asd[i]['mean_cge']:8.4f}  {results_ddd[i]['mean_cge']:8.4f}")

# %%
