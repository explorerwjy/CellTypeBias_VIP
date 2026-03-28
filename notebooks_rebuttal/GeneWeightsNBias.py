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
sys.path.insert(1, str(PROJ_DIR / "src"))
sys.path.insert(1, '/home/jw3514/Work/UNIMED/src')
from CellType_PSY import *
#from UNIMED import *
#import scanpy as sc
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# %%
HumanCT_Spec = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv", index_col=0)
HumanCT_Spec.columns = HumanCT_Spec.columns.astype(int)
HumanCT_Spec.shape


# %%
#Functions:
def plot_vip_effect_comparison(df, effect_col = "EFFECT", plot=True):
    """
    Plot comparison of effects between VIP+ and VIP- cells

    Parameters:
    df: DataFrame with columns 'VIP' and 'EFFECT'
    plot: bool, whether to produce a plot (default True)

    Returns:
    pval: float, p-value from Mann-Whitney U test
    """
    vip_pos = df[df["VIP"] >= 1]
    vip_neg = df[df["VIP"] < 1]
    
    # Prepare data for plotting or stats
    data = [vip_pos[effect_col], vip_neg[effect_col]]
    
    # Perform Mann-Whitney U test
    stat, pval = scipy.stats.mannwhitneyu(vip_pos[effect_col], 
                                          vip_neg[effect_col])
    
    if plot:
        plt.boxplot(data, labels=["VIP+", "VIP-"], showfliers=True)
        plt.plot([1] * len(vip_pos[effect_col]), 
                 vip_pos[effect_col], 'ko', alpha=0.3)
        plt.plot([2] * len(vip_neg[effect_col]), 
                 vip_neg[effect_col], 'ko', alpha=0.3)
        plt.ylabel(effect_col)
        plt.title(f"p = {pval:.2e}")
        plt.show()
    return pval

VIP_Anno = pd.read_csv(PROJ_DIR / "notebooks" / "VIP_Anno.csv", index_col=0)

# %%
GeneWeightDIR = str(PROJ_DIR / "dat" / "GeneWeights") + "/"
HIQ_GW = Fil2Dict("{}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))
LIQ_GW = Fil2Dict("{}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))
SCZ_GW = Fil2Dict("{}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(GeneWeightDIR))

# %%
SCZ_Genes = list(SCZ_GW.keys())
SCZ_Genes_valid = [g for g in SCZ_Genes if g in HumanCT_Spec.index]
SCZ_Spec = HumanCT_Spec.loc[SCZ_Genes_valid, :]

# %%
superclusters = list(Anno.Supercluster.unique())


# %%

# %%
def plot_gene_spec_across_superclusters(Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno):
    """
    Plots the specificity value for a given gene across superclusters, with boxplots sorted by median.

    Parameters:
    - Gene_Sepc: Series, specificity values indexed by cell type for the gene.
    - Gene: str, gene symbol.
    - Gene_Entrez: int or str, Entrez ID of the gene.
    - superclusters: list, unique supercluster names.
    - Anno: DataFrame, annotation with 'Supercluster' and cell type indices.
    """
    # Collect spec values and their medians for each supercluster
    supercluster_data = []
    for supercluster in superclusters:
        cts = Anno.loc[Anno.Supercluster == supercluster].index
        spec = Gene_Sepc[cts]
        median_spec = np.median(spec.values)
        supercluster_data.append((supercluster, spec.values, median_spec))

    # Sort by median
    supercluster_data_sorted = sorted(supercluster_data, key=lambda x: x[2], reverse=True)

    supercluster_labels = [item[0] for item in supercluster_data_sorted]
    supercluster_specs = [item[1] for item in supercluster_data_sorted]

    plt.figure(figsize=(10, 6))
    plt.boxplot(supercluster_specs, labels=supercluster_labels, notch=True)
    plt.ylabel(f"Spec value for {Gene} ({Gene_Entrez})")
    plt.xlabel("Supercluster")
    plt.title(f"{Gene} spec value across superclusters (sorted by median)")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()


# %%
Gene = "HCN4"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
Gene = "SCN1A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

Gene = "SCN2A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

Gene = "SCN8A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

Gene = "SCN3A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
Gene = "SCN9A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
Gene = "KCNV1"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]
plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
Gene = "GRIN2A"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
Gene = "GRIA3"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)
common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
pval =plot_vip_effect_comparison(df, plot=True)

# %%
Gene = "SOBP"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
# Align and prepare a DataFrame for plotting with Gene_Sepc (which is likely a Series of Spec values for a single gene, indexed by cell type or sample)

# Ensure the index of Gene_Sepc and VIP_Anno can be merged/joined
common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
# Create a DataFrame with the Spec values for the common indices
df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
# Add the VIP annotation as a column
df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
pval =plot_vip_effect_comparison(df, plot=False)
#plot_vip_effect_comparison(df, effect_col = "EFFECT_Adj")
# https://www.cell.com/ajhg/fulltext/S0002-9297(10)00521-5

# %%
sig_genes = []
for Gene in SCZ_Genes_valid:
    try:
        Gene_Name = Entrez2Symbol[Gene]
        Gene_Sepc = HumanCT_Spec.loc[Gene, :]
        common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
        df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
        df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
        pval = plot_vip_effect_comparison(df, plot=False)
        if pval < 0.05:
            print(f"{Gene_Name} has a significant effect on VIP status (p = {pval:.3g})")
            sig_genes.append(Gene_Name)
    except:
        print(f"{Gene_Name} has no valid Entrez ID")
# After the loop, print the significant genes
print("\nSignificant genes:", sig_genes)


# %%
for gene in sig_genes:
    print(gene)

# %%
Gene = "HCN4"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)
common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
pval =plot_vip_effect_comparison(df, plot=True)

# %%

# %%

# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

# Prepare data for heatmap
# Transpose so genes are columns (x-axis) and cell types are rows (y-axis)
heatmap_data = SCZ_Spec.T

# Cluster genes using hierarchical clustering
# Use SCZ_Spec (genes as rows) for clustering
gene_linkage = linkage(SCZ_Spec, method='ward', metric='euclidean')
gene_order = leaves_list(gene_linkage)

# Cluster cell types using hierarchical clustering
# Use heatmap_data (cell types as rows, genes as columns) for clustering
cell_type_linkage = linkage(heatmap_data, method='ward', metric='euclidean')
cell_type_order = leaves_list(cell_type_linkage)

# Get supercluster information for each cell type
cell_types = heatmap_data.index.astype(int)
supercluster_labels = Anno.loc[cell_types, 'Supercluster'].values

# Reorder cell types according to clustering
heatmap_data_clustered_ct = heatmap_data.iloc[cell_type_order, :]

# Get supercluster labels in clustered order
supercluster_labels_clustered = supercluster_labels[cell_type_order]

# Create y-axis labels: only show supercluster names at the start of each block
y_labels = [''] * len(heatmap_data_clustered_ct)
current_sc = None
for i, sc in enumerate(supercluster_labels_clustered):
    if sc != current_sc:
        y_labels[i] = sc
        current_sc = sc

# Reorder genes (columns) according to clustering
heatmap_data_clustered = heatmap_data_clustered_ct.iloc[:, gene_order]

# Convert gene Entrez IDs to symbols for x-axis labels (in clustered order)
gene_indices_clustered = [SCZ_Spec.index[i] for i in gene_order]
gene_symbols = [Entrez2Symbol.get(int(g), str(g)) for g in gene_indices_clustered]

# Create the heatmap
#plt.figure(figsize=(max(12, len(gene_symbols) * 0.3), max(8, len(cell_types) * 0.15)))
plt.figure(figsize=(15,12))
sns.heatmap(heatmap_data_clustered, 
            xticklabels=gene_symbols,
            yticklabels=y_labels,
            cmap='RdBu_r', 
            center=0,
            cbar_kws={'label': 'SCZ_Spec'},
            linewidths=0)

plt.xlabel('Genes', fontsize=12)
plt.ylabel('Cell Types (labeled by supercluster)', fontsize=12)
plt.title('SCZ_Spec Heatmap: Genes vs Cell Types', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=8)
plt.tight_layout()
plt.show()


# %%
for Gene in ["DRD1", "DRD2", "DRD3", "DRD4", "DRD5"]:
    Gene_Entrez = GeneSymbol2Entrez[Gene]
    Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

    plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)

# %%
for g in ["DRD1", "DRD2", "DRD3", "DRD4"]:
    entrez = GeneSymbol2Entrez[g]
    Gene_Spec = HumanCT_Spec.loc[entrez, :]
    
    # Get CGE interneuron subtypes from Anno
    CGE_Anno = Anno[Anno["Supercluster"] == "CGE interneuron"]
    
    # Find subtypes with specificity > 0.5
    high_spec_subtypes = []
    for ct_idx in CGE_Anno.index:
        if ct_idx in Gene_Spec.index:
            spec_value = Gene_Spec.loc[ct_idx]
            if spec_value > 0.5:
                subtype = CGE_Anno.loc[ct_idx, "Subtype auto-annotation"]
                high_spec_subtypes.append((ct_idx, subtype, spec_value))
    
    # Print results
    print(f"\n{g} (Entrez: {entrez}):")
    if high_spec_subtypes:
        # Sort by specificity value (descending)
        high_spec_subtypes.sort(key=lambda x: x[2], reverse=True)
        for ct_idx, subtype, spec_val in high_spec_subtypes:
            print(f"  Cluster ID {ct_idx} ({subtype}): {spec_val:.4f}")
    else:
        print("  No CGE interneuron subtypes with expression specificity > 0.5")


# %%
VIP_Anno

# %%
Gene = "TBX1"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)
common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
pval =plot_vip_effect_comparison(df, plot=True)

# %%
Gene = "DGCR8"
Gene_Entrez = GeneSymbol2Entrez[Gene]
Gene_Sepc = HumanCT_Spec.loc[Gene_Entrez, :]

plot_gene_spec_across_superclusters( Gene_Sepc, Gene, Gene_Entrez, superclusters, Anno)
common_indices = Gene_Sepc.index.intersection(VIP_Anno.index)
df = pd.DataFrame({'EFFECT': Gene_Sepc.loc[common_indices]})
df["VIP"] = VIP_Anno.loc[common_indices, "VIP"].values
pval =plot_vip_effect_comparison(df, plot=True)

# %%
