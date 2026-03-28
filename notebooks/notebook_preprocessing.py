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
from pathlib import Path
import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *
from tabulate import tabulate
import os
#import scanpy as sc

# %%
#Cluster_Exp_DF_log2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.csv")

# %%
Anno

# %%
import matplotlib.pyplot as plt

plt.hist(Anno["Total UMI"].dropna(), bins=50)
plt.xlabel("Total UMI")
plt.ylabel("Frequency")
plt.title("Histogram of Total UMI")
plt.show()

# %%

# %%

# %%
HumanCT_ExpR = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv", index_col=0)
HumanCT_OverallEXP = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/HumanCT.MatchDF.csv", index_col=0)


# %%
def calculate_specificity_matrix(expression_df):
    """
    Calculate cell type expression specificity directly using a matrix.
    
    Parameters:
    -----------
    expression_df : pandas DataFrame
        Expression matrix where rows are genes and columns are cell types
    
    Returns:
    --------
    pandas DataFrame
        A DataFrame with the same shape as input, containing specificity values
    """
    # Step 1: Calculate TPM normalization for each cell type (column)
    # Sum expression in each cell type
    cell_type_sums = expression_df.sum(axis=0)
    
    # Scale to TPM (multiply by 1,000,000 and divide by column sum)
    tpm_df = expression_df.div(cell_type_sums, axis=1) * 1000000
    
    # Step 2: Calculate specificity - for each gene (row), divide by the sum across cell types
    gene_total_tpm = tpm_df.sum(axis=1)
    specificity_df = tpm_df.div(gene_total_tpm, axis=0)
    
    return specificity_df

def calculate_specificity_with_filtering(expression_df, min_tpm=1):
    """
    Calculate specificity with optional TPM filtering
    
    Parameters:
    -----------
    expression_df : pandas DataFrame
        Expression matrix where rows are genes and columns are cell types
    min_tpm : float
        Minimum TPM value to keep
    
    Returns:
    --------
    tuple of (specificity_df, tpm_df)
        Both DataFrames with the same shape as input
    """
    # Calculate TPM
    cell_type_sums = expression_df.sum(axis=0)
    tpm_df = expression_df.div(cell_type_sums, axis=1) * 1000000
    
    # Create a mask for filtering (only used for specificity calculation)
    if min_tpm > 0:
        mask = tpm_df >= min_tpm
    else:
        mask = pd.DataFrame(True, index=tpm_df.index, columns=tpm_df.columns)
    
    # For specificity calculation, treat values below threshold as 0
    filtered_tpm = tpm_df.copy()
    filtered_tpm[~mask] = 0
    
    # Calculate specificity using the filtered TPM values
    gene_total_tpm = filtered_tpm.sum(axis=1)
    specificity_df = filtered_tpm.div(gene_total_tpm, axis=0)
    specificity_df = specificity_df * specificity_df.shape[1]
    
    # Handle genes that have all zeros after filtering
    specificity_df.fillna(0, inplace=True)
    
    return specificity_df, tpm_df


# %%
#specificity_df = calculate_specificity_matrix(HumanCT_ExpR)
TPM_Cut = 0.1
specificity_filt_df, tpm_df = calculate_specificity_with_filtering(HumanCT_ExpR, min_tpm=TPM_Cut)
specificity_filt_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.csv".format(TPM_Cut))
tpm_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.csv".format(TPM_Cut))

# %%
tpm_df

# %%
Mean = np.mean(specificity_filt_df.values.flatten())
#print(Mean, 1/461)
Upper = Mean * 2
specificity_filt_df_clip = specificity_filt_df.clip(lower=0, upper=Upper) # Clip max at 2x the mean
specificity_filt_df_clip.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.clip.csv".format(TPM_Cut))

# %%
fraction_zeros_per_gene = (specificity_filt_df_clip == 0).sum(axis=1) / specificity_filt_df_clip.shape[1]
fraction_zeros_per_gene

# %%
Fraction_Cut = 0.8
RobustlyExpressedGenes = fraction_zeros_per_gene[fraction_zeros_per_gene < Fraction_Cut]


# %%
RobustlyExpressedGenes

# %%
ExpL = []
Frac_Zero = []
for g in HumanCT_OverallEXP.index:
    ExpL.append(HumanCT_OverallEXP.loc[g, "Exp"])
    Frac_Zero.append((specificity_filt_df_clip.loc[g, :] == 0).sum() / specificity_filt_df_clip.shape[1])

plt.figure(figsize=(8, 5), dpi=150)
plt.scatter(ExpL, Frac_Zero, s=10, alpha=0.5, color='royalblue', edgecolor='none')
plt.xscale("log")
plt.axvline(x=1e4, color='crimson', linestyle='--', linewidth=2, label='Exp = 10,000')
plt.xlabel("Gene Expression (log scale)", fontsize=12)
plt.ylabel("Fraction of Zeros", fontsize=12)
plt.title("Gene Expression vs. Fraction of Zeros", fontsize=14, pad=12)
plt.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.4)
plt.legend(loc='upper right', fontsize=10, frameon=False)
plt.tight_layout()
plt.show()

# %%
LowExpCut = 10000
LowExpCutLog = np.log10(LowExpCut+1)
print(LowExpCutLog)
LowExpGenes = HumanCT_OverallEXP[HumanCT_OverallEXP["Exp"] < LowExpCut].index.tolist()
NGenesLowExp = len(LowExpGenes)
print(f"Number of genes with expression < {LowExpCut}: {NGenesLowExp}")

# %%
0.1 * 3000000

# %%
plt.figure(figsize=(10, 6))
plt.hist(np.log10(HumanCT_OverallEXP["Exp"]+1), bins=100, color='cornflowerblue', alpha=0.7, edgecolor='black')
# add a line at the LowExpCutLog
plt.axvline(x=LowExpCutLog, color='red', linestyle='--', linewidth=2)
plt.xlabel('Log10(Expression + 1)', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title('Distribution of Gene Expression Values', fontsize=14, pad=15)
plt.grid(True, alpha=0.3)
plt.tight_layout()

# %%
common_genes = list(set(HumanCT_OverallEXP.index.values) & set(specificity_filt_df_clip.index))
GeneExp = np.log10(HumanCT_OverallEXP.loc[common_genes, "Exp"]+1)
GeneAvgSpec = specificity_filt_df_clip.loc[common_genes, :].mean(axis=1)
print(len(GeneExp), len(GeneAvgSpec))

# %%
fig = plt.figure(dpi=200)
plt.scatter(GeneExp, GeneAvgSpec, s=0.1)
plt.xlabel("Average Expression")
plt.ylabel("Average Specificity")
plt.title("Average Expression vs Average Specificity")
#plt.ylim(0.0015, 0.0022)
plt.show()


# %%
def CellTypeMeanCorrection(SpecMat):
    for c in SpecMat.columns.values:
        mean_spec = SpecMat[c].mean()
        SpecMat[c] = SpecMat[c] - mean_spec
    return SpecMat
def CellTypeMeanCorrection(SpecMat):
    return SpecMat.subtract(SpecMat.mean(axis=0), axis=1)



# %%
print(TPM_Cut)
specificity_filt_df_clip.columns = specificity_filt_df_clip.columns.astype(int)
specificity_filt_df_clip = specificity_filt_df_clip.loc[~specificity_filt_df_clip.index.isin(LowExpGenes)]
specificity_filt_df_clip.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.clip.lowexp.cut1e4.csv".format(TPM_Cut))
specificity_filt_df_clip_centered = CellTypeMeanCorrection(specificity_filt_df_clip)
specificity_filt_df_clip_centered.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv".format(TPM_Cut))

# %%
specificity_filt_df_clip.shape

# %%
specificity_filt_df_dropLowExp = specificity_filt_df.loc[~specificity_filt_df.index.isin(LowExpGenes)]
specificity_filt_df_dropLowExp = specificity_filt_df_dropLowExp.apply(lambda x: pd.Series.rank(x, pct=True))
specificity_filt_df_dropLowExp.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.lowexp.cut1e4.Percentile.csv".format(TPM_Cut))

# %%
specificity_df_percentile = specificity_df.apply(lambda x: pd.Series.rank(x, pct=True))
specificity_filt_df_percentile = specificity_filt_df.apply(lambda x: pd.Series.rank(x, pct=True))

# %%

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
for col in specificity_filt_df_clip.columns[:100]:
    specificity_filt_df_clip[col].plot(kind='density', alpha=0.3, linewidth=1)
plt.title('Density plot for all columns')
plt.xlabel('Specificity')
plt.ylabel('Density')
plt.show()

# %%
# Group the index by Supercluster and convert to a dictionary
grouped = {name: group for name, group in Anno.groupby('Supercluster')}

# %%
grouped.keys()

# %%
grouped["Astrocyte"].index.values

# %%
avg_spec_df = pd.DataFrame({
    k: specificity_filt_df_clip.loc[:, v.index.values].mean(axis=1)
    for k, v in grouped.items()
})


# %%
avg_spec_df


# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(16, 10))  # Increased figure size
for col in avg_spec_df.columns:
    avg_spec_df[col].plot(kind='density', alpha=0.3, linewidth=1)

avg_spec_df["Astrocyte"].plot(kind='density', alpha=1, linewidth=2, label="Astrocyte")
avg_spec_df["Lower rhombic lip"].plot(kind='density', alpha=1, linewidth=2, label="Lower rhombic lip")
plt.legend(loc='upper right')
plt.title('Density plot for all columns')
plt.xlabel('Specificity')
plt.ylabel('Density')
plt.show()

# %%
plt.figure(figsize=(16, 10)) 
avg_spec_df["Astrocyte"].plot(kind='density', alpha=1, linewidth=2, label="Astrocyte")
avg_spec_df["Lower rhombic lip"].plot(kind='density', alpha=1, linewidth=2, label="Lower rhombic lip")
avg_spec_df["Choroid plexus"].plot(kind='density', alpha=1, linewidth=2, label="Choroid plexus")

for col in ["Amygdala excitatory", "MGE interneuron", "CGE interneuron"]:
    avg_spec_df[col].plot(kind='density', alpha=0.5, linewidth=1, label=col, color="grey")
plt.legend(loc='upper right')
plt.title('Density plot for all columns')
plt.xlabel('Specificity')
plt.ylabel('Density')
plt.show()

# %%
## top_spec vs UMI
UMIS = []
TopSpec = []
N_genes_gt_threshold = []
threshold = 1.5
for c in specificity_filt_df_clip.columns:
    umi = Anno.loc[int(c), "Total UMI"]
    top_spec = specificity_filt_df_clip.loc[:, c].sort_values(ascending=False).head(8000).mean()
    n_genes_gt_threshold = (specificity_filt_df_clip.loc[:, c] > threshold).sum()
    UMIS.append(umi)
    TopSpec.append(top_spec)
    N_genes_gt_threshold.append(n_genes_gt_threshold)

# %%
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

axs[0].scatter(UMIS, TopSpec)
axs[0].set_xlabel("Total UMI")
axs[0].set_ylabel("Top 8000 Specificity")
axs[0].set_title("Top 8000 Specificity vs Total UMI")

axs[1].scatter(UMIS, N_genes_gt_threshold)
axs[1].set_xlabel("Total UMI")
axs[1].set_ylabel("N genes > 1.5")
axs[1].set_title("N genes > 1.5 vs Total UMI")

plt.tight_layout()
plt.show()


# %%
def quantileNormalize(df_input):
    df = df_input.copy()
    # Sort each column and get the order of ranks for each column
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0),
                             index=df.index,
                             columns=df.columns)
    # Compute mean of each row across columns (i.e., mean for each rank)
    rank_means = sorted_df.mean(axis=1).values

    # Get the ranks for each value in each column (1-based)
    ranks = df.rank(method="min", axis=0).astype(int)

    # Create an output DataFrame of the same shape
    df_qn = df.copy()
    # Assign the mean value for each rank
    for col in df.columns:
        # ranks in this column
        col_ranks = ranks[col].values
        # assign mean for each rank (subtract 1 for 0-based index)
        df_qn[col] = [rank_means[r-1] for r in col_ranks]

    return df_qn

avg_spec_df_qn = quantileNormalize(avg_spec_df)

# %%
avg_spec_df_qn

# %%
plt.figure(figsize=(16, 10)) 
avg_spec_df_qn["Astrocyte"].plot(kind='density', alpha=1, linewidth=2, label="Astrocyte")
avg_spec_df_qn["Lower rhombic lip"].plot(kind='density', alpha=1, linewidth=2, label="Lower rhombic lip")
avg_spec_df_qn["Choroid plexus"].plot(kind='density', alpha=1, linewidth=2, label="Choroid plexus")

for col in ["Amygdala excitatory", "MGE interneuron", "CGE interneuron"]:
    avg_spec_df_qn[col].plot(kind='density', alpha=0.5, linewidth=1, label=col, color="grey")
plt.legend(loc='upper right')
plt.title('Density plot for all columns')
plt.xlabel('Specificity')
plt.ylabel('Density')
plt.show()

# %%
specificity_filt_df_clip_qn = quantileNormalize(specificity_filt_df_clip)
specificity_filt_df_clip_qn.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.{}.Filt.Spec.clip.lowexp.cut1e4.qn.csv".format(TPM_Cut))

# %%

# %%
# TPM Simulation

# %%
Anno

# %%

# %%

# %%

# %% [markdown]
# ### Subcluster

# %%
HumanCT_Subcluster = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/subcluster_MeanLogUMI.csv", index_col=0)
HumanCT_Subcluster = np.expm1(HumanCT_Subcluster)

# %%
TPM_Cut = 0.1
specificity_filt_df, tpm_df = calculate_specificity_with_filtering(HumanCT_Subcluster, min_tpm=TPM_Cut)
#specificity_filt_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Subcluster.TPM.{}.Filt.Spec.csv".format(TPM_Cut))
tpm_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.Subcluster.TPM.{}.Filt.csv".format(TPM_Cut))

# %%
Mean = np.mean(specificity_filt_df.values.flatten())
#print(Mean, 1/461)
Upper = Mean * 2
specificity_filt_df_clip = specificity_filt_df.clip(lower=0, upper=Upper) # Clip max at 2x the mean
specificity_filt_df_clip.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.Subcluster.TPM.{}.Filt.Spec.clip.csv".format(TPM_Cut))

# %%
LowExpCut = 10000
LowExpCutLog = np.log10(LowExpCut+1)
print(LowExpCutLog)
LowExpGenes = HumanCT_OverallEXP[HumanCT_OverallEXP["Exp"] < LowExpCut].index.tolist()
NGenesLowExp = len(LowExpGenes)
print(f"Number of genes with expression < {LowExpCut}: {NGenesLowExp}")

# %%
print(TPM_Cut)
specificity_filt_df_clip = specificity_filt_df_clip.loc[~specificity_filt_df_clip.index.isin(LowExpGenes)]
specificity_filt_df_clip = CellTypeMeanCorrection(specificity_filt_df_clip)
specificity_filt_df_clip.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.Subcluster.TPM.{}.Filt.Spec.clip.lowexp.cut1e4.csv".format(TPM_Cut))

# %%

# %%

# %% [markdown] heading_collapsed=true
# ## File Downloads

# %% [markdown] heading_collapsed=true
# ## Example: Make subclass Exp Mat

# %% hidden=true
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()
HumanCodingGenes = HGNC["symbol"].values
HGNC_symbols = HGNC["symbol"].values

# %% hidden=true
SaveDIR = "/home/jw3514/Work/data/Human_Brain_Cell_Atlas/Subcluster_GeneXCell/"

# %% hidden=true
Subcluster_anno = pd.read_excel("/home/jw3514/Work/data/HumanBrainCellType/subcluster_annotation.xlsx",
                               index_col="Subcluster")

# %% hidden=true
Subcluster_anno.shape

# %% hidden=true
Subcluster_anno.head(2)

# %% hidden=true
h5ad = sc.read_h5ad("/home/jw3514/Work/data/HumanBrainCellType/SuperTypeRawDat/Supercluster_CGE_interneuron.h5ad")

# %% hidden=true
#def dumpsubcluster()
Supercluster = "CGE interneuron"

Supercluster_fname = re.sub(r'\W+', '_', Supercluster)
Supercluster_fname = "Supercluster_{}.h5ad".format(Supercluster_fname)
#h5ad = sc.read_h5ad("/home/jw3514/Work/data/HumanBrainCellType/SuperTypeRawDat/{}".format(Supercluster_fname))

# %% hidden=true
_subcluster_id = 585
_subcluster_ann = Subcluster_anno.loc[_subcluster_id, :]
_subcluster_name = "{}-{}-{}".format(_subcluster_id, _subcluster_ann["Cluster"], _subcluster_ann["Supercluster"])
print(_subcluster_name)

# %% hidden=true
subcluster_h5ad = h5ad[h5ad.obs["subcluster_id"]==_subcluster_id]

# %% hidden=true
subcluster_h5ad

# %% hidden=true
subcluster_h5ad.var = subcluster_h5ad.var.set_index("Gene")

# %% hidden=true
DF = subcluster_h5ad.to_df()

# %% hidden=true
# Only keep protein coding genes
genes = np.array(list(set(DF.columns.values).intersection(set(HGNC_symbols))))
DF = DF[genes]
DF = DF.transpose()
DF.index = [int(GeneSymbol2Entrez[x]) for x in DF.index.values]
DF = DF.sort_index()

# %% hidden=true
#SaveDIR = "/home/jw3514/Work/data/Human_Brain_Cell_Atlas/Subcluster_GeneXCell/"
SaveDIR = "/home/jw3514/Work/data/HumanBrainCellType/Subcluster_GeneXCell/"
DF.to_csv("{}/{}.csv.gz".format(SaveDIR, _subcluster_id), compression="gzip")

# %% hidden=true
## Aggregate splited files into one matrix

Indv_cluster_means = []
Indv_cluster_total_UMIs = []
Gene_Total_Exp = []
Total_N_Cells = 0
for _subcluster_id, _subcluster_ann in Subcluster_anno.iterrows():
    _subcluster_name = "{}-{}-{}".format(_subcluster_id, _subcluster_ann["Cluster"], _subcluster_ann["Supercluster"])
    gene_mean_logUMI = pd.read_csv(("/home/jw3514/Work/data/HumanBrainCellType/Subcluster_GeneXCell/SplitCTs/{}.csv".format(_subcluster_id)),
                                  index_col=0)
    gene_mean_logUMI.rename(columns={"0": _subcluster_name}, inplace=True)
    #if len(gene_mean_logUMI) != 17938:
    #    print(cluster, cluster_clean_name, len(gene_mean_logUMI))
    #Total_N_Cells += cluster_v3_cell_counts
    #print(cluster_v3_cell_counts)
    #if len(Gene_Total_Exp) == 0:
    #    Gene_Total_Exp = gene_mean_logUMI * cluster_v3_cell_counts
    #else:
    #    Gene_Total_Exp += gene_mean_logUMI * cluster_v3_cell_counts
    #gene_mean_logUMI.name = cluster
    Indv_cluster_means.append(gene_mean_logUMI)

# Make and save cluster Exp Mat
Cluster_Exp_DF = pd.concat(Indv_cluster_means, axis=1)
Cluster_Exp_DF.to_csv("/home/jw3514/Work/data/HumanBrainCellType/subcluster_MeanLogUMI.csv")

# %% hidden=true
Cluster_Exp_DF

# %% [markdown] heading_collapsed=true
# ## Make Z2 Mat Cluster, from published cluster level aggregation data

# %% hidden=true
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()
HumanCodingGenes = HGNC["symbol"].values
HGNC_symbols = HGNC["symbol"].values

# %% hidden=true
fil = "/home/jw3514/Work/data/HumanBrainCellType/adult_human_20221007.agg.loom"
adata = sc.read_loom(fil)

# %% hidden=true
adata.var_names_make_unique()

# %% hidden=true
CT_ExpMat = pd.DataFrame(data=adata.X.toarray(), index=adata.obs.index.values, columns=adata.var.index.values)

# %% hidden=true
genes = np.array(list(set(CT_ExpMat.columns.values).intersection(set(HGNC_symbols))))
CT_ExpMat = CT_ExpMat[genes]
CT_ExpMat = CT_ExpMat.transpose()
CT_ExpMat.columns = [int(x) for x in CT_ExpMat.columns.values]

# %% hidden=true
Anno = pd.read_excel(str(PROJ_DIR.parent.parent / "data/HumanBrainCellType/annotation.xlsx"), index_col="Cluster")
Anno.drop(Anno.tail(1).index,inplace=True) # drop last n rows
Anno.index = [int(x) for x in Anno.index.values]
Neur_idx = [str(int(x)) for x in Anno[Anno["Class auto-annotation"]=="NEUR"].index]
Neur_idx = list(set(CT_ExpMat.columns.values).intersection(set(Neur_idx)))

# %% hidden=true
CT_ExpMat.index = [int(GeneSymbol2Entrez[x]) for x in CT_ExpMat.index.values]
CT_ExpMat.to_csv(str(PROJ_DIR / "dat/Human.CT.Exp.Entrez.csv"))

# %% [markdown] hidden=true
# ##### [Start] calculate Overall expression level for each gene

# %% hidden=true
#CT_ExpMat_log2 = np.log2(CT_ExpMat+1)
#CT_ExpMat_log2.to_csv(str(PROJ_DIR / "dat/Human.CT.Exp.Entrez.log2.csv"))

# %% hidden=true
CT_ExpMat_log2 = pd.read_csv(str(PROJ_DIR / "dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.csv"), index_col=0)

# %% hidden=true
plt.hist(CT_ExpMat_log2.values.flatten())


# %% hidden=true

# %% hidden=true
def Z1Conversion(ExpMat, outname="test.z1.mat"):
    Z_mat = []
    for g, row in ExpMat.iterrows():
        tmp = ZscoreConverting(row.values)
        Z_mat.append(tmp)
    Z_mat = np.array(Z_mat)
    CT_Z1_DF = pd.DataFrame(data=Z_mat, index=ExpMat.index.values, 
                            columns=ExpMat.columns.values)
    CT_Z1_DF.to_csv(outname)
    return CT_Z1_DF


# %% hidden=true
ExpDF_Z1 = Z1Conversion(CT_ExpMat_log2, str(PROJ_DIR / "dat/Human.CT.Exp.Entrez.log2.Z1.csv"))
max_Z, min_Z = 5, -5
ExpDF_Z1_clipped = ExpDF_Z1.clip(upper=max_Z, lower=min_Z)
ExpDF_Z1_clipped.to_csv(str(PROJ_DIR / "dat/Human.CT.Exp.Entrez.log2.Z1.clip.csv"))

# %% hidden=true
max_Z, min_Z = 3, -3
ExpDF_Z1_clipped = ExpDF_Z1.clip(upper=max_Z, lower=min_Z)
ExpDF_Z1_clipped.to_csv(str(PROJ_DIR / "dat/Human.CT.Exp.Entrez.log2.Z1.clip3.csv"))

# %% hidden=true
# Combine Z2 mat

# %% hidden=true
Z2_split_dir = "/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Z2.ExpL.HumanCTMatch/"
DFs = []
for file in os.listdir(Z2_split_dir):
    df = pd.read_csv(Z2_split_dir + file, index_col = 0)
    DFs.append(df)
Z2_MAt = pd.concat(DFs)
Z2_MAt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.csv")

# %% hidden=true
Z2_split_dir = "/home/jw3514/Work/CellType_Psy/dat/Z2.Split/Z2.ExpL.BrainSpanMatch/"
DFs = []
for file in os.listdir(Z2_split_dir):
    df = pd.read_csv(Z2_split_dir + file, index_col = 0)
    DFs.append(df)
Z2_MAt = pd.concat(DFs)
Z2_MAt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.BSP.csv")

# %% hidden=true
Z2_split_dir = "/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.ExpL.HumanCTMatch_z1clip3/"
DFs = []
for file in os.listdir(Z2_split_dir):
    df = pd.read_csv(Z2_split_dir + file, index_col = 0)
    DFs.append(df)
Z2_MAt = pd.concat(DFs)
Z2_MAt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv")

# %% [markdown] hidden=true
# #### Compare 2 Z2

# %% hidden=true
Z2_Z1clip5 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.csv", index_col=0)
Z2_Z1clip3 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv", index_col=0)

# %% hidden=true
genes = Z2_Z1clip3.index.values

# %% hidden=true
X = Z2_Z1clip5.loc[genes, "1"]
Y = Z2_Z1clip3.loc[genes, "1"]

# %% hidden=true
plt.scatter(X,Y,s=0.5)
plt.plot([-5,8], [-5, 8], color="grey")
plt.grid(True)

# %% hidden=true
sns.kdeplot(X)
sns.kdeplot(Y)

# %% [markdown] heading_collapsed=true
# ## Make Z2 Subcluster

# %% hidden=true
Cluster_Exp_DF = pd.read_csv(str(PROJ_DIR / "dat/HumanCTExpressionMats/subcluster_MeanLogUMI.csv"), index_col=0)

# %% hidden=true
plt.hist(Cluster_Exp_DF.values.flatten())

# %% hidden=true
SubCluster_Z1 = Z1Conversion(Cluster_Exp_DF, str(PROJ_DIR / "dat/Human.Subcluster.log2.Z1.csv"))
max_Z, min_Z = 3, -3
SubCluster_Z1_clipped = SubCluster_Z1.clip(upper=max_Z, lower=min_Z)
SubCluster_Z1_clipped.to_csv(str(PROJ_DIR / "dat/HumanCTExpressionMats/Human.Subcluster.log2.Z1.clip3.csv"))

# %% hidden=true
SubCluster_Z1_clipped.head(2)

# %% hidden=true
for X in SubCluster_Z1_clipped.columns.values:
    print(X)

# %% hidden=true
SubCluster_Z1_clipped = SubCluster_Z1_clipped[~SubCluster_Z1_clipped.index.duplicated(keep='first')]

# %% hidden=true
SubCluster_Z1_clipped.to_csv(str(PROJ_DIR / "dat/HumanCTExpressionMats/Human.Subcluster.log2.Z1.clip3.csv"))

# %% hidden=true
SubCluster_Z1_clipped.loc[9782, "0-297-Upper rhombic lip"]

# %% hidden=true

# %% hidden=true
Z2_split_dir = "/home/jw3514/Work/CellType_Psy/dat2/Z2.Split/Z2.Subcluster.HumanCTMatch.clip3/"
DFs = []
for file in os.listdir(Z2_split_dir):
    df = pd.read_csv(Z2_split_dir + file, index_col = 0)
    DFs.append(df)
Z2_MAt = pd.concat(DFs)
Z2_MAt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Subcluster.log2.Z2.clip3.csv")

# %% hidden=true
Z2_MAt.head(2)

# %% hidden=true
max_Z, min_Z = 3, -3
Z2_MAt = Z2_MAt.clip(upper=max_Z, lower=min_Z)
Z2_MAt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Subcluster.log2.Z2.clip3.3.csv")

# %% [markdown]
# # Re processing Cluster level data

# %%
Subcluster_anno = pd.read_excel("/home/jw3514/Work/data/HumanBrainCellType/subcluster_annotation.xlsx", index_col="Subcluster")
DatDIR = "/home/jw3514/Work/data/HumanBrainCellType/Subcluster_GeneXCell/"


# %%
def processHumanCT_Cluster(cluster, annotation, DatDIR):
    print("Processing {}".format(cluster))
    subclusters = annotation[annotation["Cluster"]==cluster].index.values
    print(subclusters)
    
    dfs = []
    for subcluster in subclusters:
        df = pd.read_csv("{}/{}.csv.gz".format(DatDIR, subcluster), index_col=0)
        dfs.append(df)
    
    # Concatenate all dataframes along columns
    combined_df = pd.concat(dfs, axis=1)
    combined_df.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/{}.csv.gz".format(cluster))
    
    # Calculate plain mean for comparison
    gene_means = combined_df.values.mean(axis=1) 
    gene_means = pd.Series(data=gene_means, index=combined_df.index)
    gene_means.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_mean.csv".format(cluster))
    
    # Apply log2(x+1) transformation and take mean across cells
    gene_log2_means = np.log2(combined_df.values + 1).mean(axis=1)
    gene_log2_means = pd.Series(data=gene_log2_means, index=combined_df.index)
    gene_log2_means.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_log2mean.csv".format(cluster))

        #Gene_Cluster_Mean = pd.Series(data=gene_dat, index=gene_index)
        #Gene_Cluster_Mean.to_csv("/home/jw3514/Work/data/HumanBrainCellType/Subcluster_GeneXCell/SplitCTs/{}.csv".format(cluster))


# %%
cluster = 308
processHumanCT_Cluster(cluster, Subcluster_anno, DatDIR)

# %%
subcluster = 1
df = pd.read_csv("{}/{}.csv.gz".format(DatDIR, subcluster), index_col=0)
df.shape



# %%
Clusters = Subcluster_anno["Cluster"].unique()
Clusters.sort()

# %%
Anno.head(2)

# %%
## Aggregate splited files into one matrix

Indv_cluster_means = []
Indv_cluster_log2means = []
Missing_Clusters = []
for cluster, row in Anno.iterrows():
    #_subcluster_name = "{}-{}-{}".format(_subcluster_id, _subcluster_ann["Cluster"], _subcluster_ann["Supercluster"])
    try:
        gene_mean_UMI = pd.read_csv(("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_mean.csv".format(cluster)), index_col=0)
        gene_mean_logUMI = pd.read_csv(("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_log2mean.csv".format(cluster)), index_col=0)
        gene_mean_UMI.rename(columns={"0": cluster}, inplace=True)
        gene_mean_logUMI.rename(columns={"0": cluster}, inplace=True)

        Indv_cluster_means.append(gene_mean_UMI)
        Indv_cluster_log2means.append(gene_mean_logUMI)
    except:
        print("{} not found".format(cluster))
        Missing_Clusters.append(cluster)
# Make and save cluster Exp Mat
Cluster_Exp_DF = pd.concat(Indv_cluster_means, axis=1)
Cluster_Exp_DF.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_MeanUMI.csv")
Cluster_Exp_DF_log2 = pd.concat(Indv_cluster_log2means, axis=1)
Cluster_Exp_DF_log2.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_MeanLogUMI.csv")


# %%
import gzip
def processHumanCT_Cluster_LargeFile(cluster, annotation, DatDIR):
    print("Processing {}".format(cluster))
    subclusters = annotation[annotation["Cluster"]==cluster].index.values
    print(subclusters)
    
    # Initialize dictionaries to store sums and counts
    gene_sums = {}
    gene_counts = {}
    gene_log2_sums = {}
    first_file = True
    
    # Process each subcluster file
    for subcluster in subclusters:
        with gzip.open("{}/{}.csv.gz".format(DatDIR, subcluster), 'rt') as f:
            reader = csv.reader(f)
            header = next(reader) # Skip header
            
            # For first file, initialize gene names
            if first_file:
                gene_names = []
                for row in reader:
                    gene_name = row[0]
                    gene_names.append(gene_name)
                    values = [float(x) for x in row[1:]]
                    gene_sums[gene_name] = sum(values)
                    gene_counts[gene_name] = len(values)
                    gene_log2_sums[gene_name] = sum(np.log2(np.array(values) + 1))
                first_file = False
            else:
                for row in reader:
                    gene_name = row[0]
                    values = [float(x) for x in row[1:]]
                    gene_sums[gene_name] += sum(values)
                    gene_counts[gene_name] += len(values)
                    gene_log2_sums[gene_name] += sum(np.log2(np.array(values) + 1))
    
    # Calculate means and save results
    with open("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_mean.csv".format(cluster), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['gene','0'])
        for gene in gene_names:
            mean = gene_sums[gene] / gene_counts[gene]
            writer.writerow([gene, mean])
            
    with open("/home/jw3514/Work/data/HumanBrainCellType/cluster_GeneXCell/SplitCTs/{}_log2mean.csv".format(cluster), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['gene','0']) 
        for gene in gene_names:
            log2_mean = gene_log2_sums[gene] / gene_counts[gene]
            writer.writerow([gene, log2_mean])

# %%
#Cluster_Exp_DF.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_MeanUMI.csv")
#Cluster_Exp_DF_log2.to_csv("/home/jw3514/Work/data/HumanBrainCellType/cluster_MeanLogUMI.csv")


# %%
Cluster_Exp_DF.head(2)

# %%
Cluster_Exp_DF_log2.head(2)


# %%
genes = Cluster_Exp_DF.index.values
ct = 10
X = Cluster_Exp_DF_log2.loc[genes, ct]
Y = Cluster_Exp_DF.loc[genes, ct]
plt.scatter(X,Y,s=0.5)
#plt.plot([-5,8], [-5, 8], color="grey")
plt.grid(True)
plt.xlabel("Log2(UMI+1) mean")
plt.ylabel("UMI mean")
plt.show()


# %%
HumanCTExp_Entrez = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv", index_col=0)
HumanCTExp_Entrez.head(2)
# make column int
HumanCTExp_Entrez.columns = [int(x) for x in HumanCTExp_Entrez.columns.values]
HumanCTExp_Entrez = HumanCTExp_Entrez.loc[~HumanCTExp_Entrez.index.duplicated(keep='first')]

# %%

Cluster_Exp_DF = Cluster_Exp_DF.loc[~Cluster_Exp_DF.index.duplicated(keep='first')]


# %%
print(len(X), len(Y))

# %%
gene1 = Cluster_Exp_DF.index.values
gene2 = HumanCTExp_Entrez.index.values
common_genes = list(set(gene1).intersection(set(gene2)))
print(len(common_genes))
ct = 10
X = HumanCTExp_Entrez.loc[common_genes, ct].values  # Add .values to ensure numpy array
Y = Cluster_Exp_DF.loc[common_genes, ct].values     # Add .values to ensure numpy array
# compute correlation
print(np.corrcoef(X, Y))
print(len(X), len(Y))
plt.figure()  # Create new figure
plt.scatter(X, Y, s=0.5)
plt.grid(True)
plt.xlabel("Log2(UMI+1) mean")
plt.ylabel("UMI mean")
plt.show()


# %%
Cluster_Exp_DF_log2 = Cluster_Exp_DF_log2.loc[~Cluster_Exp_DF_log2.index.duplicated(keep='first')]

# %%
Cluster_Exp_DF_log2.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.csv")

# %%
plt.hist(Cluster_Exp_DF_log2[300], bins=100)

# %%
plt.hist(Cluster_Exp_DF[400], bins=100)


# %%
def ZscoreConverting_V2(values, mean=np.nan, std=np.nan, low_exp = 0, min_z=-5): 
    """
    Convert values to z-scores with special handling for zeros:
    - Build distribution using only non-zero values
    - Set minimum z-score to min_z (default -5)
    - Set all zero expressions to min_z (default -5)
    
    Args:
        values: Array-like input values
        mean: Optional pre-computed mean
        std: Optional pre-computed standard deviation
        min_z: Minimum z-score (default -5)
    Returns:
        numpy array of z-scores
    """
    # Convert to numpy array and identify non-zero values
    values = np.array(values)
    non_zero_mask = values >= low_exp
    non_zero_values = values[non_zero_mask]
    
    # If no non-zero values, return array of -5
    if len(non_zero_values) == 0:
        return np.full_like(values, min_z)
    
    # Calculate mean and std from non-zero values if not provided
    if mean != mean:  # Check for nan
        mean = np.mean(non_zero_values)
    if std != std:    # Check for nan
        std = np.std(non_zero_values)
        # Handle case where std is 0
        if std == 0:
            std = 1
    
    # Calculate z-scores
    zscores = np.full_like(values, min_z)  # Initialize with -5
    non_zero_zscores = (non_zero_values - mean) / std
    
    # Apply minimum threshold 
    non_zero_zscores = np.maximum(non_zero_zscores, min_z)
    
    # Put non-zero z-scores back in original positions
    zscores[non_zero_mask] = non_zero_zscores
    
    return zscores

def Z1Conversion(ExpMat, outname="test.z1.mat"):
    Z_mat = []
    for g, row in ExpMat.iterrows():
        tmp = ZscoreConverting(row.values)
        Z_mat.append(tmp)
    Z_mat = np.array(Z_mat)
    CT_Z1_DF = pd.DataFrame(data=Z_mat, index=ExpMat.index.values, 
                            columns=ExpMat.columns.values)
    CT_Z1_DF.to_csv(outname)
    return CT_Z1_DF

def Z1Conversion_V2(ExpMat, outname="test.z1.mat", low_exp = 0, min_z=-5):
    """
    Convert expression matrix to z-scores with zero handling
    
    Args:
        ExpMat: pandas DataFrame with genes as rows and cell types as columns
        outname: output file name for saving results
        
    Returns:
        pandas DataFrame with z-scores
    """
    Z_mat = []
    for g, row in ExpMat.iterrows():
        tmp = ZscoreConverting_V2(row.values, min_z=min_z)
        Z_mat.append(tmp)
    
    Z_mat = np.array(Z_mat)
    CT_Z1_DF = pd.DataFrame(data=Z_mat, 
                           index=ExpMat.index.values,
                           columns=ExpMat.columns.values)
    
    CT_Z1_DF.to_csv(outname)
    return CT_Z1_DF


# %%
Cluster_Exp_DF_log2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.csv", index_col=0)

# %%
# Check how 0-expressed genes specificity score 
HumanCT_Z1 = Z1Conversion_V2(Cluster_Exp_DF_log2, min_z=-5)
max_Z, min_Z = 5, -5
HumanCT_Z1 = HumanCT_Z1.clip(upper=max_Z, lower=min_Z)
HumanCT_Z1.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip5.csv")

# %%
# Check how 0-expressed genes specificity score 
HumanCT_Z1 = Z1Conversion_V2(Cluster_Exp_DF_log2, min_z=-3)
max_Z, min_Z = 3, -3
HumanCT_Z1 = HumanCT_Z1.clip(upper=max_Z, lower=min_Z)
HumanCT_Z1.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip3.csv")

# %%
HumanCT_Z1_2 = Z1Conversion_V2(Cluster_Exp_DF_log2, low_exp=1e-3, min_z=-5)
max_Z, min_Z = 5, -5
HumanCT_Z1_2 = HumanCT_Z1_2.clip(upper=max_Z, lower=min_Z)
HumanCT_Z1_2.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.Z1.clip5.lowexp_1e-3.csv")

# %%

# %%

# %%
HumanCT_Z1.columns = [int(x) for x in HumanCT_Z1.columns.values]
Cluster_Exp_DF_log2.columns = [int(x) for x in Cluster_Exp_DF_log2.columns.values]

# %%
HumanCT_Z1_old = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z1.csv", index_col=0)
HumanCT_Z1_old.columns = [int(x) for x in HumanCT_Z1_old.columns.values]

# %%
# Get gene lists from both matrices
gene1 = HumanCT_Z1.index.values
gene2 = HumanCT_Z1_old.index.values
common_genes = list(set(gene1).intersection(set(gene2)))
print(len(common_genes))

# Convert column index to string since that appears to be the issue
ct = 400  # Convert to string to match column names

# Extract values for common genes
X = HumanCT_Z1.loc[common_genes, ct].values
Y = HumanCT_Z1_old.loc[common_genes, ct].values
Z = Cluster_Exp_DF_log2.loc[common_genes, ct].values


# Create scatter plot
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

ax1.scatter(X, Y, s=0.5)
ax1.grid(True)
ax1.set_xlabel("New Z1")
ax1.set_ylabel("Old Z1")

ax2.scatter(X, Z, s=0.5) 
ax2.grid(True)
ax2.set_xlabel("New Z1")
ax2.set_ylabel("Log2 Expression")

ax3.scatter(Y, Z, s=0.5)
ax3.grid(True)
ax3.set_xlabel("Old Z1")
ax3.set_ylabel("Log2 Expression")

plt.tight_layout()
plt.show()


# %%
#HumanCT_Z2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv", index_col=0)
HumanCT_Z2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z2.csv", index_col=0)
HumanCT_Z2_Res = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Residue.Z2.csv", index_col=0)
HumanCT_Z2.columns = [int(x) for x in HumanCT_Z2.columns.values]
HumanCT_Z2_Res.columns = [int(x) for x in HumanCT_Z2_Res.columns.values]
# set max min
max_Z, min_Z = 3, -3
HumanCT_Z2_Res = HumanCT_Z2_Res.clip(upper=max_Z, lower=min_Z)

# %%
#HumanCT_Z2_Res.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Residue.Z2.clip3.3.csv")

# %%
# Create figure and axes
CT_index = 400 

fig = plt.figure(figsize=(10, 10))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

# Plot histograms
zero_expr_genes = Cluster_Exp_DF_log2[Cluster_Exp_DF_log2[CT_index] == 0].index
ax1.hist(Cluster_Exp_DF_log2[CT_index], bins=100)
ax1.set_title('Cluster_Exp_DF_log2 (Zero Expression Genes)')

ax2.hist(HumanCT_Z1[CT_index], bins=100, density=True)
#ax2.hist(HumanCT_Z1.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax2.set_title('HumanCT_Z1')

ax3.hist(HumanCT_Z2[CT_index], bins=100, density=True) 
#ax3.hist(HumanCT_Z2.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax3.set_title('HumanCT_Z2')

ax4.hist(HumanCT_Z2_Res[CT_index], bins=100, density=True)
ax4.hist(HumanCT_Z2_Res.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax4.set_title('HumanCT_Z2_Res')

plt.tight_layout()
plt.show()

# %%

# %%
HumanCT_Z2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.csv", index_col=0)
HumanCT_ExpL = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.csv", index_col=0)

# %%
#Expression Level vs Specificity 
ExpLevels = []
TopNBias = []
CommonGenes = list(set(HumanCT_Z2.index.values).intersection(set(HumanCT_ExpL.index.values)))
for g in CommonGenes:
    row_sorted = HumanCT_Z2.loc[g].sort_values(ascending=False)
    top_N = row_sorted.head(10)
    top_N_spec = np.mean(top_N)
    ExpLevels.append(HumanCT_ExpL.loc[g, :].mean())
    TopNBias.append(top_N_spec)
ExpLevels = np.array(ExpLevels)
ExpLevels = np.log10(np.exp2(ExpLevels) )
TopNBias = np.array(TopNBias)


# %%
plt.figure(dpi=120)
plt.scatter(ExpLevels, TopNBias, s=0.1)
plt.xscale('log')
plt.xlabel("Expression Level")
plt.ylabel("Top {} Specificity Score".format(10))
plt.show()


# %%

# %%

# %%

# %% [markdown]
# ### Test remove zero expressed genes from Z1 matrix 
#

# %%
def compute_specificity_matrix(expr_matrix):
    """
    Compute specificity scores for zero-inflated, small-valued expression data
    
    Parameters:
    expr_matrix: genes x cell_types matrix with values >= 0
    """
    n_genes, n_celltypes = expr_matrix.shape
    specificity_matrix = pd.DataFrame(0, index=expr_matrix.index, columns=expr_matrix.columns)

    for gene_idx, expr in expr_matrix.iterrows():
        non_zero_mask = expr > 0
        non_zero_expr = expr[non_zero_mask]
        
        if len(non_zero_expr) == 0:
            specificity_matrix.loc[gene_idx,:] = 0
            continue
            
        expr_mean = np.mean(non_zero_expr)
        expr_std = np.std(non_zero_expr) if len(non_zero_expr) > 1 else expr_mean
        
        # Calculate specificity only for expressed values
        specificity = np.zeros_like(expr)
        specificity[non_zero_mask] = (expr[non_zero_mask] - expr_mean) / (expr_std + 1e-10)
        
        # Zero positions get minimum specificity
        min_spec = np.min(specificity[non_zero_mask])
        specificity[~non_zero_mask] = min_spec
        
        specificity_matrix.loc[gene_idx,:] = specificity
    
    
    return specificity_matrix

# %%
#Cluster_Spec_Test = compute_specificity_matrix(Cluster_Exp_DF_log2)
Cluster_Spec_Test = compute_specificity_matrix(Cluster_Exp_DF_log2)

# %%
Cluster_Spec_Test.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Spec.Test.Nov14.csv")

# %%
Cluster_Spec_Test

# %%
# Plot Expression vs Specificity of cell types
CT_index = 200

fig = plt.figure(figsize=(10, 10))
plt.scatter(Cluster_Exp_DF_log2[CT_index], Cluster_Spec_Test[CT_index], s=0.5)
plt.xlabel("Log2(UMI+1) mean")
plt.ylabel("Specificity")
plt.show()


# %%
# Create figure and axes
CT_index = 1 

fig = plt.figure(figsize=(10, 10))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

# Plot histograms
zero_expr_genes = Cluster_Exp_DF_log2[Cluster_Exp_DF_log2[CT_index] == 0].index
ax1.hist(Cluster_Exp_DF_log2[CT_index], bins=100)
ax1.set_title('Cluster_Exp_DF_log2 (Zero Expression Genes)')

ax2.hist(HumanCT_Z1[CT_index], bins=100, density=True)
#ax2.hist(HumanCT_Z1.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax2.set_title('HumanCT_Z1')

ax3.hist(Cluster_Spec_Test[CT_index], bins=100, density=True) 
ax3.hist(Cluster_Spec_Test.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax3.set_title('HumanCT_Z2')

ax4.hist(HumanCT_Z2_Res[CT_index], bins=100, density=True)
ax4.hist(HumanCT_Z2_Res.loc[zero_expr_genes, CT_index], bins=100, alpha=0.5, color='red', density=True)
ax4.set_title('HumanCT_Z2_Res')

plt.tight_layout()
plt.show()

# %%

# %%

# %%

# %%
expr_matrix = Cluster_Exp_DF_log2
specificity_matrix = pd.DataFrame(0, index=expr_matrix.index, columns=expr_matrix.columns)

for gene_idx, expr in expr_matrix.iterrows():
    non_zero_mask = expr > 0
    non_zero_expr = expr[non_zero_mask]
    
    if len(non_zero_expr) == 0:
        continue
        
    # Calculate specificity incorporating both components
    expr_mean = np.mean(non_zero_expr)
    expr_std = np.std(non_zero_expr) if len(non_zero_expr) > 1 else expr_mean
    
    # Modified score combining expression level and zero-fraction
    zero_fraction = 1 - (np.sum(non_zero_mask) / len(expr))
    scaled_expr = (expr - expr_mean) / (expr_std + 1e-10)
    specificity_matrix.loc[gene_idx,:] = scaled_expr * (1 - zero_fraction)

# %%
specificity_matrix

# %%
Cluster_Exp_DF_log2.shape

# %%
log2_umi_matrix = Cluster_Exp_DF_log2
count_matrix = np.exp2(log2_umi_matrix) - 1

n_genes, n_celltypes = count_matrix.shape
specificity_matrix = np.zeros_like(log2_umi_matrix)

for gene_idx, row in count_matrix.iterrows():
    counts = row
    
    # Fit ZINB model
    exog = np.ones_like(counts)
    try:
        model = sm.ZeroInflatedNegativeBinomialP(counts, exog)
        results = model.fit(disp=0)
        
        # Get probabilities of observing counts
        log_probs = results.predict(exog, which='prob')
        
        # Convert to specificity score
        specificity = -np.log10(log_probs + 1e-10)
        specificity_matrix[gene_idx,:] = specificity
        
    except:
        # Fallback if ZINB fit fails
        specificity_matrix[gene_idx,:] = stats.zscore(counts)

    break 

# %%
HumanCTExp_Entrez

# %%
count_matrix

# %%
#0 expressed genes 

# %%
#Expression level vs specificity 


# %%
HumanCT_OverallEXP = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/HumanCT.MatchDF.csv", index_col=0)

# %%
np.log10(HumanCT_OverallEXP["Exp"] + 1).hist(bins=100)

# %%
HumanCT_Z2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv", index_col=0)

# %%
# for each gene, calculate its specificity score of top N cell types 
topNCT = 1
Gene_Spec_Score = []
for g, row in HumanCT_Z2.iterrows():
    row_sorted = row.sort_values(ascending=False)
    top_N = row_sorted.head(topNCT)
    top_N_spec = np.mean(top_N)
    Gene_Spec_Score.append(top_N_spec)
Gene_Spec_Score = np.array(Gene_Spec_Score)



# %%
plt.scatter(np.log10(HumanCT_OverallEXP["Exp"] + 1), Gene_Spec_Score, s=0.1)
plt.xlabel("Log10(Expression + 1)")
plt.ylabel("top {} Specificity Score".format(topNCT))
plt.show()

# %%
topNCT = 1
Gene_Spec_Score = []
Gene_Exp = []
for g, row in HumanCT_Z2.iterrows():
    row_sorted = row.sort_values(ascending=False)
    top_N = row_sorted.head(topNCT)
    top_N_spec = np.mean(top_N)
    Gene_Spec_Score.append(top_N_spec)
    Gene_Exp.append(HumanCT_OverallEXP.loc[g, "Exp"])
    
Gene_Spec_Score = np.array(Gene_Spec_Score)
Gene_Exp = np.array(Gene_Exp)


# %%
plt.scatter(np.log10(Gene_Exp + 1), Gene_Spec_Score, s=0.1)
plt.xlabel("Log10(Expression + 1)")
plt.ylabel("top {} Specificity Score".format(topNCT))
plt.show()

# %%
HumanCT_Z2_test = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z2.csv", index_col=0)

# %%
topNCT = 10
Gene_Spec_Score = []
for g, row in HumanCT_Z2.iterrows():
    row_sorted = row.sort_values(ascending=False)
    top_N = row_sorted.head(topNCT)
    top_N_spec = np.mean(top_N)
    Gene_Spec_Score.append(top_N_spec)
Gene_Spec_Score = np.array(Gene_Spec_Score)

# %%
plt.scatter(np.log10(HumanCT_OverallEXP["Exp"] + 1), Gene_Spec_Score, s=0.1)
plt.xlabel("Log10(Expression + 1)")
plt.ylabel("top {} Specificity Score".format(topNCT))
plt.show()

# %%

# %%
# Try specificity score of each cell type 
Cluster_Exp_DF_log2 = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Exp.csv", index_col=0)

# %%
Cluster_Exp_DF_log2.head(10)

# %%
# Get column sums and mean of sums
col_sums = Cluster_Exp_DF_log2.sum(axis=0)
sum_mean = col_sums.mean()

# Normalize each column to have the same sum (1000 * mean)
target_sum = 1000 * sum_mean
normalized_expr = Cluster_Exp_DF_log2.div(col_sums, axis=1) * target_sum

# %%
normalized_expr.head(10)

# %%
# Calculate total expression sum across cell types for each gene
total_exp_per_gene = normalized_expr.sum(axis=1)

# Divide each expression value by the total for that gene to get fraction
specificity_scores = normalized_expr.div(total_exp_per_gene, axis=0)



# %%
specificity_scores.head(10)


# %%
specificity_scores.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Spec.Score.csv")

# %% [markdown]
# # Try Not log2 Z2 

# %%
Cluster_Exp_DF_log2.head(2)

# %%
Cluster_Exp_DF_original = np.expm1(Cluster_Exp_DF_log2)

# %%
Cluster_Exp_DF_original

# %%
expression_df = Cluster_Exp_DF_original

# %%
cell_type_sums = expression_df.sum(axis=0)
tpm_df = expression_df.div(cell_type_sums, axis=1) * 1000000

min_tpm = 1
# Create a mask for filtering (only used for specificity calculation)
if min_tpm > 0:
    mask = tpm_df >= min_tpm
else:
    mask = pd.DataFrame(True, index=tpm_df.index, columns=tpm_df.columns)

# For specificity calculation, treat values below threshold as 0
filtered_tpm = tpm_df.copy()
filtered_tpm[~mask] = 0
expression_filtered = expression_df.copy()
expression_filtered[~mask] = 0

# Calculate specificity using the filtered TPM values
gene_total_tpm = filtered_tpm.sum(axis=1)
specificity_df = filtered_tpm.div(gene_total_tpm, axis=0)

# Handle genes that have all zeros after filtering
specificity_df.fillna(0, inplace=True)


# %%
expression_filtered

# %%
HumanCT_Z1 = Z1Conversion_V2(expression_df, min_z=-5)
max_Z, min_Z = 5, -5
HumanCT_Z1 = HumanCT_Z1.clip(upper=max_Z, lower=min_Z)
HumanCT_Z1.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.UMI.Exp.Z1.clip5.Apr18csv")

# %%
HumanCT_Z1_withFilt = Z1Conversion_V2(expression_filtered, min_z=-5)
max_Z, min_Z = 5, -5
HumanCT_Z1_withFilt = HumanCT_Z1_withFilt.clip(upper=max_Z, lower=min_Z)
HumanCT_Z1_withFilt.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.UMI.Exp.Z1.clip5.withFilt.Apr18.csv")

# %%
expression_filtered.to_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.UMI.Exp.withFilt.Apr18.csv")

# %%
