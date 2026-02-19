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
from CellType_PSY import *
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
#Functions:
def plot_vip_effect_comparison(df, effect_col = "EFFECT"):
    """
    Plot comparison of effects between VIP+ and VIP- cells
    
    Parameters:
    df: DataFrame with columns 'VIP' and 'EFFECT'
    """
    vip_pos = df[df["VIP"] >= 1]
    vip_neg = df[df["VIP"] < 1]
    
    # plot effect of VIP+ vs VIP-
    data = [vip_pos[effect_col], vip_neg[effect_col]]
    
    # Perform Mann-Whitney U test
    stat, pval = scipy.stats.mannwhitneyu(vip_pos[effect_col], 
                                        vip_neg[effect_col])
    
    plt.boxplot(data, labels=["VIP+", "VIP-"], showfliers=True)
    plt.plot([1] * len(vip_pos[effect_col]), 
             vip_pos[effect_col], 'ko', alpha=0.3)
    plt.plot([2] * len(vip_neg[effect_col]), 
             vip_neg[effect_col], 'ko', alpha=0.3)
    plt.ylabel(effect_col)
    plt.title(f"p = {pval:.2e}")
    plt.show()
VIP_Anno = pd.read_csv("VIP_Anno.csv", index_col=0)

# %% [markdown]
# # Bias Cal

# %%
#HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/Test.BiasMat/HumanCT.Spec.clip.adj.csv", index_col=0)
#HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/Test.BiasMat/HumanCT.Spec.clip.noLowExp.cut1e4.csv", index_col=0)
#HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.1.Filt.Spec.clip.lowexp.cut1e4.csv", index_col=0)
HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv", index_col=0)
#HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv", index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
HumanCT_Z2_HCT.shape

# %% [markdown] heading_collapsed=true
# ## Load GW

# %%
GeneWeightDIR = "../dat/GeneWeights/"
#Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul07/"
Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul24_Centered/"
if not os.path.exists(Bias_Save_Dir): # make dir if not exists
    os.makedirs(Bias_Save_Dir)

# %% [markdown] heading_collapsed=true
# #### ASD Bias

# %% [markdown] hidden=true
# ### HIQ & LIQ ASD

# %%
HIQ_GW = Fil2Dict("{}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))
LIQ_GW = Fil2Dict("{}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw".format(GeneWeightDIR))
print(len(HIQ_GW), len(LIQ_GW))

# %% hidden=true
HIQ_Z2_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, HIQ_GW)
#HIQ_Z2_Bias = AdjustClusterMean(HIQ_Z2_Bias, HumanCT_Z2_HCT_ColMean)
HIQ_Z2_Bias = AnnotateCTDat(HIQ_Z2_Bias, Anno)
HIQ_Z2_Bias.to_csv("{}/HCT.ASD61.HIQ.csv".format(Bias_Save_Dir))

LIQ_Z2_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, LIQ_GW)
#LIQ_Z2_Bias = AdjustClusterMean(LIQ_Z2_Bias, HumanCT_Z2_HCT_ColMean)
LIQ_Z2_Bias = AnnotateCTDat(LIQ_Z2_Bias, Anno)
LIQ_Z2_Bias.to_csv("{}/HCT.ASD61.LIQ.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(HIQ_Z2_Bias, "HIQ ASD HCT")
SuperClusterBias_BoxPlot(LIQ_Z2_Bias, "LIQ ASD HCT")

# %%
# SuperClusterBias_BoxPlot(HIQ_Z2_Bias, "HIQ ASD HCT", EffectCol="EFFECT_Adj")
# SuperClusterBias_BoxPlot(LIQ_Z2_Bias, "LIQ ASD HCT", EffectCol="EFFECT_Adj")

# %% [markdown] heading_collapsed=true
# ## SCZ Bias

# %%
SCZ_GW = Fil2Dict("{}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(GeneWeightDIR))

# %% hidden=true
SCZ_Z2_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW)
#SCZ_Z2_Bias = AdjustClusterMean(SCZ_Z2_Bias, HumanCT_Z2_HCT_ColMean)
SCZ_Z2_Bias = AnnotateCTDat(SCZ_Z2_Bias, Anno)
SuperClusterBias_BoxPlot(SCZ_Z2_Bias, "SCZ")
SCZ_Z2_Bias.to_csv("{}/HCT.SCZ61.csv".format(Bias_Save_Dir))

# %%
SCZ_GW_Ex_NDD61 = Fil2Dict("{}/SCZ.top61.ExlNDD61.gw".format(GeneWeightDIR))
SCZ_GW_Ex_NDD297 = Fil2Dict("{}/SCZ.top61.ExlNDD297.gw".format(GeneWeightDIR))

# Calculate and save bias for SCZ Ex NDD61
SCZ_Z2_Bias_Ex_NDD61 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_Ex_NDD61)
SCZ_Z2_Bias_Ex_NDD61 = AnnotateCTDat(SCZ_Z2_Bias_Ex_NDD61, Anno)
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_Ex_NDD61, "SCZ Ex NDD61")
SCZ_Z2_Bias_Ex_NDD61.to_csv("{}/HCT.SCZ61.ExNDD61.csv".format(Bias_Save_Dir))

# Calculate and save bias for SCZ Ex NDD297
SCZ_Z2_Bias_Ex_NDD297 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_Ex_NDD297)
SCZ_Z2_Bias_Ex_NDD297 = AnnotateCTDat(SCZ_Z2_Bias_Ex_NDD297, Anno)
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_Ex_NDD297, "SCZ Ex NDD297")
SCZ_Z2_Bias_Ex_NDD297.to_csv("{}/HCT.SCZ61.ExNDD297.csv".format(Bias_Save_Dir))

# %%
common_indices = SCZ_Z2_Bias.index.intersection(VIP_Anno.index)
SCZ_Z2_Bias_sub = SCZ_Z2_Bias.loc[common_indices].copy()
SCZ_Z2_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(SCZ_Z2_Bias_sub)

# %%

# %%
SCZ_GW_200 = Fil2Dict("{}/SCZ.top200.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(GeneWeightDIR))
print(len(SCZ_GW_200))
SCZ_Z2_Bias_200 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_200)
SCZ_Z2_Bias_200 = AnnotateCTDat(SCZ_Z2_Bias_200, Anno)
SCZ_Z2_Bias_200.to_csv("{}/HCT.SCZ200.Z2.Spec.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_200, "SCZ HCT",)



# %%
SCZ_GW_500 = Fil2Dict("{}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(GeneWeightDIR))
print(len(SCZ_GW_500))
SCZ_Z2_Bias_500 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_500)
SCZ_Z2_Bias_500 = AnnotateCTDat(SCZ_Z2_Bias_500, Anno)
SCZ_Z2_Bias_500.to_csv("{}/HCT.SCZ500.Z2.Spec.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_500, "SCZ HCT")

# %%
SCZ_GW_Protect = Fil2Dict("{}/SCZ.top61.protect.gw".format(GeneWeightDIR))
print(len(SCZ_GW_Protect))
SCZ_Z2_Bias_Protect = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_Protect)
SCZ_Z2_Bias_Protect = AnnotateCTDat(SCZ_Z2_Bias_Protect, Anno)
SCZ_Z2_Bias_Protect.to_csv("{}/HCT.SCZ500.protect.Z2.Spec.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_Protect, "SCZ HCT")

# Plain weight version (all weights = 1)
SCZ_GW_Protect_plain = {k: 1 for k in SCZ_GW_Protect}
SCZ_Z2_Bias_Protect_plain = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_Protect_plain)
SCZ_Z2_Bias_Protect_plain = AnnotateCTDat(SCZ_Z2_Bias_Protect_plain, Anno)
SCZ_Z2_Bias_Protect_plain.to_csv("{}/HCT.SCZ500.protect.Z2.Spec.plain.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(SCZ_Z2_Bias_Protect_plain, "SCZ HCT (plain weight)")

# %%
SCZ_Z2_Bias_Protect

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
# Keep as Series (don't use .values) so plot_correlation can filter by Neur_idx/NonNeur_idx
# Align indices to ensure both Series have matching indices
values1_series = SCZ_Z2_Bias.sort_index()[EFFECT_Col]
values2_series = SCZ_Z2_Bias_Protect.sort_index()[EFFECT_Col]
common_idx = values1_series.index.intersection(values2_series.index)
values1 = values1_series.loc[common_idx]
values2 = values2_series.loc[common_idx]
plot_correlation(values1, values2, "SCZ", "HIQ ASD HCT", dpi=80)

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
# Keep as Series (don't use .values) so plot_correlation can filter by Neur_idx/NonNeur_idx
# Align indices to ensure both Series have matching indices
values1_series = SCZ_Z2_Bias.sort_index()[EFFECT_Col]
values2_series = SCZ_Z2_Bias_Ex_NDD61.sort_index()[EFFECT_Col]
common_idx = values1_series.index.intersection(values2_series.index)
values1 = values1_series.loc[common_idx]
values2 = values2_series.loc[common_idx] 
plot_correlation(values1, values2, "SCZ", "SCZ Ex NDD61", dpi=80)

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
# Keep as Series (don't use .values) so plot_correlation can filter by Neur_idx/NonNeur_idx
# Align indices to ensure both Series have matching indices
values1_series = SCZ_Z2_Bias.sort_index()[EFFECT_Col]
values2_series = SCZ_Z2_Bias_Ex_NDD297.sort_index()[EFFECT_Col]
common_idx = values1_series.index.intersection(values2_series.index)
values1 = values1_series.loc[common_idx]
values2 = values2_series.loc[common_idx]
plot_correlation(values1, values2, "SCZ", "HIQ ASD HCT", dpi=80)

# %% [markdown]
# ### ASD SCZ correlation after remove common genes

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
# Keep as Series (don't use .values) so plot_correlation can filter by Neur_idx/NonNeur_idx
# Align indices to ensure both Series have matching indices
values1_series = SCZ_Z2_Bias.sort_index()[EFFECT_Col]
values2_series = HIQ_Z2_Bias.sort_index()[EFFECT_Col]
common_idx = values1_series.index.intersection(values2_series.index)
values1 = values1_series.loc[common_idx]
values2 = values2_series.loc[common_idx]
plot_correlation(values1, values2, "SCZ", "HIQ ASD HCT", dpi=80)

# %%
common_keys = set(HIQ_GW.keys()) & set(SCZ_GW.keys())

# Create new dicts excluding common keys
HIQ_GW_unique = {k:v for k,v in HIQ_GW.items() if k not in common_keys}
SCZ_GW_unique = {k:v for k,v in SCZ_GW.items() if k not in common_keys}
print(len(HIQ_GW), len(SCZ_GW))
print(len(HIQ_GW_unique), len(SCZ_GW_unique))

# %% hidden=true
SCZ_Z2_Bias_unique = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, SCZ_GW_unique)
SCZ_Z2_Bias_unique = AnnotateCTDat(SCZ_Z2_Bias_unique, Anno)
ASD_HIQ_Bias_unique = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, HIQ_GW_unique)
ASD_HIQ_Bias_unique = AnnotateCTDat(ASD_HIQ_Bias_unique, Anno)
SCZ_Z2_Bias_unique.to_csv("{}/HCT.SCZ61.unique.csv".format(Bias_Save_Dir))
ASD_HIQ_Bias_unique.to_csv("{}/HCT.ASD_HIQ.unique.csv".format(Bias_Save_Dir))

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
values1 = SCZ_Z2_Bias_unique.sort_index()[EFFECT_Col].values
values2 = ASD_HIQ_Bias_unique.sort_index()[EFFECT_Col].values
plot_correlation(values2, values1, "SCZ", "HIQ ASD HCT", dpi=80)

# %% [markdown] heading_collapsed=true
# ## 22q.11 Bias

# %%
X22q_GW = Fil2Dict("../dat/GeneWeights/X22q.gw.csv")
# Exclude key 6899 from X22q_GW
#X22q_GW = {k: v for k, v in X22q_GW.items() if k != 6899}
X22q_GW_Mouse = Fil2Dict("{}/X22q.mousemodel.gw.csv".format(GeneWeightDIR))

# %% hidden=true
X22q_Z2_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, X22q_GW)
X22q_Z2_Bias = AnnotateCTDat(X22q_Z2_Bias, Anno)
X22q_Z2_Bias.to_csv("{}/HCT.X22q.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(X22q_Z2_Bias, "22q11.2 Complete deletion")


# %%
# Get subset of rows from VIP_Anno that exist in X22q_Z2_Bias
common_indices = X22q_Z2_Bias.index.intersection(VIP_Anno.index)

X22q_Z2_Bias_sub = X22q_Z2_Bias.loc[common_indices].copy()
X22q_Z2_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_Z2_Bias_sub)

# %% hidden=true
X22q_Z2_Mouse_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, X22q_GW_Mouse)
X22q_Z2_Mouse_Bias = AnnotateCTDat(X22q_Z2_Mouse_Bias, Anno)
X22q_Z2_Mouse_Bias.to_csv("{}/HCT.X22q.mousemodel.csv".format(Bias_Save_Dir))


# %%
SuperClusterBias_BoxPlot(X22q_Z2_Mouse_Bias, "22q11.2 Mouse deletion")
# SuperClusterBias_BoxPlot(X22q_Z2_Mouse_Bias_PT, "22q11.2 Mouse deletion")

# %%
common_indices = X22q_Z2_Mouse_Bias.index.intersection(VIP_Anno.index)
X22q_Z2_Mouse_Bias_sub = X22q_Z2_Mouse_Bias.loc[common_indices].copy()
X22q_Z2_Mouse_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_Z2_Mouse_Bias_sub)
#plot_vip_effect_comparison(X22q_Z2_Mouse_Bias_sub, effect_col = "EFFECT_Adj")

# %%
X22q_Genes = list(X22q_GW.keys())
X22q_Genes = [g for g in X22q_Genes if g in HumanCT_Z2_HCT.index.values]
#X22q_Genes = [g for g in X22q_Genes if g in HumanCT_Z2_HCT.index.values if g != 6899]
print(len(X22q_Genes))

# %%
HumanCT_OverallEXP = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/HumanCT.MatchDF.csv", index_col=0)
HumanCT_Exp_22q = HumanCT_OverallEXP.loc[X22q_Genes]

# %%
HumanCT_Exp_22q.head(2)

# %%
HumanCT_Exp_22q.loc[6899,:]

# %%
Anno["Supercluster"].unique()

# %%
#HumanCT_Exp_22q.hist(bins=10)
SuperCluster = "Choroid plexus"
SuperCluster_idx = Anno[Anno["Supercluster"] == SuperCluster].index.values

# Calculate expression and specificity values for each gene
log_exp_values = []
spec_values = []
for g in X22q_Genes:
    LogExp = np.log10(HumanCT_Exp_22q.loc[g, "Exp"] + 1)
    Spec = HumanCT_Z2_HCT.loc[g, SuperCluster_idx].mean()
    log_exp_values.append(LogExp)
    spec_values.append(Spec)

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(log_exp_values, spec_values, alpha=0.5)
plt.xlabel('Log10 Expression')
plt.ylabel(f'{SuperCluster} Specificity')
plt.title(f'Expression vs {SuperCluster} Specificity for 22q Genes')
plt.tight_layout()
plt.show()

# %%
#HumanCT_Exp_22q.hist(bins=10)
SuperCluster = "CGE interneuron"
SuperCluster_idx = Anno[Anno["Supercluster"] == SuperCluster].index.values

# Calculate expression and specificity values for each gene
log_exp_values = []
spec_values = []
for g in X22q_Genes:
    LogExp = np.log10(HumanCT_Exp_22q.loc[g, "Exp"] + 1)
    Spec = HumanCT_Z2_HCT.loc[g, SuperCluster_idx].mean()
    log_exp_values.append(LogExp)
    spec_values.append(Spec)

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(log_exp_values, spec_values, alpha=0.5)
plt.xlabel('Log10 Expression')
plt.ylabel(f'{SuperCluster} Specificity')
plt.title(f'Expression vs {SuperCluster} Specificity for 22q Genes')
plt.tight_layout()
plt.show()

# %%
#HumanCT_Exp_22q.hist(bins=10)
SuperCluster = "CGE interneuron"
SuperCluster_idx = Anno[Anno["Supercluster"] == SuperCluster].index.values

target_gene = 6899
# Calculate expression and specificity values for each gene

Spec = HumanCT_Z2_HCT.loc[target_gene, SuperCluster_idx]
print(Spec)


# %%
ExpL_Cut = 100000
#ExpL_Cut = 33000
HighExp_22q_Genes = HumanCT_Exp_22q[HumanCT_Exp_22q["Exp"] > ExpL_Cut].index.values
print(len(HighExp_22q_Genes))
HighExp_22q_GW = dict(zip(HighExp_22q_Genes, np.ones(len(HighExp_22q_Genes))))
HighExp_22q_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, HighExp_22q_GW)
HighExp_22q_Bias = AnnotateCTDat(HighExp_22q_Bias, Anno)
#HighExp_22q_Bias = AdjustClusterMean(HighExp_22q_Bias, HumanCT_Z2_HCT_ColMean)
HighExp_22q_Bias.to_csv("{}/HCT.22q.HighExp.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(HighExp_22q_Bias, "22q High Exp")
#SuperClusterBias_BoxPlot(HighExp_22q_Bias, "22q High Exp", EffectCol="EFFECT_Adj")


# %%
common_indices = HighExp_22q_Bias.index.intersection(VIP_Anno.index)
HighExp_22q_Bias_sub = HighExp_22q_Bias.loc[common_indices].copy()
HighExp_22q_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(HighExp_22q_Bias_sub)

# %%
VIP_ExpCut = 1
VIP_Pos = VIP_Anno[VIP_Anno["VIP"]>VIP_ExpCut].index.values
VIP_Neg = VIP_Anno[VIP_Anno["VIP"]<VIP_ExpCut].index.values
print(len(VIP_Pos), len(VIP_Neg))


dat1 = HumanCT_Z2_HCT.loc[HighExp_22q_Genes, VIP_Pos].mean(axis=1)
dat2 = HumanCT_Z2_HCT.loc[HighExp_22q_Genes, VIP_Neg].mean(axis=1)

print(len(dat1), len(dat2))

# Create boxplot with individual points
plt.figure(figsize=(8,5))
data = [dat1, dat2]
bp = plt.boxplot(data, labels=["VIP+", "VIP-"], showfliers=True)

# Add individual points
for i, d in enumerate([dat1, dat2]):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.3, s=20)

plt.ylabel("Expression Z-score") 
plt.title("Expression in VIP+ vs VIP- cells")

# Perform Mann-Whitney U test
stat, pval = scipy.stats.mannwhitneyu(dat1, dat2, alternative='two-sided')
print(f"Mann-Whitney U test p-value: {pval:.2e}")


# %%
VIP_ExpCut = 1
VIP_Pos = VIP_Anno[VIP_Anno["VIP"]>VIP_ExpCut].index.values
VIP_Neg = VIP_Anno[VIP_Anno["VIP"]<VIP_ExpCut].index.values
print(len(VIP_Pos), len(VIP_Neg))

# Get expression for each gene in both VIP+ and VIP- cells
gene_expr_vip_pos = HumanCT_Z2_HCT.loc[HighExp_22q_Genes, VIP_Pos].mean(axis=1)
gene_expr_vip_neg = HumanCT_Z2_HCT.loc[HighExp_22q_Genes, VIP_Neg].mean(axis=1)

print(len(gene_expr_vip_pos), len(gene_expr_vip_neg))

# Create boxplot with individual points
plt.figure(figsize=(8,5))
data = [gene_expr_vip_pos, gene_expr_vip_neg]

# Use tick_labels instead of deprecated labels parameter
bp = plt.boxplot(data, tick_labels=["VIP+", "VIP-"], showfliers=True)

# Add individual points with lines connecting paired values
for i in range(len(gene_expr_vip_pos)):
    # Convert gene expression values to numpy arrays to avoid potential issues
    pos_val = np.array(gene_expr_vip_pos)[i]
    neg_val = np.array(gene_expr_vip_neg)[i]
    
    # Check for NaN values before plotting
    if not (np.isnan(pos_val) or np.isnan(neg_val)):
        plt.plot([1, 2], [pos_val, neg_val], 
                 color='gray', alpha=0.3, linewidth=0.5)
        plt.scatter([1, 2], [pos_val, neg_val], 
                    alpha=0.3, s=20)

plt.ylabel("Expression Z-score") 
plt.title("Expression in VIP+ vs VIP- cells")

# Remove NaN values before statistical test
mask = ~(np.isnan(gene_expr_vip_pos) | np.isnan(gene_expr_vip_neg))
pos_clean = np.array(gene_expr_vip_pos)[mask]
neg_clean = np.array(gene_expr_vip_neg)[mask]

# Perform Wilcoxon signed-rank test (paired, one-sided)
if len(pos_clean) > 0 and len(neg_clean) > 0:
    stat, pval = scipy.stats.wilcoxon(pos_clean, neg_clean, alternative='greater')
    print(f"Wilcoxon signed-rank test p-value (VIP+ > VIP-): {pval:.2e}")
else:
    print("Not enough valid data for statistical test")


# %%
merged_loeuf_df = pd.read_csv("../dat/gnomad.LOEUF.merged.csv", index_col=0)
merged_loeuf_df.index = merged_loeuf_df.index.astype(int)
merged_loeuf_df.head(2)

# %%
LOEUF_22q_Genes = merged_loeuf_df.loc[merged_loeuf_df.index.intersection(X22q_Genes)]
print(len(LOEUF_22q_Genes))

# %%
LOEUF_cut = 1
#ExpL_Cut = 33000
high_loeuf_22q_genes = LOEUF_22q_Genes[LOEUF_22q_Genes["LOEUF"] < LOEUF_cut].index.values
high_loeuf_22q_genes_rev = LOEUF_22q_Genes[LOEUF_22q_Genes["LOEUF"] > LOEUF_cut].index.values
print(len(high_loeuf_22q_genes))
high_loeuf_22q_gw = dict(zip(high_loeuf_22q_genes, np.ones(len(high_loeuf_22q_genes))))
high_loeuf_22q_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, high_loeuf_22q_gw)
high_loeuf_22q_bias = AnnotateCTDat(high_loeuf_22q_bias, Anno)
#high_loeuf_22q_bias = AdjustClusterMean(high_loeuf_22q_bias, HumanCT_Z2_HCT_ColMean)
high_loeuf_22q_bias.to_csv("{}/HCT.22q.HighLOEUF.csv".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(high_loeuf_22q_bias, "22q High LOEUF")
#SuperClusterBias_BoxPlot(high_loeuf_22q_bias, "22q High Exp", EffectCol="EFFECT_Adj")


# %%
print([Entrez2Symbol[gid] for gid in high_loeuf_22q_genes_rev])
print([Entrez2Symbol[gid] for gid in high_loeuf_22q_genes])

# %%
common_indices = high_loeuf_22q_bias.index.intersection(VIP_Anno.index)
high_loeuf_22q_bias_sub = high_loeuf_22q_bias.loc[common_indices].copy()
high_loeuf_22q_bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(high_loeuf_22q_bias_sub)

# %% [markdown] heading_collapsed=true
# ## DDD

# %%
# DDD_top61_GW = Fil2Dict("{}/DDD.top61.gw.csv".format(GeneWeightDIR))
# DDD_hc_Bias_top61 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DDD_top61_GW)
# DDD_hc_Bias_top61 = AnnotateCTDat(DDD_hc_Bias_top61, Anno)
# DDD_hc_Bias_top61.to_csv("{}/HCT.DDDHC.Spec.top61.csv".format(Bias_Save_Dir))

# %%
DDD_top61_GW = Fil2Dict("{}/DDD.top61.gw.bgmr.csv".format(GeneWeightDIR))
DDD_hc_Bias_top61 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DDD_top61_GW)
DDD_hc_Bias_top61 = AnnotateCTDat(DDD_hc_Bias_top61, Anno)
DDD_hc_Bias_top61.to_csv("{}/HCT.DDDHC.Spec.top61".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(DDD_hc_Bias_top61, "DDD HCT top61")

# %%
DDD_top61_GW = Fil2Dict("{}/DDD.hc.gw".format(GeneWeightDIR))
DDD_hc_Bias_top61 = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DDD_top61_GW)
DDD_hc_Bias_top61 = AnnotateCTDat(DDD_hc_Bias_top61, Anno)
DDD_hc_Bias_top61.to_csv("{}/HCT.DDDHC.Spec.hc".format(Bias_Save_Dir))
SuperClusterBias_BoxPlot(DDD_hc_Bias_top61, "DDD HCT top61")

# %%
common_indices = DDD_hc_Bias_top61.index.intersection(VIP_Anno.index)
DDD_hc_Bias_top61_sub = DDD_hc_Bias_top61.loc[common_indices].copy()
DDD_hc_Bias_top61_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
#plot_vip_effect_comparison(DDD_hc_Bias_top61_sub, effect_col = "EFFECT_Adj")
plot_vip_effect_comparison(DDD_hc_Bias_top61_sub, effect_col = "EFFECT")

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# # Correlation with something

# %%
EFFECT_Col = "EFFECT"
#EFFECT_Adj_Col = "EFFECT_Adj"
values1 = Anno.sort_index()["Total UMI"].values
values2 = HIQ_Z2_Bias.sort_index()[EFFECT_Col].values
plot_correlation(values2, values1, "UMI", "HIQ ASD HCT", dpi=80)

values1 = Anno.sort_index()["Total UMI"].values
values2 = SCZ_Z2_Bias.sort_index()[EFFECT_Col].values
plot_correlation(values2, values1, "UMI", "SCZ HCT", dpi=80)

values1 = Anno.sort_index()["Total UMI"].values
values2 = LIQ_Z2_Bias.sort_index()[EFFECT_Col].values
plot_correlation(values2, values1, "UMI", "LIQ ASD HCT", dpi=80)

values1 = Anno.sort_index()["Total UMI"].values
values2 = DDD_hc_Bias_top61.sort_index()[EFFECT_Col].values
plot_correlation(values2, values1, "UMI", "DDD HCT top61", dpi=80)

# %%
Jon_PSD_Spec_percentile = pd.read_csv("/home/jw3514/Work/NeuralP/dat/Bias/Other/Jon_PSD_Spec_percentile.csv", index_col=0)
Jon_PSD_Spec_percentile.head(2)

# %%
Jon_PSD_list = "88 102 107 114 118 119 120 150 375 476 491 538 575 577 613 659 664 747 775 815 832 1000 1020 1128 1136 1499 1500 1501 1612 1627 1630 1739 1740 1741 1742 1785 1821 1838 1855 1948 1949 2039 2043 2059 2066 2171 2185 2332 2534 2596 2785 2852 2890 2891 2892 2893 2894 2895 2897 2898 2899 2900 2901 2902 2903 2904 2905 2906 2911 2913 2915 3184 3188 3337 3631 3646 3708 3756 3836 4038 4131 4218 4355 4744 4804 4842 4897 4905 4915 4985 5028 5062 5063 5064 5071 5093 5094 5142 5170 5582 5590 5621 5662 5728 5800 5802 6009 6128 6129 6132 6136 6156 6169 6175 6188 6207 6208 6222 6223 6230 6232 6252 6334 6457 6536 6543 6546 6547 6695 6711 6790 6792 6801 6853 6854 7074 7248 7249 7428 7732 7779 8087 8224 8440 8502 8516 8661 8777 8787 8825 8831 8851 8898 8927 8997 9026 9045 9101 9148 9162 9194 9201 9228 9229 9231 9419 9454 9455 9456 9463 9478 9495 9513 9743 9746 9762 9829 9863 9867 9890 9912 9921 9922 10006 10142 10243 10280 10313 10368 10369 10458 10486 10505 10509 10611 10636 11122 11178 11331 11346 22829 22849 22865 22866 22871 22883 22941 22986 22997 23043 23208 23229 23237 23362 23380 23413 23426 23513 23542 23623 23705 23767 25945 25978 26012 26037 26045 26052 27020 27091 27092 27185 27445 28964 28988 29102 29904 50488 50944 51104 51201 51225 54477 54487 54583 54910 55327 55450 55607 55737 56899 56924 57120 57142 57479 57502 57537 57554 57622 57679 57689 58489 58512 59283 59284 64084 64101 64130 64506 66000 78999 79414 79870 79953 80315 80725 80758 80852 80863 81831 81832 81926 84062 84435 84687 85358 85461 94030 112476 114798 116443 116444 145581 145773 146395 152404 158866 160622 200933 201191 254263 339451 347731 373509 388135 388336 392862 400745 401190 440073 440829 642938 729956 729993 100131897"
Jon_PSD_list = [int(x) for x in Jon_PSD_list.split(" ")]
Jon_PSD_GW = dict(zip(Jon_PSD_list, np.ones(len(Jon_PSD_list))))
Jon_PSD_Spec = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, Jon_PSD_GW)

# %%
values1 = Jon_PSD_Spec.sort_index()['EFFECT'].values
values2 = Jon_PSD_Spec_percentile.sort_index()['EFFECT'].values
plot_correlation(values2, values1, "PSD Spec", "PSD Spec Percentile", dpi=80)

# %%
#Jon_PSD_Spec_percentile

# %%
values1 = Anno.sort_index()["Total UMI"].values
values2 = Jon_PSD_Spec_percentile.sort_index()['EFFECT'].values
plot_correlation(values2, values1, "UMI", "PSD", dpi=80)

# %%
#Eff_col = "EFFECT_Adj"
Eff_col = "EFFECT"
values1 = Jon_PSD_Spec.sort_index()['EFFECT'].values
values2 = SCZ_Z2_Bias.sort_index()[Eff_col].values
plot_correlation(values1, values2, "PSD Spec", "SCZ Spec", None, dpi=80)

values1 = Jon_PSD_Spec.sort_index()['EFFECT'].values
values2 = HIQ_Z2_Bias.sort_index()[Eff_col].values
plot_correlation(values1, values2, "PSD Spec", "HIQ ASD Spec", None, dpi=80)

values1 = Jon_PSD_Spec.sort_index()['EFFECT'].values
values2 = LIQ_Z2_Bias.sort_index()[Eff_col].values
plot_correlation(values1, values2, "PSD Spec", "LIQ ASD Spec", None, dpi=80)

values1 = Jon_PSD_Spec.sort_index()['EFFECT'].values
values2 = DDD_hc_Bias_top61.sort_index()[Eff_col].values
plot_correlation(values1, values2, "PSD Spec", "DDD HCT Spec", None, dpi=80)

# %%
Eff_col = "EFFECT"
values1 = SCZ_Z2_Bias.sort_index()[Eff_col].values
values2 = HIQ_Z2_Bias.sort_index()[Eff_col].values
plot_correlation(values1, values2, "SCZ HCT", "HIQ ASD HCT", None, dpi=80)

# %%
NC_Bias = pd.read_csv("/home/jw3514/Work/NeuralP/notebooks/data/NC_Bias_Beta.csv", index_col=0)
#NC_Bias = pd.read_csv("/home/jw3514/Work/NeuralP/notebooks/NC_Bias_Beta.csv", index_col=0)
NC_Bias.head(2)

# %%
#NC_Bias.columns

# %%
values1 = SCZ_Z2_Bias.sort_index()["EFFECT"].values
values2 = NC_Bias.sort_index()["scz2022"].values
plot_correlation(values1, values2, "SCZ Rare", "SCZ GWAS", None, dpi=80)

# %%
values1 = NC_Bias.sort_index()["educational_attainment"].values
values2 = NC_Bias.sort_index()["scz2022"].values
plot_correlation(values1, values2, "EDU GWAS", "SCZ GWAS", None, dpi=80)


# %%
def compute_trait_correlations(values1, nc_bias_df):
    """
    Compute Spearman correlations between values1 and each column in nc_bias_df.
    Optionally, also compute correlations restricted to indices neur_idx.
    Returns a DataFrame sorted by Spearman_r descending.
    """
    results = []
    for effect in nc_bias_df.columns:
        values2 = nc_bias_df.sort_index()[effect].values
        r, p = spearmanr(values1, values2)
        r_neur, p_neur = spearmanr(values1[Neur_idx], values2[Neur_idx])
        results.append({"Trait": effect, "r_all": r, "p_value": p, "r_neur":r_neur, "p_value_neur":p_neur})
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="r_all", ascending=False)
    return results_df

values1 = SCZ_Z2_Bias.sort_index()["EFFECT"].values
scz_results_df = compute_trait_correlations(values1, NC_Bias)

# %%
scz_results_df

# %%
values1 = HIQ_Z2_Bias.sort_index()["EFFECT"].values
HIQ_results_df = compute_trait_correlations(values1, NC_Bias)

# %%
HIQ_results_df.sort_values(by="r_neur", ascending=False)

# %%
values1 = LIQ_Z2_Bias.sort_index()["EFFECT"].values
results = []
for effect in NC_Bias.columns:
    values2 = NC_Bias.sort_index()[effect].values
    r, p = spearmanr(values1, values2)
    results.append({"Trait": effect, "Spearman_r": r, "p_value": p})
LIQ_results_df = pd.DataFrame(results)
LIQ_results_df = LIQ_results_df.sort_values(by="Spearman_r", ascending=False)

# %%
LIQ_results_df

# %%
values1 = X22q_Z2_Bias.sort_index()["EFFECT"].values
results = []
for effect in NC_Bias.columns:
    values2 = NC_Bias.sort_index()[effect].values
    r, p = spearmanr(values1, values2)
    results.append({"Trait": effect, "Spearman_r": r, "p_value": p})
X22q_results_df = pd.DataFrame(results)
X22q_results_df = X22q_results_df.sort_values(by="Spearman_r", ascending=False)

# %%
X22q_results_df

# %%
values1 = DDD_hc_Bias_top61.sort_index()["EFFECT"].values
results = []
for effect in NC_Bias.columns:
    values2 = NC_Bias.sort_index()[effect].values
    r, p = spearmanr(values1, values2)
    results.append({"Trait": effect, "Spearman_r": r, "p_value": p})
DDD_results_df = pd.DataFrame(results)
DDD_results_df = DDD_results_df.sort_values(by="Spearman_r", ascending=False)

# %%
DDD_results_df

# %%

# %% [markdown]
# # Other data sets

# %% [markdown] heading_collapsed=true
# #### VNR

# %% hidden=true
topN = 61
VNR_Pos_GW = Fil2Dict("{}/UKBB_VNR_Pos_GW_{}.csv".format(GeneWeightDIR, topN))
VNR_Neg_GW = Fil2Dict("{}/UKBB_VNR_Neg_GW_{}.csv".format(GeneWeightDIR, topN))
VNR_Pos_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, VNR_Pos_GW)
VNR_Pos_Bias = AnnotateCTDat(VNR_Pos_Bias, Anno)
VNR_Neg_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, VNR_Neg_GW)
VNR_Neg_Bias = AnnotateCTDat(VNR_Neg_Bias, Anno)
VNR_Pos_Bias.to_csv("{}/HCT.VNR.Pos.top61.csv".format(Bias_Save_Dir))
VNR_Neg_Bias.to_csv("{}/HCT.VNR.Neg.top61.csv".format(Bias_Save_Dir))



# %%
SuperClusterBias_BoxPlot(VNR_Pos_Bias, "VNR_Pos_Bias")
SuperClusterBias_BoxPlot(VNR_Neg_Bias, "VNR_Neg_Bias")

# %%
topN = 61
EDU_Pos_GW = Fil2Dict("{}/UKBB_EDU_Pos_GW_{}.csv".format(GeneWeightDIR, topN))
EDU_Neg_GW = Fil2Dict("{}/UKBB_EDU_Neg_GW_{}.csv".format(GeneWeightDIR, topN))
EDU_Pos_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, EDU_Pos_GW)
EDU_Pos_Bias = AnnotateCTDat(EDU_Pos_Bias, Anno)
EDU_Neg_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, EDU_Neg_GW)
EDU_Neg_Bias = AnnotateCTDat(EDU_Neg_Bias, Anno)
EDU_Pos_Bias.to_csv("{}/HCT.EDU.Pos.top61.csv".format(Bias_Save_Dir))
EDU_Neg_Bias.to_csv("{}/HCT.EDU.Neg.top61.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(EDU_Pos_Bias, "EDU_Pos_Bias")
SuperClusterBias_BoxPlot(EDU_Neg_Bias, "EDU_Neg_Bias")

# %%
topN = 61
RT_Pos_GW = Fil2Dict("{}/UKBB_RT_Pos_GW_{}.csv".format(GeneWeightDIR, topN))
RT_Neg_GW = Fil2Dict("{}/UKBB_RT_Neg_GW_{}.csv".format(GeneWeightDIR, topN))
RT_Pos_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, RT_Pos_GW)
RT_Pos_Bias = AnnotateCTDat(RT_Pos_Bias, Anno)
RT_Neg_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, RT_Neg_GW)
RT_Neg_Bias = AnnotateCTDat(RT_Neg_Bias, Anno)
RT_Pos_Bias.to_csv("{}/HCT.RT.Pos.top61.csv".format(Bias_Save_Dir))
RT_Neg_Bias.to_csv("{}/HCT.RT.Neg.top61.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(RT_Pos_Bias, "RT_Pos_Bias")
SuperClusterBias_BoxPlot(RT_Neg_Bias, "RT_Neg_Bias")

# %%
VNR_Neg_Bias

# %%
common_indices = VNR_Neg_Bias.index.intersection(VIP_Anno.index)
VNR_Neg_Bias_sub = VNR_Neg_Bias.loc[common_indices].copy()
VNR_Neg_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(VNR_Neg_Bias_sub)
#plot_vip_effect_comparison(X22q_Z2_Mouse_Bias_sub, effect_col = "EFFECT_Adj")

# %%
common_indices = EDU_Neg_Bias.index.intersection(VIP_Anno.index)
EDU_Neg_Bias_sub = EDU_Neg_Bias.loc[common_indices].copy()
EDU_Neg_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(EDU_Neg_Bias_sub)

# %%
common_indices = RT_Pos_Bias.index.intersection(VIP_Anno.index)
RT_Pos_Bias_sub = RT_Pos_Bias.loc[common_indices].copy()
RT_Pos_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(RT_Pos_Bias_sub)

# %% [markdown]
# ## BP

# %%
BP_GW = Fil2Dict("/home/jw3514/Work/CellType_Psy/dat/Archive/Bipolar.top61.gw.csv")
BP_Z2_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, BP_GW)
BP_Z2_Bias = AnnotateCTDat(BP_Z2_Bias, Anno)
#BP_Z2_Bias.to_csv("../dat/Bias/CTBias.BP.Z2.HCT.csv")
BP_Z2_Bias.to_csv("{}/HCT.BP.Z2.HCT.csv".format(Bias_Save_Dir))


# %%
SuperClusterBias_BoxPlot(BP_Z2_Bias, "BP HCT")

# %% [markdown]
# ## CNVs

# %%
genes_3q29 = [
    "BDH1",       # 3-hydroxybutyrate dehydrogenase 1
    "FBXO45",     # F-box protein 45
    "SENP5",      # SUMO specific peptidase 5
    "UBXN7",      # UBX domain protein 7
    "WDR53",      # WD repeat domain 53
    "CEP19",      # centrosomal protein 19
    "DLG1",       # discs large MAGUK scaffold protein 1
    "DYNLT2B",    # dynein light chain Tctex-type 2B
    "MELTF",      # melanotransferrin
    "NRROS",      # negative regulator of reactive oxygen species
    "NCBP2",      # nuclear cap binding protein subunit 2
    "PAK2",       # p21 (RAC1) activated kinase 2
    "PCYT1A",     # phosphate cytidylyltransferase 1A, choline
    "PIGX",       # phosphatidylinositol glycan anchor biosynthesis class X
    "PIGZ",       # phosphatidylinositol glycan anchor biosynthesis class Z
    "RNF168",     # ring finger protein 168
    "SMCO1",      # single-pass membrane protein with coiled-coil domains 1
    "SLC51A",     # solute carrier family 51 subunit alpha
    "TFRC",       # transferrin receptor
    "TM4SF19",    # transmembrane 4 L six family member 19
    "ZDHHC19"     # zinc finger DHHC-type palmitoyltransferase 19
]
genes_3q29_entrez = [int(GeneSymbol2Entrez[x]) for x in genes_3q29 if x in GeneSymbol2Entrez]
genes_3q29_GW = dict(zip(genes_3q29_entrez, [1]*len(genes_3q29_entrez)))
Dict2Fil(genes_3q29_GW, "/home/jw3514/Work/CellType_Psy/dat/GeneWeights/3q29.gw.csv")
genes_3q29_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, genes_3q29_GW)
genes_3q29_Bias = AnnotateCTDat(genes_3q29_Bias, Anno)
#genes_3q29_Bias.to_csv("../dat/Bias/CTBias.3q29.Z2.HCT.csv")
genes_3q29_Bias.to_csv("{}/HCT.3q29.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(genes_3q29_Bias, "3q29 HCT")

# %%
genes_3q29_Bias.head(30)

# %%
X16p_df = pd.read_excel("~/Work/data/13073_2021_972_MOESM2_ESM.xlsx", sheet_name="16p11.2")
Gene16p11 = X16p_df[X16p_df["Type"]=="protein_coding"]["Gene_name"].values
Gene16p11 = [int(GeneSymbol2Entrez[x]) for x in Gene16p11 if x in GeneSymbol2Entrez.keys()]
print(len(Gene16p11))
X16pGW = dict(zip(Gene16p11, [1]*len(Gene16p11)))
Dict2Fil(X16pGW, "/home/jw3514/Work/CellType_Psy/dat/GeneWeights/16p11.gw.csv")

# %%
genes_16p11_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, X16pGW)
genes_16p11_Bias = AnnotateCTDat(genes_16p11_Bias, Anno)
#genes_16p11_Bias.to_csv("../dat/Bias/CTBias.16p11.Z2.HCT.csv")
genes_16p11_Bias.to_csv("{}/HCT.16p11.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(genes_16p11_Bias, "16p11.2 HCT")

# %%
common_indices = genes_16p11_Bias.index.intersection(VIP_Anno.index)
genes_16p11_Bias_sub = genes_16p11_Bias.loc[common_indices].copy()
genes_16p11_Bias_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(genes_16p11_Bias_sub)

# %%
print(len(genes_3q29_GW))
print(len(X16pGW))
print(len(X22q_GW))

# %%
ADHD_gene = [
    "MAP1A",
    "ANO8",
    "ANK2",
    "SPTBN1",
    "LRP1",
    "PIK3CD",
    "NAA15",
    "TNRC6C",
    "TSC22D1",
    "WNT1",
    "KCNA2",
    "GIGYF1",
    "CNOT3",
    "TRIM46",
    "KDR",
    "EIF3G",
    "CACNA1D",
    "PAX5",
    "UBE2S",
    "ALX1"
]

ADHD_entrez = [int(GeneSymbol2Entrez[x]) for x in ADHD_gene if x in GeneSymbol2Entrez]
ADHD_GW = dict(zip(ADHD_entrez, [1]*len(ADHD_entrez)))
Dict2Fil(ADHD_GW, "/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ADHD.gw.csv")
ADHD_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ADHD_GW)
ADHD_Bias = AnnotateCTDat(ADHD_Bias, Anno)
#ADHD_Bias.to_csv("../dat/Bias/CTBias.ADHD.Z2.HCT.csv")
ADHD_Bias.to_csv("{}/HCT.ADHD.csv".format(Bias_Save_Dir))

# %%
SuperClusterBias_BoxPlot(ADHD_Bias, "ADHD HCT")

# %%
ADHD_Bias.head(20)

# %% [markdown]
# # Neg Ctrl

# %%

HDL_C_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_continuous-30760-both_sexes--irnt_2025_11_25_15_49_42.csv", index_col=0)
T2D_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_categorical-T2D_custom-both_sexes--custom_2025_11_25_16_14_20.csv", index_col=0)
IBD_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_categorical-IBD_custom2-both_sexes--custom_2025_11_25_16_13_49.csv", index_col=0)
hba1c_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_continuous-30750-both_sexes--irnt_2025_11_25_16_32_28.csv", index_col=0)
Parkinson_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_icd_first_occurrence-131022-both_sexes--_2025_11_25_16_56_06.csv", index_col=0)
Alzheimer_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_icd_first_occurrence-131036-both_sexes--_2025_11_25_17_00_04.csv", index_col=0)

Alanine_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_continuous-30620-both_sexes--irnt_2025_12_16_20_18_07.csv", index_col=0)
redbloodcell_GeneTest = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/gene-burden-results-exomes_pLoF_continuous-30070-both_sexes--irnt_2025_12_16_20_41_11.csv", index_col=0)

# %%
df_annot = HDL_C_GeneTest.copy()
df_annot["EntrezID"] = df_annot.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)

# %%
top_entrez_ids = df_annot.head(61)["EntrezID"]
top_entrez_ids_no_nan = top_entrez_ids[top_entrez_ids.notna()].astype(int)[:200]
HDL_C_top_gw = dict(zip(list(top_entrez_ids_no_nan), [1]*len(top_entrez_ids_no_nan)))
#Dict2Fil(dict(chdl_top_gw), "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/CTRL.top61.gw.csv")
HDL_C_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, HDL_C_top_gw)
HDL_C_top_bias = AnnotateCTDat(HDL_C_top_bias, Anno)
#chdl_top_bias.to_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/CTRL/CTRL.top61.csv")
SuperClusterBias_BoxPlot(HDL_C_top_bias, "CTRL HCT")


# %%
# T2D cell type bias
df_annot_t2d = T2D_GeneTest.copy()
df_annot_t2d["EntrezID"] = df_annot_t2d.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)


# %%
T2D_top_entrez = df_annot_t2d[df_annot_t2d["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
T2D_top_gw = dict(zip(list(T2D_top_entrez), [1]*len(T2D_top_entrez))) 
T2D_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, T2D_top_gw)
T2D_top_bias = AnnotateCTDat(T2D_top_bias, Anno)
SuperClusterBias_BoxPlot(T2D_top_bias, "T2D HCT")


# %%
# IBD cell type bias
df_annot_ibd = IBD_GeneTest.copy()
df_annot_ibd["EntrezID"] = df_annot_ibd.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)


# %%
IBD_top_entrez = df_annot_ibd[df_annot_ibd["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
IBD_top_gw = dict(zip(list(IBD_top_entrez), [1]*len(IBD_top_entrez))) 
IBD_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, IBD_top_gw)
IBD_top_bias = AnnotateCTDat(IBD_top_bias, Anno)
SuperClusterBias_BoxPlot(IBD_top_bias, "IBD HCT")


# %%
# hba1c cell type bias
df_annot_hba1c = hba1c_GeneTest.copy()
df_annot_hba1c["EntrezID"] = df_annot_hba1c.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)


# %%
hba1c_top_entrez = df_annot_hba1c[df_annot_hba1c["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
hba1c_top_gw = dict(zip(list(hba1c_top_entrez), [1]*len(hba1c_top_entrez))) 
hba1c_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, hba1c_top_gw)
hba1c_top_bias = AnnotateCTDat(hba1c_top_bias, Anno)
SuperClusterBias_BoxPlot(hba1c_top_bias, "hba1c HCT")

# %%
# Parkinson cell type bias
df_annot_parkinson = Parkinson_GeneTest.copy()
df_annot_parkinson["EntrezID"] = df_annot_parkinson.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)


# %%
Parkinson_top_entrez = df_annot_parkinson[df_annot_parkinson["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
Parkinson_top_gw = dict(zip(list(Parkinson_top_entrez), [1]*len(Parkinson_top_entrez))) 
Parkinson_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, Parkinson_top_gw)
Parkinson_top_bias = AnnotateCTDat(Parkinson_top_bias, Anno)
SuperClusterBias_BoxPlot(Parkinson_top_bias, "Parkinson HCT")


# %%
# Alanine cell type bias
df_annot_alanine = Alanine_GeneTest.copy()
df_annot_alanine["EntrezID"] = df_annot_alanine.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)

Alanine_top_entrez = df_annot_alanine[df_annot_alanine["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
Alanine_top_gw = dict(zip(list(Alanine_top_entrez), [1]*len(Alanine_top_entrez))) 
Alanine_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, Alanine_top_gw)
Alanine_top_bias = AnnotateCTDat(Alanine_top_bias, Anno)
SuperClusterBias_BoxPlot(Alanine_top_bias, "Alanine HCT")


# %%
# RedBloodCell cell type bias
df_annot_redbloodcell = redbloodcell_GeneTest.copy()
df_annot_redbloodcell["EntrezID"] = df_annot_redbloodcell.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)

RedBloodCell_top_entrez = df_annot_redbloodcell[df_annot_redbloodcell["EntrezID"].notna()].head(61)["EntrezID"].astype(int)
RedBloodCell_top_gw = dict(zip(list(RedBloodCell_top_entrez), [1]*len(RedBloodCell_top_entrez))) 
RedBloodCell_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, RedBloodCell_top_gw)
RedBloodCell_top_bias = AnnotateCTDat(RedBloodCell_top_bias, Anno)
SuperClusterBias_BoxPlot(RedBloodCell_top_bias, "RedBloodCell HCT")


# %%
# ALZ (Alzheimer's) cell type bias
df_annot_alz = Alzheimer_GeneTest.copy()
df_annot_alz["EntrezID"] = df_annot_alz.index.map(lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez[x]) else None)


# %%
df_annot_alz

# %%
ALZ_top_entrez = df_annot_alz[df_annot_alz["EntrezID"].notna()].head(60)["EntrezID"].astype(int)
ALZ_top_gw = dict(zip(list(ALZ_top_entrez), [1]*len(ALZ_top_entrez))) 
ALZ_top_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ALZ_top_gw)
ALZ_top_bias = AnnotateCTDat(ALZ_top_bias, Anno)
SuperClusterBias_BoxPlot(ALZ_top_bias, "ALZ HCT")


# %%
# Save all negative control gene weights to 
Dict2Fil(HDL_C_top_gw, "{}/HDL_C.top61.gw.csv".format(GeneWeightDIR))
Dict2Fil(T2D_top_gw, "{}/T2D.top61.gw.csv".format(GeneWeightDIR))
Dict2Fil(IBD_top_gw, "{}/IBD.top61.gw.csv".format(GeneWeightDIR))
Dict2Fil(hba1c_top_gw, "{}/hba1c.top61.gw.csv".format(GeneWeightDIR))
Dict2Fil(Parkinson_top_gw, "{}/Parkinson.top61.gw.csv".format(GeneWeightDIR))
Dict2Fil(ALZ_top_gw, "{}/ALZ.top60.gw.csv".format(GeneWeightDIR))
Dict2Fil(Alanine_top_gw, "{}/Alanine.top61.gw.csv".format(GeneWeightDIR))
print("All negative control gene weights saved to GeneWeightDIR")


# %%

# %%
# for g in hba1c_top_gw.keys():
#     print(Entrez2Symbol[g])

# %%

# %%
# Negative control biases
if 'HDL_C_top_bias' in globals():
    HDL_C_top_bias.to_csv("{}/HCT.HDL_C.top61.csv".format(Bias_Save_Dir))
    print("Saved: HDL_C_top_bias")
if 'T2D_top_bias' in globals():
    T2D_top_bias.to_csv("{}/HCT.T2D.top61.csv".format(Bias_Save_Dir))
    print("Saved: T2D_top_bias")
if 'IBD_top_bias' in globals():
    IBD_top_bias.to_csv("{}/HCT.IBD.top61.csv".format(Bias_Save_Dir))
    print("Saved: IBD_top_bias")
if 'hba1c_top_bias' in globals():
    hba1c_top_bias.to_csv("{}/HCT.hba1c.top61.csv".format(Bias_Save_Dir))
    print("Saved: hba1c_top_bias")
if 'Parkinson_top_bias' in globals():
    Parkinson_top_bias.to_csv("{}/HCT.Parkinson.top61.csv".format(Bias_Save_Dir))
    print("Saved: Parkinson_top_bias")
if 'ALZ_top_bias' in globals():
    ALZ_top_bias.to_csv("{}/HCT.ALZ.top60.csv".format(Bias_Save_Dir))
    print("Saved: ALZ_top_bias")
if 'Alanine_top_bias' in globals():
    Alanine_top_bias.to_csv("{}/HCT.Alanine.top61.csv".format(Bias_Save_Dir))
    print("Saved: Alanine_top_bias")
print(f"\nAll bias files saved to: {Bias_Save_Dir}")


# %%

# %%
Bias_Save_Dir2 = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/results/main_results/random/Centering/"

VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir2 + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir2 + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir2 + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir2 + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)


# %%
VNR_Contrast = compare_biases(VNR_Neg_Bias, VNR_Pos_Bias, name1="VNR-", name2="VNR+", neurons=ALL_CTs)
VNR_Contrast_Neurons = VNR_Contrast[VNR_Contrast.index.isin(Neurons)]
VNR_Contrast_Neurons = VNR_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_PKS_Contrast = compare_biases(VNR_Neg_Bias, Parkinson_top_bias, name1="VNR-", name2="PKS+", neurons=ALL_CTs)
VNR_PKS_Contrast_Neurons = VNR_PKS_Contrast[VNR_PKS_Contrast.index.isin(Neurons)]
VNR_PKS_Contrast_Neurons = VNR_PKS_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_IBD_Contrast = compare_biases(VNR_Neg_Bias, IBD_top_bias, name1="VNR-", name2="IBD+", neurons=ALL_CTs)
VNR_IBD_Contrast_Neurons = VNR_IBD_Contrast[VNR_IBD_Contrast.index.isin(Neurons)]
VNR_IBD_Contrast_Neurons = VNR_IBD_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_T2D_Contrast = compare_biases(VNR_Neg_Bias, T2D_top_bias, name1="VNR-", name2="T2D+", neurons=ALL_CTs)
VNR_T2D_Contrast_Neurons = VNR_T2D_Contrast[VNR_T2D_Contrast.index.isin(Neurons)]
VNR_T2D_Contrast_Neurons = VNR_T2D_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_HDL_C_Contrast = compare_biases(VNR_Neg_Bias, HDL_C_top_bias, name1="VNR-", name2="HDL_C+", neurons=ALL_CTs)
VNR_HDL_C_Contrast_Neurons = VNR_HDL_C_Contrast[VNR_HDL_C_Contrast.index.isin(Neurons)]
VNR_HDL_C_Contrast_Neurons = VNR_HDL_C_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_alanine_Contrast = compare_biases(VNR_Neg_Bias, Alanine_top_bias, name1="VNR-", name2="alanine+", neurons=ALL_CTs)
VNR_alanine_Contrast_Neurons = VNR_alanine_Contrast[VNR_alanine_Contrast.index.isin(Neurons)]
VNR_alanine_Contrast_Neurons = VNR_alanine_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)

VNR_ALZ_Contrast = compare_biases(VNR_Neg_Bias, ALZ_top_bias, name1="VNR-", name2="ALZ+", neurons=ALL_CTs)
VNR_ALZ_Contrast_Neurons = VNR_ALZ_Contrast[VNR_ALZ_Contrast.index.isin(Neurons)]
VNR_ALZ_Contrast_Neurons = VNR_ALZ_Contrast_Neurons.sort_values(by="Bias_Diff", ascending=False)



# %%
# Compare VNR_Neg with each control gene set for CGE interneuron
ctrl_biases = []
ctrl_contrasts = []
ctrl_names = []

# Build lists only for variables that exist
if 'Parkinson_top_bias' in globals() and 'VNR_PKS_Contrast' in globals():
    ctrl_biases.append(Parkinson_top_bias)
    ctrl_contrasts.append(VNR_PKS_Contrast)
    ctrl_names.append("Parkinson")

if 'IBD_top_bias' in globals() and 'VNR_IBD_Contrast' in globals():
    ctrl_biases.append(IBD_top_bias)
    ctrl_contrasts.append(VNR_IBD_Contrast)
    ctrl_names.append("IBD")

if 'T2D_top_bias' in globals() and 'VNR_T2D_Contrast' in globals():
    ctrl_biases.append(T2D_top_bias)
    ctrl_contrasts.append(VNR_T2D_Contrast)
    ctrl_names.append("T2D")

if 'HDL_C_top_bias' in globals() and 'VNR_HDL_C_Contrast' in globals():
    ctrl_biases.append(HDL_C_top_bias)
    ctrl_contrasts.append(VNR_HDL_C_Contrast)
    ctrl_names.append("HDL_C")

if 'Alanine_top_bias' in globals() and 'VNR_alanine_Contrast_Neurons' in globals():
    ctrl_biases.append(Alanine_top_bias)
    ctrl_contrasts.append(VNR_alanine_Contrast)
    ctrl_names.append("alanine")

if 'ALZ_top_bias' in globals() and 'VNR_ALZ_Contrast' in globals():
    ctrl_biases.append(ALZ_top_bias)
    ctrl_contrasts.append(VNR_ALZ_Contrast)
    ctrl_names.append("ALZ")

# Run comparisons if VNR_Neg_Bias exists and we have at least one control
if 'VNR_Neg_Bias' in globals() and len(ctrl_biases) > 0:
    for ctrl_bias, ctrl_contrast, ctrl_name in zip(ctrl_biases, ctrl_contrasts, ctrl_names):
        CompareSingleCT(ctrl_bias, VNR_Neg_Bias, "CGE interneuron", ctrl_contrast, 
                       f"{ctrl_name} Control Bias", "VNR - Mutation Bias", loc=(0.0, 0))
else:
    print("Warning: VNR_Neg_Bias or control biases not found. Please run previous cells first.")
    
    

# %%

# %%

# %%
