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
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

try:
    os.chdir(f"{ProjDIR}/notebooks/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:
    print(f"Unexpected error: {e}")

#print(f"Current working directory: {os.getcwd()}")

# %% [markdown]
# ### Try 22q DEGs from https://www.nature.com/articles/s41467-022-31436-8#Abs1

# %%
DEG_22q_day0 = pd.read_excel("../dat/suppl.data/41467_2022_31436_MOESM4_ESM.xlsx", skiprows=1)
DEG_22q_day0 = DEG_22q_day0[DEG_22q_day0["gene_biotype"] == "protein_coding"]
DEG_22q_day0["EntrezID"] = DEG_22q_day0["external_gene_name"].map(GeneSymbol2Entrez)
DEG_22q_day0 = DEG_22q_day0[DEG_22q_day0["EntrezID"].notnull()].copy()
DEG_22q_day0["EntrezID"] = DEG_22q_day0["EntrezID"].astype(int)
DEG_22q_day0.set_index("EntrezID", inplace=True)
DEG_22q_day0 = DEG_22q_day0.sort_values(by="padj", ascending=True)

# %%
DEG_22q_day4 = pd.read_excel("../dat/suppl.data/41467_2022_31436_MOESM5_ESM.xlsx", skiprows=1)
DEG_22q_day4 = DEG_22q_day4[DEG_22q_day4["gene_biotype"] == "protein_coding"]
DEG_22q_day4["EntrezID"] = DEG_22q_day4["external_gene_name"].map(GeneSymbol2Entrez)
DEG_22q_day4 = DEG_22q_day4[DEG_22q_day4["EntrezID"].notnull()].copy()
DEG_22q_day4["EntrezID"] = DEG_22q_day4["EntrezID"].astype(int)
DEG_22q_day4.set_index("EntrezID", inplace=True)
DEG_22q_day4 = DEG_22q_day4.sort_values(by="padj", ascending=True)

# %%
DEG_22q_day28 = pd.read_excel("../dat/suppl.data/41467_2022_31436_MOESM6_ESM.xlsx", skiprows=1)
DEG_22q_day28 = DEG_22q_day28[DEG_22q_day28["gene_biotype"] == "protein_coding"]
DEG_22q_day28["EntrezID"] = DEG_22q_day28["external_gene_name"].map(GeneSymbol2Entrez)
DEG_22q_day28 = DEG_22q_day28[DEG_22q_day28["EntrezID"].notnull()].copy()
DEG_22q_day28["EntrezID"] = DEG_22q_day28["EntrezID"].astype(int)
DEG_22q_day28.set_index("EntrezID", inplace=True)
DEG_22q_day28 = DEG_22q_day28.sort_values(by="padj", ascending=True)

# %%
DEG_22q_day0.head(5)

# %%
DEG_22q_day4.head(5)

# %%
DEG_22q_day28.head(5)

# %%
# these are day0
FDR_cut = 0.05
DEG_22q_day0_sig = DEG_22q_day0[DEG_22q_day0["padj"] < FDR_cut]
DEG_22q_day0_sig_up = DEG_22q_day0_sig[DEG_22q_day0_sig["log2FoldChange"] > 0]
DEG_22q_day0_sig_down = DEG_22q_day0_sig[DEG_22q_day0_sig["log2FoldChange"] < 0]
print(f"Day 0: {len(DEG_22q_day0_sig)} significant genes (padj<{FDR_cut})")
print(f"    Upregulated: {len(DEG_22q_day0_sig_up)}")
print(f"    Downregulated: {len(DEG_22q_day0_sig_down)}\n")

# these are day4
DEG_22q_day4_sig = DEG_22q_day4[DEG_22q_day4["padj"] < FDR_cut]
DEG_22q_day4_sig_up = DEG_22q_day4_sig[DEG_22q_day4_sig["log2FoldChange"] > 0]
DEG_22q_day4_sig_down = DEG_22q_day4_sig[DEG_22q_day4_sig["log2FoldChange"] < 0]
print(f"Day 4: {len(DEG_22q_day4_sig)} significant genes (padj<{FDR_cut})")
print(f"    Upregulated: {len(DEG_22q_day4_sig_up)}")
print(f"    Downregulated: {len(DEG_22q_day4_sig_down)}\n")

# these are day28
DEG_22q_day28_sig = DEG_22q_day28[DEG_22q_day28["padj"] < FDR_cut]
DEG_22q_day28_sig_up = DEG_22q_day28_sig[DEG_22q_day28_sig["log2FoldChange"] > 0]
DEG_22q_day28_sig_down = DEG_22q_day28_sig[DEG_22q_day28_sig["log2FoldChange"] < 0]
print(f"Day 28: {len(DEG_22q_day28_sig)} significant genes (padj<  {FDR_cut})")
print(f"    Upregulated: {len(DEG_22q_day28_sig_up)}")
print(f"    Downregulated: {len(DEG_22q_day28_sig_down)}")


# %%
X22q_GW = Fil2Dict("../dat/GeneWeights/X22q.gw.csv")
X22q_Genes = list(X22q_GW.keys())

# %%
import matplotlib.pyplot as plt
import numpy as np

# Prepare lists for the DESeq2 results and titles
deg_list = [
    (DEG_22q_day0, "Day 0"),
    (DEG_22q_day4, "Day 4"),
    (DEG_22q_day28, "Day 28"),
]

fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
for idx, (ax, (deg_df, label)) in enumerate(zip(axes, deg_list)):
    x = deg_df['log2FoldChange']
    y = deg_df['padj']

    # Replace padj=0 with a small value to avoid log(0)
    y_log = -1 * np.log10(np.where(y > 0, y, 1e-300))

    # Identify 22q genes for red highlight
    is_22q = deg_df.index.isin(X22q_Genes)

    # Plot all genes in grey
    ax.scatter(x[~is_22q], y_log[~is_22q], c='grey', alpha=0.7, edgecolor='none', s=15, label='Other genes')
    # Plot 22q genes in red
    ax.scatter(x[is_22q], y_log[is_22q], c='red', alpha=0.8, edgecolor='none', s=18, label='22q genes')

    ax.axvline(0, color='black', lw=1, ls='--')
    ax.set_xlabel('log2(Fold Change)')
    ax.set_title(label)
    if idx == 0:
        ax.set_ylabel('-log10(padj)')
    if idx == len(axes) - 1:
        # Place legend on the last subplot only
        ax.legend()

fig.suptitle('Volcano Plots for 22q DEGs (Day 0, Day 4, Day 28)', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np

# Prepare lists for the DESeq2 results and titles, excluding 22q genes in each case
deg_list = [
    (DEG_22q_day0, "Day 0"),
    (DEG_22q_day4, "Day 4"),
    (DEG_22q_day28, "Day 28"),
]

fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
for idx, (ax, (deg_df, label)) in enumerate(zip(axes, deg_list)):
    # Exclude 22q genes
    not_22q = ~deg_df.index.isin(X22q_Genes)
    x = deg_df.loc[not_22q, 'log2FoldChange']
    y = deg_df.loc[not_22q, 'padj']

    # Replace padj=0 with a small value to avoid log(0)
    y_log = -1 * np.log10(np.where(y > 0, y, 1e-300))

    # Plot only non-22q genes
    ax.scatter(x, y_log, c='grey', alpha=0.7, edgecolor='none', s=15, label='Other genes (non-22q)')

    ax.axvline(0, color='black', lw=1, ls='--')
    ax.set_xlabel('log2(Fold Change)')
    ax.set_title(label)
    if idx == 0:
        ax.set_ylabel('-log10(padj)')
    if idx == len(axes) - 1:
        # Place legend on the last subplot only
        ax.legend()

fig.suptitle('Volcano Plots for 22q DEGs EXCLUDING 22q genes (Day 0, Day 4, Day 28)', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# %%
# Create DEG lists for Day 0, Day 4, Day 28, each EXCLUDING 22q genes
DEG_22q_sig_exl22q_day0 = DEG_22q_day0_sig[~DEG_22q_day0_sig.index.isin(X22q_Genes)]
DEG_22q_sig_up_exl22q_day0 = DEG_22q_day0_sig_up[~DEG_22q_day0_sig_up.index.isin(X22q_Genes)]
DEG_22q_sig_down_exl22q_day0 = DEG_22q_day0_sig_down[~DEG_22q_day0_sig_down.index.isin(X22q_Genes)]

DEG_22q_sig_exl22q_day4 = DEG_22q_day4_sig[~DEG_22q_day4_sig.index.isin(X22q_Genes)]
DEG_22q_sig_up_exl22q_day4 = DEG_22q_day4_sig_up[~DEG_22q_day4_sig_up.index.isin(X22q_Genes)]
DEG_22q_sig_down_exl22q_day4 = DEG_22q_day4_sig_down[~DEG_22q_day4_sig_down.index.isin(X22q_Genes)]

DEG_22q_sig_exl22q_day28 = DEG_22q_day28_sig[~DEG_22q_day28_sig.index.isin(X22q_Genes)]
DEG_22q_sig_up_exl22q_day28 = DEG_22q_day28_sig_up[~DEG_22q_day28_sig_up.index.isin(X22q_Genes)]
DEG_22q_sig_down_exl22q_day28 = DEG_22q_day28_sig_down[~DEG_22q_day28_sig_down.index.isin(X22q_Genes)]

(len(DEG_22q_sig_exl22q_day0), len(DEG_22q_sig_up_exl22q_day0), len(DEG_22q_sig_down_exl22q_day0)), \
(len(DEG_22q_sig_exl22q_day4), len(DEG_22q_sig_up_exl22q_day4), len(DEG_22q_sig_down_exl22q_day4)), \
(len(DEG_22q_sig_exl22q_day28), len(DEG_22q_sig_up_exl22q_day28), len(DEG_22q_sig_down_exl22q_day28))


# %%
HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv", index_col=0)
#HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv", index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
HumanCT_Z2_HCT.shape

# %%
#DF = DEG_22q_sig_exl22q_day0
DF = DEG_22q_sig_down_exl22q_day0
#DEG_day0_X22q_GW = dict(zip(DF.index, abs(DF["log2FoldChange"])))
DEG_day0_X22q_GW = dict(zip(DF.index, np.ones(len(DF))))

DF = DEG_22q_sig_down_exl22q_day4
#DEG_day4_X22q_GW = dict(zip(DF.index, abs(DF["log2FoldChange"])))
DEG_day4_X22q_GW = dict(zip(DF.index, np.ones(len(DF))))

DF = DEG_22q_sig_down_exl22q_day28
#DEG_day28_X22q_GW = dict(zip(DF.index, abs(DF["log2FoldChange"])))
DEG_day28_X22q_GW = dict(zip(DF.index, np.ones(len(DF))))
#DEG_X22q_GW = dict(zip(DEG_22q_sig_down.index, DEG_22q_sig_down["log2FoldChange"]))
#X22q_DEG_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_X22q_GW)

# %%
DEG_day0_X22q_GW

# %%
X22q_DEG_day0_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_day0_X22q_GW, verbose=True)
X22q_DEG_day0_Bias = AnnotateCTDat(X22q_DEG_day0_Bias, Anno)

X22q_DEG_day4_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_day4_X22q_GW, verbose=True)
X22q_DEG_day4_Bias = AnnotateCTDat(X22q_DEG_day4_Bias, Anno)

X22q_DEG_day28_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_day28_X22q_GW, verbose=True)
X22q_DEG_day28_Bias = AnnotateCTDat(X22q_DEG_day28_Bias, Anno)
#X22q_DEG_Bias.to_csv("{}/HCT.X22q.csv".format(Bias_Save_Dir))
#fdr_cut = 0.1

# Create a 2x3 subplot grid
fig, axes = plt.subplots(1, 3, figsize=(22, 10), dpi=120, facecolor='none')
fig.patch.set_alpha(0.0)
plt.style.use('seaborn-v0_8-whitegrid')

# Flatten axes array for easier indexing
axes = axes.flatten()

# Plot each dataset on its own subplot
SuperClusterBias_BoxPlot(X22q_DEG_day0_Bias, "22q day0 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[0])
SuperClusterBias_BoxPlot(X22q_DEG_day4_Bias, "22q day4 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[1])
SuperClusterBias_BoxPlot(X22q_DEG_day28_Bias, "22q day28 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[2])

plt.tight_layout()
plt.show()


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

# %%
common_indices = X22q_DEG_day0_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day0_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

common_indices = X22q_DEG_day4_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day4_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

common_indices = X22q_DEG_day28_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day28_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

# %% [markdown]
# ### Try 22q DEGs from https://www.nature.com/articles/s41591-020-1043-9#Abs1

# %%
DEG_22q_day25 = pd.read_excel("../dat/suppl.data/41591_2020_1043_MOESM3_ESM.xlsx", sheet_name="Day25")
#DEG_22q_day25 = DEG_22q_day25[DEG_22q_day25["gene_biotype"] == "protein_coding"]
DEG_22q_day25["EntrezID"] = DEG_22q_day25["Gene"].map(GeneSymbol2Entrez)
DEG_22q_day25 = DEG_22q_day25[DEG_22q_day25["EntrezID"].notnull()].copy()
DEG_22q_day25["EntrezID"] = DEG_22q_day25["EntrezID"].astype(int)
DEG_22q_day25.set_index("EntrezID", inplace=True)
DEG_22q_day25 = DEG_22q_day25.sort_values(by="FDR", ascending=True)

# %%
DEG_22q_day50 = pd.read_excel("../dat/suppl.data/41591_2020_1043_MOESM3_ESM.xlsx", sheet_name="Day50")
#DEG_22q_day50 = DEG_22q_day50[DEG_22q_day50["gene_biotype"] == "protein_coding"]
DEG_22q_day50["EntrezID"] = DEG_22q_day50["Gene"].map(GeneSymbol2Entrez)
DEG_22q_day50 = DEG_22q_day50[DEG_22q_day50["EntrezID"].notnull()].copy()
DEG_22q_day50["EntrezID"] = DEG_22q_day50["EntrezID"].astype(int)
DEG_22q_day50.set_index("EntrezID", inplace=True)
DEG_22q_day50 = DEG_22q_day50.sort_values(by="FDR", ascending=True)

# %%
DEG_22q_day75 = pd.read_excel("../dat/suppl.data/41591_2020_1043_MOESM3_ESM.xlsx", sheet_name="Day75")
#DEG_22q_day75 = DEG_22q_day75[DEG_22q_day75["gene_biotype"] == "protein_coding"]
DEG_22q_day75["EntrezID"] = DEG_22q_day75["Gene"].map(GeneSymbol2Entrez)
DEG_22q_day75 = DEG_22q_day75[DEG_22q_day75["EntrezID"].notnull()].copy()
DEG_22q_day75["EntrezID"] = DEG_22q_day75["EntrezID"].astype(int)
DEG_22q_day75.set_index("EntrezID", inplace=True)
DEG_22q_day75 = DEG_22q_day75.sort_values(by="FDR", ascending=True)

# %%
DEG_22q_day100 = pd.read_excel("../dat/suppl.data/41591_2020_1043_MOESM3_ESM.xlsx", sheet_name="Day100")
#DEG_22q_day100 = DEG_22q_day100[DEG_22q_day100["gene_biotype"] == "protein_coding"]
DEG_22q_day100["EntrezID"] = DEG_22q_day100["Gene"].map(GeneSymbol2Entrez)
DEG_22q_day100 = DEG_22q_day100[DEG_22q_day100["EntrezID"].notnull()].copy()
DEG_22q_day100["EntrezID"] = DEG_22q_day100["EntrezID"].astype(int)
DEG_22q_day100.set_index("EntrezID", inplace=True)
DEG_22q_day100 = DEG_22q_day100.sort_values(by="FDR", ascending=True)


# %%
def GeneFilt(DF):
    DF = DF[DF["FDR"] < 0.1]
    DF = DF[DF["beta"] < 0]
    #DF = DF[~DF.index.isin(X22q_Genes)]
    return DF



# %%
DEG_22q_day25_filt = GeneFilt(DEG_22q_day25)
DEG_22q_day50_filt = GeneFilt(DEG_22q_day50)
DEG_22q_day75_filt = GeneFilt(DEG_22q_day75)
DEG_22q_day100_filt = GeneFilt(DEG_22q_day100)

# %%
# Get gene weights (GW) for 22q DEGs at days 25, 50, 75, and 100
DF_25 = DEG_22q_day25_filt
DEG_22q_day25_GW = dict(zip(DF_25.index, np.ones(len(DF_25))))
#DEG_22q_day25_GW = dict(zip(DF_25.index, abs(DF_25["beta"])))
print(f"Number of genes in DEG_22q_day25_GW: {len(DEG_22q_day25_GW)}")

DF_50 = DEG_22q_day50_filt
DEG_22q_day50_GW = dict(zip(DF_50.index, np.ones(len(DF_50))))
#DEG_22q_day50_GW = dict(zip(DF_50.index, abs(DF_50["beta"])))
print(f"Number of genes in DEG_22q_day50_GW: {len(DEG_22q_day50_GW)}")

DF_75 = DEG_22q_day75_filt
DEG_22q_day75_GW = dict(zip(DF_75.index, np.ones(len(DF_75))))
#DEG_22q_day75_GW = dict(zip(DF_75.index, abs(DF_75["beta"])))
print(f"Number of genes in DEG_22q_day75_GW: {len(DEG_22q_day75_GW)}")

DF_100 = DEG_22q_day100_filt
DEG_22q_day100_GW = dict(zip(DF_100.index, np.ones(len(DF_100))))
#DEG_22q_day100_GW = dict(zip(DF_100.index, abs(DF_100["beta"])))

print(f"Number of genes in DEG_22q_day100_GW: {len(DEG_22q_day100_GW)}")

# %%
# Compute HumanCT bias for 22q DEGs at days 25, 50, 75, and 100
X22q_DEG_day25_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_22q_day25_GW, verbose=True)
X22q_DEG_day25_Bias = AnnotateCTDat(X22q_DEG_day25_Bias, Anno)

try:
    X22q_DEG_day50_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_22q_day50_GW, verbose=True)
    X22q_DEG_day50_Bias = AnnotateCTDat(X22q_DEG_day50_Bias, Anno)
except:
    X22q_DEG_day50_Bias = None

X22q_DEG_day75_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_22q_day75_GW, verbose=True)
X22q_DEG_day75_Bias = AnnotateCTDat(X22q_DEG_day75_Bias, Anno)

try:
    X22q_DEG_day100_Bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, DEG_22q_day100_GW, verbose=True)
    X22q_DEG_day100_Bias = AnnotateCTDat(X22q_DEG_day100_Bias, Anno)
except:
    X22q_DEG_day100_Bias = None

Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul24_Centered/"
X22q_DEG_day25_Bias.to_csv("{}/DEG_22q_day25_Bias.csv".format(Bias_Save_Dir))
X22q_DEG_day75_Bias.to_csv("{}/DEG_22q_day75_Bias.csv".format(Bias_Save_Dir))
X22q_DEG_day50_Bias.to_csv("{}/DEG_22q_day50_Bias.csv".format(Bias_Save_Dir))
X22q_DEG_day100_Bias.to_csv("{}/DEG_22q_day100_Bias.csv".format(Bias_Save_Dir))

#Create a 1x4 subplot grid for easier visual comparison
fig, axes = plt.subplots(1, 4, figsize=(28, 10), dpi=120, facecolor='none')
fig.patch.set_alpha(0.0)
plt.style.use('seaborn-v0_8-whitegrid')
axes = axes.flatten()
SuperClusterBias_BoxPlot(X22q_DEG_day25_Bias, "22q day25 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[0])
SuperClusterBias_BoxPlot(X22q_DEG_day50_Bias, "22q day50 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[1])
SuperClusterBias_BoxPlot(X22q_DEG_day75_Bias, "22q day75 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[2])
SuperClusterBias_BoxPlot(X22q_DEG_day100_Bias, "22q day100 DEGs", NeuroOnly=False, sortby="mean", EffectCol="EFFECT", ax=axes[3])
plt.tight_layout()
plt.show()


# %%
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

plt.figure(figsize=(6,6))
venn2([set(DEG_22q_day25_GW.keys()), set(DEG_22q_day75_GW.keys())], set_labels=("Day 25 DEGs", "Day 75 DEGs"))
plt.title("Venn Diagram of DEG Keys: Day 25 vs Day 75")
plt.show()

# %%
SuperClusterBias_BoxPlot(X22q_DEG_day25_Bias, "22q day25 DEGs", NeuroOnly=False, sortby="mean")
SuperClusterBias_BoxPlot(X22q_DEG_day75_Bias, "22q day75 DEGs", NeuroOnly=False, sortby="mean")

# %%
common_indices = X22q_DEG_day25_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day25_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

common_indices = X22q_DEG_day50_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day50_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

common_indices = X22q_DEG_day75_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day75_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

common_indices = X22q_DEG_day100_Bias.index.intersection(VIP_Anno.index)
X22q_DEG_Test = X22q_DEG_day100_Bias.loc[common_indices].copy()
X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
plot_vip_effect_comparison(X22q_DEG_Test)

# %%
# common_indices = X22q_DEG_day25_Bias.index.intersection(VIP_Anno.index)
# X22q_DEG_Test = X22q_DEG_day25_Bias.loc[common_indices].copy()
# X22q_DEG_Test["VIP"] = VIP_Anno.loc[common_indices, "VIP"]
# X22q_DEG_Test

# %%
DEG_22q_day25

# %%
Test_28 = DEG_22q_day28[DEG_22q_day28["padj"] < 0.05]["external_gene_name"]
Test_28

# %%
Test_25 = DEG_22q_day25[DEG_22q_day25["FDR"] < 0.05]["Gene"]
Test_25

# %%
X22q_Genes_symbol = [Entrez2Symbol.get(x, None) for x in X22q_Genes]
X22q_Genes_symbol

# %%
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Convert the Series to sets for Venn diagram
set_25 = set(Test_25)
set_28 = set(Test_28)

# Compute the intersection
common_genes = set_25 & set_28

# Find out how many of these are 22q11.2 genes
# Assume you have a set called genes_22q11 that contains all 22q11.2 gene symbols
try:
    count_22q_in_common = len(common_genes & set(X22q_Genes_symbol))
except NameError:
    print("Define `genes_22q11` as a set of 22q11.2 gene symbols to count overlap.")
    count_22q_in_common = None

plt.figure(figsize=(6, 6))
venn2([set_25, set_28], set_labels=("Day 25 DEGs", "Day 28 DEGs"))
plt.title(f"Venn diagram of Day 25 and Day 28 DEGs\nCommon in 22q11.2: {count_22q_in_common}")
plt.show()

if count_22q_in_common is not None:
    print(f"Number of overlapping genes that are 22q11.2: {count_22q_in_common}")

# %%

# %%
