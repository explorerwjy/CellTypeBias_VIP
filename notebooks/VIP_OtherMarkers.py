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


import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection, multipletests

# %%
#HumanCT_res_df_GeneL = pd.read_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.cluster.June10.csv", index_col=0)
#HumanCT_res_df_GeneL.head(5)

# %%
ExpL = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv", index_col=0)
ExpL.columns = [int(x) for x in ExpL.columns.values]
ExpL = ExpL.applymap(lambda x: np.log2(x + 1))

# %%
TPM = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.Filt.csv", index_col=0)
TPM = np.log1p(TPM)
TPM.columns = [int(x) for x in TPM.columns.values ]

# %%
CGE_idx = Anno[Anno["Supercluster"]=="CGE interneuron"].index.values
MGE_idx = Anno[Anno["Supercluster"]=="MGE interneuron"].index.values

# Define markers to plot
markers = {
    'SST': GeneSymbol2Entrez["SST"],
    'PVALB': GeneSymbol2Entrez["PVALB"],
    "CCK": GeneSymbol2Entrez["CCK"],
    "CALB2": GeneSymbol2Entrez["CALB2"],
    "SNCG": GeneSymbol2Entrez["SNCG"],
    "VIP": GeneSymbol2Entrez["VIP"],
    "M2R": GeneSymbol2Entrez["CHRM2"],
    "CNR1": GeneSymbol2Entrez["CNR1"]
}

cutoffs = {
    "VIP": 1,
    "CCK": 5,
    "CR": 1.5,
    "SNCG": 0.5,
    "M2R": 2,
    "CNR1": 1
}

# Calculate number of rows needed (2 markers per row)
n_markers = len(markers)
n_cols = 2
n_rows = math.ceil(n_markers / n_cols)

# Create subplot layout
fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows), dpi=120)
fig.suptitle("CGE vs MGE Expression Comparison", fontsize=14, y=1.02)

# Flatten axes array for easier indexing
if n_rows == 1:
    axes = axes.reshape(1, -1)
axes = axes.flatten()

for idx, (marker_name, entrez_id) in enumerate(markers.items()):
    cge_expr = ExpL.loc[entrez_id, CGE_idx]
    mge_expr = ExpL.loc[entrez_id, MGE_idx]
    ax = axes[idx]
    sns.histplot(cge_expr, bins=20, ax=ax, color='blue', alpha=0.7, label='CGE')
    sns.histplot(mge_expr, bins=20, ax=ax, color='red', alpha=0.7, label='MGE')
    ax.set_title(f"{marker_name}")
    ax.set_xlabel("Log2(UMI + 1)")
    ax.set_ylabel("Frequency")
    ax.legend()

    # Add vertical line for cutoff if available
    cutoff_val = None
    # We match marker_name to key in cutoffs (marker_name or e.g. 'CR' for 'CALB2')
    cutoff_key = marker_name
    if marker_name == "CALB2":
        cutoff_key = "CR"  # CALB2 is the gene, CR is the marker label for cutoff
    if cutoff_key in cutoffs:
        cutoff_val = cutoffs[cutoff_key]
    if cutoff_val is not None:
        ax.axvline(cutoff_val, color="k", linestyle="--", linewidth=1.5, label="Cutoff")
        # To avoid duplicate legends, only add cutoff to one handle per ax
        handles, labels = ax.get_legend_handles_labels()
        # Remove duplicate "Cutoff" if already present
        if labels.count("Cutoff") > 1:
            cutoff_indices = [i for i, label in enumerate(labels) if label == "Cutoff"]
            for i in cutoff_indices[1:]:
                del handles[i]
                del labels[i]
        ax.legend(handles, labels)

# Hide any unused subplots
for idx in range(n_markers, len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.show()


# %%
#Bias_Save_Dict = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jun09/"
#Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul07/"
Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/results/random/Centering/"

SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
X22q_mouse_Bias = pd.read_csv(Bias_Save_Dir + "22q_small_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)

# VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
# VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_Bias_addP.csv", index_col=0)
# EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_Bias_addP.csv", index_col=0)
# EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_Bias_addP.csv", index_col=0)

#VNR_Pos_Bias = pd.read_csv(os.path.join(Bias_Save_Dict, "CTBias.VNR.Pos.top61.csv"), index_col=0)
#VNR_Neg_Bias = pd.read_csv(os.path.join(Bias_Save_Dict, "CTBias.VNR.Neg.top61.csv"), index_col=0)
Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul24_Centered/"
X22q_Hexp_Bias = pd.read_csv("{}/HCT.22q.HighExp.csv".format(Bias_Save_Dir), index_col=0)
X22q_DEG_day25_Bias = pd.read_csv("{}/DEG_22q_day25_Bias.csv".format(Bias_Save_Dir), index_col=0)
X22q_DEG_day50_Bias = pd.read_csv("{}/DEG_22q_day50_Bias.csv".format(Bias_Save_Dir), index_col=0)
X22q_DEG_day75_Bias = pd.read_csv("{}/DEG_22q_day75_Bias.csv".format(Bias_Save_Dir), index_col=0)
X22q_DEG_day100_Bias = pd.read_csv("{}/DEG_22q_day100_Bias.csv".format(Bias_Save_Dir), index_col=0)

# %%

# %%
#BiasDFs = [X22q_Bias, HighIQ_ASD_Bias, LowIQ_ASD_Bias, VNR_Pos_Bias, VNR_Neg_Bias, DDD_Bias, SCZ_Bias]
BiasDFs = [X22q_Bias, X22q_mouse_Bias,  X22q_Hexp_Bias, X22q_DEG_day25_Bias, X22q_DEG_day50_Bias, X22q_DEG_day75_Bias, X22q_DEG_day100_Bias]
BiasDF_Names = ["X22q", "X22q_mouse", "X22q_Hexp", "X22q_DEG_day25", "X22q_DEG_day50", "X22q_DEG_day75", "X22q_DEG_day100"]
for i, row in X22q_Bias.iterrows():
    for BiasDF, name in zip(BiasDFs, BiasDF_Names):
        X22q_Bias.loc[i, name+"_Bias"] = BiasDF.loc[i, "EFFECT"]

# %%
#TmpDF = X22q_Bias[["Class","Supercluster","X22q_Bias", "X22q_mouse_Bias", "LowIQ_ASD_Bias", "VNR_Pos_Bias", "VNR_Neg_Bias", "DDD_Bias", "SCZ_Bias"]]
TmpDF = X22q_Bias[["Class","Supercluster","X22q_Bias", "X22q_mouse_Bias", "X22q_Hexp_Bias", "X22q_DEG_day25_Bias", "X22q_DEG_day50_Bias", "X22q_DEG_day75_Bias", "X22q_DEG_day100_Bias"]]

# %%
# More efficient, avoids chained assignment warnings and uses .loc for columns in bulk

Biases_CGE_ExpL = TmpDF[TmpDF["Supercluster"] == "CGE interneuron"].copy()

# List marker gene names and their column suffixes
marker_genes = {
    "VIP_ExpL": "VIP",
    "CCK_ExpL": "CCK",
    "CR_ExpL": "CALB2",
    "SNCG_ExpL": "SNCG",
    "M2R_ExpL": "CHRM2",
    "CNR1_ExpL": "CNR1"
}

for colname, symbol in marker_genes.items():
    entrez_id = GeneSymbol2Entrez[symbol]
    # Extract matching columns in proper order
    Biases_CGE_ExpL.loc[:, colname] = ExpL.loc[entrez_id, Biases_CGE_ExpL.index].values


# %%
Biases_CGE_ExpL.sort_values(by="X22q_Bias", ascending=False, inplace=True)
#Biases_CGE_ExpL

# %%
CGE_idx = Anno[Anno["Supercluster"]=="CGE interneuron"].index.values
MGE_idx = Anno[Anno["Supercluster"]=="MGE interneuron"].index.values

# Define markers to plot
markers = {
    'SST': GeneSymbol2Entrez["SST"],
    'PVALB': GeneSymbol2Entrez["PVALB"],
    "CCK": GeneSymbol2Entrez["CCK"],
    "CALB2": GeneSymbol2Entrez["CALB2"],
    "SNCG": GeneSymbol2Entrez["SNCG"],
    "VIP": GeneSymbol2Entrez["VIP"],
    "M2R": GeneSymbol2Entrez["CHRM2"],
    "CNR1": GeneSymbol2Entrez["CNR1"]
}

cutoffs = {
    "VIP": 1,
    "CCK": 2,
    "CR": 1.5,
    "SNCG": 0.5,
    "M2R": 2,
    "CNR1": 1
}

# Calculate number of rows needed (2 markers per row)
n_markers = len(markers)
n_cols = 2
n_rows = math.ceil(n_markers / n_cols)

# Create subplot layout
fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows), dpi=120)
fig.suptitle("CGE vs MGE Expression Comparison", fontsize=14, y=1.02)

# Flatten axes array for easier indexing
if n_rows == 1:
    axes = axes.reshape(1, -1)
axes = axes.flatten()

for idx, (marker_name, entrez_id) in enumerate(markers.items()):
    cge_expr = ExpL.loc[entrez_id, CGE_idx]
    mge_expr = ExpL.loc[entrez_id, MGE_idx]
    ax = axes[idx]
    sns.histplot(cge_expr, bins=20, ax=ax, color='blue', alpha=0.7, label='CGE')
    sns.histplot(mge_expr, bins=20, ax=ax, color='red', alpha=0.7, label='MGE')
    ax.set_title(f"{marker_name}")
    ax.set_xlabel("Log2(UMI + 1)")
    ax.set_ylabel("Frequency")
    ax.legend()

    # Add vertical line for cutoff if available
    cutoff_val = None
    # We match marker_name to key in cutoffs (marker_name or e.g. 'CR' for 'CALB2')
    cutoff_key = marker_name
    if marker_name == "CALB2":
        cutoff_key = "CR"  # CALB2 is the gene, CR is the marker label for cutoff
    if cutoff_key in cutoffs:
        cutoff_val = cutoffs[cutoff_key]
    if cutoff_val is not None:
        ax.axvline(cutoff_val, color="k", linestyle="--", linewidth=1.5, label="Cutoff")
        # To avoid duplicate legends, only add cutoff to one handle per ax
        handles, labels = ax.get_legend_handles_labels()
        # Remove duplicate "Cutoff" if already present
        if labels.count("Cutoff") > 1:
            cutoff_indices = [i for i, label in enumerate(labels) if label == "Cutoff"]
            for i in cutoff_indices[1:]:
                del handles[i]
                del labels[i]
        ax.legend(handles, labels)

# Hide any unused subplots
for idx in range(n_markers, len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.show()


# %%
# Use cutoffs from the dict defined above (keyed by feature/column names)
VIP_Pos = Biases_CGE_ExpL[Biases_CGE_ExpL["VIP_ExpL"] > cutoffs["VIP"]]
VIP_Neg = Biases_CGE_ExpL[Biases_CGE_ExpL["VIP_ExpL"] < cutoffs["VIP"]]

# Access cutoffs for each marker
CCK_Cutoff = cutoffs["CCK"]
CR_Cutoff = cutoffs["CR"]
SNCG_Cutoff = cutoffs["SNCG"]
M2R_Cutoff = cutoffs["M2R"]

# Main subpopulation selections
CCKBC = Biases_CGE_ExpL[
    (Biases_CGE_ExpL["CR_ExpL"] < CR_Cutoff)
    & (Biases_CGE_ExpL["SNCG_ExpL"] > SNCG_Cutoff)
    & (Biases_CGE_ExpL["CCK_ExpL"] > CCK_Cutoff)
    & (Biases_CGE_ExpL["M2R_ExpL"] < M2R_Cutoff)
]
ISI2 = VIP_Pos[
    (VIP_Pos["CR_ExpL"] < CR_Cutoff)
    & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)
    & (VIP_Pos["SNCG_ExpL"] < SNCG_Cutoff)
]
ISI3 = VIP_Pos[
    (VIP_Pos["CR_ExpL"] > CR_Cutoff)
    & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)
    & (VIP_Pos["SNCG_ExpL"] < SNCG_Cutoff)
]

# Combine indices and get the "Other" group
combined_indices = np.concatenate([CCKBC.index.values, ISI2.index.values, ISI3.index.values])
Other = VIP_Pos.loc[~VIP_Pos.index.isin(combined_indices)]

print(len(CCKBC), len(ISI2), len(ISI3), len(Other))

# %%
CCKBC

# %%
ISI2

# %%
ISI3

# %%
Other

# %%
VIP_Neg


# %%
def plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias"):
    plt.figure(figsize=(8,6))

    # Create boxplot
    bp = plt.boxplot([CCKBC[label], ISI2[label], ISI3[label], Other[label], VIP_Neg[label]],
                     tick_labels=['CCKBC', 'ISI2', 'ISI3', 'Other', 'VIP-'])

    # Add individual points with jitter
    for i, data in enumerate([CCKBC[label], ISI2[label], ISI3[label], Other[label], VIP_Neg[label]], 1):
        x = np.random.normal(i, 0.04, size=len(data))
        plt.plot(x, data, 'o', alpha=0.5, color='black', markersize=4)

    plt.ylabel(label)
    plt.title(label)

    # Perform Kruskal-Wallis H-test
    h_stat, p_val = stats.kruskal(CCKBC[label], ISI2[label], ISI3[label], Other[label], VIP_Neg[label])
    plt.text(0.02, 0.98, f'Kruskal-Wallis p={p_val:.2e}', 
             transform=plt.gca().transAxes, verticalalignment='top')

    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def plot_beta_distribution(
    CCKBC, ISI2, ISI3, Other, VIP_Neg,
    label="X22q_Bias",
    figsize=(4.8, 3.6),
    jitter=0.06):
    fig, ax = plt.subplots(figsize=figsize)

    groups = [
        CCKBC[label],
        ISI2[label],
        ISI3[label],
        Other[label],
        VIP_Neg[label]
    ]
    group_labels = ['CCKBC', 'ISI2', 'ISI3', 'Other', 'VIP−']

    # --- Boxplot (lighter, transparent) ---
    bp = ax.boxplot(
        groups,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='black', linewidth=1.4),
        boxprops=dict(
            facecolor='lightgray',
            edgecolor='black',
            linewidth=1.2,
            alpha=0.35
        ),
        whiskerprops=dict(color='black', linewidth=1.1),
        capprops=dict(color='black', linewidth=1.1)
    )

    # --- Overlay individual points (larger, prominent) ---
    for i, data in enumerate(groups, 1):
        x = np.random.normal(i, jitter, size=len(data))
        ax.scatter(
            x, data,
            s=28,                # ⬅ bigger dots
            color='black',
            alpha=0.75,
            edgecolors='black',  # ⬅ outline for visibility
            linewidths=0.3,
            zorder=3
        )

    # --- Axes formatting ---
    ax.set_xticks(range(1, len(group_labels) + 1))
    ax.set_xticklabels(group_labels, fontsize=9)
    ax.set_ylabel(label, fontsize=9)
    ax.set_title(label, fontsize=10, pad=6)

    # Clean spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)

    # --- Statistics ---
    h_stat, p_kw = stats.kruskal(*groups)
    _, p_isi2 = stats.mannwhitneyu(
        ISI2[label], CCKBC[label], alternative='greater'
    )
    _, p_isi3 = stats.mannwhitneyu(
        ISI3[label], CCKBC[label], alternative='greater'
    )
    _, p_isi23_cck = stats.mannwhitneyu(
        list(ISI2[label]) + list(ISI3[label]), list(CCKBC[label]), alternative='greater'
    )

    _, p_isi23_vip_neg = stats.mannwhitneyu(
        list(ISI2[label]) + list(ISI3[label]), list(VIP_Neg[label]), alternative='greater'
    )
    # stat_text = (
    #     f"Kruskal–Wallis p = {p_kw:.2e}\n"
    #     f"ISI2 > CCKBC: p = {p_isi2:.2e}\n"
    #     f"ISI3 > CCKBC: p = {p_isi3:.2e}"
    # )

    stat_text = (
        f"Kruskal–Wallis p = {p_kw:.2e}\n"
        f"ISI2+ISI3 > CCKBC: p = {p_isi23_cck:.2e}\n"
        f"ISI2+ISI3 > VIP-: p = {p_isi23_vip_neg:.2e}\n"
    )

    ax.text(
        0.02, 0.98,
        stat_text,
        transform=ax.transAxes,
        fontsize=8,
        va='top'
    )

    plt.tight_layout()
    plt.show()



# %%
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_mouse_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Hexp_Bias")

plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day25_Bias")
#plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day50_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day75_Bias")
#plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day100_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="LowIQ_ASD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Pos_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Neg_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="DDD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="SCZ_Bias")

# %%
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_mouse_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Hexp_Bias")

plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day25_Bias")
#plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day50_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day75_Bias")
#plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_DEG_day100_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="LowIQ_ASD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Pos_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Neg_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="DDD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="SCZ_Bias")

# %% [markdown]
# ### Try CPM instead of LogUMI

# %%
Biases_CGE = TmpDF[TmpDF["Supercluster"]=="CGE interneuron"]
Biases_CGE_TPM = Biases_CGE.copy()
Biases_CGE_TPM.loc[:, "VIP_ExpL"]  = TPM.loc[GeneSymbol2Entrez["VIP"],  Biases_CGE_TPM.index].values
Biases_CGE_TPM.loc[:, "CCK_ExpL"]  = TPM.loc[GeneSymbol2Entrez["CCK"],  Biases_CGE_TPM.index].values
Biases_CGE_TPM.loc[:, "CR_ExpL"]   = TPM.loc[GeneSymbol2Entrez["CALB2"], Biases_CGE_TPM.index].values
Biases_CGE_TPM.loc[:, "SNCG_ExpL"] = TPM.loc[GeneSymbol2Entrez["SNCG"], Biases_CGE_TPM  .index].values
Biases_CGE_TPM.loc[:, "M2R_ExpL"]  = TPM.loc[GeneSymbol2Entrez["CHRM2"], Biases_CGE_TPM.index].values
Biases_CGE_TPM.loc[:, "CNR1_ExpL"] = TPM.loc[GeneSymbol2Entrez["CNR1"], Biases_CGE_TPM.index].values

# %%
# -- Fixed: Dynamically compute rows/cols based on number of features to avoid IndexError --

import math

# Only include unique features for plotting (remove duplicate "CNR1_ExpL")
plot_feature_names = []
seen_cols = set()
for col, label in feature_names:
    if col not in seen_cols:
        plot_feature_names.append((col, label))
        seen_cols.add(col)

nplots = len(plot_feature_names)
ncols = 3
nrows = math.ceil(nplots / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
if nrows * ncols > 1:
    axes = axes.flatten()
else:
    axes = [axes]
fig.suptitle("Expression Value Distributions (CGE interneurons)", fontsize=16)

for idx, (col, label) in enumerate(plot_feature_names):
    ax = axes[idx]
    if col in Biases_CGE_TPM.columns:
        ax.hist(Biases_CGE_TPM[col].dropna(), bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_title(label)
        ax.set_xlabel("CPM")
        ax.set_ylabel("Count")
        # Draw cutoff line if exists
        if cutoffs.get(col, None) is not None:
            ax.axvline(cutoffs[col], color='red', linestyle='--', linewidth=2, label=f"Cutoff = {cutoffs[col]}")
            ax.legend()
    else:
        ax.axis('off')

# Hide any unused axes (if feature count doesn't fill grid)
for j in range(idx + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# %%
Biases_CGE_TPM.sort_values(by="SNCG_ExpL", ascending=False, inplace=True)
Biases_CGE_TPM[Biases_CGE_TPM["VIP_ExpL"]>VIP_Cutoff]

# %%
VIP_Pos = Biases_CGE_TPM[Biases_CGE_TPM["VIP_ExpL"] > cutoffs["VIP_ExpL"]]
VIP_Neg = Biases_CGE_TPM[Biases_CGE_TPM["VIP_ExpL"] < cutoffs["VIP_ExpL"]]
# VIP_Cutoff = 100
# CCK_Cutoff = 1000
# CR_Cutoff = 200
# M2R_Cutoff = 100
# SNCG_Cutoff = 10
# CCKBC = VIP_Pos[(VIP_Pos["CCK_ExpL"] > CCK_Cutoff) & (VIP_Pos["CR_ExpL"] < CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]
# ISI2 = VIP_Pos[(VIP_Pos["CCK_ExpL"] < CCK_Cutoff) & (VIP_Pos["CR_ExpL"] < CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]
# ISI3 = VIP_Pos[(VIP_Pos["CCK_ExpL"] < CCK_Cutoff) & (VIP_Pos["CR_ExpL"] > CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]

CCKBC = Biases_CGE_TPM[(Biases_CGE_TPM["CR_ExpL"] < cutoffs["CR_ExpL"]) 
    & (Biases_CGE_TPM["SNCG_ExpL"] > cutoffs["SNCG_ExpL"]) & (Biases_CGE_TPM["CCK_ExpL"] > cutoffs["CCK_ExpL"])]
ISI2 = VIP_Pos[(VIP_Pos["CR_ExpL"] < cutoffs["CR_ExpL"]) & (VIP_Pos["M2R_ExpL"] < cutoffs["M2R_ExpL"]) & (VIP_Pos["SNCG_ExpL"] < cutoffs["SNCG_ExpL"])] #& (VIP_Pos["SNCG_ExpL"] < SNCG_Cutoff)
ISI3 = VIP_Pos[(VIP_Pos["CR_ExpL"] > cutoffs["CR_ExpL"]) & (VIP_Pos["M2R_ExpL"] < cutoffs["M2R_ExpL"]) & (VIP_Pos["SNCG_ExpL"] < cutoffs["SNCG_ExpL"])]

# Fix: Concatenate index values using np.concatenate instead of + operator
combined_indices = np.concatenate([CCKBC.index.values, ISI2.index.values, ISI3.index.values])
Other = VIP_Pos.loc[~VIP_Pos.index.isin(combined_indices)]

print(len(CCKBC), len(ISI2), len(ISI3), len(Other))

# %%
CCKBC

# %%
ISI2

# %%
ISI3

# %%
Other

# %%
VIP_Neg

# %%
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias")
plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_mouse_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="LowIQ_ASD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Pos_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="VNR_Neg_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="DDD_Bias")
# plot_beta_distribution(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="SCZ_Bias")

# %%
import matplotlib.pyplot as plt

def plot_gene_scatter(df, x_gene, y_gene, x_label=None, y_label=None, title=None, alpha=0.7, figsize=(4,4)):
    """
    Plots a scatter plot between two gene expression columns.

    Parameters:
    - df: DataFrame containing expression columns
    - x_gene: str, column name for x-axis
    - y_gene: str, column name for y-axis
    - x_label: str, label for x-axis (default: x_gene)
    - y_label: str, label for y-axis (default: y_gene)
    - title: str, plot title (default: 'Scatter plot of {x_gene} vs {y_gene} Expression')
    - alpha: float, transparency of points
    - figsize: tuple, figure size
    """
    plt.figure(figsize=figsize)
    plt.scatter(df[x_gene], df[y_gene], alpha=alpha)
    plt.xlabel(x_label if x_label is not None else f"{x_gene} Expression Level")
    plt.ylabel(y_label if y_label is not None else f"{y_gene} Expression Level")
    if title is not None:
        plt.title(title)
    else:
        plt.title(f"Scatter plot of {x_gene} vs {y_gene} Expression")
    plt.grid(True)
    plt.show()

# Example usage:
plot_gene_scatter(Biases_CGE_TPM, "CR_ExpL", "SNCG_ExpL")
plot_gene_scatter(Biases_CGE_TPM, "VIP_ExpL", "SNCG_ExpL")
plot_gene_scatter(Biases_CGE_TPM, "VIP_ExpL", "CCK_ExpL")
plot_gene_scatter(Biases_CGE_TPM, "VIP_ExpL", "CR_ExpL")
#plot_gene_scatter(Biases_CGE, "CR_ExpL", "SNCG_ExpL")

# %%
test = Biases_CGE[["VIP_ExpL", "CCK_ExpL", "CNR1_ExpL", "CR_ExpL", "M2R_ExpL"]]

# %%
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(test.fillna(0))  # Fill NaNs with 0 for PCA

# Plot the first two principal components
plt.figure(figsize=(6, 5))
plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7)
plt.xlabel('PC1 (%.2f%%)' % (pca.explained_variance_ratio_[0]*100))
plt.ylabel('PC2 (%.2f%%)' % (pca.explained_variance_ratio_[1]*100))
plt.title('PCA of VIP/CCK/CNR1/CR/M2R Expression')
plt.grid(True)
plt.show()

# %%
test

# %%
from sklearn.manifold import TSNE

# t-SNE on test (fill NaNs with 0)
tsne = TSNE(n_components=2, random_state=42, perplexity=5)
X_tsne = tsne.fit_transform(test.fillna(0))

def plot_tsne_mask(X_tsne, mask, mask_label='Masked', unmask_label='Unmasked', 
                   mask_color='red', unmask_color='blue', title='t-SNE of VIP/CCK/CNR1/CR/M2R Expression'):
    """
    Plots t-SNE results with a boolean mask to highlight points.

    Parameters:
    - X_tsne: np.ndarray, shape (n_samples, 2)
    - mask: boolean array, shape (n_samples,)
    - mask_label: label for masked points
    - unmask_label: label for unmasked points
    - mask_color: color for masked points
    - unmask_color: color for unmasked points
    - title: plot title
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6, 5))
    plt.scatter(X_tsne[~mask, 0], X_tsne[~mask, 1], alpha=0.7, color=unmask_color, label=unmask_label)
    plt.scatter(X_tsne[mask, 0], X_tsne[mask, 1], alpha=0.9, color=mask_color, label=mask_label)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()

# Example usage:
# Highlight points with CCK_ExpL > 3000 and CNR1_ExpL > 2000
mask = (test["CCK_ExpL"] > 3000) & (test["CNR1_ExpL"] > 2000 & (test["VIP_ExpL"] > 100))
plot_tsne_mask(
    X_tsne, 
    mask, 
    mask_label='CCK > 3000 & CNR1 > 2000', 
    unmask_label='CCK ≤ 3000 or CNR1 ≤ 2000'
)


# %%
from sklearn.manifold import TSNE

# t-SNE on test (fill NaNs with 0)
tsne = TSNE(n_components=2, random_state=42, perplexity=5)
X_tsne = tsne.fit_transform(test.fillna(0))

def plot_tsne_multi_mask(X_tsne, masks, labels=None, colors=None, title='t-SNE of VIP/CCK/CNR1/CR/M2R Expression'):
    """
    Plots t-SNE results with multiple boolean masks to highlight groups of points.

    Parameters:
    - X_tsne: np.ndarray, shape (n_samples, 2)
    - masks: list of boolean arrays, each shape (n_samples,)
    - labels: list of str, label for each mask
    - colors: list of str, color for each mask
    - title: plot title
    """
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(masks)
    if labels is None:
        labels = [f"Mask {i+1}" for i in range(n)]
    if colors is None:
        # Use some default matplotlib colors
        default_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        colors = default_colors[:n]

    plt.figure(figsize=(6, 5), dpi=120)
    already_masked = np.zeros(X_tsne.shape[0], dtype=bool)
    for mask, label, color in zip(masks, labels, colors):
        # Only plot points not already assigned to a previous mask
        mask_to_plot = mask & (~already_masked)
        plt.scatter(X_tsne[mask_to_plot, 0], X_tsne[mask_to_plot, 1], alpha=0.8, color=color, label=label)
        already_masked = already_masked | mask
    # Optionally plot the rest as "Other"
    if not np.all(already_masked):
        plt.scatter(X_tsne[~already_masked, 0], X_tsne[~already_masked, 1], alpha=0.3, color='lightgray', label='Other')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.grid(True)
    plt.show()



# %%

# Example usage:
# Define multiple masks
mask1 = (test["CCK_ExpL"] > 3000) & (test["CNR1_ExpL"] > 2000)
mask2 = (test["VIP_ExpL"] < 100)
mask3 = (test["CR_ExpL"] > 200)
masks = [mask1, mask2, mask3]
labels = ['CCK > 3000 & CNR1 > 2000', 'VIP-', 'CR > 200']
colors = ['red', 'blue', 'green']

plot_tsne_multi_mask(
    X_tsne,
    masks,
    labels=labels,
    colors=colors,
    title='t-SNE: Multiple Masks'
)

# %%
PBS_CGE = HumanCT_res_df_GeneL[HumanCT_res_df_GeneL["Supercluster"]=="CGE interneuron"].copy()
PBS_CGE["VIP_ExpL"] = TPM.loc[GeneSymbol2Entrez["VIP"], PBS_CGE.index]
PBS_CGE["SST_ExpL"] = TPM.loc[GeneSymbol2Entrez["SST"], PBS_CGE.index]
PBS_CGE["PVALB_ExpL"] = TPM.loc[GeneSymbol2Entrez["PVALB"], PBS_CGE.index]
PBS_CGE["NPY_ExpL"] = TPM.loc[GeneSymbol2Entrez["NPY"], PBS_CGE.index]
PBS_CGE["CCK_ExpL"] = TPM.loc[GeneSymbol2Entrez["CCK"], PBS_CGE.index]
PBS_CGE["CNR1_ExpL"] = TPM.loc[GeneSymbol2Entrez["CNR1"], PBS_CGE.index]
PBS_CGE["CR_ExpL"] = TPM.loc[GeneSymbol2Entrez["CALB2"], PBS_CGE.index]
PBS_CGE["M2R_ExpL"] = TPM.loc[GeneSymbol2Entrez["CHRM2"], PBS_CGE.index]

PBS_OtherIN = HumanCT_res_df_GeneL[
    (HumanCT_res_df_GeneL["Supercluster"] == "MGE interneuron") | 
    (HumanCT_res_df_GeneL["Supercluster"] == "LAMP5-LHX6 and Chandelier")
].copy()
PBS_OtherIN["VIP_ExpL"] = TPM.loc[GeneSymbol2Entrez["VIP"], PBS_OtherIN.index]
PBS_OtherIN["SST_ExpL"] = TPM.loc[GeneSymbol2Entrez["SST"], PBS_OtherIN.index]
PBS_OtherIN["PVALB_ExpL"] = TPM.loc[GeneSymbol2Entrez["PVALB"], PBS_OtherIN.index]
PBS_OtherIN["NPY_ExpL"] = TPM.loc[GeneSymbol2Entrez["NPY"], PBS_OtherIN.index]
PBS_OtherIN["CCK_ExpL"] = TPM.loc[GeneSymbol2Entrez["CCK"], PBS_OtherIN.index]
PBS_OtherIN["CNR1_ExpL"] = TPM.loc[GeneSymbol2Entrez["CNR1"], PBS_OtherIN.index]
PBS_OtherIN["CR_ExpL"] = TPM.loc[GeneSymbol2Entrez["CALB2"], PBS_OtherIN.index]
PBS_OtherIN["M2R_ExpL"] = TPM.loc[GeneSymbol2Entrez["CHRM2"], PBS_OtherIN.index]


# %%
fig, axes = plt.subplots(3, 3, figsize=(18, 12))
genes = ["VIP_ExpL", "SST_ExpL", "PVALB_ExpL", "NPY_ExpL", "CCK_ExpL", "CNR1_ExpL", "CR_ExpL", "M2R_ExpL"]
titles = ["VIP", "SST", "PVALB", "NPY", "CCK", "CNR1", "CR", "M2R"]

for ax, gene, title in zip(axes.flat, genes, titles):
    ax.hist(
        [np.log10(PBS_CGE[gene]), np.log10(PBS_OtherIN[gene])],
        label=["CGE", "Not CGE"],
        alpha=0.7,
        bins=20,
        density=True
    )
    ax.set_title(title)
    ax.set_xlabel('log10(TPM)')
    ax.set_ylabel('Count')
    ax.legend()

plt.tight_layout()
plt.show()

# %%
VIP_Cutoff = 100
CCK_Cutoff = 3000
CR_Cutoff = 100
M2R_Cutoff = 100

VIP_Pos = PBS_CGE[PBS_CGE["VIP_ExpL"] > VIP_Cutoff]
VIP_Neg = PBS_CGE[PBS_CGE["VIP_ExpL"] < VIP_Cutoff]

CCKBC = VIP_Pos[(VIP_Pos["CCK_ExpL"] > CCK_Cutoff) & (VIP_Pos["CR_ExpL"] < CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]
ISI2 = VIP_Pos[(VIP_Pos["CCK_ExpL"] < CCK_Cutoff) & (VIP_Pos["CR_ExpL"] < CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]
ISI3 = VIP_Pos[(VIP_Pos["CCK_ExpL"] < CCK_Cutoff) & (VIP_Pos["CR_ExpL"] > CR_Cutoff) & (VIP_Pos["M2R_ExpL"] < M2R_Cutoff)]

# Fix: Concatenate index values using np.concatenate instead of + operator
combined_indices = np.concatenate([CCKBC.index.values, ISI2.index.values, ISI3.index.values])
Other = VIP_Pos.loc[~VIP_Pos.index.isin(combined_indices)]

print(len(CCKBC), len(ISI2), len(ISI3), len(Other))

# %%
PBS_CGE[PBS_CGE["CCK_ExpL"] > 5000]

# %%
PBS_CGE[PBS_CGE["p_beta_perm_FDR"]>0.05]

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(6, 4))
Marker1 = "CCK"
Marker2 = "CR"
plt.scatter(PBS_CGE[f"{Marker1}_ExpL"], PBS_CGE[f"{Marker2}_ExpL"], alpha=0.7)
plt.xlabel(f"{Marker1} TPM")
plt.ylabel(f"{Marker2} TPM")
plt.grid(True)
plt.show()

# %%
