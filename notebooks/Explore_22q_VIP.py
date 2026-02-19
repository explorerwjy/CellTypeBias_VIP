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
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Exploratory: 22q11.2 and VIP Analysis
#
# CGE marker expression, VIP+/VIP- bias across disorders,
# multi-marker comparison, SCZ rm NDD.
# This is an exploratory notebook — not for final figures.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os

ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/"
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *

os.chdir(f"{ProjDIR}/notebooks/")
plt.style.use('seaborn-v0_8-whitegrid')

# %% [markdown]
# ## Load data

# %%
Bias_Save_Dir = os.path.join(ProjDIR, "results/main_results/random/Centering/")

HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
X22q_mousemodel = pd.read_csv(Bias_Save_Dir + "22q_small_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)

# Pre-computed CGE annotations
CONTRAST_DIR = os.path.join(ProjDIR, "results/main_results/contrasts/")
CGE_anno = pd.read_csv(CONTRAST_DIR + "CGE_VIP_annotation.csv", index_col=0)
CGE_anno_22q = pd.read_csv(CONTRAST_DIR + "CGE_VIP_annotation_22q.csv", index_col=0)

# %%
# Expression matrix for marker KDEs
ExpL = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv", index_col=0)
ExpL.columns = [int(x) for x in ExpL.columns.values]
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# %% [markdown]
# ## CGE Marker Expression KDE

# %%
CGE_idx = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values

CCK_expl = ExpL.loc[GeneSymbol2Entrez["CCK"], CGE_idx]
CALB2_expl = ExpL.loc[GeneSymbol2Entrez["CALB2"], CGE_idx]
VIP_expl = ExpL.loc[GeneSymbol2Entrez["VIP"], CGE_idx]
SNCG_expl = ExpL.loc[GeneSymbol2Entrez["SNCG"], CGE_idx]
LAMP5_expl = ExpL.loc[GeneSymbol2Entrez["LAMP5"], CGE_idx]

plt.figure(dpi=120)
sns.kdeplot(CCK_expl, label="CCK")
sns.kdeplot(CALB2_expl, label="CALB2")
sns.kdeplot(VIP_expl, label="VIP")
sns.kdeplot(SNCG_expl, label="SNCG")
sns.kdeplot(LAMP5_expl, label="LAMP5")
plt.xlabel("ExpL")
plt.legend()
plt.show()

# %%
plt.hist(CCK_expl, label="CCK")
plt.hist(CALB2_expl, label="CALB2")
plt.hist(VIP_expl, label="VIP")
plt.xlabel("ExpL")
plt.legend()
plt.show()

# %% [markdown]
# ## VIP+ vs VIP- bias across ALL disorders

# %%
from scipy.stats import mannwhitneyu

tmpdf = CGE_anno
cutoff = 1.0

VIP1 = [tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_22q.11"],
        tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_SCZ"],
        tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_LIQ"],
        tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_DD"],
        tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_VNR"],
        tmpdf[tmpdf["VIP"] > cutoff]["EFFECT_EDU"]]

VIP2 = [tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_22q.11"],
        tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_SCZ"],
        tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_LIQ"],
        tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_DD"],
        tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_VNR"],
        tmpdf[tmpdf["VIP"] < cutoff]["EFFECT_EDU"]]

ticks = ['22q.11', 'SCZ', 'ASD/ID', 'DD/ID', "VNR", "EDU"]

plt.figure(dpi=480)
summer_rain_plot = plt.boxplot(VIP1, positions=np.array(np.arange(len(VIP1)))*2.0-0.35, widths=0.6)
winter_rain_plot = plt.boxplot(VIP2, positions=np.array(np.arange(len(VIP2)))*2.0+0.35, widths=0.6)

def define_box_properties(plot_name, color_code, label):
    for k, v in plot_name.items():
        plt.setp(plot_name.get(k), color=color_code)
    plt.plot([], c=color_code, label=label)
    plt.legend()

define_box_properties(summer_rain_plot, '#D7191C', 'VIP+')
define_box_properties(winter_rain_plot, '#2C7BB6', 'VIP-')

plt.xticks(np.arange(0, len(ticks) * 2, 2), ticks)
plt.ylabel("Cell Type Bias")
plt.title('Cell Type Bias for CGE Interneurons')
plt.show()

# %% [markdown]
# ## Multi-marker comparison (VIP, CR, CCK, GABA)

# %%
cutoff = 1.0

VIP1 = [CGE_anno_22q[CGE_anno_22q["VIP"] > cutoff]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["CR"] > cutoff]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["CCK"] > 20]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["GABA"] > 0.7]["EFFECT_22q.11"]]

VIP2 = [CGE_anno_22q[CGE_anno_22q["VIP"] < cutoff]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["CR"] < cutoff]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["CCK"] < 20]["EFFECT_22q.11"],
        CGE_anno_22q[CGE_anno_22q["GABA"] < 0.7]["EFFECT_22q.11"]]

ticks = ['VIP', 'CR', 'CCK', 'GABA']

plt.figure(dpi=600, figsize=(8, 6))

summer_rain_plot = plt.boxplot(VIP1, positions=np.array(np.arange(len(VIP1)))*2.0-0.35,
                               widths=0.5, patch_artist=False,
                               boxprops=dict(color="red", linewidth=1.5),
                               medianprops=dict(color="darkred", linewidth=2),
                               whiskerprops=dict(color="red", linewidth=1.5),
                               capprops=dict(color="red", linewidth=1.5))

winter_rain_plot = plt.boxplot(VIP2, positions=np.array(np.arange(len(VIP2)))*2.0+0.35,
                               widths=0.5, patch_artist=False,
                               boxprops=dict(color="blue", linewidth=1.5),
                               medianprops=dict(color="darkblue", linewidth=2),
                               whiskerprops=dict(color="blue", linewidth=1.5),
                               capprops=dict(color="blue", linewidth=1.5))

plt.plot([], c='red', label='Marker+')
plt.plot([], c='blue', label='Marker-')
plt.legend(frameon=False, fontsize=15, loc='upper left')

for i in range(len(VIP1)):
    x1 = np.random.normal(i*2.0-0.35, 0.04, size=len(VIP1[i]))
    x2 = np.random.normal(i*2.0+0.35, 0.04, size=len(VIP2[i]))
    plt.scatter(x1, VIP1[i], color='#D7191C', alpha=0.6, s=20, edgecolor='k', linewidth=0.3)
    plt.scatter(x2, VIP2[i], color='#2C7BB6', alpha=0.6, s=20, edgecolor='k', linewidth=0.3)

    mean_diff = np.mean(VIP1[i]) - np.mean(VIP2[i])
    stat, p = mannwhitneyu(VIP1[i], VIP2[i])
    x1_pos, x2_pos = i*2.0-0.35, i*2.0+0.35
    y, h, col = max(max(VIP1[i]), max(VIP2[i])) + 0.1, max(max(VIP1[i]), max(VIP2[i])) * 0.05, 'k'
    plt.plot([x1_pos, x1_pos, x2_pos, x2_pos], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1_pos + x2_pos)*.5, y+h, f"P = {p:.2e}\n\u0394 = {mean_diff:.3f}",
             ha='center', va='bottom', color=col, fontsize=10, fontweight='bold')

plt.xticks(np.arange(0, len(ticks) * 2, 2), ticks, fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.ylabel("CGE Cluster Cell Type Bias", fontsize=14, fontweight='bold')

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## SCZ rm NDD (exploratory)

# %%
# These files are from an older analysis directory — only load if they exist
old_bias_dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Spec_Bias_Jul24_Centered/"
scz_rm_ndd61_path = old_bias_dir + "HCT.SCZ61.ExNDD61.csv"
scz_rm_ndd297_path = old_bias_dir + "HCT.SCZ61.ExNDD297.csv"

if os.path.exists(scz_rm_ndd61_path) and os.path.exists(scz_rm_ndd297_path):
    SCZ_rm_NDD61_Bias = pd.read_csv(scz_rm_ndd61_path, index_col=0)
    SCZ_rm_NDD297_Bias = pd.read_csv(scz_rm_ndd297_path, index_col=0)

    datasets_rm = {
        'ASD w/o ID': HighIQ_ASD_Bias,
        'SCZ rm NDD297': SCZ_rm_NDD297_Bias,
        'SCZ rm NDD61': SCZ_rm_NDD61_Bias,
        'SCZ': SCZ_Bias
    }
    TestPairs_rm = [("SCZ rm NDD297", "SCZ"), ("SCZ rm NDD61", "SCZ"),
                    ("SCZ rm NDD297", "ASD w/o ID"), ("SCZ rm NDD61", "ASD w/o ID")]
    plot_mutation_bias_comparison("CGE interneuron", datasets_rm, Anno, TestPairs=TestPairs_rm)
else:
    print("SCZ rm NDD files not found — skipping this analysis.")
