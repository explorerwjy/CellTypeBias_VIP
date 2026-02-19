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
import matplotlib.font_manager as fm
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)  # Only if you're adding a new font file
fm._load_fontmanager(try_read_cache=False)

# %% [markdown]
# # HumanCT IQ Phenotype
#

# %%
HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv", index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)

# %%
CogDF = pd.read_excel("../dat/41588_2023_1398_MOESM3_ESM.xlsx", sheet_name="Table S4")
CogDF = CogDF[CogDF["POPULATION"]=="EUR"]
for i, row in CogDF.iterrows():
    genesymbol = row["GENE"]
    if genesymbol in GeneSymbol2Entrez:
        CogDF.loc[i, "Entrez"] = GeneSymbol2Entrez[genesymbol]
    else:
        CogDF.loc[i, "Entrez"] = -1
CogDF = CogDF[CogDF["Entrez"] != -1]
# Entrez contert to int 
CogDF["Entrez"] = CogDF["Entrez"].astype(int)


# %%
VNR_DF = CogDF[CogDF["PHENOTYPE"]=="VNR"]
VNR_DF = VNR_DF[VNR_DF["Entrez"]!=-1]

RT_DF = CogDF[CogDF["PHENOTYPE"]=="RT"]
RT_DF = RT_DF[RT_DF["Entrez"]!=-1]

EDU_DF = CogDF[CogDF["PHENOTYPE"]=="EDU"]
EDU_DF = EDU_DF[EDU_DF["Entrez"]!=-1]

# %%
from scipy.stats import spearmanr

def BiasVsPheno_UKBB(MutPhenoDF, BiasMat, STR, label):
    biases = []
    IQs = []
    for i, row in MutPhenoDF.iterrows():
        if label == 'label':
            gene = row["HGNC"]
        else:
            gene = int(row["Entrez"])
        if gene in BiasMat.index.values:
            bias = BiasMat.loc[gene, STR]
            if bias == bias:
                IQ = row["BETA"]
                biases.append(bias)
                IQs.append(IQ)
    return biases, IQs

def Plot_Bias_vs_UKBBPheno(STR, Mut_n_IQ_conf, ylabel, HCT_Z2_MAT_HCT = HumanCT_Z2_HCT):
    biases, IQs = BiasVsPheno_UKBB(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , STR, 'XX')
    pho, p = spearmanr(biases, IQs)
    plt.figure(dpi=150, figsize=(5, 4))
    plt.scatter(biases, IQs, s=50, color="#2c7bb6", edgecolor="black", alpha=0.8, zorder=10)

    b, a = np.polyfit(biases, IQs, deg=1)
    xseq = np.linspace(min(biases), max(biases), num=100)
    plt.plot(xseq, a + b * xseq, color="#d7191c", lw=2.5, linestyle='--', zorder=5)
    _SuperCluster = Anno.loc[STR, "Supercluster"]
    # Add title with improved formatting
    plt.title(f'{_SuperCluster} - {STR} \nSpearman ρ = {pho:.2f}, p = {p:.2e}', fontsize=14, fontweight='bold')

    # Labeling axes
    plt.xlabel("Cell Type Bias", fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')

    # Grid lines for better readability
    plt.grid(True, linestyle='--', alpha=0.5)

    # Adjust tick parameters
    plt.xticks(fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    # Tight layout for better spacing
    plt.tight_layout()

    # Show the plot
    plt.show()


# %%
VNR_DF_highConf = VNR_DF[VNR_DF["P"]<0.05]
VNR_DF_highConf_Pos = VNR_DF_highConf[VNR_DF_highConf["BETA"]>0]
VNR_DF_highConf_Neg = VNR_DF_highConf[VNR_DF_highConf["BETA"]<0]

# %%
#VNR_DF_highConf = VNR_DF[VNR_DF["P"]<0.05]
VNR_DF_highConf = VNR_DF
VNR_DF_highConf_Pos = VNR_DF_highConf[VNR_DF_highConf["BETA"]>0]
VNR_DF_highConf_Neg = VNR_DF_highConf[VNR_DF_highConf["BETA"]<0]

# %%
Neuron_idx = Anno[Anno["Supercluster"].isin(Neurons)].index.values

# %%
from scipy.stats import pearsonr

def linear_fit(biases, IQs, alpha=0.05):
    model = sm.OLS(IQs, sm.add_constant(biases))
    results = model.fit()
    
    intercept = results.params[0]
    beta = results.params[1]
    # get CI of beta
    ci = results.conf_int(alpha=alpha)

    ci_low = ci[1][0]
    ci_high = ci[1][1]     
    r_value = results.rsquared
    p_value = results.pvalues[1]
    std_err = results.bse[1]
    rho, p_rho = spearmanr(biases, IQs)
    r, p_r = pearsonr(biases, IQs)
    
    return intercept, beta, ci_low, ci_high, r_value, p_value, std_err, rho, p_rho, r, p_r

def Make_HumanCT_DF(Mut_n_IQ_conf, HCT_Z2_MAT_HCT, output_file=None, alpha=0.05):
    names, supercluster, spearmanr, spearmanp, pearsonr_list, pearsonp, beta_values, beta_ci_low, beta_ci_high, intercept_values, r_value_values, p_value_values, std_err_values = [],[],[],[],[],[],[],[],[],[],[],[],[]
    for Idx in HCT_Z2_MAT_HCT.columns.values:
        biases, IQs = BiasVsPheno_UKBB(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , Idx, 'xx')
        intercept, beta, ci_low, ci_high, r_value, p_value, std_err, rho, p_rho, r, p_r = linear_fit(biases, IQs, alpha=0.05)
        
        names.append("{}".format(Idx))
        supercluster.append(Anno.loc[Idx, "Supercluster"])
        spearmanr.append(rho)
        spearmanp.append(p_rho)
        pearsonr_list.append(r)
        pearsonp.append(p_r)
        beta_values.append(beta)
        beta_ci_low.append(ci_low)
        beta_ci_high.append(ci_high)
        intercept_values.append(intercept)
        r_value_values.append(r_value)
        p_value_values.append(p_value)
        std_err_values.append(std_err)

    str_res_df = pd.DataFrame(data={"CT":names, "Supercluster":supercluster, "SpearmanR":spearmanr, "SpearmanP":spearmanp, 
                                            "PearsonR":pearsonr_list, "PearsonP":pearsonp, "beta":beta_values, "CI_low":beta_ci_low, "CI_high":beta_ci_high, "intercept":intercept_values, "r_value":r_value_values, 
                                            "p_value":p_value_values, "std_err":std_err_values})
    str_res_df = str_res_df.sort_values("SpearmanR")
    #str_res_df = ADJ_P(str_res_df)
    if output_file is not None:
        str_res_df.to_csv(output_file)
    return str_res_df

#human_ct_res_df_mutL.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.mutL.csv")


# %%
VNR_DF_highConf_Neg.head(2)

# %%
human_ct_df_VNR = Make_HumanCT_DF(VNR_DF_highConf, HumanCT_Z2_HCT)
human_ct_df_VNR_neg = Make_HumanCT_DF(VNR_DF_highConf_Neg, HumanCT_Z2_HCT)
human_ct_df_VNR_pos = Make_HumanCT_DF(VNR_DF_highConf_Pos, HumanCT_Z2_HCT)

# %%
human_ct_df_VNR.head(10)

# %%
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR, "VNR Spec", plot_metric="beta", sortby="median")

# %%
# def SuperClusterBias_BoxPlot_CorrIQ(DF1, flip_axis=True, plot_metric="beta", figsize=(8, 10), xlabel="", sortby="mean", FDR_label="FDR", Pval_label="p_beta_perm_Log"):
#     # Prepare data
#     dat_Z2 = []
#     mean_Z2 = []
#     for _CT in Neurons:
#         tmp = DF1[DF1["Supercluster"] == _CT]
#         dat_Z2.append(tmp[plot_metric].values)
#         if sortby == "mean":
#             mean_Z2.append(np.mean(tmp[plot_metric].values))
#         elif sortby == "median":
#             mean_Z2.append(np.median(tmp[plot_metric].values))
#     mean_Z2 = np.array(mean_Z2)
#     # Sort data###
#     if flip_axis:
#         sort_idx = np.argsort(mean_Z2)[::-1]
#     else:
#         sort_idx = np.argsort(mean_Z2)
#     show_dat_Z2 = [dat_Z2[x] for x in sort_idx]
#     show_CTs = [ALL_CTs[x] for x in sort_idx]

#     # Set publication-quality style
#     plt.style.use('seaborn-v0_8-paper')
#     sns.set_context("paper", font_scale=1.4)
    
#     # Create figure
#     fig, ax = plt.subplots(dpi=600, figsize=figsize)

#     # Define colors and styles
#     box_color = '#2E5A88'  # Darker blue
#     point_color = '#1f77b4'  # Seaborn blue
    
#     # Customize boxplot appearance
#     boxprops = dict(linestyle='-', linewidth=1.5, color=box_color, alpha=0.8)
#     medianprops = dict(linestyle='-', linewidth=2.5, color='#D62728')  # Red median line
#     whiskerprops = dict(color=box_color, linewidth=1.5, alpha=0.8)
#     capprops = dict(color=box_color, linewidth=1.5)
#     flierprops = dict(marker='', color=box_color, alpha=0)  # Hide outliers

#     # Draw boxplot
#     bp = ax.boxplot(show_dat_Z2, labels=show_CTs, vert=False, patch_artist=True,
#                     boxprops=boxprops, medianprops=medianprops,
#                     whiskerprops=whiskerprops, capprops=capprops, flierprops=flierprops)
    
#     # Fill boxes with lighter color
#     for patch in bp['boxes']:
#         patch.set_facecolor('#A8C8E4')  # Light blue
#         patch.set_alpha(0.3)

#     # Add individual points with jitter
#     for i in range(len(show_dat_Z2)):
#         y = np.random.normal(i + 1, 0.08, size=len(show_dat_Z2[i]))
#         ax.scatter(show_dat_Z2[i], y, s=20, color=point_color, alpha=0.4, 
#                   edgecolor='none', zorder=2)

#     # Customize grid
#     ax.grid(True, axis='x', linestyle='--', alpha=0.3, color='gray')
#     ax.grid(False, axis='y')

#     # Labels and title
#     if plot_metric == "beta":
#         ax.set_xlabel("PBS", fontsize=20, fontweight='normal')
#     elif plot_metric == "SpearmanR":
#         ax.set_xlabel("Spearman Correlation", fontsize=20, fontweight='normal')
#         #ax.set_ylabel("Superclusters", fontsize=20, fontweight='normal')
#     elif plot_metric == "p_beta_perm_Log" or plot_metric == "-log10(p)":
#         #ax.set_xlabel("PBS -log10(P)", fontsize=20, fontweight='normal')
#         ax.set_xlabel(xlabel, fontsize=20, fontweight='normal')
#         FDR_cut = GetFDRCut_PBS(DF1, FDR=0.05, FDR_label=FDR_label, Pval_label=Pval_label)
#         ax.axvline(x=FDR_cut, color='red', linestyle='--', linewidth=1, alpha=0.8)
#         #ax.text(FDR_cut, 0.95, "FDR=0.05", fontsize=12, fontweight='normal', ha='right', va='top', color='red', alpha=0.5)
#         ax.text(FDR_cut, 0.95, "FDR=0.05", fontsize=12, fontweight='normal', va='top', color='red', alpha=0.8)
#     else:
#         ax.set_xlabel(xlabel, fontsize=20, fontweight='normal')

#     # Customize axis appearance
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.spines['left'].set_visible(False)
    
#     # Adjust tick parameters
#     ax.tick_params(axis='both', which='major', labelsize=12)
#     ax.tick_params(axis='y', length=0)  # Remove y-axis ticks
    
#     # Flip axis if needed
#     if flip_axis:
#         ax.invert_xaxis()
    
#     if plot_metric == "beta":
#         #ax.set_xlim(7.5, -10)
#         pass



#     # Add zero line
#     ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
#     plt.tight_layout()
#     plt.show()

# %%
human_ct_df_VNR['-log10(p)'] = -np.log10(human_ct_df_VNR['p_value'])
human_ct_df_VNR["FDR"] = fdrcorrection(human_ct_df_VNR["p_value"])[1]
human_ct_df_VNR.head(20)
human_ct_df_VNR[human_ct_df_VNR["FDR"]<0.05].head(50)
human_ct_df_VNR.head(20)
human_ct_df_VNR['-log10(p)'] = -np.log10(human_ct_df_VNR['p_value'])
human_ct_df_VNR["FDR"] = fdrcorrection(human_ct_df_VNR["p_value"])[1]
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR, "VNR Spec", FDR_label="FDR", plot_metric="-log10(p)", sortby="median")

# %%
human_ct_df_VNR[human_ct_df_VNR["FDR"]<0.05]

# %%
#SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR_neg, "VNR Spec Neg", plot_metric="beta")
#SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR_pos, "VNR Spec Pos", plot_metric="beta")

# %%
EDU_DF_highConf = EDU_DF[EDU_DF["P"]<0.05]
RT_DF_highConf = RT_DF[RT_DF["P"]<0.05]
EDU_DF_highConf_Pos = EDU_DF_highConf[EDU_DF_highConf["BETA"]>0]
EDU_DF_highConf_Neg = EDU_DF_highConf[EDU_DF_highConf["BETA"]<0]

RT_DF_highConf_Pos = RT_DF_highConf[RT_DF_highConf["BETA"]>0]
RT_DF_highConf_Neg = RT_DF_highConf[RT_DF_highConf["BETA"]<0]


human_ct_df_EDU = Make_HumanCT_DF(EDU_DF_highConf, HumanCT_Z2_HCT)
human_ct_df_EDU_neg = Make_HumanCT_DF(EDU_DF_highConf_Neg, HumanCT_Z2_HCT)
human_ct_df_EDU_pos = Make_HumanCT_DF(EDU_DF_highConf_Pos, HumanCT_Z2_HCT)

human_ct_df_RT = Make_HumanCT_DF(RT_DF_highConf, HumanCT_Z2_HCT)
human_ct_df_RT_neg = Make_HumanCT_DF(RT_DF_highConf_Neg, HumanCT_Z2_HCT)
human_ct_df_RT_pos = Make_HumanCT_DF(RT_DF_highConf_Pos, HumanCT_Z2_HCT)

#SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR_neg, "VNR Spec Neg", plot_metric="beta")
#SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR_pos, "VNR Spec Pos", plot_metric="beta")

# %%
human_ct_df_EDU["FDR"] = fdrcorrection(human_ct_df_EDU["p_value"])[1]
human_ct_df_VNR["FDR"] = fdrcorrection(human_ct_df_VNR["p_value"])[1]
human_ct_df_RT["FDR"] = fdrcorrection(human_ct_df_RT["p_value"])[1]

# %%
human_ct_df_EDU.head(20)

# %%
human_ct_df_EDU[human_ct_df_EDU["FDR"]<0.05].head(50)

# %%
human_ct_df_VNR[human_ct_df_VNR["FDR"]<0.05].head(50)

# %%
human_ct_df_VNR.head(20)

# %%
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_EDU, "EDU Spec", plot_metric="beta", sortby="median")
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_RT, "RT Spec", plot_metric="beta", sortby="median")

# %%
human_ct_df_EDU['-log10(p)'] = -np.log10(human_ct_df_EDU['p_value'])
human_ct_df_VNR['-log10(p)'] = -np.log10(human_ct_df_VNR['p_value'])
human_ct_df_RT['-log10(p)'] = -np.log10(human_ct_df_RT['p_value'])


# %%
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_EDU, xlabel="EDU PBS -log10(P)", plot_metric="-log10(p)", sortby="median", flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_VNR, xlabel="VNR PBS -log10(P)", plot_metric="-log10(p)", sortby="median", flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")
SuperClusterBias_BoxPlot_CorrIQ(human_ct_df_RT, xlabel="RT PBS -log10(P)", plot_metric="-log10(p)", sortby="median", flip_axis=False, FDR_label="FDR", Pval_label="-log10(p)")

# %%
human_ct_df_EDU.head(10)


# %%
def CleanUpUKBBPheno(DF):
    res = DF.copy(deep=True)
    res = res[["CT", "Supercluster", "beta", "CI_low", "CI_high", "p_value", "FDR", "r_value"]]
    res = res.rename(columns={"CT": "Cluster", "Supercluster": "SuperCluster"})
    return res


# %%
human_ct_df_EDU_clean = CleanUpUKBBPheno(human_ct_df_EDU)
human_ct_df_VNR_clean = CleanUpUKBBPheno(human_ct_df_VNR)
human_ct_df_RT_clean = CleanUpUKBBPheno(human_ct_df_RT)
human_ct_df_EDU_clean.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.EDU.csv")
human_ct_df_VNR_clean.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.VNR.csv")
human_ct_df_RT_clean.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.RT.csv")

# %%
VIP_Anno = pd.read_csv("VIP_Anno.csv", index_col=0)

# %%
common_indices = human_ct_df_VNR_clean.index.intersection(VIP_Anno.index)
human_ct_df_VNR_clean_sub = human_ct_df_VNR_clean.loc[common_indices].copy()
human_ct_df_VNR_clean_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# %%
# Adjust code for VNR
common_indices = human_ct_df_VNR_clean.index.intersection(VIP_Anno.index)
human_ct_df_VNR_clean_sub = human_ct_df_VNR_clean.loc[common_indices].copy()
human_ct_df_VNR_clean_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# Filter for VNR phenotype
VNR_df = human_ct_df_VNR_clean_sub.copy()

# Split into VIP+ and VIP- groups
VNR_VIP_pos = VNR_df[VNR_df["VIP"] >= 1]
VNR_VIP_neg = VNR_df[VNR_df["VIP"] < 1]

# Prepare data for plotting
data = [VNR_VIP_pos["beta"], VNR_VIP_neg["beta"]]
from scipy.stats import mannwhitneyu

stat, pval = mannwhitneyu(VNR_VIP_pos["beta"], VNR_VIP_neg["beta"], alternative="greater")

# Boxplot with scatter
plt.figure(figsize=(5, 6))
plt.boxplot(data, labels=["VIP+", "VIP-"])
for i, d in enumerate(data):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect size (beta) for VNR")
plt.title(f"VNR: Mann-Whitney U p = {pval:.2e}")
plt.show()


# %%
# Adjust code for EDU
common_indices = human_ct_df_EDU_clean.index.intersection(VIP_Anno.index)
human_ct_df_EDU_clean_sub = human_ct_df_EDU_clean.loc[common_indices].copy()
human_ct_df_EDU_clean_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# Filter for EDU phenotype
EDU_df = human_ct_df_EDU_clean_sub.copy()

# Split into VIP+ and VIP- groups
EDU_VIP_pos = EDU_df[EDU_df["VIP"] >= 1]
EDU_VIP_neg = EDU_df[EDU_df["VIP"] < 1]

# Prepare data for plotting
data = [EDU_VIP_pos["beta"], EDU_VIP_neg["beta"]]

stat, pval = mannwhitneyu(EDU_VIP_pos["beta"], EDU_VIP_neg["beta"], alternative="greater")

# Boxplot with scatter
plt.figure(figsize=(5, 6))
plt.boxplot(data, labels=["VIP+", "VIP-"])
for i, d in enumerate(data):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect size (beta) for EDU")
plt.title(f"EDU: Mann-Whitney U p = {pval:.2e}")
plt.show()


# %%
# Adjust code for RT (Reaction Time)
common_indices = human_ct_df_RT_clean.index.intersection(VIP_Anno.index)
human_ct_df_RT_clean_sub = human_ct_df_RT_clean.loc[common_indices].copy()
human_ct_df_RT_clean_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# Filter for RT phenotype
RT_df = human_ct_df_RT_clean_sub.copy()

# Split into VIP+ and VIP- groups
RT_VIP_pos = RT_df[RT_df["VIP"] >= 1]
RT_VIP_neg = RT_df[RT_df["VIP"] < 1]

# Prepare data for plotting
data = [RT_VIP_pos["beta"], RT_VIP_neg["beta"]]

stat, pval = mannwhitneyu(RT_VIP_pos["beta"], RT_VIP_neg["beta"])

# Boxplot with scatter
plt.figure(figsize=(5, 6))
plt.boxplot(data, labels=["VIP+", "VIP-"])
for i, d in enumerate(data):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect size (beta) for RT")
plt.title(f"RT: Mann-Whitney U p = {pval:.2e}")
plt.show()


# %%

# %%

# %%

# %%
