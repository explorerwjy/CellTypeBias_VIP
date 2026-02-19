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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import sys

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
Mut_n_IQ = pd.read_csv("../dat/ASD_IQ_Mut.csv")
Mut_n_IQ.columns.values

# %%
Spark_Denovo = pd.read_excel("../dat/41588_2022_1148_MOESM4_ESM.xlsx",
                           skiprows=2, sheet_name="Table S7")
Spark_Denovo = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]!="."]
Spark_Denovo_ExomeWide = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]<=1.3e-6]
Spark_Denovo_ExomeWide.shape

# %%
top_Genes = Spark_Denovo.head(61)["HGNC"].values
Mut_n_IQ_conf = Mut_n_IQ[Mut_n_IQ["HGNC"].isin(top_Genes)]
Mut_n_IQ_conf.to_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv", index=False)
Mut_n_IQ_conf.shape

# %%
#### Gene Level
Genes = list(set(Mut_n_IQ_conf["Entrez"].values))
data = []
for g in Genes:
    tmp_df = Mut_n_IQ_conf[Mut_n_IQ_conf["Entrez"]==g]
    avg_IQ = tmp_df["IQ"].mean()
    row = [g, avg_IQ]
    data.append(row)
columns = ["Entrez", "IQ"]
Avg_Gene_IQ_DF = pd.DataFrame(data=data, columns=columns)
Avg_Gene_IQ_DF.to_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.GeneL.csv")

# %%
STR = 276
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HumanCT_Z2_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax2)
plt.tight_layout()
plt.show()

# %%
STR = 276
biases, IQs = BiasVsPheno(Avg_Gene_IQ_DF, HumanCT_Z2_HCT , STR, 'XX')
# save biases, IQs to file
np.savetxt("biases.txt", biases)
np.savetxt("IQs.txt", IQs)

# load biases, IQs from file 
biases = np.loadtxt("biases.txt")
IQs = np.loadtxt("IQs.txt")


# %%
STR = 276
plt.style.use('seaborn-v0_8-paper')
fig, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
#Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax1)
plt.tight_layout()
plt.show()

# %%
STR = 295
plt.style.use('seaborn-v0_8-paper')
fig, (ax1) = plt.subplots(1, 1, figsize=(7, 6))
#Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax1, Pval=0.037)
plt.tight_layout()
plt.show()

# %%
STR = 294
plt.style.use('seaborn-v0_8-paper')
fig, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
#Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HumanCT_Z2_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax1, Pval=0.0001)
plt.tight_layout()
plt.show()

# %%
STR = 276
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HumanCT_Z2_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax2)
plt.tight_layout()
plt.show()

# %%
STR = 212
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HumanCT_Z2_HCT, ax=ax1)
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax2)
plt.tight_layout()
plt.show()

# %%
Neuron_idx = Anno[Anno["Supercluster"].isin(Neurons)].index.values


# %%
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

def Make_HumanCT_DF(Mut_n_IQ_conf, HCT_Z2_MAT_HCT, output_file, alpha=0.05):
    names, supercluster, spearmanr, spearmanp, pearsonr, pearsonp, beta_values, beta_ci_low, beta_ci_high, intercept_values, r_value_values, p_value_values, std_err_values = [],[],[],[],[],[],[],[],[],[],[],[],[]
    for Idx in HCT_Z2_MAT_HCT.columns.values:
        biases, IQs = BiasVsPheno(Mut_n_IQ_conf, HCT_Z2_MAT_HCT , Idx, 'xx')
        intercept, beta, ci_low, ci_high, r_value, p_value, std_err, rho, p_rho, r, p_r = linear_fit(biases, IQs, alpha=0.05)
        
        names.append("{}".format(Idx))
        supercluster.append(Anno.loc[Idx, "Supercluster"])
        spearmanr.append(rho)
        spearmanp.append(p_rho)
        pearsonr.append(r)
        pearsonp.append(p_r)
        beta_values.append(beta)
        beta_ci_low.append(ci_low)
        beta_ci_high.append(ci_high)
        intercept_values.append(intercept)
        r_value_values.append(r_value)
        p_value_values.append(p_value)
        std_err_values.append(std_err)

    str_res_df = pd.DataFrame(data={"CT":names, "Supercluster":supercluster, "SpearmanR":spearmanr, "SpearmanP":spearmanp, 
                                            "PearsonR":pearsonr, "PearsonP":pearsonp, "beta":beta_values, "CI_low":beta_ci_low, "CI_high":beta_ci_high, "intercept":intercept_values, "r_value":r_value_values, 
                                            "p_value":p_value_values, "std_err":std_err_values})
    str_res_df = str_res_df.sort_values("SpearmanR")
    #str_res_df = ADJ_P(str_res_df)
    str_res_df.to_csv(output_file)
    return str_res_df


# %%
HumanCT_res_df_MutL = Make_HumanCT_DF(Mut_n_IQ_conf, HumanCT_Z2_HCT, "../dat/Pheno_Bias_vs_IQ/HumanCT.spec.MutL.csv")
HumanCT_res_df_GeneL = Make_HumanCT_DF(Avg_Gene_IQ_DF, HumanCT_Z2_HCT, "../dat/Pheno_Bias_vs_IQ/HumanCT.spec.GeneL.csv")


# %%
Mut_n_IQ_conf

# %%
HumanCT_res_df_GeneL.sort_values("SpearmanR").head(20)

# %%
HumanCT_res_df_GeneL.sort_values("beta").head(20)

# %%
HumanCT_res_df_GeneL.sort_values("p_value").head(20)

# %%
print(np.mean(HumanCT_res_df_GeneL[HumanCT_res_df_GeneL["Supercluster"]=="CGE interneuron"]["beta"]))
print(np.mean(HumanCT_res_df_GeneL[HumanCT_res_df_GeneL["Supercluster"]=="LAMP5-LHX6 and Chandelier"]["beta"]))
print(np.mean(HumanCT_res_df_GeneL[HumanCT_res_df_GeneL["Supercluster"]=="MGE interneuron"]["beta"]))

# %%
HumanCT_res_df_GeneL.head(20)

# %%
SuperClusterBias_BoxPlot_CorrIQ(HumanCT_res_df_GeneL, flip_axis=True, figsize=(6, 8), plot_metric="SpearmanR")

# %%
SuperClusterBias_BoxPlot_CorrIQ(HumanCT_res_df_GeneL, flip_axis=True, figsize=(6, 8), plot_metric="beta")

# %%
VIP_Anno = pd.read_csv("VIP_Anno.csv", index_col=0)

# %%
# Get subset of rows from VIP_Anno that exist in X22q_Z2_Bias_PT
common_indices = HumanCT_res_df_GeneL.index.intersection(VIP_Anno.index)
HumanCT_res_df_GeneL_sub = HumanCT_res_df_GeneL.loc[common_indices].copy()
HumanCT_res_df_GeneL_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# %%
HumanCT_res_df_GeneL_sub_VIP_pos = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] >= 1]
HumanCT_res_df_GeneL_sub_VIP_neg = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] < 1]
# plot effect of VIP+ vs VIP-
data = [HumanCT_res_df_GeneL_sub_VIP_pos["beta"], HumanCT_res_df_GeneL_sub_VIP_neg["beta"]]
# Perform Mann-Whitney U test
stat, pval = scipy.stats.mannwhitneyu(HumanCT_res_df_GeneL_sub_VIP_pos["beta"], 
                                    HumanCT_res_df_GeneL_sub_VIP_neg["beta"], alternative="less")
# Create boxplot with individual points
bp = plt.boxplot(data, labels=["VIP+", "VIP-"])
# Add scatter points
for i, d in enumerate([HumanCT_res_df_GeneL_sub_VIP_pos["beta"], HumanCT_res_df_GeneL_sub_VIP_neg["beta"]]):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect")
plt.title(f"p = {pval:.2e}")
plt.show()


# %%
Effect = "SpearmanR"
HumanCT_res_df_GeneL_sub_VIP_pos = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] >= 1]
HumanCT_res_df_GeneL_sub_VIP_neg = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] < 1]
# plot effect of VIP+ vs VIP-
data = [HumanCT_res_df_GeneL_sub_VIP_pos[Effect], HumanCT_res_df_GeneL_sub_VIP_neg[Effect]]
# Perform Mann-Whitney U test
stat, pval = scipy.stats.mannwhitneyu(HumanCT_res_df_GeneL_sub_VIP_pos[Effect], 
                                    HumanCT_res_df_GeneL_sub_VIP_neg[Effect])
# Create boxplot with individual points
bp = plt.boxplot(data, labels=["VIP+", "VIP-"])
# Add scatter points
for i, d in enumerate([HumanCT_res_df_GeneL_sub_VIP_pos[Effect], HumanCT_res_df_GeneL_sub_VIP_neg[Effect]]):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect")
plt.title(f"p = {pval:.2e}")
plt.show()

# %%

# %% [markdown]
# ### Mut IQ permutation

# %%
#Perm_DIR = "/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT_Feb25/ALL/"
#Perm_DIR = "/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT_Feb25_Seed42/ALL/"
Perm_DIR = "/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT_June13/ALL/"
Perm_DFs = []
for df in os.listdir(Perm_DIR):
    df = pd.read_csv(f"{Perm_DIR}/{df}", index_col=0)
    df.index = df.index.astype(int)
    Perm_DFs.append(df)
print(len(Perm_DFs))


# %%
def plot_null_distributions(HumanCT_res_df_GeneL, Perm_DFs, CT=276, plot=False):
    """
    Plot null distributions of Rho and Beta values compared to observed values for a given cell type.
    
    Parameters:
    -----------
    HumanCT_res_df_GeneL : pandas DataFrame
        DataFrame containing observed results
    Perm_DFs : list
        List of DataFrames containing permutation results
    CT : int
        Cell type ID to analyze (default: 276)
    """
    Obs_Rho = HumanCT_res_df_GeneL.loc[CT, "SpearmanR"]
    Obs_Beta = HumanCT_res_df_GeneL.loc[CT, "beta"]
    
    Null_Rho = []
    Null_Beta = []
    for df in Perm_DFs:
        Null_Rho.append(df.loc[CT, "SpearmanR"])
        Null_Beta.append(df.loc[CT, "beta"])
    Null_Rho = np.array(Null_Rho)
    Null_Beta = np.array(Null_Beta)
    # Calculate p-values
    p_rho = (Null_Rho <= Obs_Rho).mean()
    p_beta = (Null_Beta <= Obs_Beta).mean()
    
    if plot:
        # Create plots comparing observed values to null distributions
        plt.figure(figsize=(10,5))

        # Plot Rho distribution and observed value
        plt.subplot(121)
        plt.hist(Null_Rho, bins=20, alpha=0.5)
        plt.axvline(x=Obs_Rho, color='red', linestyle='--', label='Observed')
        plt.xlabel('Rho')
        plt.ylabel('Count')
        plt.text(0.60, 0.60, f'p_rho = {p_rho:.2e}', transform=plt.gca().transAxes, fontsize=12, ha='right', va='top')
        plt.title('Rho: Null Distribution vs Observed')
        plt.legend()

        # Plot Beta distribution and observed value  
        plt.subplot(122)
        plt.hist(Null_Beta, bins=20, alpha=0.5)
        plt.axvline(x=Obs_Beta, color='red', linestyle='--', label='Observed')
        plt.xlabel('Beta')
        plt.ylabel('Count')
        plt.title('Beta: Null Distribution vs Observed')
        plt.legend()
    # annotate p-values

        plt.text(0.60, 0.60, f'p_beta = {p_beta:.2e}', transform=plt.gca().transAxes, fontsize=12, ha='right', va='top')
        plt.tight_layout()
        plt.show()
    
    return p_rho, p_beta

def plot_null_suptercluster_distributions(ClusterIdx, HumanCT_res_df_GeneL, Perm_DFs, plot=False):

    Obs_Rho = HumanCT_res_df_GeneL.loc[ClusterIdx, "SpearmanR"].mean()
    Obs_Beta = HumanCT_res_df_GeneL.loc[ClusterIdx, "beta"].mean()
    

    Null_Rho = []
    Null_Beta = []
    for df in Perm_DFs:
        Null_Rho.append(df.loc[ClusterIdx, "SpearmanR"].mean())
        Null_Beta.append(df.loc[ClusterIdx, "beta"].mean())
    Null_Rho = np.array(Null_Rho)
    Null_Beta = np.array(Null_Beta)
    # Calculate p-values
    p_rho = (Null_Rho <= Obs_Rho).mean()
    p_beta = (Null_Beta <= Obs_Beta).mean()

    
    if plot:
        # Create plots comparing observed values to null distributions
        plt.figure(figsize=(10,5))

        # Plot Rho distribution and observed value
        plt.subplot(121)
        plt.hist(Null_Rho, bins=20, alpha=0.5)
        plt.axvline(x=Obs_Rho, color='red', linestyle='--', label='Observed')
        plt.xlabel('Rho')
        plt.ylabel('Count')
        plt.text(0.60, 0.60, f'p_rho = {p_rho:.2e}', transform=plt.gca().transAxes, fontsize=12, ha='right', va='top')
        plt.title('Rho: Null Distribution vs Observed')
        plt.legend()

        # Plot Beta distribution and observed value  
        plt.subplot(122)
        plt.hist(Null_Beta, bins=20, alpha=0.5)
        plt.axvline(x=Obs_Beta, color='red', linestyle='--', label='Observed')
        plt.xlabel('Beta')
        plt.ylabel('Count')
        plt.title('Beta: Null Distribution vs Observed')
        plt.legend()
    # annotate p-values

        plt.text(0.60, 0.60, f'p_beta = {p_beta:.2e}', transform=plt.gca().transAxes, fontsize=12, ha='right', va='top')
        plt.tight_layout()
        plt.show()
    
    return p_rho, p_beta

# %%
plot_null_distributions(HumanCT_res_df_GeneL, Perm_DFs, CT=295, plot=True)

# %%
Supercluster = "CGE interneuron"
ClusterIdx = Anno[Anno["Supercluster"]==Supercluster].index.values

p_rho, p_beta = plot_null_suptercluster_distributions(ClusterIdx, HumanCT_res_df_GeneL, Perm_DFs, plot=True)

# %%
#Supercluster = "CGE interneuron"
# Create DataFrame to store results
results_df = pd.DataFrame(columns=['Supercluster', 'p_rho', 'p_beta'])

for Supercluster in Neurons:
    ClusterIdx = Anno[Anno["Supercluster"]==Supercluster].index.values
    p_rho, p_beta = plot_null_suptercluster_distributions(ClusterIdx, HumanCT_res_df_GeneL, Perm_DFs, plot=False)
    
    # Add results to DataFrame
    results_df = results_df.append({
        'Supercluster': Supercluster,
        "mean_PBS": HumanCT_res_df_GeneL.loc[ClusterIdx, "beta"].mean(),
        'p_rho': p_rho,
        'p_beta': p_beta
    }, ignore_index=True)

# %%
# Add FDR corrected p-values
results_df['p_rho_FDR'] = multipletests(results_df['p_rho'], method='fdr_bh')[1]
results_df['p_beta_FDR'] = multipletests(results_df['p_beta'], method='fdr_bh')[1]

results_df.sort_values("p_beta")

# %%

# %%
for i, row in HumanCT_res_df_GeneL.iterrows():
    ct = int(row["CT"])
    p_rho, p_beta = plot_null_distributions(HumanCT_res_df_GeneL, Perm_DFs, CT=ct, plot=False)
    HumanCT_res_df_GeneL.loc[i, "p_rho_perm"] = p_rho
    HumanCT_res_df_GeneL.loc[i, "p_beta_perm"] = p_beta
#HumanCT_res_df_GeneL.sort_values("p_beta_FDR")

# %%
HumanCT_res_df_GeneL_neuron = HumanCT_res_df_GeneL[HumanCT_res_df_GeneL["Supercluster"].isin(Neurons)]

# %%
HumanCT_res_df_GeneL = HumanCT_res_df_GeneL.sort_values("p_beta_perm")
HumanCT_res_df_GeneL_neuron = HumanCT_res_df_GeneL_neuron.sort_values("p_beta_perm")
print(HumanCT_res_df_GeneL.shape, HumanCT_res_df_GeneL_neuron.shape)


# %%
HumanCT_res_df_GeneL["p_beta_perm_Log"] = HumanCT_res_df_GeneL["p_beta_perm"].apply(lambda x: -np.log10(x))

# %%
HumanCT_res_df_GeneL = HumanCT_res_df_GeneL.sort_values("p_beta_perm")
HumanCT_res_df_GeneL.head(5)

# %%
#umanCT_res_df_GeneL
SuperClusterBias_BoxPlot_CorrIQ(HumanCT_res_df_GeneL, flip_axis=True, figsize=(6, 8), plot_metric="p_beta_perm_Log")

# %%
HumanCT_res_df_GeneL['p_rho_perm_FDR'] = multipletests(HumanCT_res_df_GeneL['p_rho_perm'], method='fdr_bh')[1]
HumanCT_res_df_GeneL['p_beta_perm_FDR'] = multipletests(HumanCT_res_df_GeneL['p_beta_perm'], method='fdr_bh')[1]

HumanCT_res_df_GeneL_neuron['p_rho_perm_FDR'] = multipletests(HumanCT_res_df_GeneL_neuron['p_rho_perm'], method='fdr_bh')[1]
HumanCT_res_df_GeneL_neuron['p_beta_perm_FDR'] = multipletests(HumanCT_res_df_GeneL_neuron['p_beta_perm'], method='fdr_bh')[1]

# %%
HumanCT_res_df_GeneL_neuron.head(10)

# %%
HumanCT_res_df_GeneL = HumanCT_res_df_GeneL.sort_values("beta")

# %%
HumanCT_res_df_GeneL.head(10)

# %%
HumanCT_res_df_GeneL.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.cluster.June10.csv")
results_df.to_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.Supercluster.June10.csv")

# %%
HumanCT_res_df_GeneL.head(2)

# %%
MGE interneuron

# %%
STR = 295
plt.style.use('seaborn-v0_8-paper')
fig, (ax1) = plt.subplots(1, 1, figsize=(8, 6))
#Plot_Bias_vs_IQ_HumanCT(STR, Mut_n_IQ_conf, HCT_Z2_MAT_HCT, ax=ax1)
#Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax1, Pval=0.015367, textPos=(0.1, 0.2))
Plot_Bias_vs_IQ_HumanCT(STR, Avg_Gene_IQ_DF, HumanCT_Z2_HCT, ax=ax1, Pval=0.0002, textPos=(0.1, 0.2))
plt.tight_layout()
plt.show()

# %%
VIP_Anno = pd.read_csv("VIP_Anno.csv", index_col=0)
common_indices = HumanCT_res_df_GeneL.index.intersection(VIP_Anno.index)
HumanCT_res_df_GeneL_sub = HumanCT_res_df_GeneL.loc[common_indices].copy()
HumanCT_res_df_GeneL_sub["VIP"] = VIP_Anno.loc[common_indices, "VIP"]

# %%
HumanCT_res_df_GeneL_sub_VIP_pos = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] >= 1]
HumanCT_res_df_GeneL_sub_VIP_neg = HumanCT_res_df_GeneL_sub[HumanCT_res_df_GeneL_sub["VIP"] < 1]
# plot effect of VIP+ vs VIP-
EFFECT = "p_beta_perm"
data = [HumanCT_res_df_GeneL_sub_VIP_pos[EFFECT], HumanCT_res_df_GeneL_sub_VIP_neg[EFFECT]]
# Perform Mann-Whitney U test
stat, pval = scipy.stats.mannwhitneyu(HumanCT_res_df_GeneL_sub_VIP_pos[EFFECT], 
                                    HumanCT_res_df_GeneL_sub_VIP_neg[EFFECT])
# Create boxplot with individual points
bp = plt.boxplot(data, labels=["VIP+", "VIP-"])
# Add scatter points
for i, d in enumerate([HumanCT_res_df_GeneL_sub_VIP_pos[EFFECT], HumanCT_res_df_GeneL_sub_VIP_neg[EFFECT]]):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)
plt.ylabel("Effect")
plt.title(f"p = {pval:.2e}")
plt.show()

# %%
HumanCT_res_df_GeneL_sub.sort_values(by="p_beta_perm_FDR", ascending=True).head(15)

# %%
HumanCT_res_df_GeneL_sub
FDR_cutoff = 0.05
N_VIP_Sig = HumanCT_res_df_GeneL_sub_VIP_pos[HumanCT_res_df_GeneL_sub_VIP_pos["p_beta_perm_FDR"]< FDR_cutoff].shape[0]
N_nonVIP_Sig = HumanCT_res_df_GeneL_sub_VIP_neg[HumanCT_res_df_GeneL_sub_VIP_neg["p_beta_perm_FDR"]< FDR_cutoff].shape[0]
N_total_VIP = HumanCT_res_df_GeneL_sub_VIP_pos.shape[0]
N_total_nonVIP = HumanCT_res_df_GeneL_sub_VIP_neg.shape[0]
print(f"VIP+: {N_VIP_Sig}/{N_total_VIP}, VIP-: {N_nonVIP_Sig}/{N_total_nonVIP}")
print(f"VIP+: {N_VIP_Sig/N_total_VIP}, VIP-: {N_nonVIP_Sig/N_total_nonVIP}")

# %%
# Test if N_VIP_Sig is significantly more than expected by chance using a binomial test
from scipy.stats import binom_test

# Assume the expected proportion of significant genes is the same as in VIP- group
expected_prop = N_nonVIP_Sig / N_total_nonVIP if N_total_nonVIP > 0 else 0
#expected_prop = (2+10) / (14 + 7)
print(expected_prop)

# Perform one-sided binomial test: is N_VIP_Sig greater than expected under null?
p_binom = binom_test(N_VIP_Sig, N_total_VIP, expected_prop, alternative='greater')

print(f"Binomial test p-value (VIP+ > expected): {p_binom:.3g}")
N_VIP_Sig, N_nonVIP_Sig, N_total_VIP, N_total_nonVIP

# %%
# Test if VIP+ is associated with being significant using Fisher's exact test
from scipy.stats import fisher_exact

# Construct contingency table:
#            Significant   Not Significant
# VIP+       N_VIP_Sig     N_total_VIP - N_VIP_Sig
# VIP-       N_nonVIP_Sig  N_total_nonVIP - N_nonVIP_Sig
contingency_table = [
    [N_VIP_Sig, N_total_VIP - N_VIP_Sig],
    [N_nonVIP_Sig, N_total_nonVIP - N_nonVIP_Sig]
]

# Perform Fisher's exact test (test for enrichment of significance in VIP+)
oddsratio, p_fisher = fisher_exact(contingency_table, alternative='greater')

print(f"Fisher's exact test p-value (VIP+ enrichment for significance): {p_fisher:.3g}")
print("Contingency table:")
print(f"         Significant  Not Significant")
print(f"VIP+     {N_VIP_Sig:<12} {N_total_VIP - N_VIP_Sig}")
print(f"VIP-     {N_nonVIP_Sig:<12} {N_total_nonVIP - N_nonVIP_Sig}")
print(f"Odds ratio: {oddsratio:.3g}")
N_VIP_Sig, N_nonVIP_Sig, N_total_VIP, N_total_nonVIP

# %%
# Is Fisher's exact test the most powerful here?
# For small sample sizes or when expected counts are low, Fisher's exact test is appropriate and provides an exact p-value.
# For larger sample sizes, a chi-squared test of independence may have more power, but can be inaccurate if expected counts are small.
# Here, we use both tests for comparison.

from scipy.stats import fisher_exact, chi2_contingency

# Construct contingency table:
#            Significant   Not Significant
# VIP+       N_VIP_Sig     N_total_VIP - N_VIP_Sig
# VIP-       N_nonVIP_Sig  N_total_nonVIP - N_nonVIP_Sig
contingency_table = [
    [N_VIP_Sig, N_total_VIP - N_VIP_Sig],
    [N_nonVIP_Sig, N_total_nonVIP - N_nonVIP_Sig]
]

# Fisher's exact test (exact, best for small samples or low expected counts)
oddsratio, p_fisher = fisher_exact(contingency_table, alternative='greater')

# Chi-squared test (approximate, more powerful for large samples)
chi2, p_chi2, dof, expected = chi2_contingency(contingency_table, correction=False)

print(f"Fisher's exact test p-value (VIP+ enrichment for significance): {p_fisher:.3g}")
print(f"Chi-squared test p-value: {p_chi2:.3g}")
print("Contingency table:")
print(f"         Significant  Not Significant")
print(f"VIP+     {N_VIP_Sig:<12} {N_total_VIP - N_VIP_Sig}")
print(f"VIP-     {N_nonVIP_Sig:<12} {N_total_nonVIP - N_nonVIP_Sig}")
print(f"Odds ratio (Fisher): {oddsratio:.3g}")
print(f"Expected counts (chi2):\n{expected}")
N_VIP_Sig, N_nonVIP_Sig, N_total_VIP, N_total_nonVIP

# %%
# Barnard's Exact Test for 2x2 tables (more powerful than Fisher's for some cases)
# Note: scipy does not implement Barnard's test as of 2024. Use 'barnard_exact' from statsmodels if available.

from scipy.stats import fisher_exact
try:
    from statsmodels.stats.contingency_tables import barnard_exact
    barnard_available = True
except ImportError:
    barnard_available = False

# Construct contingency table:
#            Significant   Not Significant
# VIP+       N_VIP_Sig     N_total_VIP - N_VIP_Sig
# VIP-       N_nonVIP_Sig  N_total_nonVIP - N_nonVIP_Sig
contingency_table = [
    [N_VIP_Sig, N_total_VIP - N_VIP_Sig],
    [N_nonVIP_Sig, N_total_nonVIP - N_nonVIP_Sig]
]

# Fisher's exact test (for comparison)
oddsratio, p_fisher = fisher_exact(contingency_table, alternative='greater')

print(f"Fisher's exact test p-value (VIP+ enrichment for significance): {p_fisher:.3g}")
print(f"Contingency table:")
print(f"         Significant  Not Significant")
print(f"VIP+     {N_VIP_Sig:<12} {N_total_VIP - N_VIP_Sig}")
print(f"VIP-     {N_nonVIP_Sig:<12} {N_total_nonVIP - N_nonVIP_Sig}")
print(f"Odds ratio (Fisher): {oddsratio:.3g}")

# Barnard's exact test (if available)
if barnard_available:
    res = barnard_exact(contingency_table, alternative='greater')
    print(f"Barnard's exact test p-value (VIP+ enrichment for significance): {res.pvalue:.3g}")
    print(f"Barnard's test statistic: {res.statistic:.3g}")
else:
    print("Barnard's exact test is not available (requires statsmodels >= 0.13.0).")

N_VIP_Sig, N_nonVIP_Sig, N_total_VIP, N_total_nonVIP

# %%
import scipy.stats as stats
res = stats.fisher_exact([[11, 3], [2, 5]], alternative="greater")
print(res.pvalue)

# %%
import scipy.stats as stats
res = stats.barnard_exact([[10, 4], [2, 5]], alternative="greater")
print(res.pvalue)

# %%
import scipy.stats as stats
res = stats.boschloo_exact([[10, 4], [2, 5]], alternative="greater")
print(res.pvalue)


# %%
def permutation_test(group1_success, group1_total, group2_success, group2_total, n_permutations=100000):
    """
    Permutation test for difference in proportions.
    """
    # Original data
    group1_data = [1] * group1_success + [0] * (group1_total - group1_success)
    group2_data = [1] * group2_success + [0] * (group2_total - group2_success)
    
    # Combine all data
    all_data = group1_data + group2_data
    n1, n2 = len(group1_data), len(group2_data)
    
    # Observed difference
    obs_diff = (group1_success / group1_total) - (group2_success / group2_total)
    
    print(f"Observed difference: {obs_diff:.4f}")
    print(f"Running {n_permutations:,} permutations...")
    
    # Permutation test
    extreme_count = 0
    random.seed(42)  # For reproducibility
    
    for i in range(n_permutations):
        # Randomly shuffle the combined data
        shuffled = random.sample(all_data, len(all_data))
        
        # Split back into two groups
        perm_group1 = shuffled[:n1]
        perm_group2 = shuffled[n1:]
        
        # Calculate difference for this permutation
        perm_diff = (sum(perm_group1) / n1) - (sum(perm_group2) / n2)
        
        # Count if as extreme or more extreme than observed
        if perm_diff >= obs_diff:  # One-tailed test
            extreme_count += 1
    
    p_value_perm = extreme_count / n_permutations
    
    return p_value_perm, obs_diff

# Run permutation test
p_perm, obs_diff = permutation_test(10, 14, 2, 7, n_permutations=1000000)
print(p_perm)

# %%

# %% [markdown]
# #### Test VIP+ and VIP-

# %%
CGE_VIP_Pos_list = np.loadtxt("../dat/Other/CGE_VIP_Pos.txt", dtype=int)
CGE_VIP_Neg_list = np.loadtxt("../dat/Other/CGE_VIP_Neg.txt", dtype=int)

# %%
VIP_pos_meanPBS = HumanCT_res_df_GeneL[HumanCT_res_df_GeneL.index.isin(CGE_VIP_Pos_list)]["beta"].mean()
VIP_neg_meanPBS = HumanCT_res_df_GeneL[HumanCT_res_df_GeneL.index.isin(CGE_VIP_Neg_list)]["beta"].mean()
DIFF_obs = VIP_pos_meanPBS - VIP_neg_meanPBS
print(VIP_pos_meanPBS, VIP_neg_meanPBS, VIP_pos_meanPBS - VIP_neg_meanPBS)

# %%
null_diffs = []
for i in range(10000):
    tmpDF = Perm_DFs[i]
    tmp_VIP_pos_meanPBS = tmpDF[tmpDF.index.isin(CGE_VIP_Pos_list)]["beta"].mean()
    tmp_VIP_neg_meanPBS = tmpDF[tmpDF.index.isin(CGE_VIP_Neg_list)]["beta"].mean()
    null_diffs.append(tmp_VIP_pos_meanPBS - tmp_VIP_neg_meanPBS)
null_diffs = np.array(null_diffs)

# %%
# Plot observed difference vs null distribution
plt.figure(figsize=(8,6))
plt.hist(null_diffs, bins=30, density=True, alpha=0.5, label='Null distribution')
plt.axvline(DIFF_obs, color='red', linestyle='dashed', label='Observed difference')
plt.xlabel('VIP+ vs VIP- mean bias difference')
plt.ylabel('Density')
plt.legend()

#compute 1-sided p-value
pval = np.mean(null_diffs <= DIFF_obs)
print(f"Observed difference: {DIFF_obs:.3f}")
print(f"P-value: {pval:.3e}")

# %%
HumanCT_res_df_GeneL = pd.read_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.cluster.June10.csv", index_col=0)
results_df = pd.read_csv("../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.Supercluster.June10.csv")

# %%
results_df

# %%
CGE_VIP_Pos_list

# %%
HumanCT_res_df_GeneL_VIP = HumanCT_res_df_GeneL.loc[CGE_VIP_Pos_list, :]

# Get all interneuron clusters except the VIP+ ones
other_in_mask = HumanCT_res_df_GeneL["Supercluster"].isin([
    "LAMP5-LHX6 and Chandelier", "CGE interneuron", "MGE interneuron"
])
# Exclude VIP+ clusters
other_in_indices = HumanCT_res_df_GeneL.index[other_in_mask & ~HumanCT_res_df_GeneL.index.isin(CGE_VIP_Pos_list)]
HumanCT_res_df_GeneL_Other_IN = HumanCT_res_df_GeneL.loc[other_in_indices, :]
print(len(HumanCT_res_df_GeneL_VIP))
print(len(HumanCT_res_df_GeneL_Other_IN))

# %%
HumanCT_res_df_GeneL_VIP

# %%

EFFECT = "p_beta_perm_Log"
data = [HumanCT_res_df_GeneL_VIP[EFFECT], HumanCT_res_df_GeneL_Other_IN[EFFECT]]

# Perform Mann-Whitney U test
stat, pval = scipy.stats.mannwhitneyu(
    HumanCT_res_df_GeneL_VIP[EFFECT], 
    HumanCT_res_df_GeneL_Other_IN[EFFECT]
)

# Create boxplot with individual points
bp = plt.boxplot(data, labels=["VIP+", "VIP-"])

# Add scatter points for each group
for i, d in enumerate(data):
    x = np.random.normal(i+1, 0.04, size=len(d))
    plt.scatter(x, d, alpha=0.4, s=20)

plt.ylabel("Effect")
plt.title(f"p = {pval:.2e}")
plt.show()

# %%

# %%

# %% [markdown]
# ### LGD vs Dmis

# %%
Mut_n_IQ_conf_LGD = Mut_n_IQ_conf[Mut_n_IQ_conf["GeneEff"]!="missense"]
Mut_n_IQ_conf_Dmis = Mut_n_IQ_conf[Mut_n_IQ_conf["GeneEff"]=="missense"]
print(Mut_n_IQ_conf.shape, Mut_n_IQ_conf_LGD.shape, Mut_n_IQ_conf_Dmis.shape)

# %%
# Calculate gene level for Mut_n_IQ_conf
Genes_LGD = list(set(Mut_n_IQ_conf_LGD["Entrez"].values))
data_LGD = []
for g in Genes_LGD:
    tmp_df = Mut_n_IQ_conf_LGD[Mut_n_IQ_conf_LGD["Entrez"]==g]
    avg_IQ = tmp_df["IQ"].mean()
    row = [g, avg_IQ]
    data_LGD.append(row)
columns = ["Entrez", "IQ"]
Avg_Gene_IQ_DF_LGD = pd.DataFrame(data=data_LGD, columns=columns)

# Calculate gene level for Mut_n_IQ_Dmis
Genes_Dmis = list(set(Mut_n_IQ_conf_Dmis["Entrez"].values))
data_Dmis = []
for g in Genes_Dmis:
    tmp_df = Mut_n_IQ_conf_Dmis[Mut_n_IQ_conf_Dmis["Entrez"]==g]
    avg_IQ = tmp_df["IQ"].mean()
    row = [g, avg_IQ]
    data_Dmis.append(row)
Avg_Gene_IQ_DF_Dmis = pd.DataFrame(data=data_Dmis, columns=columns)



# %%
HumanCT_res_df_GeneL_LGD = Make_HumanCT_DF(Avg_Gene_IQ_DF_LGD, HumanCT_Z2_HCT, "../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.LGD.csv")
HumanCT_res_df_GeneL_Dmis = Make_HumanCT_DF(Avg_Gene_IQ_DF_Dmis, HumanCT_Z2_HCT, "../dat/Pheno_Bias_vs_IQ/HumanCT.GeneL.Dmis.csv")

# %%
SuperClusterBias_BoxPlot_CorrIQ(HumanCT_res_df_GeneL_LGD, flip_axis=True, figsize=(6, 8))

# %%
SuperClusterBias_BoxPlot_CorrIQ(HumanCT_res_df_GeneL_Dmis, flip_axis=True, figsize=(6, 8))

# %%
Perm_DIR_LGD = "/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT_Feb25/LGD/"
Perm_DFs_LGD = []

max_count = 100000

i = 0
for df in os.listdir(Perm_DIR_LGD):
    df = pd.read_csv(f"{Perm_DIR_LGD}/{df}", index_col=0)
    df.index = df.index.astype(int)
    Perm_DFs_LGD.append(df)
    i += 1
    if i > max_count:
        break
print(len(Perm_DFs_LGD))

Perm_DIR_Dmis = "/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT_Feb25/Dmis/"
Perm_DFs_Dmis = []

max_count = 100000

i = 0
for df in os.listdir(Perm_DIR_Dmis):
    df = pd.read_csv(f"{Perm_DIR_Dmis}/{df}", index_col=0)
    df.index = df.index.astype(int)
    Perm_DFs_Dmis.append(df)
    i += 1
    if i > max_count:
        break
print(len(Perm_DFs_Dmis))

# %%
HumanCT_res_df_GeneL_LGD

# %%
results_df_LGD = pd.DataFrame(columns=['Supercluster', 'p_rho', 'p_beta'])

for Supercluster in Neurons:
    ClusterIdx = Anno[Anno["Supercluster"]==Supercluster].index.values
    p_rho_LGD, p_beta_LGD = plot_null_suptercluster_distributions(ClusterIdx, HumanCT_res_df_GeneL_LGD, Perm_DFs_LGD, plot=False)
    
    # Add results to DataFrame
    results_df_LGD = results_df_LGD.append({
        'Supercluster': Supercluster,
        'p_rho': p_rho_LGD,
        'p_beta': p_beta_LGD
    }, ignore_index=True)

# %%
results_df_Dmis = pd.DataFrame(columns=['Supercluster', 'p_rho', 'p_beta'])

for Supercluster in Neurons:
    ClusterIdx = Anno[Anno["Supercluster"]==Supercluster].index.values
    p_rho_Dmis, p_beta_Dmis = plot_null_suptercluster_distributions(ClusterIdx, HumanCT_res_df_GeneL_Dmis, Perm_DFs_Dmis, plot=False)
    
    # Add results to DataFrame
    results_df_Dmis = results_df_Dmis.append({
        'Supercluster': Supercluster,
        'p_rho': p_rho_Dmis,
        'p_beta': p_beta_Dmis
    }, ignore_index=True)

# %%
results_df_LGD

# %%
results_df_Dmis

# %% [markdown]
# # UKBB

# %%
