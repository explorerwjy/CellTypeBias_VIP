# Compute permutation p-values for IQ-bias association for MouseCT / MouseSTR
# Example: python script_IQ_Permutation.py --permutation_type MouseCT


import sys
sys.path.insert(1, '/home/jw3514/Work/CellType_Psy/src')
from CellType_PSY import *

from statsmodels.stats.multitest import fdrcorrection, multipletests
import statsmodels.api as sm
from multiprocessing import Pool, cpu_count

import argparse

# Load the data

# Function to compute gene level IQ
def compute_gene_level_IQ(mut_level_df):
    Genes = list(set(mut_level_df["Entrez"].values))
    data = []
    for g in Genes:
        tmp_df = mut_level_df[mut_level_df["Entrez"]==g]
        avg_IQ = tmp_df["IQ"].mean()
        row = [g, avg_IQ]
        data.append(row)
    columns = ["Entrez", "IQ"]
    Avg_Gene_IQ_DF = pd.DataFrame(data=data, columns=columns)
    return Avg_Gene_IQ_DF

# Function to compute bias vs phenotype
def BiasVsPheno(MutPhenoDF, BiasMat, STR, label):
    biases = []
    IQs = []
    for i, row in MutPhenoDF.iterrows():
        if label == 'label':
            gene = row["HGNC"]
        else:
            gene = int(row["Entrez"])
        if gene in BiasMat.index.values:
            bias = BiasMat.loc[gene, STR]
            if bias == bias:  # This is checking for NaN in a non-ideal way
                IQ = row["IQ"]
                biases.append(bias)
                IQs.append(IQ)
    return np.array(biases), np.array(IQs)

# Function to perform linear fit
def linear_fit(biases, IQs, alpha=0.05):
    if len(biases) < 2:  # Add check for minimum required samples
        raise ValueError("Need at least 2 samples for linear regression")
        
    model = sm.OLS(IQs, sm.add_constant(biases))
    results = model.fit()
    
    intercept = results.params[0]
    beta = results.params[1]
    ci = results.conf_int(alpha=alpha)
    ci_low = ci[1][0]
    ci_high = ci[1][1]     
    r_value = results.rsquared
    p_value = results.pvalues[1]
    std_err = results.bse[1]
    pho, p = spearmanr(biases, IQs)
    
    return intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p

# Function to create Mouse STR DataFrame
def Make_Mouse_STR_DF(Mut_n_IQ_conf, Mouse_STR_Z2_Mat, STR_Ann, output_file, alpha=0.05):
    names, supercluster, correlation, pvalues, beta_values, beta_ci_low, beta_ci_high, intercept_values, r_value_values, p_value_values, std_err_values = [],[],[],[],[],[],[],[],[],[],[]
    for STR in Mouse_STR_Z2_Mat.columns.values:
        try:
            biases, IQs = BiasVsPheno(Mut_n_IQ_conf, Mouse_STR_Z2_Mat , STR, 'xx')
            if len(biases) < 2:
                continue
            intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p = linear_fit(biases, IQs, alpha=0.05)
            
            names.append("{}".format(STR))
            supercluster.append(STR_Ann.loc[STR, "REG"])
            correlation.append(pho)
            pvalues.append(p)
            beta_values.append(beta)
            beta_ci_low.append(ci_low)
            beta_ci_high.append(ci_high)
            intercept_values.append(intercept)
            r_value_values.append(r_value)
            p_value_values.append(p_value)
            std_err_values.append(std_err)
        except Exception as e:
            print(f"Error processing STR {STR}: {str(e)}")
            continue

    str_res_df = pd.DataFrame(data={"STR":names, "Region":supercluster, "SpearmanR":correlation, "SpearmanP":pvalues, 
                                            "beta":beta_values, "CI_low":beta_ci_low, "CI_high":beta_ci_high, "intercept":intercept_values, "r_value":r_value_values, 
                                            "p_value":p_value_values, "std_err":std_err_values})
    str_res_df = str_res_df.sort_values("SpearmanR")
    str_res_df.to_csv(output_file, compression='gzip')
    return str_res_df
    
def process_permutation_MouseSTR(idx, Mut_n_IQ_conf, Mouse_STR_Z2_Mat, STR_Ann, ASD_Circuits, seed):
    try:
        np.random.seed(seed)
        permuted_mut_df = Mut_n_IQ_conf.copy(deep=True)
        permuted_mut_df['IQ'] = np.random.permutation(permuted_mut_df['IQ'].values)
        permuted_mut_df_LGD = permuted_mut_df[permuted_mut_df["GeneEff"]!="missense"]
        permuted_mut_df_Dmis = permuted_mut_df[permuted_mut_df["GeneEff"]=="missense"]

        permuted_gene_df = compute_gene_level_IQ(permuted_mut_df)
        output_path = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT/MouseSTR.GeneL_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL = Make_Mouse_STR_DF(permuted_gene_df, Mouse_STR_Z2_Mat, STR_Ann, output_path)

        permuted_gene_df_LGD = compute_gene_level_IQ(permuted_mut_df_LGD)
        output_path_LGD = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT_LGD/MouseSTR.GeneL_LGD_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL_LGD = Make_Mouse_STR_DF(permuted_gene_df_LGD, Mouse_STR_Z2_Mat, STR_Ann, output_path_LGD)

        permuted_gene_df_Dmis = compute_gene_level_IQ(permuted_mut_df_Dmis)
        output_path_Dmis = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT_Dmis/MouseSTR.GeneL_Dmis_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL_Dmis = Make_Mouse_STR_DF(permuted_gene_df_Dmis, Mouse_STR_Z2_Mat, STR_Ann, output_path_Dmis)

        mouse_str_res_df_geneL_ASD = mouse_str_res_df_GeneL[mouse_str_res_df_GeneL["STR"].isin(ASD_Circuits)]
        mouse_str_res_df_geneL_ASD = mouse_str_res_df_geneL_ASD.set_index("STR")
        mouse_str_res_df_geneL_ASD.loc["Bed_nuclei_of_the_stria_terminalis", "Region"] = "Amygdala"
        if "Anterior_pretectal_nucleus" in mouse_str_res_df_geneL_ASD.index:
            mouse_str_res_df_geneL_ASD.drop("Anterior_pretectal_nucleus", inplace=True)

        mean_beta = mouse_str_res_df_geneL_ASD["beta"].mean()
        return mean_beta
    except Exception as e:
        print(f"Error in permutation {idx}: {str(e)}")
        return None

def MouseSTR_Permutation(Mut_n_IQ_conf=None):  # Add default None
    n_permutations = 10000
    n_threads = 20

    # Calculate bias
    str2reg = STR2Region()
    ASD_CircuitsSet = pd.read_csv(
        "/home/jw3514/Work/ASD_Circuits/notebooks/ASD.SA.Circuits.Size46.csv",
        index_col="idx")
    ASD_Circuits = ASD_CircuitsSet.loc[3, "STRs"].split(";")

    CIR_REGIONS = ['Isocortex', 'Hippocampus', 'Cortical_subplate',
         'Amygdala', 'Striatum', 'Thalamus', 'Olfactory_areas']
    CIR_REGIONS_Dict = {}
    for i in range(len(CIR_REGIONS)):
        CIR_REGIONS_Dict[CIR_REGIONS[i]] = []
    for _str in ASD_Circuits:
        for i in range(len(CIR_REGIONS)):
            if str2reg[_str] == CIR_REGIONS[i]:
                CIR_REGIONS_Dict[CIR_REGIONS[i]].append(_str)
                break
    CIR_REGIONS_Dict["Amygdala"].append("Bed_nuclei_of_the_stria_terminalis")

    Mouse_STR_Z2_Mat = pd.read_csv("../../ASD_Circuits/dat/allen-mouse-exp/AllenMouseBrain_Z2bias.csv", index_col=0)
    STR_Ann = pd.read_csv("/home/jw3514/Work/ASD_Circuits/dat/structure2region.tsv", delimiter="\t", index_col=0)

    if Mut_n_IQ_conf is None:
        Mut_n_IQ_conf = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv")

    # Compute IQ_bias association for each permutation in parallel
    seeds = np.random.randint(0, 1000000, size=n_permutations)
    with Pool(n_threads) as pool:
        mean_betas = pool.starmap(process_permutation_MouseSTR, [(idx, Mut_n_IQ_conf, Mouse_STR_Z2_Mat, STR_Ann, ASD_Circuits, seed) for idx, seed in enumerate(seeds)])

    # Filter out None values and save mean_betas to a file
    mean_betas = [beta for beta in mean_betas if beta is not None]
    mean_betas_df = pd.DataFrame(mean_betas, columns=["mean_beta"])
    mean_betas_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/mean_betas.csv", index=False)



def process_permutation_MouseCT(idx, Mut_n_IQ_conf, Mouse_CT_Z2_Mat, ASD_CTs, seed, outDIR):
    try:
        np.random.seed(seed)
        permuted_mut_df = Mut_n_IQ_conf.copy(deep=True)
        permuted_mut_df['IQ'] = np.random.permutation(permuted_mut_df['IQ'].values)
        permuted_mut_df_LGD = permuted_mut_df[permuted_mut_df["GeneEff"]!="missense"]
        permuted_mut_df_Dmis = permuted_mut_df[permuted_mut_df["GeneEff"]=="missense"]


        permuted_gene_df = compute_gene_level_IQ(permuted_mut_df)
        output_path = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT/MouseSTR.GeneL_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL = Make_Mouse_CT_DF(permuted_gene_df, Mouse_CT_Z2_Mat, ASD_CTs, output_path)

        permuted_gene_df_LGD = compute_gene_level_IQ(permuted_mut_df_LGD)
        output_path_LGD = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT_LGD/MouseSTR.GeneL_LGD_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL_LGD = Make_Mouse_CT_DF(permuted_gene_df_LGD, Mouse_CT_Z2_Mat, ASD_CTs, output_path_LGD)    

        permuted_gene_df_Dmis = compute_gene_level_IQ(permuted_mut_df_Dmis)
        output_path_Dmis = f"/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT_Dmis/MouseSTR.GeneL_Dmis_perm_{idx}.csv.gz"
        mouse_str_res_df_GeneL_Dmis = Make_Mouse_CT_DF(permuted_gene_df_Dmis, Mouse_CT_Z2_Mat, ASD_CTs, output_path_Dmis)


        # permuted_gene_df = compute_gene_level_IQ(permuted_mut_df)
        # output_path = f"{outDIR}/MouseCT.GeneL_perm_{idx}.csv.gz"
        # mouse_str_res_df_GeneL = Make_Mouse_CT_DF(permuted_gene_df, Mouse_CT_Z2_Mat, ASD_CTs, output_path)

        mean_beta = mouse_str_res_df_GeneL["beta"].mean()
        return mean_beta
    except Exception as e:
        print(f"Error in permutation {idx}: {str(e)}")
        return None

def Make_Mouse_CT_DF(Mut_n_IQ_conf, Mouse_CT_Z2_Mat, ClusterAnn, output_file, alpha=0.05):
    names, supercluster, correlation, pvalues, beta_values, beta_ci_low, beta_ci_high, intercept_values, r_value_values, p_value_values, std_err_values = [],[],[],[],[],[],[],[],[],[],[]
    for CT in ClusterAnn.index.values:
        try:
            if CT not in Mouse_CT_Z2_Mat.columns.values:    
                continue
            biases, IQs = BiasVsPheno(Mut_n_IQ_conf, Mouse_CT_Z2_Mat , CT, 'xx')
            if len(biases) < 2:
                continue
            intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p = linear_fit(biases, IQs, alpha=0.05)
            
            names.append("{}".format(CT))
            supercluster.append(ClusterAnn.loc[CT, "class_label"])
            correlation.append(pho)
            pvalues.append(p)
            beta_values.append(beta)
            beta_ci_low.append(ci_low)
            beta_ci_high.append(ci_high)
            intercept_values.append(intercept)
            r_value_values.append(r_value)
            p_value_values.append(p_value)
            std_err_values.append(std_err)
        except Exception as e:
            print(f"Error processing CT {CT}: {str(e)}")
            continue

    str_res_df = pd.DataFrame(data={"cluster_id":names, "Class":supercluster, "SpearmanR":correlation, "SpearmanP":pvalues, 
                                            "beta":beta_values, "CI_low":beta_ci_low, "CI_high":beta_ci_high, "intercept":intercept_values, "r_value":r_value_values, 
                                            "p_value":p_value_values, "std_err":std_err_values})
    str_res_df = str_res_df.sort_values("SpearmanR")
    str_res_df.to_csv(output_file, compression='gzip')
    return str_res_df

def MouseCT_Permutation(Mut_n_IQ_conf, outDIR):
    n_permutations = 10000
    n_threads = 20

    Mouse_CT_Z2_Mat = pd.read_csv("../AllenBrainCellAtlas/dat/SC_UMI_Mats/MouseCT.Z2Mat.ASD.csv", index_col=0)
    ASD_CTs = pd.read_csv("/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/MouseCT.Selected_ASD_CTs.csv", index_col=0)
    
    if Mut_n_IQ_conf is None:
        Mut_n_IQ_conf = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv")
    else:
        Mut_n_IQ_conf = pd.read_csv(Mut_n_IQ_conf)

    seeds = np.random.randint(0, 1000000, size=n_permutations)
    with Pool(n_threads) as pool:
        mean_betas = pool.starmap(process_permutation_MouseCT, [(idx, Mut_n_IQ_conf, Mouse_CT_Z2_Mat, ASD_CTs, seed, outDIR) for idx, seed in enumerate(seeds)])

    # Filter out None values and save mean_betas to a file
    mean_betas = [beta for beta in mean_betas if beta is not None]
    mean_betas_df = pd.DataFrame(mean_betas, columns=["mean_beta"])
    mean_betas_df.to_csv("/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/mean_betas.ct.csv", index=False)

def Make_Human_CT_DF(Mut_n_IQ_conf, Human_CT_Z2_Mat, output_file, alpha=0.05):
    names, supercluster, correlation, pvalues, beta_values, beta_ci_low, beta_ci_high, intercept_values, r_value_values, p_value_values, std_err_values = [],[],[],[],[],[],[],[],[],[],[]
    for CT in Anno.index.values:
        try:
            if CT not in Human_CT_Z2_Mat.columns.values:    
                continue
            biases, IQs = BiasVsPheno(Mut_n_IQ_conf, Human_CT_Z2_Mat , CT, 'xx')
            if len(biases) < 2:
                continue
            intercept, beta, ci_low, ci_high, r_value, p_value, std_err, pho, p = linear_fit(biases, IQs, alpha=0.05)
            
            names.append("{}".format(CT))
            supercluster.append(Anno.loc[CT, "Supercluster"])
            correlation.append(pho)
            pvalues.append(p)
            beta_values.append(beta)
            beta_ci_low.append(ci_low)
            beta_ci_high.append(ci_high)
            intercept_values.append(intercept)
            r_value_values.append(r_value)
            p_value_values.append(p_value)
            std_err_values.append(std_err)
        except Exception as e:
            print(f"Error processing CT {CT}: {str(e)}")
            continue

    str_res_df = pd.DataFrame(data={"cluster_id":names, "Supercluster":supercluster, "SpearmanR":correlation, "SpearmanP":pvalues, 
                                            "beta":beta_values, "CI_low":beta_ci_low, "CI_high":beta_ci_high, "intercept":intercept_values, "r_value":r_value_values, 
                                            "p_value":p_value_values, "std_err":std_err_values})
    str_res_df = str_res_df.sort_values("SpearmanR")
    str_res_df.to_csv(output_file, compression='gzip')
    return str_res_df

def process_permutation_HumanCT(idx, Mut_n_IQ_conf, Human_CT_Z2_Mat, seed, outDIR):
    try:
        np.random.seed(seed)
        permuted_mut_df = Mut_n_IQ_conf.copy(deep=True)
        permuted_mut_df['IQ'] = np.random.permutation(permuted_mut_df['IQ'].values)
        permuted_mut_df_LGD = permuted_mut_df[permuted_mut_df["GeneEff"]!="missense"]
        permuted_mut_df_Dmis = permuted_mut_df[permuted_mut_df["GeneEff"]=="missense"]


        permuted_gene_df = compute_gene_level_IQ(permuted_mut_df)
        output_path = f"{outDIR}/ALL/HumanCT.GeneL_perm_{idx}.csv.gz"
        res_df_GeneL = Make_Human_CT_DF(permuted_gene_df, Human_CT_Z2_Mat, output_path)

        permuted_gene_df_LGD = compute_gene_level_IQ(permuted_mut_df_LGD)
        output_path_LGD = f"{outDIR}/LGD/HumanCT.GeneL_LGD_perm_{idx}.csv.gz"
        res_df_GeneL_LGD = Make_Human_CT_DF(permuted_gene_df_LGD, Human_CT_Z2_Mat, output_path_LGD)    

        permuted_gene_df_Dmis = compute_gene_level_IQ(permuted_mut_df_Dmis)
        output_path_Dmis = f"{outDIR}/Dmis/HumanCT.GeneL_Dmis_perm_{idx}.csv.gz"
        res_df_GeneL_Dmis = Make_Human_CT_DF(permuted_gene_df_Dmis, Human_CT_Z2_Mat, output_path_Dmis)

    except Exception as e:
        print(f"Error in permutation {idx}: {str(e)}")
        return None

def HumanCT_Permutation(Mut_n_IQ_conf, outDIR, n_permutations):
    n_threads = 20

    Human_CT_Z2_Mat = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Jan21.csv", index_col=0)
    Human_CT_Z2_Mat.columns = Human_CT_Z2_Mat.columns.astype(int)
    
    if Mut_n_IQ_conf is None:
        Mut_n_IQ_conf = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv")
    else:
        Mut_n_IQ_conf = pd.read_csv(Mut_n_IQ_conf)
    try:
        os.makedirs(f"{outDIR}/ALL", exist_ok=True)
        os.makedirs(f"{outDIR}/LGD", exist_ok=True) 
        os.makedirs(f"{outDIR}/Dmis", exist_ok=True)
    except Exception as e:
        print(f"Error creating directories: {str(e)}")
        raise
          
    seeds = np.random.randint(0, 1000000, size=n_permutations)
    with Pool(n_threads) as pool:
        pool.starmap(process_permutation_HumanCT, [(idx, Mut_n_IQ_conf, Human_CT_Z2_Mat, seed, outDIR) for idx, seed in enumerate(seeds)])


def main():
    parser = argparse.ArgumentParser(description="Run MouseSTR or MouseCT permutation")
    parser.add_argument("--permutation_type", type=str, choices=["MouseSTR", "MouseCT", "HumanCT"], default="MouseSTR", help="Type of permutation to run: MouseSTR or MouseCT or HumanCT")
    parser.add_argument("--Mut_n_IQ", type=str, help="File path to Mut_n_IQ_conf.csv")
    parser.add_argument("--n_permutations", type=int, default=10000, help="Number of permutations to run")
    parser.add_argument("--outDIR", type=str, help="Output directory")
    args = parser.parse_args()

    if args.permutation_type == "MouseSTR":
        MouseSTR_Permutation(args.Mut_n_IQ)  # Pass Mut_n_IQ argument
    elif args.permutation_type == "MouseCT":
        MouseCT_Permutation(args.Mut_n_IQ, args.outDIR)
    elif args.permutation_type == "HumanCT":
        HumanCT_Permutation(args.Mut_n_IQ, args.outDIR, args.n_permutations)
    else:
        raise ValueError("Invalid permutation type. Choose 'MouseSTR' or 'MouseCT' or 'HumanCT'.")

if __name__ == "__main__":
    main()

# Example: python script_IQ_Permutation.py --permutation_type MouseCT --Mut_n_IQ /home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv --outDIR /home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/MouseCT
# Example: python script_IQ_Permutation.py --permutation_type HumanCT --Mut_n_IQ /home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/Mut_n_IQ_conf.csv --outDIR /home/jw3514/Work/CellType_Psy/dat/Pheno_Bias_vs_IQ/IQ_Permuts/HumanCT
