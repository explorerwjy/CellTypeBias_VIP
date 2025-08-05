import sys
import os
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *
import argparse
import pandas as pd



def ReadNullBias(Bias_Null):
    # Bias_Null is a DF, index are cell type (cluster_id), columns are bias from a simulation of random genes 
    Bias_Null_DF = pd.read_csv(Bias_Null, index_col=0)
    return Bias_Null_DF

def AddPvalue_vec(df, Null_Bias_values, method='vectorized', **kwargs):
#def AddPvalue_optimized(df, control_dfs, method='vectorized', **kwargs):
    """Optimized version of AddPvalue with multiple acceleration strategies."""
    # Convert to numpy arrays for faster processing
    observed_vals = df["EFFECT"].values
    cell_type_ids = df.index.tolist()
    
    # Strategy 1: Pure vectorization
    null_matrix = Null_Bias_values.loc[cell_type_ids].values.T
    #print(null_matrix.shape)
    z_scores, p_values, obs_adjs = GetPermutationP_vectorized(null_matrix, observed_vals)
        
    # Add results to dataframe
    df = df.copy()
    df["P-value"] = p_values
    df["Z-score"] = z_scores  
    df["EFFECT_adj"] = obs_adjs

    # Calculate FDR-corrected q-values
    from statsmodels.stats.multitest import multipletests
    _, q_values = multipletests(df["P-value"].values, alpha=0.1, method="fdr_i")[0:2]
    df["q-value"] = q_values
    df["-logP"] = -np.log10(df["P-value"])
    
    return df

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate bias for a gene set using provided weights and expression matrix.")
    parser.add_argument('--SpecMat', required=True, help='Path to expression matrix (CSV)')
    parser.add_argument('--gw', required=True, help='Path to gene weights file (CSV)')
    #parser.add_argument('--geneset', required=True, help='Name of the gene set')
    parser.add_argument('--Bias_Out', required=True, help='Output file for bias results (TSV)')
    parser.add_argument('--Bias_Null', required=True, help='Output file for bias results (TSV)')
    return parser.parse_args()

def main():
    args = parse_args()
    SpecMat = pd.read_csv(args.SpecMat, index_col=0)
    GW = Fil2Dict(args.gw)
    Bias = HumanCT_AvgZ_Weighted(SpecMat, GW) # compute bias
    Bias = AnnotateCTDat(Bias, Anno)

    Null_Bias_values = pd.read_csv(args.Bias_Null, index_col=0)

    Bias_add_P = AddPvalue_vec(Bias, Null_Bias_values) # 
    Bias_add_P.to_csv(args.Bias_Out)
    

if __name__ == "__main__":
    main()
