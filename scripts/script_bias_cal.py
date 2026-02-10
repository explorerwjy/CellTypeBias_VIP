import sys
import os

# Get project directory from script location (works regardless of where script is called from)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ProjDIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(1, os.path.join(ProjDIR, 'src'))

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

def AddPvalue_Supercluster(df, Null_Bias_values):
    """
    Calculate supercluster-level bias and p-values.

    Parameters:
    -----------
    df : DataFrame
        Cluster-level bias results with 'EFFECT' and 'Supercluster' columns
    Null_Bias_values : DataFrame
        Null distribution matrix (clusters x simulations)

    Returns:
    --------
    DataFrame with supercluster-level results including p-values and FDR correction
    """
    from statsmodels.stats.multitest import multipletests

    # Check that df has Supercluster column
    if 'Supercluster' not in df.columns:
        raise ValueError("Input dataframe must have 'Supercluster' column. Run AnnotateCTDat first.")

    # Group by supercluster and calculate mean EFFECT (observed)
    supercluster_obs = df.groupby('Supercluster')['EFFECT'].mean()

    # For null distribution, group clusters by supercluster and calculate mean for each simulation
    # First, create mapping of cluster_id to supercluster
    cluster_to_supercluster = df['Supercluster'].to_dict()

    # Get list of superclusters
    superclusters = supercluster_obs.index.tolist()
    n_sims = Null_Bias_values.shape[1]

    # Initialize null matrix for superclusters (n_sims x n_superclusters)
    null_supercluster_matrix = np.zeros((n_sims, len(superclusters)))

    # For each simulation, calculate mean EFFECT per supercluster
    for sim_idx in range(n_sims):
        sim_col = Null_Bias_values.columns[sim_idx]

        # Get null bias for this simulation
        sim_bias = Null_Bias_values[sim_col]

        # Add supercluster annotation
        sim_df = pd.DataFrame({'EFFECT': sim_bias})
        sim_df['Supercluster'] = sim_df.index.map(cluster_to_supercluster)

        # Calculate mean per supercluster
        supercluster_means = sim_df.groupby('Supercluster')['EFFECT'].mean()

        # Store in matrix (align with superclusters order)
        for sc_idx, sc_name in enumerate(superclusters):
            if sc_name in supercluster_means.index:
                null_supercluster_matrix[sim_idx, sc_idx] = supercluster_means.loc[sc_name]
            else:
                null_supercluster_matrix[sim_idx, sc_idx] = np.nan

    # Calculate p-values using vectorized function
    observed_vals = supercluster_obs.values
    z_scores, p_values, obs_adjs = GetPermutationP_vectorized(null_supercluster_matrix, observed_vals)

    # Create results dataframe
    results = pd.DataFrame({
        'Supercluster': superclusters,
        'EFFECT': observed_vals,
        'P-value': p_values,
        'Z-score': z_scores,
        'EFFECT_adj': obs_adjs
    })
    results.set_index('Supercluster', inplace=True)

    # Calculate FDR-corrected q-values
    _, q_values = multipletests(results["P-value"].values, alpha=0.1, method="fdr_i")[0:2]
    results["q-value"] = q_values
    results["-logP"] = -np.log10(results["P-value"])

    # Add number of clusters per supercluster
    cluster_counts = df.groupby('Supercluster').size()
    results['n_clusters'] = results.index.map(cluster_counts)

    return results

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate bias for a gene set using provided weights and expression matrix.")
    parser.add_argument('--SpecMat', required=True, help='Path to expression matrix (CSV)')
    parser.add_argument('--gw', required=True, help='Path to gene weights file (CSV)')
    parser.add_argument('--Bias_Out', required=True, help='Output file for cluster-level bias results (CSV)')
    parser.add_argument('--Bias_Null', required=True, help='Path to null bias distribution (CSV)')
    parser.add_argument('--Bias_Out_Supercluster', default=None,
                        help='Output file for supercluster-level results (CSV). If not provided, derived from --Bias_Out')
    return parser.parse_args()

def main():
    args = parse_args()
    SpecMat = pd.read_csv(args.SpecMat, index_col=0)
    GW = Fil2Dict(args.gw)
    Bias = HumanCT_AvgZ_Weighted(SpecMat, GW) # compute bias
    Bias = AnnotateCTDat(Bias, Anno)

    Null_Bias_values = pd.read_csv(args.Bias_Null, index_col=0)

    # Cluster-level p-values
    Bias_add_P = AddPvalue_vec(Bias, Null_Bias_values)
    Bias_add_P.to_csv(args.Bias_Out)
    print(f"Cluster-level results saved to: {args.Bias_Out}")

    # Supercluster-level p-values
    Bias_Supercluster = AddPvalue_Supercluster(Bias, Null_Bias_values)

    # Use explicit path if provided, otherwise derive from Bias_Out
    if args.Bias_Out_Supercluster:
        supercluster_out = args.Bias_Out_Supercluster
    elif args.Bias_Out.endswith('_bias_addP.csv'):
        supercluster_out = args.Bias_Out.replace('_bias_addP.csv', '_bias_addP_supercluster.csv')
    elif args.Bias_Out.endswith('.csv'):
        supercluster_out = args.Bias_Out.replace('.csv', '_supercluster.csv')
    else:
        supercluster_out = args.Bias_Out + '_supercluster.csv'

    Bias_Supercluster.to_csv(supercluster_out)
    print(f"Supercluster-level results saved to: {supercluster_out}")
    

if __name__ == "__main__":
    main()
