# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_run_ctrl_sim.v3.py
# Run control simutations for expression biases (Cell Type; Structures)
# ========================================================================================================

import argparse
import sys
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *



def RandomGenes(ExpMat, WeightDF, outfile, GeneProb, n_sims=10000):
    import pandas as pd
    import numpy as np
    import os

    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    print(len(Gene_Weights))

    # Prepare simulation matrix: rows = genes, columns = [weight, sim0, sim1, ...]
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)

    # Adjust Prob
    if GeneProb is not None:
        Gene2Prob = pd.read_csv(GeneProb, index_col=0)
        # Only keep genes in Gene2Prob that are also in valid_genes
        filtered_genes = [g for g in Gene2Prob.index.values if g in valid_genes]
        Gene2Prob = Gene2Prob.loc[filtered_genes, :]
        probs = Gene2Prob["Prob"].values
        total = np.sum(probs)
        probs = probs / total
        # Ensure sum to 1
        if len(probs) > 1:
            probs[-1] = 1 - np.sum(probs[:-1])
        Gene2Prob["Prob"] = probs

        gene_pool = Gene2Prob.index.values
        gene_probs = Gene2Prob["Prob"].values

        for i in range(n_sims):
            Genes = np.random.choice(gene_pool, size=len(Gene_Weights), p=gene_probs, replace=False)
            sim_matrix[:, i] = Genes
    else:
        for i in range(n_sims):
            Genes = np.random.choice(valid_genes, size=len(Gene_Weights), replace=False)
            sim_matrix[:, i] = Genes

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} simulations to {outfile}")


###########################################################################
## Args and Main Functions
###########################################################################
def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', type=str, help='Output file')
    #parser.add_argument('-o', '--outdir', type=str, help='Output directory')
    parser.add_argument('-w', '--WeightDF', type=str, help='Weight DF for control geneset')
    parser.add_argument('-p', '--GeneProb', default=None, help='GeneProb Filname or None if dont use')
    parser.add_argument('--n_sims', type=int, default=10000, help='Number of ctrl simulations')
    parser.add_argument('--GW_Dir', type=str, help="dirctory of ctrl gene weights")
    parser.add_argument('--SpecMat', type=str, help="Filename of bias matrix")

    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    SpecMat = args.SpecMat
    WeightDF = args.WeightDF
    outfile = args.outfile
    GeneProb = args.GeneProb
    n_sims = args.n_sims
    RandomGenes(SpecMat, WeightDF, outfile, GeneProb, n_sims)


    return

if __name__ == '__main__':
    main()
