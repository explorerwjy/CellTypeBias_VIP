import sys
sys.path.insert(1, '../src/')
from ASD_Circuits import *
import argparse

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help = 'Match Feature Table')
    parser.add_argument('-o', '--output', required=True, type=str, help = 'Output Dir')
    args = parser.parse_args()
    return args

def Match_a_Gene(Gene_ID, MatchFeature, sample_size=10000, outfil = None, kernel="uniform", 
                 interval=0.05):
    gene_rank = MatchFeature.loc[Gene_ID, "Rank"]
    quantile = MatchFeature.loc[Gene_ID, "quantile"]
    min_quantile = max(0, quantile - interval)
    max_quantile = min(1, quantile + interval)
    Interval_genes = MatchFeature[(MatchFeature["quantile"]>=min_quantile)&
                               (MatchFeature["quantile"]<=max_quantile)&
                               (MatchFeature.index != Gene_ID)].index.values
    if kernel == "uniform":
        match_genes = np.random.choice(Interval_genes, size=sample_size, replace=True)
    elif kernel == "tricubic":
        print("To Be Implemented")
    if outfil:
        List2Fil(match_genes, outfil)

def main():
    args = GetOptions()
    MatchFeature = pd.read_csv(args.input, index_col=0)
    OutDir = args.output
    for g, row in MatchFeature.iterrows():
        Match_a_Gene(g, MatchFeature, outfil="{}/{}.csv".format(OutDir, g))

if __name__=='__main__':
    main()

