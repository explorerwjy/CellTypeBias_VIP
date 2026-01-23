# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_run_ctrl_sim.py
# Run control simutations for expression biases (Cell Type; Structures)
# ========================================================================================================

import argparse
import sys
import os
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
sys.path.insert(1, '/home/jw3514/Work/UNIMED/src')
from CellType_PSY import *
from UNIMED import *
import time
import json
import pickle
import multiprocessing
from multiprocessing import Pool

def SaveDict(Dict, fname):
    with open(fname, 'wb') as hand:
        pickle.dump(Dict, hand)
    return

def LoadDict(fname):
    with open(fname, 'rb') as hand:
        b = pickle.load(hand)
        return b

def GeneWeightSimulation(GW):
    return

def parallel_SubSampleSibling():
    return

def SubSampleSibling(WeightDF, outdir, GeneProb, n_sims=10000):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    SibWeightDF = pd.read_csv("/home/jw3514/Work/ASD_Circuits/dat/Unionize_bias/sibling_weights_LGD_Dmis.csv", header=None)
    with open('/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ClusterV3_Z2_Mat_Genes.pkl', 'rb') as file:
        MouseCT_Genes = pickle.load(file)
    SibGenes = SibWeightDF[0].values
    SibGenes = [g for g in SibGenes if g in MouseCT_Genes]
    WeightDF = pd.read_csv(WeightDF, header=None)
    #Gene_Weights = [1] * len(WeightDF[1].values)
    Gene_Weights = WeightDF[1].values

    # Adjust Prob
    if GeneProb != None:
        Gene2Prob = pd.read_csv(GeneProb, index_col=0)
        SibGenes = [g for g in Gene2Prob.index.values if g in SibGenes]
        Gene2Prob = Gene2Prob.loc[SibGenes, :]
        probs = Gene2Prob["Prob"].values
        total = np.sum(probs)
        probs = probs/total
        probs[-1] = 1 - np.sum(probs[:-1])
        Gene2Prob["Prob"] = probs

        for i in range(n_sims):
            Genes = np.random.choice(Gene2Prob.index.values, size=len(Gene_Weights), p=Gene2Prob["Prob"].values, replace=False)
            tmp_dict = dict(zip(Genes, Gene_Weights))
            #print(len(Genes), len(Gene_Weights), len(tmp_dict))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))
    else:
        for i in range(n_sims):
            Genes = np.random.choice(SibGenes, size=len(Gene_Weights))
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))

def RandomGenes(WeightDF, outdir, GeneProb, n_sims=10000):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    GeneSetDF = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCT.cluster.Exp.Var.CV.csv", index_col=0)
    ALL_Genes = GeneSetDF.index.values
    HighCV_Genes = GeneSetDF[GeneSetDF["CV"]>0.1].index.values
    WeightDF = pd.read_csv(WeightDF, header=None)
    #Gene_Weights = [1] * len(WeightDF[1].values)
    Gene_Weights = WeightDF[1].values

    print(len(Gene_Weights))
    # Adjust Prob
    if GeneProb != None:
        Gene2Prob = pd.read_csv(GeneProb, index_col=0)
        SibGenes = [g for g in Gene2Prob.index.values if g in SibGenes]
        Gene2Prob = Gene2Prob.loc[SibGenes, :]
        probs = Gene2Prob["Prob"].values
        total = np.sum(probs)
        probs = probs/total
        probs[-1] = 1 - np.sum(probs[:-1])
        Gene2Prob["Prob"] = probs

        for i in range(n_sims):
            Genes = np.random.choice(Gene2Prob.index.values, size=len(Gene_Weights), p=Gene2Prob["Prob"].values)
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))
    else:
        for i in range(n_sims):
            #Genes = np.random.choice(ALL_Genes, size=len(Gene_Weights))
            Genes = np.random.choice(HighCV_Genes, size=len(Gene_Weights))
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))


def ExpressionMatch_RandGenes(WeightDF, outdir, ExpMatchDir, n_sims=10000):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir+"/weighted/"):
        os.makedirs(outdir+"/weighted/")
    if not os.path.exists(outdir+"/uniform/"):
        os.makedirs(outdir+"/uniform/")

    WeightDF = pd.read_csv(WeightDF, header=None)
    Gene_Weights = WeightDF[1].values
#
    # Load Exp Match genes for each gene
    MatchGenes = []
    for g in WeightDF[0].values:
        match_genes = loadgenelist(ExpMatchDir+"/{}.csv".format(g), toint=True)
        MatchGenes.append(match_genes)
    MatchGenes = np.array(MatchGenes)
    print(MatchGenes.shape)

    for i in range(n_sims):
        rand_genes = MatchGenes[:,i]
        rand_gw = dict(zip(rand_genes, Gene_Weights))
        rand_uniform = dict(zip(rand_genes, np.ones(len(rand_genes))))
        Dict2Fil(rand_gw, "{}/weighted/cont.gw.{}.csv".format(outdir, i))
        Dict2Fil(rand_uniform, "{}/uniform/cont.gw.{}.csv".format(outdir, i))


###########################################################################
## Human Cell Type CTRL Generation
###########################################################################
def AdjustClusterMean(BiasDF, HumanCT_Z2_HCT_ColMean):
    for i, row in BiasDF.iterrows():
        BiasDF.loc[i,"EFFECT"] = row["EFFECT"] - HumanCT_Z2_HCT_ColMean[i]
    return BiasDF
    
def BiasCal_SingleCtrl_HumanCT(GW_Fil, SpecMat, outDir):
    GW = Fil2Dict(GW_Fil)
    idx = GW_Fil.split("/")[-1].split(".")[2]

    Ctrl_BiasDF = HumanCT_AvgZ_Weighted(SpecMat, GW)

    #SpecMat_ColMean = SpecMat.mean(axis=0)
    #Ctrl_BiasDF = AdjustClusterMean(Ctrl_BiasDF, SpecMat_ColMean)
    #Ctrl_BiasDF = AnnotateCTDat(Ctrl_BiasDF, Anno)
    Ctrl_BiasDF.to_csv("{}/cont.bias.{}.csv.gz".format(outDir, idx), compression="gzip")
    return

def CtrlBiasCal_HumanCT(GW_Dir, SpecMat, outDir, n_processes=20): 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    Ctrl_GW_Fils = []
    for root, dirs, file_names in os.walk(GW_Dir):
        for file_name in file_names:
            Ctrl_GW_Fils.append(os.path.join(root, file_name))
    #print(Ctrl_GW_Fils)
    #print(len(Ctrl_GW_Fils))
    SpecMat = pd.read_csv(SpecMat, index_col=0)
    #SpecMat = SpecMat * 461
    SpecMat.columns = SpecMat.columns.astype(int)
    #max_Z, min_Z = 3, -3
    #SpecMat = SpecMat.clip(upper=max_Z, lower=min_Z)
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(BiasCal_SingleCtrl_HumanCT, [(GW_Fil, SpecMat, outDir) for GW_Fil in Ctrl_GW_Fils])
    pool.close()
    pool.join()

###########################################################################
## Mouse Cell Type CTRL Generation
###########################################################################
def BiasCal_SingleCtrl_MouseCT(GW_Fil, SpecMat, outDir):
    GW = Fil2Dict(GW_Fil)
    idx = GW_Fil.split("/")[-1].split(".")[2]
    Ctrl_BiasDF = ABC_AvgCTZ_Weighted(SpecMat, GW, csv_fil="{}/cont.bias.{}.csv".format(outDir, idx))
    return

def BiasCal_SingleCtrl_MouseCT_DN(GW_Fil, SpecMat, outDir, ISH_SC_CorrDF):
    GW = Fil2Dict(GW_Fil)
    GW_adj = {}
    for k,v in GW.items():
        if k in ISH_SC_CorrDF.index.values:
            GW_adj[k] = v * (ISH_SC_CorrDF.loc[k, "V2_V3_CT_Corr"]**2)
    idx = GW_Fil.split("/")[-1].split(".")[2]
    Ctrl_BiasDF = ABC_AvgCTZ_Weighted(SpecMat, GW_adj, csv_fil="{}/cont.bias.{}.csv".format(outDir, idx))
    return

def CtrlBiasCal_MouseCT(GW_Dir, SpecMat, outDir, n_processes=20, DN=True): 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    Ctrl_GW_Fils = []
    for root, dirs, file_names in os.walk(GW_Dir):
        for file_name in file_names:
            Ctrl_GW_Fils.append(os.path.join(root, file_name))
    SpecMat = pd.read_csv(SpecMat, index_col=0)
    pool = multiprocessing.Pool(processes=n_processes)
    print(DN)
    if DN:
        print("Use Denoise")
        ISH_SC_CorrDF = pd.read_csv("/home/jw3514/Work/CellType_Psy/AllenBrainCellAtlas/dat/ISH_MERFISH_Gene_CorssSTR_Corr.v2.csv", index_col=0)
        results = pool.starmap(BiasCal_SingleCtrl_MouseCT_DN, [(GW_Fil, SpecMat, outDir, ISH_SC_CorrDF) for GW_Fil in Ctrl_GW_Fils])
    else:
        print("Dont use Denoise")
        results = pool.starmap(BiasCal_SingleCtrl_MouseCT, [(GW_Fil, SpecMat, outDir) for GW_Fil in Ctrl_GW_Fils])

    pool.close()
    pool.join()

###########################################################################
## Mouse Structure CTRL Generation
###########################################################################
def BiasCal_SingleCtrl(GW_Fil, SpecMat, outDir):
    GW = Fil2Dict(GW_Fil)
    idx = GW_Fil.split("/")[-1].split(".")[2]
    Ctrl_BiasDF = AvgSTRZ_Weighted(SpecMat, GW, csv_fil="{}/cont.bias.{}.csv".format(outDir, idx))
    return

def CtrlBiasCal(GW_Dir, SpecMat, outDir, n_processes=20): 
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    Ctrl_GW_Fils = []
    for root, dirs, file_names in os.walk(GW_Dir):
        for file_name in file_names:
            Ctrl_GW_Fils.append(os.path.join(root, file_name))
    #print(Ctrl_GW_Fils)
    #print(len(Ctrl_GW_Fils))
    SpecMat = pd.read_csv(SpecMat, index_col=0)
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(BiasCal_SingleCtrl, [(GW_Fil, SpecMat, outDir) for GW_Fil in Ctrl_GW_Fils])
    pool.close()
    pool.join()

###########################################################################
## Args and Main Functions
###########################################################################
def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=str, help='Mode of program 1: gene weight generateion; 2. bias calculaton')
    parser.add_argument('-o', '--outdir', type=str, help='Output directory')
    parser.add_argument('-w', '--WeightDF', type=str, help='Weight DF for control geneset')
    parser.add_argument('-p', '--GeneProb', default=None, help='GeneProb Filname or None if dont use')
    parser.add_argument('--n_sims', type=int, default=10000, help='Number of ctrl simulations')

    parser.add_argument('--GW_Dir', type=str, help="dirctory of ctrl gene weights")
    parser.add_argument('--SpecMat', type=str, help="Filename of bias matrix")
    parser.add_argument('--MatchDir', type=str, help="Dir of exp match genes")
    parser.add_argument('--n_processes', type=int, default=20, help="Filename of bias matrix")
    parser.add_argument('--DN', type=bool, help="Use Gene Denoise or not")
    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    mode = args.mode
    if mode == 'gw':
        WeightDF = args.WeightDF
        outdir = args.outdir
        GeneProb = args.GeneProb
        n_sims = args.n_sims
        #SubSampleSibling(WeightDF, outdir, GeneProb)
        #RandomGenes(WeightDF, outdir, GeneProb, n_sims)
        ExpressionMatch_RandGenes(WeightDF, outdir, args.MatchDir, n_sims)
    if mode == 'bias':
        GW_Dir = args.GW_Dir
        SpecMat = args.SpecMat
        outDir = args.outdir
        n_processes = args.n_processes
        CtrlBiasCal(GW_Dir, SpecMat, outDir, n_processes)
    if mode == 'human_ct_bias':
        GW_Dir = args.GW_Dir
        SpecMat = args.SpecMat
        outDir = args.outdir
        n_processes = args.n_processes
        CtrlBiasCal_HumanCT(GW_Dir, SpecMat, outDir, n_processes)
    if mode == "mouse_ct_bias":
        DN=args.DN
        print(DN)
        #DN = 0
        GW_Dir = args.GW_Dir
        SpecMat = args.SpecMat
        outDir = args.outdir
        n_processes = args.n_processes
        CtrlBiasCal_MouseCT(GW_Dir, SpecMat, outDir, n_processes, DN)

    return

if __name__ == '__main__':
    main()
