# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_run_ctrl_sim.py
# Run control simutations for expression biases (Cell Type; Structures)
# ========================================================================================================

import argparse
import sys
sys.path.insert(1, '/home/jw3514/Work/CellType_Psy/src')
from CellType_PSY import *
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
    SibGenes = SibWeightDF[0].values
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
            Genes = np.random.choice(Gene2Prob.index.values, size=len(Gene_Weights), p=Gene2Prob["Prob"].values)
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))
    else:
        for i in range(n_sims):
            Genes = np.random.choice(SibGenes, size=len(Gene_Weights))
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))

def RandomGenes(WeightDF, outdir, GeneProb, n_sims=10000):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    GeneSetDF = pd.read_csv("/home/jw3514/Work/ASD_Circuits/dat/allen-mouse-exp/AllenMouseBrain_Z2bias.csv", header=None)
    Genes = Genes.values
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
            Genes = np.random.choice(Gene2Prob.index.values, size=len(Gene_Weights), p=Gene2Prob["Prob"].values)
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))
    else:
        for i in range(n_sims):
            Genes = np.random.choice(Genes, size=len(Gene_Weights))
            tmp_dict = dict(zip(Genes, Gene_Weights))
            Dict2Fil(tmp_dict, "{}/cont.gw.{}.csv".format(outdir, i))

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

def process_Similarity_ASD_SCZ(idx, CT_Z2_MAT_HC, ASD_GeneLofZ, SCZ_GeneLofZ):
    DF1 = ASD_GeneLofZ.sample(frac=1, random_state=idx)
    DF2 = SCZ_GeneLofZ.sample(frac=1, random_state=idx)
    yy = []
    for i in range(0, 31, 1):
        tmp_ASD_GW = dict(zip(DF1["Entrez"].values[i:], DF1["GW"].values[i:]))
        tmp_SCZ_GW = dict(zip(DF2["Entrez"].values[i:], DF2["GW"].values[i:]))
        tmp_ASD_Bias = AvgCTZ_Weighted(CT_Z2_MAT_HC, tmp_ASD_GW, Method = 1)
        tmp_ASD_Bias = AnnotateCTDat(tmp_ASD_Bias, Anno)
        tmp_SCZ_Bias = AvgCTZ_Weighted(CT_Z2_MAT_HC, tmp_SCZ_GW, Method = 1)
        tmp_SCZ_Bias = AnnotateCTDat(tmp_SCZ_Bias, Anno)
        r,p = GetSingeCellBiasCorr(tmp_ASD_Bias, tmp_SCZ_Bias)
        yy.append(r)
    #print(yy)
    return yy 

def Similarity_ASD_SCZ(n_processes=20):
    HCT_Z2_MAT_HCT = pd.read_csv("../dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv", index_col=0)
    max_Z, min_Z = 3, -3
    HCT_Z2_MAT_HCT = HCT_Z2_MAT_HCT.clip(upper=max_Z, lower=min_Z)
    #CT_Z2_MAT_HC = pd.read_csv("../dat/HumanCellType.AllCell.HCT.Z2bias.entrez.csv", index_col=0)
    CT_Z2_MAT_HC = HCT_Z2_MAT_HCT
    ASD_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/ASD_GeneLofZ.csv")
    SCZ_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/SCZ_GeneLofZ.csv")
    JobArrays = np.arange(1000)
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(process_Similarity_ASD_SCZ, [(idx, CT_Z2_MAT_HC, ASD_GeneLofZ, SCZ_GeneLofZ) for idx in JobArrays])
    pool.close()
    pool.join()
    ALL_RES = []
    for List in results:
        ALL_RES.append(List)
    ALL_RES = np.array(ALL_RES)
    with open("ASD_SCZ_BiasCorrRandomGeneRemove.npy", 'wb') as f:
        np.save(f, ALL_RES)

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=str, help='Mode of program 1: gene weight generateion; 2. bias calculaton')
    parser.add_argument('-o', '--outdir', type=str, help='Output directory')
    parser.add_argument('-w', '--WeightDF', type=str, help='Weight DF for control geneset')
    parser.add_argument('-p', '--GeneProb', default=None, help='GeneProb Filname or None if dont use')
    parser.add_argument('--n_sims', type=int, default=10000, help='Number of ctrl simulations')

    parser.add_argument('--GW_Dir', type=str, help="dirctory of ctrl gene weights")
    parser.add_argument('--SpecMat', type=str, help="Filename of bias matrix")
    parser.add_argument('--n_processes', type=int, default=20, help="Filename of bias matrix")
    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    Similarity_ASD_SCZ()
    return

if __name__ == '__main__':
    main()
