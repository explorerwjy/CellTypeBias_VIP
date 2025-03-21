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

def process_Similarity_ASD_SCZ_HumanCT(idx, CT_Z2_MAT_HC, ASD_GeneLofZ, SCZ_GeneLofZ):
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
    return yy 

def Similarity_ASD_SCZ_HumanCT(n_processes=20):
    HCT_Z2_MAT_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.Cluster.Log2Mean.Z1clip5.Z2.clip3.Dec30.csv", index_col=0)
    HCT_Z2_MAT_HCT.columns = HCT_Z2_MAT_HCT.columns.astype(int)
    CT_Z2_MAT_HC = HCT_Z2_MAT_HCT
    #ASD_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/ASD_GeneLofZ.csv")
    #SCZ_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/SCZ_GeneLofZ.csv")
    ASD_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Other/ASD_GeneLofZ.LGD_Dmis_SameWeight.csv")
    SCZ_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Other/SCZ_GeneLofZ.LGD_Dmis_SameWeight.csv")
    JobArrays = np.arange(1000)
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(process_Similarity_ASD_SCZ_HumanCT, [(idx, CT_Z2_MAT_HC, ASD_GeneLofZ, SCZ_GeneLofZ) for idx in JobArrays])
    pool.close()
    pool.join()
    ALL_RES = []
    for List in results:
        ALL_RES.append(List)
    ALL_RES = np.array(ALL_RES)
    #with open("../dat/Other/ASD_SCZ_HumanCT_BiasCorrRandomGeneRemove.npy", 'wb') as f:
    with open("../dat/Other/ASD_SCZ_MouseSTR_BiasCorrRandomGeneRemove.LGD_Dmis_SameWeight.npy", 'wb') as f:
        np.save(f, ALL_RES)

def process_Similarity_ASD_SCZ_MouseSTR(idx, CT_Z2_MAT_HC, ASD_GeneLofZ, SCZ_GeneLofZ):
    DF1 = ASD_GeneLofZ.sample(frac=1, random_state=idx)
    DF2 = SCZ_GeneLofZ.sample(frac=1, random_state=idx)
    yy = []
    for i in range(0, 31, 1):
        tmp_ASD_GW = dict(zip(DF1["Entrez"].values[i:], DF1["GW"].values[i:]))
        tmp_SCZ_GW = dict(zip(DF2["Entrez"].values[i:], DF2["GW"].values[i:]))
        tmp_ASD_Bias = AvgSTRZ_Weighted(CT_Z2_MAT_HC, tmp_ASD_GW, Method = 1)
        tmp_SCZ_Bias = AvgSTRZ_Weighted(CT_Z2_MAT_HC, tmp_SCZ_GW, Method = 1)
        r,p = GetSingeCellBiasCorr(tmp_ASD_Bias, tmp_SCZ_Bias)
        yy.append(r)
    return yy 

def Similarity_ASD_SCZ_MouseSTR(n_processes=20):
    Mouse_STR_Z2_Mat = pd.read_csv("../../ASD_Circuits/dat/allen-mouse-exp/AllenMouseBrain_Z2bias.csv", index_col=0)
    #ASD_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/ASD_GeneLofZ.csv")
    #SCZ_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/notebooks3/SCZ_GeneLofZ.csv")
    ASD_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Other/ASD_GeneLofZ.LGD_Dmis_SameWeight.csv")
    SCZ_GeneLofZ = pd.read_csv("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/Other/SCZ_GeneLofZ.LGD_Dmis_SameWeight.csv")
    JobArrays = np.arange(1000)
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(process_Similarity_ASD_SCZ_MouseSTR, [(idx, Mouse_STR_Z2_Mat, ASD_GeneLofZ, SCZ_GeneLofZ) for idx in JobArrays])
    pool.close()
    pool.join()
    ALL_RES = []
    for List in results:
        ALL_RES.append(List)
    ALL_RES = np.array(ALL_RES)
    #with open("ASD_SCZ_MouseSTR_BiasCorrRandomGeneRemove.npy", 'wb') as f:
    with open("../dat/Other/ASD_SCZ_MouseSTR_BiasCorrRandomGeneRemove.LGD_Dmis_SameWeight.npy", 'wb') as f:
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
    Similarity_ASD_SCZ_HumanCT()
    #Similarity_ASD_SCZ_MouseSTR()
    return

if __name__ == '__main__':
    main()
