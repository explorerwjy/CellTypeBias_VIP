# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_run_ctrl_sim.v3.py
# Run control simutations for expression biases (Cell Type; Structures)
# ========================================================================================================

import argparse
import sys
import os
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *
import time
import json
import pickle
import multiprocessing
from multiprocessing import Pool

###########################################################################
## Human Cell Type CTRL Generation
###########################################################################

def BiasCal_SingleCtrl_HumanCT(Random_GW_DF, SpecMat, idx):
    GW = dict(zip(Random_GW_DF[idx].values, Random_GW_DF["GeneWeight"].values))
    Ctrl_BiasDF = HumanCT_AvgZ_Weighted(SpecMat, GW)
    return Ctrl_BiasDF

def CtrlBiasCal_HumanCT(Random_GW_Fil, SpecMat, outfile, n_processes=20): 
    print("HumanCT Null Bias Simulation -", Random_GW_Fil)

    Random_GW_DF = pd.read_csv(Random_GW_Fil, index_col=0) # index is real gene,first column is gene weights, rest are random gene entrez
    Task_idx = Random_GW_DF.columns.values[1:]
    #print(Task_idx)

    SpecMat = pd.read_csv(SpecMat, index_col=0)
    SpecMat.columns = SpecMat.columns.astype(int)

    pool = multiprocessing.Pool(processes=n_processes)
    result_dfs = pool.starmap(BiasCal_SingleCtrl_HumanCT, [(Random_GW_DF, SpecMat, idx) for idx in Task_idx])

    pool.close()
    pool.join()
    print("end multiprocessing, merging results")
    str_ids = sorted(result_dfs[0].index.values)
    effects_matrix = np.stack([df.loc[str_ids, "EFFECT"].values for df in result_dfs], axis=0)
    effects_df = pd.DataFrame(effects_matrix.T, index=str_ids, columns=[str(i) for i in range(effects_matrix.shape[0])])
    effects_df.to_csv(outfile)
    return
    

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
    parser.add_argument('-o', '--outfile', type=str, help='Output file')
    parser.add_argument('--Ctrl_Genes_Fil', type=str, required=True, help="Filename of ctrl genes")
    parser.add_argument('--SpecMat', type=str, help="Filename of bias matrix")
    parser.add_argument('--n_processes', type=int, default=20, help="Filename of bias matrix")
    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    mode = args.mode
    if mode == 'bias':
        SpecMat = args.SpecMat
        Ctrl_Genes_Fil = args.Ctrl_Genes_Fil
        outfile = args.outfile
        n_processes = args.n_processes
        CtrlBiasCal(Ctrl_Genes_Fil, SpecMat, outfile, n_processes) # 
    if mode == 'human_ct_bias':
        SpecMat = args.SpecMat
        Ctrl_Genes_Fil = args.Ctrl_Genes_Fil
        outfile = args.outfile
        n_processes = args.n_processes
        CtrlBiasCal_HumanCT(Ctrl_Genes_Fil, SpecMat, outfile, n_processes)
    if mode == "mouse_ct_bias":
        DN=args.DN
        print(DN)
        #DN = 0
        GW_Dir = args.GW_Dir
        SpecMat = args.SpecMat
        outDir = args.outdir
        n_processes = args.n_processes
        CtrlBiasCal_MouseCT(Ctrl_Genes_Fil, SpecMat, outfile, n_processes, DN)

    return

if __name__ == '__main__':
    main()
