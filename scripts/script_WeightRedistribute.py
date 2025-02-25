#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_WeightRedistribute.py
# Generate weights according to the distribution of mutations
# Includde De novo mutations / Case control mutations  
# ========================================================================================================

import argparse
import sys
sys.path.insert(1, '/home/jw3514/Work/CellType_Psy/src')
from CellType_PSY import *
import multiprocessing
from multiprocessing import Pool


HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

def normalize(probs):
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]
def Reduce(List, Counts):
    res = []
    Counts = list(Counts)
    for L in List:
        res.append(Counts.count(L))
    return res

##########################
# ASD ReWeights
##########################
def Table2GW_ASD(Mutable):
    gene2MutN = {}
    for g, row in Mutable.iterrows():
        gene2MutN[g] = row["NLGD"] * 0.554 + row["NDmis"] * 0.333
    return gene2MutN

# strategy 1: Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
# strategy 2: Randomly redistribute the mutations, coording to background mutation rate
# strategy 3: Randomly redistribute the mutations, and assign new weights to random genes 
# strategy 4: Randomly redistribute the mutations, and assign new weights to expression matched genes
def SingleSet_ReWeights_ASD(idx, DZ_Muts, type=1, BGMR=None, BiasMat=None, ExpMatchDir=None, outDir=None):
    np.random.seed(idx)
    N_LGD = DZ_Muts[DZ_Muts["GeneEff"]!="missense"].shape[0]
    N_Dmis = DZ_Muts[DZ_Muts["GeneEff"]=="missense"].shape[0]
    
    Genes = DZ_Muts["Entrez"].unique()
    if type == 1: # Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
        Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 2: # Randomly redistribute the mutations, coording to background mutation rate
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 3: # Randomly redistribute the mutations, and assign new weights to random genes 
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        tmp_Mutable.index = np.random.choice(BiasMat.index.values, size=tmp_Mutable.shape[0], replace=False)
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 4: # Randomly redistribute the mutations, and assign new weights to expression matched genes
        Permed_LGD = np.random.choice(Genes, size=N_LGD, replace=True)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, replace=True)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for g, row in tmp_Mutable.iterrows():
            match_genes = loadgenelist(ExpMatchDir+"/{}.csv".format(g), toint=True)
            match_gene = np.random.choice(match_genes, 1)
            tmp_Mutable.loc[g, "MatchG"] = match_gene
        tmp_Mutable = tmp_Mutable.set_index("MatchG")
        perm_GW = Table2GW_ASD(tmp_Mutable)

    Dict2Fil(perm_GW, "{}/{}.perm.GW.txt".format(outDir, idx))
    Perm_Z2_Bias = AvgCTZ_Weighted(BiasMat, perm_GW, Method = 1)
    Perm_Z2_Bias = AnnotateCTDat(Perm_Z2_Bias, Anno)
    Perm_Z2_Bias.to_csv("{}/{}.perm.Z2.csv.gz".format(outDir, idx))
    return perm_GW


def SingleSet_ReWeights_ASD_ALL(idx, ASD_GeneDF, type=1, BGMR=None, BiasMat=None, ExpMatchDir=None, outDir=None):
    np.random.seed(idx)
    N_LGD = ASD_GeneDF["AutismMerged_LoF"].sum().astype('int')
    N_Dmis = ASD_GeneDF["AutismMerged_Dmis_REVEL0.5"].sum().astype('int')
    Genes = ASD_GeneDF["EntrezID"].unique()

    if type == 1: # Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
        Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 2: # Randomly redistribute the mutations, coording to background mutation rate
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 3: # Randomly redistribute the mutations, and assign new weights to random genes 
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        tmp_Mutable.index = np.random.choice(BiasMat.index.values, size=tmp_Mutable.shape[0], replace=False)
        perm_GW = Table2GW_ASD(tmp_Mutable)
    if type == 4: # Randomly redistribute the mutations, and assign new weights to expression matched genes
        Permed_LGD = np.random.choice(Genes, size=N_LGD, replace=True)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, replace=True)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for g, row in tmp_Mutable.iterrows():
            match_genes = loadgenelist(ExpMatchDir+"/{}.csv".format(g), toint=True)
            match_gene = np.random.choice(match_genes, 1)
            tmp_Mutable.loc[g, "MatchG"] = match_gene
        tmp_Mutable = tmp_Mutable.set_index("MatchG")
        perm_GW = Table2GW_ASD(tmp_Mutable)

    Dict2Fil(perm_GW, "{}/{}.perm.GW.txt".format(outDir, idx))
    Perm_Z2_Bias = AvgCTZ_Weighted(BiasMat, perm_GW, Method = 1)
    Perm_Z2_Bias = AnnotateCTDat(Perm_Z2_Bias, Anno)
    Perm_Z2_Bias.to_csv("{}/{}.perm.Z2.csv.gz".format(outDir, idx))
    return perm_GW

##########################
# SCZ ReWeights
##########################
def Table2GW_SCZ(Mutable):
    gene2MutN = {}
    for g, row in Mutable.iterrows():
        gene2MutN[g] = row["nLGD"] * 0.26 + row["nMis3"] * 0.25 + row["nMis2"] * 0.06
    return gene2MutN

def ModifyMutCount(CaseCount, ContCount, dnvCount, CaseN=24248, ContN=97322):
    if dnvCount!=dnvCount:
        dnvCount = 0
    return max(CaseCount - ContCount / ContN * CaseN + dnvCount, 0)

# Stratgy 1: Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
# Fixed genes (all genes are still SCZ genes)
def SingleSet_ReWeight_SCZ(idx, SCZ_MutDF, type=1, BiasMat=None, ExpMatchDir=None, outDir=None):
    #Mut_types = ["PTV", "mis3", "mis2"]
    np.random.seed(idx)
    TotalCase = 24248
    TotalCtrl = 97322
    N_PTV = SCZ_MutDF["Case PTV"].sum().astype('int')
    N_mis3 = SCZ_MutDF["Case mis3"].sum().astype('int')
    N_mis2 = SCZ_MutDF["Case mis2"].sum().astype('int')
    
    Genes = SCZ_MutDF.index.unique()
    if type == 1: # Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
        Permed_PTV = np.random.choice(Genes, size=N_PTV, replace=True)
        Permed_PTV_Count = Reduce(Genes, Permed_PTV)
        Permed_Mis3 = np.random.choice(Genes, size=N_mis3, replace=True)
        Permed_Mis3_Count = Reduce(Genes, Permed_Mis3)
        Permed_Mis2 = np.random.choice(Genes, size=N_mis2, replace=True)
        Permed_Mis2_Count = Reduce(Genes, Permed_Mis2)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "Case PTV":Permed_PTV_Count,
                                        "Case mis3":Permed_Mis3_Count, "Case mis2":Permed_Mis2_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for i, row in tmp_Mutable.iterrows():
            tmp_Mutable.loc[i, "Ctrl PTV"] = SCZ_MutDF.loc[i, "Ctrl PTV"]
            tmp_Mutable.loc[i, "Ctrl mis3"] = SCZ_MutDF.loc[i, "Ctrl mis3"]
            tmp_Mutable.loc[i, "Ctrl mis2"] = SCZ_MutDF.loc[i, "Ctrl mis2"]
            #tmp_Mutable.loc[i, "nLGD"] = ModifyMutCount(tmp_Mutable.loc[i, "Case PTV"], tmp_Mutable.loc[i, "Ctrl PTV"], 0)
            #tmp_Mutable.loc[i, "nMis3"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis3"], tmp_Mutable.loc[i, "Ctrl mis3"], 0)
            #tmp_Mutable.loc[i, "nMis2"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis2"], tmp_Mutable.loc[i, "Ctrl mis2"], 0)
            tmp_Mutable.loc[i, "nLGD"] = tmp_Mutable.loc[i, "Case PTV"]
            tmp_Mutable.loc[i, "nMis3"] = tmp_Mutable.loc[i, "Case mis3"]
            tmp_Mutable.loc[i, "nMis2"] = tmp_Mutable.loc[i, "Case mis2"]
        perm_GW = Table2GW_SCZ(tmp_Mutable)
    if type == 2: # Randomly redistribute the mutations, coording to case syn counts 
        PTV_P, mis3_P, mis2_P = [],[],[],
        for g in Genes:
            PTV_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
            mis3_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
            mis2_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
        Permed_PTV = np.random.choice(Genes, size=N_PTV, p=normalize(PTV_P), replace=True)
        Permed_PTV_Count = Reduce(Genes, Permed_PTV)
        Permed_mis3 = np.random.choice(Genes, size=N_mis3, p=normalize(mis3_P), replace=True)
        Permed_mis3_Count = Reduce(Genes, Permed_mis3)
        Permed_mis2 = np.random.choice(Genes, size=N_mis2, p=normalize(mis2_P), replace=True)
        Permed_mis2_Count = Reduce(Genes, Permed_mis2)
        
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "Case PTV":Permed_PTV_Count,
                                        "Case mis3":Permed_mis3_Count, "Case mis2":Permed_mis2_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for i, row in tmp_Mutable.iterrows():
            tmp_Mutable.loc[i, "Ctrl PTV"] = SCZ_MutDF.loc[i, "Ctrl PTV"]
            tmp_Mutable.loc[i, "Ctrl mis3"] = SCZ_MutDF.loc[i, "Ctrl mis3"]
            tmp_Mutable.loc[i, "Ctrl mis2"] = SCZ_MutDF.loc[i, "Ctrl mis2"]
            tmp_Mutable.loc[i, "nLGD"] = ModifyMutCount(tmp_Mutable.loc[i, "Case PTV"], tmp_Mutable.loc[i, "Ctrl PTV"], 0)
            tmp_Mutable.loc[i, "nMis3"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis3"], tmp_Mutable.loc[i, "Ctrl mis3"], 0)
            tmp_Mutable.loc[i, "nMis2"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis2"], tmp_Mutable.loc[i, "Ctrl mis2"], 0)
        perm_GW = Table2GW_SCZ(tmp_Mutable)

    if type == 3: # Randomly redistribute the mutations, coording to control counts of each mutation type
        PTV_P, mis3_P, mis2_P = [],[],[],
        for g in Genes:
            PTV_P.append(SCZ_MutDF.loc[g, "Ctrl PTV"]+1)
            mis3_P.append(SCZ_MutDF.loc[g, "Ctrl mis3"]+1)
            mis2_P.append(SCZ_MutDF.loc[g, "Ctrl mis2"]+1)
        Permed_PTV = np.random.choice(Genes, size=N_PTV, p=normalize(PTV_P), replace=True)  
        Permed_PTV_Count = Reduce(Genes, Permed_PTV)
        Permed_mis3 = np.random.choice(Genes, size=N_mis3, p=normalize(mis3_P), replace=True)
        Permed_mis3_Count = Reduce(Genes, Permed_mis3)
        Permed_mis2 = np.random.choice(Genes, size=N_mis2, p=normalize(mis2_P), replace=True)
        Permed_mis2_Count = Reduce(Genes, Permed_mis2)
        
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "Case PTV":Permed_PTV_Count,
                                        "Case mis3":Permed_mis3_Count, "Case mis2":Permed_mis2_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for i, row in tmp_Mutable.iterrows():
            tmp_Mutable.loc[i, "Ctrl PTV"] = SCZ_MutDF.loc[i, "Ctrl PTV"]
            tmp_Mutable.loc[i, "Ctrl mis3"] = SCZ_MutDF.loc[i, "Ctrl mis3"]
            tmp_Mutable.loc[i, "Ctrl mis2"] = SCZ_MutDF.loc[i, "Ctrl mis2"]
            tmp_Mutable.loc[i, "nLGD"] = tmp_Mutable.loc[i, "Case PTV"]
            tmp_Mutable.loc[i, "nMis3"] = tmp_Mutable.loc[i, "Case mis3"]
            tmp_Mutable.loc[i, "nMis2"] = tmp_Mutable.loc[i, "Case mis2"]
            tmp_Mutable.loc[i, "nLGD"] = ModifyMutCount(tmp_Mutable.loc[i, "Case PTV"], tmp_Mutable.loc[i, "Ctrl PTV"], 0)
            tmp_Mutable.loc[i, "nMis3"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis3"], tmp_Mutable.loc[i, "Ctrl mis3"], 0)
            tmp_Mutable.loc[i, "nMis2"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis2"], tmp_Mutable.loc[i, "Ctrl mis2"], 0)
        perm_GW = Table2GW_SCZ(tmp_Mutable)
    if type == 4: # Randomly redistribute the mutations, and assign new weights to random genes
        PTV_P, mis3_P, mis2_P = [],[],[],
        for g in Genes:
            PTV_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
            mis3_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
            mis2_P.append(SCZ_MutDF.loc[g, "Ctrl syn"]+1)
        Permed_PTV = np.random.choice(Genes, size=N_PTV, p=normalize(PTV_P), replace=True)  
        Permed_PTV_Count = Reduce(Genes, Permed_PTV)
        Permed_mis3 = np.random.choice(Genes, size=N_mis3, p=normalize(mis3_P), replace=True)
        Permed_mis3_Count = Reduce(Genes, Permed_mis3)
        Permed_mis2 = np.random.choice(Genes, size=N_mis2, p=normalize(mis2_P), replace=True)
        Permed_mis2_Count = Reduce(Genes, Permed_mis2)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "Case PTV":Permed_PTV_Count,
                                        "Case mis3":Permed_mis3_Count, "Case mis2":Permed_mis2_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")   
        for i, row in tmp_Mutable.iterrows():
            tmp_Mutable.loc[i, "Ctrl PTV"] = SCZ_MutDF.loc[i, "Ctrl PTV"]
            tmp_Mutable.loc[i, "Ctrl mis3"] = SCZ_MutDF.loc[i, "Ctrl mis3"]
            tmp_Mutable.loc[i, "Ctrl mis2"] = SCZ_MutDF.loc[i, "Ctrl mis2"]
            tmp_Mutable.loc[i, "nLGD"] = ModifyMutCount(tmp_Mutable.loc[i, "Case PTV"], tmp_Mutable.loc[i, "Ctrl PTV"], 0)
            tmp_Mutable.loc[i, "nMis3"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis3"], tmp_Mutable.loc[i, "Ctrl mis3"], 0)
            tmp_Mutable.loc[i, "nMis2"] = ModifyMutCount(tmp_Mutable.loc[i, "Case mis2"], tmp_Mutable.loc[i, "Ctrl mis2"], 0)

        tmp_Mutable.index = np.random.choice(BiasMat.index.values, size=tmp_Mutable.shape[0], replace=False)
        perm_GW = Table2GW_SCZ(tmp_Mutable)
    Dict2Fil(perm_GW, "{}/{}.perm.GW.txt".format(outDir, idx))
    Perm_Z2_Bias = AvgCTZ_Weighted(BiasMat, perm_GW, Method = 1)
    Perm_Z2_Bias = AnnotateCTDat(Perm_Z2_Bias, Anno)
    Perm_Z2_Bias.to_csv("{}/{}.perm.Z2.csv.gz".format(outDir, idx))
    return perm_GW


def Table2GW_NDD(Mutable):
    gene2MutN = {}
    for g, row in Mutable.iterrows():
        gene2MutN[g] = row["NLGD"] *  0.347 + row["NDmis"] * 0.194
    return gene2MutN
def SingleSet_ReWeights_NDD(idx, NDD_GeneDF, type=1, BGMR=None, BiasMat=None, ExpMatchDir=None, outDir=None):
    np.random.seed(idx)
    N_LGD = NDD_GeneDF["frameshift_variant"].sum().astype('int') + NDD_GeneDF["splice_acceptor_variant"].sum().astype('int') + NDD_GeneDF["splice_donor_variant"].sum().astype('int') + NDD_GeneDF["stop_gained"].sum().astype('int') + NDD_GeneDF["stop_lost"].sum().astype('int')
    N_Dmis = NDD_GeneDF["missense_variant"].sum().astype('int') 

    Genes = NDD_GeneDF["EntrezID"].unique()

    if type == 1: # Randomly redistribute the mutations, coording to nothing (each mutation is independent, property of genes doesnt play any role)
        Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_NDD(tmp_Mutable)
    if type == 2: # Randomly redistribute the mutations, coording to background mutation rate
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        perm_GW = Table2GW_NDD(tmp_Mutable)
    if type == 3: # Randomly redistribute the mutations, and assign new weights to random genes 
        LGD_P = []
        Dmis_P = []
        for g in Genes:
            try:
                LGD_P.append(BGMR.loc[g, "p_LGD"])
                Dmis_P.append(BGMR.loc[g, "prevel_0.5"])
            except:
                LGD_P.append(BGMR["p_LGD"].mean())
                Dmis_P.append(BGMR["prevel_0.5"].mean())
        Permed_LGD = np.random.choice(Genes, size=N_LGD, p=normalize(LGD_P), replace=True)
        #Permed_LGD = np.random.choice(Genes, size=N_LGD)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, p=normalize(Dmis_P), replace=True)
        #Permed_Dmis = np.random.choice(Genes, size=N_Dmis)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        tmp_Mutable.index = np.random.choice(BiasMat.index.values, size=tmp_Mutable.shape[0], replace=False)
        perm_GW = Table2GW_NDD(tmp_Mutable)
    if type == 4: # Randomly redistribute the mutations, and assign new weights to expression matched genes
        Permed_LGD = np.random.choice(Genes, size=N_LGD, replace=True)
        Permed_LGD_Count = Reduce(Genes, Permed_LGD)
        Permed_Dmis = np.random.choice(Genes, size=N_Dmis, replace=True)
        Permed_Dmis_Count = Reduce(Genes, Permed_Dmis)
        tmp_Mutable = pd.DataFrame(data={"Entrez":Genes, "NLGD":Permed_LGD_Count,
                                        "NDmis":Permed_Dmis_Count})
        tmp_Mutable = tmp_Mutable.set_index("Entrez")
        for g, row in tmp_Mutable.iterrows():
            match_genes = loadgenelist(ExpMatchDir+"/{}.csv".format(g), toint=True)
            match_gene = np.random.choice(match_genes, 1)
            tmp_Mutable.loc[g, "MatchG"] = match_gene
        tmp_Mutable = tmp_Mutable.set_index("MatchG")
        perm_GW = Table2GW_NDD(tmp_Mutable)

    Dict2Fil(perm_GW, "{}/{}.perm.GW.txt".format(outDir, idx))
    Perm_Z2_Bias = AvgCTZ_Weighted(BiasMat, perm_GW, Method = 1)
    Perm_Z2_Bias = AnnotateCTDat(Perm_Z2_Bias, Anno)
    Perm_Z2_Bias.to_csv("{}/{}.perm.Z2.csv.gz".format(outDir, idx))
    return perm_GW

def SinglePerm_TwoSetReWeights(idx, ASD_MutationDat, SCZ_MutationDat, HCT_Z2_MAT_HCT, BGMR, outDir):
    ASD_RW_GW = SingleSet_ReWeights_ASD(ASD_MutationDat, BGMR)
    RW_ASD_Bias = AvgCTZ_Weighted(HCT_Z2_MAT_HCT, ASD_RW_GW)
    RW_ASD_Bias = AnnotateCTDat(RW_ASD_Bias, Anno)

    SCZ_RW_GW = SingleSet_ReWeight_SCZ(SCZ_MutationDat)
    #SCZ_RW_GW = SingleSet_ReWeight_SCZ_v2(SCZ_MutationDat, BGMR)
    RW_SCZ_Bias = AvgCTZ_Weighted(HCT_Z2_MAT_HCT, SCZ_RW_GW)
    RW_SCZ_Bias = AnnotateCTDat(RW_SCZ_Bias, Anno)

    RW_ASD_Bias.to_csv("{}/{}.ASD.perm.bias.csv".format(outDir, idx))
    RW_SCZ_Bias.to_csv("{}/{}.SCZ.perm.bias.csv".format(outDir, idx))
    return

def TwoSetReWeights(outDir, N_perms=10000, n_processes=20): 
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    HCT_Z2_MAT_HCT = pd.read_csv("../dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.log2.Z2.HCT.z1clip3.csv", index_col=0)
    max_Z, min_Z = 3, -3
    HCT_Z2_MAT_HCT = HCT_Z2_MAT_HCT.clip(upper=max_Z, lower=min_Z)

    BGMR = pd.read_csv("~/Work/Resources/MutationRate_20170710_rate.txt", delimiter="\t")
    BGMR["Entrez"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in BGMR["GeneName"].values]
    BGMR = BGMR[BGMR["Entrez"].isin(HCT_Z2_MAT_HCT.index.values)]
    BGMR.index = BGMR["Entrez"].values

    
    ASD_MutationDat = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/GeneWeights/ASD.HighIQ.Mutable.csv", index_col=0)
    SCZ_MutationDat = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/GeneWeights/SCZ.top61.Mutable.csv", index_col=0)

    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(SinglePerm_TwoSetReWeights, 
        [(idx, ASD_MutationDat, SCZ_MutationDat, HCT_Z2_MAT_HCT, BGMR, outDir) for idx in np.arange(N_perms)])
    pool.close()
    pool.join()

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('--BiasMat', type=str,
                        help='Cell Type Bias matrix')
    parser.add_argument('--dataset', type=str, choices=['ASD', 'SCZ', 'NDD', 'ASD_ALL'],
                        help='dataset aim for permutation')
    parser.add_argument('--mutable', type=str,
                        help='mutation table used for permutation')
    parser.add_argument('--sim_type', type=int,
                        help='type of simulation strategy')
    parser.add_argument('-n', '--N_perms', type=int, default = 1000, help='Number of permutations')
    parser.add_argument('--n_processes', type=int, default = 10, help='number of processes')
    parser.add_argument('-d', '--dir', type=str, required=True, help='Save folder')
    args = parser.parse_args()

    return args

def LoadBGMR(BiasMat):
    BGMR = pd.read_csv("~/Work/Resources/MutationRate_20170710_rate.txt", delimiter="\t")
    BGMR["Entrez"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in BGMR["GeneName"].values]
    BGMR = BGMR[BGMR["Entrez"].isin(BiasMat.index.values)]
    BGMR.index = BGMR["Entrez"].values
    return BGMR

def main():
    args = GetOptions()
    #n_processes = args.n_processes
    BiasMat = pd.read_csv(args.BiasMat, index_col=0)
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    
    pool = multiprocessing.Pool(processes=args.n_processes)
    

    if args.dataset == 'ASD':
        DZ_Muts = pd.read_csv(args.mutable, index_col=0)
        if args.sim_type == 1: # random mutation
            results = pool.starmap(SingleSet_ReWeights_ASD, 
                                    [(idx, DZ_Muts, args.sim_type, None, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 2: # background mutation rate
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap(SingleSet_ReWeights_ASD, 
                                    [(idx, DZ_Muts, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 3: # random weight
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap(SingleSet_ReWeights_ASD, 
                                    [(idx, DZ_Muts, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 4: # weight matched
            BGMR = LoadBGMR(BiasMat)
            SingleSet_ReWeights_ASD(args.mutable, type=args.sim_type, BGMR=BGMR, HCT_Z2_MAT_HCT=None, ExpMatchDir=None, outDir=args.dir)
    elif args.dataset == 'ASD_ALL':
        ASD_Gene_FILE = pd.read_csv(args.mutable, index_col=0)
        if args.sim_type == 1: # random mutation
            results = pool.starmap( SingleSet_ReWeights_ASD_ALL, 
                                    [(idx, ASD_Gene_FILE, args.sim_type, None, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 2: # background mutation rate
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap( SingleSet_ReWeights_ASD_ALL, 
                                    [(idx, ASD_Gene_FILE, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 3: # random weight
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap( SingleSet_ReWeights_ASD_ALL, 
                                    [(idx, ASD_Gene_FILE, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 4: # weight matched
            BGMR = LoadBGMR(BiasMat)
            SingleSet_ReWeights_ASD_ALL(args.mutable, type=args.sim_type, BGMR=BGMR, HCT_Z2_MAT_HCT=None, ExpMatchDir=None, outDir=args.dir)
    elif args.dataset == 'SCZ':
        DZ_Muts = pd.read_csv(args.mutable, index_col=0)
        if args.sim_type == 1: # random mutation
            results = pool.starmap(SingleSet_ReWeight_SCZ, 
                                    [(idx, DZ_Muts, args.sim_type, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 2: # sim with syn rate
            results = pool.starmap(SingleSet_ReWeight_SCZ, 
                                    [(idx, DZ_Muts, args.sim_type, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 3: # sim with control mutation counts
            results = pool.starmap(SingleSet_ReWeight_SCZ, 
                                    [(idx, DZ_Muts, args.sim_type, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 4: # weight matched
            results = pool.starmap(SingleSet_ReWeight_SCZ, 
                                    [(idx, DZ_Muts, args.sim_type, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])

    elif args.dataset == 'NDD':
        NDD_Gene_FILE = pd.read_csv(args.mutable, index_col=0)
        if args.sim_type == 1: # random mutation
            results = pool.starmap( SingleSet_ReWeights_NDD, 
                                    [(idx, NDD_Gene_FILE, args.sim_type, None, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 2: # background mutation rate
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap( SingleSet_ReWeights_NDD, 
                                    [(idx, NDD_Gene_FILE, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 3: # random weight
            BGMR = LoadBGMR(BiasMat)
            results = pool.starmap( SingleSet_ReWeights_NDD, 
                                    [(idx, NDD_Gene_FILE, args.sim_type, BGMR, BiasMat, None, args.dir) for idx in np.arange(args.N_perms)])
        elif args.sim_type == 4: # weight matched
            BGMR = LoadBGMR(BiasMat)
            SingleSet_ReWeights_ASD_ALL(args.mutable, type=args.sim_type, BGMR=BGMR, HCT_Z2_MAT_HCT=None, ExpMatchDir=None, outDir=args.dir)

    pool.close()
    pool.join()
    return


if __name__ == '__main__':
    main()
