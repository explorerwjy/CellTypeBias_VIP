# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (gencic)
#     language: python
#     name: gencic
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import sys
ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *
import os 
#import scanpy as sc
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

try:
    os.chdir(f"{ProjDIR}/notebooks/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:
    print(f"Unexpected error: {e}")

# %%
HumanCT_Z2_HCT = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv", index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
HumanCT_Z2_HCT.shape

# %% [markdown] heading_collapsed=true
# # SCZ

# %% hidden=true
GeneDF = pd.read_excel("../dat/41586_2022_4556_MOESM3_ESM.xlsx",
                    sheet_name="Table S5 - Gene Results")
ExAC_pLI = pd.read_csv("/home/jw3514/Work/Resources/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t",
                      index_col="gene")

# %%
Ncase = 24248
Nctrl = 97322

# %%
GeneDF.head(2)


# %% [markdown]
# #### Modify Mutation Counts

# %% hidden=true
def oddsratio(NcaseMut, NctrlMut, dnvCount, Ncase = 24248, Nctrl=97322):
    if dnvCount!=dnvCount:
        dnvCount = 0
    NcaseMut += 1 # add pseudo count to prevent inf 
    NctrlMut += 1
    AD = (NcaseMut) * (Nctrl-NctrlMut) 
    BC = (NctrlMut) * (Ncase-NcaseMut)
    return AD/BC + dnvCount

def Penetrance(NcaseMut, NctrlMut, Ncase = 24248, Nctrl=97322, prevelence=0.45/100):
    NcaseMut += 1
    NctrlMut += 1
    Ntotal_ =  Ncase/prevelence
    Nctrl_ = Ntotal_ - Ncase
    NctrlMut_ = Nctrl_*(NctrlMut/Nctrl)
    #print(Nctrl_, NctrlMut_)
    p = NcaseMut/(NcaseMut+NctrlMut_)
    return p


# %% hidden=true
print(Penetrance(15, 3))
print(Penetrance(8, 1))
print(Penetrance(12, 1))
print(Penetrance(129, 466))

# %% hidden=true
print(oddsratio(15, 3, 0))
print(oddsratio(8, 1, 0))
print(oddsratio(12, 1, 0))


# %% hidden=true
# Effective mutation count in case  
def ModifyMutCount(CaseCount, ContCount, dnvCount, CaseN=24248, ContN=97322):
    if isinstance(dnvCount, float) and np.isnan(dnvCount):
        dnvCount = 0
    #return max(CaseCount - ContCount / ContN * CaseN + dnvCount, 0)
    return CaseCount - ContCount / ContN * CaseN + dnvCount
    
for i, row in GeneDF.iterrows():
    symbol = row["Gene Symbol"]
    #print(i, symbol)
    
    try:
        GeneDF.loc[i, "Entrez"] = int(GeneSymbol2Entrez[symbol])
    except:
        GeneDF.loc[i, "Entrez"] = None
        
    try:
        GeneDF.loc[i, "pLI"] = ExAC_pLI.loc[symbol, "pLI"]
    except:
        GeneDF.loc[i, "pLI"] = 0
        
    case_ptv = float(row["Case PTV"]) if not pd.isna(row["Case PTV"]) else 0
    ctrl_ptv = float(row["Ctrl PTV"]) if not pd.isna(row["Ctrl PTV"]) else 0
    case_mis3 = float(row["Case mis3"]) if not pd.isna(row["Case mis3"]) else 0 
    ctrl_mis3 = float(row["Ctrl mis3"]) if not pd.isna(row["Ctrl mis3"]) else 0
    case_mis2 = float(row["Case mis2"]) if not pd.isna(row["Case mis2"]) else 0
    ctrl_mis2 = float(row["Ctrl mis2"]) if not pd.isna(row["Ctrl mis2"]) else 0
    
    dnv_ptv = float(row["De novo PTV"]) if not pd.isna(row["De novo PTV"]) else 0
    dnv_mis3 = float(row["De novo mis3"]) if not pd.isna(row["De novo mis3"]) else 0
    dnv_mis2 = float(row["De novo mis2"]) if not pd.isna(row["De novo mis2"]) else 0

    GeneDF.loc[i, "nLGD"] = ModifyMutCount(case_ptv, ctrl_ptv, 0)
    GeneDF.loc[i, "nMis3"] = ModifyMutCount(case_mis3, ctrl_mis3, 0)
    GeneDF.loc[i, "nMis2"] = ModifyMutCount(case_mis2, ctrl_mis2, 0)

    GeneDF.loc[i, "LGD_OR"] = oddsratio(case_ptv, ctrl_ptv, dnv_ptv)
    GeneDF.loc[i, "Mis3_OR"] = oddsratio(case_mis3, ctrl_mis3, dnv_mis3)
    GeneDF.loc[i, "Mis2_OR"] = oddsratio(case_mis2, ctrl_mis2, dnv_mis2)
    
    GeneDF.loc[i, "LGD_pen"] = Penetrance(case_ptv, ctrl_ptv)
    GeneDF.loc[i, "Mis3_pen"] = Penetrance(case_mis3, ctrl_mis3)
    GeneDF.loc[i, "Mis2_pen"] = Penetrance(case_mis2, ctrl_mis2)
    
GeneDF = GeneDF.dropna(subset=["Entrez"])
GeneDF = GeneDF.set_index("Entrez")
#GeneDF = GeneDF[GeneDF.index.isin(HumanCT_Z2_HCT.index.values)]
GeneDF.to_csv("../dat/SCZ.ALLGENE.MutCountModified.csv")
GeneDF.shape

# %%
GeneDF = pd.read_csv("../dat/SCZ.ALLGENE.MutCountModified.csv", index_col=0)

# %% hidden=true
GeneDF.head(5)

# %%
SCZ_Protext = GeneDF[GeneDF["nLGD"]<0].head(61)

# %%
SCZ_Protext

# %%
SCZ_61_protect_GW = Aggregate_Gene_Weights_SCZ_Daly(SCZ_Protext, HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
SCZ_61_protect_GW = { k: abs(v) for k,v in SCZ_61_protect_GW.items() }
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61_protect_GW, "{}/SCZ.top61.protect.gw".format(OutDIR))


# %%

# %%
def Aggregate_Gene_Weights_SCZ_Daly(MutFil, allen_mouse_genes, usepLI=False, Bmis=False, out=None, mode="MC", 
                                  lgd_weight=0.33, mis3_weight=0.27, mis2_weight=0.12):
    print("New")
    assert mode in ["OR", "MC", "ORMC"]
    print(mode)
    gene2MutN = {}
    for i, row in MutFil.iterrows():
        try:
            g = int(i)
            if g not in allen_mouse_genes:
                print(g, "not in Expression dataset")
                continue
        except:
            print(g, "Error converting Entrez ID")
        if usepLI:
            try:
                pLI = float(row["pLI"])
            except:
                print(g, "don't have pLI score on file, set to 0")
                pLI = 0.0
            if pLI >= 0.5:
                gene2MutN[g] = row["nLGD"] * 0.26 + row["nMis3"] * 0.25 + row["nMis2"] * 0.06  
            else:
                gene2MutN[g] = row["nLGD"] * 0.01 + row["nMis3"] * 0.01 + row["nMis2"] * 0 
        else:
            if mode == "OR":
                gene2MutN[g] = row["LGD_OR"] * lgd_weight + row["Mis3_OR"] * mis3_weight + row["Mis2_OR"] * mis2_weight
            elif mode == "MC":
                gene2MutN[g] = row["nLGD"] * lgd_weight + row["nMis3"] * mis3_weight + row["nMis2"] * mis2_weight
            ##elif mode == "ORMC":
            #    gene2MutN[g] = row["nLGD"] * lgd_weight + row["nMis3"] * mis3_weight + row["nMis2"] * mis2_weight
            #gene2MutN[g] = row["LGD_OR"] * 0.26 + row["Mis3_OR"] * 0.25 + row["nMis2_OR"] * 0.06
            #gene2MutN[g] = row["LGD_OR"] * 0.26 + row["Mis3_OR"] * 0.25 + row["nMis2_OR"] * 0.06
    if out != None:
        writer = csv.writer(open(out, 'wt'))
        for k,v in sorted(gene2MutN.items(), key=lambda x:x[1], reverse=True):
           writer.writerow([k,v]) 
    return gene2MutN



# %% hidden=true
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF.head(61), HumanCT_Z2_HCT.index.values, mode="MC")
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.nopLI.MC.gw".format(OutDIR))

SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF.head(61), HumanCT_Z2_HCT.index.values, mode="OR")
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.nopLI.OR.gw".format(OutDIR))

# %%
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF.head(61), HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.gw".format(OutDIR))

# %% hidden=true
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF.head(61), HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(OutDIR))

# %%
GeneDF_pos = GeneDF[GeneDF["nLGD"]>1]
SCZ_100GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF_pos.head(100), HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_100GW, "{}/SCZ.top100.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(OutDIR))

# %%
GeneDF_pos = GeneDF[GeneDF["nLGD"]>1]
SCZ_200GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF_pos.head(200), HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_200GW, "{}/SCZ.top200.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(OutDIR))

# %%
SCZ_500GW = Aggregate_Gene_Weights_SCZ_Daly(GeneDF_pos.head(500), HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_500GW, "{}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw".format(OutDIR))

# %% hidden=true
NLGDs = []
NDmis = []
Bias_CT = {}
for i, row in GeneDF.head(61).iterrows():
    try:
        gene = int(i)
        NLGDs.append(row["nLGD"])
        NDmis.append(row["nMis3"] + row["nMis2"])
        for CT in Neurons:
            CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
            if CT not in Bias_CT:
                Bias_CT[CT] = []
            Bias_CT[CT].append(CT_Z2_MAT_HC.loc[gene, CT_Idx].mean())
    except:
        print(i, "not in exp data")

# %% hidden=true
NLGDs = np.array(NLGDs)

# %% hidden=true
print(len(NLGDs))
print(len(Bias_CT["Amygdala excitatory"]))

# %% hidden=true
CT = "CGE interneuron"
CT = "MGE interneuron"
CT = "Thalamic excitatory"
_CT_Bias = np.array(Bias_CT[CT])
_CT_Bias > 0
NLGD_Pos = NLGDs[_CT_Bias > 0]
NLGD_Neg = NLGDs[_CT_Bias <= 0]
Total_G = len(NLGD_Pos) + len(NLGD_Neg)
print(len(NLGD_Pos), len(NLGD_Neg))
N_Mut_Pos = int(NLGD_Pos.sum())
Total_Mut = int(N_Mut_Pos + int(NLGD_Neg.sum()))
print(N_Mut_Pos, Total_Mut, N_Mut_Pos/Total_Mut / (len(NLGD_Pos)/Total_G))
print(scipy.stats.binom_test(N_Mut_Pos, Total_Mut, p=len(NLGD_Pos)/Total_G))

# %% [markdown] heading_collapsed=true
# # ASD

# %% [markdown]
# ## Complete ASD group (Not used in the paper)

# %% hidden=true
Spark_Meta_1stage = pd.read_excel("~/Work/SPARK2020/TabS_DenovoWEST_Stage1.xlsx",
                           skiprows=1, sheet_name="AllGenes")

# %% hidden=true
Spark_Meta_1stage.columns.values

# %% hidden=true
Spark_Meta_1stage.head(5)

# %% hidden=true
_, ASD_top200_GW = Aggregate_Gene_Weights2(Spark_Meta_1stage.head(200), ExpZ2Mat.index.values)
Dict2Fil(ASD_top200_GW, "../dat3/ASD.top200.gw.csv")

# %% hidden=true
len(ASD_top200_GW)

# %% hidden=true
Spark_Meta_2stage = pd.read_excel("../../ASD_Circuits//dat/genes/asd/TabS_DenovoWEST_Stage1+2.xlsx",
                           skiprows=2, sheet_name="TopDnEnrich")
Spark_Meta_2stage.shape

# %%
Spark_Meta_2stage.to_csv("../dat/Spark_TabS_DenovoWEST_Stage1+2.csv")

# %% hidden=true
Spark_Meta_2stage.columns.values

# %% code_folding=[] hidden=true
NLGDs = []
NDmis = []
Bias_CT = {}
for i, row in Spark_Meta_2stage.head(159).iterrows():
    gene = int(row["EntrezID"])
    NLGDs.append(row["AutismMerged_LoF"])
    NDmis.append(row["AutismMerged_Dmis_REVEL0.5"])
    for CT in Neurons:
        CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
        if CT not in Bias_CT:
            Bias_CT[CT] = []
        Bias_CT[CT].append(CT_Z2_MAT_HC.loc[gene, CT_Idx].mean())
NLGDs = np.array(NLGDs)
NDmis = np.array(NDmis)

# %% hidden=true
CT = "MGE interneuron"
CT = "Medium spiny neuron"
_CT_Bias = np.array(Bias_CT[CT])
_CT_Bias > 0
NLGD_Pos = NLGDs[_CT_Bias > 0]
NLGD_Neg = NLGDs[_CT_Bias <= 0]
print(len(NLGD_Pos), len(NLGD_Neg))
Total_G = len(NLGD_Pos) + len(NLGD_Neg)
N_Mut_Pos = int(NLGD_Pos.sum())
Total_Mut = int(N_Mut_Pos + int(NLGD_Neg.sum()))
print(N_Mut_Pos, Total_Mut, N_Mut_Pos/Total_Mut / (len(NLGD_Pos)/Total_G))
print(scipy.stats.binom_test(N_Mut_Pos, Total_Mut, p=len(NLGD_Pos)/len(NLGDs)))
#print("P_binom=")

# %% hidden=true
## Rank vs Bias
CT = "Medium spiny neuron"
CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
Gene_Idx = Spark_Meta_2stage.head(61)["EntrezID"].values
dat = CT_Z2_MAT_HC.loc[Gene_Idx, CT_Idx].mean(axis=1)
plt.plot(np.arange(len(Gene_Idx)), dat,)

# %% hidden=true
CT = "MGE interneuron"
CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
Gene_Idx = Spark_Meta_2stage.head(61)["EntrezID"].values
dat = CT_Z2_MAT_HC.loc[Gene_Idx, CT_Idx].mean(axis=1)
plt.plot(np.arange(len(Gene_Idx)), dat,)

# %% hidden=true
Z2Mat = pd.read_csv("../../ASD_Circuits/dat/allen-mouse-exp//AllenMouseBrain_Z2bias.csv", index_col=0)
allen_mouse_genes = Z2Mat.index.values
BGMR = pd.read_csv("~/Work/Resources/MutationRate_20170710_rate.txt", delimiter="\t")
BGMR["Entrez"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in BGMR["GeneName"].values]
BGMR = BGMR[BGMR["Entrez"].isin(allen_mouse_genes)]
BGMR.index = BGMR["Entrez"].values

# %% hidden=true
NLGDs = []
NDmis = []
Bias_CT = {}
for i, row in Spark_Meta_2stage.head(60).iterrows():
    gene = int(row["EntrezID"])
    if gene in BGMR.index.values:
        rate1 = BGMR.loc[gene, "p_LGD"]*1e6
        rate2 = BGMR.loc[gene, "prevel_0.5"]*1e6
    else:
        rate1 = BGMR["p_LGD"].mean()*1e6
        rate2 = BGMR["prevel_0.5"].mean()*1e6
    NLGDs.append(row["AutismMerged_LoF"]/rate1)
    NDmis.append(row["AutismMerged_Dmis_REVEL0.5"]/rate2)
    for CT in Neurons:
        CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
        if CT not in Bias_CT:
            Bias_CT[CT] = []
        Bias_CT[CT].append(CT_Z2_MAT_HC.loc[gene, CT_Idx].mean())
NLGDs = np.array(NLGDs)
NDmis = np.array(NDmis)

# %% hidden=true
NLGDs = []
NDmis = []
Bias_CT = {}
for i, row in Spark_Meta_2stage.head(159).iterrows():
    gene = int(row["EntrezID"])

    NLGDs.append(row["AutismMerged_LoF"])
    NDmis.append(row["AutismMerged_Dmis_REVEL0.5"])
    for CT in Neurons:
        CT_Idx = Anno[Anno["Supercluster"]==CT].index.values
        if CT not in Bias_CT:
            Bias_CT[CT] = []
        Bias_CT[CT].append(CT_Z2_MAT_HC.loc[gene, CT_Idx].mean())
NLGDs = np.array(NLGDs)
NDmis = np.array(NDmis)

# %% hidden=true
CT = "Medium spiny neuron"
plt.scatter(Bias_CT[CT], NLGDs, label="LGD", color="red", zorder=2)
plt.scatter(Bias_CT[CT], NDmis, label="Dmis", color="blue", zorder=1)
plt.xlabel("Cell Type Bias")
plt.ylabel("Number of Mutations")
r1, p1 = spearmanr(Bias_CT[CT], NLGDs)
r2, p2 = spearmanr(Bias_CT[CT], NDmis)
plt.title("%s \nLGD Spearman's R=%.2f P=%.1e\nDmis Spearman's R=%.2f P=%.1e"%(CT, r1, p1, r2, p2))
plt.legend()
plt.show()

# %% hidden=true
for CT in Neurons:
    print(CT)
    plt.scatter(Bias_CT[CT], NLGDs, label="LGD", color="red", zorder=2)
    plt.scatter(Bias_CT[CT], NDmis, label="Dmis", color="blue", zorder=1)
    plt.xlabel("Cell Type Bias")
    plt.ylabel("Number of Mutations")
    r1, p1 = spearmanr(Bias_CT[CT], NLGDs)
    r2, p2 = spearmanr(Bias_CT[CT], NDmis)
    plt.title("%s \nLGD Spearman's R=%.2f P=%.1e\nDmis Spearman's R=%.2f P=%.1e"%(CT, r1, p1, r2, p2))
    plt.legend()
    plt.show()

# %% hidden=true

# %% [markdown] heading_collapsed=true hidden=true
# ## High IQ vs Low IQ

# %% hidden=true
Spark_Denovo = pd.read_excel("../dat/41588_2022_1148_MOESM4_ESM.xlsx",
                           skiprows=2, sheet_name="Table S7")
Spark_Denovo = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]!="."]
#Spark_Denovo_ExomeWide = Spark_Denovo[Spark_Denovo[
#    "pDenovoWEST_Meta"]<=1.3e-6]
Spark_Denovo_ExomeWide = Spark_Denovo[Spark_Denovo[
    "pDenovoWEST_Meta"]<=0.01]
Spark_Denovo_ExomeWide.shape

# %% hidden=true
Mut_n_IQ = pd.read_csv("../dat/ASD_IQ_Mut.csv")
Mut_n_IQ = Mut_n_IQ[Mut_n_IQ["Entrez"].isin(HumanCT_Z2_HCT.index.values)]

# %% hidden=true
top_Genes = Spark_Denovo_ExomeWide["HGNC"].values
Mut_n_IQ_conf = Mut_n_IQ[Mut_n_IQ["HGNC"].isin(top_Genes)]
Mut_n_IQ_conf.shape

# %% hidden=true
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"]>70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"]<=70]

HighIQMuts = Mut_n_IQ_conf[Mut_n_IQ_conf["IQ"]>70]
LowIQMuts = Mut_n_IQ_conf[Mut_n_IQ_conf["IQ"]<=70]

# %%
LowIQMuts


# %%
def CountMut(DF):
    N_LGD, N_mis, N_Dmis, N_syn = 0,0,0,0
    for i, row in DF.iterrows():
        GeneEff = row["GeneEff"].split(";")[0]
        if GeneEff in ["frameshift", "splice_acceptor", "splice_donor", "start_lost", "stop_gained", "stop_lost"]:
            N_LGD += 1
        elif GeneEff == "missense":
            N_mis += 1
            row["REVEL"] = row["REVEL"].split(";")[0]
            if row["REVEL"] != ".":
                if float(row["REVEL"]) > 0.5:
                    N_Dmis += 1
        elif GeneEff == "synonymous":
            N_syn += 1
    return N_LGD, N_mis, N_Dmis, N_syn
def Mut2GeneDF(MutDF, PPVs = [0.554, 0.333, 0.138, 0.130], LGD=True, Dmis=True):
    genes = np.array(list(set(MutDF["HGNC"].values)))
    dat = []
    gene2MutN = {}
    for g in genes:
        try:
            Entrez = int(GeneSymbol2Entrez[g])
        except:
            Entrez = -1
            continue
        Muts = MutDF[MutDF["HGNC"]==g]
        try:
            pLI = float(Muts["ExACpLI"].values[0])
        except:
            pLI = 0
        N_LGD, N_Mis, N_Dmis, N_Syn = CountMut(Muts)
        if not LGD:
            N_LGD = 0
        if not Dmis:
            N_Dmis = 0
        if pLI > 0.5:
            gene2MutN[Entrez] = N_LGD * PPVs[0] + N_Dmis * PPVs[1]
        else:
            gene2MutN[Entrez] = N_LGD * PPVs[2] + N_Dmis * PPVs[3]
    return gene2MutN


# %%
def Mut2GeneDF(MutDF, GeneSymbol2Entrez, pLI=True, LGD_weight_high=0.554, Dmis_weight_high=0.333, LGD_weight_low=0.138, Dmis_weight_low=0.130, LGD_weight_nopli=0.357, Dmis_weight_nopli=0.231):
    Select_Genes = np.array(list(set(MutDF["HGNC"].values)))
    dat = []
    gene2MutN = {}
    for g in Select_Genes:
        try:
            Entrez = int(GeneSymbol2Entrez[g])
        except:
            Entrez = -1
            continue
        Muts = MutDF[MutDF["HGNC"]==g]
        N_LGD, N_Mis, N_Dmis, N_Syn = CountMut(Muts)
        if pLI:
            try:
                pLI = float(Muts["ExACpLI"].values[0])
            except:
                pLI = 0.0
            if pLI >= 0.5:
                gene2MutN[Entrez] = N_LGD * LGD_weight_high + N_Dmis * Dmis_weight_high
            else:
                gene2MutN[Entrez] = N_LGD * LGD_weight_low + N_Dmis * Dmis_weight_low
        else:
            gene2MutN[Entrez] = N_LGD * LGD_weight_nopli + N_Dmis * Dmis_weight_nopli
    return gene2MutN


# %% hidden=true
HIQ_GW = Mut2GeneDF(HighIQMuts, GeneSymbol2Entrez)
LIQ_GW = Mut2GeneDF(LowIQMuts, GeneSymbol2Entrez)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(HIQ_GW, "{}/HIQ.top61.nopLI.2.gw".format(OutDIR))
Dict2Fil(LIQ_GW, "{}/LIQ.top61.nopLI.2.gw".format(OutDIR))

# %%
HIQ_GW = Mut2GeneDF(HighIQMuts, GeneSymbol2Entrez, pLI=False, LGD_weight_nopli=1, Dmis_weight_nopli=1)
LIQ_GW = Mut2GeneDF(LowIQMuts, GeneSymbol2Entrez, pLI=False, LGD_weight_nopli=1, Dmis_weight_nopli=1)
print(len(HIQ_GW), len(LIQ_GW))
OutDIR = "../dat/GeneWeights/"
Dict2Fil(HIQ_GW, "{}/HIQ.top100.nopLI.LGD_Dmis_SameWeight.gw".format(OutDIR))
Dict2Fil(LIQ_GW, "{}/LIQ.top100.nopLI.LGD_Dmis_SameWeight.gw".format(OutDIR))

# %%

# %%
# Get genes with mutations in HighIQMuts
genes_with_mutations = set(HighIQMuts_ALL["HGNC"].unique())

# Sort Spark_Denovo by pDenovoWEST_Meta and get top genes that have mutations
top_genes = Spark_Denovo.sort_values("pDenovoWEST_Meta")[
    Spark_Denovo["HGNC"].isin(genes_with_mutations)
].head(61)

print(f"Number of top genes: {len(top_genes)}")
print("\nTop genes sorted by pDenovoWEST_Meta:")
print(top_genes[["HGNC", "pDenovoWEST_Meta"]].to_string())

# %%
top_genes

# %%
(Spark_Denovo["pDenovoWEST_Meta"]<1e-3).sum()

# %% [markdown]
# # 22q11.2 

# %%
Gene22q = [8214, 5625, 9993, 23617, 2928, 6576, 8218, 7290, 64976, 128977,
           8318, 7122, 2812, 6899, 54584, 10587, 1312, 421, 128989, 54487,
           27037, 5902, 29801, 388849, 65078, 85359, 373856, 7625, 91179, 84861,
           51586, 5297, 3053, 9342, 1399, 150209, 8216, 80764, 9127, 6545, 400891,
           8220, 7353, 5413, 79680, 728418]
X22q_GW = dict(zip(Gene22q, [1]*len(Gene22q)))

X22q_mouse_model_genes_corrected = ["DGCR6", "PRODH", "DGCR2", "TSSK2", "ESS2", "GSC2", "SLC25A1",
                         "CLTCL1", "HIRA", "MRPL40", "UFD1", "CDC45", "CLDN5", "SEPTIN5",
                         "GP1BB", "TBX1", "GNB1L", "TXNRD2", "COMT", "ARVCF", "CD38", 
                         "DGCR8", "TRMT2A", "RANBP1", "ZDHHC8", "RTN4R", "PRODH", "DGCR6L"]
X22q_mouse_model_entrez = [GeneSymbol2Entrez[x] for x in X22q_mouse_model_genes_corrected]
X22q_V2_GW = dict(zip(X22q_mouse_model_entrez, [1]*len(X22q_mouse_model_entrez)))

Dict2Fil(X22q_GW, "../dat/GeneWeights/X22q.gw.csv")
Dict2Fil(X22q_V2_GW, "../dat/GeneWeights/X22q.mousemodel.gw.csv")


# %%
print(len(Gene22q), len(X22q_mouse_model_entrez))

# %% [markdown] heading_collapsed=true
# # Bipolar

# %% hidden=true
df = pd.read_csv("../../data/Bipolar/Bipolar_gene_res.tsv", sep="\t")
df1 = df[df["group"]=="Bipolar Disorder"]
df2 = df[df["group"]=="Bipolar Disorder (including Schizoaffective)"]
print(df.shape, df2.shape)

# %% hidden=true
# Annotate genes and combine p from dmis and PTV
for i, row in df1.iterrows():
    gene_id = row["gene_id"]
    gene_entrez = ENSID2Entrez.get(gene_id, 0)
    gene_symbol = Entrez2Symbol.get(gene_entrez, ".")
    ncase = row["n_cases"]
    ncont = row["n_controls"]
    df1.loc[i, "EntrezID"] = gene_entrez
    df1.loc[i, "Symbol"] = gene_symbol
    #print(gene_entrez, gene_symbol)
    _, fisherP = scipy.stats.combine_pvalues([row["damaging_missense_fisher_gnom_non_psych_pval"], 
                                    row["ptv_fisher_gnom_non_psych_pval"]])
    #print(fisherP)
    Eff_LGD = row["ptv_case_count"] - row["ptv_control_count"]/ncont*ncase
    Eff_Dmis = row["damaging_missense_case_count"] - row["damaging_missense_control_count"]/ncont*ncase
    df1.loc[i, "Pcombine"] = fisherP
    df1.loc[i, "Effect.LGD"] = max(0, Eff_LGD)
    df1.loc[i, "Effect.Dmis"] = max(0, Eff_Dmis)
df1 = df1.sort_values("Pcombine", ascending=True)
df1 = df1[(df1["Effect.LGD"]>0) | (df1["Effect.Dmis"]>0)]
df1["EntrezID"] = [int(x) for x in df1["EntrezID"].values]
df1.to_csv("../dat/bipolar.genes.modified.csv")

# %% hidden=true
# Annotate genes and combine p from dmis and PTV
for i, row in df2.iterrows():
    gene_id = row["gene_id"]
    gene_entrez = ENSID2Entrez.get(gene_id, 0)
    gene_symbol = Entrez2Symbol.get(gene_entrez, ".")
    ncase = row["n_cases"]
    ncont = row["n_controls"]
    df2.loc[i, "EntrezID"] = gene_entrez
    df2.loc[i, "Symbol"] = gene_symbol
    #print(gene_entrez, gene_symbol)
    _, fisherP = scipy.stats.combine_pvalues([row["damaging_missense_fisher_gnom_non_psych_pval"], 
                                    row["ptv_fisher_gnom_non_psych_pval"]])
    #print(fisherP)
    Eff_LGD = row["ptv_case_count"] - row["ptv_control_count"]/ncont*ncase
    Eff_Dmis = row["damaging_missense_case_count"] - row["damaging_missense_control_count"]/ncont*ncase
    df2.loc[i, "Pcombine"] = fisherP
    df2.loc[i, "Effect.LGD"] = max(0, Eff_LGD)
    df2.loc[i, "Effect.Dmis"] = max(0, Eff_Dmis)
df2 = df2.sort_values("Pcombine", ascending=True)
df2 = df2[(df2["Effect.LGD"]>0) | (df2["Effect.Dmis"]>0)]
df2["EntrezID"] = [int(x) for x in df2["EntrezID"].values]
df2.to_csv("../dat/bipolar.wscz.genes.modified.csv")

# %% hidden=true
BP_topG = df2.head(61)["Symbol"].values
BP_top100G = df2.head(100)["Symbol"].values
BP_top200G = df2.head(200)["Symbol"].values

# %% hidden=true
BP_topG

# %% hidden=true

# %% hidden=true
ASD_SCZ_DF = pd.read_csv("GeneFuncs_ASD_SCZ-NetBag.tsv", sep="\t")

# %% hidden=true
SCZ_Symbols = ASD_SCZ_DF[ASD_SCZ_DF["SCZ"]==1]["Gene"].values
ASD_Symbols = ASD_SCZ_DF[ASD_SCZ_DF["ASD"]==1]["Gene"].values

# %% hidden=true
set(SCZ_Symbols).intersection(set(BP_topG))

# %% hidden=true
set(ASD_Symbols).intersection(set(BP_topG))

# %% hidden=true
BP_topG

# %% hidden=true

# %% hidden=true
CT_Z2_MAT_HC = pd.read_csv("../dat/HumanCellType.AllCell.HCT.Z2bias.entrez.csv", index_col=0)
CT_Z2_MAT_HC.columns = [int(x) for x in CT_Z2_MAT_HC.columns.values]

# %% hidden=true
Bipolar_GW1 = Aggregate_Gene_Weights_Bipolar(df1.head(61), CT_Z2_MAT_HC.index.values)
Bipolar_GW2 = Aggregate_Gene_Weights_Bipolar(df2.head(61), CT_Z2_MAT_HC.index.values)

Bipolar_100GW = Aggregate_Gene_Weights_Bipolar(df1.head(100), CT_Z2_MAT_HC.index.values)
Bipolar_200GW2 = Aggregate_Gene_Weights_Bipolar(df1.head(200), CT_Z2_MAT_HC.index.values)

# %% hidden=true
Dict2Fil(Bipolar_GW1, "../dat/Bipolar.top61.gw.csv")
Dict2Fil(Bipolar_GW2, "../dat/Bipolar.2.top61.gw.csv")
Dict2Fil(Bipolar_100GW, "../dat/Bipolar.top100.gw.csv")
Dict2Fil(Bipolar_200GW2, "../dat/Bipolar.top200.gw.csv")

# %% hidden=true

# %% hidden=true

# %% hidden=true

# %% [markdown] heading_collapsed=true
# # NDD

# %%
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv")
BGMR.head(2)

# %% hidden=true
df = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
df = df.sort_values("denovoWEST_p_full")
hc_df = df[df["denovoWEST_p_full"]<=0.05/18762]
entrez_ids = [int(GeneSymbol2Entrez.get(x, -1)) for x in hc_df["symbol"].values]
hc_df["EntrezID"] = entrez_ids
hc_df.shape

# %%
#DDD_hc_GW = dict(zip(DDD_hc_entrez, [1]*len(DDD_hc_entrez)))
# DDD_hc_GW = Aggregate_Gene_Weights_NDD(hc_df, out="../dat/GeneWeights/DDD.hc.gw.csv")
# NDD_top61_DF = hc_df.head(61)
# DDD_top61_GW = Aggregate_Gene_Weights_NDD(NDD_top61_DF, out="../dat/GeneWeights/DDD.top61.gw.csv")
#Dict2Fil(DDD_top61_GW, "../dat/GeneWeights/DDD.top61.gw.csv")

DDD_hc_GW = Aggregate_Gene_Weights_NDD(hc_df, BGMR, wLGD=1, wMis=1, out="../dat/GeneWeights/DDD.hc.gw")
NDD_top61_DF = hc_df.head(61)
DDD_top61_GW = Aggregate_Gene_Weights_NDD(NDD_top61_DF, BGMR, wLGD=1, wMis=1, out="../dat/GeneWeights/DDD.top61.gw")

# %%
NDD_top61_genes = hc_df.head(61)["EntrezID"].values
NDD_top297_genes = hc_df.head(297)["EntrezID"].values

# %%
NDD_top61_genes

# %%
SCZ_top61 = set(int(x) for x in GeneDF.head(61).index.values)  # SCZ top 61 Entrez gene IDs
NDD_top61_genes_set = set(int(x) for x in NDD_top61_genes)
NDD_top297_genes_set = set(int(x) for x in NDD_top297_genes)

overlap_61 = SCZ_top61.intersection(NDD_top61_genes_set)
overlap_297 = SCZ_top61.intersection(NDD_top297_genes_set)

print(f"Number of SCZ top 61 overlapping NDD top 61: {len(overlap_61)}")
print(f"Number of SCZ top 61 overlapping NDD top 297: {len(overlap_297)}")

# Create list(s) for SCZ top 61 excluding overlaps with NDD top 61 and NDD top 297
SCZ_top61_ex_NDD61 = SCZ_top61 - NDD_top61_genes_set
SCZ_top61_ex_NDD297 = SCZ_top61 - NDD_top297_genes_set

print(f"Number of SCZ top 61 excluding NDD top 61: {len(SCZ_top61_ex_NDD61)}")
print(f"Number of SCZ top 61 excluding NDD top 297: {len(SCZ_top61_ex_NDD297)}")


# %%
tmp_GeneDF = GeneDF.loc[list(SCZ_top61_ex_NDD61), :]
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(tmp_GeneDF, HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.ExlNDD61.gw".format(OutDIR))

tmp_GeneDF = GeneDF.loc[list(SCZ_top61_ex_NDD297), :]
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(tmp_GeneDF, HumanCT_Z2_HCT.index.values, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0)
OutDIR = "../dat/GeneWeights/"
Dict2Fil(SCZ_61GW, "{}/SCZ.top61.ExlNDD297.gw".format(OutDIR))

# %% [markdown] heading_collapsed=true
# # Epi25

# %% hidden=true
Epi25DF = pd.read_excel("../../data/Epilepsy/1-s2.0-S0002929721001403-mmc2.xlsx", 
                        sheet_name="S17 All Epi Deleterious")

# %% hidden=true
#Epi25DF.columns.values

# %% hidden=true
# 'Missense.Case', 'PTV.Case', 'Missense.Ctrl', 'PTV.Ctrl'

# %% hidden=true
Epi25DF = Epi25DF[Epi25DF["estimate"]>1]
for i, row in Epi25DF.iterrows():
    symbol = Epi25DF.loc[i, "gene"].strip("''")
    Epi25DF.loc[i, "gene"] = symbol
    try:
        Epi25DF.loc[i, "Entrez"] = int(GeneSymbol2Entrez[symbol])
    except:
        Epi25DF.loc[i, "Entrez"] = None
    try:
        Epi25DF.loc[i, "pLI"] = ExAC_pLI.loc[symbol, "pLI"]
    except:
        Epi25DF.loc[i, "pLI"] = 0
    Epi25DF.loc[i, "nLGD"] = ModifyMutCount(Epi25DF.loc[i, "PTV.Case"], 
                                    Epi25DF.loc[i, "PTV.Ctrl"], 0, CaseN=13487, ContN=15678) 
    Epi25DF.loc[i, "nMis"] = ModifyMutCount(Epi25DF.loc[i, "Missense.Case"], 
                                    Epi25DF.loc[i, "Missense.Ctrl"], 0, CaseN=13487, ContN=15678) 
    
Epi25DF = Epi25DF.dropna(subset="Entrez")
Epi25DF_Filt = Epi25DF[["gene", "Entrez", "pLI", "nLGD", "nMis"]]
Epi25DF_Filt["Entrez"] = [int(x) for x in Epi25DF_Filt["Entrez"].values]
Epi25DF_Filt = Epi25DF_Filt.set_index("Entrez")
Epi25DF_Filt.shape

# %% hidden=true
Epi25DF_Filt.head(10)

# %% hidden=true
Epi25_GW = Aggregate_Gene_Weights_Epi25(Epi25DF_Filt, ExpZ2Mat.index.values)

# %% hidden=true
Dict2Fil(Epi25_GW, "../dat3/Epi25.gw.csv")

# %% hidden=true
Epi25DF = pd.read_excel("../../data/Epilepsy/1-s2.0-S0002929721001403-mmc2.xlsx", 
                        sheet_name="S17 All Epi Deleterious")
Epi25DF_Neg = Epi25DF[Epi25DF["estimate"]<1]
for i, row in Epi25DF_Neg.iterrows():
    symbol = Epi25DF_Neg.loc[i, "gene"].strip("''")
    Epi25DF_Neg.loc[i, "gene"] = symbol
    try:
        Epi25DF_Neg.loc[i, "Entrez"] = int(GeneSymbol2Entrez[symbol])
    except:
        Epi25DF_Neg.loc[i, "Entrez"] = None
    try:
        Epi25DF_Neg.loc[i, "pLI"] = ExAC_pLI.loc[symbol, "pLI"]
    except:
        Epi25DF_Neg.loc[i, "pLI"] = 0
    Epi25DF_Neg.loc[i, "nLGD"] = ModifyMutCount(Epi25DF_Neg.loc[i, "PTV.Case"], 
                                    Epi25DF.loc[i, "PTV.Ctrl"], -1000, CaseN=13487, ContN=15678) 
    Epi25DF_Neg.loc[i, "nMis"] = ModifyMutCount(Epi25DF_Neg.loc[i, "Missense.Case"], 
                                    Epi25DF.loc[i, "Missense.Ctrl"], -1000, CaseN=13487, ContN=15678) 
    
Epi25DF_Neg = Epi25DF_Neg.dropna(subset="Entrez")
Epi25DF_Neg_Filt = Epi25DF_Neg[["gene", "Entrez", "pLI", "nLGD", "nMis"]]
Epi25DF_Neg_Filt["Entrez"] = [int(x) for x in Epi25DF_Neg_Filt["Entrez"].values]
Epi25DF_Neg_Filt = Epi25DF_Neg_Filt.set_index("Entrez")
Epi25DF_Neg_Filt.shape

# %% hidden=true
Epi25DF_Neg_Filt

# %% hidden=true
Epi25DF_Protect_gw = dict(zip(Epi25DF_Neg_Filt.index.values, [1]*len(Epi25DF_Neg_Filt.index.values)))
Dict2Fil(Epi25DF_Protect_gw, "../dat3/Epi25.protect.gw.csv")

# %% hidden=true
Epi25DF_Protect_gw

# %% [markdown] heading_collapsed=true hidden=true
# #### subtypes

# %% hidden=true
DEE_DF = pd.read_excel("../../data/Epilepsy/1-s2.0-S0002929721001403-mmc2.xlsx", 
                        sheet_name="S09 DEE Deleterious")
GGE_DF = pd.read_excel("../../data/Epilepsy/1-s2.0-S0002929721001403-mmc2.xlsx", 
                        sheet_name="S12 GGE Deleterious")
NAFE_DF = pd.read_excel("../../data/Epilepsy/1-s2.0-S0002929721001403-mmc2.xlsx", 
                        sheet_name="S14 NAFE Deleterious")


# %% hidden=true
def Filt_Epilepsy_DF(Epi25DF, Ncase, Ncontrol):
    Epi25DF = Epi25DF[Epi25DF["estimate"]>1]
    for i, row in Epi25DF.iterrows():
        symbol = Epi25DF.loc[i, "gene"].strip("''")
        Epi25DF.loc[i, "gene"] = symbol
        try:
            Epi25DF.loc[i, "Entrez"] = int(GeneSymbol2Entrez[symbol])
        except:
            Epi25DF.loc[i, "Entrez"] = None
        try:
            Epi25DF.loc[i, "pLI"] = ExAC_pLI.loc[symbol, "pLI"]
        except:
            Epi25DF.loc[i, "pLI"] = 0
        Epi25DF.loc[i, "nLGD"] = ModifyMutCount(Epi25DF.loc[i, "PTV.Case"], 
                                        Epi25DF.loc[i, "PTV.Ctrl"], 0, CaseN=Ncase, ContN=Ncontrol) 
        Epi25DF.loc[i, "nMis"] = ModifyMutCount(Epi25DF.loc[i, "Missense.Case"], 
                                        Epi25DF.loc[i, "Missense.Ctrl"], 0, CaseN=Ncase, ContN=Ncontrol) 

    Epi25DF = Epi25DF.dropna(subset="Entrez")
    Epi25DF_Filt = Epi25DF[["gene", "Entrez", "pLI", "nLGD", "nMis"]]
    Epi25DF_Filt["Entrez"] = [int(x) for x in Epi25DF_Filt["Entrez"].values]
    Epi25DF_Filt = Epi25DF_Filt.set_index("Entrez")
    return Epi25DF_Filt


# %% hidden=true
DEE_DF_Filt = Filt_Epilepsy_DF(DEE_DF, 1835, 13978)
GEE_DF_Filt = Filt_Epilepsy_DF(GGE_DF, 5303, 15677)
NAFE_DF_Filt = Filt_Epilepsy_DF(NAFE_DF, 6439, 15678)

# %% hidden=true
DEE_GW = Aggregate_Gene_Weights_Epi25(DEE_DF_Filt, ExpZ2Mat.index.values)
Dict2Fil(DEE_GW, "../dat3/Epi25.DEE.gw.csv")

GEE_GW = Aggregate_Gene_Weights_Epi25(GEE_DF_Filt, ExpZ2Mat.index.values)
Dict2Fil(GEE_GW, "../dat3/Epi25.GEE.gw.csv")

NAFE_GW = Aggregate_Gene_Weights_Epi25(NAFE_DF_Filt, ExpZ2Mat.index.values)
Dict2Fil(NAFE_GW, "../dat3/Epi25.NAFE.gw.csv")

# %% hidden=true
DEE_GW

# %% hidden=true
GEE_GW

# %% [markdown]
# # UKBB Cog function

# %% [markdown]
# Paper: https://www.nature.com/articles/s41588-023-01398-8
#
# educational attainment (EDU); 
#
# reaction time (RT); 
#
# verbal-numerical reasoning (VNR)

# %%
CogDF = pd.read_excel("../dat/41588_2023_1398_MOESM3_ESM.xlsx", sheet_name="Table S4")
CogDF = CogDF[CogDF["POPULATION"]=="EUR"]


# %%
def DF2GW_UKBB(DF, Ncarrier=False):
    res = {}
    for i, row in DF.iterrows():
        entrez = int(GeneSymbol2Entrez.get(row["GENE"], 0))
        if entrez == 0 :
            continue
        if Ncarrier:
            res[entrez] = abs(row["BETA"]) * row["NCARRIER"]
        else:
            res[entrez] = abs(row["BETA"]) 
    return res


# %%
topN = 61
VNR_DF = CogDF[CogDF["PHENOTYPE"]=="VNR"].sort_values("P")
VNR_Pos_DF = VNR_DF[VNR_DF["BETA"]>0].head(topN)
VNR_Neg_DF = VNR_DF[VNR_DF["BETA"]<0].head(topN)

VNR_Pos_GW = DF2GW_UKBB(VNR_Pos_DF)
VNR_Neg_GW = DF2GW_UKBB(VNR_Neg_DF)
Dict2Fil(VNR_Pos_GW, "{}/UKBB_VNR_Pos_GW_{}.csv".format(OutDIR, topN))
Dict2Fil(VNR_Neg_GW, "{}/UKBB_VNR_Neg_GW_{}.csv".format(OutDIR, topN))

EDU_DF = CogDF[CogDF["PHENOTYPE"]=="EDU"].sort_values("P")
EDU_Pos_DF = EDU_DF[EDU_DF["BETA"]>0].head(topN)
EDU_Neg_DF = EDU_DF[EDU_DF["BETA"]<0].head(topN)

EDU_Pos_GW = DF2GW_UKBB(EDU_Pos_DF)
EDU_Neg_GW = DF2GW_UKBB(EDU_Neg_DF)
Dict2Fil(EDU_Pos_GW, "{}/UKBB_EDU_Pos_GW_{}.csv".format(OutDIR, topN))
Dict2Fil(EDU_Neg_GW, "{}/UKBB_EDU_Neg_GW_{}.csv".format(OutDIR, topN))

RT_DF = CogDF[CogDF["PHENOTYPE"]=="RT"].sort_values("P")
RT_Pos_DF = RT_DF[RT_DF["BETA"]>0].head(topN)
RT_Neg_DF = RT_DF[RT_DF["BETA"]<0].head(topN)

RT_Pos_GW = DF2GW_UKBB(RT_Pos_DF)
RT_Neg_GW = DF2GW_UKBB(RT_Neg_DF)
Dict2Fil(RT_Pos_GW, "{}/UKBB_RT_Pos_GW_{}.csv".format(OutDIR, topN))
Dict2Fil(RT_Neg_GW, "{}/UKBB_RT_Neg_GW_{}.csv".format(OutDIR, topN))

# %%
VNR_NoEff_DF = VNR_DF.tail(61)
VNR_NoEff_GW = DF2GW_UKBB(VNR_NoEff_DF)
Dict2Fil(VNR_NoEff_GW, "{}/UKBB_VNR_NoEff_GW_{}.csv".format(OutDIR, topN))

# %%
for k, v in VNR_NoEff_GW.items():
    print(Entrez2Symbol[k])

# %%
