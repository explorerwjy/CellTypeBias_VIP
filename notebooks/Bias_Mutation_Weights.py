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

# %% [markdown]
# # Gene Weight Generation
#
# This notebook generates ALL gene weight files used by the bias pipeline.
# Each section produces `.gw` or `.csv` files in `dat/GeneWeights/`.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
from pathlib import Path
from collections import OrderedDict
import yaml

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *
import os

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# Expression matrix (for valid gene filtering)
HumanCT_Z2_HCT = pd.read_csv(
    "/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv",
    index_col=0,
)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
valid_genes = HumanCT_Z2_HCT.index.values

# BGMR — background mutation rate, indexed by Entrez ID
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv")
BGMR = BGMR.set_index("entrez_id")
BGMR.index = BGMR.index.astype(int)

# Output directory
GW_DIR = PROJ_DIR / "dat" / "GeneWeights"

print(f"Expression matrix: {HumanCT_Z2_HCT.shape}")
print(f"BGMR genes: {len(BGMR)}")
print(f"Output dir: {GW_DIR}")

# %% [markdown]
# # SCZ

# %%
# Load SCZ gene results from Singh et al. 2022
GeneDF = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41586_2022_4556_MOESM3_ESM.xlsx"),
    sheet_name="Table S5 - Gene Results",
)
ExAC_pLI = pd.read_csv(
    "/home/jw3514/Work/Resources/gnomad.v2.1.1.lof_metrics.by_gene.txt",
    sep="\t",
    index_col="gene",
)

Ncase = 24248
Nctrl = 97322

# %% [markdown]
# ## Modify Mutation Counts

# %%
def oddsratio(NcaseMut, NctrlMut, dnvCount, Ncase=24248, Nctrl=97322):
    if dnvCount != dnvCount:
        dnvCount = 0
    NcaseMut += 1  # pseudo count to prevent inf
    NctrlMut += 1
    AD = (NcaseMut) * (Nctrl - NctrlMut)
    BC = (NctrlMut) * (Ncase - NcaseMut)
    return AD / BC + dnvCount


def Penetrance(NcaseMut, NctrlMut, Ncase=24248, Nctrl=97322, prevelence=0.45 / 100):
    NcaseMut += 1
    NctrlMut += 1
    Ntotal_ = Ncase / prevelence
    Nctrl_ = Ntotal_ - Ncase
    NctrlMut_ = Nctrl_ * (NctrlMut / Nctrl)
    p = NcaseMut / (NcaseMut + NctrlMut_)
    return p


def ModifyMutCount(CaseCount, ContCount, dnvCount, CaseN=24248, ContN=97322):
    if isinstance(dnvCount, float) and np.isnan(dnvCount):
        dnvCount = 0
    return CaseCount - ContCount / ContN * CaseN + dnvCount

# %%
for i, row in GeneDF.iterrows():
    symbol = row["Gene Symbol"]

    try:
        GeneDF.loc[i, "Entrez"] = int(GeneSymbol2Entrez[symbol])
    except Exception:
        GeneDF.loc[i, "Entrez"] = None

    try:
        GeneDF.loc[i, "pLI"] = ExAC_pLI.loc[symbol, "pLI"]
    except Exception:
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
GeneDF.to_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"))
print(f"SCZ GeneDF: {GeneDF.shape}")

# %%
# Reload to ensure clean state
GeneDF = pd.read_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"), index_col=0)

# %% [markdown]
# ## SCZ Protective Gene Weights

# %%
SCZ_Protect = GeneDF[GeneDF["nLGD"] < 0].head(61)
SCZ_61_protect_GW = Aggregate_Gene_Weights_SCZ_Daly(
    SCZ_Protect, valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
SCZ_61_protect_GW = {k: abs(v) for k, v in SCZ_61_protect_GW.items()}
Dict2Fil(SCZ_61_protect_GW, f"{GW_DIR}/SCZ.top61.protect.gw")
print(f"SCZ protective: {len(SCZ_61_protect_GW)} genes")

# %% [markdown]
# ## SCZ Risk Gene Weights

# %%
# Main pipeline file: top 61, LGD+Dmis same weight, exclude Mis2
SCZ_61GW = Aggregate_Gene_Weights_SCZ_Daly(
    GeneDF.head(61), valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_61GW, f"{GW_DIR}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw")

# Also save with Mis2 included (same weights, just a different filename)
Dict2Fil(SCZ_61GW, f"{GW_DIR}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.gw")

print(f"SCZ top 61: {len(SCZ_61GW)} genes")

# %%
# Larger gene sets for sensitivity analysis
GeneDF_pos = GeneDF[GeneDF["nLGD"] > 1]

SCZ_100GW = Aggregate_Gene_Weights_SCZ_Daly(
    GeneDF_pos.head(100), valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_100GW, f"{GW_DIR}/SCZ.top100.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw")

SCZ_200GW = Aggregate_Gene_Weights_SCZ_Daly(
    GeneDF_pos.head(200), valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_200GW, f"{GW_DIR}/SCZ.top200.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw")

SCZ_500GW = Aggregate_Gene_Weights_SCZ_Daly(
    GeneDF_pos.head(500), valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_500GW, f"{GW_DIR}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw")

print(f"SCZ top 100/200/500: {len(SCZ_100GW)}/{len(SCZ_200GW)}/{len(SCZ_500GW)} genes")

# Store SCZ top 61 Entrez IDs for later exclusion analysis
SCZ_top61 = set(int(x) for x in GeneDF.head(61).index.values)

# %% [markdown]
# # ASD

# %% [markdown]
# ## Load ASD data

# %%
# Load Spark DenovoWEST meta-analysis data (sorted by p-value)
Spark_Denovo = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41588_2022_1148_MOESM4_ESM.xlsx"),
    skiprows=2,
    sheet_name="Table S7",
)
Spark_Denovo = Spark_Denovo[Spark_Denovo["pDenovoWEST_Meta"] != "."]
Spark_Denovo["pDenovoWEST_Meta"] = Spark_Denovo["pDenovoWEST_Meta"].astype(float)
Spark_Denovo = Spark_Denovo.sort_values("pDenovoWEST_Meta").reset_index(drop=True)
print(f"Spark Denovo genes: {Spark_Denovo.shape[0]}")

# %%
# Load IQ mutation data
Mut_n_IQ = pd.read_csv(str(PROJ_DIR / "dat/ASD_IQ_Mut.csv"))
Mut_n_IQ = Mut_n_IQ[Mut_n_IQ["Entrez"].isin(valid_genes)]
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] > 70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] <= 70]
print(f"HIQ mutations: {len(HighIQMuts_ALL)}, LIQ mutations: {len(LowIQMuts_ALL)}")

# %% [markdown]
# ## BGMR-corrected ASD gene weights

# %%
def compute_asd_bgmr_weights(spark_sub, count_col, BGMR, N_proband, n_top=61,
                              lof_col="AutismMerged_LoF", dmis_col="AutismMerged_Dmis_REVEL0.5"):
    """
    Compute BGMR-corrected ASD gene weights using LoF + Dmis (REVEL>0.5).

    1. Take top n_top genes by p-value (spark_sub must be pre-sorted)
    2. Filter to genes with mutations in the specified IQ group (count_col > 0)
    3. Compute weight = (observed_LoF - expected_LoF) + (observed_Dmis - expected_Dmis)
       Uses BGMR prevel_0.5 for expected Dmis rate (matches REVEL>0.5 filtering)
    4. Clip negative weights to 0 (edge cases where gene is significant by pattern not count)

    Parameters
    ----------
    spark_sub : DataFrame indexed by EntrezID, sorted by p-value
    count_col : str or None — column for filtering (e.g. 'HIQ_counts'); None = no filter
    BGMR : DataFrame indexed by entrez_id with p_LGD, prevel_0.5 columns
    N_proband : int — cohort size for expected mutation calculation
    n_top : int — number of top genes by p-value to consider
    lof_col : str — column for observed LoF counts
    dmis_col : str — column for observed Dmis (REVEL>0.5) counts

    Returns OrderedDict in p-value order.
    """
    top = spark_sub.head(n_top)
    if count_col is not None:
        top = top[top[count_col] > 0]
    gw = OrderedDict()
    for g in top.index:
        if g not in BGMR.index:
            continue
        nLGD = float(top.loc[g, lof_col])
        nDmis = float(top.loc[g, dmis_col])
        exp_LGD = BGMR.loc[g, "p_LGD"] * 2 * N_proband
        exp_Dmis = BGMR.loc[g, "prevel_0.5"] * 2 * N_proband
        weight = (nLGD - exp_LGD) + (nDmis - exp_Dmis)
        gw[g] = max(weight, 0)  # clip edge cases
    return gw

# %%
# Cohort sizes for BGMR correction (WES probands with IQ measured)
# Derived from cross-referencing:
#   - ASC/SSC: 1-s2.0-S0092867419313984-mmc4.xlsx Phenotype sheet
#   - SPARK: core_descriptive_variables.csv
#   - WES sample list: ASD_Discov_Trios.txt + ASD_Rep_Trios.txt
N_ASD_ALL = 42607   # Full merged cohort (Fu et al. 2022)
N_ASD_HIQ = 3610    # WES probands with IQ > 70 (SSC/ASC: 3009, SPARK: 601)
N_ASD_LIQ = 1792    # WES probands with IQ <= 70 (SSC/ASC: 1494, SPARK: 298)
print(f"Cohort sizes: ASD All={N_ASD_ALL}, HIQ={N_ASD_HIQ}, LIQ={N_ASD_LIQ}")

# %%
# Build Spark_sub with IQ-specific mutation counts per gene
Spark_sub = Spark_Denovo[
    [
        "EntrezID",
        "HGNC",
        "pDenovoWEST_Meta",
        "AutismMerged_LoF",
        "AutismMerged_Dmis_REVEL0.5",
    ]
].copy()
Spark_sub = Spark_sub.set_index("EntrezID")

# Count IQ-specific mutations per gene (LGD and Dmis separately)
# Dmis = missense with REVEL > 0.5 (matching BGMR prevel_0.5 rate)
for i in Spark_sub.index:
    hiq_muts = HighIQMuts_ALL[HighIQMuts_ALL["Entrez"] == i]
    liq_muts = LowIQMuts_ALL[LowIQMuts_ALL["Entrez"] == i]
    Spark_sub.loc[i, "HIQ_counts"] = len(hiq_muts)
    Spark_sub.loc[i, "LIQ_counts"] = len(liq_muts)
    # IQ-specific LGD and Dmis (REVEL>0.5) counts for BGMR correction
    lgd_types = ["frameshift", "splice_acceptor", "splice_donor", "start_lost", "stop_gained", "stop_lost"]
    for prefix, muts in [("HIQ", hiq_muts), ("LIQ", liq_muts)]:
        eff = muts["GeneEff"].str.split(";").str[0]
        n_lgd = (eff.isin(lgd_types)).sum()
        # Dmis: missense with REVEL > 0.5
        is_mis = eff == "missense"
        revel = muts.loc[is_mis, "REVEL"].str.split(";").str[0]
        n_dmis = (revel.apply(lambda x: float(x) > 0.5 if x != "." else False)).sum()
        Spark_sub.loc[i, f"{prefix}_LGD"] = n_lgd
        Spark_sub.loc[i, f"{prefix}_Dmis"] = n_dmis

print(f"Spark_sub: {Spark_sub.shape}")
print(f"HIQ genes with mutations: {(Spark_sub['HIQ_counts'] > 0).sum()}")
print(f"LIQ genes with mutations: {(Spark_sub['LIQ_counts'] > 0).sum()}")

# %%
# ASD All (Spark Meta) — BGMR corrected (LoF + Dmis REVEL>0.5, full cohort)
asd_all_gw = compute_asd_bgmr_weights(
    Spark_sub, count_col=None, BGMR=BGMR, N_proband=N_ASD_ALL, n_top=61,
    lof_col="AutismMerged_LoF", dmis_col="AutismMerged_Dmis_REVEL0.5",
)
Dict2Fil(asd_all_gw, f"{GW_DIR}/Spark_Meta_EWS.GeneWeight.bgmr.csv")

# ASD All — raw (no BGMR)
asd_all_raw = OrderedDict()
for g in Spark_sub.head(61).index:
    nLGD = float(Spark_sub.loc[g, "AutismMerged_LoF"])
    nDmis = float(Spark_sub.loc[g, "AutismMerged_Dmis_REVEL0.5"])
    asd_all_raw[g] = nLGD + nDmis
Dict2Fil(asd_all_raw, f"{GW_DIR}/Spark_Meta_EWS.GeneWeight.csv")

# HIQ — BGMR corrected with IQ-specific LoF+Dmis counts and cohort size
hiq_gw = compute_asd_bgmr_weights(
    Spark_sub, count_col="HIQ_counts", BGMR=BGMR, N_proband=N_ASD_HIQ, n_top=61,
    lof_col="HIQ_LGD", dmis_col="HIQ_Dmis",
)
Dict2Fil(hiq_gw, f"{GW_DIR}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw")

# LIQ — BGMR corrected with IQ-specific LoF+Dmis counts and cohort size
liq_gw = compute_asd_bgmr_weights(
    Spark_sub, count_col="LIQ_counts", BGMR=BGMR, N_proband=N_ASD_LIQ, n_top=61,
    lof_col="LIQ_LGD", dmis_col="LIQ_Dmis",
)
Dict2Fil(liq_gw, f"{GW_DIR}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw")

# Print summaries
for label, gw in [("ASD All bgmr", asd_all_gw), ("ASD All raw", asd_all_raw),
                  ("HIQ bgmr", hiq_gw), ("LIQ bgmr", liq_gw)]:
    print(f"  {label}: {len(gw)} genes, first={Entrez2Symbol.get(list(gw.keys())[0], '?')}")

# %% [markdown]
# ## Raw HIQ/LIQ gene weights (PPV-weighted, no BGMR)

# %%
def CountMut(DF):
    N_LGD, N_mis, N_Dmis, N_syn = 0, 0, 0, 0
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


def Mut2GeneDF(MutDF, GeneSymbol2Entrez, pLI=True,
               LGD_weight_high=0.554, Dmis_weight_high=0.333,
               LGD_weight_low=0.138, Dmis_weight_low=0.130,
               LGD_weight_nopli=0.357, Dmis_weight_nopli=0.231):
    Select_Genes = np.array(list(set(MutDF["HGNC"].values)))
    gene2MutN = {}
    for g in Select_Genes:
        try:
            Entrez = int(GeneSymbol2Entrez[g])
        except Exception:
            continue
        Muts = MutDF[MutDF["HGNC"] == g]
        N_LGD, N_Mis, N_Dmis, N_Syn = CountMut(Muts)
        if pLI:
            try:
                pLI_val = float(Muts["ExACpLI"].values[0])
            except Exception:
                pLI_val = 0.0
            if pLI_val >= 0.5:
                gene2MutN[Entrez] = N_LGD * LGD_weight_high + N_Dmis * Dmis_weight_high
            else:
                gene2MutN[Entrez] = N_LGD * LGD_weight_low + N_Dmis * Dmis_weight_low
        else:
            gene2MutN[Entrez] = N_LGD * LGD_weight_nopli + N_Dmis * Dmis_weight_nopli
    return gene2MutN

# %%
# Filter to exome-wide significant genes for raw HIQ/LIQ
Spark_Denovo_ExomeWide = Spark_Denovo[Spark_Denovo["pDenovoWEST_Meta"] <= 0.01]
top_Genes = Spark_Denovo_ExomeWide["HGNC"].values
HighIQMuts = HighIQMuts_ALL[HighIQMuts_ALL["HGNC"].isin(top_Genes)]
LowIQMuts = LowIQMuts_ALL[LowIQMuts_ALL["HGNC"].isin(top_Genes)]

# Raw LGD+Dmis same weight (no pLI)
HIQ_GW = Mut2GeneDF(HighIQMuts, GeneSymbol2Entrez, pLI=False, LGD_weight_nopli=1, Dmis_weight_nopli=1)
LIQ_GW = Mut2GeneDF(LowIQMuts, GeneSymbol2Entrez, pLI=False, LGD_weight_nopli=1, Dmis_weight_nopli=1)
print(f"Raw HIQ: {len(HIQ_GW)} genes, Raw LIQ: {len(LIQ_GW)} genes")
Dict2Fil(HIQ_GW, f"{GW_DIR}/HIQ.top100.nopLI.LGD_Dmis_SameWeight.gw")
Dict2Fil(LIQ_GW, f"{GW_DIR}/LIQ.top100.nopLI.LGD_Dmis_SameWeight.gw")

# %% [markdown]
# # 22q11.2

# %%
Gene22q = [
    8214, 5625, 9993, 23617, 2928, 6576, 8218, 7290, 64976, 128977,
    8318, 7122, 2812, 6899, 54584, 10587, 1312, 421, 128989, 54487,
    27037, 5902, 29801, 388849, 65078, 85359, 373856, 7625, 91179, 84861,
    51586, 5297, 3053, 9342, 1399, 150209, 8216, 80764, 9127, 6545, 400891,
    8220, 7353, 5413, 79680, 728418,
]
X22q_GW = dict(zip(Gene22q, [1] * len(Gene22q)))

X22q_mouse_model_genes_corrected = [
    "DGCR6", "PRODH", "DGCR2", "TSSK2", "ESS2", "GSC2", "SLC25A1",
    "CLTCL1", "HIRA", "MRPL40", "UFD1", "CDC45", "CLDN5", "SEPTIN5",
    "GP1BB", "TBX1", "GNB1L", "TXNRD2", "COMT", "ARVCF", "CD38",
    "DGCR8", "TRMT2A", "RANBP1", "ZDHHC8", "RTN4R", "PRODH", "DGCR6L",
]
X22q_mouse_model_entrez = [GeneSymbol2Entrez[x] for x in X22q_mouse_model_genes_corrected]
X22q_V2_GW = dict(zip(X22q_mouse_model_entrez, [1] * len(X22q_mouse_model_entrez)))

Dict2Fil(X22q_GW, f"{GW_DIR}/X22q.gw.csv")
Dict2Fil(X22q_V2_GW, f"{GW_DIR}/X22q.mousemodel.gw.csv")

print(f"22q full: {len(X22q_GW)} genes, 22q mouse model: {len(X22q_V2_GW)} genes")

# %% [markdown]
# # NDD / DDD

# %%
# Load DDD data
df = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
df = df.sort_values("denovoWEST_p_full")
hc_df = df[df["denovoWEST_p_full"] <= 0.05 / 18762]
entrez_ids = [int(GeneSymbol2Entrez.get(x, -1)) for x in hc_df["symbol"].values]
hc_df["EntrezID"] = entrez_ids
# Filter out unmapped genes and genes not in BGMR
hc_df = hc_df[hc_df["EntrezID"] > 0]
hc_df = hc_df[hc_df["EntrezID"].isin(BGMR.index)]
print(f"DDD high-confidence genes (in BGMR): {hc_df.shape[0]}")

# %%
# DDD all high-confidence — BGMR corrected (Nproband=31058)
DDD_hc_GW = Aggregate_Gene_Weights_NDD(
    hc_df, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1,
    out=f"{GW_DIR}/DDD.hc.gw",
)

# DDD top 61 — BGMR corrected
NDD_top61_DF = hc_df.head(61)
DDD_top61_GW = Aggregate_Gene_Weights_NDD(
    NDD_top61_DF, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1,
)
Dict2Fil(DDD_top61_GW, f"{GW_DIR}/DDD.top61.gw.bgmr.csv")

# DDD top 285 — BGMR corrected
NDD_top285_DF = hc_df.head(285)
DDD_top285_GW = Aggregate_Gene_Weights_NDD(
    NDD_top285_DF, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1,
)
Dict2Fil(DDD_top285_GW, f"{GW_DIR}/DDD.top285.gw.bgmr.csv")

print(f"DDD hc: {len(DDD_hc_GW)} genes, top 61: {len(DDD_top61_GW)} genes, top 285: {len(DDD_top285_GW)} genes")

# %% [markdown]
# ## SCZ excluding NDD genes

# %%
NDD_top61_genes = hc_df.head(61)["EntrezID"].values
NDD_top297_genes = hc_df.head(297)["EntrezID"].values
NDD_top61_genes_set = set(int(x) for x in NDD_top61_genes)
NDD_top297_genes_set = set(int(x) for x in NDD_top297_genes)

SCZ_top61_ex_NDD61 = SCZ_top61 - NDD_top61_genes_set
SCZ_top61_ex_NDD297 = SCZ_top61 - NDD_top297_genes_set

print(f"SCZ top 61 overlap NDD top 61: {len(SCZ_top61.intersection(NDD_top61_genes_set))}")
print(f"SCZ top 61 overlap NDD top 297: {len(SCZ_top61.intersection(NDD_top297_genes_set))}")
print(f"SCZ top 61 ex NDD61: {len(SCZ_top61_ex_NDD61)} genes")
print(f"SCZ top 61 ex NDD297: {len(SCZ_top61_ex_NDD297)} genes")

# %%
# SCZ top 61 excluding NDD top 61
tmp_GeneDF = GeneDF.loc[list(SCZ_top61_ex_NDD61), :]
SCZ_exNDD61_GW = Aggregate_Gene_Weights_SCZ_Daly(
    tmp_GeneDF, valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_exNDD61_GW, f"{GW_DIR}/SCZ.top61.ExlNDD61.gw")

# SCZ top 61 excluding NDD top 297
tmp_GeneDF = GeneDF.loc[list(SCZ_top61_ex_NDD297), :]
SCZ_exNDD297_GW = Aggregate_Gene_Weights_SCZ_Daly(
    tmp_GeneDF, valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
Dict2Fil(SCZ_exNDD297_GW, f"{GW_DIR}/SCZ.top61.ExlNDD297.gw")

print(f"SCZ exNDD61: {len(SCZ_exNDD61_GW)} genes, SCZ exNDD297: {len(SCZ_exNDD297_GW)} genes")

# %% [markdown]
# # UKBB Cognitive Function
#
# Paper: https://www.nature.com/articles/s41588-023-01398-8
#
# - EDU: educational attainment
# - RT: reaction time
# - VNR: verbal-numerical reasoning

# %%
CogDF = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41588_2023_1398_MOESM3_ESM.xlsx"),
    sheet_name="Table S4",
)
CogDF = CogDF[CogDF["POPULATION"] == "EUR"]


def DF2GW_UKBB(DF, Ncarrier=False):
    res = {}
    for i, row in DF.iterrows():
        entrez = int(GeneSymbol2Entrez.get(row["GENE"], 0))
        if entrez == 0:
            continue
        if Ncarrier:
            res[entrez] = abs(row["BETA"]) * row["NCARRIER"]
        else:
            res[entrez] = abs(row["BETA"])
    return res

# %%
topN = 61

# VNR
VNR_DF = CogDF[CogDF["PHENOTYPE"] == "VNR"].sort_values("P")
VNR_Pos_DF = VNR_DF[VNR_DF["BETA"] > 0].head(topN)
VNR_Neg_DF = VNR_DF[VNR_DF["BETA"] < 0].head(topN)
VNR_Pos_GW = DF2GW_UKBB(VNR_Pos_DF)
VNR_Neg_GW = DF2GW_UKBB(VNR_Neg_DF)
Dict2Fil(VNR_Pos_GW, f"{GW_DIR}/UKBB_VNR_Pos_GW_{topN}.csv")
Dict2Fil(VNR_Neg_GW, f"{GW_DIR}/UKBB_VNR_Neg_GW_{topN}.csv")

# EDU
EDU_DF = CogDF[CogDF["PHENOTYPE"] == "EDU"].sort_values("P")
EDU_Pos_DF = EDU_DF[EDU_DF["BETA"] > 0].head(topN)
EDU_Neg_DF = EDU_DF[EDU_DF["BETA"] < 0].head(topN)
EDU_Pos_GW = DF2GW_UKBB(EDU_Pos_DF)
EDU_Neg_GW = DF2GW_UKBB(EDU_Neg_DF)
Dict2Fil(EDU_Pos_GW, f"{GW_DIR}/UKBB_EDU_Pos_GW_{topN}.csv")
Dict2Fil(EDU_Neg_GW, f"{GW_DIR}/UKBB_EDU_Neg_GW_{topN}.csv")

# RT
RT_DF = CogDF[CogDF["PHENOTYPE"] == "RT"].sort_values("P")
RT_Pos_DF = RT_DF[RT_DF["BETA"] > 0].head(topN)
RT_Neg_DF = RT_DF[RT_DF["BETA"] < 0].head(topN)
RT_Pos_GW = DF2GW_UKBB(RT_Pos_DF)
RT_Neg_GW = DF2GW_UKBB(RT_Neg_DF)
Dict2Fil(RT_Pos_GW, f"{GW_DIR}/UKBB_RT_Pos_GW_{topN}.csv")
Dict2Fil(RT_Neg_GW, f"{GW_DIR}/UKBB_RT_Neg_GW_{topN}.csv")

# VNR no-effect control
VNR_NoEff_DF = VNR_DF.tail(61)
VNR_NoEff_GW = DF2GW_UKBB(VNR_NoEff_DF)
Dict2Fil(VNR_NoEff_GW, f"{GW_DIR}/UKBB_VNR_NoEff_GW_{topN}.csv")

print(f"UKBB VNR pos/neg: {len(VNR_Pos_GW)}/{len(VNR_Neg_GW)}")
print(f"UKBB EDU pos/neg: {len(EDU_Pos_GW)}/{len(EDU_Neg_GW)}")
print(f"UKBB RT pos/neg: {len(RT_Pos_GW)}/{len(RT_Neg_GW)}")

# %% [markdown]
# # Negative Control Traits (genebass pLoF burden)
#
# Non-brain traits from UKBB exome sequencing (genebass).
# Top 61 genes by p-value, uniform weight=1.
# BETA direction is mixed for non-brain traits, so abs(BETA) weighting is not meaningful.
# Used to demonstrate CGE interneuron enrichment is specific to psychiatric/cognitive disorders.

# %%
CTRL_DIR = PROJ_DIR / "dat" / "CTRL"

neg_ctrl_traits = {
    # Non-brain metabolic/hematological
    "NegCtrl_HDL": CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30760-both_sexes--irnt_2025_11_25_15_49_42.csv",
    "NegCtrl_Alanine": CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30620-both_sexes--irnt_2025_12_16_20_18_07.csv",
    "NegCtrl_HbA1c": CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30750-both_sexes--irnt_2025_11_25_16_32_28.csv",
    "NegCtrl_T2D": CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-T2D_custom-both_sexes--custom_2025_11_25_16_14_20.csv",
    "NegCtrl_IBD": CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-IBD_custom2-both_sexes--custom_2025_11_25_16_13_49.csv",
    "NegCtrl_RBC": CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30070-both_sexes--irnt_2025_12_16_20_41_11.csv",
    # Neurological (positive controls — expected to show some neuronal signal)
    "NegCtrl_Parkinson": CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131022-both_sexes--_2025_11_25_16_56_06.csv",
    "NegCtrl_Alzheimer": CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131036-both_sexes--_2025_11_25_17_00_04.csv",
}

for trait_name, filepath in neg_ctrl_traits.items():
    df_trait = pd.read_csv(filepath, index_col=0)
    # Annotate with Entrez IDs
    df_trait["EntrezID"] = df_trait.index.map(
        lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez.get(x)) else None
    )
    # Top 61 genes by p-value, uniform weight=1
    # (BETA direction is mixed for non-brain traits — abs(BETA) weighting not meaningful)
    top = df_trait[df_trait["EntrezID"].notna()].head(61)
    gw = dict(zip(
        top["EntrezID"].astype(int).tolist(),
        [1] * len(top),
    ))
    Dict2Fil(gw, f"{GW_DIR}/{trait_name}.gw.csv")
    print(f"  {trait_name}: {len(gw)} genes")

# %% [markdown]
# # Verification

# %%
# Verify all pipeline gene weight files
for fname in [
    "Spark_Meta_EWS.GeneWeight.bgmr.csv",
    "Spark_Meta_EWS.GeneWeight.csv",
    "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "SCZ.top61.nopLI.LGD_Dmis_SameWeight.gw",
    "SCZ.top61.protect.gw",
    "SCZ.top100.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "SCZ.top200.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "SCZ.top61.ExlNDD61.gw",
    "SCZ.top61.ExlNDD297.gw",
    "DDD.top61.gw.bgmr.csv",
    "DDD.top285.gw.bgmr.csv",
    "DDD.hc.gw",
    "X22q.gw.csv",
    "X22q.mousemodel.gw.csv",
    "HIQ.top100.nopLI.LGD_Dmis_SameWeight.gw",
    "LIQ.top100.nopLI.LGD_Dmis_SameWeight.gw",
    "UKBB_VNR_Neg_GW_61.csv",
    "UKBB_VNR_Pos_GW_61.csv",
    "UKBB_EDU_Neg_GW_61.csv",
    "UKBB_EDU_Pos_GW_61.csv",
    "UKBB_RT_Neg_GW_61.csv",
    "UKBB_RT_Pos_GW_61.csv",
    "UKBB_VNR_NoEff_GW_61.csv",
    "NegCtrl_HDL.gw.csv",
    "NegCtrl_Alanine.gw.csv",
    "NegCtrl_HbA1c.gw.csv",
    "NegCtrl_T2D.gw.csv",
    "NegCtrl_IBD.gw.csv",
    "NegCtrl_RBC.gw.csv",
    "NegCtrl_Parkinson.gw.csv",
    "NegCtrl_Alzheimer.gw.csv",
]:
    fpath = GW_DIR / fname
    if fpath.exists():
        df_check = pd.read_csv(fpath, header=None)
        n_neg = (df_check[1] < 0).sum()
        first = Entrez2Symbol.get(int(df_check.iloc[0, 0]), "?")
        print(f"  {fname}: {len(df_check)} genes, first={first}, neg={n_neg}")
    else:
        print(f"  MISSING: {fname}")

# %%
