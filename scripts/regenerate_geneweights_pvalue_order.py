#!/usr/bin/env python3
"""
Regenerate gene weight files in p-value order (not sorted by weight).

Bug fix: Aggregate_Gene_Weights_NDD previously sorted output by weight descending.
Gene weight files should preserve the original p-value ranking from gene discovery.

This script regenerates ALL gene weight files that were affected by the sorting bug.
"""
import sys
import os
import csv
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# Load project config
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    cfg = yaml.safe_load(f)
PROJ_DIR = Path(cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import Aggregate_Gene_Weights_NDD, Aggregate_Gene_Weights_SCZ_Daly
from UNIMED import LoadGeneINFO

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

GW_DIR = PROJ_DIR / "dat" / "GeneWeights"

# Load expression matrix for valid gene filtering
HumanCT = pd.read_csv(
    PROJ_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    index_col=0
)
valid_genes = HumanCT.index.values

# Load BGMR
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv", low_memory=False)
BGMR["entrez_id"] = BGMR["entrez_id"].astype(int)
BGMR = BGMR.set_index("entrez_id")


def dict2fil_ordered(gene_dict, filepath):
    """Write gene weight dict to CSV, preserving insertion order."""
    with open(filepath, 'wt') as f:
        writer = csv.writer(f)
        for k, v in gene_dict.items():
            writer.writerow([k, v])
    print(f"  Wrote {len(gene_dict)} genes to {filepath}")


def verify_order(filepath, expected_first_gene, label=""):
    """Verify the first gene in the file matches expected."""
    df = pd.read_csv(filepath, header=None)
    actual_first = df.iloc[0, 0]
    symbol = Entrez2Symbol.get(actual_first, "?")
    status = "OK" if actual_first == expected_first_gene else "MISMATCH"
    print(f"  [{status}] {label}: first gene = {actual_first} ({symbol})")


# =====================================================================
# 1. DDD gene weights (sorted by denovoWEST_p_full)
# =====================================================================
print("\n" + "=" * 60)
print("1. DDD Gene Weights")
print("=" * 60)

DDD_Genes = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
DDD_Genes = DDD_Genes.sort_values("denovoWEST_p_full").reset_index(drop=True)
print(f"DDD genes loaded: {len(DDD_Genes)}, sorted by denovoWEST_p_full")
print(f"Top gene: {DDD_Genes.iloc[0]['symbol']} (p={DDD_Genes.iloc[0]['denovoWEST_p_full']:.2e})")

# Add EntrezID
DDD_Genes["EntrezID"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_Genes["symbol"]]

# Prepare mutation columns needed by Aggregate_Gene_Weights_NDD
# Fill missing mutation columns with 0
for col in ["frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "stop_lost", "missense_variant"]:
    if col not in DDD_Genes.columns:
        DDD_Genes[col] = 0
    DDD_Genes[col] = DDD_Genes[col].fillna(0)

# DDD top 61 (bgmr-corrected)
top61 = DDD_Genes.head(61)
top61_valid = top61[top61["EntrezID"].isin(BGMR.index)]
DDD_top61_gw = Aggregate_Gene_Weights_NDD(top61_valid, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1)
dict2fil_ordered(DDD_top61_gw, GW_DIR / "DDD.top61.gw.bgmr.csv")
verify_order(GW_DIR / "DDD.top61.gw.bgmr.csv", top61_valid.iloc[0]["EntrezID"], "DDD.top61.bgmr")

# DDD top 61 (no bgmr)
DDD_top61_gw_raw = Aggregate_Gene_Weights_NDD(top61_valid, Nproband=31058, BGMR=None, wLGD=1, wMis=1)
dict2fil_ordered(DDD_top61_gw_raw, GW_DIR / "DDD.top61.gw.csv")

# DDD top 61 (bgmr, written to .gw file too for compatibility)
Aggregate_Gene_Weights_NDD(top61_valid, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1,
                           out=str(GW_DIR / "DDD.top61.gw"))

# DDD top 285 (bgmr-corrected)
top285 = DDD_Genes.head(285)
top285_valid = top285[top285["EntrezID"].isin(BGMR.index)]
DDD_top285_gw = Aggregate_Gene_Weights_NDD(top285_valid, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1)
dict2fil_ordered(DDD_top285_gw, GW_DIR / "DDD.top285.gw.bgmr.csv")
verify_order(GW_DIR / "DDD.top285.gw.bgmr.csv", top285_valid.iloc[0]["EntrezID"], "DDD.top285.bgmr")

# DDD all high-confidence (hc)
hc = DDD_Genes[DDD_Genes["denovoWEST_p_full"] < 0.05]
hc_valid = hc[hc["EntrezID"].isin(BGMR.index)]
DDD_hc_gw = Aggregate_Gene_Weights_NDD(hc_valid, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1)
dict2fil_ordered(DDD_hc_gw, GW_DIR / "DDD.hc.gw.csv")
Aggregate_Gene_Weights_NDD(hc_valid, Nproband=31058, BGMR=BGMR, wLGD=1, wMis=1,
                           out=str(GW_DIR / "DDD.hc.gw"))

# =====================================================================
# 2. ASD gene weights (sorted by pDenovoWEST_Meta)
# =====================================================================
print("\n" + "=" * 60)
print("2. ASD Gene Weights")
print("=" * 60)

Spark_Denovo = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41588_2022_1148_MOESM4_ESM.xlsx"),
    skiprows=2, sheet_name="Table S7"
)
Spark_Denovo = Spark_Denovo[Spark_Denovo["pDenovoWEST_Meta"] != "."]
Spark_Denovo["pDenovoWEST_Meta"] = Spark_Denovo["pDenovoWEST_Meta"].astype(float)
Spark_Denovo = Spark_Denovo.sort_values("pDenovoWEST_Meta").reset_index(drop=True)
print(f"Spark genes loaded: {len(Spark_Denovo)}, sorted by pDenovoWEST_Meta")
print(f"Top gene: {Spark_Denovo.iloc[0]['HGNC']} (p={Spark_Denovo.iloc[0]['pDenovoWEST_Meta']:.2e})")

# Load IQ-stratified mutation data
Mut_n_IQ = pd.read_csv(str(PROJ_DIR / "dat/ASD_IQ_Mut.csv"))
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] > 70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] <= 70]

# Build per-gene mutation counts by IQ group, in p-value order
Spark_sub = Spark_Denovo[["EntrezID", "HGNC", "pDenovoWEST_Meta",
                          "AutismMerged_LoF", "AutismMerged_Dmis_REVEL0.5"]].copy()
Spark_sub = Spark_sub.set_index("EntrezID")

for i in Spark_sub.index:
    Spark_sub.loc[i, "HIQ_counts"] = HighIQMuts_ALL[HighIQMuts_ALL["Entrez"] == i].shape[0]
    Spark_sub.loc[i, "LIQ_counts"] = LowIQMuts_ALL[LowIQMuts_ALL["Entrez"] == i].shape[0]

# BGMR-correct: compute weights for top-61 genes that have HIQ/LIQ mutations
# HIQ: top genes by p-value that have HIQ mutations
hiq_genes_with_muts = Spark_sub[Spark_sub["HIQ_counts"] > 0]
liq_genes_with_muts = Spark_sub[Spark_sub["LIQ_counts"] > 0]

def compute_asd_bgmr_weights(gene_sub, count_col, n_top=61):
    """Compute BGMR-corrected ASD gene weights preserving p-value order."""
    # Take top N genes (by p-value, already sorted) that have mutations
    top = gene_sub.head(n_top)
    gw = {}
    for g in top.index:
        if g not in BGMR.index:
            continue
        nLGD = top.loc[g, "AutismMerged_LoF"] if "AutismMerged_LoF" in top.columns else 0
        nMis = top.loc[g, "AutismMerged_Dmis_REVEL0.5"] if "AutismMerged_Dmis_REVEL0.5" in top.columns else 0
        # BGMR correction
        exp_LGD = BGMR.loc[g, "p_LGD"] * 2 * 42607
        exp_Mis = BGMR.loc[g, "p_misense"] * 2 * 42607
        weight = (float(nLGD) - exp_LGD) * 1 + (float(nMis) - exp_Mis) * 1
        gw[g] = weight
    return gw

HIQ_top61_gw = compute_asd_bgmr_weights(hiq_genes_with_muts, "HIQ_counts", n_top=61)
LIQ_top61_gw = compute_asd_bgmr_weights(liq_genes_with_muts, "LIQ_counts", n_top=61)

dict2fil_ordered(HIQ_top61_gw, GW_DIR / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw")
dict2fil_ordered(LIQ_top61_gw, GW_DIR / "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw")
verify_order(GW_DIR / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
             list(HIQ_top61_gw.keys())[0], "HIQ.bgmr")
verify_order(GW_DIR / "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
             list(LIQ_top61_gw.keys())[0], "LIQ.bgmr")

# Spark Meta (ASD All) — top 61 by p-value, BGMR-corrected
all_asd_top61 = Spark_sub.head(61)
asd_all_gw = {}
for g in all_asd_top61.index:
    if g not in BGMR.index:
        continue
    nLGD = float(all_asd_top61.loc[g, "AutismMerged_LoF"])
    nMis = float(all_asd_top61.loc[g, "AutismMerged_Dmis_REVEL0.5"])
    exp_LGD = BGMR.loc[g, "p_LGD"] * 2 * 42607
    exp_Mis = BGMR.loc[g, "p_misense"] * 2 * 42607
    asd_all_gw[g] = (nLGD - exp_LGD) * 1 + (nMis - exp_Mis) * 1
dict2fil_ordered(asd_all_gw, GW_DIR / "Spark_Meta_EWS.GeneWeight.bgmr.csv")

# Non-bgmr version too
asd_all_gw_raw = {}
for g in all_asd_top61.index:
    nLGD = float(all_asd_top61.loc[g, "AutismMerged_LoF"])
    nMis = float(all_asd_top61.loc[g, "AutismMerged_Dmis_REVEL0.5"])
    asd_all_gw_raw[g] = nLGD * 1 + nMis * 1
dict2fil_ordered(asd_all_gw_raw, GW_DIR / "Spark_Meta_EWS.GeneWeight.csv")

print(f"ASD All (Spark Meta): {len(asd_all_gw)} genes")
print(f"ASD HIQ bgmr: {len(HIQ_top61_gw)} genes")
print(f"ASD LIQ bgmr: {len(LIQ_top61_gw)} genes")

# =====================================================================
# 3. SCZ ExNDD gene weights (bgmr-corrected, preserving p-value order)
# =====================================================================
print("\n" + "=" * 60)
print("3. SCZ ExNDD Gene Weights")
print("=" * 60)

GeneDF = pd.read_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"), index_col=0)
# GeneDF is already sorted by "P meta" (ascending)

# Load DDD gene lists for exclusion
DDD_top61_entrez = set(DDD_Genes.head(61)["EntrezID"].values)
DDD_top285_entrez = set(DDD_Genes.head(285)["EntrezID"].values)

# SCZ top 61 gene list
SCZ_top61 = GeneDF.head(61)
SCZ_top61_genes = set(SCZ_top61.index.values)

# SCZ minus NDD61
SCZ_ex_NDD61 = SCZ_top61_genes - DDD_top61_entrez
tmp_df = GeneDF.loc[[g for g in GeneDF.index if g in SCZ_ex_NDD61]].head(len(SCZ_ex_NDD61))
SCZ_ex_NDD61_gw = Aggregate_Gene_Weights_SCZ_Daly(
    tmp_df, valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
dict2fil_ordered(SCZ_ex_NDD61_gw, GW_DIR / "SCZ.top61.ExlNDD61.bgmr.gw")
print(f"SCZ ex NDD61: {len(SCZ_ex_NDD61_gw)} genes")

# SCZ minus NDD285
SCZ_ex_NDD285 = SCZ_top61_genes - DDD_top285_entrez
tmp_df = GeneDF.loc[[g for g in GeneDF.index if g in SCZ_ex_NDD285]].head(len(SCZ_ex_NDD285))
SCZ_ex_NDD285_gw = Aggregate_Gene_Weights_SCZ_Daly(
    tmp_df, valid_genes, mode="MC", lgd_weight=1, mis3_weight=1, mis2_weight=0
)
dict2fil_ordered(SCZ_ex_NDD285_gw, GW_DIR / "SCZ.top61.ExlNDD285.bgmr.gw")
print(f"SCZ ex NDD285: {len(SCZ_ex_NDD285_gw)} genes")

# =====================================================================
# Summary: Verify top gene is correct for all critical files
# =====================================================================
print("\n" + "=" * 60)
print("VERIFICATION: First gene in each file (should be top by p-value)")
print("=" * 60)

for fname, expected_symbol in [
    ("DDD.top61.gw.bgmr.csv", "CTCF"),
    ("DDD.top285.gw.bgmr.csv", "CTCF"),
    ("Spark_Meta_EWS.GeneWeight.bgmr.csv", "SCN2A"),
]:
    fpath = GW_DIR / fname
    if fpath.exists():
        df = pd.read_csv(fpath, header=None)
        first_gene = df.iloc[0, 0]
        symbol = Entrez2Symbol.get(first_gene, "?")
        status = "PASS" if symbol == expected_symbol else "FAIL"
        print(f"  [{status}] {fname}: {symbol} (expected {expected_symbol})")

print("\nDone.")
