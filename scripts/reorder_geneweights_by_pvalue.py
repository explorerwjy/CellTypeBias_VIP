#!/usr/bin/env python3
"""
Reorder gene weight files by p-value (not weight).

Reads existing gene weight files and reorders genes by their discovery p-value
from the original source data. Does NOT change which genes are included or their weights.
Only fixes the ordering.
"""
import sys
import csv
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    cfg = yaml.safe_load(f)
PROJ_DIR = Path(cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from UNIMED import LoadGeneINFO
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

GW_DIR = PROJ_DIR / "dat" / "GeneWeights"


def load_gw(filepath):
    """Load gene weight file as OrderedDict-like list of (entrez, weight)."""
    df = pd.read_csv(filepath, header=None, names=["EntrezID", "Weight"])
    return df


def write_gw(df, filepath):
    """Write gene weight DataFrame preserving order, no header."""
    df.to_csv(filepath, header=False, index=False)
    print(f"  Wrote {len(df)} genes to {filepath}")


def reorder_by_pval(gw_df, pval_rank_dict, label=""):
    """Reorder gene weight DataFrame by p-value rank. Genes not in rank dict go to end."""
    gw_df = gw_df.copy()
    gw_df["pval_rank"] = gw_df["EntrezID"].map(pval_rank_dict).fillna(99999)
    gw_df = gw_df.sort_values("pval_rank").drop(columns="pval_rank")
    first_gene = gw_df.iloc[0]["EntrezID"]
    symbol = Entrez2Symbol.get(int(first_gene), "?")
    print(f"  [{label}] First gene after reorder: {int(first_gene)} ({symbol})")
    return gw_df


# =====================================================================
# Build p-value rankings from source data
# =====================================================================

# DDD p-value ranking
DDD_Genes = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
DDD_Genes = DDD_Genes.sort_values("denovoWEST_p_full").reset_index(drop=True)
DDD_Genes["EntrezID"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_Genes["symbol"]]
ddd_pval_rank = {row["EntrezID"]: idx for idx, row in DDD_Genes.iterrows()}
print(f"DDD: Top gene = {DDD_Genes.iloc[0]['symbol']} (p={DDD_Genes.iloc[0]['denovoWEST_p_full']:.2e})")

# ASD (Spark) p-value ranking
Spark = pd.read_excel(
    str(PROJ_DIR / "dat/suppl.data/41588_2022_1148_MOESM4_ESM.xlsx"),
    skiprows=2, sheet_name="Table S7"
)
Spark = Spark[Spark["pDenovoWEST_Meta"] != "."]
Spark["pDenovoWEST_Meta"] = Spark["pDenovoWEST_Meta"].astype(float)
Spark = Spark.sort_values("pDenovoWEST_Meta").reset_index(drop=True)
asd_pval_rank = {int(row["EntrezID"]): idx for idx, row in Spark.iterrows()}
print(f"ASD: Top gene = {Spark.iloc[0]['HGNC']} (p={Spark.iloc[0]['pDenovoWEST_Meta']:.2e})")

# SCZ p-value ranking
GeneDF = pd.read_csv(str(PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"), index_col=0)
# GeneDF index order IS the p-value ranking (sorted by "P meta")
scz_pval_rank = {int(g): idx for idx, g in enumerate(GeneDF.index)}
top_scz_gene = GeneDF.index[0]
print(f"SCZ: Top gene = {int(top_scz_gene)} ({Entrez2Symbol.get(int(top_scz_gene), '?')})")

# =====================================================================
# Reorder DDD files
# =====================================================================
print("\n" + "=" * 60)
print("DDD files")
print("=" * 60)

for fname in [
    "DDD.top61.gw.bgmr.csv",
    "DDD.top61.gw.csv",
    "DDD.top61.gw",
    "DDD.top285.gw.bgmr.csv",
    "DDD.hc.gw.csv",
    "DDD.hc.gw",
]:
    fpath = GW_DIR / fname
    if fpath.exists():
        gw = load_gw(fpath)
        gw_new = reorder_by_pval(gw, ddd_pval_rank, label=fname)
        write_gw(gw_new, fpath)
    else:
        print(f"  SKIP (not found): {fname}")

# =====================================================================
# Reorder ASD files
# =====================================================================
print("\n" + "=" * 60)
print("ASD files")
print("=" * 60)

for fname in [
    "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "Spark_Meta_EWS.GeneWeight.bgmr.csv",
    "Spark_Meta_EWS.GeneWeight.csv",
]:
    fpath = GW_DIR / fname
    if fpath.exists():
        gw = load_gw(fpath)
        gw_new = reorder_by_pval(gw, asd_pval_rank, label=fname)
        write_gw(gw_new, fpath)
    else:
        print(f"  SKIP (not found): {fname}")

# =====================================================================
# Reorder SCZ ExNDD files
# =====================================================================
print("\n" + "=" * 60)
print("SCZ ExNDD files")
print("=" * 60)

for fname in [
    "SCZ.top61.ExlNDD61.bgmr.gw",
    "SCZ.top61.ExlNDD285.bgmr.gw",
]:
    fpath = GW_DIR / fname
    if fpath.exists():
        gw = load_gw(fpath)
        gw_new = reorder_by_pval(gw, scz_pval_rank, label=fname)
        write_gw(gw_new, fpath)
    else:
        print(f"  SKIP (not found): {fname}")

# =====================================================================
# Verification: compare old vs new first gene
# =====================================================================
print("\n" + "=" * 60)
print("FINAL VERIFICATION")
print("=" * 60)

checks = [
    ("DDD.top61.gw.bgmr.csv", "CTCF", ddd_pval_rank),
    ("DDD.top285.gw.bgmr.csv", "CTCF", ddd_pval_rank),
    ("Spark_Meta_EWS.GeneWeight.bgmr.csv", None, asd_pval_rank),
    ("HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw", None, asd_pval_rank),
    ("LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw", None, asd_pval_rank),
]

for fname, expected, rank_dict in checks:
    fpath = GW_DIR / fname
    if fpath.exists():
        gw = load_gw(fpath)
        first = int(gw.iloc[0]["EntrezID"])
        sym = Entrez2Symbol.get(first, "?")
        n_genes = len(gw)
        # Check monotonically increasing rank
        ranks = [rank_dict.get(int(g), 99999) for g in gw["EntrezID"]]
        is_sorted = all(ranks[i] <= ranks[i+1] for i in range(len(ranks)-1))
        status = "PASS" if is_sorted else "FAIL"
        if expected and sym != expected:
            status = "FAIL"
        print(f"  [{status}] {fname}: {n_genes} genes, first={sym}, p-val sorted={is_sorted}")

print("\nDone.")
