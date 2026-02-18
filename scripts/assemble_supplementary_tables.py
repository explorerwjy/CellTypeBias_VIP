#!/usr/bin/env python3
"""
Assemble supplementary tables for the CellTypeBias_VIP manuscript.

Outputs:
    results/supplementary/Table_S_bias_cluster.csv       - Cluster-level bias for all traits
    results/supplementary/Table_S_bias_supercluster.csv  - Supercluster-level bias for all traits
    results/supplementary/Table_S_gene_lists.csv         - Annotated gene lists with mutation counts

Usage:
    conda activate gencic
    python scripts/assemble_supplementary_tables.py
"""

import os
import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJ_DIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP"
RESULTS_DIR = os.path.join(
    PROJ_DIR,
    "results/conservation_model_WB_mean_phastCons_n_CDS_bases_Best1000/Centering",
)
# Allow override from the worktree when results are symlinked there
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WORKTREE_DIR = os.path.dirname(SCRIPT_DIR)
WORKTREE_RESULTS = os.path.join(
    WORKTREE_DIR,
    "results/conservation_model_WB_mean_phastCons_n_CDS_bases_Best1000/Centering",
)
if os.path.isdir(WORKTREE_RESULTS):
    RESULTS_DIR = WORKTREE_RESULTS

OUTPUT_DIR = os.path.join(WORKTREE_DIR, "results", "supplementary")

# Source data
ASD_ALL_EXCEL = os.path.join(PROJ_DIR, "dat/41588_2022_1148_MOESM4_ESM.xlsx")
ASD_IQ_CSV = os.path.join(PROJ_DIR, "dat/ASD_IQ_Mut.csv")
DDD_EXCEL = "/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx"
SCZ_CSV = os.path.join(PROJ_DIR, "dat/SCZ.ALLGENE.MutCountModified.csv")
BGMR_CSV = "/home/jw3514/Work/Resources/BGMR.withEntrez.csv"

# Gene weight directory
GW_DIR = os.path.join(PROJ_DIR, "dat/GeneWeights")

# Trait definitions
TRAITS = [
    "ASD_All",
    "ASD_HIQ",
    "ASD_LIQ",
    "SCZ",
    "SCZ_Neg",
    "DDD_61",
    "DDD_297",
    "22q_del",
    "NegCtrl_HDL",
    "NegCtrl_Alanine",
]

# Gene weight file mapping (must match config.yaml)
GW_FILES = {
    "ASD_All": os.path.join(GW_DIR, "Spark_Meta_EWS.GeneWeight.bgmr.csv"),
    "ASD_HIQ": os.path.join(GW_DIR, "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    "ASD_LIQ": os.path.join(GW_DIR, "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    "SCZ": os.path.join(GW_DIR, "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw"),
    "SCZ_Neg": os.path.join(GW_DIR, "SCZ.top61.protect.gw"),
    "DDD_61": os.path.join(GW_DIR, "DDD.top61.gw.bgmr.csv"),
    "DDD_297": os.path.join(GW_DIR, "DDD.hc.gw.csv"),
    "22q_del": os.path.join(GW_DIR, "X22q.gw.csv"),
    "NegCtrl_HDL": os.path.join(GW_DIR, "HDL_C.top61.gw.csv"),
    "NegCtrl_Alanine": os.path.join(GW_DIR, "Alanine.top61.gw.csv"),
}

# Cohort sizes for BGMR expected count calculation
COHORT_N = {
    "ASD_All": 42607,
    "ASD_HIQ": 4876,
    "ASD_LIQ": 2619,
    "DDD_61": 31058,
    "DDD_297": 31058,
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def load_gw(path):
    """Load a 2-column gene weight file (Entrez, Weight), no header."""
    df = pd.read_csv(path, header=None, names=["Entrez_ID", "Weight"])
    df["Entrez_ID"] = df["Entrez_ID"].astype(int)
    # Filter out invalid Entrez IDs (e.g., -1 from unmapped genes)
    df = df[df["Entrez_ID"] > 0].reset_index(drop=True)
    return df


def load_bgmr():
    """Load BGMR background mutation rate table, indexed by entrez_id."""
    bgmr = pd.read_csv(BGMR_CSV, low_memory=False)
    bgmr["entrez_id"] = bgmr["entrez_id"].astype(int)
    bgmr = bgmr.set_index("entrez_id")
    return bgmr


def build_entrez2symbol(bgmr):
    """Build a comprehensive Entrez -> Gene Symbol mapping from multiple sources."""
    e2s = {}
    # Primary: BGMR
    for eid, row in bgmr.iterrows():
        e2s[int(eid)] = row["GeneName"]
    # Fallback: SCZ mutation data
    scz = pd.read_csv(SCZ_CSV)
    for _, row in scz.iterrows():
        eid = int(row["Entrez"])
        if eid not in e2s:
            e2s[eid] = row["Gene Symbol"]
    # Fallback: ASD source data
    try:
        asd = pd.read_excel(ASD_ALL_EXCEL, sheet_name="Table S7", header=2)
        for _, row in asd.iterrows():
            eid = int(row["EntrezID"])
            if eid not in e2s:
                e2s[eid] = row["HGNC"]
    except Exception:
        pass
    # Fallback: ASD IQ mutation data
    try:
        iq = pd.read_csv(ASD_IQ_CSV)
        for _, row in iq.iterrows():
            eid = int(row["Entrez"])
            if eid > 0 and eid not in e2s:
                e2s[eid] = row["HGNC"]
    except Exception:
        pass
    # Hardcoded fallback for genes not in any source
    # (verified against NCBI Gene database)
    ncbi_fallback = {
        388849: "DGCR5",       # DiGeorge critical region 5 (22q11.21)
        400891: "LRRC74B",     # leucine rich repeat containing 74B (22q11.21)
        79680: "RTL10",        # retrotransposon Gag like 10 (22q11.21)
        728418: "POM121L7P",   # POM121 nucleoporin like 7 pseudogene (22q11.21)
        340533: "NEXMIF",      # neurite extension and migration factor
        10476: "ATP5PD",       # ATP synthase peripheral stalk subunit d
        55908: "ANGPTL8",      # angiopoietin like 8
        29919: "RMC1",         # regulator of MON1-CCZ1
        57082: "KNL1",         # kinetochore scaffold 1
        55248: "PACC1",        # proton activated chloride channel 1
        5413: "SEPTIN5",       # septin 5 (22q11.21)
        7353: "UFD1",          # ubiquitin recognition factor in ER degradation 1 (22q11.21)
        373856: "USP41P",      # ubiquitin specific peptidase 41 pseudogene (22q11.21)
    }
    for eid, sym in ncbi_fallback.items():
        if eid not in e2s:
            e2s[eid] = sym
    return e2s


# ---------------------------------------------------------------------------
# Table S_bias: cluster-level
# ---------------------------------------------------------------------------
def assemble_bias_cluster():
    """Concatenate all per-trait *_bias_addP.csv into a long-format table."""
    frames = []
    for trait in TRAITS:
        path = os.path.join(RESULTS_DIR, f"{trait}_bias_addP.csv")
        if not os.path.exists(path):
            print(f"WARNING: missing {path}, skipping {trait}")
            continue
        df = pd.read_csv(path, index_col=0)
        df.insert(0, "Trait", trait)
        df.index.name = "Cluster"
        frames.append(df)

    result = pd.concat(frames, axis=0)
    # Rename EFFECT -> Bias for clarity
    result = result.rename(columns={"EFFECT": "Bias"})
    # Select and reorder columns
    keep_cols = [
        "Trait",
        "Bias",
        "P-value",
        "q-value",
        "Z-score",
        "EFFECT_adj",
        "Class",
        "Supercluster",
        "Subtype",
        "Neurotransmitter",
        "Top three regions",
        "Top three dissections",
        "Number of cells",
        "Neuropeptide auto-annotation",
        "-logP",
        "Rank",
    ]
    result = result[[c for c in keep_cols if c in result.columns]]
    result = result.sort_values(["Trait", "P-value", "Bias"], ascending=[True, True, False])
    return result


# ---------------------------------------------------------------------------
# Table S_bias: supercluster-level
# ---------------------------------------------------------------------------
def assemble_bias_supercluster():
    """Concatenate all per-trait *_bias_addP_supercluster.csv."""
    frames = []
    for trait in TRAITS:
        path = os.path.join(RESULTS_DIR, f"{trait}_bias_addP_supercluster.csv")
        if not os.path.exists(path):
            print(f"WARNING: missing {path}, skipping {trait}")
            continue
        df = pd.read_csv(path, index_col=0)
        df.insert(0, "Trait", trait)
        df.index.name = "Supercluster"
        frames.append(df)

    result = pd.concat(frames, axis=0)
    result = result.rename(columns={"EFFECT": "Bias"})
    keep_cols = [
        "Trait",
        "Bias",
        "P-value",
        "q-value",
        "Z-score",
        "EFFECT_adj",
        "n_clusters",
        "-logP",
    ]
    result = result[[c for c in keep_cols if c in result.columns]]
    result = result.sort_values(["Trait", "P-value", "Bias"], ascending=[True, True, False])
    return result


# ---------------------------------------------------------------------------
# Table S_genes: annotated gene lists with mutation counts
# ---------------------------------------------------------------------------
def _build_asd_all_genes(bgmr, entrez2symbol):
    """ASD_All: de novo study with pLI-stratified weights, BGMR correction."""
    gw = load_gw(GW_FILES["ASD_All"])
    asd = pd.read_excel(ASD_ALL_EXCEL, sheet_name="Table S7", header=2)
    asd["EntrezID"] = asd["EntrezID"].astype(int)
    asd = asd.set_index("EntrezID")

    N = COHORT_N["ASD_All"]
    rows = []
    for _, gw_row in gw.iterrows():
        eid = int(gw_row["Entrez_ID"])
        weight = gw_row["Weight"]
        gene_symbol = entrez2symbol.get(eid, "")

        obs_lgd = asd.loc[eid, "AutismMerged_LoF"] if eid in asd.index else np.nan
        obs_dmis = (
            asd.loc[eid, "AutismMerged_Dmis_REVEL0.5"] if eid in asd.index else np.nan
        )
        pli = asd.loc[eid, "ExACpLI"] if eid in asd.index else np.nan

        # Expected counts from BGMR
        if eid in bgmr.index:
            exp_lgd = bgmr.loc[eid, "p_LGD"] * 2 * N
            exp_dmis = bgmr.loc[eid, "prevel_0.5"] * 2 * N
        else:
            exp_lgd = np.nan
            exp_dmis = np.nan

        # pLI-stratified weight multipliers
        try:
            pli_val = float(pli)
        except (ValueError, TypeError):
            pli_val = 0.0
        if pli_val >= 0.5:
            wLGD, wDmis = 0.554, 0.333
        else:
            wLGD, wDmis = 0.138, 0.130

        rows.append(
            {
                "Trait": "ASD_All",
                "Gene_Symbol": gene_symbol,
                "Entrez_ID": eid,
                "obs_LGD": obs_lgd,
                "obs_Dmis": obs_dmis,
                "exp_LGD": round(exp_lgd, 4) if not np.isnan(exp_lgd) else np.nan,
                "exp_Dmis": round(exp_dmis, 4) if not np.isnan(exp_dmis) else np.nan,
                "pLI": pli,
                "wLGD": wLGD,
                "wDmis": wDmis,
                "Weight": round(weight, 4),
            }
        )
    return pd.DataFrame(rows)


def _count_iq_mutations(iq_df, eid):
    """Count LGD and Dmis mutations for a gene in the ASD IQ data."""
    gene_muts = iq_df[iq_df["Entrez"] == eid]
    n_lgd = 0
    n_dmis = 0
    for _, row in gene_muts.iterrows():
        eff = str(row["GeneEff"]).split(";")[0]
        if eff in [
            "frameshift",
            "splice_acceptor",
            "splice_donor",
            "start_lost",
            "stop_gained",
            "stop_lost",
        ]:
            n_lgd += 1
        elif eff == "missense":
            revel = str(row["REVEL"]).split(";")[0]
            if revel != "." and revel != "nan":
                try:
                    if float(revel) > 0.5:
                        n_dmis += 1
                except ValueError:
                    pass
    return n_lgd, n_dmis


def _build_asd_iq_genes(trait, iq_threshold, bgmr, entrez2symbol):
    """ASD_HIQ or ASD_LIQ: de novo with equal LGD/Dmis weights, BGMR correction."""
    gw = load_gw(GW_FILES[trait])
    iq_df = pd.read_csv(ASD_IQ_CSV)

    if trait == "ASD_HIQ":
        iq_df = iq_df[iq_df["IQ"] >= iq_threshold]
        N = COHORT_N["ASD_HIQ"]
    else:
        iq_df = iq_df[iq_df["IQ"] < iq_threshold]
        N = COHORT_N["ASD_LIQ"]

    rows = []
    for _, gw_row in gw.iterrows():
        eid = int(gw_row["Entrez_ID"])
        weight = gw_row["Weight"]
        gene_symbol = entrez2symbol.get(eid, "")

        obs_lgd, obs_dmis = _count_iq_mutations(iq_df, eid)

        if eid in bgmr.index:
            exp_lgd = bgmr.loc[eid, "p_LGD"] * 2 * N
            exp_dmis = bgmr.loc[eid, "prevel_0.5"] * 2 * N
        else:
            exp_lgd = np.nan
            exp_dmis = np.nan

        pli_vals = iq_df[iq_df["Entrez"] == eid]["ExACpLI"]
        pli = float(pli_vals.iloc[0]) if len(pli_vals) > 0 else np.nan

        rows.append(
            {
                "Trait": trait,
                "Gene_Symbol": gene_symbol,
                "Entrez_ID": eid,
                "obs_LGD": obs_lgd,
                "obs_Dmis": obs_dmis,
                "exp_LGD": round(exp_lgd, 4) if not pd.isna(exp_lgd) else np.nan,
                "exp_Dmis": round(exp_dmis, 4) if not pd.isna(exp_dmis) else np.nan,
                "pLI": pli,
                "wLGD": 1.0,
                "wDmis": 1.0,
                "Weight": round(weight, 4),
            }
        )
    return pd.DataFrame(rows)


def _build_ddd_genes(trait, bgmr, entrez2symbol):
    """DDD_61 or DDD_297: de novo with equal weights, BGMR correction."""
    gw = load_gw(GW_FILES[trait])
    ddd = pd.read_excel(DDD_EXCEL, sheet_name="kaplanis_samocha_denovoWEST_res")

    # Build symbol -> entrez mapping from entrez2symbol (reverse)
    symbol2entrez = {}
    for eid, sym in entrez2symbol.items():
        symbol2entrez[sym] = eid

    # Build entrez -> DDD row mapping
    entrez2ddd = {}
    for _, row in ddd.iterrows():
        sym = row["symbol"]
        if sym in symbol2entrez:
            entrez2ddd[symbol2entrez[sym]] = row

    N = COHORT_N.get(trait, COHORT_N["DDD_61"])

    rows = []
    for _, gw_row in gw.iterrows():
        eid = int(gw_row["Entrez_ID"])
        weight = gw_row["Weight"]
        gene_symbol = entrez2symbol.get(eid, "")

        if eid in entrez2ddd:
            ddd_row = entrez2ddd[eid]
            obs_lgd = (
                ddd_row["frameshift_variant"]
                + ddd_row["splice_acceptor_variant"]
                + ddd_row["splice_donor_variant"]
                + ddd_row["stop_gained"]
                + ddd_row["stop_lost"]
            )
            obs_dmis = ddd_row["missense_variant"]
        else:
            obs_lgd = np.nan
            obs_dmis = np.nan

        if eid in bgmr.index:
            exp_lgd = bgmr.loc[eid, "p_LGD"] * 2 * N
            exp_dmis = bgmr.loc[eid, "p_misense"] * 2 * N
        else:
            exp_lgd = np.nan
            exp_dmis = np.nan

        rows.append(
            {
                "Trait": trait,
                "Gene_Symbol": gene_symbol,
                "Entrez_ID": eid,
                "obs_LGD": obs_lgd,
                "obs_Dmis": obs_dmis,
                "exp_LGD": round(exp_lgd, 4) if not pd.isna(exp_lgd) else np.nan,
                "exp_Dmis": round(exp_dmis, 4) if not pd.isna(exp_dmis) else np.nan,
                "pLI": np.nan,
                "wLGD": 1.0,
                "wDmis": 1.0,
                "Weight": round(weight, 4),
            }
        )
    return pd.DataFrame(rows)


def _build_scz_genes(trait, bgmr, entrez2symbol):
    """SCZ or SCZ_Neg: case-control excess."""
    gw = load_gw(GW_FILES[trait])
    scz = pd.read_csv(SCZ_CSV, index_col="Entrez")

    rows = []
    for _, gw_row in gw.iterrows():
        eid = int(gw_row["Entrez_ID"])
        weight = gw_row["Weight"]
        gene_symbol = entrez2symbol.get(eid, "")

        if eid in scz.index:
            scz_row = scz.loc[eid]
            # Handle potential duplicate index
            if isinstance(scz_row, pd.DataFrame):
                scz_row = scz_row.iloc[0]
            case_ptv = scz_row["Case PTV"]
            ctrl_ptv = scz_row["Ctrl PTV"]
            case_mis3 = scz_row["Case mis3"]
            ctrl_mis3 = scz_row["Ctrl mis3"]
            nlgd = scz_row["nLGD"]
            nmis3 = scz_row["nMis3"]
        else:
            case_ptv = np.nan
            ctrl_ptv = np.nan
            case_mis3 = np.nan
            ctrl_mis3 = np.nan
            nlgd = np.nan
            nmis3 = np.nan

        rows.append(
            {
                "Trait": trait,
                "Gene_Symbol": gene_symbol,
                "Entrez_ID": eid,
                "Case_PTV": case_ptv,
                "Ctrl_PTV": ctrl_ptv,
                "Case_mis3": case_mis3,
                "Ctrl_mis3": ctrl_mis3,
                "nLGD_excess": nlgd,
                "nMis3_excess": nmis3,
                "Weight": round(weight, 4),
            }
        )
    return pd.DataFrame(rows)


def _build_simple_genes(trait, entrez2symbol):
    """22q_del, NegCtrl_HDL, NegCtrl_Alanine: just gene symbol + entrez + weight."""
    gw = load_gw(GW_FILES[trait])
    rows = []
    for _, gw_row in gw.iterrows():
        eid = int(gw_row["Entrez_ID"])
        weight = gw_row["Weight"]
        gene_symbol = entrez2symbol.get(eid, "")
        rows.append(
            {
                "Trait": trait,
                "Gene_Symbol": gene_symbol,
                "Entrez_ID": eid,
                "Weight": round(weight, 4),
            }
        )
    return pd.DataFrame(rows)


def assemble_gene_lists():
    """Build gene list table with per-gene mutation details."""
    bgmr = load_bgmr()
    entrez2symbol = build_entrez2symbol(bgmr)

    print("  Building ASD_All gene list...")
    asd_all = _build_asd_all_genes(bgmr, entrez2symbol)

    print("  Building ASD_HIQ gene list...")
    asd_hiq = _build_asd_iq_genes("ASD_HIQ", 70, bgmr, entrez2symbol)

    print("  Building ASD_LIQ gene list...")
    asd_liq = _build_asd_iq_genes("ASD_LIQ", 70, bgmr, entrez2symbol)

    print("  Building DDD_61 gene list...")
    ddd_61 = _build_ddd_genes("DDD_61", bgmr, entrez2symbol)

    print("  Building DDD_297 gene list...")
    ddd_297 = _build_ddd_genes("DDD_297", bgmr, entrez2symbol)

    print("  Building SCZ gene list...")
    scz = _build_scz_genes("SCZ", bgmr, entrez2symbol)

    print("  Building SCZ_Neg gene list...")
    scz_neg = _build_scz_genes("SCZ_Neg", bgmr, entrez2symbol)

    print("  Building 22q_del gene list...")
    del_22q = _build_simple_genes("22q_del", entrez2symbol)

    print("  Building NegCtrl_HDL gene list...")
    negctrl_hdl = _build_simple_genes("NegCtrl_HDL", entrez2symbol)

    print("  Building NegCtrl_Alanine gene list...")
    negctrl_ala = _build_simple_genes("NegCtrl_Alanine", entrez2symbol)

    # Concatenate all — columns that don't apply will be NaN
    all_genes = pd.concat(
        [asd_all, asd_hiq, asd_liq, ddd_61, ddd_297, scz, scz_neg, del_22q, negctrl_hdl, negctrl_ala],
        axis=0,
        ignore_index=True,
    )

    # Reorder: common columns first, then de novo, then case-control
    col_order = [
        "Trait",
        "Gene_Symbol",
        "Entrez_ID",
        "obs_LGD",
        "obs_Dmis",
        "exp_LGD",
        "exp_Dmis",
        "pLI",
        "wLGD",
        "wDmis",
        "Case_PTV",
        "Ctrl_PTV",
        "Case_mis3",
        "Ctrl_mis3",
        "nLGD_excess",
        "nMis3_excess",
        "Weight",
    ]
    col_order = [c for c in col_order if c in all_genes.columns]
    all_genes = all_genes[col_order]
    all_genes = all_genes.sort_values(["Trait", "Weight"], ascending=[True, False])
    return all_genes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Assembling cluster-level bias table...")
    bias_cluster = assemble_bias_cluster()
    out_cluster = os.path.join(OUTPUT_DIR, "Table_S_bias_cluster.csv")
    bias_cluster.to_csv(out_cluster)
    print(f"  -> {out_cluster}  ({len(bias_cluster)} rows)")

    print("Assembling supercluster-level bias table...")
    bias_sc = assemble_bias_supercluster()
    out_sc = os.path.join(OUTPUT_DIR, "Table_S_bias_supercluster.csv")
    bias_sc.to_csv(out_sc)
    print(f"  -> {out_sc}  ({len(bias_sc)} rows)")

    print("Assembling gene list table...")
    gene_lists = assemble_gene_lists()
    out_genes = os.path.join(OUTPUT_DIR, "Table_S_gene_lists.csv")
    gene_lists.to_csv(out_genes, index=False)
    print(f"  -> {out_genes}  ({len(gene_lists)} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
