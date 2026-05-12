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
# # Assemble Supplementary Table for Supplementary Table: Mutation Bias Results
#
# Reviewer 3 minor 1: *"Please provide full mutation bias estimates for all
# traits tested, not just selected highlights. Preferably provide all mutational
# bias results for all traits tested in a 'long' format, with additional column
# describing the tested trait. Gene lists and mutation counts used for each
# trait should be published in supplementary tables."*
#
# This notebook reads `dat/SuppTable_MutationBias_manifest.yaml` and assembles
# a single multi-sheet Excel workbook containing:
#
# 1. **Trait_Index** — one row per trait with display name, category, source,
#    figure references, and gene/cell type counts
# 2. **GeneWeights_DeNovo** — de novo mutation gene weights with raw counts,
#    expected rates, and BGMR-corrected weights (ASD, DDD traits)
# 3. **GeneWeights_CaseControl** — case-control gene weights with case/ctrl
#    counts and effective excess weights (SCZ traits)
# 4. **GeneWeights_Membership** — membership gene weights (UKBB, CNV, controls)
# 5. **Bias_Long** — long-format bias results from BOTH null models:
#    (trait, cell_type, EFFECT, P_random, q_random, P_matched, q_matched)
# 6. **Bias_Cluster_Wide** — pivot of EFFECT (cluster x trait)
# 7. **Pvalue_Wide** — pivot of P_random (cluster x trait)
# 8. **Pvalue_Matched_Wide** — pivot of P_matched (cluster x trait)
# 9. **Qvalue_Wide** — pivot of q_random (cluster x trait)
# 10. **Qvalue_Matched_Wide** — pivot of q_matched (cluster x trait)
# 11. **Bias_Supercluster** — long-format supercluster-level bias with both nulls
#
# Output: `results/SuppTable_MutationBias.xlsx`

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import warnings; warnings.filterwarnings('ignore')
import yaml
from pathlib import Path
from collections import OrderedDict
import numpy as np
import pandas as pd

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import LoadGeneINFO

_, _, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

MANIFEST_FILE = PROJ_DIR / "dat/SuppTable_MutationBias_manifest.yaml"
OUTPUT_FILE = PROJ_DIR / "results/SuppTable_MutationBias.xlsx"

# ---- Source data for de novo and case-control gene weight sheets ----
# BGMR — background mutation rate, indexed by Entrez ID
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv")
BGMR = BGMR.set_index("entrez_id")
BGMR.index = BGMR.index.astype(int)

# ASD SPARK denovo meta-analysis (Table S7)
ASD_SPARK_FILE = PROJ_DIR / "dat/suppl.data/41588_2022_1148_MOESM4_ESM.xlsx"
# DDD Kaplanis et al. 2020
DDD_FILE = Path("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
# ASD IQ mutation data (for HIQ/LIQ split)
ASD_IQ_FILE = PROJ_DIR / "dat/ASD_IQ_Mut.csv"
# SCZ case-control mutation counts (modified)
SCZ_MUT_FILE = PROJ_DIR / "dat/SCZ.ALLGENE.MutCountModified.csv"

# Cohort sizes
N_ASD_ALL = 42607
N_ASD_HIQ = 3610
N_ASD_LIQ = 1792
N_DDD = 31058
N_SCZ_CASE = 24248
N_SCZ_CTRL = 97322

# %% [markdown]
# ## Load manifest

# %%
with open(MANIFEST_FILE) as f:
    manifest = yaml.safe_load(f)

print(f"Loaded manifest with {len(manifest)} traits")
for trait, info in manifest.items():
    print(f"  {trait:20s}  {info['category']:20s}  {info['display_name']}")

# %% [markdown]
# ## Load gene weights and bias for each trait

# %%
def load_geneweight_file(path):
    """Gene weight files are 2-column CSVs with no header: Entrez_ID, Weight."""
    df = pd.read_csv(path, header=None, names=["Entrez", "Weight"])
    df["Entrez"] = df["Entrez"].astype(int)
    df["Symbol"] = df["Entrez"].map(lambda e: Entrez2Symbol.get(e, ""))
    return df


# Validate all required files exist before processing
missing = []
for trait, info in manifest.items():
    for key in ("weight_file", "bias_file", "bias_supercluster",
                "bias_file_matched", "bias_supercluster_matched"):
        p = PROJ_DIR / info[key]
        if not p.exists():
            missing.append((trait, key, str(p)))
if missing:
    print("Missing files:")
    for m in missing:
        print(f"  {m}")
    raise FileNotFoundError(f"{len(missing)} required files missing")
print("All manifest files exist.")

# %%
# Load gene weights (basic summary for Trait_Index)
gw_summary = {}
for trait, info in manifest.items():
    gw = load_geneweight_file(PROJ_DIR / info["weight_file"])
    gw_summary[trait] = {
        "n_genes": len(gw),
        "min_weight": gw["Weight"].min(),
        "max_weight": gw["Weight"].max(),
        "mean_weight": gw["Weight"].mean(),
    }

# %% [markdown]
# ### GeneWeights_DeNovo (ASD, DDD traits)

# %%
# --- Load source data for de novo gene weights ---
# ASD SPARK Table S7
Spark_Denovo = pd.read_excel(
    str(ASD_SPARK_FILE), skiprows=2, sheet_name="Table S7",
)
Spark_Denovo = Spark_Denovo[Spark_Denovo["pDenovoWEST_Meta"] != "."]
Spark_Denovo["pDenovoWEST_Meta"] = Spark_Denovo["pDenovoWEST_Meta"].astype(float)
Spark_Denovo = Spark_Denovo.sort_values("pDenovoWEST_Meta").reset_index(drop=True)
Spark_Denovo = Spark_Denovo.set_index("EntrezID")

# ASD IQ mutation data for HIQ/LIQ split
Mut_n_IQ = pd.read_csv(str(ASD_IQ_FILE))
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] > 70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] <= 70]

# Count IQ-specific LGD and Dmis per gene (same logic as Bias_Mutation_Weights.py)
lgd_types = ["frameshift", "splice_acceptor", "splice_donor", "start_lost", "stop_gained", "stop_lost"]
for entrez_id in Spark_Denovo.index:
    for prefix, muts in [("HIQ", HighIQMuts_ALL), ("LIQ", LowIQMuts_ALL)]:
        gene_muts = muts[muts["Entrez"] == entrez_id]
        if len(gene_muts) == 0:
            Spark_Denovo.loc[entrez_id, f"{prefix}_LGD"] = 0
            Spark_Denovo.loc[entrez_id, f"{prefix}_Dmis"] = 0
            continue
        eff = gene_muts["GeneEff"].str.split(";").str[0]
        n_lgd = (eff.isin(lgd_types)).sum()
        is_mis = eff == "missense"
        revel = gene_muts.loc[is_mis, "REVEL"].str.split(";").str[0]
        n_dmis = (revel.apply(lambda x: float(x) > 0.5 if x != "." else False)).sum()
        Spark_Denovo.loc[entrez_id, f"{prefix}_LGD"] = n_lgd
        Spark_Denovo.loc[entrez_id, f"{prefix}_Dmis"] = n_dmis

# DDD Kaplanis data
DDD_df = pd.read_excel(str(DDD_FILE))
DDD_df = DDD_df.sort_values("denovoWEST_p_full")
DDD_df["EntrezID"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_df["symbol"].values]
DDD_df = DDD_df[DDD_df["EntrezID"] > 0]
DDD_df = DDD_df[DDD_df["EntrezID"].isin(BGMR.index)]
DDD_df = DDD_df.set_index("EntrezID")
print(f"Spark Denovo: {len(Spark_Denovo)} genes, DDD: {len(DDD_df)} genes")

# %%
# Build GeneWeights_DeNovo rows
denovo_rows = []

# ASD trait configs: (trait_key, n_top, lof_col, dmis_col, N_proband, mis_filter, count_filter_col)
asd_configs = {
    "ASD_All": (61, "AutismMerged_LoF", "AutismMerged_Dmis_REVEL0.5", N_ASD_ALL, "Dmis REVEL>0.5", None),
    "ASD_HIQ": (61, "HIQ_LGD", "HIQ_Dmis", N_ASD_HIQ, "Dmis REVEL>0.5", "HIQ_LGD"),
    "ASD_LIQ": (61, "LIQ_LGD", "LIQ_Dmis", N_ASD_LIQ, "Dmis REVEL>0.5", "LIQ_LGD"),
}

for trait_key, (n_top, lof_col, dmis_col, n_proband, mis_filter, count_col) in asd_configs.items():
    top = Spark_Denovo.head(n_top).copy()
    # For HIQ/LIQ, filter to genes with any mutation in that IQ group
    if count_col is not None:
        # Total mutations = LGD + Dmis for this IQ group
        lgd_c = count_col
        dmis_c = count_col.replace("_LGD", "_Dmis")
        top = top[(top[lgd_c].fillna(0) + top[dmis_c].fillna(0)) > 0]
    for g in top.index:
        if g not in BGMR.index:
            continue
        nLGD = float(top.loc[g, lof_col])
        nMis = float(top.loc[g, dmis_col])
        exp_LGD = BGMR.loc[g, "p_LGD"] * 2 * n_proband
        exp_Mis = BGMR.loc[g, "prevel_0.5"] * 2 * n_proband
        weight = max((nLGD - exp_LGD) + (nMis - exp_Mis), 0)
        denovo_rows.append({
            "Trait": trait_key,
            "Gene_Entrez": int(g),
            "Gene_Symbol": Entrez2Symbol.get(int(g), ""),
            "nLGD": nLGD,
            "nMis": nMis,
            "Mis_Filter": mis_filter,
            "N_proband": n_proband,
            "exp_LGD": round(exp_LGD, 6),
            "exp_Mis": round(exp_Mis, 6),
            "Weight": round(weight, 6),
        })

# DDD configs
ddd_configs = {
    "DDD_61": 61,
    "DDD_285": 285,
}

for trait_key, n_top in ddd_configs.items():
    top = DDD_df.head(n_top)
    for g in top.index:
        if g not in BGMR.index:
            continue
        nLGD = (top.loc[g, "frameshift_variant"] + top.loc[g, "splice_acceptor_variant"]
                + top.loc[g, "splice_donor_variant"] + top.loc[g, "stop_gained"]
                + top.loc[g, "stop_lost"])
        nMis = top.loc[g, "missense_variant"]
        exp_LGD = BGMR.loc[g, "p_LGD"] * 2 * N_DDD
        exp_Mis = BGMR.loc[g, "p_misense"] * 2 * N_DDD
        weight = (nLGD - exp_LGD) * 1 + (nMis - exp_Mis) * 1
        denovo_rows.append({
            "Trait": trait_key,
            "Gene_Entrez": int(g),
            "Gene_Symbol": Entrez2Symbol.get(int(g), ""),
            "nLGD": nLGD,
            "nMis": nMis,
            "Mis_Filter": "All missense",
            "N_proband": N_DDD,
            "exp_LGD": round(exp_LGD, 6),
            "exp_Mis": round(exp_Mis, 6),
            "Weight": round(weight, 6),
        })

GeneWeights_DeNovo = pd.DataFrame(denovo_rows)
print(f"GeneWeights_DeNovo: {len(GeneWeights_DeNovo)} rows")
for t in GeneWeights_DeNovo["Trait"].unique():
    n = (GeneWeights_DeNovo["Trait"] == t).sum()
    print(f"  {t}: {n} genes")

# %% [markdown]
# ### GeneWeights_CaseControl (SCZ traits)

# %%
# Load SCZ case-control data
SCZ_GeneDF = pd.read_csv(str(SCZ_MUT_FILE), index_col=0)

cc_traits = ["SCZ", "SCZ_protect", "SCZ_rmNDD61", "SCZ_rmNDD285"]
cc_rows = []

for trait_key in cc_traits:
    info = manifest[trait_key]
    # Load which genes are in this trait's gene weight file
    gw = load_geneweight_file(PROJ_DIR / info["weight_file"])
    trait_entrez = set(gw["Entrez"].values)

    for g in trait_entrez:
        if g not in SCZ_GeneDF.index:
            continue
        row = SCZ_GeneDF.loc[g]
        case_ptv = float(row["Case PTV"]) if not pd.isna(row.get("Case PTV")) else 0
        ctrl_ptv = float(row["Ctrl PTV"]) if not pd.isna(row.get("Ctrl PTV")) else 0
        case_mis3 = float(row["Case mis3"]) if not pd.isna(row.get("Case mis3")) else 0
        ctrl_mis3 = float(row["Ctrl mis3"]) if not pd.isna(row.get("Ctrl mis3")) else 0
        # Weight from the gene weight file (already computed)
        weight = float(gw.loc[gw["Entrez"] == g, "Weight"].values[0])
        cc_rows.append({
            "Trait": trait_key,
            "Gene_Entrez": int(g),
            "Gene_Symbol": Entrez2Symbol.get(int(g), ""),
            "Case_PTV": case_ptv,
            "Ctrl_PTV": ctrl_ptv,
            "Case_Mis3": case_mis3,
            "Ctrl_Mis3": ctrl_mis3,
            "N_case": N_SCZ_CASE,
            "N_ctrl": N_SCZ_CTRL,
            "Weight": round(weight, 6),
        })

GeneWeights_CaseControl = pd.DataFrame(cc_rows)
print(f"GeneWeights_CaseControl: {len(GeneWeights_CaseControl)} rows")
for t in GeneWeights_CaseControl["Trait"].unique():
    n = (GeneWeights_CaseControl["Trait"] == t).sum()
    print(f"  {t}: {n} genes")

# %% [markdown]
# ### GeneWeights_Membership (all other traits)

# %%
membership_rows = []
for trait, info in manifest.items():
    if info.get("gene_weight_type") != "membership":
        continue
    gw = load_geneweight_file(PROJ_DIR / info["weight_file"])
    for _, r in gw.iterrows():
        membership_rows.append({
            "Trait": trait,
            "Gene_Entrez": int(r["Entrez"]),
            "Gene_Symbol": r["Symbol"],
            "Weight": round(float(r["Weight"]), 6),
        })

GeneWeights_Membership = pd.DataFrame(membership_rows)
print(f"GeneWeights_Membership: {len(GeneWeights_Membership)} rows")
for t in GeneWeights_Membership["Trait"].unique():
    n = (GeneWeights_Membership["Trait"] == t).sum()
    print(f"  {t}: {n} genes")

# %%
# Load cluster-level bias from BOTH null models and merge
bias_long_rows = []
bias_summary = {}
for trait, info in manifest.items():
    # Random null
    bias_r = pd.read_csv(PROJ_DIR / info["bias_file"], index_col=0)
    bias_r = bias_r.reset_index().rename(columns={"index": "CellType"})

    # Matched null
    bias_m = pd.read_csv(PROJ_DIR / info["bias_file_matched"], index_col=0)
    bias_m = bias_m.reset_index().rename(columns={"index": "CellType"})

    # Merge: keep EFFECT and annotation from random, add P/q from both
    # Annotation columns (same in both files)
    anno_cols_in = ["CellType", "Class", "Supercluster", "Subtype",
                    "Neurotransmitter", "Top three regions", "Top three dissections",
                    "Number of cells", "Neuropeptide auto-annotation"]
    anno_cols_present = [c for c in anno_cols_in if c in bias_r.columns]

    merged = bias_r[anno_cols_present + ["EFFECT"]].copy()
    merged["P_random"] = bias_r["P-value"].values
    merged["q_random"] = bias_r["q-value"].values
    merged["P_matched"] = bias_m["P-value"].values
    merged["q_matched"] = bias_m["q-value"].values
    merged["Trait"] = trait
    merged["Trait_DisplayName"] = info["display_name"]

    bias_long_rows.append(merged)
    bias_summary[trait] = {
        "n_cell_types": len(merged),
        "n_sig_q01": (merged["q_random"] < 0.1).sum(),
        "n_sig_q005": (merged["q_random"] < 0.05).sum(),
    }

Bias_Long = pd.concat(bias_long_rows, ignore_index=True)
# Reorder columns: Trait info first, then bias columns
front_cols = ["Trait", "Trait_DisplayName", "CellType", "Class", "Supercluster", "Subtype"]
back_cols = [c for c in Bias_Long.columns if c not in front_cols]
Bias_Long = Bias_Long[front_cols + back_cols]
print(f"Bias_Long: {len(Bias_Long)} rows ({Bias_Long['Trait'].nunique()} traits × {Bias_Long['CellType'].nunique()} clusters)")

# %%
# Load supercluster-level bias from BOTH null models and merge
sc_long_rows = []
for trait, info in manifest.items():
    sc_r = pd.read_csv(PROJ_DIR / info["bias_supercluster"], index_col=0)
    sc_r = sc_r.reset_index().rename(columns={"index": "Supercluster"})

    sc_m = pd.read_csv(PROJ_DIR / info["bias_supercluster_matched"], index_col=0)
    sc_m = sc_m.reset_index().rename(columns={"index": "Supercluster"})

    merged = sc_r[["Supercluster", "EFFECT", "n_clusters"]].copy()
    merged["P_random"] = sc_r["P-value"].values
    merged["q_random"] = sc_r["q-value"].values
    merged["P_matched"] = sc_m["P-value"].values
    merged["q_matched"] = sc_m["q-value"].values
    merged["Trait"] = trait
    merged["Trait_DisplayName"] = info["display_name"]
    sc_long_rows.append(merged)

Bias_Supercluster = pd.concat(sc_long_rows, ignore_index=True)
front_cols_sc = ["Trait", "Trait_DisplayName", "Supercluster"]
back_cols_sc = [c for c in Bias_Supercluster.columns if c not in front_cols_sc]
Bias_Supercluster = Bias_Supercluster[front_cols_sc + back_cols_sc]
print(f"Bias_Supercluster: {len(Bias_Supercluster)} rows")

# %% [markdown]
# ## Build Trait Index sheet

# %%
trait_index_rows = []
for trait, info in manifest.items():
    row = {
        "Trait": trait,
        "Display_Name": info["display_name"],
        "Category": info["category"],
        "Gene_Weight_Type": info.get("gene_weight_type", "membership"),
        "N_Genes": gw_summary[trait]["n_genes"],
        "Weight_Format": info["weight_format"],
        "Mean_Weight": round(gw_summary[trait]["mean_weight"], 4),
        "Min_Weight": round(gw_summary[trait]["min_weight"], 4),
        "Max_Weight": round(gw_summary[trait]["max_weight"], 4),
        "N_CellTypes": bias_summary[trait]["n_cell_types"],
        "N_Sig_q05_random": bias_summary[trait]["n_sig_q005"],
        "N_Sig_q10_random": bias_summary[trait]["n_sig_q01"],
        "Source": info["source"],
        "Figure_Reference": info["figure_ref"],
    }
    trait_index_rows.append(row)

Trait_Index = pd.DataFrame(trait_index_rows)
print(Trait_Index[["Trait", "Category", "N_Genes", "N_Sig_q10_random"]].to_string(index=False))

# %% [markdown]
# ## Build Bias_Cluster_Wide pivot (cluster × trait, EFFECT values)

# %%
# Pivot EFFECT values: rows = cell type cluster, columns = trait
Bias_Cluster_Wide = Bias_Long.pivot(index="CellType", columns="Trait", values="EFFECT")

# Add cluster annotation columns from one of the bias files (they all have the same annotation)
first_trait = list(manifest.keys())[0]
ref_bias = pd.read_csv(PROJ_DIR / manifest[first_trait]["bias_file"], index_col=0)
anno_cols = ["Class", "Supercluster", "Subtype", "Neurotransmitter",
             "Top three regions", "Top three dissections", "Number of cells"]
anno_present = [c for c in anno_cols if c in ref_bias.columns]
anno_df = ref_bias[anno_present]

Bias_Cluster_Wide = anno_df.join(Bias_Cluster_Wide, how="right")
Bias_Cluster_Wide.index.name = "CellType"
Bias_Cluster_Wide = Bias_Cluster_Wide.reset_index()

# Order trait columns by manifest order
trait_order = list(manifest.keys())
final_cols = ["CellType"] + anno_present + trait_order
Bias_Cluster_Wide = Bias_Cluster_Wide[final_cols]
print(f"Bias_Cluster_Wide: {Bias_Cluster_Wide.shape}")

# %% [markdown]
# ## Build P-value and q-value Wide pivots (random + matched)

# %%
def build_wide_pivot(bias_long, value_col, anno_df, anno_present, trait_order):
    """Build a wide pivot table from Bias_Long for the given value column."""
    wide = bias_long.pivot(index="CellType", columns="Trait", values=value_col)
    wide = anno_df.join(wide, how="right")
    wide.index.name = "CellType"
    wide = wide.reset_index()
    wide = wide[["CellType"] + anno_present + trait_order]
    return wide

Pvalue_Wide = build_wide_pivot(Bias_Long, "P_random", anno_df, anno_present, trait_order)
Pvalue_Matched_Wide = build_wide_pivot(Bias_Long, "P_matched", anno_df, anno_present, trait_order)
Qvalue_Wide = build_wide_pivot(Bias_Long, "q_random", anno_df, anno_present, trait_order)
Qvalue_Matched_Wide = build_wide_pivot(Bias_Long, "q_matched", anno_df, anno_present, trait_order)

print(f"Pvalue_Wide: {Pvalue_Wide.shape}")
print(f"Pvalue_Matched_Wide: {Pvalue_Matched_Wide.shape}")
print(f"Qvalue_Wide: {Qvalue_Wide.shape}")
print(f"Qvalue_Matched_Wide: {Qvalue_Matched_Wide.shape}")

# %% [markdown]
# ## README sheet

# %%
readme_lines = [
    ("Supplementary Table: Mutation Bias Results — Supplementary Table", ""),
    ("", ""),
    ("This workbook provides full mutation bias estimates and gene weights for", ""),
    ("all traits tested in the manuscript, in long format with explicit trait", ""),
    ("identifiers for downstream analysis.", ""),
    ("", ""),
    ("SHEETS", ""),
    ("------", ""),
    ("Trait_Index", "One row per trait. Includes display name, category, gene weight type, source, figure reference, and gene/cell type counts."),
    ("GeneWeights_DeNovo", "De novo mutation gene weights (ASD, DDD). Columns: Trait, Gene_Entrez, Gene_Symbol, nLGD, nMis, Mis_Filter, N_proband, exp_LGD, exp_Mis, Weight. Weight = (nLGD - exp_LGD) + (nMis - exp_Mis) after BGMR correction."),
    ("GeneWeights_CaseControl", "Case-control gene weights (SCZ). Columns: Trait, Gene_Entrez, Gene_Symbol, Case_PTV, Ctrl_PTV, Case_Mis3, Ctrl_Mis3, N_case, N_ctrl, Weight. Weight = effective case excess (LGD + Dmis MPC>3)."),
    ("GeneWeights_Membership", "Membership gene weights (UKBB, CNV, controls). Columns: Trait, Gene_Entrez, Gene_Symbol, Weight. Weight is absolute beta (UKBB) or uniform=1 (CNV, controls)."),
    ("Bias_Long", "Cell-type-level mutation bias for each trait with BOTH null models. Columns: Trait, CellType, Class, Supercluster, Subtype, EFFECT, P_random, q_random, P_matched, q_matched. Each (trait, cell type cluster) is one row."),
    ("Bias_Cluster_Wide", "Pivot of EFFECT values (cell type cluster x trait). Includes cell type annotation columns."),
    ("Pvalue_Wide", "Same layout as Bias_Cluster_Wide but with P_random (empirical P-values from random null)."),
    ("Pvalue_Matched_Wide", "Same layout but with P_matched (empirical P-values from matched null)."),
    ("Qvalue_Wide", "Same layout but with q_random (FDR-corrected q-values from random null, Benjamini-Hochberg)."),
    ("Qvalue_Matched_Wide", "Same layout but with q_matched (FDR-corrected q-values from matched null)."),
    ("Bias_Supercluster", "Long-format supercluster-level bias summaries with both null models (P_random, q_random, P_matched, q_matched)."),
    ("", ""),
    ("NULL MODELS", ""),
    ("-----------", ""),
    ("Random null:", "10,000 random gene sets of matched size. Tests whether bias is greater than expected by chance for any gene set of that size."),
    ("Matched null:", "10,000 gene sets matched on whole-brain expression (WB), evolutionary conservation (mean_phastCons), and coding sequence length (n_CDS_bases) using best-of-1000 sampling. Tests whether bias is greater than expected for gene sets with similar genomic properties."),
    ("", ""),
    ("GENE WEIGHT METHODS", ""),
    ("-------------------", ""),
    ("ASD:", "BGMR-corrected mutation excess: weight = (obs_LGD - exp_LGD) + (obs_Dmis - exp_Dmis). Dmis filtered at REVEL>0.5. Expected counts from BGMR prevel_0.5 and p_LGD rates. ASD_HIQ/LIQ split by proband IQ (>70 / <=70)."),
    ("DDD:", "BGMR-corrected mutation excess: weight = (obs_LGD - exp_LGD) + (obs_Mis - exp_Mis). All missense counted. Expected counts from BGMR p_misense and p_LGD rates."),
    ("SCZ:", "Effective case-control excess for LGD + Dmis MPC>3 (mis2 excluded). Weight = nLGD + nMis3 where nLGD = Case_PTV - Ctrl_PTV * (N_case/N_ctrl)."),
    ("UKBB:", "Absolute beta coefficients from rare-variant burden analysis (Chen et al. 2023 Nat Genet)."),
    ("CNV/Controls:", "Uniform weights (weight = 1 per gene)."),
    ("", ""),
    ("Cell type clusters: 461 cell types from Siletti et al. 2023 (Nature),", ""),
    ("Human Brain Cell Atlas; cluster IDs match the Siletti taxonomy.", ""),
]
README = pd.DataFrame(readme_lines, columns=["Field", "Description"])

# %% [markdown]
# ## Write to Excel

# %%
print(f"Writing to {OUTPUT_FILE}...")
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False)
    Trait_Index.to_excel(writer, sheet_name="Trait_Index", index=False)
    GeneWeights_DeNovo.to_excel(writer, sheet_name="GeneWeights_DeNovo", index=False)
    GeneWeights_CaseControl.to_excel(writer, sheet_name="GeneWeights_CaseControl", index=False)
    GeneWeights_Membership.to_excel(writer, sheet_name="GeneWeights_Membership", index=False)
    Bias_Long.to_excel(writer, sheet_name="Bias_Long", index=False)
    Bias_Cluster_Wide.to_excel(writer, sheet_name="Bias_Cluster_Wide", index=False)
    Pvalue_Wide.to_excel(writer, sheet_name="Pvalue_Wide", index=False)
    Pvalue_Matched_Wide.to_excel(writer, sheet_name="Pvalue_Matched_Wide", index=False)
    Qvalue_Wide.to_excel(writer, sheet_name="Qvalue_Wide", index=False)
    Qvalue_Matched_Wide.to_excel(writer, sheet_name="Qvalue_Matched_Wide", index=False)
    Bias_Supercluster.to_excel(writer, sheet_name="Bias_Supercluster", index=False)

print(f"Done. Saved to {OUTPUT_FILE}")
print(f"  Size: {OUTPUT_FILE.stat().st_size / 1e6:.1f} MB")

# %% [markdown]
# ## Sanity checks

# %%
print("=== Sanity checks ===")
print(f"Trait count: {len(Trait_Index)}")

# Gene weight sheets cover all traits
gw_denovo_traits = set(GeneWeights_DeNovo["Trait"].unique())
gw_cc_traits = set(GeneWeights_CaseControl["Trait"].unique())
gw_mem_traits = set(GeneWeights_Membership["Trait"].unique())
all_gw_traits = gw_denovo_traits | gw_cc_traits | gw_mem_traits
assert all_gw_traits == set(manifest.keys()), f"Gene weight traits mismatch: {all_gw_traits ^ set(manifest.keys())}"
print(f"Gene weight sheets cover all {len(manifest)} traits: OK")
print(f"  DeNovo: {len(GeneWeights_DeNovo)} rows ({len(gw_denovo_traits)} traits)")
print(f"  CaseControl: {len(GeneWeights_CaseControl)} rows ({len(gw_cc_traits)} traits)")
print(f"  Membership: {len(GeneWeights_Membership)} rows ({len(gw_mem_traits)} traits)")

print(f"Bias rows: {len(Bias_Long)}")
expected_bias = len(manifest) * 461
assert len(Bias_Long) == expected_bias, f"Mismatch: {len(Bias_Long)} vs {expected_bias}"
print(f"  Matches {len(manifest)} traits x 461 clusters ({expected_bias}): OK")

# Check no NaN in EFFECT
n_nan = Bias_Long["EFFECT"].isna().sum()
print(f"NaN EFFECT values: {n_nan}")
assert n_nan == 0, "Found NaN EFFECT values"

# Check P-values are reasonable (both null models)
for pcol in ["P_random", "P_matched"]:
    print(f"{pcol} range: [{Bias_Long[pcol].min():.2g}, {Bias_Long[pcol].max():.2g}]")
    assert (Bias_Long[pcol] >= 0).all() and (Bias_Long[pcol] <= 1).all()
print("All checks passed.")
