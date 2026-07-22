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
# # Assemble Supplementary Table: Mutation Bias Results
#
# Reviewer 3 minor 1: *"Please provide full mutation bias estimates for all
# traits tested, not just selected highlights. Preferably provide all mutational
# bias results for all traits tested in a 'long' format, with additional column
# describing the tested trait. Gene lists and mutation counts used for each
# trait should be published in supplementary tables."*
#
# This notebook reads `dat/SuppTable_MutationBias_manifest.yaml` and assembles
# a single multi-sheet Excel workbook with the following sheets:
#
# 1. **README** — sheet descriptions, null-model definitions, weighting methods
# 2. **Trait_Index** — one row per trait with display name, category, source,
#    figure references, gene count, # significant cell types
# 3. **GeneWeights_Master** — wide table indexed by Gene_Symbol. Per-trait
#    blocks of (membership, mutation counts / β, weight) columns. NA when the
#    gene is not in the trait's gene set.
# 4. **Bias_AllTraits** — wide table indexed by CellType (cluster_id) with
#    hierarchical columns: top level = "CellType_Annotation" or trait name;
#    sub level = annotation fields once, then (bias, Rank, P_random, q_random,
#    P_matched, q_matched) per trait.
# 5. **Contrast_***— 5 supercluster-level bias-contrast sheets:
#    SCZ vs ASD w/ ID, ASD w/o ID vs SCZ, ASD w/o ID vs ASD w/ ID,
#    ASD w/o ID vs DDD, VNR- vs VNR+.
#
# Output: `results/SuppTable_MutationBias.xlsx`

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import warnings; warnings.filterwarnings('ignore')
import yaml
from pathlib import Path
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

# Expression matrix used for bias calc (to flag In_ExpMat)
EXPMAT_FILE = PROJ_DIR / _cfg["analysis_types"]["Centering"]

# %% [markdown]
# ## Load manifest

# %%
with open(MANIFEST_FILE) as f:
    manifest = yaml.safe_load(f)

# Preserve manifest order
trait_order = list(manifest.keys())
print(f"Loaded manifest with {len(manifest)} traits")
for trait, info in manifest.items():
    print(f"  {trait:20s}  {info['category']:25s}  {info['display_name']}")

# %% [markdown]
# ## Helpers

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
    for key in ("weight_file", "bias_file", "bias_file_matched"):
        p = PROJ_DIR / info[key]
        if not p.exists():
            missing.append((trait, key, str(p)))
if missing:
    print("Missing files:")
    for m in missing:
        print(f"  {m}")
    raise FileNotFoundError(f"{len(missing)} required files missing")
print("All manifest files exist.")

# %% [markdown]
# ## Load expression-matrix gene set for the In_ExpMat flag
#
# We only need the row index of the expression matrix to know which genes
# are eligible to contribute to bias.

# %%
expmat_genes = pd.read_csv(EXPMAT_FILE, usecols=[0], index_col=0).index
expmat_entrez_set = set(int(g) for g in expmat_genes if str(g).isdigit() or
                        (isinstance(g, (int, float)) and not pd.isna(g)))
# Some indices may be strings; coerce robustly
expmat_entrez_set = set()
for g in expmat_genes:
    try:
        expmat_entrez_set.add(int(g))
    except (ValueError, TypeError):
        continue
print(f"Expression matrix has {len(expmat_entrez_set)} unique Entrez IDs (gene rows)")

# %% [markdown]
# ## Load source mutation tables (ASD SPARK, DDD, SCZ case-control)

# %%
# ASD SPARK Table S7 (full sorted gene list with per-gene LoF / Dmis counts)
Spark_Denovo = pd.read_excel(
    str(ASD_SPARK_FILE), skiprows=2, sheet_name="Table S7",
)
Spark_Denovo = Spark_Denovo[Spark_Denovo["pDenovoWEST_Meta"] != "."]
Spark_Denovo["pDenovoWEST_Meta"] = Spark_Denovo["pDenovoWEST_Meta"].astype(float)
Spark_Denovo = Spark_Denovo.sort_values("pDenovoWEST_Meta").reset_index(drop=True)
Spark_Denovo = Spark_Denovo.set_index("EntrezID")

# IQ-stratified counts for HIQ/LIQ
Mut_n_IQ = pd.read_csv(str(ASD_IQ_FILE))
HighIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] > 70]
LowIQMuts_ALL = Mut_n_IQ[Mut_n_IQ["IQ"] <= 70]

lgd_types = ["frameshift", "splice_acceptor", "splice_donor",
             "start_lost", "stop_gained", "stop_lost"]
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

# DDD (Kaplanis 2020) full sorted gene list
DDD_df = pd.read_excel(str(DDD_FILE))
DDD_df = DDD_df.sort_values("denovoWEST_p_full")
DDD_df["EntrezID"] = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_df["symbol"].values]
DDD_df = DDD_df[DDD_df["EntrezID"] > 0]
DDD_df = DDD_df.set_index("EntrezID")

# SCZ SCHEMA case-control gene counts
SCZ_GeneDF = pd.read_csv(str(SCZ_MUT_FILE), index_col=0)

print(f"Source loaded:  Spark={len(Spark_Denovo)},  DDD={len(DDD_df)},  SCZ={len(SCZ_GeneDF)}")

# %% [markdown]
# ## Build GeneWeights_Master (wide; one row per gene)
#
# - Index: Gene_Symbol (union across all 16 traits' gene sets)
# - Per-trait block conventions:
#   * `{trait}_in`              — 1 if gene in this trait's set, NA otherwise
#   * `{trait}_nLGD/nDmis/nMis` — observed counts (de novo)
#   * `{trait}_Case_PTV` etc.   — case/ctrl counts (SCZ)
#   * `{trait}_weight`          — weight used in bias calc (from the gw file)
# - `In_ExpMat` (1/0) flags whether the gene is present in the expression matrix
#   used for bias calculation. Genes with In_ExpMat=0 contribute 0 to bias even
#   if they appear in a trait's source gene set (e.g. ~4/46 22q11.2 genes).

# %%
# Helper: build a single-trait row block for a given gene
def block_columns_for(trait):
    """Return ordered list of (col_suffix, label) for this trait's block."""
    cat = manifest[trait]["category"]
    gwtype = manifest[trait].get("gene_weight_type", "membership")
    cols = ["in"]
    if trait in ("ASD_All", "ASD_HIQ", "ASD_LIQ"):
        cols += ["nLGD", "nDmis"]
    elif trait in ("DDD_61", "DDD_285"):
        cols += ["nLGD", "nMis"]
    elif trait == "SCZ":
        cols += ["Case_PTV", "Ctrl_PTV", "Case_Mis3", "Ctrl_Mis3"]
    cols += ["weight"]
    return cols


# Build per-trait records keyed by Entrez ID
per_trait_records = {}  # trait -> {entrez_id: {field: value}}
trait_gene_lists = {}   # trait -> set(entrez)

for trait in trait_order:
    gw = load_geneweight_file(PROJ_DIR / manifest[trait]["weight_file"])
    trait_gene_lists[trait] = set(int(e) for e in gw["Entrez"].values)
    cols = block_columns_for(trait)
    rec = {}
    for _, r in gw.iterrows():
        eid = int(r["Entrez"])
        row = {"in": 1}
        # Mutation count columns
        if trait == "ASD_All":
            if eid in Spark_Denovo.index:
                row["nLGD"] = float(Spark_Denovo.loc[eid, "AutismMerged_LoF"])
                row["nDmis"] = float(Spark_Denovo.loc[eid, "AutismMerged_Dmis_REVEL0.5"])
            else:
                row["nLGD"] = np.nan
                row["nDmis"] = np.nan
        elif trait == "ASD_HIQ":
            if eid in Spark_Denovo.index:
                row["nLGD"] = float(Spark_Denovo.loc[eid, "HIQ_LGD"])
                row["nDmis"] = float(Spark_Denovo.loc[eid, "HIQ_Dmis"])
            else:
                row["nLGD"] = np.nan
                row["nDmis"] = np.nan
        elif trait == "ASD_LIQ":
            if eid in Spark_Denovo.index:
                row["nLGD"] = float(Spark_Denovo.loc[eid, "LIQ_LGD"])
                row["nDmis"] = float(Spark_Denovo.loc[eid, "LIQ_Dmis"])
            else:
                row["nLGD"] = np.nan
                row["nDmis"] = np.nan
        elif trait in ("DDD_61", "DDD_285"):
            if eid in DDD_df.index:
                nLGD = (DDD_df.loc[eid, "frameshift_variant"]
                        + DDD_df.loc[eid, "splice_acceptor_variant"]
                        + DDD_df.loc[eid, "splice_donor_variant"]
                        + DDD_df.loc[eid, "stop_gained"]
                        + DDD_df.loc[eid, "stop_lost"])
                row["nLGD"] = float(nLGD)
                row["nMis"] = float(DDD_df.loc[eid, "missense_variant"])
            else:
                row["nLGD"] = np.nan
                row["nMis"] = np.nan
        elif trait == "SCZ":
            if eid in SCZ_GeneDF.index:
                src = SCZ_GeneDF.loc[eid]
                row["Case_PTV"] = float(src.get("Case PTV", 0) or 0)
                row["Ctrl_PTV"] = float(src.get("Ctrl PTV", 0) or 0)
                row["Case_Mis3"] = float(src.get("Case mis3", 0) or 0)
                row["Ctrl_Mis3"] = float(src.get("Ctrl mis3", 0) or 0)
            else:
                row["Case_PTV"] = np.nan
                row["Ctrl_PTV"] = np.nan
                row["Case_Mis3"] = np.nan
                row["Ctrl_Mis3"] = np.nan
        row["weight"] = float(r["Weight"])
        rec[eid] = row
    per_trait_records[trait] = rec

# Union of Entrez IDs across all traits
all_entrez = set()
for s in trait_gene_lists.values():
    all_entrez |= s
all_entrez = sorted(all_entrez)
print(f"Union of gene sets: {len(all_entrez)} unique Entrez IDs")

# %%
# Assemble the master wide table
rows = []
for eid in all_entrez:
    row = {
        "Gene_Symbol": Entrez2Symbol.get(eid, ""),
        "Entrez_ID": eid,
        "In_ExpMat": 1 if eid in expmat_entrez_set else 0,
    }
    for trait in trait_order:
        cols = block_columns_for(trait)
        if eid in per_trait_records[trait]:
            for c in cols:
                row[f"{trait}_{c}"] = per_trait_records[trait][eid].get(c, np.nan)
        else:
            for c in cols:
                row[f"{trait}_{c}"] = np.nan
    rows.append(row)

GeneWeights_Master = pd.DataFrame(rows)
# Sort: in_ExpMat first (so dropped genes are visible at the bottom),
# then by gene symbol
GeneWeights_Master = GeneWeights_Master.sort_values(
    ["In_ExpMat", "Gene_Symbol"], ascending=[False, True]
).reset_index(drop=True)
print(f"GeneWeights_Master shape: {GeneWeights_Master.shape}")
print(f"  Genes in ExpMat: {(GeneWeights_Master['In_ExpMat']==1).sum()}")
print(f"  Genes NOT in ExpMat (dropped from bias calc): {(GeneWeights_Master['In_ExpMat']==0).sum()}")

# %% [markdown]
# ## Build Bias_AllTraits (hierarchical wide)
#
# - Row index: `CellType` (cluster_id from Siletti taxonomy)
# - Column MultiIndex level 0: "CellType_Annotation" or trait name
# - Level 1 (annotation): Class, Supercluster, Subtype, Neurotransmitter,
#   Top three regions, Top three dissections, Number of cells,
#   Neuropeptide auto-annotation
# - Level 1 (per trait): bias, Rank, P_random, q_random, P_matched, q_matched

# %%
anno_cols_in_order = [
    "Class", "Supercluster", "Subtype", "Neurotransmitter",
    "Top three regions", "Top three dissections", "Number of cells",
    "Neuropeptide auto-annotation",
]

# Per-trait bias CSVs are sorted by EFFECT (descending), so cluster order
# differs across traits. We reindex everything on cluster_id (ascending) for
# a stable canonical row order.
ref_trait = "ASD_All"
ref_bias = pd.read_csv(PROJ_DIR / manifest[ref_trait]["bias_file"], index_col=0)
ref_bias.index.name = "CellType"
# Canonical cluster index = sorted by cluster_id (numeric ascending)
celltype_idx = sorted(ref_bias.index.tolist())

anno_present = [c for c in anno_cols_in_order if c in ref_bias.columns]
anno_df = ref_bias[anno_present].reindex(celltype_idx).copy()

# Per-trait block columns
trait_block_cols = ["bias", "Rank", "P_random", "q_random", "P_matched", "q_matched"]

# Collect per-trait DataFrames aligned on canonical index
trait_frames = {}
for trait in trait_order:
    bias_r = pd.read_csv(PROJ_DIR / manifest[trait]["bias_file"], index_col=0)
    bias_m = pd.read_csv(PROJ_DIR / manifest[trait]["bias_file_matched"], index_col=0)
    if sorted(bias_r.index) != celltype_idx:
        raise ValueError(f"Unexpected cluster set in {trait} (random null)")
    if sorted(bias_m.index) != celltype_idx:
        raise ValueError(f"Unexpected cluster set in {trait} (matched null)")
    bias_r = bias_r.reindex(celltype_idx)
    bias_m = bias_m.reindex(celltype_idx)
    df = pd.DataFrame(index=celltype_idx)
    df.index.name = "CellType"
    df["bias"] = bias_r["EFFECT"].values
    df["Rank"] = bias_r["Rank"].values
    df["P_random"] = bias_r["P-value"].values
    df["q_random"] = bias_r["q-value"].values
    df["P_matched"] = bias_m["P-value"].values
    df["q_matched"] = bias_m["q-value"].values
    trait_frames[trait] = df

# Build hierarchical-column DataFrame
# Annotation block first
anno_part = anno_df.copy()
anno_part.columns = pd.MultiIndex.from_tuples(
    [("CellType_Annotation", c) for c in anno_part.columns]
)

trait_parts = []
for trait in trait_order:
    df = trait_frames[trait]
    df = df.copy()
    df.columns = pd.MultiIndex.from_tuples([(trait, c) for c in df.columns])
    trait_parts.append(df)

Bias_AllTraits = pd.concat([anno_part] + trait_parts, axis=1)
Bias_AllTraits.index.name = "CellType"
print(f"Bias_AllTraits shape: {Bias_AllTraits.shape}")
print(f"  Columns: {len(anno_present)} annotation + "
      f"{len(trait_order)} × {len(trait_block_cols)} per-trait "
      f"= {len(Bias_AllTraits.columns)}")

# %% [markdown]
# ## Load supercluster-level bias contrasts (5 sheets)

# %%
CONTRAST_DIR = PROJ_DIR / "results/main_results/contrasts"
contrast_specs = [
    ("Contrast_SCZ_vs_ASDwID",       "SCZ_vs_ASD_wID_contrast.csv"),
    ("Contrast_ASDwoID_vs_SCZ",      "ASD_woID_vs_SCZ_contrast.csv"),
    ("Contrast_ASDwoID_vs_ASDwID",   "ASD_woID_vs_ASD_wID_contrast.csv"),
    ("Contrast_ASDwoID_vs_DDD",      "ASD_woID_vs_DDD_contrast.csv"),
    ("Contrast_VNRneg_vs_VNRpos",    "VNR_neg_vs_pos_contrast.csv"),
]
contrast_frames = {}
for sheet, fname in contrast_specs:
    p = CONTRAST_DIR / fname
    if not p.exists():
        raise FileNotFoundError(f"Missing contrast file: {p}")
    contrast_frames[sheet] = pd.read_csv(p)
    print(f"  {sheet:35s} {p.name:45s} {contrast_frames[sheet].shape}")

# %% [markdown]
# ## Build Trait_Index

# %%
trait_index_rows = []
for trait in trait_order:
    info = manifest[trait]
    gw = load_geneweight_file(PROJ_DIR / info["weight_file"])
    bias_r = pd.read_csv(PROJ_DIR / info["bias_file"], index_col=0)
    trait_index_rows.append({
        "Trait": trait,
        "Display_Name": info["display_name"],
        "Category": info["category"],
        "Weight_Format": info["weight_format"],
        "Gene_Weight_Type": info.get("gene_weight_type", "membership"),
        "N_Genes_Source": len(gw),
        "N_Genes_In_ExpMat": int(sum(1 for e in gw["Entrez"].values
                                     if int(e) in expmat_entrez_set)),
        "Weight_Min": round(gw["Weight"].min(), 6),
        "Weight_Max": round(gw["Weight"].max(), 6),
        "Weight_Mean": round(gw["Weight"].mean(), 6),
        "N_CellTypes": len(bias_r),
        "N_Sig_q05_random": int((bias_r["q-value"] < 0.05).sum()),
        "N_Sig_q10_random": int((bias_r["q-value"] < 0.10).sum()),
        "Source": info["source"],
        "Figure_Reference": info["figure_ref"],
    })

Trait_Index = pd.DataFrame(trait_index_rows)
print(Trait_Index[["Trait", "Category", "N_Genes_Source",
                   "N_Genes_In_ExpMat", "N_Sig_q05_random"]].to_string(index=False))

# %% [markdown]
# ## README sheet

# %%
readme_lines = [
    ("Supplementary Table: Mutation Bias Results", ""),
    ("", ""),
    ("Reviewer 3 minor 1: this workbook provides full mutation-bias estimates", ""),
    ("and gene weights for all traits tested in the manuscript.", ""),
    ("", ""),
    ("SHEETS", ""),
    ("------", ""),
    ("Trait_Index",
     "One row per trait. Display name, category, weighting method, gene "
     "counts, # significant cell types (random null)."),
    ("GeneWeights_Master",
     "Wide table indexed by Gene_Symbol with per-trait blocks. "
     "Columns: Entrez_ID, In_ExpMat, then per trait: {trait}_in (1 if gene in "
     "set, NA otherwise), {trait}_nLGD/nDmis/nMis or {trait}_Case_PTV etc., "
     "{trait}_weight (weight used in bias calc)."),
    ("Bias_AllTraits",
     "Cell-type-level mutation bias for every trait under both null models. "
     "Index = CellType (Siletti cluster_id). Two-row column header: top row "
     "is 'CellType_Annotation' or a trait name; sub row is the field. "
     "Annotation appears once (Class, Supercluster, Subtype, Neurotransmitter, "
     "Top three regions/dissections, Number of cells, "
     "Neuropeptide auto-annotation). Per trait: bias (=EFFECT), Rank, "
     "P_random, q_random, P_matched, q_matched."),
    ("Contrast_SCZ_vs_ASDwID",      "Supercluster-level bias contrast: SCZ - ASD with ID."),
    ("Contrast_ASDwoID_vs_SCZ",     "Supercluster-level bias contrast: ASD w/o ID - SCZ."),
    ("Contrast_ASDwoID_vs_ASDwID",  "Supercluster-level bias contrast: ASD w/o ID - ASD with ID."),
    ("Contrast_ASDwoID_vs_DDD",     "Supercluster-level bias contrast: ASD w/o ID - DDD."),
    ("Contrast_VNRneg_vs_VNRpos",   "Supercluster-level bias contrast: VNR negative vs VNR positive."),
    ("", ""),
    ("NULL MODELS", ""),
    ("-----------", ""),
    ("Random null",
     "10,000 random gene sets of matched size. Tests whether bias is greater "
     "than expected by chance for any gene set of that size."),
    ("Matched null",
     "10,000 gene sets matched on whole-brain expression (WB), evolutionary "
     "conservation (mean_phastCons), and coding sequence length (n_CDS_bases) "
     "using best-of-1000 sampling. Tests whether bias is greater than "
     "expected for gene sets with similar genomic properties."),
    ("", ""),
    ("WEIGHTING SCHEMES (per trait family)", ""),
    ("------------------------------------", ""),
    ("ASD (All / HIQ / LIQ)",
     "BGMR-corrected de novo excess: weight = (nLGD - exp_LGD) + "
     "(nDmis - exp_Dmis), floored at 0. Dmis filtered at REVEL > 0.5. "
     "HIQ = SPARK probands with FSIQ > 70 (N=3,610); "
     "LIQ = FSIQ <= 70 (N=1,792)."),
    ("DDD (61 / 285)",
     "BGMR-corrected de novo excess: weight = (nLGD - exp_LGD) + "
     "(nMis - exp_Mis). All missense counted (Kaplanis et al. did not use a "
     "REVEL threshold; per-individual mutation data not available to apply one)."),
    ("SCZ",
     "Effective case excess for LGD + Dmis MPC>3 (Mis2 excluded). "
     "weight = (Case_PTV + Case_Mis3) - (Ctrl_PTV + Ctrl_Mis3) * "
     "(N_case / N_ctrl). Source: SCHEMA, Singh et al. 2022."),
    ("UKBB cognitive (VNR/EDU/RT, ±)",
     "weight = abs(beta) from rare-variant burden analysis "
     "(Chen et al. 2023, Nat Genet). Top 61 genes by association direction."),
    ("22q11.2 (del / small_del)",
     "Uniform weight = 1 per gene. del = 46 genes in 3 Mb interval; "
     "small_del = 27 genes corresponding to mouse Df(16)A homolog."),
    ("NegCtrl (HDL/Alanine/RBC/IBD)",
     "Uniform weight = 1 per gene (top 61 genes by burden in UKBB Genebass). "
     "Non-brain negative controls."),
    ("", ""),
    ("In_ExpMat", "1 if Entrez ID present in the expression matrix used for "
                   "bias calculation; 0 if absent. Genes with In_ExpMat=0 are "
                   "listed in GeneWeights_Master for completeness but do not "
                   "contribute to the bias score."),
    ("", ""),
    ("Cell-type clusters: 461 cell types from Siletti et al. 2023 (Nature), "
     "Human Brain Cell Atlas; cluster IDs match the Siletti taxonomy.", ""),
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
    GeneWeights_Master.to_excel(writer, sheet_name="GeneWeights_Master", index=False)
    # Hierarchical-column sheet: index=True so CellType cluster_id is shown
    Bias_AllTraits.to_excel(writer, sheet_name="Bias_AllTraits", index=True)
    for sheet, _ in contrast_specs:
        contrast_frames[sheet].to_excel(writer, sheet_name=sheet, index=False)

print(f"Done. Saved to {OUTPUT_FILE}")
print(f"  Size: {OUTPUT_FILE.stat().st_size / 1e6:.1f} MB")

# %% [markdown]
# ## Sanity checks

# %%
print("=== Sanity checks ===")
print(f"Trait count: {len(Trait_Index)} (expected 18)")
assert len(Trait_Index) == 18

# GeneWeights_Master spot checks
print(f"\nGeneWeights_Master: {GeneWeights_Master.shape}")
for trait in trait_order:
    n = GeneWeights_Master[f"{trait}_in"].sum()
    n_expected = len(load_geneweight_file(PROJ_DIR / manifest[trait]["weight_file"]))
    print(f"  {trait:20s} membership = {int(n):4d}  (expected {n_expected})")
    assert int(n) == n_expected, f"Membership mismatch for {trait}"

# 22q size check
n22 = GeneWeights_Master["22q_del_in"].sum()
print(f"\n  22q_del genes in master table: {int(n22)} (expected 46)")
assert int(n22) == 46

# Bias_AllTraits checks
print(f"\nBias_AllTraits: {Bias_AllTraits.shape}")
n_anno = sum(1 for c in Bias_AllTraits.columns if c[0] == "CellType_Annotation")
n_per_trait = {t: sum(1 for c in Bias_AllTraits.columns if c[0] == t)
               for t in trait_order}
print(f"  Annotation columns: {n_anno} (expected {len(anno_present)})")
print(f"  Per-trait column counts: {set(n_per_trait.values())} (expected {{{len(trait_block_cols)}}})")
assert n_anno == len(anno_present)
assert all(n == len(trait_block_cols) for n in n_per_trait.values())

# bias column range checks
all_bias_vals = []
for trait in trait_order:
    vals = Bias_AllTraits[(trait, "bias")].values
    all_bias_vals.extend(vals)
    assert not np.isnan(vals).any(), f"NaN bias values in {trait}"
print(f"  Bias range across all traits: [{min(all_bias_vals):.3f}, {max(all_bias_vals):.3f}]")

# P-value range checks
for pcol in ["P_random", "P_matched"]:
    for trait in trait_order:
        v = Bias_AllTraits[(trait, pcol)].values
        assert (v >= 0).all() and (v <= 1).all(), f"Bad {pcol} in {trait}"
print(f"  P_random / P_matched all in [0, 1]: OK")

print("\nAll checks passed.")
