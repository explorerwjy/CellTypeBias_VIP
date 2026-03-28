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
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Compute Bias Contrasts
#
# Pre-computes all pairwise bias contrasts, VIP classification, and BrainSpan
# expression data used by the figure notebooks.
#
# **Outputs** are saved to `results/main_results/contrasts/`.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os

from pathlib import Path
import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

# %% [markdown]
# ## 1. Load all bias DataFrames

# %%
Bias_Save_Dir = str(PROJ_DIR / "results/main_results/random/Centering/") + "/"

ASD_All_Bias = pd.read_csv(Bias_Save_Dir + "ASD_All_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)

VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)

X22q_mousemodel = pd.read_csv(Bias_Save_Dir + "22q_small_del_bias_addP.csv", index_col=0)

print(f"Loaded {11} bias DataFrames")

# %% [markdown]
# ## 2. Pairwise contrasts

# %%
OUT_DIR = str(PROJ_DIR / "results/main_results/contrasts/") + "/"
os.makedirs(OUT_DIR, exist_ok=True)

EffLabel = "EFFECT"

# ASD w/o ID vs SCZ
ASD_SCZ_Contrast = compare_biases(HighIQ_ASD_Bias, SCZ_Bias, name1="ASD w/o ID", name2="SCZ", efflabel=EffLabel)
ASD_SCZ_Contrast.to_csv(OUT_DIR + "ASD_woID_vs_SCZ_contrast.csv")
print("Saved ASD_woID_vs_SCZ_contrast.csv")

# ASD with ID vs SCZ
ASD_wID_SCZ_Contrast = compare_biases(LowIQ_ASD_Bias, SCZ_Bias, name1="ASD with ID", name2="SCZ", efflabel=EffLabel, neurons=Neurons)
ASD_wID_SCZ_Contrast.to_csv(OUT_DIR + "ASD_wID_vs_SCZ_contrast.csv")
print("Saved ASD_wID_vs_SCZ_contrast.csv")

# ASD w/o ID vs ASD with ID
ASD_woID_wID_Contrast = compare_biases(HighIQ_ASD_Bias, LowIQ_ASD_Bias, name1="ASD w/o ID", name2="ASD with ID", efflabel=EffLabel, neurons=ALL_CTs)
ASD_woID_wID_Contrast.to_csv(OUT_DIR + "ASD_woID_vs_ASD_wID_contrast.csv")
print("Saved ASD_woID_vs_ASD_wID_contrast.csv")

# ASD w/o ID vs DDD
ASD_woID_DDD_Contrast = compare_biases(HighIQ_ASD_Bias, DDD_Bias, name1="ASD w/o ID", name2="DD/ID", efflabel=EffLabel, neurons=Neurons)
ASD_woID_DDD_Contrast.to_csv(OUT_DIR + "ASD_woID_vs_DDD_contrast.csv")
print("Saved ASD_woID_vs_DDD_contrast.csv")

# SCZ vs ASD with ID
SCZ_ASD_wID_Contrast = compare_biases(SCZ_Bias, LowIQ_ASD_Bias, name1="SCZ", name2="ASD with ID", efflabel=EffLabel, neurons=Neurons)
SCZ_ASD_wID_Contrast.to_csv(OUT_DIR + "SCZ_vs_ASD_wID_contrast.csv")
print("Saved SCZ_vs_ASD_wID_contrast.csv")

# VNR- vs VNR+
VNR_Contrast = compare_biases(VNR_Neg_Bias, VNR_Pos_Bias, name1="VNR-", name2="VNR+", neurons=ALL_CTs)
VNR_Contrast.to_csv(OUT_DIR + "VNR_neg_vs_pos_contrast.csv")
print("Saved VNR_neg_vs_pos_contrast.csv")

# EDU- vs EDU+
EDU_Contrast = compare_biases(EDU_Neg_Bias, EDU_Pos_Bias, name1="EDU-", name2="EDU+", neurons=ALL_CTs)
EDU_Contrast.to_csv(OUT_DIR + "EDU_neg_vs_pos_contrast.csv")
print("Saved EDU_neg_vs_pos_contrast.csv")

# %% [markdown]
# ## 3. Combined contrasts + global FDR

# %%
# Filter to neurons for contrasts that need it
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
ASD_woID_wID_Contrast_Neurons = ASD_woID_wID_Contrast[ASD_woID_wID_Contrast.index.isin(Neurons)]
SCZ_ASD_wID_Contrast_Neurons = SCZ_ASD_wID_Contrast[SCZ_ASD_wID_Contrast.index.isin(Neurons)]
ASD_woID_DDD_Contrast_Neurons = ASD_woID_DDD_Contrast[ASD_woID_DDD_Contrast.index.isin(Neurons)]

Contrast_List = [VNR_Contrast, ASD_SCZ_Contrast, ASD_woID_wID_Contrast_Neurons,
                 SCZ_ASD_wID_Contrast_Neurons, ASD_woID_DDD_Contrast_Neurons]
Test = [("VNR+", "VNR-"), ("ASD w/o ID", "SCZ"), ("ASD w/o ID", "ASD with ID"),
        ("SCZ", "ASD with ID"), ("ASD w/o ID", "DD/ID")]

all_contrasts_rows = []
for pair, DF in zip(Test, Contrast_List):
    Pairname = pair[0] + " - " + pair[1]
    for supercluster, row in DF.iterrows():
        all_contrasts_rows.append([
            Pairname, supercluster,
            row["Bias_" + pair[0]], row["Bias_" + pair[1]],
            row["Bias_Diff"], row["Mann_Whitney_P"], row["Mann_Whitney_FDR"]
        ])

all_contrasts_df = pd.DataFrame(
    all_contrasts_rows,
    columns=["Pair", "SuperCluster", "Bias1", "Bias2", "BiasDiff", "Pval", "MWU_FDR"]
)
all_contrasts_df["ALL_FDR"] = fdrcorrection(all_contrasts_df["Pval"])[1]

all_contrasts_df.to_csv(OUT_DIR + "all_contrasts_fdr.csv", index=False)
print(f"Saved all_contrasts_fdr.csv ({len(all_contrasts_df)} rows)")

# %% [markdown]
# ## 4. Cluster-level bias diffs for Fig 3B

# %%
from scipy.stats import mannwhitneyu

datasets = {
    'ASD w/o ID': HighIQ_ASD_Bias,
    'ASD with ID': LowIQ_ASD_Bias,
    'VNR+': VNR_Pos_Bias,
    'VNR-': VNR_Neg_Bias,
    'DD/ID': DDD_Bias,
    'SCZ': SCZ_Bias
}

superclusters = ["CGE interneuron", "MGE interneuron", "LAMP5-LHX6 and Chandelier"]
disorder_pairs = [
    ("SCZ", "ASD with ID"),
    ("SCZ", "ASD w/o ID"),
    ("ASD with ID", "ASD w/o ID"),
    ("DD/ID", "ASD w/o ID"),
    ("VNR-", "VNR+"),
]

def disorderpair_label(d1, d2):
    return f"{d1} to {d2}"

all_biasdiffs = []
for d1, d2 in disorder_pairs:
    for sc in superclusters:
        clusters1 = datasets[d1].loc[datasets[d1]["Supercluster"] == sc].index.values
        clusters2 = datasets[d2].loc[datasets[d2]["Supercluster"] == sc].index.values
        common_clusters = set(clusters1).intersection(set(clusters2))
        for clust in common_clusters:
            val1 = datasets[d1].loc[(datasets[d1]["Supercluster"] == sc)].loc[[clust], "EFFECT"]
            val2 = datasets[d2].loc[(datasets[d2]["Supercluster"] == sc)].loc[[clust], "EFFECT"]
            if not val1.empty and not val2.empty:
                biasdiff = val1.values[0] - val2.values[0]
                all_biasdiffs.append({
                    "DisorderPair": disorderpair_label(d1, d2),
                    "Supercluster": sc,
                    "Cluster": clust,
                    "BiasDiff": biasdiff
                })

plot_df = pd.DataFrame(all_biasdiffs)
plot_df.to_csv(OUT_DIR + "cluster_biasdiff.csv", index=False)
print(f"Saved cluster_biasdiff.csv ({len(plot_df)} rows)")

# %%
# MWU tests between superclusters for each disorder pair
disorderpair_order = [disorderpair_label(d1, d2) for d1, d2 in disorder_pairs]

test_results = []
for dp in disorderpair_order:
    cge_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "CGE interneuron")]["BiasDiff"].values
    mge_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "MGE interneuron")]["BiasDiff"].values
    lamp_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "LAMP5-LHX6 and Chandelier")]["BiasDiff"].values

    # CGE vs MGE (two-tailed)
    if len(cge_diffs) > 0 and len(mge_diffs) > 0:
        stat, pval = mannwhitneyu(cge_diffs, mge_diffs, alternative='two-sided')
        test_results.append({
            "DisorderPair": dp, "Test": "CGE != MGE",
            "CGE_median": np.median(cge_diffs), "MGE_median": np.median(mge_diffs),
            "MWU_stat": stat, "MWU_pval": pval
        })
        # One-tailed MGE > CGE for SCZ to ASD with ID
        if "SCZ" in dp and "ASD with ID" in dp:
            stat_gt, pval_gt = mannwhitneyu(mge_diffs, cge_diffs, alternative='greater')
            test_results.append({
                "DisorderPair": dp, "Test": "MGE > CGE",
                "MGE_median": np.median(mge_diffs), "CGE_median": np.median(cge_diffs),
                "MWU_stat": stat_gt, "MWU_pval": pval_gt
            })

    # CGE vs LAMP5 (two-tailed)
    if len(cge_diffs) > 0 and len(lamp_diffs) > 0:
        stat, pval = mannwhitneyu(cge_diffs, lamp_diffs, alternative='two-sided')
        test_results.append({
            "DisorderPair": dp, "Test": "CGE != LAMP5",
            "CGE_median": np.median(cge_diffs), "LAMP5_median": np.median(lamp_diffs),
            "MWU_stat": stat, "MWU_pval": pval
        })

test_df = pd.DataFrame(test_results)
test_df.to_csv(OUT_DIR + "cluster_biasdiff_tests.csv", index=False)
print(f"Saved cluster_biasdiff_tests.csv ({len(test_df)} rows)")

# %% [markdown]
# ## 5. VIP classification — CGE annotation DataFrames

# %%
ExpL = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv", index_col=0)
ExpL.columns = [int(x) for x in ExpL.columns.values]

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

CGE_idx = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values

# Build CGE annotation with all disorder biases
records = []
for idx in CGE_idx:
    records.append({
        "Index": idx,
        "PV": ExpL.loc[GeneSymbol2Entrez["PVALB"], idx],
        "SST": ExpL.loc[GeneSymbol2Entrez["SST"], idx],
        "VIP": ExpL.loc[GeneSymbol2Entrez["VIP"], idx],
        "CCK": ExpL.loc[GeneSymbol2Entrez["CCK"], idx],
        "GABA": ExpL.loc[GeneSymbol2Entrez["SLC32A1"], idx],
        "CR": ExpL.loc[GeneSymbol2Entrez["CALB2"], idx],
        "EFFECT_HIQ ASD": HighIQ_ASD_Bias.loc[idx, "EFFECT"],
        "EFFECT_22q.11": X22q_Bias.loc[idx, "EFFECT"],
        "EFFECT_SCZ": SCZ_Bias.loc[idx, "EFFECT"],
        "EFFECT_VNR": VNR_Neg_Bias.loc[idx, "EFFECT"],
        "EFFECT_DD": DDD_Bias.loc[idx, "EFFECT"],
        "EFFECT_LIQ": LowIQ_ASD_Bias.loc[idx, "EFFECT"],
        "EFFECT_EDU": EDU_Neg_Bias.loc[idx, "EFFECT"],
    })
CGE_anno = pd.DataFrame(records).set_index("Index")
CGE_anno.to_csv(OUT_DIR + "CGE_VIP_annotation.csv")
print(f"Saved CGE_VIP_annotation.csv ({len(CGE_anno)} rows)")

# %%
# Build 22q-specific CGE annotation (human + mouse)
records_22q = []
for idx in CGE_idx:
    records_22q.append({
        "Index": idx,
        "PV": ExpL.loc[GeneSymbol2Entrez["PVALB"], idx],
        "SST": ExpL.loc[GeneSymbol2Entrez["SST"], idx],
        "VIP": ExpL.loc[GeneSymbol2Entrez["VIP"], idx],
        "CCK": ExpL.loc[GeneSymbol2Entrez["CCK"], idx],
        "GABA": ExpL.loc[GeneSymbol2Entrez["SLC32A1"], idx],
        "CR": ExpL.loc[GeneSymbol2Entrez["CALB2"], idx],
        "EFFECT_22q.11": X22q_Bias.loc[idx, "EFFECT"],
        "EFFECT_22q.11_mouse_gene": X22q_mousemodel.loc[idx, "EFFECT"],
    })
CGE_anno_22q = pd.DataFrame(records_22q).set_index("Index")

# Add VIP+/- flag
cutoff = 1.0
CGE_anno_22q["VIP_group"] = np.where(CGE_anno_22q["VIP"] >= cutoff, "VIP+", "VIP-")
CGE_anno_22q.to_csv(OUT_DIR + "CGE_VIP_annotation_22q.csv")
print(f"Saved CGE_VIP_annotation_22q.csv ({len(CGE_anno_22q)} rows)")

# %%
# Write VIP+/VIP- cluster lists (side effect kept for compatibility)
CGE_VIP_Pos = CGE_anno_22q[CGE_anno_22q["VIP"] >= cutoff].index.values
CGE_VIP_Neg = CGE_anno_22q[CGE_anno_22q["VIP"] < cutoff].index.values

os.makedirs(str(PROJ_DIR / "dat/Other"), exist_ok=True)
with open(str(PROJ_DIR / "dat/Other/CGE_VIP_Pos.txt"), "w") as f:
    for gene in CGE_VIP_Pos:
        f.write(str(gene) + "\n")
with open(str(PROJ_DIR / "dat/Other/CGE_VIP_Neg.txt"), "w") as f:
    for gene in CGE_VIP_Neg:
        f.write(str(gene) + "\n")
print(f"VIP+ clusters: {len(CGE_VIP_Pos)}, VIP- clusters: {len(CGE_VIP_Neg)}")

# %% [markdown]
# ## 6. BrainSpan expression data

# %%
ExpMat = pd.read_csv("/home/jw3514/Work/BrainDisorders/data/expression/brainspan/gene_matrix/gene_exp_avg2time.csv", index_col=0)
LogExpMat = np.log2(ExpMat + 1)
qnLogExpMat = quantileNormalize(LogExpMat)
Meta = pd.read_csv("/home/jw3514/Work/BrainDisorders/data/expression/brainspan/gene_matrix/rows_metadata.csv", index_col=0)
Meta = Meta[~Meta["entrez_id"].isna()]
Meta["entrez_id"] = [int(x) for x in Meta["entrez_id"].values]

Time = ['mean_2A', 'mean_2B', 'mean_3A', 'mean_3B', 'mean_4', 'mean_5',
        'mean_6', 'mean_7', 'mean_8', 'mean_9', 'mean_10', 'mean_11']

# Load gene sets
GeneWeightDIR = str(PROJ_DIR / "dat/GeneWeights/") + "/"
SCZ_Genes = pd.read_csv(f"{GeneWeightDIR}/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw", header=None)[0].values
HIQ_ASD_Genes = pd.read_csv(f"{GeneWeightDIR}/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw", header=None)[0].values
LIQ_ASD_Genes = pd.read_csv(f"{GeneWeightDIR}/LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw", header=None)[0].values

# Filter to gene sets
SCZ_Idx = Meta[Meta["entrez_id"].isin(SCZ_Genes)].index
HIQASD_Idx = Meta[Meta["entrez_id"].isin(HIQ_ASD_Genes)].index
LIQASD_Idx = Meta[Meta["entrez_id"].isin(LIQ_ASD_Genes)].index

SCZ_Dat = qnLogExpMat.loc[SCZ_Idx, Time]
HIQASD_Dat = qnLogExpMat.loc[HIQASD_Idx, Time]
LIQASD_Dat = qnLogExpMat.loc[LIQASD_Idx, Time]
ALL_Dat = qnLogExpMat.loc[:, Time]

# Build summary table: mean and SEM per gene set per developmental stage
brainspan_records = []
for gene_set_name, dat in [("ASD w/o ID", HIQASD_Dat), ("ASD with ID", LIQASD_Dat),
                            ("SCZ", SCZ_Dat), ("ALL", ALL_Dat)]:
    for col in Time:
        brainspan_records.append({
            "gene_set": gene_set_name,
            "stage": col,
            "mean": dat[col].mean(),
            "sem": dat[col].std() / math.sqrt(dat.shape[0]),
        })

brainspan_df = pd.DataFrame(brainspan_records)
brainspan_df.to_csv(OUT_DIR + "brainspan_expression.csv", index=False)
print(f"Saved brainspan_expression.csv ({len(brainspan_df)} rows)")

# %% [markdown]
# ## Summary
#
# All pre-computed files are now in `results/main_results/contrasts/`.

# %%
import glob
saved_files = sorted(glob.glob(OUT_DIR + "*.csv"))
print(f"\n{len(saved_files)} files in {OUT_DIR}:")
for f in saved_files:
    print(f"  {os.path.basename(f)}")
