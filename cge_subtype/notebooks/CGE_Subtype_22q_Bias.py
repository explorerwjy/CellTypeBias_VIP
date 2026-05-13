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
# # CGE Subtype 22q Mutation Bias
#
# Reviewers 2 and 3 ask whether the CCKBC vs ISI2/ISI3/VIPo difference observed
# in mouse 22q electrophysiology is recapitulated in the human mutation bias analysis.
#
# This notebook:
# 1. Classifies 21 CGE clusters into subtypes using marker expression (VIP, CCK, CR, SNCG, M2R)
# 2. Loads Patch-seq electrophysiology from scANVI-mapped human interneurons
# 3. Compares marker-based and ephys-based cluster assignments
# 4. Tests whether CCKBC shows lower 22q bias than ISI/VIP subtypes
# 5. Extends the comparison to other cognitive disorders

# %% [markdown]
# ## 1. Setup

# %%
# %load_ext autoreload
# %autoreload 2

import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/"
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

try:
    os.chdir(f"{ProjDIR}/notebooks/")
    print(f"Working directory: {os.getcwd()}")
except Exception as e:
    print(f"Error: {e}")

# Paths
BIAS_DIR = Path(ProjDIR) / "results" / "main_results" / "random" / "Centering"
CONTRAST_DIR = Path(ProjDIR) / "results" / "main_results" / "contrasts"
EPHYS_FX_PATH = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")
MAPPING_PATH = Path("/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/human_interneuron/scvi_mapped/mapping_results.csv")
EXPR_PATH = "/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv"

# Plotting style
sns.set_style("ticks")
plt.rcParams["figure.dpi"] = 120

# %% [markdown]
# ## 2. Load Data

# %%
# CGE annotation tables (pre-computed marker + bias per cluster)
cge_anno = pd.read_csv(CONTRAST_DIR / "CGE_VIP_annotation.csv", index_col="Index")
cge_anno_22q = pd.read_csv(CONTRAST_DIR / "CGE_VIP_annotation_22q.csv", index_col="Index")

print(f"CGE clusters: {len(cge_anno)}")
print(f"Columns (main): {list(cge_anno.columns)}")
print(f"Columns (22q):  {list(cge_anno_22q.columns)}")

# %%
# Load expression matrix and log2-transform
ExpL = pd.read_csv(EXPR_PATH, index_col=0)
ExpL.columns = [int(x) for x in ExpL.columns.values]
ExpL = ExpL.map(lambda x: np.log2(x + 1))
print(f"Expression matrix: {ExpL.shape[0]} genes x {ExpL.shape[1]} cell types")

# %%
# CGE cluster indices
CGE_idx = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values
MGE_idx = Anno[Anno["Supercluster"] == "MGE interneuron"].index.values
print(f"CGE clusters: {len(CGE_idx)}, MGE clusters: {len(MGE_idx)}")

# %% [markdown]
# ## 3. Marker Expression Profiles

# %%
# Define markers and cutoffs (log2(UMI+1) scale)
markers = {
    "SST": GeneSymbol2Entrez["SST"],
    "PVALB": GeneSymbol2Entrez["PVALB"],
    "CCK": GeneSymbol2Entrez["CCK"],
    "CALB2": GeneSymbol2Entrez["CALB2"],
    "SNCG": GeneSymbol2Entrez["SNCG"],
    "VIP": GeneSymbol2Entrez["VIP"],
    "M2R": GeneSymbol2Entrez["CHRM2"],
    "CNR1": GeneSymbol2Entrez["CNR1"],
}

cutoffs = {
    "VIP": 1,
    "CCK": 2,
    "CR": 1.5,
    "SNCG": 0.5,
    "M2R": 2,
    "CNR1": 1,
}

# %%
# Histogram: CGE vs MGE expression for each marker
n_markers = len(markers)
n_cols = 2
n_rows = math.ceil(n_markers / n_cols)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows))
fig.patch.set_alpha(0)
fig.suptitle("CGE vs MGE Marker Expression (log2 UMI+1)", fontsize=14, y=1.02)
axes = axes.flatten()

for idx, (marker_name, entrez_id) in enumerate(markers.items()):
    cge_expr = ExpL.loc[entrez_id, CGE_idx]
    mge_expr = ExpL.loc[entrez_id, MGE_idx]
    ax = axes[idx]
    ax.patch.set_alpha(0)
    sns.histplot(cge_expr, bins=20, ax=ax, color="blue", alpha=0.7, label="CGE")
    sns.histplot(mge_expr, bins=20, ax=ax, color="red", alpha=0.7, label="MGE")
    ax.set_title(marker_name)
    ax.set_xlabel("Log2(UMI + 1)")
    ax.set_ylabel("Frequency")
    ax.legend()

    # Cutoff line
    cutoff_key = "CR" if marker_name == "CALB2" else marker_name
    if cutoff_key in cutoffs:
        ax.axvline(cutoffs[cutoff_key], color="k", linestyle="--", linewidth=1.5, label="Cutoff")

for idx in range(n_markers, len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Marker-Based Subtype Classification

# %%
# Build combined DataFrame: bias + marker expression for all 21 CGE clusters
# Start with multi-disorder bias from the main annotation
Biases_CGE = cge_anno.copy()

# Rename EFFECT columns for clarity
bias_rename = {
    "EFFECT_HIQ ASD": "ASD_HIQ_Bias",
    "EFFECT_22q.11": "X22q_Bias",
    "EFFECT_SCZ": "SCZ_Bias",
    "EFFECT_VNR": "VNR_Bias",
    "EFFECT_DD": "DDD_Bias",
    "EFFECT_LIQ": "ASD_LIQ_Bias",
    "EFFECT_EDU": "EDU_Bias",
}
Biases_CGE = Biases_CGE.rename(columns=bias_rename)

# Add 22q mouse-gene bias from the 22q-specific annotation
Biases_CGE["X22q_mouse_Bias"] = cge_anno_22q["EFFECT_22q.11_mouse_gene"]

# Add log2(UMI+1) marker expression per cluster
marker_genes = {
    "VIP_ExpL": "VIP",
    "CCK_ExpL": "CCK",
    "CR_ExpL": "CALB2",
    "SNCG_ExpL": "SNCG",
    "M2R_ExpL": "CHRM2",
    "CNR1_ExpL": "CNR1",
}
for colname, symbol in marker_genes.items():
    entrez_id = GeneSymbol2Entrez[symbol]
    Biases_CGE.loc[:, colname] = ExpL.loc[entrez_id, Biases_CGE.index].values

print(f"Biases_CGE shape: {Biases_CGE.shape}")
Biases_CGE[["VIP_ExpL", "CCK_ExpL", "CR_ExpL", "SNCG_ExpL", "M2R_ExpL", "CNR1_ExpL"]].round(2)

# %%
# Apply marker cutoffs to classify subtypes
# CCKBC from ALL CGE (not just VIP+): CR-, SNCG+, CCK+, M2R-
CCKBC = Biases_CGE[
    (Biases_CGE["CR_ExpL"] < cutoffs["CR"])
    & (Biases_CGE["SNCG_ExpL"] > cutoffs["SNCG"])
    & (Biases_CGE["CCK_ExpL"] > cutoffs["CCK"])
    & (Biases_CGE["M2R_ExpL"] < cutoffs["M2R"])
]

VIP_Pos = Biases_CGE[Biases_CGE["VIP_ExpL"] > cutoffs["VIP"]]
VIP_Neg = Biases_CGE[Biases_CGE["VIP_ExpL"] < cutoffs["VIP"]]

# ISI2: VIP+, CR-, M2R-, SNCG-
ISI2 = VIP_Pos[
    (VIP_Pos["CR_ExpL"] < cutoffs["CR"])
    & (VIP_Pos["M2R_ExpL"] < cutoffs["M2R"])
    & (VIP_Pos["SNCG_ExpL"] < cutoffs["SNCG"])
]

# ISI3: VIP+, CR+, M2R-, SNCG-
ISI3 = VIP_Pos[
    (VIP_Pos["CR_ExpL"] > cutoffs["CR"])
    & (VIP_Pos["M2R_ExpL"] < cutoffs["M2R"])
    & (VIP_Pos["SNCG_ExpL"] < cutoffs["SNCG"])
]

# Other VIP+: remaining VIP+ not in named groups
combined_indices = np.concatenate([CCKBC.index.values, ISI2.index.values, ISI3.index.values])
Other = VIP_Pos.loc[~VIP_Pos.index.isin(combined_indices)]

print(f"CCKBC: {len(CCKBC)}, ISI2: {len(ISI2)}, ISI3: {len(ISI3)}, Other VIP+: {len(Other)}, VIP-: {len(VIP_Neg)}")

# %%
# Assign subtype labels to every cluster
Biases_CGE["marker_subtype"] = "Unassigned"
Biases_CGE.loc[CCKBC.index, "marker_subtype"] = "CCKBC"
Biases_CGE.loc[ISI2.index, "marker_subtype"] = "ISI2"
Biases_CGE.loc[ISI3.index, "marker_subtype"] = "ISI3"
Biases_CGE.loc[Other.index, "marker_subtype"] = "Other VIP+"
Biases_CGE.loc[VIP_Neg.index, "marker_subtype"] = "VIP-"

print("\nClassification table:")
display_cols = ["marker_subtype", "VIP_ExpL", "CCK_ExpL", "CR_ExpL", "SNCG_ExpL", "M2R_ExpL", "X22q_Bias"]
Biases_CGE[display_cols].sort_values("marker_subtype").round(3)

# %% [markdown]
# ## 5. Load Ephys Data

# %%
# Ephys cluster groupings (from TransEphys cge_cckbc_ephys_analysis)
EPHYS_GROUPS = {
    "CCKBC": [284, 289],
    "High-rate VIP": [290],
    "Irregular VIP": [291, 292],
    "Regular VIP": [279, 280],
    "LAMP5": [278, 287, 288],
}
EPHYS_GROUP_ORDER = ["CCKBC", "High-rate VIP", "Irregular VIP", "Regular VIP", "LAMP5"]
EPHYS_COLORS = {
    "CCKBC": "#d62728",
    "High-rate VIP": "#ff7f0e",
    "Irregular VIP": "#2ca02c",
    "Regular VIP": "#1f77b4",
    "LAMP5": "#9467bd",
}

# %%
# Load scANVI mapping results
mapping = pd.read_csv(MAPPING_PATH, index_col=0)
mapping.index = mapping.index.astype(int)
print(f"Mapping results: {mapping.shape[0]} cells")

# Load ephys features
ephys = pd.read_csv(EPHYS_FX_PATH, index_col=0)
ephys.index = ephys.index.astype(int)
print(f"Ephys features: {ephys.shape[0]} cells, {ephys.shape[1]} features")

# Filter to CGE interneurons
cge_mask = mapping["assigned_supercluster"] == "CGE interneuron"
cge_mapping = mapping[cge_mask]
print(f"CGE interneurons: {len(cge_mapping)}")

# Join with ephys
shared_ids = cge_mapping.index.intersection(ephys.index)
print(f"CGE cells with ephys features: {len(shared_ids)}")

df_ephys = cge_mapping.loc[shared_ids].copy()
df_ephys = df_ephys.join(ephys.loc[shared_ids])

# Add ephys group labels
cluster_to_group = {}
for group, clusters in EPHYS_GROUPS.items():
    for c in clusters:
        cluster_to_group[c] = group

df_ephys["ephys_group"] = df_ephys["assigned_cluster"].map(cluster_to_group)
df_grouped = df_ephys.dropna(subset=["ephys_group"])

print(f"\nCGE cells in defined ephys groups: {len(df_grouped)}")
print("Group distribution:")
print(df_grouped["ephys_group"].value_counts().reindex(EPHYS_GROUP_ORDER))

# %% [markdown]
# ## 6. Ephys Profile Comparison

# %%
# Key ephys features to compare
plot_features = {
    "avg_rate_hero": "Firing Rate (Hz)",
    "upstroke_downstroke_ratio_ramp": "UD Ratio",
    "rheobase_i": "Rheobase (pA)",
    "isi_cv_hero": "ISI CV",
}

fig, axes = plt.subplots(1, 4, figsize=(16, 4))
fig.patch.set_alpha(0)

for ax, (col, label) in zip(axes, plot_features.items()):
    ax.patch.set_alpha(0)
    data = []
    positions = []
    colors = []
    for i, group in enumerate(EPHYS_GROUP_ORDER):
        vals = df_grouped.loc[df_grouped["ephys_group"] == group, col].dropna().values
        if len(vals) > 0:
            data.append(vals)
            positions.append(i)
            colors.append(EPHYS_COLORS[group])

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.4),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Overlay individual points
    rng = np.random.default_rng(42)
    for i, (group, pos) in enumerate(zip(EPHYS_GROUP_ORDER, positions)):
        vals = df_grouped.loc[df_grouped["ephys_group"] == group, col].dropna().values
        if len(vals) > 0:
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(pos + jitter, vals, c=colors[i], s=10, alpha=0.5, zorder=3, edgecolors="none")

    ax.set_xticks(range(len(EPHYS_GROUP_ORDER)))
    ax.set_xticklabels([g.replace(" ", "\n") for g in EPHYS_GROUP_ORDER], fontsize=8)
    ax.set_ylabel(label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle("CGE Interneuron Subtypes -- Electrophysiology Features", fontsize=12, y=1.02)
plt.tight_layout()
plt.show()

# %%
# Mann-Whitney U tests: CCKBC vs each other group
test_features = ["avg_rate_hero", "upstroke_downstroke_ratio_ramp", "rheobase_i", "isi_cv_hero"]
cckbc_ephys = df_grouped[df_grouped["ephys_group"] == "CCKBC"]

results = []
for feat in test_features:
    cckbc_vals = cckbc_ephys[feat].dropna().values
    if len(cckbc_vals) < 2:
        continue
    for other_group in EPHYS_GROUP_ORDER[1:]:
        other_vals = df_grouped.loc[df_grouped["ephys_group"] == other_group, feat].dropna().values
        if len(other_vals) < 2:
            continue
        stat, pval = stats.mannwhitneyu(cckbc_vals, other_vals, alternative="two-sided")
        results.append({
            "Feature": feat,
            "Comparison": f"CCKBC vs {other_group}",
            "CCKBC median": np.median(cckbc_vals),
            "Other median": np.median(other_vals),
            "U": stat,
            "p-value": pval,
        })

if results:
    stats_df = pd.DataFrame(results)
    stats_df
else:
    print("Insufficient CCKBC data for statistical tests.")

# %% [markdown]
# ## 7. Convergent Classification
#
# Compare marker-based subtype assignments with ephys-based cluster groupings.

# %%
# Compute median ephys features per cluster (for clusters with ephys data)
ephys_per_cluster = df_grouped.groupby("assigned_cluster").agg(
    n_cells=("ephys_group", "size"),
    ephys_group=("ephys_group", "first"),
    firing_rate_median=("avg_rate_hero", "median"),
    isi_cv_median=("isi_cv_hero", "median"),
    ud_ratio_median=("upstroke_downstroke_ratio_ramp", "median"),
    rheobase_median=("rheobase_i", "median"),
).round(2)

# Build convergent summary table
convergent = Biases_CGE[["marker_subtype", "VIP_ExpL", "CCK_ExpL", "CR_ExpL", "SNCG_ExpL", "X22q_Bias"]].copy()
convergent = convergent.join(ephys_per_cluster[["ephys_group", "n_cells", "firing_rate_median", "isi_cv_median"]], how="left")
convergent["ephys_group"] = convergent["ephys_group"].fillna("No ephys")
convergent = convergent.round(3)

print("Convergent Classification: Marker-based vs Ephys-based")
convergent.sort_values(["marker_subtype", "ephys_group"])

# %%
# Highlight agreement/disagreement
agree = convergent[
    (convergent["marker_subtype"] == "CCKBC") & (convergent["ephys_group"] == "CCKBC")
]
disagree = convergent[
    ((convergent["marker_subtype"] == "CCKBC") & (convergent["ephys_group"] != "CCKBC"))
    | ((convergent["marker_subtype"] != "CCKBC") & (convergent["ephys_group"] == "CCKBC"))
]

print(f"\nAgreement (both marker + ephys = CCKBC): clusters {list(agree.index)}")
if len(disagree) > 0:
    print(f"Disagreement: clusters {list(disagree.index)}")
    print(disagree[["marker_subtype", "ephys_group", "VIP_ExpL", "CCK_ExpL", "CR_ExpL"]])
else:
    print("No disagreements between marker and ephys assignments.")

# %% [markdown]
# ## 8. 22q Bias Comparison Across Subtypes

# %%
def plot_subtype_bias(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias",
                      figsize=(4.8, 3.6), jitter=0.06):
    """Boxplot of bias across marker-based subtypes with statistical tests."""
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    groups = [CCKBC[label], ISI2[label], ISI3[label], Other[label], VIP_Neg[label]]
    group_labels = ["CCKBC", "ISI2", "ISI3", "Other", "VIP\u2212"]

    bp = ax.boxplot(
        groups,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.4),
        boxprops=dict(facecolor="lightgray", edgecolor="black", linewidth=1.2, alpha=0.35),
        whiskerprops=dict(color="black", linewidth=1.1),
        capprops=dict(color="black", linewidth=1.1),
    )

    for i, data in enumerate(groups, 1):
        x = np.random.normal(i, jitter, size=len(data))
        ax.scatter(x, data, s=28, color="black", alpha=0.75, edgecolors="black",
                   linewidths=0.3, zorder=3)

    ax.set_xticks(range(1, len(group_labels) + 1))
    ax.set_xticklabels(group_labels, fontsize=9)
    ax.set_ylabel(label, fontsize=9)
    ax.set_title(label, fontsize=10, pad=6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=8)

    # Statistical tests
    h_stat, p_kw = stats.kruskal(*groups)
    _, p_isi23_cck = stats.mannwhitneyu(
        list(ISI2[label]) + list(ISI3[label]), list(CCKBC[label]), alternative="greater"
    )
    _, p_isi23_vip_neg = stats.mannwhitneyu(
        list(ISI2[label]) + list(ISI3[label]), list(VIP_Neg[label]), alternative="greater"
    )

    stat_text = (
        f"Kruskal\u2013Wallis p = {p_kw:.2e}\n"
        f"ISI2+ISI3 > CCKBC: p = {p_isi23_cck:.2e}\n"
        f"ISI2+ISI3 > VIP\u2212: p = {p_isi23_vip_neg:.2e}"
    )
    ax.text(0.02, 0.98, stat_text, transform=ax.transAxes, fontsize=8, va="top")

    plt.tight_layout()
    plt.show()

    return {"KW_p": p_kw, "ISI23_vs_CCKBC_p": p_isi23_cck, "ISI23_vs_VIPneg_p": p_isi23_vip_neg}


# %%
# 22q human-gene bias
print("=== 22q Human-Gene Bias ===")
res_22q = plot_subtype_bias(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_Bias")

# %%
# 22q mouse-gene bias
print("=== 22q Mouse-Gene Bias ===")
res_22q_mouse = plot_subtype_bias(CCKBC, ISI2, ISI3, Other, VIP_Neg, label="X22q_mouse_Bias")

# %%
# Spearmans' R: VIP expression vs 22q bias (uses all 21 clusters)
rho, pval = stats.spearmanr(Biases_CGE["VIP_ExpL"], Biases_CGE["X22q_Bias"])

fig, ax = plt.subplots(figsize=(5, 4))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

colors = {"CCKBC": "#d62728", "ISI2": "#ff7f0e", "ISI3": "#2ca02c",
          "Other VIP+": "#1f77b4", "VIP-": "#9467bd"}
for subtype, color in colors.items():
    mask = Biases_CGE["marker_subtype"] == subtype
    ax.scatter(Biases_CGE.loc[mask, "VIP_ExpL"], Biases_CGE.loc[mask, "X22q_Bias"],
               c=color, label=subtype, s=40, edgecolors="black", linewidths=0.3, zorder=3)

ax.set_xlabel("VIP Expression (log2 UMI+1)", fontsize=10)
ax.set_ylabel("22q Bias", fontsize=10)
ax.set_title(f"VIP Expression vs 22q Bias\nSpearmans' R={rho:.3f}, p={pval:.3e}", fontsize=10)
ax.legend(fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 9. Multi-Disorder Comparison

# %%
# Load full bias DataFrames and extract CGE cluster values
disorder_files = {
    "ASD (w/o ID)": "ASD_HIQ_bias_addP.csv",
    "ASD (with ID)": "ASD_LIQ_bias_addP.csv",
    "SCZ": "SCZ_bias_addP.csv",
    "DDD": "DDD_61_bias_addP.csv",
    "22q": "22q_del_bias_addP.csv",
    "22q (mouse)": "22q_small_del_bias_addP.csv",
    "VNR": "UKBB_VNR_Pos_bias_addP.csv",
    "EDU": "UKBB_EDU_Pos_bias_addP.csv",
}

for name, fname in disorder_files.items():
    fpath = BIAS_DIR / fname
    if fpath.exists():
        bias_df = pd.read_csv(fpath, index_col=0)
        col = f"{name}_Bias"
        for idx in Biases_CGE.index:
            if idx in bias_df.index:
                Biases_CGE.loc[idx, col] = bias_df.loc[idx, "EFFECT"]
    else:
        print(f"Warning: {fpath} not found")

# Refresh subgroup DataFrames with new columns
CCKBC = Biases_CGE.loc[Biases_CGE["marker_subtype"] == "CCKBC"]
ISI2 = Biases_CGE.loc[Biases_CGE["marker_subtype"] == "ISI2"]
ISI3 = Biases_CGE.loc[Biases_CGE["marker_subtype"] == "ISI3"]
Other = Biases_CGE.loc[Biases_CGE["marker_subtype"] == "Other VIP+"]
VIP_Neg = Biases_CGE.loc[Biases_CGE["marker_subtype"] == "VIP-"]

# %%
# Generate boxplots for each disorder
multi_results = {}
for name in disorder_files.keys():
    col = f"{name}_Bias"
    if col in Biases_CGE.columns:
        print(f"\n=== {name} ===")
        res = plot_subtype_bias(CCKBC, ISI2, ISI3, Other, VIP_Neg, label=col)
        multi_results[name] = res

# %%
# Summary table of multi-disorder p-values
if multi_results:
    summary = pd.DataFrame(multi_results).T
    summary.columns = ["KW p-value", "ISI2+3 > CCKBC p", "ISI2+3 > VIP- p"]
    print("\nMulti-disorder statistical summary:")
    summary

# %% [markdown]
# ## 10. Mouse CCKBC Markers (SNCG, SERPINF1) in Human
#
# In mouse M1 patch-seq (Gouwens et al.), the Sncg subclass (mostly VIP- and
# strongly CCK+) likely corresponds to CCK basket cells. Among VIP+ interneurons,
# the Vip Sncg and Vip Serpinf1 types are also putative CCKBCs.
#
# **Question**: Do the human ephys-CCKBC clusters (284, 289) show higher SNCG and
# SERPINF1 expression, as the mouse analogy predicts?
#
# We check at two levels:
# 1. **Atlas level** -- cluster-average expression from the specificity matrix
# 2. **Single-cell level** -- individual patch-seq cells from the scANVI mapping

# %% [markdown]
# ### 10a. Atlas-Level SNCG & SERPINF1

# %%
# Add SNCG and SERPINF1 to the atlas-level CGE table
sncg_entrez = GeneSymbol2Entrez["SNCG"]
serpinf1_entrez = GeneSymbol2Entrez["SERPINF1"]

atlas_markers = pd.DataFrame({
    "SNCG": ExpL.loc[sncg_entrez, CGE_idx].values,
    "SERPINF1": ExpL.loc[serpinf1_entrez, CGE_idx].values,
    "VIP": ExpL.loc[GeneSymbol2Entrez["VIP"], CGE_idx].values,
    "CCK": ExpL.loc[GeneSymbol2Entrez["CCK"], CGE_idx].values,
}, index=CGE_idx)

# Map clusters to ephys groups
_c2g = {}
for _g, _cs in EPHYS_GROUPS.items():
    for _c in _cs:
        _c2g[_c] = _g
atlas_markers["ephys_group"] = [_c2g.get(x, "Other/No ephys") for x in CGE_idx]

# Display per cluster
print("Atlas-level marker expression (log2 UMI+1) per CGE cluster:")
atlas_markers.sort_values("ephys_group").round(3)

# %%
# Group summary
print("Group medians (atlas):")
atlas_group_med = atlas_markers.groupby("ephys_group")[["SNCG", "SERPINF1"]].agg(["median", "mean"]).round(3)
atlas_group_med

# %%
# Barplot comparing SNCG and SERPINF1 across ephys groups
_order = ["CCKBC", "High-rate VIP", "Irregular VIP", "Regular VIP", "LAMP5", "Other/No ephys"]
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
fig.patch.set_alpha(0)

for ax, gene in zip(axes, ["SNCG", "SERPINF1"]):
    ax.patch.set_alpha(0)
    medians = []
    group_labels = []
    for g in _order:
        sub = atlas_markers[atlas_markers["ephys_group"] == g]
        medians.append(sub[gene].median())
        group_labels.append(g)
        # Overlay individual cluster points
        jitter = np.random.normal(len(medians) - 1, 0.08, size=len(sub))
        ax.scatter(jitter, sub[gene].values, s=30, color=EPHYS_COLORS.get(g, "gray"),
                   edgecolors="black", linewidths=0.3, zorder=3, alpha=0.8)

    bars = ax.bar(range(len(medians)), medians, color="lightgray", edgecolor="black",
                  alpha=0.4, zorder=1)
    ax.set_xticks(range(len(group_labels)))
    ax.set_xticklabels([g.replace(" ", "\n") for g in group_labels], fontsize=8)
    ax.set_ylabel(f"{gene} (log2 UMI+1)")
    ax.set_title(f"{gene} by Ephys Group (Atlas)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.show()

# %% [markdown]
# ### 10b. Single-Cell Patch-Seq SNCG
#
# The scANVI h5ad contains ~3,000 HVGs including SNCG but not SERPINF1.
# We use this for single-cell SNCG analysis.

# %%
import scanpy as sc

adata = sc.read_h5ad(str(MAPPING_PATH).replace("mapping_results.csv", "query_mapped.h5ad"))

# Filter to CGE
cge_sc = adata[adata.obs["assigned_supercluster"] == "CGE interneuron"].copy()
cge_sc.obs["ephys_group"] = cge_sc.obs["assigned_cluster"].astype(str).map(
    {str(k): v for k, v in _c2g.items()}
).fillna("Other/No ephys")

# Extract SNCG counts
sncg_sc_idx = list(cge_sc.var_names).index("SNCG")
cge_sc.obs["SNCG_raw"] = np.array(cge_sc.X[:, sncg_sc_idx].todense()).flatten()
cge_sc.obs["SNCG_log1p"] = np.log1p(cge_sc.obs["SNCG_raw"])

print(f"CGE patch-seq cells: {cge_sc.shape[0]}")
print(f"SNCG in HVG: True, SERPINF1 in HVG: {'SERPINF1' in cge_sc.var_names}")

# %%
# SNCG per ephys group (single-cell)
print("Single-cell SNCG (log1p counts) by ephys group:")
sc_summary = []
for g in ["CCKBC", "High-rate VIP", "Irregular VIP", "Regular VIP", "LAMP5", "Other/No ephys"]:
    sub = cge_sc.obs[cge_sc.obs["ephys_group"] == g]
    if len(sub) > 0:
        sc_summary.append({
            "Group": g,
            "n_cells": len(sub),
            "SNCG_median": sub["SNCG_log1p"].median(),
            "SNCG_mean": sub["SNCG_log1p"].mean(),
            "pct_SNCG_positive": 100 * (sub["SNCG_raw"] > 0).mean(),
        })
sc_summary_df = pd.DataFrame(sc_summary).set_index("Group").round(2)
sc_summary_df

# %%
# Per-cluster SNCG in patch-seq cells
print("SNCG per cluster (patch-seq cells):")
sc_cluster = cge_sc.obs.groupby("assigned_cluster").agg(
    n=("SNCG_raw", "count"),
    ephys_grp=("ephys_group", "first"),
    SNCG_mean=("SNCG_log1p", "mean"),
    SNCG_median=("SNCG_log1p", "median"),
    pct_positive=("SNCG_raw", lambda x: round(100 * (x > 0).mean(), 1)),
).sort_values("SNCG_mean", ascending=False).round(2)
sc_cluster

# %%
# Boxplot: single-cell SNCG by ephys group
fig, ax = plt.subplots(figsize=(6, 4))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

_plot_order = ["CCKBC", "High-rate VIP", "Irregular VIP", "Regular VIP", "LAMP5"]
data_sc = []
colors_sc = []
for g in _plot_order:
    vals = cge_sc.obs.loc[cge_sc.obs["ephys_group"] == g, "SNCG_log1p"].values
    data_sc.append(vals)
    colors_sc.append(EPHYS_COLORS[g])

bp = ax.boxplot(data_sc, widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color="black", linewidth=1.4))
for patch, color in zip(bp["boxes"], colors_sc):
    patch.set_facecolor(color)
    patch.set_alpha(0.5)

rng = np.random.default_rng(42)
for i, (g, vals) in enumerate(zip(_plot_order, data_sc)):
    jitter = rng.uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(i + 1 + jitter, vals, c=colors_sc[i], s=12, alpha=0.5, zorder=3, edgecolors="none")

ax.set_xticks(range(1, len(_plot_order) + 1))
ax.set_xticklabels([g.replace(" ", "\n") for g in _plot_order], fontsize=8)
ax.set_ylabel("SNCG (log1p counts)")
ax.set_title("Single-Cell SNCG Expression by Ephys Group")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# MWU: CCKBC vs Irregular VIP (highest SNCG group)
if len(data_sc[0]) >= 2 and len(data_sc[2]) >= 2:
    _, p_irr = stats.mannwhitneyu(data_sc[0], data_sc[2], alternative="two-sided")
    ax.text(0.02, 0.98, f"CCKBC vs Irregular VIP: p={p_irr:.3f}\n(CCKBC n={len(data_sc[0])})",
            transform=ax.transAxes, fontsize=8, va="top")

plt.tight_layout()
plt.show()

# %% [markdown]
# ### 10c. Interpretation
#
# In mouse M1, the Sncg subclass (VIP-, CCK+) and Vip Sncg / Vip Serpinf1 types
# are putative CCK basket cells. If human ephys-CCKBC clusters (284, 289) are
# homologous, we would expect elevated SNCG and SERPINF1.
#
# **Findings:**
# - **SNCG**: Ephys-CCKBC (284, 289) have the **lowest** SNCG expression among all
#   CGE groups, at both atlas and single-cell levels. Irregular VIP (cl. 291) has
#   the highest. This is the **opposite** of the mouse prediction.
# - **SERPINF1**: At the atlas level, CCKBC clusters have the **highest** SERPINF1
#   (median=0.53 vs <0.16 for all other groups). This is driven mainly by cluster
#   284. However, overall expression is low (all values <0.65).
# - **Caveat**: Only 5 patch-seq cells in the ephys-CCKBC group (4 in cl.284, 1 in
#   cl.289), severely limiting statistical power.
#
# The marker-based CCKBC clusters (276, 279) -- defined by CR-, SNCG+, CCK+ --
# do express SNCG (by definition), but they map to "Regular VIP" or have no ephys.
# The ephys-CCKBC clusters (284, 289) are actually CR+ and SNCG-, which is why
# they were classified as ISI3 by markers. This suggests that the mouse Sncg/CCKBC
# equivalence may not directly translate to human, or that the human CCKBC identity
# is defined by different molecular markers.

# %% [markdown]
# ## 11. Summary
#
# ### Key findings
#
# **Marker-based classification** divides 21 CGE clusters into 5 subtypes:
# - **CCKBC**: CR-, SNCG+, CCK+, M2R- (from all CGE, not just VIP+)
# - **ISI2**: VIP+, CR-, M2R-, SNCG-
# - **ISI3**: VIP+, CR+, M2R-, SNCG-
# - **Other VIP+**: VIP+ clusters not fitting above
# - **VIP-**: Low VIP expression
#
# **Ephys validation**: scANVI-mapped Patch-seq data confirms that clusters 284/289
# (CCKBC by ephys) show distinct electrophysiological profiles from VIP subtypes.
#
# **Marker vs Ephys mismatch**: Ephys-CCKBC clusters (284, 289) are classified as
# ISI3 by markers (CR+, SNCG-), while marker-CCKBC clusters (276, 279) are CR-,
# SNCG+. No cluster is called CCKBC by both methods.
#
# **Mouse CCKBC markers**: SNCG is **not** elevated in human ephys-CCKBC (opposite
# of mouse). SERPINF1 shows a weak signal in cluster 284 only. The mouse Sncg
# subclass -> CCKBC mapping does not directly translate to human.
#
# **22q bias pattern**: The key reviewer question -- do CCKBC show lower 22q bias
# than ISI/VIP subtypes? -- is tested by:
# 1. Kruskal-Wallis across all 5 groups
# 2. One-sided Mann-Whitney U: ISI2+ISI3 > CCKBC
# 3. Spearmans' R: VIP expression vs 22q bias (all 21 clusters)
#
# **Caveats**: N=21 clusters total with groups of 2-8 clusters limits statistical
# power. The Spearmans' R uses all data points and may be more informative
# than group comparisons. Only 5 patch-seq cells map to ephys-CCKBC clusters.
