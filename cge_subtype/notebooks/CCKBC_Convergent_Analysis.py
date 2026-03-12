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
# # Convergent CCKBC Analysis: Cross-Species Mapping + Ephys Transfer + 22q Bias
#
# **Question**: Do imputed human CCK basket cells (CCKBCs) show different 22q11.2
# deletion bias compared to VIP+ irregular-spiking (ISI) interneurons?
#
# **Key challenge**: CCKBC is a mouse-defined cell type (Sncg subclass, VIP−, CCK+).
# No discrete CCKBC cluster exists in human transcriptomic atlases (Bakken 2021:
# AUROC=0.50 for Sncg cross-species classification; Darmanis 2015: CCK mRNA
# ubiquitous in human neurons).
#
# **Approach**: We impute human CCKBC identity using two independent modalities:
# 1. **Transcriptomic mapping** (Harmony): map mouse patch-seq CCKBCs to human CGE clusters
# 2. **Ephys signature transfer**: define mouse CCKBC ephys fingerprint, find matching human cells
#
# Then compare 22q bias between imputed CCKBCs and remaining VIP+ (ISI) clusters.

# %% [markdown]
# ## 1. Setup

# %%
# %load_ext autoreload
# %autoreload 2

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from matplotlib.lines import Line2D
from pathlib import Path

ProjDIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
CGE_DIR = ProjDIR / "cge_subtype"
sys.path.insert(1, str(ProjDIR / "src"))
from CellType_PSY import AnnotateCTDat

# Paths
BIAS_DIR = ProjDIR / "results" / "main_results" / "random" / "Centering"
HARMONY_PATH = CGE_DIR / "results" / "harmony_patchseq_mapping_results.csv"
SCVI_PATH = CGE_DIR / "results" / "patchseq_mapping_results.csv"
HUMAN_EPHYS_PATH = CGE_DIR / "dat" / "LeeDalley_ephys_fx.csv"
HUMAN_MAPPING_PATH = CGE_DIR / "results" / "human_scvi_mapping_results.csv"
MOUSE_EPHYS_PATH = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_ephys_features.csv")
MOUSE_META_PATH = Path("/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv")
ANNO_PATH = ProjDIR / "dat" / "annotation.xlsx"
OLDER_BIAS_DIR = ProjDIR / "dat" / "Spec_Bias_Jul24_Centered"

sns.set_style("ticks")
plt.rcParams["figure.dpi"] = 120

# %% [markdown]
# ## 2. Load Annotation & 22q-Related Bias Sets

# %%
# Cell type annotation (462 types)
anno = pd.read_excel(ANNO_PATH, index_col=0)
anno = anno[anno.index.notna()]  # drop NaN rows
anno.index = anno.index.astype(int)
CGE_idx = anno[anno["Supercluster"] == "CGE interneuron"].index.values

# Identify VIP+ (ISI) clusters within CGE
VIP_idx = anno[(anno["Supercluster"] == "CGE interneuron") &
               (anno["Subtype auto-annotation"] == "INT-VIP")].index.values
print(f"CGE clusters: {len(CGE_idx)} (indices: {sorted(CGE_idx)})")
print(f"  INT-VIP: {len(VIP_idx)} — {sorted(VIP_idx)}")

# %%
# Load 22q-related bias gene sets
bias_sets = {}

# 22q deletion (main)
bias_sets["22q_del"] = pd.read_csv(BIAS_DIR / "22q_del_bias_addP.csv", index_col=0)
# 22q mouse model deletion
bias_sets["22q_mouse"] = pd.read_csv(BIAS_DIR / "22q_small_del_bias_addP.csv", index_col=0)
# 22q highly expressed genes
bias_sets["22q_HighExp"] = pd.read_csv(OLDER_BIAS_DIR / "HCT.22q.HighExp.csv", index_col=0)
# 22q DEGs (day 75 — representative timepoint)
bias_sets["22q_DEG_d75"] = pd.read_csv(OLDER_BIAS_DIR / "DEG_22q_day75_Bias.csv", index_col=0)

for name, df in bias_sets.items():
    print(f"  {name}: {df.shape}")

# Extract CGE bias for each gene set
cge_bias = pd.DataFrame(index=CGE_idx)
for name, df in bias_sets.items():
    cge_bias[f"bias_{name}"] = df.loc[CGE_idx, "EFFECT"]
cge_bias["cluster_name"] = anno.loc[CGE_idx, "Subtype auto-annotation"]
print(f"\nCGE bias: {cge_bias.shape}")
cge_bias.head()

# %% [markdown]
# ## 3. Cross-Species Transcriptomic Mapping (Harmony)
#
# Map mouse patch-seq CCKBCs → human CGE clusters using Harmony-corrected
# PCA space with k-NN (k=30) distance-weighted majority vote.

# %%
harmony = pd.read_csv(HARMONY_PATH, index_col=0)
harmony_cge = harmony[harmony["assigned_supercluster"] == "CGE interneuron"]

print(f"Total mapped cells: {len(harmony)}")
print(f"  → CGE: {len(harmony_cge)}")
print(f"  → CCKBCs total: {harmony['is_cckbc'].sum()}")
print(f"  → CCKBCs in CGE: {harmony_cge['is_cckbc'].sum()}")

# Per-cluster CCKBC fraction (Harmony)
harm_stats = []
for cl in sorted(harmony_cge["assigned_cluster"].unique()):
    cl_cells = harmony_cge[harmony_cge["assigned_cluster"] == cl]
    n_total = len(cl_cells)
    n_cckbc = cl_cells["is_cckbc"].sum()
    harm_stats.append({
        "cluster_id": int(cl), "n_harmony": n_total, "n_cckbc_harmony": n_cckbc,
        "cckbc_frac_harmony": n_cckbc / n_total if n_total > 0 else 0,
        "mean_conf_harmony": cl_cells["conf_cluster"].mean(),
    })
harm_df = pd.DataFrame(harm_stats).set_index("cluster_id")

# Also load scVI for comparison
scvi = pd.read_csv(SCVI_PATH, index_col=0)
scvi_cge = scvi[scvi["assigned_human_supercluster"] == "CGE interneuron"]
scvi_stats = []
for cl in sorted(scvi_cge["assigned_human_cluster"].unique()):
    cl_cells = scvi_cge[scvi_cge["assigned_human_cluster"] == cl]
    n_total = len(cl_cells)
    n_cckbc = cl_cells["is_cckbc"].sum()
    scvi_stats.append({
        "cluster_id": int(cl), "n_scvi": n_total, "n_cckbc_scvi": n_cckbc,
        "cckbc_frac_scvi": n_cckbc / n_total if n_total > 0 else 0,
    })
scvi_df = pd.DataFrame(scvi_stats).set_index("cluster_id")

# Combine
mapping_df = harm_df.join(scvi_df, how="outer").fillna(0)
print(f"\n=== Harmony vs scVI CCKBC fraction (CGE clusters) ===")
print(mapping_df[["cckbc_frac_harmony", "cckbc_frac_scvi", "n_harmony", "n_scvi"]]
      .sort_values("cckbc_frac_harmony", ascending=False)
      .head(15).to_string())

# %%
# Concordance between Harmony and scVI
both = mapping_df[(mapping_df["n_harmony"] > 0) & (mapping_df["n_scvi"] > 0)]
rho, p = spearmanr(both["cckbc_frac_harmony"], both["cckbc_frac_scvi"])
print(f"Harmony-scVI concordance: Spearman rho={rho:.3f}, p={p:.2e}, N={len(both)} clusters")

fig, ax = plt.subplots(figsize=(5, 5))
fig.patch.set_alpha(0); ax.patch.set_alpha(0)
ax.scatter(both["cckbc_frac_scvi"], both["cckbc_frac_harmony"], s=40, alpha=0.7)
for cl in both.index:
    if both.loc[cl, "cckbc_frac_harmony"] > 0.3 or both.loc[cl, "cckbc_frac_scvi"] > 0.3:
        ax.annotate(str(cl), (both.loc[cl, "cckbc_frac_scvi"], both.loc[cl, "cckbc_frac_harmony"]),
                    fontsize=8, ha="left")
ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
ax.set_xlabel("CCKBC fraction (scVI)")
ax.set_ylabel("CCKBC fraction (Harmony)")
ax.set_title(f"Cross-species mapping concordance\nSpearman rho={rho:.2f}, p={p:.2e}")
plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_harmony_scvi_concordance.pdf"), transparent=True)
plt.show()

# %% [markdown]
# ## 4. Ephys Signature Transfer
#
# Define mouse CCKBC electrophysiology fingerprint relative to other VIP cells
# (within-mouse z-score), then score human CGE cells for CCKBC-likeness.

# %%
# Load mouse M1 ephys + metadata
m_ephys = pd.read_csv(MOUSE_EPHYS_PATH).rename(columns={"cell id": "Cell"})
m_meta = pd.read_csv(MOUSE_META_PATH, sep="\t")
m_merged = m_ephys.merge(m_meta[["Cell", "RNA family", "RNA type"]], on="Cell")

# Define CCKBC (Sncg family + Vip Sncg/Serpinf1)
cckbc_types = ["Vip Sncg", "Vip Serpinf1_1", "Vip Serpinf1_2",
               "Sncg Npy2r", "Sncg Col14a1", "Sncg Calb1_1", "Sncg Calb1_2"]
vip_other_types = ["Vip Mybpc1_1", "Vip Mybpc1_2", "Vip Mybpc1_3", "Vip Gpc3",
                   "Vip Chat_1", "Vip C1ql1", "Vip Htr1f", "Vip Col15a1"]
m_merged["is_cckbc"] = m_merged["RNA type"].isin(cckbc_types)
m_merged["is_vip_other"] = m_merged["RNA type"].isin(vip_other_types)
mouse_cge = m_merged[m_merged["is_cckbc"] | m_merged["is_vip_other"]].copy()

print(f"Mouse CGE pool: {len(mouse_cge)} (CCKBC={mouse_cge.is_cckbc.sum()}, VIP-other={mouse_cge.is_vip_other.sum()})")

# %%
# Shared feature mapping (mouse → human, conceptually equivalent)
feature_map = {
    "AP width": "width_ramp",
    "AP threshold": "threshold_v_ramp",
    "ISI CV": "isi_cv_hero",
    "R_input": "input_resistance",
    "tau": "tau",
}

# Compute CCKBC direction vector (within-mouse z-scores)
cckbc_directions = {}
print("Mouse CCKBC ephys signature (within-mouse z-scores):")
print(f"  {'Feature':18s} {'CCKBC_z':>8} {'VIP_z':>8} {'Direction':>10}")
print("  " + "-" * 50)
for m_feat in feature_map:
    vals = mouse_cge[m_feat].dropna()
    mu_all, sd_all = vals.mean(), vals.std()
    if sd_all == 0:
        continue
    cckbc_z = (mouse_cge.loc[mouse_cge["is_cckbc"], m_feat].dropna().mean() - mu_all) / sd_all
    vip_z = (mouse_cge.loc[mouse_cge["is_vip_other"], m_feat].dropna().mean() - mu_all) / sd_all
    cckbc_directions[m_feat] = cckbc_z
    direction = "higher" if cckbc_z > vip_z else "lower"
    print(f"  {m_feat:18s} {cckbc_z:>+8.3f} {vip_z:>+8.3f} {'CCKBC ' + direction:>10}")

# %%
# Load human ephys + cluster mapping
h_ephys = pd.read_csv(HUMAN_EPHYS_PATH)
h_map = pd.read_csv(HUMAN_MAPPING_PATH, index_col=0)
h_merged = h_ephys.set_index("specimen_id").join(h_map, how="inner")
h_cge = h_merged[h_merged["assigned_supercluster"] == "CGE interneuron"].copy()
print(f"Human CGE cells with ephys: {len(h_cge)}")

# Z-score human features within human CGE pool
human_z = pd.DataFrame(index=h_cge.index)
for m_feat, h_feat in feature_map.items():
    if m_feat not in cckbc_directions:
        continue
    vals = pd.to_numeric(h_cge[h_feat], errors="coerce")
    mu, sd = vals.mean(), vals.std()
    if sd > 0:
        human_z[m_feat] = (vals - mu) / sd

# Project onto CCKBC direction (higher score = more CCKBC-like)
direction_vec = np.array([cckbc_directions[f] for f in human_z.columns])
direction_norm = direction_vec / np.linalg.norm(direction_vec)
h_cge["cckbc_ephys_score"] = human_z.values @ direction_norm

# Per-cluster ephys CCKBC score
ephys_stats = h_cge.groupby("assigned_cluster").agg(
    n_ephys=("cckbc_ephys_score", "count"),
    ephys_score=("cckbc_ephys_score", "mean"),
).rename_axis("cluster_id")
ephys_stats.index = ephys_stats.index.astype(int)

# %%
# Correlation: ephys score vs Harmony CCKBC fraction
combined = mapping_df.join(ephys_stats, how="inner")
valid = combined[combined["n_ephys"] >= 3]
rho_eh, p_eh = spearmanr(valid["ephys_score"], valid["cckbc_frac_harmony"])
print(f"Ephys score vs Harmony CCKBC frac: Spearman rho={rho_eh:.3f}, p={p_eh:.4f}, N={len(valid)} clusters")

fig, ax = plt.subplots(figsize=(5, 5))
fig.patch.set_alpha(0); ax.patch.set_alpha(0)
ax.scatter(valid["cckbc_frac_harmony"], valid["ephys_score"], s=40, alpha=0.7)
for cl in valid.index:
    ax.annotate(str(cl), (valid.loc[cl, "cckbc_frac_harmony"], valid.loc[cl, "ephys_score"]),
                fontsize=8, ha="left")
ax.set_xlabel("CCKBC fraction (Harmony transcriptomic mapping)")
ax.set_ylabel("CCKBC ephys score (signature transfer)")
ax.set_title(f"Transcriptomic-ephys convergence\nSpearman rho={rho_eh:.2f}, p={p_eh:.3f}")
plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_ephys_harmony_convergence.pdf"), transparent=True)
plt.show()

# %% [markdown]
# ## 5. Define Imputed Human CCKBC vs VIP+ (ISI) Clusters
#
# Use Harmony CCKBC fraction as primary criterion (transcriptomic mapping),
# supported by ephys concordance. Threshold: CCKBC fraction > 0.5.
#
# Compare against remaining INT-VIP clusters (the VIP+/ISI type), excluding
# LAMP5 and unannotated clusters — since we already show non-VIP CGE has
# lower 22q bias independently.

# %%
# Merge all data into master CGE table
master = pd.DataFrame(index=CGE_idx)
master["cluster_name"] = anno.loc[CGE_idx, "Subtype auto-annotation"]

# Harmony mapping
for col in ["n_harmony", "n_cckbc_harmony", "cckbc_frac_harmony", "mean_conf_harmony"]:
    master[col] = mapping_df[col].reindex(CGE_idx).fillna(0)

# scVI mapping
for col in ["n_scvi", "n_cckbc_scvi", "cckbc_frac_scvi"]:
    master[col] = mapping_df[col].reindex(CGE_idx).fillna(0)

# Ephys score
master["ephys_score"] = ephys_stats["ephys_score"].reindex(CGE_idx)
master["n_ephys"] = ephys_stats["n_ephys"].reindex(CGE_idx).fillna(0).astype(int)

# All 22q-related bias sets
for name in bias_sets:
    master[f"bias_{name}"] = cge_bias[f"bias_{name}"]

# Classify: CCKBC if Harmony fraction > 0.5
master["is_imputed_cckbc"] = master["cckbc_frac_harmony"] > 0.5
cckbc_clusters = master[master["is_imputed_cckbc"]].index.tolist()

# VIP+ (ISI) = INT-VIP annotated clusters that are NOT imputed CCKBC
vip_isi_clusters = [c for c in VIP_idx if c not in cckbc_clusters]

print(f"Imputed CCKBC clusters (Harmony frac > 0.5): {sorted(cckbc_clusters)}")
print(f"VIP+ (ISI) clusters: {sorted(vip_isi_clusters)}")
print(f"\nImputed CCKBC details:")
print(master.loc[cckbc_clusters, ["cluster_name", "cckbc_frac_harmony", "cckbc_frac_scvi",
                                    "ephys_score", "n_harmony"]].to_string())
print(f"\nVIP+ (ISI) details:")
print(master.loc[vip_isi_clusters, ["cluster_name", "cckbc_frac_harmony", "cckbc_frac_scvi",
                                      "ephys_score", "n_harmony"]].to_string())

# %% [markdown]
# ## 6. 22q-Related Bias Comparison: Imputed CCKBC vs VIP+ (ISI)

# %%
# Compare bias across all 22q-related gene sets: CCKBC vs VIP+ (ISI)
print("=== CCKBC vs VIP+ (ISI): 22q-Related Bias ===\n")
print(f"  CCKBC clusters (n={len(cckbc_clusters)}): {sorted(cckbc_clusters)}")
print(f"  VIP+ ISI clusters (n={len(vip_isi_clusters)}): {sorted(vip_isi_clusters)}")

gene_set_labels = {
    "22q_del": "22q Deletion",
    "22q_mouse": "22q Mouse Model",
    "22q_HighExp": "22q High Expressed",
    "22q_DEG_d75": "22q DEG (day 75)",
}

results = []
for name, label in gene_set_labels.items():
    col = f"bias_{name}"
    cckbc_vals = master.loc[cckbc_clusters, col].dropna()
    vip_vals = master.loc[vip_isi_clusters, col].dropna()

    # One-tailed: test if CCKBC < VIP+ (ISI), consistent with mouse prediction
    u_stat, mwu_p = mannwhitneyu(cckbc_vals, vip_vals, alternative="less")
    diff = cckbc_vals.mean() - vip_vals.mean()
    results.append({
        "gene_set": label,
        "CCKBC_mean": cckbc_vals.mean(),
        "VIP_mean": vip_vals.mean(),
        "diff": diff,
        "direction": "CCKBC higher" if diff > 0 else "CCKBC lower",
        "MWU_p": mwu_p,
    })

results_df = pd.DataFrame(results)
print("\n" + results_df.to_string(index=False))

# %%
# Box plots: CCKBC vs VIP+ (ISI) across 22q-related gene sets
fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
fig.patch.set_alpha(0)

for idx, (name, label) in enumerate(gene_set_labels.items()):
    ax = axes[idx]
    ax.patch.set_alpha(0)

    col = f"bias_{name}"
    data_cckbc = master.loc[cckbc_clusters, col].dropna()
    data_vip = master.loc[vip_isi_clusters, col].dropna()

    plot_data = pd.DataFrame({
        "Bias": pd.concat([data_cckbc, data_vip]),
        "Group": (["Imputed\nCCKBC"] * len(data_cckbc) +
                  ["VIP+ (ISI)"] * len(data_vip)),
    })

    sns.boxplot(data=plot_data, x="Group", y="Bias", ax=ax,
                palette=["#E74C3C", "#3498DB"], width=0.5, fliersize=3)
    sns.stripplot(data=plot_data, x="Group", y="Bias", ax=ax,
                  color="black", size=4, alpha=0.6)

    row = results_df[results_df["gene_set"] == label].iloc[0]
    p_str = f"p={row['MWU_p']:.3f}" if row["MWU_p"] >= 0.001 else f"p={row['MWU_p']:.1e}"
    ax.set_title(f"{label}\n{row['direction']}, {p_str}", fontsize=10)
    ax.set_ylabel("Bias (EFFECT)" if idx == 0 else "")
    ax.set_xlabel("")

plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_cckbc_vs_vip_22q_bias.pdf"), transparent=True)
plt.show()

# %% [markdown]
# ## 6b. Three-Way Comparison: VIP- CCKBC vs VIP+ CCKBC vs VIP+ (ISI)
#
# The 5 imputed CCKBC clusters split by VIP expression:
# - **VIP- CCKBC** (277, 278): VIP-negative, low 22q bias
# - **VIP+ CCKBC** (279, 280, 281): VIP-positive, high 22q bias
#
# Do VIP+ CCKBCs differ from VIP+ ISI, or is VIP status the real predictor?

# %%
# Define the three groups
vip_neg_cckbc = [c for c in cckbc_clusters if c not in VIP_idx]
vip_pos_cckbc = [c for c in cckbc_clusters if c in VIP_idx]

print(f"VIP- CCKBC: {sorted(vip_neg_cckbc)} (n={len(vip_neg_cckbc)})")
print(f"VIP+ CCKBC: {sorted(vip_pos_cckbc)} (n={len(vip_pos_cckbc)})")
print(f"VIP+ (ISI): {sorted(vip_isi_clusters)} (n={len(vip_isi_clusters)})")

# %%
# Load expression data for marker comparison
ExpL = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv",
                   index_col=0)
ExpL.columns = [int(x) for x in ExpL.columns.values]
ExpL = ExpL.map(lambda x: np.log2(x + 1))

sys.path.insert(1, str(ProjDIR / "src"))
from CellType_PSY import LoadGeneINFO
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

marker_genes = {
    "VIP": "VIP", "CCK": "CCK", "CNR1": "CNR1", "CALB2": "CALB2",
    "SNCG": "SNCG", "CHRM2": "CHRM2", "TAC3": "TAC3", "HTR3A": "HTR3A",
    "RELN": "RELN", "SP8": "SP8", "PROX1": "PROX1", "ADARB2": "ADARB2",
}

groups = {"VIP- CCKBC": vip_neg_cckbc, "VIP+ CCKBC": vip_pos_cckbc, "VIP+ (ISI)": vip_isi_clusters}

# Marker expression table
marker_rows = []
for name, symbol in marker_genes.items():
    eid = GeneSymbol2Entrez.get(symbol)
    if eid is None or eid not in ExpL.index:
        continue
    row = {"Marker": name}
    for grp_name, grp_idx in groups.items():
        row[grp_name] = np.mean([ExpL.loc[eid, c] for c in grp_idx])
    marker_rows.append(row)
marker_df = pd.DataFrame(marker_rows).set_index("Marker")
print("=== Marker Expression (log2) ===")
print(marker_df.to_string(float_format="%.2f"))

# %%
# Ephys comparison (human)
ephys_features = ["input_resistance", "tau", "isi_cv_hero", "sag",
                  "fi_fit_slope", "upstroke_downstroke_ratio_ramp", "avg_rate_hero"]

ephys_rows = []
for feat in ephys_features:
    vals = pd.to_numeric(h_cge[feat], errors="coerce")
    row = {"Feature": feat}
    for grp_name, grp_idx in groups.items():
        grp_vals = vals[h_cge["assigned_cluster"].isin(grp_idx)].dropna()
        row[grp_name] = grp_vals.mean() if len(grp_vals) > 0 else np.nan
        row[f"n_{grp_name}"] = len(grp_vals)
    ephys_rows.append(row)
ephys_df = pd.DataFrame(ephys_rows).set_index("Feature")
print("\n=== Human Ephys ===")
print(ephys_df[[c for c in ephys_df.columns if not c.startswith("n_")]].to_string(float_format="%.3f"))
print(f"\nN cells: VIP-CCKBC={ephys_df['n_VIP- CCKBC'].iloc[0]:.0f}, "
      f"VIP+CCKBC={ephys_df['n_VIP+ CCKBC'].iloc[0]:.0f}, "
      f"VIP+(ISI)={ephys_df['n_VIP+ (ISI)'].iloc[0]:.0f}")

# %%
# Mouse ground-truth ephys: CCKBC vs VIP-other firing properties
mouse_cge_cckbc = mouse_cge[mouse_cge["is_cckbc"]]
mouse_cge_vip = mouse_cge[mouse_cge["is_vip_other"]]

mouse_feats = ["AP count", "AP count 2nd half", "ISI CV", "AP width", "R_input", "tau"]
print("\n=== Mouse M1 Patch-Seq Ephys (ground-truth labels) ===")
print(f"  {'Feature':>20}  {'CCKBC':>8}  {'VIP-other':>10}  {'MWU_p':>8}")
print("  " + "-" * 55)
for feat in mouse_feats:
    cv = pd.to_numeric(mouse_cge_cckbc[feat], errors="coerce").dropna()
    vv = pd.to_numeric(mouse_cge_vip[feat], errors="coerce").dropna()
    if len(cv) >= 2 and len(vv) >= 2:
        _, p = mannwhitneyu(cv, vv, alternative="two-sided")
        print(f"  {feat:>20}  {cv.mean():>8.3f}  {vv.mean():>10.3f}  {p:>8.4f}")

# %%
# Three-way bias comparison
print("\n=== 22q Bias: Three-Way Comparison ===")
print(f"  {'Gene Set':>20}  {'VIP-CCKBC':>10}  {'VIP+CCKBC':>10}  {'VIP+(ISI)':>10}")
print("  " + "-" * 55)
for name, label in gene_set_labels.items():
    col = f"bias_{name}"
    vals = {}
    for grp_name, grp_idx in groups.items():
        vals[grp_name] = master.loc[grp_idx, col].mean()
    print(f"  {label:>20}  {vals['VIP- CCKBC']:>10.4f}  {vals['VIP+ CCKBC']:>10.4f}  {vals['VIP+ (ISI)']:>10.4f}")

# %%
# Box plot: three-way comparison for 22q deletion bias
fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
fig.patch.set_alpha(0)

group_colors = {"VIP-\nCCKBC": "#2ECC71", "VIP+\nCCKBC": "#E74C3C", "VIP+\n(ISI)": "#3498DB"}

for idx, (name, label) in enumerate(gene_set_labels.items()):
    ax = axes[idx]
    ax.patch.set_alpha(0)
    col = f"bias_{name}"

    plot_parts = []
    for grp_label, grp_idx, short in [("VIP-\nCCKBC", vip_neg_cckbc, "VIP-CCKBC"),
                                       ("VIP+\nCCKBC", vip_pos_cckbc, "VIP+CCKBC"),
                                       ("VIP+\n(ISI)", vip_isi_clusters, "VIP+ISI")]:
        vals = master.loc[grp_idx, col].dropna()
        plot_parts.append(pd.DataFrame({"Bias": vals, "Group": grp_label}))
    plot_data = pd.concat(plot_parts)

    sns.boxplot(data=plot_data, x="Group", y="Bias", ax=ax,
                palette=list(group_colors.values()), width=0.5, fliersize=3,
                order=list(group_colors.keys()))
    sns.stripplot(data=plot_data, x="Group", y="Bias", ax=ax,
                  color="black", size=4, alpha=0.6, order=list(group_colors.keys()))
    ax.set_title(label, fontsize=10)
    ax.set_ylabel("Bias (EFFECT)" if idx == 0 else "")
    ax.set_xlabel("")

plt.suptitle("22q Bias: VIP- CCKBC vs VIP+ CCKBC vs VIP+ (ISI)", fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_threeway_22q_bias.pdf"), transparent=True)
plt.show()

# %% [markdown]
# ## 6c. Cross-Species Ephys Alignment: Human Cells in Mouse Ephys Space
#
# Validate human cluster identity by embedding human CGE cells alongside
# mouse CCKBC and VIP-other cells in shared electrophysiology feature space.
# Uses 16 matched feature pairs for PCA/UMAP (all available shared features),
# with within-species z-scoring to remove cross-species baseline differences.

# %%
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
import umap

# All 16 matched mouse → human feature pairs
full_feature_map = {
    "AP count":             "avg_rate_hero",              # firing rate
    "AP count 2nd half":    "adapt_hero",                 # sustained firing
    "ISI CV":               "isi_cv_hero",                # firing regularity
    "AP average amp adapt": "upstroke_adapt_ratio",       # waveform adaptation
    "latency":              "latency_hero",               # first-spike latency
    "R_input":              "input_resistance",            # input resistance
    "tau":                  "tau",                         # membrane time constant
    "AP width":             "width_ramp",                 # spike width
    "AP threshold":         "threshold_v_ramp",           # spike threshold
    "AP amplitude":         "peak_deltav_ramp",           # spike amplitude
    "AHP":                  "fast_trough_deltav_ramp",    # afterhyperpolarization
    "AP amp adapt":         "peak_v_adapt_ratio",         # amplitude adaptation
    "3rd AP width":         "width_hero",                 # sustained spike width
    "3rd AP threshold":     "threshold_v_hero",           # sustained threshold
    "3rd AP amplitude":     "peak_deltav_hero",           # sustained amplitude
    "3rd AHP":              "fast_trough_deltav_hero",    # sustained AHP
}

# Discriminative subset (p < 0.01 in mouse CCKBC vs VIP-other) — used for embedding
discrim_feature_map = {
    "AP count":             "avg_rate_hero",              # p=0.0001
    "AP count 2nd half":    "adapt_hero",                 # p=0.0001
    "ISI CV":               "isi_cv_hero",                # p=0.0002
    "AP average amp adapt": "upstroke_adapt_ratio",       # p<0.0001
    "latency":              "latency_hero",               # p=0.008
    "AP count 1st half":    "first_isi_inv_hero",         # p=0.0014
}

# Display name lookup
feature_display = {v: k for k, v in full_feature_map.items()}
feature_display.update({"first_isi_inv_hero": "AP count 1st half"})

# Define human cell groups for cross-species alignment
def classify_human_group(cluster_id):
    if cluster_id in vip_neg_cckbc:
        return "Human VIP- CCKBC"
    elif cluster_id in vip_pos_cckbc:
        return "Human VIP+ CCKBC"
    elif cluster_id in vip_isi_clusters:
        return "Human VIP+ ISI"
    return None

# --- Prepare data using ALL 16 features (for violin plots) ---
all_mouse_cols = list(full_feature_map.keys())
all_human_cols = list(full_feature_map.values())

mouse_feat_all = mouse_cge[all_mouse_cols].copy()
mouse_feat_all.columns = all_human_cols
mouse_feat_all.index = mouse_cge["Cell"].values
mouse_labels_all = pd.Series(
    ["Mouse CCKBC" if x else "Mouse VIP-other" for x in mouse_cge["is_cckbc"].values],
    index=mouse_cge["Cell"].values,
)

h_cge["human_group"] = h_cge["assigned_cluster"].map(classify_human_group)
h_of_interest = h_cge[h_cge["human_group"].notna()].copy()
human_feat_all = h_of_interest[all_human_cols].copy()
human_labels_all = h_of_interest["human_group"]

# --- Prepare data using discriminative features only (for embedding/kNN/classifier) ---
discrim_mouse_cols = list(discrim_feature_map.keys())
discrim_human_cols = list(discrim_feature_map.values())

mouse_feat_d = mouse_cge[discrim_mouse_cols].copy()
mouse_feat_d.columns = discrim_human_cols
mouse_feat_d.index = mouse_cge["Cell"].values

human_feat_d = h_of_interest[discrim_human_cols].copy()

# Drop NaN — discriminative features for embedding
mouse_feat_d_clean = mouse_feat_d.dropna()
human_feat_d_clean = human_feat_d.dropna()
mouse_labels_clean = mouse_labels_all.loc[mouse_feat_d_clean.index]
human_labels_clean = human_labels_all.loc[human_feat_d_clean.index]

print(f"Embedding features (discriminative, p<0.01): {len(discrim_feature_map)}")
print(f"Violin plot features (all matched pairs): {len(full_feature_map)}")
print(f"\nMouse: {len(mouse_feat_d_clean)} cells (CCKBC={(mouse_labels_clean=='Mouse CCKBC').sum()}, VIP-other={(mouse_labels_clean=='Mouse VIP-other').sum()})")
print(f"Human: {len(human_feat_d_clean)} cells")
for g in ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]:
    print(f"  {g}: {(human_labels_clean == g).sum()}")

# Z-score within species (discriminative features)
scaler_m = StandardScaler()
scaler_h = StandardScaler()
mouse_z = pd.DataFrame(scaler_m.fit_transform(mouse_feat_d_clean),
                        index=mouse_feat_d_clean.index, columns=mouse_feat_d_clean.columns)
human_z_cs = pd.DataFrame(scaler_h.fit_transform(human_feat_d_clean),
                           index=human_feat_d_clean.index, columns=human_feat_d_clean.columns)

# Combine for PCA/UMAP
combined_z = pd.concat([mouse_z, human_z_cs])
combined_labels = pd.concat([mouse_labels_clean, human_labels_clean])

# %%
# PCA + UMAP on discriminative features
pca_cs = PCA(n_components=min(6, combined_z.shape[1]))
pca_coords = pca_cs.fit_transform(combined_z)

reducer = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42)
umap_coords = reducer.fit_transform(combined_z)

print(f"PCA explained variance: {pca_cs.explained_variance_ratio_.round(3)}")
print(f"Cumulative: {np.cumsum(pca_cs.explained_variance_ratio_).round(3)}")
print(f"\nPC1 loadings:")
for j, feat in enumerate(combined_z.columns):
    print(f"  {feature_display.get(feat, feat):>25}: {pca_cs.components_[0, j]:+.3f}")

# %%
# Figure: PCA + UMAP side by side
GROUP_COLORS = {
    "Mouse CCKBC": "#1f77b4", "Mouse VIP-other": "#ff7f0e",
    "Human VIP- CCKBC": "#2ca02c", "Human VIP+ CCKBC": "#d62728", "Human VIP+ ISI": "#9467bd",
}
GROUP_MARKERS = {
    "Mouse CCKBC": "o", "Mouse VIP-other": "s",
    "Human VIP- CCKBC": "^", "Human VIP+ CCKBC": "D", "Human VIP+ ISI": "v",
}
PLOT_ORDER = ["Mouse VIP-other", "Mouse CCKBC", "Human VIP+ ISI", "Human VIP- CCKBC", "Human VIP+ CCKBC"]

fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
fig.patch.set_alpha(0)

n_d = len(discrim_feature_map)
for ax, coords, title, xlab, ylab in [
    (axes[0], pca_coords[:, :2], f"PCA: Cross-species Ephys ({n_d} discriminative features)",
     f"PC1 ({pca_cs.explained_variance_ratio_[0]*100:.1f}%)",
     f"PC2 ({pca_cs.explained_variance_ratio_[1]*100:.1f}%)"),
    (axes[1], umap_coords, f"UMAP: Cross-species Ephys ({n_d} discriminative features)", "UMAP1", "UMAP2"),
]:
    ax.patch.set_alpha(0)
    for grp in PLOT_ORDER:
        mask = combined_labels.values == grp
        if mask.sum() == 0:
            continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=GROUP_COLORS[grp], marker=GROUP_MARKERS[grp],
                   s=40 if "Mouse" in grp else 60,
                   alpha=0.5 if "Mouse" in grp else 0.85,
                   edgecolors="white", linewidths=0.3,
                   label=f"{grp} (n={mask.sum()})",
                   zorder=2 if "Mouse" in grp else 3)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab); ax.set_title(title, fontweight="bold")

handles = [Line2D([0], [0], marker=GROUP_MARKERS[g], color="w",
                   markerfacecolor=GROUP_COLORS[g], markersize=9, label=g,
                   markeredgecolor="grey", linewidth=0) for g in PLOT_ORDER]
fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=10,
           frameon=True, framealpha=0.8)
plt.tight_layout(rect=[0, 0.12, 1, 1])
plt.savefig(str(CGE_DIR / "results" / "fig_cross_species_ephys_pca_umap.pdf"), transparent=True, dpi=150)
plt.show()

# %%
# k-NN validation: for each human cell, what fraction of k=5 nearest mouse neighbors are CCKBC?
K = 5
nn = NearestNeighbors(n_neighbors=K, metric="euclidean")
nn.fit(mouse_z.values)

distances, indices = nn.kneighbors(human_z_cs.values)
mouse_labels_arr = mouse_labels_clean.values

cckbc_fracs = [np.mean(mouse_labels_arr[indices[i]] == "Mouse CCKBC") for i in range(len(human_z_cs))]
nn_results = pd.DataFrame({
    "group": human_labels_clean.values,
    "frac_cckbc_neighbors": cckbc_fracs,
    "mean_neighbor_dist": distances.mean(axis=1),
}, index=human_z_cs.index)

print("Fraction of nearest mouse neighbors that are CCKBC:")
human_grp_order = ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]
for grp in human_grp_order:
    sub = nn_results[nn_results["group"] == grp]
    if len(sub) == 0:
        continue
    print(f"  {grp} (n={len(sub)}): mean={sub['frac_cckbc_neighbors'].mean():.3f}, "
          f">50% CCKBC: {(sub['frac_cckbc_neighbors'] > 0.5).sum()}/{len(sub)}")

# Statistical tests
vip_neg_nn = nn_results.loc[nn_results["group"] == "Human VIP- CCKBC", "frac_cckbc_neighbors"]
vip_pos_cckbc_nn = nn_results.loc[nn_results["group"] == "Human VIP+ CCKBC", "frac_cckbc_neighbors"]
vip_isi_nn = nn_results.loc[nn_results["group"] == "Human VIP+ ISI", "frac_cckbc_neighbors"]

all_cckbc_nn = pd.concat([vip_neg_nn, vip_pos_cckbc_nn])
if len(all_cckbc_nn) > 0 and len(vip_isi_nn) > 0:
    u, p_mwu = mannwhitneyu(all_cckbc_nn, vip_isi_nn, alternative="two-sided")
    print(f"\nAll CCKBC vs VIP+ ISI: MWU p={p_mwu:.4e}")

if len(vip_neg_nn) > 0 and len(vip_pos_cckbc_nn) > 0 and len(vip_isi_nn) > 0:
    h_kw, p_kw = kruskal(vip_neg_nn, vip_pos_cckbc_nn, vip_isi_nn)
    print(f"Kruskal-Wallis (3-group): H={h_kw:.2f}, p={p_kw:.4e}")

# %%
# Figure: k-NN CCKBC neighbor fraction box plot
fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_alpha(0); ax.patch.set_alpha(0)

group_data = [nn_results.loc[nn_results["group"] == g, "frac_cckbc_neighbors"].values
              for g in human_grp_order]
nn_colors = [GROUP_COLORS[g] for g in human_grp_order]

bp = ax.boxplot(group_data, patch_artist=True, widths=0.5,
                medianprops=dict(color="black", linewidth=2))
for patch, color in zip(bp["boxes"], nn_colors):
    patch.set_facecolor(color); patch.set_alpha(0.6)

rng = np.random.default_rng(42)
for i, (data, color) in enumerate(zip(group_data, nn_colors)):
    jitter = rng.uniform(-0.12, 0.12, size=len(data))
    ax.scatter(np.full(len(data), i + 1) + jitter, data,
               c=color, alpha=0.5, s=20, edgecolors="white", linewidths=0.3, zorder=3)

ax.set_xticklabels([f"{g.replace('Human ', '')}\n(n={len(d)})"
                     for g, d in zip(human_grp_order, group_data)], fontsize=10)
ax.set_ylabel("Fraction CCKBC among k=5 nearest mouse neighbors", fontsize=11)
ax.set_title("Cross-Species Ephys Validation of Human Cluster Identity", fontsize=12, fontweight="bold")

baseline = (mouse_labels_clean == "Mouse CCKBC").sum() / len(mouse_labels_clean)
ax.axhline(baseline, ls=":", color="darkblue", alpha=0.5, label=f"Mouse CCKBC base rate ({baseline:.2f})")
ax.axhline(0.5, ls="--", color="grey", alpha=0.5, label="50% threshold")
ax.legend(fontsize=9); ax.set_ylim(-0.05, 1.05)

plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_cross_species_nn_validation.pdf"), transparent=True, dpi=150)
plt.show()

# %%
# Violin plots: all 16 matched features, raw values (mouse and human shown separately)
# Mouse CCKBC vs VIP-other p-values for annotation
mouse_pvals = {}
for mfeat, hfeat in full_feature_map.items():
    cv = pd.to_numeric(mouse_cge.loc[mouse_cge["is_cckbc"], mfeat], errors="coerce").dropna()
    vv = pd.to_numeric(mouse_cge.loc[mouse_cge["is_vip_other"], mfeat], errors="coerce").dropna()
    if len(cv) >= 3 and len(vv) >= 3:
        _, p = mannwhitneyu(cv, vv)
        mouse_pvals[hfeat] = p

# Build raw-value DataFrames for mouse and human (no z-scoring)
mouse_raw = mouse_feat_all.copy()
mouse_raw["group"] = mouse_labels_all
human_raw = human_feat_all.copy()
human_raw["group"] = human_labels_all

# 4x4 grid: left 2 violins = mouse (CCKBC, VIP-other), right 3 = human groups
ncols, nrows = 4, 4
fig, axes = plt.subplots(nrows, ncols, figsize=(18, 16))
fig.patch.set_alpha(0)

mouse_grps = ["Mouse CCKBC", "Mouse VIP-other"]
human_grps = ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]
xlabels = ["M CCKBC", "M VIP-o", "H VIP-\nCCKBC", "H VIP+\nCCKBC", "H VIP+\nISI"]
all_grps = mouse_grps + human_grps
grp_colors = [GROUP_COLORS[g] for g in all_grps]

for i, (mfeat, hfeat) in enumerate(full_feature_map.items()):
    ax = axes[i // ncols, i % ncols]; ax.patch.set_alpha(0)

    # Gather raw data per group
    data_by_group = []
    for g in mouse_grps:
        vals = pd.to_numeric(mouse_raw.loc[mouse_raw["group"] == g, hfeat], errors="coerce").dropna().values
        data_by_group.append(vals)
    for g in human_grps:
        vals = pd.to_numeric(human_raw.loc[human_raw["group"] == g, hfeat], errors="coerce").dropna().values
        data_by_group.append(vals)

    # Skip if any group is empty
    if any(len(d) == 0 for d in data_by_group):
        ax.set_visible(False)
        continue

    parts = ax.violinplot(data_by_group, positions=range(len(all_grps)),
                          showmeans=True, showextrema=False)
    for j, (pc, color) in enumerate(zip(parts["bodies"], grp_colors)):
        pc.set_facecolor(color); pc.set_alpha(0.6)
    parts["cmeans"].set_color("black")

    # Vertical separator between mouse and human
    ax.axvline(1.5, ls=":", color="grey", alpha=0.4)

    ax.set_xticks(range(len(all_grps)))
    ax.set_xticklabels(xlabels, fontsize=6, rotation=0)

    # Title with mouse p-value
    p = mouse_pvals.get(hfeat, 1.0)
    p_str = f"p={p:.1e}" if p < 0.01 else f"p={p:.2f}"
    star = " ***" if p < 0.001 else " **" if p < 0.01 else " *" if p < 0.05 else ""
    ax.set_title(f"{mfeat}\n(mouse {p_str}{star})", fontsize=8,
                 fontweight="bold" if p < 0.01 else "normal")

fig.suptitle("All 16 Matched Ephys Features (raw values)\nBold = discriminative in mouse CCKBC vs VIP-other (p<0.01)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_cross_species_feature_violins.pdf"), transparent=True, dpi=150)
plt.show()

# %%
# Classifier: train on mouse, predict CCKBC probability for human cells
mouse_y = (mouse_labels_clean == "Mouse CCKBC").astype(int).values

lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(mouse_z.values, mouse_y)
human_prob = lr.predict_proba(human_z_cs.values)[:, 1]
human_pred = lr.predict(human_z_cs.values)

print("Logistic Regression (trained on mouse, predicted on human):")
for grp in human_grp_order:
    mask = human_labels_clean.values == grp
    if mask.sum() == 0:
        continue
    rate = human_pred[mask].mean()
    prob = human_prob[mask].mean()
    print(f"  {grp} (n={mask.sum()}): {rate*100:.1f}% classified CCKBC, mean P(CCKBC)={prob:.3f}")

# %% [markdown]
# **Cross-species ephys alignment summary:**
# - 16 matched mouse–human feature pairs identified; 6 are discriminative (p<0.01)
#   between mouse CCKBC and VIP-other: firing rate, sustained firing, ISI CV,
#   waveform adaptation, latency, early firing rate.
# - PCA/UMAP embedding uses the 6 discriminative features; violin plots show all 16.
# - VIP- CCKBC (277, 278) show the **strongest CCKBC electrophysiological signature**,
#   confirming they are "classic" CCKBCs.
# - VIP+ CCKBC (279, 280, 281) are intermediate, consistent with a transitional/novel subtype.
# - VIP+ ISI clusters have the lowest CCKBC neighbor fraction, near the base rate.

# %% [markdown]
# ## 7. Continuous Analysis: CCKBC Fraction vs 22q-Related Bias

# %%
# Test whether CCKBC fraction correlates with bias across VIP clusters
# (restrict to INT-VIP + imputed CCKBC clusters only)
vip_all = sorted(set(list(VIP_idx) + cckbc_clusters))
vip_master = master.loc[vip_all]

fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
fig.patch.set_alpha(0)

corr_results = []
for idx, (name, label) in enumerate(gene_set_labels.items()):
    ax = axes[idx]
    ax.patch.set_alpha(0)

    col = f"bias_{name}"
    valid_mask = vip_master["cckbc_frac_harmony"].notna() & vip_master[col].notna()
    x = vip_master.loc[valid_mask, "cckbc_frac_harmony"]
    y = vip_master.loc[valid_mask, col]
    rho, p = spearmanr(x, y)
    corr_results.append({"gene_set": label, "rho": rho, "p": p, "n": len(x)})

    colors = ["#E74C3C" if c in cckbc_clusters else "#3498DB" for c in x.index]
    ax.scatter(x, y, s=50, alpha=0.7, c=colors, edgecolors="k", linewidths=0.5)
    for cl in x.index:
        ax.annotate(str(cl), (x[cl], y[cl]), fontsize=7, ha="left", alpha=0.7)

    ax.set_xlabel("CCKBC fraction (Harmony)")
    ax.set_ylabel(f"Bias (EFFECT)" if idx == 0 else "")
    ax.set_title(f"{label}\nrho={rho:.2f}, p={p:.3f}")

    if idx == 0:
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', markersize=8, label='Imputed CCKBC'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498DB', markersize=8, label='VIP+ (ISI)'),
        ]
        ax.legend(handles=legend_elements, loc="best", framealpha=0.5, fontsize=8)

plt.tight_layout()
plt.savefig(str(CGE_DIR / "results" / "fig_cckbc_frac_vs_22q_bias.pdf"), transparent=True)
plt.show()

corr_df = pd.DataFrame(corr_results)
print(corr_df.to_string(index=False))

# %% [markdown]
# ## 8. Summary Table

# %%
# Final summary: key columns for VIP + CCKBC clusters
summary_cols = ["cluster_name", "cckbc_frac_harmony", "cckbc_frac_scvi",
                "ephys_score", "n_ephys", "is_imputed_cckbc"]
summary_cols += [f"bias_{name}" for name in bias_sets]

vip_all = sorted(set(list(VIP_idx) + cckbc_clusters))
summary = master.loc[vip_all, summary_cols].sort_values("cckbc_frac_harmony", ascending=False)
summary.to_csv(CGE_DIR / "results" / "cckbc_convergent_bias_summary.csv")
print("Saved: cckbc_convergent_bias_summary.csv")
print()
print(summary.to_string())

# %% [markdown]
# ## 9. Conclusion
#
# **Summary**:
# - Mouse CCKBCs map consistently to human clusters 277–281 via both
#   Harmony and scVI cross-species integration.
# - Ephys signature transfer independently confirms these clusters as most CCKBC-like
#   (Spearman rho with Harmony fraction, p<0.05).
# - 22q11.2 deletion bias: imputed CCKBC vs VIP+ (ISI) clusters show
#   [result from section 6 above].
#
# **Caveats**:
# - CCKBC is a mouse-defined cell type with no discrete human transcriptomic homolog
#   (Bakken et al. 2021; Darmanis et al. 2015).
# - Sample sizes are small (5 imputed CCKBC clusters vs 11 VIP+ ISI).
# - Results should be treated as exploratory, reviewer-only evidence.
