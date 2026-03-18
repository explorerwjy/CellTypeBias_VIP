"""
Module E: Multi-modal CCKBC classification + updated 22q bias.

Steps:
1. Load MetaNeighbor AUROC → find best-matching mouse cluster per human CGE cluster (276-296)
2. Determine if best-matching mouse clusters are Sncg subclass (module_a_is_sncg)
3. Load existing classification data for marker evidence (CCK+, SNCG+)
4. Build evidence table and classify via classify_clusters()
5. Compute 22q bias per cluster per gene set using HumanCT_AvgZ_Weighted
6. Group comparisons (Mann-Whitney U) and Spearman correlation
7. Save results to cge_subtype/results/updated_22q_bias/
"""

import sys
import os
import warnings
warnings.filterwarnings("ignore")

# Add project root to path
REPO_ROOT = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP"
sys.path.insert(0, REPO_ROOT)

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr

from src.CellType_PSY import HumanCT_AvgZ_Weighted
from cge_subtype.src.multimodal_classification import classify_clusters, EVIDENCE_COLUMNS

# Paths
AUROC_PATH = os.path.join(REPO_ROOT, "cge_subtype/results/cluster_bridge/metaneighbor_auroc_cge.csv")
ALLEN_META_PATH = "/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv"
BIAS_SUMMARY_PATH = os.path.join(REPO_ROOT, "cge_subtype/results/cckbc_convergent_bias_summary.csv")
CLASSIFICATION_PATH = os.path.join(REPO_ROOT, "cge_subtype/results/convergent_cckbc_classification.csv")
PSEUDOBULK_HUMAN_PATH = os.path.join(REPO_ROOT, "cge_subtype/results/pseudobulk/human_pseudobulk.csv")
SPEC_MAT_PATH = os.path.join(REPO_ROOT, "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv")
GW_DIR = os.path.join(REPO_ROOT, "dat/GeneWeights")
OUT_DIR = os.path.join(REPO_ROOT, "cge_subtype/results/updated_22q_bias")

os.makedirs(OUT_DIR, exist_ok=True)

CGE_CLUSTER_IDS = list(range(276, 297))  # 276-296 inclusive

print("=" * 70)
print("MODULE E: Multi-modal CCKBC Classification + Updated 22q Bias")
print("=" * 70)

# =============================================================================
# STEP 1: Load MetaNeighbor AUROC and find best-matching mouse cluster
# =============================================================================
print("\n--- Step 1: MetaNeighbor AUROC → best mouse cluster per human CGE cluster ---")

auroc = pd.read_csv(AUROC_PATH, index_col=0)
print(f"AUROC matrix: {auroc.shape[0]} mouse clusters x {auroc.shape[1]} human CGE clusters")
print(f"Human CGE clusters in AUROC: {sorted([int(c) for c in auroc.columns])}")

# Best-matching mouse cluster per human CGE cluster (column argmax)
best_mouse_per_human = auroc.idxmax(axis=0).astype(float)
best_auroc_per_human = auroc.max(axis=0)

print("\nBest mouse cluster per human CGE cluster:")
for hc in [str(c) for c in CGE_CLUSTER_IDS]:
    if hc in best_mouse_per_human.index:
        mc = best_mouse_per_human[hc]
        auc = best_auroc_per_human[hc]
        print(f"  Human {hc} → Mouse {int(mc)} (AUROC={auc:.3f})")

# =============================================================================
# STEP 2: Look up subclass for best-matching mouse clusters
# =============================================================================
print("\n--- Step 2: Allen metadata subclass lookup ---")

# Load only needed columns from large CSV
print("Loading Allen metadata (this may take a moment)...")
allen_meta = pd.read_csv(
    ALLEN_META_PATH,
    usecols=["cluster_alias", "subclass", "cluster"],
)
# Build per-cluster lookup (cluster_alias → subclass, cluster name)
cluster_info = allen_meta.groupby("cluster_alias")[["subclass", "cluster"]].first()
print(f"Loaded {len(cluster_info)} unique mouse clusters from Allen metadata")

# Map best mouse clusters to subclass
best_cluster_subclass = {}
module_a_is_sncg = {}

for hc in [str(c) for c in CGE_CLUSTER_IDS]:
    if hc not in best_mouse_per_human.index:
        module_a_is_sncg[int(hc)] = False
        best_cluster_subclass[int(hc)] = "N/A"
        continue
    mc_id = int(best_mouse_per_human[hc])
    if mc_id in cluster_info.index:
        subclass = cluster_info.loc[mc_id, "subclass"]
        cluster_name = cluster_info.loc[mc_id, "cluster"]
        is_sncg = "Sncg" in str(subclass) or "sncg" in str(subclass).lower()
        best_cluster_subclass[int(hc)] = f"{subclass} ({cluster_name})"
        module_a_is_sncg[int(hc)] = is_sncg
    else:
        best_cluster_subclass[int(hc)] = "not_found"
        module_a_is_sncg[int(hc)] = False

print("\nSubclass assignments (module_a_is_sncg):")
for hc in CGE_CLUSTER_IDS:
    is_sncg = module_a_is_sncg.get(hc, False)
    subclass = best_cluster_subclass.get(hc, "N/A")
    flag = "YES" if is_sncg else "no"
    print(f"  Human {hc}: {subclass} → Sncg={flag}")

# =============================================================================
# STEP 3: Load existing classification data for marker evidence
# =============================================================================
print("\n--- Step 3: Load existing classification data + pseudobulk markers ---")

# Load convergent classification (has existing subtypes, ephys groups)
existing_cls = pd.read_csv(CLASSIFICATION_PATH)
existing_cls = existing_cls.set_index("cluster_id")

# Load bias summary (has cckbc_frac_harmony, is_imputed_cckbc)
bias_summary = pd.read_csv(BIAS_SUMMARY_PATH, index_col=0)
# The unnamed column is the cluster_id
bias_summary.index = bias_summary.index.astype(int)

# Load human pseudobulk for CCK and SNCG expression
pb_human = pd.read_csv(PSEUDOBULK_HUMAN_PATH, index_col=0)
pb_human.index = pb_human.index.astype(int)

# Define marker thresholds
# CCK: very low expression; cluster 276 is the only one with detectable CCK
# Use >0 (any expression) as CCK positive threshold at pseudobulk level
# SNCG: use median as threshold (above-median = relatively positive)
cck_values = pb_human.loc[
    [i for i in CGE_CLUSTER_IDS if i in pb_human.index], "CCK"
]
sncg_values = pb_human.loc[
    [i for i in CGE_CLUSTER_IDS if i in pb_human.index], "SNCG"
]

# Thresholds: CCK > 0 (any detectable); SNCG > median across CGE clusters
CCK_THRESH = 0.0  # strictly > 0 means any expression
SNCG_THRESH = sncg_values.median()

print(f"CCK threshold: > {CCK_THRESH} (any detectable expression)")
print(f"SNCG threshold: > {SNCG_THRESH:.4f} (median across CGE clusters)")
print(f"\nCCK-positive clusters: {sorted([i for i in CGE_CLUSTER_IDS if i in pb_human.index and pb_human.loc[i, 'CCK'] > CCK_THRESH])}")
print(f"SNCG-positive clusters (above median): {sorted([i for i in CGE_CLUSTER_IDS if i in pb_human.index and pb_human.loc[i, 'SNCG'] > SNCG_THRESH])}")

# =============================================================================
# STEP 4: Build evidence table
# =============================================================================
print("\n--- Step 4: Build multi-modal evidence table ---")

evidence_rows = []
for hc in CGE_CLUSTER_IDS:
    row = {"cluster_id": hc}

    # Module A: best matching mouse cluster is Sncg subclass
    row["module_a_is_sncg"] = module_a_is_sncg.get(hc, False)

    # Marker evidence from pseudobulk
    if hc in pb_human.index:
        row["marker_cck_positive"] = bool(pb_human.loc[hc, "CCK"] > CCK_THRESH)
        row["marker_sncg_positive"] = bool(pb_human.loc[hc, "SNCG"] > SNCG_THRESH)
    else:
        row["marker_cck_positive"] = False
        row["marker_sncg_positive"] = False

    # Module C, D: not yet computed → False
    row["module_c_direct_cckbc"] = False
    row["module_c_indirect_cckbc"] = False
    row["module_d_fast_spiking"] = False

    # Supplementary info from existing analyses
    if hc in existing_cls.index:
        ec = existing_cls.loc[hc]
        row["subtype"] = ec.get("subtype", np.nan)
        row["ephys_group"] = ec.get("ephys_group", np.nan)
        row["mouse_cckbc_fraction_harmony"] = ec.get("mouse_cckbc_fraction", np.nan)
    else:
        row["subtype"] = np.nan
        row["ephys_group"] = np.nan
        row["mouse_cckbc_fraction_harmony"] = np.nan

    if hc in bias_summary.index:
        bs = bias_summary.loc[hc]
        row["cckbc_frac_harmony_summary"] = bs.get("cckbc_frac_harmony", np.nan)
        row["is_imputed_cckbc_prev"] = bs.get("is_imputed_cckbc", np.nan)
    else:
        row["cckbc_frac_harmony_summary"] = np.nan
        row["is_imputed_cckbc_prev"] = np.nan

    row["best_mouse_cluster"] = int(best_mouse_per_human.get(str(hc), np.nan)) if str(hc) in best_mouse_per_human.index else np.nan
    row["best_mouse_subclass"] = best_cluster_subclass.get(hc, "N/A")
    row["best_auroc"] = float(best_auroc_per_human.get(str(hc), np.nan)) if str(hc) in best_auroc_per_human.index else np.nan

    evidence_rows.append(row)

evidence_df = pd.DataFrame(evidence_rows).set_index("cluster_id")
print(f"Evidence table: {evidence_df.shape[0]} clusters x {evidence_df.shape[1]} columns")

# =============================================================================
# STEP 5: Classify clusters
# =============================================================================
print("\n--- Step 5: Multi-modal classification ---")

classified_df = classify_clusters(evidence_df)

print("\nClassification results:")
print(f"{'Cluster':>8} {'module_a':>9} {'mrkr_cck':>9} {'mrkr_sncg':>10} {'confidence':>11} {'tier':>25} {'subtype':>15}")
print("-" * 90)
for cid, row in classified_df.iterrows():
    print(f"{cid:>8} {str(row['module_a_is_sncg']):>9} {str(row['marker_cck_positive']):>9} "
          f"{str(row['marker_sncg_positive']):>10} {row['cckbc_confidence']:>11} "
          f"{row['cckbc_tier']:>25} {str(row.get('subtype', 'N/A')):>15}")

tier_counts = classified_df["cckbc_tier"].value_counts()
print(f"\nTier distribution:")
for tier, count in tier_counts.items():
    print(f"  {tier}: {count} clusters")

# =============================================================================
# STEP 6: Compute 22q bias per cluster per gene set
# =============================================================================
print("\n--- Step 6: Compute 22q bias per CGE cluster per gene set ---")

# Load specificity matrix (genes x cell types)
print("Loading specificity matrix...")
spec_mat = pd.read_csv(SPEC_MAT_PATH, index_col=0)
print(f"Spec matrix: {spec_mat.shape[0]} genes x {spec_mat.shape[1]} cell types")

# Subset to CGE clusters 276-296
cge_cols = [str(c) for c in CGE_CLUSTER_IDS if str(c) in spec_mat.columns]
spec_cge = spec_mat[cge_cols]
print(f"CGE-subset spec matrix: {spec_cge.shape[0]} genes x {spec_cge.shape[1]} clusters")

# Find 22q gene weight files
gw_22q_files = {
    "22q_del": os.path.join(GW_DIR, "X22q.gw.csv"),
    "22q_mousemodel": os.path.join(GW_DIR, "X22q.mousemodel.gw.csv"),
}

bias_records = {}

for gs_name, gw_path in gw_22q_files.items():
    print(f"\n  Gene set: {gs_name}")
    gw_df = pd.read_csv(gw_path, header=None)
    gw_df.columns = ["gene_id", "weight"]

    # Build Gene2Weights dict; gene IDs are Entrez, spec_mat index is likely gene symbols
    # Check if spec_mat index matches Entrez or gene symbols
    gene2weights = dict(zip(gw_df["gene_id"].astype(int).astype(str), gw_df["weight"].astype(float)))

    # Try string Entrez IDs as index
    # Spec matrix index dtype:
    spec_index_sample = spec_cge.index[:5].tolist()

    # Match: try both int and str Entrez
    valid_genes_str = [g for g in gene2weights.keys() if g in spec_cge.index.astype(str).values]
    gene2weights_matched = {g: gene2weights[str(g)] for g in spec_cge.index.astype(str).values if str(g) in gene2weights}

    # Reindex spec matrix with string index
    spec_str = spec_cge.copy()
    spec_str.index = spec_str.index.astype(str)

    print(f"    Gene set size: {len(gw_df)}, matched to spec matrix: {len(gene2weights_matched)}")

    if len(gene2weights_matched) == 0:
        print(f"    WARNING: No genes matched! Skipping {gs_name}")
        continue

    # Compute bias for each CGE cluster
    bias_result = HumanCT_AvgZ_Weighted(spec_str, gene2weights_matched)

    # Extract bias scores for CGE clusters only
    cge_bias = {}
    for hc in CGE_CLUSTER_IDS:
        ct_key = str(hc)
        if ct_key in bias_result.index:
            cge_bias[hc] = bias_result.loc[ct_key, "EFFECT"]
        else:
            cge_bias[hc] = np.nan

    bias_records[gs_name] = cge_bias
    print(f"    Bias range across CGE clusters: [{min(cge_bias.values()):.4f}, {max(cge_bias.values()):.4f}]")

# Build bias dataframe
bias_df = pd.DataFrame(bias_records)
bias_df.index.name = "cluster_id"
print(f"\nBias table shape: {bias_df.shape}")
print("\nBias per cluster (22q_del):")
print(bias_df["22q_del"].sort_values(ascending=False).to_string())

# =============================================================================
# STEP 7: Group comparisons and correlations
# =============================================================================
print("\n--- Step 7: Group comparisons and correlations ---")

# Merge classification with bias
merged = classified_df.join(bias_df, how="left")

# Identify groups
high_conf_cckbc = merged[merged["cckbc_tier"] == "high-confidence CCKBC"]
tentative_cckbc = merged[merged["cckbc_tier"] == "tentative CCKBC"]
isi_vip = merged[merged["cckbc_tier"] == "ISI VIP"]

print(f"\nGroup sizes:")
print(f"  High-confidence CCKBC: {len(high_conf_cckbc)} clusters → {sorted(high_conf_cckbc.index.tolist())}")
print(f"  Tentative CCKBC: {len(tentative_cckbc)} clusters → {sorted(tentative_cckbc.index.tolist())}")
print(f"  ISI VIP: {len(isi_vip)} clusters → {sorted(isi_vip.index.tolist())}")

stats_records = []

for gs_name in gw_22q_files.keys():
    if gs_name not in bias_df.columns:
        continue

    print(f"\n  === Gene set: {gs_name} ===")

    hc_bias = high_conf_cckbc[gs_name].dropna().values
    tent_bias = tentative_cckbc[gs_name].dropna().values
    isi_bias = isi_vip[gs_name].dropna().values
    all_bias = merged[gs_name].dropna().values
    all_conf = merged.loc[merged[gs_name].notna(), "cckbc_confidence"].values

    # Print group means
    if len(hc_bias) > 0:
        print(f"  High-conf CCKBC: n={len(hc_bias)}, mean={np.mean(hc_bias):.4f}, median={np.median(hc_bias):.4f}")
    if len(tent_bias) > 0:
        print(f"  Tentative CCKBC: n={len(tent_bias)}, mean={np.mean(tent_bias):.4f}, median={np.median(tent_bias):.4f}")
    if len(isi_bias) > 0:
        print(f"  ISI VIP:         n={len(isi_bias)}, mean={np.mean(isi_bias):.4f}, median={np.median(isi_bias):.4f}")

    # 2-way: high-confidence CCKBC vs ISI VIP
    if len(hc_bias) >= 2 and len(isi_bias) >= 2:
        stat, p = mannwhitneyu(hc_bias, isi_bias, alternative="two-sided")
        print(f"\n  Mann-Whitney U (high-conf CCKBC vs ISI VIP): U={stat:.0f}, p={p:.4f}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "MWU_highconf_vs_ISI",
            "group1": "high-confidence CCKBC",
            "group2": "ISI VIP",
            "n1": len(hc_bias),
            "n2": len(isi_bias),
            "mean1": np.mean(hc_bias),
            "mean2": np.mean(isi_bias),
            "U_stat": stat,
            "p_value": p,
        })
    else:
        print(f"\n  Mann-Whitney U (high-conf vs ISI): insufficient data (n_hc={len(hc_bias)}, n_isi={len(isi_bias)})")
        stats_records.append({
            "gene_set": gs_name,
            "test": "MWU_highconf_vs_ISI",
            "group1": "high-confidence CCKBC",
            "group2": "ISI VIP",
            "n1": len(hc_bias),
            "n2": len(isi_bias),
            "mean1": np.mean(hc_bias) if len(hc_bias) > 0 else np.nan,
            "mean2": np.mean(isi_bias) if len(isi_bias) > 0 else np.nan,
            "U_stat": np.nan,
            "p_value": np.nan,
        })

    # 3-way ANOVA-like: compare all tiers using Kruskal-Wallis
    groups_3way = []
    labels_3way = []
    if len(hc_bias) >= 2:
        groups_3way.append(hc_bias)
        labels_3way.append("high-conf CCKBC")
    if len(tent_bias) >= 2:
        groups_3way.append(tent_bias)
        labels_3way.append("tentative CCKBC")
    if len(isi_bias) >= 2:
        groups_3way.append(isi_bias)
        labels_3way.append("ISI VIP")

    if len(groups_3way) >= 2:
        kw_stat, kw_p = stats.kruskal(*groups_3way)
        print(f"  Kruskal-Wallis ({' vs '.join(labels_3way)}): H={kw_stat:.3f}, p={kw_p:.4f}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "KruskalWallis_3way",
            "group1": " + ".join(labels_3way),
            "group2": "",
            "n1": sum(len(g) for g in groups_3way),
            "n2": np.nan,
            "mean1": np.nan,
            "mean2": np.nan,
            "U_stat": kw_stat,
            "p_value": kw_p,
        })

    # Pairwise MWU: tentative vs ISI
    if len(tent_bias) >= 2 and len(isi_bias) >= 2:
        stat_t, p_t = mannwhitneyu(tent_bias, isi_bias, alternative="two-sided")
        print(f"  Mann-Whitney U (tentative CCKBC vs ISI VIP): U={stat_t:.0f}, p={p_t:.4f}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "MWU_tentative_vs_ISI",
            "group1": "tentative CCKBC",
            "group2": "ISI VIP",
            "n1": len(tent_bias),
            "n2": len(isi_bias),
            "mean1": np.mean(tent_bias),
            "mean2": np.mean(isi_bias),
            "U_stat": stat_t,
            "p_value": p_t,
        })

    # High-conf + tentative vs ISI
    cckbc_combined = np.concatenate([hc_bias, tent_bias])
    if len(cckbc_combined) >= 2 and len(isi_bias) >= 2:
        stat_c, p_c = mannwhitneyu(cckbc_combined, isi_bias, alternative="two-sided")
        print(f"  Mann-Whitney U (all CCKBC tiers vs ISI VIP): U={stat_c:.0f}, p={p_c:.4f}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "MWU_allCCKBC_vs_ISI",
            "group1": "high-conf+tentative CCKBC",
            "group2": "ISI VIP",
            "n1": len(cckbc_combined),
            "n2": len(isi_bias),
            "mean1": np.mean(cckbc_combined),
            "mean2": np.mean(isi_bias),
            "U_stat": stat_c,
            "p_value": p_c,
        })

    # Spearman: confidence score vs 22q bias
    valid_mask = ~np.isnan(all_bias) & ~np.isnan(all_conf)
    if valid_mask.sum() >= 4:
        rho, pval = spearmanr(all_conf[valid_mask], all_bias[valid_mask])
        print(f"  Spearman (CCKBC confidence vs 22q bias): rho={rho:.3f}, p={pval:.4f}, N={valid_mask.sum()}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "Spearman_confidence_vs_bias",
            "group1": "CCKBC confidence score",
            "group2": "22q bias",
            "n1": int(valid_mask.sum()),
            "n2": np.nan,
            "mean1": np.nan,
            "mean2": np.nan,
            "U_stat": rho,
            "p_value": pval,
        })

    # Also correlation with harmony cckbc fraction
    harm_frac = merged.loc[merged[gs_name].notna(), "cckbc_frac_harmony_summary"].values
    bias_vals = merged.loc[merged[gs_name].notna(), gs_name].values
    valid_mask2 = ~np.isnan(harm_frac) & ~np.isnan(bias_vals)
    if valid_mask2.sum() >= 4:
        rho2, pval2 = spearmanr(harm_frac[valid_mask2], bias_vals[valid_mask2])
        print(f"  Spearman (Harmony CCKBC frac vs 22q bias): rho={rho2:.3f}, p={pval2:.4f}, N={valid_mask2.sum()}")
        stats_records.append({
            "gene_set": gs_name,
            "test": "Spearman_harmonyFrac_vs_bias",
            "group1": "Harmony CCKBC fraction",
            "group2": "22q bias",
            "n1": int(valid_mask2.sum()),
            "n2": np.nan,
            "mean1": np.nan,
            "mean2": np.nan,
            "U_stat": rho2,
            "p_value": pval2,
        })

stats_df = pd.DataFrame(stats_records)

# =============================================================================
# STEP 8: Save results
# =============================================================================
print("\n--- Step 8: Saving results ---")

# 1. Multimodal classification
cls_out = classified_df.copy()
cls_out.to_csv(os.path.join(OUT_DIR, "multimodal_classification.csv"))
print(f"Saved: {OUT_DIR}/multimodal_classification.csv")

# 2. 22q bias per cluster per gene set
bias_out = bias_df.copy()
bias_out.to_csv(os.path.join(OUT_DIR, "22q_bias_by_cluster.csv"))
print(f"Saved: {OUT_DIR}/22q_bias_by_cluster.csv")

# 3. Stats
stats_df.to_csv(os.path.join(OUT_DIR, "group_comparison_stats.csv"), index=False)
print(f"Saved: {OUT_DIR}/group_comparison_stats.csv")

# 4. Human-readable summary
summary_lines = []
summary_lines.append("=" * 70)
summary_lines.append("MODULE E ANALYSIS SUMMARY")
summary_lines.append(f"Date: 2026-03-18")
summary_lines.append("=" * 70)
summary_lines.append("")
summary_lines.append("CLASSIFICATION")
summary_lines.append("-" * 40)
summary_lines.append(f"  Total CGE clusters analyzed: {len(classified_df)}")
for tier, count in tier_counts.items():
    clusters_in_tier = sorted(classified_df[classified_df["cckbc_tier"] == tier].index.tolist())
    summary_lines.append(f"  {tier}: {count} clusters → {clusters_in_tier}")

summary_lines.append("")
summary_lines.append("EVIDENCE SUMMARY (per cluster)")
summary_lines.append("-" * 40)
for cid, row in classified_df.iterrows():
    summary_lines.append(
        f"  Cluster {cid}: conf={row['cckbc_confidence']}, tier={row['cckbc_tier']}, "
        f"module_a={row['module_a_is_sncg']}, cck={row['marker_cck_positive']}, "
        f"sncg={row['marker_sncg_positive']}"
    )

summary_lines.append("")
summary_lines.append("22q BIAS RESULTS")
summary_lines.append("-" * 40)
for gs_name in gw_22q_files.keys():
    if gs_name not in bias_df.columns:
        continue
    summary_lines.append(f"  Gene set: {gs_name}")
    for tier in ["high-confidence CCKBC", "tentative CCKBC", "ISI VIP"]:
        tier_clusters = classified_df[classified_df["cckbc_tier"] == tier].index
        tier_bias = bias_df.loc[tier_clusters, gs_name].dropna()
        if len(tier_bias) > 0:
            summary_lines.append(
                f"    {tier}: n={len(tier_bias)}, mean={np.mean(tier_bias):.4f}, "
                f"median={np.median(tier_bias):.4f}, range=[{np.min(tier_bias):.4f}, {np.max(tier_bias):.4f}]"
            )

summary_lines.append("")
summary_lines.append("STATISTICAL TESTS")
summary_lines.append("-" * 40)
for _, row in stats_df.iterrows():
    if pd.isna(row["p_value"]):
        summary_lines.append(f"  [{row['gene_set']}] {row['test']}: insufficient data")
    else:
        stat_val = row["U_stat"]
        summary_lines.append(
            f"  [{row['gene_set']}] {row['test']}: stat={stat_val:.3f}, p={row['p_value']:.4f}"
        )

summary_lines.append("")
summary_lines.append("NOTES")
summary_lines.append("-" * 40)
summary_lines.append(f"  CCK threshold: > {CCK_THRESH} (pseudobulk TPM-normalized)")
summary_lines.append(f"  SNCG threshold: > {SNCG_THRESH:.4f} (median across CGE clusters)")
summary_lines.append(f"  Spec matrix: {SPEC_MAT_PATH}")
summary_lines.append("  Module C (direct/indirect CCKBC path), Module D (ephys fast-spiking): not yet computed → False")
summary_lines.append("  Only Modules A + Markers available for this run.")

summary_text = "\n".join(summary_lines)
print("\n" + summary_text)

with open(os.path.join(OUT_DIR, "analysis_summary.txt"), "w") as f:
    f.write(summary_text)
print(f"\nSaved: {OUT_DIR}/analysis_summary.txt")

print("\n" + "=" * 70)
print("MODULE E COMPLETE")
print("=" * 70)
