#!/usr/bin/env python3
"""
Module D: Ephys Convergence Analysis
=====================================
Cross-species electrophysiological convergence test for CGE interneuron subtypes.

ID matching strategy
--------------------
Mouse cells:
  - Combined ephys cell_id format: NWB BIDS prefix, e.g.
      'sub_mouse-AVLEY_ses_20190515-sample-4_...'
  - patchseq_mapping_results.csv (Unnamed: 0) contains two sub-datasets:
      * M1 (dataset=='M1'): Allen Institute Patch-seq IDs like '20171204_sample_2'
        -> matched to NWB prefix via m1_patchseq_meta_data.tsv
      * V1 (dataset=='V1'): PS-style IDs like 'PS0810_E1-50_S88'
        -> no NWB ephys available; these are excluded
  - Only M1 cells have matching ephys.

Human cells:
  - Ephys cell_id format: BIDS NWB, e.g. 'sub-964127089_ses-967660413_icephys'
  - human_scvi_mapping_results.csv (Unnamed: 0): specimen_ids (Allen Institute)
  - Crosswalk: human_patchseq_query_preprocessed.h5ad obs['nwb_id']
      specimen_id -> nwb_id (= ses-XXXX in BIDS)
  - Match: map_results[specimen_id] -> h5ad[nwb_id] -> ephys[ses-XXXX]

Usage
-----
    conda run -n gencic python cge_subtype/scripts/run_ephys_convergence.py

Outputs (all in cge_subtype/results/ephys_convergence/):
  - ephys_convergence_results.csv
  - ephys_permutation_pvalue.txt
  - ephys_convergence_summary.txt
"""

import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
sys.path.insert(0, str(REPO_ROOT))

COMBINED_EPHYS = REPO_ROOT / "cge_subtype/results/ephys_convergence/combined_ephys_features.csv"
MOUSE_MAP_PATH = REPO_ROOT / "cge_subtype/results/patchseq_mapping_results.csv"
HUMAN_MAP_PATH = REPO_ROOT / "cge_subtype/results/human_scvi_mapping_results.csv"
M1_META_PATH = Path("/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv")
HUMAN_QUERY_H5AD = Path(
    "/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results"
    "/human_interneuron/preprocessed/human_patchseq_query_preprocessed.h5ad"
)
OUTDIR = REPO_ROOT / "cge_subtype/results/ephys_convergence"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Imports from project source
# ---------------------------------------------------------------------------

from cge_subtype.src.ephys_harmonization import (
    BIO_FEATURES,
    ISI_FEATURES,
    CV_ISI_FEATURES,
    ALL_FEATURES,
    zscore_within_species,
    permutation_test_cluster_similarity,
)

ALL_FEATURE_COLS = BIO_FEATURES + ISI_FEATURES + CV_ISI_FEATURES


# ---------------------------------------------------------------------------
# Step 1. Load combined ephys features
# ---------------------------------------------------------------------------

print("=" * 70)
print("Step 1: Loading combined ephys features")
print("=" * 70)

combined = pd.read_csv(COMBINED_EPHYS)
print(f"  Combined shape: {combined.shape}")
print(f"  Species: {combined['species'].value_counts().to_dict()}")

mouse_ephys = combined[combined["species"] == "mouse"].copy().reset_index(drop=True)
human_ephys = combined[combined["species"] == "human"].copy().reset_index(drop=True)
print(f"  Mouse cells: {len(mouse_ephys)}")
print(f"  Human cells: {len(human_ephys)}")


# ---------------------------------------------------------------------------
# Step 2. Load cluster assignments and merge
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 2: Loading cluster assignments")
print("=" * 70)

# ---- 2a. Mouse M1 cells: match via metadata NWB prefix --------------------

print("\n[Mouse] Loading patchseq mapping …")
mouse_map_raw = pd.read_csv(MOUSE_MAP_PATH)
mouse_map_raw = mouse_map_raw.rename(columns={"Unnamed: 0": "patchseq_id"})

m1_map = mouse_map_raw[mouse_map_raw["dataset"] == "M1"].copy()
v1_map = mouse_map_raw[mouse_map_raw["dataset"] == "V1"].copy()
print(f"  M1 cells in mapping: {len(m1_map)}")
print(f"  V1 cells in mapping: {len(v1_map)}")
print(f"  V1 cells have PS-style IDs (e.g. {v1_map['patchseq_id'].iloc[0]}); "
      f"no NWB ephys — excluded.")

print("\n[Mouse] Building NWB prefix crosswalk from M1 metadata …")
meta = pd.read_csv(M1_META_PATH, sep="\t")
meta["mouse_id"] = meta["Mouse"].str.replace("mouse_", "", regex=False)
meta["date_fmt"] = meta["Date"].apply(
    lambda d: datetime.strptime(d, "%m/%d/%y").strftime("%Y%m%d") if pd.notna(d) else None
)
meta["sample_num"] = meta["Sample"].str.extract(r"(\d+)")[0]
meta["nwb_prefix"] = (
    "sub_mouse-" + meta["mouse_id"] + "_ses_" +
    meta["date_fmt"] + "-sample-" + meta["sample_num"]
)
# meta['Cell'] column matches patchseq_id in m1_map
meta_crosswalk = meta[["Cell", "nwb_prefix"]].rename(columns={"Cell": "patchseq_id"})

m1_map = m1_map.merge(meta_crosswalk, on="patchseq_id", how="left")
print(f"  M1 cells with NWB prefix: {m1_map['nwb_prefix'].notna().sum()} / {len(m1_map)}")

# Match NWB prefix to ephys cell_id (prefix match, since NWB IDs are longer)
mouse_ephys_ids = mouse_ephys["cell_id"].tolist()

def _find_ephys_id(prefix):
    if pd.isna(prefix):
        return None
    for eid in mouse_ephys_ids:
        if eid.startswith(prefix):
            return eid
    return None

m1_map["cell_id"] = m1_map["nwb_prefix"].apply(_find_ephys_id)
n_mouse_matched = m1_map["cell_id"].notna().sum()
print(f"  M1 cells matched to ephys cell_id: {n_mouse_matched} / {len(m1_map)}")

# Build mouse cluster assignment table
mouse_cluster_df = (
    m1_map[m1_map["cell_id"].notna()][["cell_id", "assigned_human_cluster"]]
    .rename(columns={"assigned_human_cluster": "cluster"})
    .copy()
)
print(f"  Mouse cells with cluster + ephys: {len(mouse_cluster_df)}")


# ---- 2b. Human cells: specimen_id -> nwb_id (ses) -> ephys cell_id -------

print("\n[Human] Loading scVI mapping …")
human_map_raw = pd.read_csv(HUMAN_MAP_PATH)
human_map_raw = human_map_raw.rename(columns={"Unnamed: 0": "specimen_id"})
human_map_raw["specimen_id"] = human_map_raw["specimen_id"].astype(str)
print(f"  Human mapping shape: {human_map_raw.shape}")

print("[Human] Loading crosswalk from h5ad obs['nwb_id'] …")
try:
    import anndata as ad
    query_adata = ad.read_h5ad(HUMAN_QUERY_H5AD)
    crosswalk = query_adata.obs[["nwb_id"]].copy()
    crosswalk.index.name = "specimen_id"
    crosswalk = crosswalk.reset_index()
    crosswalk["specimen_id"] = crosswalk["specimen_id"].astype(str)
    crosswalk["ses_id_str"] = (
        crosswalk["nwb_id"]
        .dropna()
        .astype(int)
        .astype(str)
        .reindex(crosswalk.index)
    )
    crosswalk_valid = crosswalk[crosswalk["ses_id_str"].notna()][["specimen_id", "ses_id_str"]]
    print(f"  Crosswalk entries (non-NaN nwb_id): {len(crosswalk_valid)} / {len(crosswalk)}")
except Exception as exc:
    print(f"  WARNING: Could not load h5ad crosswalk: {exc}")
    crosswalk_valid = pd.DataFrame(columns=["specimen_id", "ses_id_str"])

# Merge mapping with crosswalk
human_map_ses = human_map_raw.merge(crosswalk_valid, on="specimen_id", how="left")
n_ses_filled = human_map_ses["ses_id_str"].notna().sum()
print(f"  Human mapping cells with ses_id: {n_ses_filled} / {len(human_map_ses)}")

# Extract ses_id from ephys cell_id
human_ephys["ses_id_str"] = human_ephys["cell_id"].str.extract(r"ses-(\d+)_icephys")[0]

# Merge on ses_id
human_merged = human_map_ses.merge(
    human_ephys[["cell_id", "ses_id_str"]],
    on="ses_id_str",
    how="inner",
)
print(f"  Human cells matched (mapping + ephys): {len(human_merged)}")

# Build human cluster assignment table
human_cluster_df = (
    human_merged[["cell_id", "assigned_cluster"]]
    .rename(columns={"assigned_cluster": "cluster"})
    .copy()
)
print(f"  Human cells with cluster + ephys: {len(human_cluster_df)}")

# ---- 2c. Summary of ID formats ----------------------------------------

print("\n[ID Format Summary]")
print(f"  Mouse ephys cell_id example: {mouse_ephys['cell_id'].iloc[0]}")
print(f"  Mouse mapping ID example (M1): {m1_map['patchseq_id'].iloc[0]}")
print(f"  Human ephys cell_id example: {human_ephys['cell_id'].iloc[0]}")
print(f"  Human mapping ID example: {human_map_raw['specimen_id'].iloc[0]}")
print(f"  Human h5ad nwb_id example: {crosswalk_valid['ses_id_str'].iloc[0] if len(crosswalk_valid) else 'N/A'}")


# ---------------------------------------------------------------------------
# Step 3. Merge ephys features with cluster assignments
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 3: Merging ephys features with cluster assignments")
print("=" * 70)

mouse_annotated = mouse_ephys.merge(mouse_cluster_df, on="cell_id", how="left")
human_annotated = human_ephys.merge(human_cluster_df, on="cell_id", how="left")

# Combined annotated dataset
combined_annotated = pd.concat([mouse_annotated, human_annotated], ignore_index=True)

n_mouse_clus = (mouse_annotated["cluster"].notna()).sum()
n_human_clus = (human_annotated["cluster"].notna()).sum()
print(f"  Mouse cells with cluster: {n_mouse_clus} / {len(mouse_annotated)}")
print(f"  Human cells with cluster: {n_human_clus} / {len(human_annotated)}")
print(f"  Combined annotated: {combined_annotated.shape}")

# Keep only cells with cluster assignments
combined_valid = combined_annotated[combined_annotated["cluster"].notna()].copy()
combined_valid = combined_valid.reset_index(drop=True)
print(f"  Cells with cluster (both species): {len(combined_valid)}")


# ---------------------------------------------------------------------------
# Step 4. Identify clusters with cells from BOTH species
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 4: Identifying dual-species clusters")
print("=" * 70)

cluster_species = (
    combined_valid.groupby("cluster")["species"]
    .apply(lambda x: set(x.unique()))
    .reset_index()
    .rename(columns={"species": "species_set"})
)
cluster_species["has_both"] = cluster_species["species_set"].apply(
    lambda s: "mouse" in s and "human" in s
)
dual_clusters = cluster_species[cluster_species["has_both"]]["cluster"].tolist()
mouse_only = cluster_species[~cluster_species["has_both"] & cluster_species["species_set"].apply(lambda s: "mouse" in s)]["cluster"].tolist()
human_only = cluster_species[~cluster_species["has_both"] & cluster_species["species_set"].apply(lambda s: "human" in s)]["cluster"].tolist()

print(f"  Clusters with BOTH species: {len(dual_clusters)}")
print(f"  Mouse-only clusters: {len(mouse_only)}")
print(f"  Human-only clusters: {len(human_only)}")
print(f"  Dual-species cluster IDs: {sorted(dual_clusters)}")

# Cells counts per dual cluster
for cl in sorted(dual_clusters):
    n_m = ((combined_valid["cluster"] == cl) & (combined_valid["species"] == "mouse")).sum()
    n_h = ((combined_valid["cluster"] == cl) & (combined_valid["species"] == "human")).sum()
    print(f"    Cluster {cl}: mouse={n_m}, human={n_h}")


# ---------------------------------------------------------------------------
# Step 5. Z-score features within species
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 5: Within-species z-scoring")
print("=" * 70)

feat_cols = [c for c in ALL_FEATURE_COLS if c in combined_valid.columns]
print(f"  Feature columns available: {len(feat_cols)}")

# Drop rows missing ALL features
feat_present = combined_valid[feat_cols].notna().any(axis=1)
n_drop = (~feat_present).sum()
if n_drop > 0:
    print(f"  Dropping {n_drop} cells with no feature data")
combined_valid = combined_valid[feat_present].copy().reset_index(drop=True)

features_raw = combined_valid[feat_cols]
species_labels = combined_valid["species"].values

features_zscored = zscore_within_species(features_raw, species_labels)
print(f"  Z-scored features shape: {features_zscored.shape}")

# Sanity check
for sp in ["mouse", "human"]:
    mask = species_labels == sp
    grand_mean = features_zscored.loc[mask].mean().mean()
    print(f"  {sp} grand mean (expect ~0): {grand_mean:.4f}")


# ---------------------------------------------------------------------------
# Step 6. Global permutation test
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 6: Global permutation test (n_perm=1000)")
print("=" * 70)

# Use only dual-species clusters for the test
dual_mask = combined_valid["cluster"].isin(dual_clusters)
perm_features = features_zscored[dual_mask]
perm_clusters = combined_valid.loc[dual_mask, "cluster"].values
perm_species = combined_valid.loc[dual_mask, "species"].values

print(f"  Cells used for permutation test: {len(perm_features)}")
print(f"  Dual-species clusters: {len(np.unique(perm_clusters))}")

N_PERM = 1000
SEED = 42

p_value = permutation_test_cluster_similarity(
    perm_features,
    perm_clusters,
    perm_species,
    n_perm=N_PERM,
    seed=SEED,
)
print(f"\n  Global permutation p-value: {p_value:.4f}")
if p_value < 0.05:
    print("  => Cross-species cluster similarity is GREATER than chance (p < 0.05)")
elif np.isnan(p_value):
    print("  => p-value is NaN (insufficient dual-species clusters)")
else:
    print("  => No significant cross-species cluster similarity (p >= 0.05)")


# ---------------------------------------------------------------------------
# Step 7. Per-cluster Spearman correlation of mean ephys profiles
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 7: Per-cluster Spearman correlation (mouse vs human mean profiles)")
print("=" * 70)

bio_cols = [c for c in BIO_FEATURES if c in feat_cols]
results = []

for cl in sorted(dual_clusters):
    mouse_mask = (combined_valid["cluster"] == cl) & (combined_valid["species"] == "mouse")
    human_mask = (combined_valid["cluster"] == cl) & (combined_valid["species"] == "human")

    n_mouse = mouse_mask.sum()
    n_human = human_mask.sum()

    # Use BIO_FEATURES centroids for Spearman correlation
    mouse_centroid = features_zscored.loc[mouse_mask, bio_cols].mean()
    human_centroid = features_zscored.loc[human_mask, bio_cols].mean()

    # Drop features where either centroid has NaN
    valid_mask = mouse_centroid.notna() & human_centroid.notna()
    if valid_mask.sum() < 3:
        rho, pval = np.nan, np.nan
    else:
        rho, pval = spearmanr(
            mouse_centroid[valid_mask].values,
            human_centroid[valid_mask].values,
        )

    # Euclidean distance in full feature space
    full_mouse = features_zscored.loc[mouse_mask, feat_cols].mean()
    full_human = features_zscored.loc[human_mask, feat_cols].mean()
    valid_full = full_mouse.notna() & full_human.notna()
    euclid_dist = float(np.sqrt(np.nansum((full_mouse - full_human) ** 2)))

    results.append({
        "cluster": cl,
        "n_mouse": n_mouse,
        "n_human": n_human,
        "spearman_rho": rho,
        "spearman_pval": pval,
        "euclidean_distance": euclid_dist,
        "n_features_used": int(valid_mask.sum()),
    })

    print(f"  Cluster {cl}: mouse={n_mouse}, human={n_human}, "
          f"rho={rho:.3f} (p={pval:.3f}), dist={euclid_dist:.3f}")

results_df = pd.DataFrame(results)

# Summary stats
n_pos = (results_df["spearman_rho"] > 0).sum()
n_sig = ((results_df["spearman_rho"] > 0) & (results_df["spearman_pval"] < 0.05)).sum()
mean_rho = results_df["spearman_rho"].mean()
mean_dist = results_df["euclidean_distance"].mean()

print(f"\n  Summary across {len(results_df)} dual-species clusters:")
print(f"    Mean Spearman rho: {mean_rho:.3f}")
print(f"    Clusters with positive rho: {n_pos} / {len(results_df)}")
print(f"    Clusters with rho > 0 and p < 0.05: {n_sig}")
print(f"    Mean Euclidean distance: {mean_dist:.3f}")


# ---------------------------------------------------------------------------
# Step 8. Save results
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Step 8: Saving results")
print("=" * 70)

# 8a. Per-cluster CSV
out_csv = OUTDIR / "ephys_convergence_results.csv"
results_df.to_csv(out_csv, index=False)
print(f"  Saved: {out_csv}")

# 8b. Permutation p-value
out_pval = OUTDIR / "ephys_permutation_pvalue.txt"
with open(out_pval, "w") as fh:
    fh.write(f"{p_value}\n")
print(f"  Saved: {out_pval}")

# 8c. Human-readable summary
out_summary = OUTDIR / "ephys_convergence_summary.txt"
with open(out_summary, "w") as fh:
    fh.write("Module D: Ephys Convergence Analysis\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    fh.write("--- Data ---\n")
    fh.write(f"Combined ephys cells: {len(combined)}\n")
    fh.write(f"  Mouse: {(combined['species']=='mouse').sum()}\n")
    fh.write(f"  Human: {(combined['species']=='human').sum()}\n\n")

    fh.write("--- ID Matching ---\n")
    fh.write(f"Mouse cells matched (M1 patch-seq -> NWB ephys): {n_mouse_matched}\n")
    fh.write(f"  Mouse M1 mapping IDs (e.g.): {m1_map['patchseq_id'].iloc[0]}\n")
    fh.write(f"  Mouse NWB ephys IDs (e.g.): {mouse_ephys['cell_id'].iloc[0][:50]}...\n")
    fh.write(f"  V1 cells excluded (PS-style IDs, no ephys data): {len(v1_map)}\n")
    fh.write(f"Human cells matched (specimen_id -> nwb_id -> BIDS): {len(human_merged)}\n")
    fh.write(f"  Human mapping IDs (e.g.): {human_map_raw['specimen_id'].iloc[0]}\n")
    fh.write(f"  Human ephys IDs (e.g.): {human_ephys['cell_id'].iloc[0]}\n\n")

    fh.write("--- Cluster Annotations ---\n")
    fh.write(f"Mouse cells with cluster assignment: {n_mouse_clus}\n")
    fh.write(f"Human cells with cluster assignment: {n_human_clus}\n")
    fh.write(f"Dual-species clusters (mouse AND human): {len(dual_clusters)}\n")
    fh.write(f"Mouse-only clusters: {len(mouse_only)}\n")
    fh.write(f"Human-only clusters: {len(human_only)}\n\n")

    fh.write("--- Global Permutation Test ---\n")
    fh.write(f"n_permutations: {N_PERM}\n")
    fh.write(f"Cells used: {len(perm_features)}\n")
    fh.write(f"p-value: {p_value:.4f}\n")
    if p_value < 0.05:
        interp = "Significant: cross-species cluster similarity > chance"
    elif np.isnan(p_value):
        interp = "NaN: insufficient dual-species clusters"
    else:
        interp = "Not significant: no detectable cross-species cluster similarity"
    fh.write(f"Interpretation: {interp}\n\n")

    fh.write("--- Per-Cluster Spearman Correlations ---\n")
    fh.write(f"Features used: BIO_FEATURES ({len(bio_cols)} features)\n")
    fh.write(f"Mean Spearman rho across {len(results_df)} dual-species clusters: {mean_rho:.3f}\n")
    fh.write(f"Clusters with rho > 0: {n_pos} / {len(results_df)}\n")
    fh.write(f"Clusters with rho > 0 AND p < 0.05: {n_sig}\n")
    fh.write(f"Mean Euclidean distance (z-scored full features): {mean_dist:.3f}\n\n")

    fh.write("--- Per-Cluster Details ---\n")
    fh.write(results_df.to_string(index=False))
    fh.write("\n")

print(f"  Saved: {out_summary}")
print("\nDone.")
