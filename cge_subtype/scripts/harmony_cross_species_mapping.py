"""
Harmony-based cross-species mapping: mouse CCKBC → human CGE clusters.

Alternative to scVI for Approach 2.
Uses the preprocessed cross-species reference (human+mouse, 2828 HVGs)
and maps mouse patch-seq query cells via k-NN in Harmony-corrected PCA space.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import harmonypy as hm
from scipy.spatial import cKDTree
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# Paths
# ============================================================
BASE = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
REF_HVG = BASE / "results/cross_species_preprocessed/cross_species_reference_hvg.h5ad"
HUMAN_REF = BASE / "results/cross_species_preprocessed/human_reference.h5ad"
QUERY = BASE / "results/cross_species_preprocessed/mouse_patchseq_query.h5ad"
OUT_DIR = BASE / "results"

# ============================================================
# 1. Load reference (already HVG-subsetted, 227K cells x 2828 genes)
# ============================================================
print("Loading reference (HVG subset)...")
ref = sc.read_h5ad(str(REF_HVG))
print(f"  Reference: {ref.shape}")

# Get cluster labels from full human reference
print("Loading human reference for cluster labels...")
human_full = sc.read_h5ad(str(HUMAN_REF), backed="r")
human_labels = human_full.obs[["supercluster_term", "cluster_id"]].copy()
human_full.file.close()

# Merge cluster labels into reference (human cells only)
# Note: HVG reference indices have "-0" suffix; human reference does not
ref.obs["supercluster_term"] = pd.Series(np.nan, index=ref.obs.index, dtype=object)
ref.obs["cluster_id"] = pd.Series(np.nan, index=ref.obs.index, dtype=object)
human_mask = ref.obs["species"] == "human"
human_idx = ref.obs.index[human_mask]
# Strip "-0" suffix for matching
idx_stripped = human_idx.str.replace(r"-\d+$", "", regex=True)
idx_map = pd.Series(human_idx, index=idx_stripped)
shared = idx_stripped.intersection(human_labels.index)
ref.obs.loc[idx_map[shared].values, "supercluster_term"] = human_labels.loc[shared, "supercluster_term"].values
ref.obs.loc[idx_map[shared].values, "cluster_id"] = human_labels.loc[shared, "cluster_id"].values
print(f"  Human cells with cluster labels: {ref.obs['cluster_id'].notna().sum()}")

# ============================================================
# 2. Load mouse patch-seq query
# ============================================================
print("Loading mouse patch-seq query...")
query = sc.read_h5ad(str(QUERY))
print(f"  Query: {query.shape}, CCKBCs: {query.obs['is_cckbc'].sum()}")

# Subset query to reference HVGs
shared_genes = ref.var_names.intersection(query.var_names)
print(f"  Shared HVGs: {len(shared_genes)}")
query_sub = query[:, shared_genes].copy()
ref_sub = ref[:, shared_genes].copy()

# ============================================================
# 3. Concatenate reference + query
# ============================================================
print("Concatenating reference + query...")
query_sub.obs["data_type"] = "query"
ref_sub.obs["data_type"] = "reference"
# Ensure is_cckbc exists in both
if "is_cckbc" not in ref_sub.obs.columns:
    ref_sub.obs["is_cckbc"] = False

combined = sc.concat([ref_sub, query_sub], join="inner")
print(f"  Combined: {combined.shape}")

# ============================================================
# 4. Normalize and PCA
# ============================================================
print("Normalizing and running PCA...")
# Log-normalize (if not already)
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, n_comps=50, svd_solver="arpack")
print(f"  PCA done: {combined.obsm['X_pca'].shape}")

# ============================================================
# 5. Run Harmony (correct for species)
# ============================================================
print("Running Harmony (batch_key=species)...")
ho = hm.run_harmony(
    combined.obsm["X_pca"],
    combined.obs,
    "species",
    max_iter_harmony=30,
    random_state=42,
)
Z_corr = np.asarray(ho.Z_corr)
if Z_corr.shape[0] != combined.n_obs:
    Z_corr = Z_corr.T  # old API returns (d, N)
combined.obsm["X_harmony"] = Z_corr
print(f"  Harmony done: {combined.obsm['X_harmony'].shape}")

# ============================================================
# 6. k-NN assignment: query → human reference
# ============================================================
print("Running k-NN assignment...")
ref_mask = combined.obs["data_type"] == "reference"
human_ref_mask = ref_mask & (combined.obs["species"] == "human")
query_mask = combined.obs["data_type"] == "query"

ref_embedding = combined.obsm["X_harmony"][human_ref_mask.values]
query_embedding = combined.obsm["X_harmony"][query_mask.values]
ref_obs = combined.obs[human_ref_mask].copy()
query_obs = combined.obs[query_mask].copy()

# Merge cluster labels into ref_obs (labels from ref.obs, which has them)
# ref.obs has supercluster_term/cluster_id; combined.obs may have lost them in concat
for col in ["supercluster_term", "cluster_id"]:
    if col not in ref_obs.columns:
        ref_obs[col] = ref.obs.loc[ref_obs.index, col].values

# Build k-d tree on human reference
k = 30
tree = cKDTree(ref_embedding)
distances, indices = tree.query(query_embedding, k=k)

# Assign supercluster and cluster by majority vote
ref_superclusters = ref_obs["supercluster_term"].values
ref_clusters = ref_obs["cluster_id"].values

assignments = []
for i in range(len(query_obs)):
    nn_idx = indices[i]
    nn_sc = ref_superclusters[nn_idx]
    nn_cl = ref_clusters[nn_idx]
    nn_dist = distances[i]

    # Majority vote for supercluster (distance-weighted)
    weights = 1.0 / (nn_dist + 1e-8)
    sc_votes = {}
    for sc_val, w in zip(nn_sc, weights):
        if pd.notna(sc_val):
            sc_votes[sc_val] = sc_votes.get(sc_val, 0) + w
    best_sc = max(sc_votes, key=sc_votes.get) if sc_votes else np.nan
    sc_conf = sc_votes.get(best_sc, 0) / sum(sc_votes.values()) if sc_votes else 0

    # Majority vote for cluster
    cl_votes = {}
    for cl_val, w in zip(nn_cl, weights):
        if pd.notna(cl_val):
            cl_votes[cl_val] = cl_votes.get(cl_val, 0) + w
    best_cl = max(cl_votes, key=cl_votes.get) if cl_votes else np.nan
    cl_conf = cl_votes.get(best_cl, 0) / sum(cl_votes.values()) if cl_votes else 0

    assignments.append({
        "assigned_supercluster": best_sc,
        "conf_supercluster": sc_conf,
        "assigned_cluster": best_cl,
        "conf_cluster": cl_conf,
        "mean_dist": np.mean(nn_dist),
    })

result_df = pd.DataFrame(assignments, index=query_obs.index)
result_df["is_cckbc"] = query_obs["is_cckbc"].values

# Add mouse metadata
for col in ["dataset", "species"]:
    if col in query_obs.columns:
        result_df[col] = query_obs[col].values

# Also transfer RNA family/type if available
for col in ["RNA family", "RNA type", "corresponding_AIT2.3.1_alias"]:
    if col in query.obs.columns:
        # Need to align indices
        result_df[f"mouse_{col}"] = query.obs.loc[result_df.index, col].values

# ============================================================
# 7. Save and print summary
# ============================================================
out_path = OUT_DIR / "harmony_patchseq_mapping_results.csv"
result_df.to_csv(out_path)
print(f"\nSaved to {out_path}")

# Summary
print("\n=== HARMONY MAPPING SUMMARY ===")
print(f"Total mapped cells: {len(result_df)}")
print(f"\nSupercluster distribution (ALL cells):")
print(result_df["assigned_supercluster"].value_counts().to_string())

cckbc = result_df[result_df["is_cckbc"] == True]
print(f"\nCCKBC cells: {len(cckbc)}")
print(f"CCKBC → supercluster:")
print(cckbc["assigned_supercluster"].value_counts().to_string())

cckbc_cge = cckbc[cckbc["assigned_supercluster"] == "CGE interneuron"]
print(f"\nCCKBC → CGE clusters ({len(cckbc_cge)}/{len(cckbc)} = {100*len(cckbc_cge)/len(cckbc):.1f}%):")
if len(cckbc_cge) > 0:
    print(cckbc_cge["assigned_cluster"].value_counts().head(15).to_string())

print(f"\nMean confidence (cluster): {result_df['conf_cluster'].mean():.3f}")
print(f"CCKBC mean confidence (cluster): {cckbc['conf_cluster'].mean():.3f}")

# Specificity: CCKBC fraction per CGE cluster
all_cge = result_df[result_df["assigned_supercluster"] == "CGE interneuron"]
cge_stats = []
for cl_id in sorted(all_cge["assigned_cluster"].unique()):
    cl_cells = all_cge[all_cge["assigned_cluster"] == cl_id]
    n_total = len(cl_cells)
    n_cckbc = cl_cells["is_cckbc"].sum()
    cge_stats.append({"cluster": int(cl_id), "n_total": n_total, "n_cckbc": n_cckbc,
                       "cckbc_frac": n_cckbc/n_total if n_total > 0 else 0})
cge_stats_df = pd.DataFrame(cge_stats).sort_values("cckbc_frac", ascending=False)
print(f"\nCCKBC fraction per CGE cluster (Harmony):")
print(cge_stats_df.to_string(index=False))

# Compare with scVI: key clusters
print(f"\n=== KEY COMPARISON: scVI vs Harmony ===")
for cl in [279, 280, 281, 284, 289, 292]:
    row = cge_stats_df[cge_stats_df["cluster"] == cl]
    if len(row) > 0:
        frac = row["cckbc_frac"].values[0]
        n = row["n_total"].values[0]
        nc = row["n_cckbc"].values[0]
        print(f"  Cluster {cl}: {nc}/{n} CCKBC = {frac:.3f}")
    else:
        print(f"  Cluster {cl}: not mapped")
