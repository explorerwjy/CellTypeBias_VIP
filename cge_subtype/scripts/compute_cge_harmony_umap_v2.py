#!/usr/bin/env python3
"""
Compute Harmony-corrected UMAP for cross-species CGE visualization.

Uses the cross-species HVG reference (227k cells) + mouse patch-seq query.
Runs PCA + Harmony with strong integration (theta=4) and computes UMAP
for visualization comparison with the scVI UMAP.

Output:
  results/harmony_umap_cge_harmony.h5ad
"""

import gc
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import harmonypy as hm
import torch
import scipy.sparse as sp

warnings.filterwarnings("ignore")

_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

BASE = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
REF_HVG = BASE / "results/cross_species_preprocessed/cross_species_reference_hvg.h5ad"
HUMAN_REF = BASE / "results/cross_species_preprocessed/human_reference.h5ad"
QUERY = BASE / "results/cross_species_preprocessed/mouse_patchseq_query.h5ad"
OUT_PATH = BASE / "results" / "harmony_umap_cge_harmony.h5ad"

N_HUMAN_CGE = 25000
N_MOUSE_REF = 8000
SEED = 42
THETA = 4.0  # Stronger integration than default 2.0

print("=" * 60)
print("Harmony UMAP for cross-species CGE")
print(f"  theta = {THETA}")
print("=" * 60)

# 1. Load HVG reference
print(f"\n[1/6] Loading HVG reference...")
ref = sc.read_h5ad(str(REF_HVG))
print(f"  Reference: {ref.shape}")

# 2. Attach human cluster labels
print("\n[2/6] Attaching cluster labels...")
human_full = sc.read_h5ad(str(HUMAN_REF), backed="r")
human_labels = human_full.obs[["supercluster_term", "cluster_id"]].copy()
human_full.file.close()
del human_full

ref.obs["supercluster_term"] = pd.Series(np.nan, index=ref.obs.index, dtype=object)
ref.obs["cluster_id"] = pd.Series(np.nan, index=ref.obs.index, dtype=object)

human_mask = ref.obs["species"] == "human"
human_idx = ref.obs.index[human_mask]
idx_stripped = human_idx.str.replace(r"-\d+$", "", regex=True)
idx_map = pd.Series(human_idx, index=idx_stripped)
shared = idx_stripped.intersection(human_labels.index)
ref.obs.loc[idx_map[shared].values, "supercluster_term"] = human_labels.loc[shared, "supercluster_term"].values
ref.obs.loc[idx_map[shared].values, "cluster_id"] = human_labels.loc[shared, "cluster_id"].values

# 3. Filter and subsample
print("\n[3/6] Filtering and subsampling...")
human_cge_mask = (ref.obs["species"] == "human") & (
    ref.obs["supercluster_term"] == "CGE interneuron"
)
mouse_ref_mask = ref.obs["species"] == "mouse"

np.random.seed(SEED)
human_cge_idx = ref.obs.index[human_cge_mask]
if len(human_cge_idx) > N_HUMAN_CGE:
    sub_obs = ref.obs.loc[human_cge_idx].copy()
    sub_obs = sub_obs[sub_obs["cluster_id"].notna()]
    sampled_parts = []
    for cid, grp in sub_obs.groupby("cluster_id"):
        n_take = max(1, int(N_HUMAN_CGE * len(grp) / len(sub_obs)))
        n_take = min(n_take, len(grp))
        sampled_parts.append(grp.sample(n=n_take, random_state=SEED))
    sampled = pd.concat(sampled_parts)
    if len(sampled) > N_HUMAN_CGE:
        sampled = sampled.sample(n=N_HUMAN_CGE, random_state=SEED)
    human_cge_idx = sampled.index

mouse_idx = ref.obs.index[mouse_ref_mask]
if len(mouse_idx) > N_MOUSE_REF:
    mouse_idx = pd.Index(
        np.random.RandomState(SEED).choice(mouse_idx, N_MOUSE_REF, replace=False)
    )
print(f"  Subsampled: human CGE={len(human_cge_idx)}, mouse ref={len(mouse_idx)}")

keep_idx = human_cge_idx.union(mouse_idx)
ref_sub = ref[keep_idx].copy()
del ref
gc.collect()

# 4. Load query and unify species, filter to M1 + CGE only
print("\n[4/6] Loading query...")
query = sc.read_h5ad(str(QUERY))
print(f"  Query (V1+M1, all): {query.shape}, CCKBCs: {int(query.obs['is_cckbc'].sum())}")

# Filter: M1 + CGE families only
CGE_FAMILIES = ["Sncg", "Vip", "Lamp5"]
m1_mask = query.obs["dataset"].astype(str) == "M1"
cge_mask = query.obs["RNA family"].astype(str).isin(CGE_FAMILIES)
keep_mask = m1_mask & cge_mask
print(f"  M1+CGE filter: {keep_mask.sum()} cells "
      f"(of which {int(query.obs.loc[keep_mask, 'is_cckbc'].sum())} CCKBC)")
query = query[keep_mask].copy()
query.obs["RNA_family"] = query.obs["RNA family"].astype(str)

query.obs["original_species"] = query.obs["species"].astype(str)
query.obs["species"] = "mouse"  # Unify with mouse reference for batch correction
query.obs["data_type"] = "query"

shared_genes = ref_sub.var_names.intersection(query.var_names)
query = query[:, shared_genes].copy()
ref_sub = ref_sub[:, shared_genes].copy()

ref_sub.obs["data_type"] = "reference"
if "is_cckbc" not in ref_sub.obs.columns:
    ref_sub.obs["is_cckbc"] = False

# Save metadata before concat
ref_obs_saved = ref_sub.obs.copy()
query_obs_saved = query.obs.copy()

# 5. Concatenate, normalize, PCA, Harmony, UMAP
print("\n[5/6] Concatenate + Harmony + UMAP...")
combined = sc.concat([ref_sub, query], join="inner")
combined.var_names_make_unique()

# Re-attach metadata after concat
all_obs_cols = set(ref_obs_saved.columns) | set(query_obs_saved.columns)
combined_obs_idx = combined.obs_names
for col in ["species", "supercluster_term", "cluster_id", "is_cckbc", "data_type", "RNA_family"]:
    if col in ref_obs_saved.columns or col in query_obs_saved.columns:
        ref_vals = ref_obs_saved[col] if col in ref_obs_saved.columns else pd.Series([np.nan]*len(ref_obs_saved), index=ref_obs_saved.index)
        q_vals = query_obs_saved[col] if col in query_obs_saved.columns else pd.Series([np.nan]*len(query_obs_saved), index=query_obs_saved.index)
        merged = pd.concat([ref_vals, q_vals])
        combined.obs[col] = merged.loc[combined_obs_idx].values

print(f"  Combined: {combined.shape}")

print("  Normalizing + PCA...")
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, n_comps=50, svd_solver="arpack")

print(f"  Running Harmony (batch=species, theta={THETA})...")
ho = hm.run_harmony(
    combined.obsm["X_pca"],
    combined.obs,
    "species",
    theta=THETA,
    max_iter_harmony=30,
    random_state=SEED,
)
Z = ho.Z_corr
if hasattr(Z, "cpu"):
    Z = Z.cpu().numpy()
Z = np.array(Z)
if Z.shape[0] == 50 and Z.shape[1] == combined.n_obs:
    Z = Z.T
combined.obsm["X_harmony"] = Z
print(f"  Harmony: {Z.shape}")

print("  Computing UMAP from Harmony embedding...")
sc.pp.neighbors(combined, use_rep="X_harmony", n_neighbors=30)
sc.tl.umap(combined, min_dist=0.3, spread=1.0, random_state=SEED)
print(f"  UMAP: {combined.obsm['X_umap'].shape}")

# 6. Save
print(f"\n[6/6] Saving to {OUT_PATH.name}...")
n = combined.n_obs
empty_X = sp.csr_matrix((n, 1), dtype=np.float32)

obs_clean = combined.obs.copy()
obs_clean["source"] = ["reference" if dt == "reference" else "query"
                       for dt in obs_clean["data_type"]]
for col in obs_clean.columns:
    if obs_clean[col].dtype == object:
        obs_clean[col] = obs_clean[col].astype(str).replace("nan", "")
    elif hasattr(obs_clean[col], "cat"):
        obs_clean[col] = obs_clean[col].astype(str)

out = ad.AnnData(X=empty_X, obs=obs_clean, var=pd.DataFrame(index=["dummy"]))
out.obsm["X_umap"] = combined.obsm["X_umap"]
out.obsm["X_harmony"] = combined.obsm["X_harmony"]
out.write(OUT_PATH)
print(f"  Saved: {OUT_PATH}")
print(f"  Size: {OUT_PATH.stat().st_size / 1e6:.1f} MB")

# Quick check: where do CCKBCs land?
query_cckbc_mask = (out.obs["source"] == "query") & (out.obs["is_cckbc"] == "True")
print(f"\nMouse CCKBCs: {query_cckbc_mask.sum()}")
if query_cckbc_mask.sum() > 0:
    cckbc_umap = out.obsm["X_umap"][query_cckbc_mask.values]
    print(f"  UMAP1 range: {cckbc_umap[:, 0].min():.1f} to {cckbc_umap[:, 0].max():.1f}")
    print(f"  UMAP2 range: {cckbc_umap[:, 1].min():.1f} to {cckbc_umap[:, 1].max():.1f}")
