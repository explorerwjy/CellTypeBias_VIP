#!/usr/bin/env python3
"""
Compute joint UMAP for cross-species CGE visualization using pre-computed
scVI latent embeddings (reference + scArches-mapped query).

Loads:
  - reference_with_latent.h5ad  (227k cells, scVI X_latent)
  - patchseq_query_mapped.h5ad  (5764 cells, scArches X_latent in same space)

Filters reference to human CGE + mouse cells, subsamples, concatenates latents
with the query, then computes a joint UMAP for visualization.

Output:
  results/harmony_umap_cge.h5ad  (small: only obs + obsm['X_umap'])
"""

import gc
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
import scipy.sparse as sp

warnings.filterwarnings("ignore")

# torch.load patch
_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

BASE = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
REF_LATENT = BASE / "results/cross_species_preprocessed/reference_with_latent.h5ad"
QUERY_MAPPED = Path("/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/cross_species/mapped/patchseq_query_mapped.h5ad")
QUERY_RAW = BASE / "results/cross_species_preprocessed/mouse_patchseq_query.h5ad"
OUT_PATH = BASE / "results" / "harmony_umap_cge.h5ad"

N_HUMAN_CGE = 25000
N_MOUSE_REF = 8000
SEED = 42

print("=" * 60)
print("Joint UMAP from scVI latents (cross-species CGE)")
print("=" * 60)

# ---------------------------------------------------------------------------
# 1. Load reference (with scVI latent)
# ---------------------------------------------------------------------------
print(f"\n[1/5] Loading reference latent from {REF_LATENT.name}...")
ref = sc.read_h5ad(str(REF_LATENT))
print(f"  Reference: {ref.shape}")
print(f"  Obsm keys: {list(ref.obsm.keys())}")
print(f"  Species: {ref.obs['species'].value_counts().to_dict()}")

# Filter: human CGE + all mouse
human_cge_mask = (ref.obs["species"] == "human") & (
    ref.obs["supercluster_term"] == "CGE interneuron"
)
mouse_mask = ref.obs["species"] == "mouse"
print(f"  Human CGE: {human_cge_mask.sum()}, Mouse ref: {mouse_mask.sum()}")

# Subsample
np.random.seed(SEED)
human_idx = ref.obs.index[human_cge_mask]
if len(human_idx) > N_HUMAN_CGE:
    sub_obs = ref.obs.loc[human_idx].copy()
    # Drop rows with NaN cluster_id, fill empty categories
    sub_obs = sub_obs[sub_obs["cluster_id"].notna()]
    sampled_parts = []
    for cid, grp in sub_obs.groupby("cluster_id"):
        n_take = max(1, int(N_HUMAN_CGE * len(grp) / len(sub_obs)))
        n_take = min(n_take, len(grp))
        sampled_parts.append(grp.sample(n=n_take, random_state=SEED))
    sampled = pd.concat(sampled_parts)
    if len(sampled) > N_HUMAN_CGE:
        sampled = sampled.sample(n=N_HUMAN_CGE, random_state=SEED)
    human_idx = sampled.index

mouse_idx = ref.obs.index[mouse_mask]
if len(mouse_idx) > N_MOUSE_REF:
    mouse_idx = pd.Index(
        np.random.RandomState(SEED).choice(mouse_idx, N_MOUSE_REF, replace=False)
    )

print(f"  Subsampled: human CGE={len(human_idx)}, mouse ref={len(mouse_idx)}")

keep_idx = human_idx.union(mouse_idx)
ref_sub = ref[keep_idx].copy()
ref_latent = ref_sub.obsm["X_latent"].copy()
del ref
gc.collect()
print(f"  Subset reference latent: {ref_latent.shape}")

# ---------------------------------------------------------------------------
# 2. Load query (with scArches latent), filter to M1 + CGE only
# ---------------------------------------------------------------------------
print(f"\n[2/5] Loading query mapped latent from {QUERY_MAPPED.name}...")
query = sc.read_h5ad(str(QUERY_MAPPED))
print(f"  Query (V1+M1, all): {query.shape}")
print(f"  Obsm keys: {list(query.obsm.keys())}")

# Need is_cckbc + RNA family + dataset, all in the raw query file
print(f"\n  Loading is_cckbc + metadata from raw query file...")
query_raw = sc.read_h5ad(str(QUERY_RAW), backed="r")
qmeta = query_raw.obs.loc[query.obs_names, ["is_cckbc", "dataset", "RNA family"]].copy()
query_raw.file.close()

# Filter: M1 dataset + CGE families only (Sncg + Vip + Lamp5)
CGE_FAMILIES = ["Sncg", "Vip", "Lamp5"]
m1_mask = qmeta["dataset"].astype(str) == "M1"
cge_mask = qmeta["RNA family"].astype(str).isin(CGE_FAMILIES)
keep_mask = m1_mask & cge_mask

print(f"  M1 cells: {m1_mask.sum()}")
print(f"  CGE-family (Sncg+Vip+Lamp5) cells: {cge_mask.sum()}")
print(f"  M1 + CGE (intersection): {keep_mask.sum()}")

query = query[keep_mask.values].copy()
query.obs["is_cckbc"] = qmeta.loc[query.obs_names, "is_cckbc"].values
query.obs["RNA_family"] = qmeta.loc[query.obs_names, "RNA family"].values
query_latent = query.obsm["X_latent"].copy()

print(f"  Filtered query: {query.shape}")
print(f"  CCKBCs in filtered set: {int(query.obs['is_cckbc'].sum())}")
print(f"  Non-CCKBC by family: {query.obs.loc[~query.obs['is_cckbc'].astype(bool), 'RNA_family'].value_counts().to_dict()}")

# ---------------------------------------------------------------------------
# 3. Combine latents and compute joint UMAP
# ---------------------------------------------------------------------------
print(f"\n[3/5] Building combined AnnData with joint latents...")

n_ref = ref_sub.n_obs
n_query = query.n_obs
n_total = n_ref + n_query

combined_latent = np.vstack([ref_latent, query_latent])
print(f"  Combined latent: {combined_latent.shape}")

# Build a minimal AnnData for UMAP
empty_X = sp.csr_matrix((n_total, 1), dtype=np.float32)
combined_obs = pd.DataFrame(
    index=list(ref_sub.obs_names) + list(query.obs_names)
)
combined_obs["source"] = ["reference"] * n_ref + ["query"] * n_query
combined_obs["species"] = (
    ref_sub.obs["species"].astype(str).tolist()
    + ["mouse_patchseq"] * n_query
)
combined_obs["supercluster_term"] = (
    ref_sub.obs["supercluster_term"].astype(str).tolist()
    + [""] * n_query
)
combined_obs["cluster_id"] = (
    ref_sub.obs["cluster_id"].astype(str).tolist()
    + [""] * n_query
)
combined_obs["is_cckbc"] = (
    [False] * n_ref + query.obs["is_cckbc"].astype(bool).tolist()
)
combined_obs["is_cckbc"] = combined_obs["is_cckbc"].astype(str)
combined_obs["RNA_family"] = (
    [""] * n_ref + query.obs["RNA_family"].astype(str).tolist()
)

combined = ad.AnnData(X=empty_X, obs=combined_obs, var=pd.DataFrame(index=["dummy"]))
combined.obsm["X_latent"] = combined_latent

# ---------------------------------------------------------------------------
# 4. Compute joint UMAP from scVI latent
# ---------------------------------------------------------------------------
print(f"\n[4/5] Computing UMAP from scVI latent ({n_total} cells)...")
sc.pp.neighbors(combined, use_rep="X_latent", n_neighbors=30)
sc.tl.umap(combined, min_dist=0.3, spread=1.0, random_state=SEED)
print(f"  UMAP: {combined.obsm['X_umap'].shape}")

# ---------------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------------
print(f"\n[5/5] Saving to {OUT_PATH.name}...")

# Clean object dtypes for h5ad write
for col in combined.obs.columns:
    if combined.obs[col].dtype == object:
        combined.obs[col] = combined.obs[col].astype(str).replace("nan", "")

combined.write(OUT_PATH)
print(f"  Saved: {OUT_PATH}")
print(f"  Size: {OUT_PATH.stat().st_size / 1e6:.1f} MB")

# Quick sanity check: where do CCKBCs land?
query_cckbc_mask = (combined.obs["source"] == "query") & (combined.obs["is_cckbc"] == "True")
print(f"\nMouse CCKBCs in UMAP: {query_cckbc_mask.sum()}")
if query_cckbc_mask.sum() > 0:
    cckbc_umap = combined.obsm["X_umap"][query_cckbc_mask.values]
    print(f"  UMAP1 range: {cckbc_umap[:, 0].min():.1f} to {cckbc_umap[:, 0].max():.1f}")
    print(f"  UMAP2 range: {cckbc_umap[:, 1].min():.1f} to {cckbc_umap[:, 1].max():.1f}")
