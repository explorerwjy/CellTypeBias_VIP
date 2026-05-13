#!/usr/bin/env python3
"""
Compute UMAP for the within-mouse validation that Sncg subclass = CCKBC.

Concatenates WMB-10Xv3 cortical GABAergic atlas cells (reference, subsampled)
with M1 patch-seq cells (query), runs PCA + Harmony (batch=source) + UMAP,
and saves a small h5ad with metadata + UMAP coordinates for downstream
plotting in Fig_Supp_CCKBC_Mapping.ipynb.

Output:
  cge_subtype/results/sncg_validation_umap.h5ad
"""

import gc
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import harmonypy as hm
import scipy.sparse as sp
import torch
import yaml

warnings.filterwarnings("ignore")

_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

PROJECT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
CONFIG = PROJECT_DIR / "configs" / "cross_species_config.yaml"
M1_META = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_metadata.csv")
M1_COUNTS = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_counts.csv.gz")
OUT_PATH = PROJECT_DIR / "results" / "sncg_validation_umap.h5ad"

# Subsampling targets
N_REF = 25000      # WMB GABAergic reference cells (cortex)
SEED = 42

# CCKBC type definitions (Gouwens 2020)
CCKBC_RNA_TYPES = {
    "Sncg Col14a1", "Sncg Slc17a8", "Sncg Calb1_1", "Sncg Calb1_2",
    "Sncg Npy2r",
    "Vip Sncg",
    "Vip Serpinf1_1", "Vip Serpinf1_2", "Vip Serpinf1_3",
}

# ---------------------------------------------------------------------------
# 1. Load WMB-10Xv3 GABAergic cortex reference (subsampled)
# ---------------------------------------------------------------------------
print("=" * 60)
print("Compute Sncg validation UMAP (within-mouse mapping)")
print("=" * 60)

with open(CONFIG) as f:
    config = yaml.safe_load(f)

print("\n[1/5] Loading WMB-10Xv3 metadata...")
meta_path = config["mouse_reference"]["metadata_csv"]
meta = pd.read_csv(meta_path, low_memory=False)
meta = meta[meta["dataset_label"] == "WMB-10Xv3"]
print(f"  WMB-10Xv3 cells: {len(meta):,}")

# Region filter (cortex)
region_prefixes = config["mouse_reference"]["region_filter"]["include_prefixes"]
roi_mask = pd.Series(False, index=meta.index)
for prefix in region_prefixes:
    roi_mask |= meta["region_of_interest_acronym"].astype(str).str.startswith(prefix)
meta = meta[roi_mask]
print(f"  After region filter: {len(meta):,}")

# GABAergic only (CGE + MGE)
gaba_classes = ["06 CTX-CGE GABA", "07 CTX-MGE GABA"]
meta = meta[meta["class"].isin(gaba_classes)]
print(f"  GABAergic: {len(meta):,}")

# Subsample stratified by subclass
np.random.seed(SEED)
if len(meta) > N_REF:
    meta = meta.groupby("subclass", group_keys=False).apply(
        lambda x: x.sample(n=max(1, int(N_REF * len(x) / len(meta))), random_state=SEED)
    )
    if len(meta) > N_REF:
        meta = meta.sample(n=N_REF, random_state=SEED)
print(f"  Subsampled: {len(meta):,}")
print(f"  Subclass distribution:")
for sc_name, cnt in meta["subclass"].value_counts().items():
    print(f"    {sc_name}: {cnt}")

# Load expression for these cells
expr_dir = Path(config["mouse_reference"]["expression_dir"])
expr_files = config["mouse_reference"]["expression_files"]
cell_labels = set(meta["cell_label"].tolist())

print(f"\n[2/5] Loading WMB expression...")
chunks = []
for fname in expr_files:
    fpath = expr_dir / fname
    if not fpath.exists():
        continue
    print(f"  Loading {fname} ...")
    a = sc.read_h5ad(fpath, backed="r")
    file_labels = a.obs_names.values if "cell_label" not in a.obs.columns else a.obs["cell_label"].values
    keep_idx = [i for i, c in enumerate(file_labels) if c in cell_labels]
    if not keep_idx:
        a.file.close()
        continue
    keep_idx = np.array(keep_idx)
    print(f"    {len(keep_idx)} matching cells")
    Xs = a.X[keep_idx]
    if not sp.issparse(Xs):
        Xs = sp.csr_matrix(Xs)
    var_df = a.var.copy()
    if "gene_symbol" in var_df.columns:
        gene_names = var_df["gene_symbol"].tolist()
    else:
        gene_names = var_df.index.tolist()
    chunk_obs = pd.DataFrame(index=[file_labels[i] for i in keep_idx])
    chunk = ad.AnnData(X=Xs.copy(), obs=chunk_obs, var=pd.DataFrame(index=gene_names))
    chunk.var_names_make_unique()
    chunks.append(chunk)
    a.file.close()
    del a
    gc.collect()

# Merge
common_genes = set(chunks[0].var_names)
for c in chunks[1:]:
    common_genes &= set(c.var_names)
common_genes = sorted(common_genes)
for i in range(len(chunks)):
    chunks[i] = chunks[i][:, common_genes].copy()
adata_ref = ad.concat(chunks, join="inner")
adata_ref.var_names_make_unique()
del chunks
gc.collect()

# Attach metadata
meta_idx = meta.set_index("cell_label", drop=False)
for col in ["class", "subclass"]:
    adata_ref.obs[col] = meta_idx.loc[adata_ref.obs_names, col].values
adata_ref.obs["source"] = "atlas"
adata_ref.obs["batch"] = "atlas"
adata_ref.layers["counts"] = adata_ref.X.copy()
print(f"  Reference: {adata_ref.shape}")

# ---------------------------------------------------------------------------
# 3. Load M1 patch-seq query
# ---------------------------------------------------------------------------
print(f"\n[3/5] Loading M1 patch-seq query...")
counts = pd.read_csv(M1_COUNTS, index_col=0, compression="gzip").T
m1_meta = pd.read_csv(M1_META, sep="\t").set_index("Cell")
common = sorted(set(counts.index) & set(m1_meta.index))
counts = counts.loc[common]
m1_meta = m1_meta.loc[common]
print(f"  M1 cells: {len(counts)}")

# Build query AnnData
X_q = sp.csr_matrix(counts.values.astype(np.float32))
adata_q = ad.AnnData(X=X_q,
                     obs=pd.DataFrame(index=counts.index),
                     var=pd.DataFrame(index=counts.columns))
adata_q.layers["counts"] = adata_q.X.copy()

# Tag CCKBC + RNA family
adata_q.obs["RNA_family"] = m1_meta["RNA family"].values
adata_q.obs["RNA_type"] = m1_meta["RNA type"].values
adata_q.obs["is_cckbc"] = (
    m1_meta["RNA type"].astype(str).isin(CCKBC_RNA_TYPES)
    | (m1_meta["RNA family"].astype(str) == "Sncg")
).values
adata_q.obs["source"] = "patchseq"
adata_q.obs["batch"] = "patchseq"
print(f"  Query: {adata_q.shape}, CCKBCs: {int(adata_q.obs['is_cckbc'].sum())}")

# ---------------------------------------------------------------------------
# 4. Concatenate, normalize, PCA, Harmony, UMAP
# ---------------------------------------------------------------------------
print(f"\n[4/5] Concatenate + Harmony + UMAP...")
shared_genes = sorted(set(adata_ref.var_names) & set(adata_q.var_names))
print(f"  Shared genes: {len(shared_genes)}")
adata_ref_s = adata_ref[:, shared_genes].copy()
adata_q_s = adata_q[:, shared_genes].copy()

# Save metadata
ref_obs_saved = adata_ref_s.obs.copy()
q_obs_saved = adata_q_s.obs.copy()

combined = sc.concat([adata_ref_s, adata_q_s], join="inner")
combined.var_names_make_unique()

# Re-attach metadata
for col in ["subclass", "class", "RNA_family", "RNA_type", "is_cckbc", "source", "batch"]:
    if col in ref_obs_saved.columns or col in q_obs_saved.columns:
        ref_vals = ref_obs_saved[col] if col in ref_obs_saved.columns else pd.Series([np.nan]*len(ref_obs_saved), index=ref_obs_saved.index)
        q_vals = q_obs_saved[col] if col in q_obs_saved.columns else pd.Series([np.nan]*len(q_obs_saved), index=q_obs_saved.index)
        merged = pd.concat([ref_vals, q_vals])
        combined.obs[col] = merged.loc[combined.obs_names].values

print(f"  Combined: {combined.shape}")

print("  Normalizing + PCA...")
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
# HVGs from atlas only
atlas_mask = combined.obs["source"] == "atlas"
atlas_tmp = combined[atlas_mask].copy()
sc.pp.highly_variable_genes(atlas_tmp, n_top_genes=3000, flavor="seurat_v3", layer=None, span=0.3)
hvgs = atlas_tmp.var_names[atlas_tmp.var.highly_variable].tolist()
del atlas_tmp; gc.collect()
print(f"  HVGs: {len(hvgs)}")
combined_hvg = combined[:, hvgs].copy()
sc.pp.scale(combined_hvg, max_value=10)
sc.tl.pca(combined_hvg, n_comps=50)
combined.obsm["X_pca"] = combined_hvg.obsm["X_pca"]
del combined_hvg; gc.collect()

print("  Harmony...")
ho = hm.run_harmony(combined.obsm["X_pca"], combined.obs, "batch", max_iter_harmony=30)
Z = ho.Z_corr
if hasattr(Z, "cpu"):
    Z = Z.cpu().numpy()
Z = np.array(Z)
if Z.shape[0] == 50 and Z.shape[1] == combined.n_obs:
    Z = Z.T
combined.obsm["X_harmony"] = Z
print(f"  Harmony: {Z.shape}")

print("  UMAP from Harmony...")
sc.pp.neighbors(combined, use_rep="X_harmony", n_neighbors=30)
sc.tl.umap(combined, min_dist=0.3, spread=1.0, random_state=SEED)

# ---------------------------------------------------------------------------
# 5. Save (small file: only obs + obsm['X_umap'])
# ---------------------------------------------------------------------------
print(f"\n[5/5] Saving to {OUT_PATH.name}...")
n = combined.n_obs
empty_X = sp.csr_matrix((n, 1), dtype=np.float32)

obs_clean = combined.obs.copy()
for col in obs_clean.columns:
    if obs_clean[col].dtype == object:
        obs_clean[col] = obs_clean[col].astype(str).replace("nan", "")
    elif hasattr(obs_clean[col], "cat"):
        obs_clean[col] = obs_clean[col].astype(str)

out = ad.AnnData(X=empty_X, obs=obs_clean, var=pd.DataFrame(index=["dummy"]))
out.obsm["X_umap"] = combined.obsm["X_umap"]
out.write(OUT_PATH)
print(f"  Saved: {OUT_PATH}")
print(f"  Size: {OUT_PATH.stat().st_size / 1e6:.1f} MB")

print("\nQuery summary:")
q_obs_final = out.obs[out.obs["source"] == "patchseq"]
print(f"  Total query: {len(q_obs_final)}")
print(f"  CCKBC (is_cckbc=True): {(q_obs_final['is_cckbc'].astype(str) == 'True').sum()}")
print(f"  By RNA family: {q_obs_final['RNA_family'].value_counts().to_dict()}")
