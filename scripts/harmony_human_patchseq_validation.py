"""
Validate human patch-seq cluster assignments by mapping to Siletti atlas
using Harmony batch correction + k-NN classification.

Compare results with existing scVI-based mapping.
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
from scipy.sparse import issparse
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
PATCHSEQ_EXPR = "/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/expression_by_cell.csv.gz"
PATCHSEQ_META = "/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/cell_metadata.csv"
ANNOTATION = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/annotation.xlsx"
SCVI_RESULTS = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/results/human_scvi_mapping_results.csv"
REF_DIR = "/mnt/data0/HumanBrainCellType/SuperTypeRawDat"
OUTFILE = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/results/harmony_human_patchseq_validation.csv"

# Reference superclusters to load (the ones that patch-seq cells map to)
REF_SUPERCLUSTERS = [
    "CGE_interneuron",
    "MGE_interneuron",
    "LAMP5_LHX6_and_Chandelier",
    "Microglia",
    "Astrocyte",
]

# ── Parameters ─────────────────────────────────────────────────────────────
N_HVG = 3000          # Number of highly variable genes to use
N_PCS = 50            # Number of PCs for Harmony
K_NEIGHBORS = 15      # k for k-NN classification
MAX_REF_PER_SC = 5000 # Subsample each supercluster reference to this max
RANDOM_STATE = 42

np.random.seed(RANDOM_STATE)

# ── 1. Load annotation ────────────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Loading annotation")
print("=" * 70)
ann = pd.read_excel(ANNOTATION)
ann = ann.dropna(subset=["Cluster"])
ann["Cluster"] = ann["Cluster"].astype(int)
cluster_to_sc = dict(zip(ann["Cluster"], ann["Supercluster"]))
cluster_to_subtype = dict(zip(ann["Cluster"], ann["Subtype auto-annotation"]))
print(f"  Loaded annotation for {len(ann)} clusters")

# ── 2. Load patch-seq expression ──────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Loading patch-seq expression data")
print("=" * 70)
patchseq_expr = pd.read_csv(PATCHSEQ_EXPR, index_col=0)
patchseq_meta = pd.read_csv(PATCHSEQ_META, index_col=0)
print(f"  Patch-seq: {patchseq_expr.shape[0]} cells x {patchseq_expr.shape[1]} genes")
print(f"  Metadata: {patchseq_meta.shape[0]} cells")

# Gene names from patch-seq (gene symbols)
patchseq_genes = set(patchseq_expr.columns)

# ── 3. Load reference data (subsample for efficiency) ────────────────────
print("\n" + "=" * 70)
print("STEP 3: Loading Siletti reference data")
print("=" * 70)

ref_dfs = []
ref_meta_list = []

for sc_name in REF_SUPERCLUSTERS:
    h5ad_path = os.path.join(REF_DIR, f"Supercluster_{sc_name}.h5ad")
    if not os.path.exists(h5ad_path):
        print(f"  WARNING: {h5ad_path} not found, skipping")
        continue

    print(f"  Loading {sc_name}...", end=" ", flush=True)
    adata = sc.read_h5ad(h5ad_path)
    print(f"{adata.shape[0]} cells x {adata.shape[1]} genes")

    # Build gene symbol mapping from var
    gene_map = adata.var["Gene"].to_dict()  # Ensembl -> Symbol

    # Find shared genes with patch-seq
    ref_symbols = set(adata.var["Gene"].values)
    shared = patchseq_genes & ref_symbols

    # Subsample if too large
    if adata.shape[0] > MAX_REF_PER_SC:
        idx = np.random.choice(adata.shape[0], MAX_REF_PER_SC, replace=False)
        adata = adata[idx].copy()
        print(f"    Subsampled to {MAX_REF_PER_SC} cells")

    # Extract expression matrix (convert to dense if sparse)
    X = adata.X
    if issparse(X):
        X = X.toarray()

    # Create DataFrame with gene symbols as columns
    gene_symbols = adata.var["Gene"].values
    df = pd.DataFrame(X, index=adata.obs.index, columns=gene_symbols)

    # Handle duplicate gene symbols (keep first occurrence)
    df = df.loc[:, ~df.columns.duplicated()]

    # Store metadata
    meta = pd.DataFrame({
        "cluster_id": adata.obs["cluster_id"].values,
        "supercluster": sc_name.replace("_", " ").replace("LAMP5 LHX6 and Chandelier", "LAMP5-LHX6 and Chandelier"),
        "source": "reference",
    }, index=adata.obs.index)

    ref_dfs.append(df)
    ref_meta_list.append(meta)

    del adata, X

# Fix supercluster naming
for meta in ref_meta_list:
    meta["supercluster"] = meta["supercluster"].replace({
        "CGE interneuron": "CGE interneuron",
        "MGE interneuron": "MGE interneuron",
        "LAMP5 LHX6 and Chandelier": "LAMP5-LHX6 and Chandelier",
    })

ref_expr = pd.concat(ref_dfs, axis=0)
ref_meta = pd.concat(ref_meta_list, axis=0)
print(f"\n  Combined reference: {ref_expr.shape[0]} cells x {ref_expr.shape[1]} genes")

# ── 4. Find shared genes and prepare combined matrix ─────────────────────
print("\n" + "=" * 70)
print("STEP 4: Finding shared genes and combining data")
print("=" * 70)

shared_genes = sorted(set(patchseq_expr.columns) & set(ref_expr.columns))
print(f"  Shared genes: {len(shared_genes)}")

# Subset both to shared genes
patchseq_sub = patchseq_expr[shared_genes].copy()
ref_sub = ref_expr[shared_genes].copy()

# Combine
combined = pd.concat([ref_sub, patchseq_sub], axis=0)
batch_labels = (
    ["reference"] * ref_sub.shape[0] + ["patchseq"] * patchseq_sub.shape[0]
)
print(f"  Combined matrix: {combined.shape[0]} cells x {combined.shape[1]} genes")
print(f"    Reference: {ref_sub.shape[0]} cells")
print(f"    Patch-seq: {patchseq_sub.shape[0]} cells")

# Free memory
del patchseq_expr, ref_expr, ref_sub, patchseq_sub, ref_dfs

# ── 5. Normalize and select HVGs ─────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Normalizing and selecting HVGs")
print("=" * 70)

# Create AnnData for scanpy processing
adata_combined = sc.AnnData(
    X=combined.values.astype(np.float32),
    obs=pd.DataFrame({"batch": batch_labels}, index=combined.index),
)
adata_combined.var_names = combined.columns

# The patch-seq data appears to be already log-transformed (values like 4.98, 7.28)
# The reference data is raw counts. We need to normalize the reference consistently.
# Strategy: normalize reference to CPM, log1p, then use the combined data.

# Check value ranges
ref_mask = adata_combined.obs["batch"] == "reference"
ps_mask = adata_combined.obs["batch"] == "patchseq"

ref_mean = adata_combined.X[ref_mask.values].mean()
ps_mean = adata_combined.X[ps_mask.values].mean()
ref_max = adata_combined.X[ref_mask.values].max()
ps_max = adata_combined.X[ps_mask.values].max()

print(f"  Reference - mean: {ref_mean:.3f}, max: {ref_max:.3f}")
print(f"  Patch-seq - mean: {ps_mean:.3f}, max: {ps_max:.3f}")

# Normalize reference and patchseq separately, then recombine
# Reference: total-count normalize + log1p
ref_adata = adata_combined[ref_mask].copy()
sc.pp.normalize_total(ref_adata, target_sum=1e4)
sc.pp.log1p(ref_adata)

# Patch-seq: already log-transformed (log2(CPM+1) style), convert to natural log scale
# The values suggest log2(CPM+1) since max ~16 which is log2(65536)
# Convert: ln(2^x - 1 + 1) = x * ln(2) for x >> 0, approximately
# Actually, let's normalize patch-seq too for consistency
ps_adata = adata_combined[ps_mask].copy()

# Check if patch-seq is already log-scale by looking at distribution
ps_nonzero = ps_adata.X[ps_adata.X > 0]
print(f"  Patch-seq nonzero values - median: {np.median(ps_nonzero):.3f}, "
      f"max: {np.max(ps_nonzero):.3f}")

# Patch-seq appears log2(CPM+1) scaled. Convert to linear, then normalize like reference.
# 2^x - 1 to reverse log2(x+1)
print("  Converting patch-seq from log2(CPM+1) to counts, then re-normalizing...")
ps_linear = np.power(2, ps_adata.X) - 1
ps_adata.X = ps_linear.astype(np.float32)
sc.pp.normalize_total(ps_adata, target_sum=1e4)
sc.pp.log1p(ps_adata)

# Recombine
adata_combined = sc.concat([ref_adata, ps_adata], axis=0)
adata_combined.obs["batch"] = batch_labels

del ref_adata, ps_adata

# Select highly variable genes
sc.pp.highly_variable_genes(
    adata_combined, n_top_genes=N_HVG, batch_key="batch", flavor="seurat_v3"
)
n_hvg = adata_combined.var["highly_variable"].sum()
print(f"  Selected {n_hvg} highly variable genes")

# Subset to HVGs
adata_hvg = adata_combined[:, adata_combined.var["highly_variable"]].copy()

# Scale
sc.pp.scale(adata_hvg, max_value=10)

# PCA
sc.tl.pca(adata_hvg, n_comps=N_PCS, random_state=RANDOM_STATE)
print(f"  PCA done: {N_PCS} components")

# ── 6. Run Harmony ───────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Running Harmony batch correction")
print("=" * 70)

pca_embeddings = adata_hvg.obsm["X_pca"]
meta_df = pd.DataFrame({"batch": batch_labels}, index=adata_hvg.obs.index)

ho = hm.run_harmony(
    pca_embeddings,
    meta_df,
    "batch",
    max_iter_harmony=30,
    random_state=RANDOM_STATE,
)

# Z_corr may be a torch tensor on GPU; convert to numpy
Z_corr = ho.Z_corr
if hasattr(Z_corr, 'cpu'):
    Z_corr = Z_corr.cpu().numpy()
elif hasattr(Z_corr, 'values'):
    Z_corr = Z_corr.values

# Z_corr shape is (n_pcs, n_cells), transpose to (n_cells, n_pcs)
if Z_corr.shape[0] == N_PCS and Z_corr.shape[1] != N_PCS:
    harmony_embeddings = Z_corr.T
else:
    harmony_embeddings = Z_corr
print(f"  Harmony embeddings shape: {harmony_embeddings.shape}")

# Store in adata
adata_hvg.obsm["X_harmony"] = harmony_embeddings

# ── 7. k-NN classification ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 7: k-NN classification in Harmony space")
print("=" * 70)

ref_mask = np.array(batch_labels) == "reference"
ps_mask = np.array(batch_labels) == "patchseq"

X_ref = harmony_embeddings[ref_mask]
X_ps = harmony_embeddings[ps_mask]

# Reference labels
ref_cluster_ids = ref_meta["cluster_id"].values
ref_superclusters = ref_meta["supercluster"].values

# Train k-NN for supercluster assignment
print(f"  Training k-NN (k={K_NEIGHBORS}) for supercluster...")
knn_sc = KNeighborsClassifier(n_neighbors=K_NEIGHBORS, weights="distance", metric="euclidean")
knn_sc.fit(X_ref, ref_superclusters)
ps_supercluster = knn_sc.predict(X_ps)
ps_supercluster_proba = knn_sc.predict_proba(X_ps)
ps_supercluster_conf = ps_supercluster_proba.max(axis=1)

# Train k-NN for cluster assignment
print(f"  Training k-NN (k={K_NEIGHBORS}) for cluster...")
knn_cl = KNeighborsClassifier(n_neighbors=K_NEIGHBORS, weights="distance", metric="euclidean")
knn_cl.fit(X_ref, ref_cluster_ids)
ps_cluster = knn_cl.predict(X_ps)
ps_cluster_proba = knn_cl.predict_proba(X_ps)
ps_cluster_conf = ps_cluster_proba.max(axis=1)

# Map cluster -> supercluster for validation
ps_cluster_supercluster = [cluster_to_sc.get(int(c), "Unknown") for c in ps_cluster]

print(f"  Classified {X_ps.shape[0]} patch-seq cells")

# ── 8. Build results DataFrame ──────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 8: Building results")
print("=" * 70)

# Get patch-seq cell indices
ps_indices = combined.index[len(ref_meta):]

results = pd.DataFrame({
    "harmony_supercluster": ps_supercluster,
    "harmony_supercluster_conf": ps_supercluster_conf,
    "harmony_cluster": ps_cluster,
    "harmony_cluster_conf": ps_cluster_conf,
    "harmony_cluster_supercluster": ps_cluster_supercluster,
}, index=ps_indices)

# Add subtype annotation
results["harmony_cluster_subtype"] = [
    cluster_to_subtype.get(int(c), "Unknown") for c in ps_cluster
]

print(f"  Results shape: {results.shape}")
print(f"\n  Harmony supercluster assignments:")
print(results["harmony_supercluster"].value_counts().to_string())

# ── 9. Compare with scVI results ────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 9: Comparing with scVI assignments")
print("=" * 70)

scvi = pd.read_csv(SCVI_RESULTS, index_col=0)
scvi.index = scvi.index.astype(str)
results.index = results.index.astype(str)

# Merge
merged = results.join(scvi[["assigned_supercluster", "conf_supercluster",
                             "assigned_cluster", "conf_cluster"]],
                       how="inner", rsuffix="_scvi")
merged.rename(columns={
    "assigned_supercluster": "scvi_supercluster",
    "conf_supercluster": "scvi_supercluster_conf",
    "assigned_cluster": "scvi_cluster",
    "conf_cluster": "scvi_cluster_conf",
}, inplace=True)

print(f"  Merged cells: {merged.shape[0]} (of {results.shape[0]} harmony, {scvi.shape[0]} scVI)")

# Overall supercluster concordance
sc_match = (merged["harmony_supercluster"] == merged["scvi_supercluster"]).sum()
sc_total = len(merged)
print(f"\n  OVERALL SUPERCLUSTER CONCORDANCE: {sc_match}/{sc_total} = {sc_match/sc_total:.1%}")

# Per-supercluster concordance
print("\n  Per-supercluster concordance:")
for sc_name in sorted(merged["scvi_supercluster"].unique()):
    mask = merged["scvi_supercluster"] == sc_name
    n = mask.sum()
    match = (merged.loc[mask, "harmony_supercluster"] == sc_name).sum()
    print(f"    {sc_name}: {match}/{n} ({match/n:.1%})")

# Confusion matrix at supercluster level
print("\n  Supercluster confusion matrix (rows=scVI, cols=Harmony):")
ct = pd.crosstab(merged["scvi_supercluster"], merged["harmony_supercluster"])
print(ct.to_string())

# ── 10. Focus on CGE interneurons ────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 10: CGE interneuron analysis")
print("=" * 70)

# Cells assigned to CGE by scVI
cge_scvi = merged[merged["scvi_supercluster"] == "CGE interneuron"]
print(f"\n  CGE cells (scVI): {len(cge_scvi)}")

# Of those, how many are also CGE by Harmony?
cge_harmony_match = (cge_scvi["harmony_supercluster"] == "CGE interneuron").sum()
print(f"  Also CGE by Harmony: {cge_harmony_match}/{len(cge_scvi)} = {cge_harmony_match/len(cge_scvi):.1%}")

# Where do the disagreements go?
cge_disagree = cge_scvi[cge_scvi["harmony_supercluster"] != "CGE interneuron"]
if len(cge_disagree) > 0:
    print(f"\n  Disagreements (scVI=CGE, Harmony!=CGE): {len(cge_disagree)}")
    print(cge_disagree["harmony_supercluster"].value_counts().to_string())

# Cells assigned to CGE by Harmony
cge_harmony = merged[merged["harmony_supercluster"] == "CGE interneuron"]
print(f"\n  CGE cells (Harmony): {len(cge_harmony)}")
harmony_only_cge = cge_harmony[cge_harmony["scvi_supercluster"] != "CGE interneuron"]
if len(harmony_only_cge) > 0:
    print(f"  Harmony=CGE but scVI!=CGE: {len(harmony_only_cge)}")
    print(harmony_only_cge["scvi_supercluster"].value_counts().to_string())

# Cluster-level concordance for CGE cells (where both agree on CGE)
cge_both = merged[
    (merged["scvi_supercluster"] == "CGE interneuron") &
    (merged["harmony_supercluster"] == "CGE interneuron")
]
print(f"\n  CGE cells agreed by both methods: {len(cge_both)}")

if len(cge_both) > 0:
    # Check cluster-level agreement
    # scVI cluster is a cluster index, harmony cluster is also a cluster index
    cl_match = (cge_both["harmony_cluster"].astype(int) == cge_both["scvi_cluster"].astype(int)).sum()
    print(f"  Exact cluster match: {cl_match}/{len(cge_both)} = {cl_match/len(cge_both):.1%}")

    # Show cluster assignment comparison
    print("\n  Cluster assignment comparison (scVI vs Harmony) for CGE cells:")
    ct_cluster = pd.crosstab(
        cge_both["scvi_cluster"].astype(int),
        cge_both["harmony_cluster"].astype(int),
    )
    # Only show non-zero entries
    ct_cluster = ct_cluster.loc[ct_cluster.sum(axis=1) > 0, ct_cluster.sum(axis=0) > 0]
    print(ct_cluster.to_string())

    # Subtype-level concordance
    print("\n  Harmony cluster subtypes for CGE cells:")
    print(cge_both["harmony_cluster_subtype"].value_counts().to_string())

# ── 11. Confidence analysis ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 11: Confidence analysis")
print("=" * 70)

# Compare confidence distributions
print("  Harmony supercluster confidence:")
print(f"    Mean: {merged['harmony_supercluster_conf'].mean():.3f}")
print(f"    Median: {merged['harmony_supercluster_conf'].median():.3f}")
print(f"    >0.5: {(merged['harmony_supercluster_conf'] > 0.5).sum()}/{len(merged)}")
print(f"    >0.8: {(merged['harmony_supercluster_conf'] > 0.8).sum()}/{len(merged)}")

print("\n  High-confidence concordance (Harmony conf > 0.8):")
hc = merged[merged["harmony_supercluster_conf"] > 0.8]
if len(hc) > 0:
    hc_match = (hc["harmony_supercluster"] == hc["scvi_supercluster"]).sum()
    print(f"    {hc_match}/{len(hc)} = {hc_match/len(hc):.1%}")

print("\n  High-confidence concordance (both methods conf > 0.5):")
hc2 = merged[
    (merged["harmony_supercluster_conf"] > 0.5) &
    (merged["scvi_supercluster_conf"] > 0.5)
]
if len(hc2) > 0:
    hc2_match = (hc2["harmony_supercluster"] == hc2["scvi_supercluster"]).sum()
    print(f"    {hc2_match}/{len(hc2)} = {hc2_match/len(hc2):.1%}")

# ── 12. Save results ────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 12: Saving results")
print("=" * 70)

# Save full merged results
merged.to_csv(OUTFILE)
print(f"  Saved to: {OUTFILE}")
print(f"  Shape: {merged.shape}")
print(f"  Columns: {merged.columns.tolist()}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
