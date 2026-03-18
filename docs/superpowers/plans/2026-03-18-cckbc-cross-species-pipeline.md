# Cross-Species Multi-Modal CCKBC/ISI VIP Pipeline — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a 5-module pipeline that integrates 4 datasets (mouse SC, mouse patch-seq, human SC, human patch-seq) to robustify CCKBC vs ISI VIP 22q bias analysis for reviewer response.

**Architecture:** Cluster-level cross-species bridge (Spearman RBH + MetaNeighbor AUROC), improved mouse patch-seq mapping (Harmony + hyperparameter sweep), path concordance validation, cross-species ephys convergence (same EphysSumStats pipeline), and updated 22q bias on multi-modal CCKBC groupings.

**Tech Stack:** Python 3.x, scanpy, anndata, harmonypy, scipy, sklearn, matplotlib, pandas, numpy. Conda env: `gencic`.

**Spec:** `docs/superpowers/specs/2026-03-17-cckbc-cross-species-pipeline-design.md`

**Prerequisite check before starting:**
```bash
conda activate gencic
python -c "import harmonypy; print('harmonypy OK')" || pip install harmonypy
python -c "import scanpy, anndata, scipy, sklearn; print('core deps OK')"
mkdir -p /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/tests
```

---

## File Structure

```
cge_subtype/
├── src/
│   ├── __init__.py                      # existing
│   ├── scvi_integration.py              # existing (unchanged)
│   ├── cluster_correspondence.py        # NEW: Module A — pseudobulk, Spearman RBH, MetaNeighbor AUROC
│   ├── harmony_mapping.py               # NEW: Module B — subsample, Harmony, kNN transfer, evaluation
│   ├── concordance.py                   # NEW: Module C — path chaining, concordance metrics
│   ├── ephys_harmonization.py           # NEW: Module D — aggregation, z-scoring, permutation test
│   └── multimodal_classification.py     # NEW: Module E — multi-modal confidence scoring
├── scripts/
│   ├── 06_compute_pseudobulk.py         # NEW: Module A — RAM-safe pseudobulk for both atlases
│   ├── 07_cluster_bridge.py             # NEW: Module A — Spearman RBH (Pass 1)
│   ├── 07b_metaneighbor_pass2.py        # NEW: Module A — MetaNeighbor AUROC on unresolved clusters (Pass 2)
│   ├── 07c_subsample_reference.py       # NEW: Module A/B — Create subsampled reference h5ad
│   ├── 08_harmony_mapping.py            # NEW: Module B — Harmony + kNN + evaluation
│   ├── 09_harmony_sweep.py              # NEW: Module B — hyperparameter sweep runner
│   └── 10_aggregate_ephys.py            # NEW: Module D — aggregate EphysSumStats analysis.csv
├── notebooks/
│   ├── Path_Concordance.py              # NEW: Module C — concordance analysis + figures
│   ├── Ephys_Convergence.py             # NEW: Module D — ephys convergence analysis + figures
│   └── Updated_22q_Bias.py              # NEW: Module E — 22q bias with multi-modal groupings
├── tests/
│   ├── __init__.py                      # NEW
│   ├── test_cluster_correspondence.py   # NEW
│   ├── test_concordance.py              # NEW
│   ├── test_ephys_harmonization.py      # NEW
│   └── test_multimodal_classification.py # NEW
└── configs/
    └── cross_species_config.yaml        # existing (may add new keys)
```

**Dependency order:**
- Tasks 1-2 (Module A library + scripts): independent, can start immediately
- Tasks 3-4 (Module B library + scripts): independent, can run in parallel with A
- Task 5 (Module C): depends on A + B outputs
- Task 6 (Module D): data prep independent; analysis depends on existing mappings
- Task 7 (Module E): depends on all previous

---

## Chunk 1: Module A — Cross-Species Cluster Bridge

### Task 1: Core Library — cluster_correspondence.py

**Files:**
- Create: `cge_subtype/src/cluster_correspondence.py`
- Create: `cge_subtype/tests/__init__.py`
- Create: `cge_subtype/tests/test_cluster_correspondence.py`

- [ ] **Step 1: Write tests for pseudobulk computation**

```python
# cge_subtype/tests/test_cluster_correspondence.py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cge_subtype.src.cluster_correspondence import compute_pseudobulk


def test_pseudobulk_dense():
    """Pseudobulk = mean expression per cluster."""
    # 6 cells, 4 genes, 2 clusters
    expr = pd.DataFrame(
        [[1, 2, 3, 4], [3, 4, 5, 6], [5, 6, 7, 8],    # cluster A
         [10, 20, 30, 40], [30, 40, 50, 60], [50, 60, 70, 80]],  # cluster B
        columns=["g1", "g2", "g3", "g4"],
    )
    labels = pd.Series(["A", "A", "A", "B", "B", "B"])
    result = compute_pseudobulk(expr, labels)
    assert result.shape == (2, 4)
    np.testing.assert_array_almost_equal(result.loc["A"], [3, 4, 5, 6])
    np.testing.assert_array_almost_equal(result.loc["B"], [30, 40, 50, 60])


def test_pseudobulk_single_cell_cluster():
    """Cluster with one cell returns that cell's expression."""
    expr = pd.DataFrame([[1, 2], [3, 4]], columns=["g1", "g2"])
    labels = pd.Series(["A", "B"])
    result = compute_pseudobulk(expr, labels)
    np.testing.assert_array_almost_equal(result.loc["A"], [1, 2])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_cluster_correspondence.py -v`
Expected: FAIL (ImportError — module doesn't exist yet)

- [ ] **Step 3: Implement compute_pseudobulk**

```python
# cge_subtype/src/cluster_correspondence.py
"""Cross-species cluster correspondence: pseudobulk, Spearman RBH, MetaNeighbor AUROC."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


def compute_pseudobulk(expr_df, cluster_labels):
    """Compute mean expression per cluster (pseudobulk centroid).

    Parameters
    ----------
    expr_df : pd.DataFrame
        Cells x genes expression matrix.
    cluster_labels : pd.Series
        Cluster label per cell, same length as expr_df rows.

    Returns
    -------
    pd.DataFrame
        Clusters x genes mean expression matrix.
    """
    expr_df = expr_df.copy()
    expr_df["_cluster"] = cluster_labels.values
    return expr_df.groupby("_cluster").mean()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_cluster_correspondence.py::test_pseudobulk_dense cge_subtype/tests/test_cluster_correspondence.py::test_pseudobulk_single_cell_cluster -v`
Expected: PASS

- [ ] **Step 5: Write tests for Spearman RBH**

Add to `cge_subtype/tests/test_cluster_correspondence.py`:

```python
from cge_subtype.src.cluster_correspondence import (
    compute_spearman_corr_matrix,
    find_reciprocal_best_hits,
    determine_rbh_threshold,
)


def test_spearman_corr_matrix_shape():
    """Correlation matrix has shape (n_mouse, n_human)."""
    mouse = pd.DataFrame(np.random.rand(5, 100), index=[f"m{i}" for i in range(5)])
    human = pd.DataFrame(np.random.rand(8, 100), index=[f"h{i}" for i in range(8)])
    corr = compute_spearman_corr_matrix(mouse, human)
    assert corr.shape == (5, 8)
    assert corr.index.tolist() == [f"m{i}" for i in range(5)]
    assert corr.columns.tolist() == [f"h{i}" for i in range(8)]


def test_spearman_corr_identical():
    """Identical vectors should have correlation ~1.0."""
    shared = np.random.rand(100)
    mouse = pd.DataFrame([shared, np.random.rand(100)], index=["m0", "m1"])
    human = pd.DataFrame([shared, np.random.rand(100)], index=["h0", "h1"])
    corr = compute_spearman_corr_matrix(mouse, human)
    assert corr.loc["m0", "h0"] > 0.99


def test_rbh_perfect_pairs():
    """Perfect reciprocal best hits with clear separation."""
    # 3 mouse, 3 human: m0<->h0, m1<->h1, m2<->h2
    corr = pd.DataFrame(
        [[0.95, 0.3, 0.2], [0.2, 0.90, 0.3], [0.1, 0.2, 0.85]],
        index=["m0", "m1", "m2"],
        columns=["h0", "h1", "h2"],
    )
    rbh = find_reciprocal_best_hits(corr, threshold=0.5)
    assert len(rbh) == 3
    assert rbh[rbh["mouse_cluster"] == "m0"]["human_cluster"].values[0] == "h0"
    assert rbh[rbh["mouse_cluster"] == "m1"]["human_cluster"].values[0] == "h1"


def test_rbh_threshold_filters():
    """RBH below threshold are excluded."""
    corr = pd.DataFrame(
        [[0.95, 0.3], [0.2, 0.45]],  # m1<->h1 is RBH but below threshold
        index=["m0", "m1"],
        columns=["h0", "h1"],
    )
    rbh = find_reciprocal_best_hits(corr, threshold=0.5)
    assert len(rbh) == 1  # only m0<->h0


def test_determine_threshold_permutation():
    """Permutation-based threshold should be > 0."""
    corr = pd.DataFrame(
        np.random.rand(10, 10),
        index=[f"m{i}" for i in range(10)],
        columns=[f"h{i}" for i in range(10)],
    )
    threshold = determine_rbh_threshold(corr, method="permutation", n_perm=100, seed=42)
    assert 0 < threshold < 1
```

- [ ] **Step 6: Implement Spearman RBH functions**

Add to `cge_subtype/src/cluster_correspondence.py`:

```python
def compute_spearman_corr_matrix(mouse_centroids, human_centroids):
    """Compute Spearman correlation between all mouse and human cluster centroids.

    Parameters
    ----------
    mouse_centroids : pd.DataFrame
        Mouse clusters x genes.
    human_centroids : pd.DataFrame
        Human clusters x genes.

    Returns
    -------
    pd.DataFrame
        Mouse clusters x human clusters correlation matrix.
    """
    shared_genes = mouse_centroids.columns.intersection(human_centroids.columns)
    mouse_mat = mouse_centroids[shared_genes].values
    human_mat = human_centroids[shared_genes].values

    n_mouse = mouse_mat.shape[0]
    n_human = human_mat.shape[0]
    corr_mat = np.zeros((n_mouse, n_human))

    for i in range(n_mouse):
        for j in range(n_human):
            corr_mat[i, j], _ = spearmanr(mouse_mat[i], human_mat[j])

    return pd.DataFrame(
        corr_mat,
        index=mouse_centroids.index,
        columns=human_centroids.index,
    )


def find_reciprocal_best_hits(corr_matrix, threshold=None):
    """Find reciprocal best hits from a correlation matrix.

    Parameters
    ----------
    corr_matrix : pd.DataFrame
        Mouse clusters (rows) x human clusters (columns).
    threshold : float or None
        Minimum correlation for a pair to be considered resolved.
        If None, all RBH pairs are returned regardless of correlation.

    Returns
    -------
    pd.DataFrame
        Columns: mouse_cluster, human_cluster, correlation, is_rbh.
    """
    mouse_best = corr_matrix.idxmax(axis=1)  # best human for each mouse
    human_best = corr_matrix.idxmax(axis=0)  # best mouse for each human

    pairs = []
    for m_clust in corr_matrix.index:
        h_best = mouse_best[m_clust]
        corr_val = corr_matrix.loc[m_clust, h_best]
        is_rbh = human_best[h_best] == m_clust
        if is_rbh:
            if threshold is None or corr_val >= threshold:
                pairs.append({
                    "mouse_cluster": m_clust,
                    "human_cluster": h_best,
                    "correlation": corr_val,
                    "is_rbh": True,
                    "method": "spearman_rbh",
                })
    return pd.DataFrame(pairs)


def determine_rbh_threshold(corr_matrix, method="permutation", n_perm=1000, seed=42):
    """Determine the RBH correlation threshold from data.

    Parameters
    ----------
    corr_matrix : pd.DataFrame
        Mouse x human correlation matrix.
    method : str
        "permutation" — 95th percentile of shuffled RBH correlations.
    n_perm : int
        Number of permutations.
    seed : int
        Random seed.

    Returns
    -------
    float
        Threshold value.
    """
    rng = np.random.default_rng(seed)
    null_corrs = []

    for _ in range(n_perm):
        shuffled = corr_matrix.copy()
        shuffled.index = rng.permutation(shuffled.index)
        rbh = find_reciprocal_best_hits(shuffled, threshold=None)
        if len(rbh) > 0:
            null_corrs.extend(rbh["correlation"].tolist())

    if len(null_corrs) == 0:
        return 0.5  # fallback

    return float(np.percentile(null_corrs, 95))
```

- [ ] **Step 7: Run all tests to verify they pass**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_cluster_correspondence.py -v`
Expected: all PASS

- [ ] **Step 8: Write tests for MetaNeighbor AUROC**

Add to `cge_subtype/tests/test_cluster_correspondence.py`:

```python
from cge_subtype.src.cluster_correspondence import compute_metaneighbor_auroc


def test_metaneighbor_perfect_separation():
    """Well-separated clusters should have high AUROC."""
    np.random.seed(42)
    # 2 clusters, clearly separated
    cluster_a = np.random.randn(20, 50) + 5
    cluster_b = np.random.randn(20, 50) - 5
    expr = pd.DataFrame(np.vstack([cluster_a, cluster_b]))
    labels_mouse = pd.Series(["A"] * 20 + ["B"] * 20)
    labels_human = pd.Series(["X"] * 20 + ["Y"] * 20)
    species = pd.Series(["mouse"] * 20 + ["human"] * 20)

    auroc = compute_metaneighbor_auroc(expr, labels_mouse, labels_human, species)
    # Cross-species AUROC for matched types should be high
    assert auroc.loc["A", "X"] > 0.8
    assert auroc.loc["B", "Y"] > 0.8


def test_metaneighbor_random_data():
    """Random data should give AUROC ~0.5."""
    np.random.seed(42)
    expr = pd.DataFrame(np.random.randn(40, 50))
    labels_mouse = pd.Series(["A"] * 10 + ["B"] * 10 + ["C"] * 10 + ["D"] * 10)
    labels_human = labels_mouse.copy()
    species = pd.Series(["mouse"] * 20 + ["human"] * 20)

    auroc = compute_metaneighbor_auroc(expr, labels_mouse, labels_human, species)
    # All AUROCs should be near 0.5
    assert auroc.mean().mean() < 0.7
```

- [ ] **Step 9: Implement MetaNeighbor AUROC**

Add to `cge_subtype/src/cluster_correspondence.py`:

```python
def compute_metaneighbor_auroc(expr_df, labels_mouse, labels_human, species):
    """Compute cross-species cluster-to-cluster AUROC (MetaNeighbor unsupervised).

    For each cell, rank all cells from the OTHER species by Spearman correlation,
    then compute AUROC: do same-label cells rank higher than different-label cells?

    Parameters
    ----------
    expr_df : pd.DataFrame
        Combined cells x genes matrix (mouse + human, subsampled).
    labels_mouse : pd.Series
        Cluster labels for mouse cells (NaN for human cells).
    labels_human : pd.Series
        Cluster labels for human cells (NaN for mouse cells).
    species : pd.Series
        "mouse" or "human" for each cell.

    Returns
    -------
    pd.DataFrame
        Mouse clusters (rows) x human clusters (columns) AUROC matrix.
    """
    mouse_mask = species == "mouse"
    human_mask = species == "human"

    mouse_expr = expr_df.loc[mouse_mask].values
    human_expr = expr_df.loc[human_mask].values
    mouse_labels = labels_mouse.loc[mouse_mask].values
    human_labels = labels_human.loc[human_mask].values

    # Spearman correlation: each mouse cell vs each human cell
    from scipy.stats import rankdata
    mouse_ranked = np.apply_along_axis(rankdata, 1, mouse_expr)
    human_ranked = np.apply_along_axis(rankdata, 1, human_expr)

    # Pearson on ranks = Spearman
    mouse_centered = mouse_ranked - mouse_ranked.mean(axis=1, keepdims=True)
    human_centered = human_ranked - human_ranked.mean(axis=1, keepdims=True)
    mouse_norm = np.sqrt((mouse_centered**2).sum(axis=1, keepdims=True))
    human_norm = np.sqrt((human_centered**2).sum(axis=1, keepdims=True))

    # Correlation matrix: (n_mouse, n_human)
    corr = (mouse_centered @ human_centered.T) / (mouse_norm @ human_norm.T + 1e-10)

    # AUROC: for each mouse cluster vs each human cluster
    mouse_clusters = np.unique(mouse_labels[~pd.isna(mouse_labels)])
    human_clusters = np.unique(human_labels[~pd.isna(human_labels)])

    auroc_matrix = pd.DataFrame(
        np.nan,
        index=mouse_clusters,
        columns=human_clusters,
    )

    for m_clust in mouse_clusters:
        m_idx = np.where(mouse_labels == m_clust)[0]
        for h_clust in human_clusters:
            # For each mouse cell in m_clust, do human cells in h_clust
            # rank higher than human cells NOT in h_clust?
            h_pos = np.where(human_labels == h_clust)[0]
            h_neg = np.where(human_labels != h_clust)[0]

            if len(h_pos) == 0 or len(h_neg) == 0:
                continue

            scores = []
            for mi in m_idx:
                pos_corrs = corr[mi, h_pos]
                neg_corrs = corr[mi, h_neg]
                y_true = np.concatenate([np.ones(len(pos_corrs)), np.zeros(len(neg_corrs))])
                y_score = np.concatenate([pos_corrs, neg_corrs])
                try:
                    scores.append(roc_auc_score(y_true, y_score))
                except ValueError:
                    continue

            if scores:
                auroc_matrix.loc[m_clust, h_clust] = np.mean(scores)

    return auroc_matrix
```

- [ ] **Step 10: Run all tests**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_cluster_correspondence.py -v`
Expected: all PASS

- [ ] **Step 11: Commit Module A library**

```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
git add cge_subtype/src/cluster_correspondence.py cge_subtype/tests/__init__.py cge_subtype/tests/test_cluster_correspondence.py
git commit -m "Add cluster correspondence library (Module A): pseudobulk, Spearman RBH, MetaNeighbor AUROC"
```

---

### Task 2: Scripts — Pseudobulk Computation + Cluster Bridge

**Files:**
- Create: `cge_subtype/scripts/06_compute_pseudobulk.py`
- Create: `cge_subtype/scripts/07_cluster_bridge.py`

- [ ] **Step 1: Write 06_compute_pseudobulk.py**

This script loads each atlas one at a time in backed mode, computes pseudobulk per cluster, and saves CSV files. RAM-safe design is critical.

```python
#!/usr/bin/env python
"""Compute pseudobulk centroids for mouse and human atlases (RAM-safe).

Loads each atlas in backed mode, iterates over clusters, computes mean
expression per cluster, and saves cluster x gene CSV files.

Usage:
    conda activate gencic
    python cge_subtype/scripts/06_compute_pseudobulk.py
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_config():
    config_path = PROJECT_ROOT / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def compute_pseudobulk_backed(h5ad_paths, metadata_path, cluster_col, gene_col="gene_symbol"):
    """Compute pseudobulk from backed h5ad files, one cluster at a time.

    Parameters
    ----------
    h5ad_paths : list[str]
        Paths to h5ad files (will be concatenated if multiple).
    metadata_path : str
        Path to metadata CSV with cluster annotations.
    cluster_col : str
        Column name for cluster identity in metadata.
    gene_col : str
        Column name for gene symbols in adata.var.

    Returns
    -------
    pd.DataFrame
        Clusters x genes mean expression.
    """
    # Load metadata first (small)
    meta = pd.read_csv(metadata_path)
    log.info(f"Loaded metadata: {len(meta)} cells, {meta[cluster_col].nunique()} clusters")

    # Load first h5ad to get gene names
    adata_ref = ad.read_h5ad(h5ad_paths[0], backed="r")
    if gene_col in adata_ref.var.columns:
        gene_names = adata_ref.var[gene_col].values
    else:
        gene_names = adata_ref.var_names.values
    n_genes = len(gene_names)
    log.info(f"Gene space: {n_genes} genes")
    adata_ref.file.close()

    # Build cell_id -> cluster mapping
    cell_to_cluster = dict(zip(meta["cell_label"], meta[cluster_col]))
    clusters = sorted(meta[cluster_col].unique())
    log.info(f"Computing pseudobulk for {len(clusters)} clusters")

    # Accumulate sums and counts per cluster
    cluster_sums = {c: np.zeros(n_genes, dtype=np.float64) for c in clusters}
    cluster_counts = {c: 0 for c in clusters}

    for h5ad_path in h5ad_paths:
        log.info(f"Processing {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path, backed="r")
        cell_ids = adata.obs_names.values

        # Process in chunks of 1000 cells
        chunk_size = 1000
        for start in range(0, len(cell_ids), chunk_size):
            end = min(start + chunk_size, len(cell_ids))
            chunk_ids = cell_ids[start:end]

            # Load chunk into memory
            chunk_expr = adata[start:end].X
            if hasattr(chunk_expr, "toarray"):
                chunk_expr = chunk_expr.toarray()
            chunk_expr = np.asarray(chunk_expr, dtype=np.float64)

            for i, cid in enumerate(chunk_ids):
                clust = cell_to_cluster.get(cid)
                if clust is not None:
                    cluster_sums[clust] += chunk_expr[i]
                    cluster_counts[clust] += 1

            if (start // chunk_size) % 50 == 0:
                log.info(f"  Processed {end}/{len(cell_ids)} cells")

        adata.file.close()

    # Compute means
    result = {}
    for c in clusters:
        if cluster_counts[c] > 0:
            result[c] = cluster_sums[c] / cluster_counts[c]
        else:
            log.warning(f"Cluster {c} has 0 cells, skipping")

    return pd.DataFrame(result, index=gene_names).T


def main():
    parser = argparse.ArgumentParser(description="Compute pseudobulk centroids")
    parser.add_argument("--species", choices=["mouse", "human", "both"], default="both")
    parser.add_argument("--outdir", default=str(PROJECT_ROOT / "results" / "pseudobulk"))
    args = parser.parse_args()

    config = load_config()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.species in ("mouse", "both"):
        log.info("=== Computing mouse pseudobulk ===")
        mouse_paths = [
            "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/"
            "WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-raw.h5ad",
            "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/"
            "WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-2-raw.h5ad",
        ]
        mouse_meta = (
            "/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/"
            "20230830/views/cell_metadata_with_cluster_annotation.csv"
        )
        mouse_pb = compute_pseudobulk_backed(mouse_paths, mouse_meta, cluster_col="cluster_alias")
        out_path = outdir / "mouse_pseudobulk.csv"
        mouse_pb.to_csv(out_path)
        log.info(f"Saved mouse pseudobulk: {mouse_pb.shape} to {out_path}")

    if args.species in ("human", "both"):
        log.info("=== Computing human pseudobulk ===")
        # Human: iterate over supercluster files
        import glob
        human_paths = sorted(glob.glob(
            "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_*.h5ad"
        ))
        # Human metadata — need to identify the metadata file/column
        # This will be determined during data exploration (Task 2, Step 2)
        # Placeholder: use adata.obs cluster column directly
        log.info(f"Found {len(human_paths)} human supercluster files")

        # For human, cluster labels are in adata.obs — process differently
        human_sums = {}
        human_counts = {}

        for hp in human_paths:
            log.info(f"Processing {hp}")
            adata = ad.read_h5ad(hp, backed="r")

            # Determine cluster column (check obs columns)
            cluster_col = None
            for candidate in ["cluster_id", "cluster", "supertype", "cell_type"]:
                if candidate in adata.obs.columns:
                    cluster_col = candidate
                    break

            if cluster_col is None:
                log.warning(f"No cluster column found in {hp}, skipping")
                adata.file.close()
                continue

            # Gene symbol column may be "gene_symbol", "Gene", or use var_names
            gene_col_found = None
            for gc in ["gene_symbol", "Gene", "gene_name"]:
                if gc in adata.var.columns:
                    gene_col_found = gc
                    break
            gene_names = adata.var[gene_col_found].values if gene_col_found else adata.var_names.values
            cell_ids = adata.obs_names.values
            labels = adata.obs[cluster_col].values

            chunk_size = 1000
            for start in range(0, len(cell_ids), chunk_size):
                end = min(start + chunk_size, len(cell_ids))
                chunk_expr = adata[start:end].X
                if hasattr(chunk_expr, "toarray"):
                    chunk_expr = chunk_expr.toarray()
                chunk_expr = np.asarray(chunk_expr, dtype=np.float64)

                for i in range(chunk_expr.shape[0]):
                    clust = str(labels[start + i])
                    if clust not in human_sums:
                        human_sums[clust] = np.zeros(len(gene_names), dtype=np.float64)
                        human_counts[clust] = 0
                    human_sums[clust] += chunk_expr[i]
                    human_counts[clust] += 1

            adata.file.close()

        human_result = {}
        for c in human_sums:
            if human_counts[c] > 0:
                human_result[c] = human_sums[c] / human_counts[c]

        human_pb = pd.DataFrame(human_result, index=gene_names).T
        out_path = outdir / "human_pseudobulk.csv"
        human_pb.to_csv(out_path)
        log.info(f"Saved human pseudobulk: {human_pb.shape} to {out_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Data exploration — verify atlas structure**

Before running, verify the h5ad structure and metadata columns. This is a manual exploration step:

Run:
```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic python -c "
import anndata as ad
# Mouse
adata = ad.read_h5ad('/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-raw.h5ad', backed='r')
print('Mouse obs columns:', list(adata.obs.columns[:10]))
print('Mouse var columns:', list(adata.var.columns[:10]))
print('Mouse shape:', adata.shape)
adata.file.close()
# Human
adata = ad.read_h5ad('/mnt/data0/HumanBrainCellType/SuperTypeRawDat/$(ls /mnt/data0/HumanBrainCellType/SuperTypeRawDat/ | head -1)', backed='r')
print('Human obs columns:', list(adata.obs.columns[:10]))
print('Human var columns:', list(adata.var.columns[:10]))
print('Human shape:', adata.shape)
adata.file.close()
"
```

Expected: column names that identify cluster identity and gene symbols. Update script if column names differ from assumptions.

- [ ] **Step 3: Write 07_cluster_bridge.py**

```python
#!/usr/bin/env python
"""Cross-species cluster bridge: Spearman RBH (Pass 1) + MetaNeighbor AUROC (Pass 2).

Reads pseudobulk centroids from 06_compute_pseudobulk.py output,
applies ortholog conversion, and computes cluster correspondences.

Usage:
    conda activate gencic
    python cge_subtype/scripts/07_cluster_bridge.py
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from cge_subtype.src.cluster_correspondence import (
    compute_metaneighbor_auroc,
    compute_spearman_corr_matrix,
    determine_rbh_threshold,
    find_reciprocal_best_hits,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_ortholog_mapping(path):
    """Load mouse->human ortholog mapping. Returns dict mouse_gene -> human_gene.

    Expected CSV columns: mouse_symbol, human_symbol, human_ensembl,
    mapping_source, is_one_to_one, many_mouse_to_one_human.
    """
    df = pd.read_csv(path)
    # Keep only 1:1 orthologs
    if "is_one_to_one" in df.columns:
        one_to_one = df[df["is_one_to_one"] == True]
    elif "ortholog_type" in df.columns:
        one_to_one = df[df["ortholog_type"] == "one2one"]
    else:
        one_to_one = df
    return dict(zip(one_to_one["mouse_symbol"], one_to_one["human_symbol"]))


def select_hvgs(combined_centroids, n_hvgs=3000):
    """Select top HVGs by variance across cluster centroids."""
    variances = combined_centroids.var(axis=0)
    top_genes = variances.nlargest(n_hvgs).index
    return top_genes.tolist()


def main():
    parser = argparse.ArgumentParser(description="Cross-species cluster bridge")
    parser.add_argument("--mouse-pb", default=str(PROJECT_ROOT / "results/pseudobulk/mouse_pseudobulk.csv"))
    parser.add_argument("--human-pb", default=str(PROJECT_ROOT / "results/pseudobulk/human_pseudobulk.csv"))
    parser.add_argument(
        "--ortholog-map",
        default="/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/cross_species/orthologs/ortholog_mapping.csv",
    )
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--outdir", default=str(PROJECT_ROOT / "results" / "cluster_bridge"))
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load pseudobulk
    log.info("Loading pseudobulk centroids")
    mouse_pb = pd.read_csv(args.mouse_pb, index_col=0)
    human_pb = pd.read_csv(args.human_pb, index_col=0)
    log.info(f"Mouse: {mouse_pb.shape}, Human: {human_pb.shape}")

    # 2. Ortholog conversion: rename mouse genes to human
    log.info("Loading ortholog mapping")
    ortho_map = load_ortholog_mapping(args.ortholog_map)
    mouse_genes_converted = {g: ortho_map.get(g, None) for g in mouse_pb.columns}
    valid = {k: v for k, v in mouse_genes_converted.items() if v is not None}
    log.info(f"Mapped {len(valid)}/{len(mouse_pb.columns)} mouse genes to human orthologs")

    mouse_pb_converted = mouse_pb[list(valid.keys())].copy()
    mouse_pb_converted.columns = [valid[g] for g in mouse_pb_converted.columns]

    # Remove duplicate gene columns (if multiple mouse genes map to same human gene)
    mouse_pb_converted = mouse_pb_converted.loc[:, ~mouse_pb_converted.columns.duplicated()]

    # 3. Shared genes
    shared_genes = mouse_pb_converted.columns.intersection(human_pb.columns)
    log.info(f"Shared genes after ortholog conversion: {len(shared_genes)}")

    mouse_shared = mouse_pb_converted[shared_genes]
    human_shared = human_pb[shared_genes]

    # 4. HVG selection
    combined = pd.concat([mouse_shared, human_shared], axis=0)
    hvgs = select_hvgs(combined, n_hvgs=args.n_hvgs)
    log.info(f"Selected {len(hvgs)} HVGs")

    mouse_hvg = mouse_shared[hvgs]
    human_hvg = human_shared[hvgs]

    # 5. Pass 1: Spearman correlation + RBH
    log.info("Computing Spearman correlation matrix")
    corr = compute_spearman_corr_matrix(mouse_hvg, human_hvg)
    corr.to_csv(outdir / "spearman_corr_matrix.csv")
    log.info(f"Correlation matrix: {corr.shape}")

    # Determine threshold from data
    log.info(f"Determining RBH threshold via permutation (n_perm={args.n_perm})")
    threshold = determine_rbh_threshold(corr, method="permutation", n_perm=args.n_perm)
    log.info(f"Data-driven RBH threshold: {threshold:.4f}")

    rbh = find_reciprocal_best_hits(corr, threshold=threshold)
    rbh.to_csv(outdir / "rbh_resolved_pairs.csv", index=False)
    log.info(f"Pass 1 resolved {len(rbh)} RBH pairs (threshold={threshold:.4f})")

    # Identify unresolved clusters
    resolved_mouse = set(rbh["mouse_cluster"])
    resolved_human = set(rbh["human_cluster"])
    unresolved_mouse = [c for c in corr.index if c not in resolved_mouse]
    unresolved_human = [c for c in corr.columns if c not in resolved_human]
    log.info(f"Unresolved: {len(unresolved_mouse)} mouse, {len(unresolved_human)} human clusters")

    # Save unresolved list for Pass 2
    pd.DataFrame({"mouse_unresolved": pd.Series(unresolved_mouse)}).to_csv(
        outdir / "unresolved_mouse.csv", index=False
    )
    pd.DataFrame({"human_unresolved": pd.Series(unresolved_human)}).to_csv(
        outdir / "unresolved_human.csv", index=False
    )

    # Save summary
    summary = {
        "total_mouse_clusters": len(corr.index),
        "total_human_clusters": len(corr.columns),
        "rbh_threshold": threshold,
        "rbh_resolved_pairs": len(rbh),
        "unresolved_mouse": len(unresolved_mouse),
        "unresolved_human": len(unresolved_human),
    }
    pd.Series(summary).to_csv(outdir / "pass1_summary.csv")
    log.info(f"Pass 1 summary: {summary}")

    log.info("Pass 1 complete. Run MetaNeighbor (Pass 2) on unresolved clusters with subsampled cells.")
    log.info("Pass 2 requires loading subsampled cells — see spec Module A Pass 2.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run pseudobulk computation (mouse first, then human)**

Run:
```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic python cge_subtype/scripts/06_compute_pseudobulk.py --species mouse
```

Expected: `results/pseudobulk/mouse_pseudobulk.csv` created. Monitor RAM usage — should stay below 30 GB.

Then:
```bash
conda run -n gencic python cge_subtype/scripts/06_compute_pseudobulk.py --species human
```

Expected: `results/pseudobulk/human_pseudobulk.csv` created.

- [ ] **Step 5: Run cluster bridge**

Run:
```bash
conda run -n gencic python cge_subtype/scripts/07_cluster_bridge.py --n-perm 1000
```

Expected: Output files in `results/cluster_bridge/`:
- `spearman_corr_matrix.csv`
- `rbh_resolved_pairs.csv`
- `pass1_summary.csv`
- `unresolved_mouse.csv`, `unresolved_human.csv`

Inspect results: how many pairs resolved? What's the threshold? Are CGE clusters in the unresolved set (expected)?

- [ ] **Step 6: Write 07b_metaneighbor_pass2.py (Pass 2 for unresolved clusters)**

```python
#!/usr/bin/env python
"""MetaNeighbor Pass 2: Resolve ambiguous clusters from Pass 1.

Subsamples cells from unresolved clusters and computes cross-species
AUROC scores using the MetaNeighbor algorithm (Python-native).

Usage:
    conda activate gencic
    python cge_subtype/scripts/07b_metaneighbor_pass2.py
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cge_subtype.src.cluster_correspondence import compute_metaneighbor_auroc
from cge_subtype.src.harmony_mapping import subsample_reference

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_ortholog_mapping(path):
    """Load mouse->human ortholog mapping."""
    df = pd.read_csv(path)
    if "is_one_to_one" in df.columns:
        one_to_one = df[df["is_one_to_one"] == True]
    else:
        one_to_one = df
    return dict(zip(one_to_one["mouse_symbol"], one_to_one["human_symbol"]))


def main():
    parser = argparse.ArgumentParser(description="MetaNeighbor Pass 2 on unresolved clusters")
    parser.add_argument("--max-per-cluster", type=int, default=200)
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument(
        "--ortholog-map",
        default="/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/cross_species/orthologs/ortholog_mapping.csv",
    )
    parser.add_argument("--bridge-dir", default=str(PROJECT_ROOT / "results" / "cluster_bridge"))
    parser.add_argument("--outdir", default=str(PROJECT_ROOT / "results" / "cluster_bridge"))
    args = parser.parse_args()

    bridge_dir = Path(args.bridge_dir)
    outdir = Path(args.outdir)

    # Load unresolved cluster lists from Pass 1
    unresolved_mouse = pd.read_csv(bridge_dir / "unresolved_mouse.csv")["mouse_unresolved"].tolist()
    unresolved_human = pd.read_csv(bridge_dir / "unresolved_human.csv")["human_unresolved"].tolist()
    log.info(f"Unresolved: {len(unresolved_mouse)} mouse, {len(unresolved_human)} human clusters")

    if len(unresolved_mouse) == 0 or len(unresolved_human) == 0:
        log.info("No unresolved clusters — Pass 2 not needed")
        return

    # Load ortholog mapping
    ortho_map = load_ortholog_mapping(args.ortholog_map)

    # Subsample mouse atlas cells from unresolved clusters only
    log.info("Loading mouse atlas cells for unresolved clusters")
    mouse_meta = pd.read_csv(
        "/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/"
        "20230830/views/cell_metadata_with_cluster_annotation.csv"
    )
    mouse_unresolved_cells = mouse_meta[mouse_meta["cluster_alias"].isin(unresolved_mouse)]
    log.info(f"Mouse unresolved cells: {len(mouse_unresolved_cells)}")

    # Subsample to max_per_cluster
    mouse_sampled = mouse_unresolved_cells.groupby("cluster_alias").apply(
        lambda x: x.sample(min(len(x), args.max_per_cluster), random_state=42)
    ).reset_index(drop=True)
    log.info(f"Mouse sampled: {len(mouse_sampled)} cells")

    # Load expression for sampled mouse cells (chunked from backed h5ad)
    mouse_cell_ids = set(mouse_sampled["cell_label"])
    mouse_expr_rows = []
    for part in [1, 2]:
        h5ad_path = (
            f"/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/"
            f"WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-{part}-raw.h5ad"
        )
        adata = ad.read_h5ad(h5ad_path, backed="r")
        mask = np.isin(adata.obs_names.values, list(mouse_cell_ids))
        if mask.sum() > 0:
            chunk = adata[mask].to_memory()
            mouse_expr_rows.append(chunk)
            log.info(f"  Part {part}: loaded {mask.sum()} cells")
        adata.file.close()

    mouse_adata = ad.concat(mouse_expr_rows) if len(mouse_expr_rows) > 1 else mouse_expr_rows[0]

    # Similarly load human cells from unresolved clusters
    log.info("Loading human atlas cells for unresolved clusters")
    import glob
    human_paths = sorted(glob.glob(
        "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_*.h5ad"
    ))
    human_expr_rows = []
    for hp in human_paths:
        adata = ad.read_h5ad(hp, backed="r")
        cluster_col = None
        for candidate in ["cluster_id", "cluster", "supertype"]:
            if candidate in adata.obs.columns:
                cluster_col = candidate
                break
        if cluster_col is None:
            adata.file.close()
            continue
        mask = adata.obs[cluster_col].astype(str).isin(unresolved_human)
        if mask.sum() > 0:
            chunk = adata[mask.values].to_memory()
            # Subsample per cluster
            sc_labels = chunk.obs[cluster_col].astype(str)
            keep = []
            for clust in sc_labels.unique():
                idx = np.where(sc_labels == clust)[0]
                if len(idx) > args.max_per_cluster:
                    idx = np.random.default_rng(42).choice(idx, args.max_per_cluster, replace=False)
                keep.extend(idx)
            chunk = chunk[sorted(keep)]
            human_expr_rows.append(chunk)
            log.info(f"  {Path(hp).name}: loaded {len(chunk)} cells")
        adata.file.close()

    human_adata = ad.concat(human_expr_rows)
    log.info(f"Human sampled: {human_adata.n_obs} cells")

    # Convert mouse genes to human orthologs
    mouse_genes = mouse_adata.var_names.values
    human_gene_map = {g: ortho_map.get(g) for g in mouse_genes if ortho_map.get(g) is not None}
    mouse_adata = mouse_adata[:, list(human_gene_map.keys())].copy()
    mouse_adata.var_names = pd.Index([human_gene_map[g] for g in mouse_adata.var_names])

    # Shared genes
    shared = mouse_adata.var_names.intersection(human_adata.var_names)
    log.info(f"Shared ortholog genes: {len(shared)}")

    mouse_expr = pd.DataFrame(
        mouse_adata[:, shared].X.toarray() if hasattr(mouse_adata.X, "toarray") else mouse_adata[:, shared].X,
        columns=shared,
    )
    human_expr = pd.DataFrame(
        human_adata[:, shared].X.toarray() if hasattr(human_adata.X, "toarray") else human_adata[:, shared].X,
        columns=shared,
    )

    # HVG selection on combined
    combined_expr = pd.concat([mouse_expr, human_expr], ignore_index=True)
    variances = combined_expr.var(axis=0)
    hvgs = variances.nlargest(args.n_hvgs).index.tolist()
    mouse_expr = mouse_expr[hvgs]
    human_expr = human_expr[hvgs]
    combined_expr = pd.concat([mouse_expr, human_expr], ignore_index=True)

    # Build labels
    mouse_labels = mouse_sampled.set_index("cell_label").loc[
        mouse_adata.obs_names, "cluster_alias"
    ].values
    human_labels = human_adata.obs[cluster_col].astype(str).values

    labels_mouse = pd.Series(
        list(mouse_labels) + [np.nan] * len(human_labels)
    )
    labels_human = pd.Series(
        [np.nan] * len(mouse_labels) + list(human_labels)
    )
    species = pd.Series(
        ["mouse"] * len(mouse_labels) + ["human"] * len(human_labels)
    )

    # Compute MetaNeighbor AUROC
    log.info("Computing MetaNeighbor AUROC")
    auroc = compute_metaneighbor_auroc(combined_expr, labels_mouse, labels_human, species)
    auroc.to_csv(outdir / "metaneighbor_auroc_matrix.csv")
    log.info(f"AUROC matrix: {auroc.shape}")

    # Find best matches from AUROC
    pass2_pairs = []
    for m_clust in auroc.index:
        best_h = auroc.loc[m_clust].idxmax()
        best_auroc = auroc.loc[m_clust, best_h]
        if not np.isnan(best_auroc):
            pass2_pairs.append({
                "mouse_cluster": m_clust,
                "human_cluster": best_h,
                "auroc": best_auroc,
                "method": "metaneighbor",
            })

    pass2_df = pd.DataFrame(pass2_pairs)
    pass2_df.to_csv(outdir / "pass2_metaneighbor_pairs.csv", index=False)
    log.info(f"Pass 2 resolved {len(pass2_df)} pairs")

    # Merge Pass 1 + Pass 2
    pass1 = pd.read_csv(bridge_dir / "rbh_resolved_pairs.csv")
    combined_pairs = pd.concat([pass1, pass2_df], ignore_index=True)
    combined_pairs.to_csv(outdir / "all_cluster_correspondences.csv", index=False)
    log.info(f"Total correspondences: {len(combined_pairs)} (Pass 1: {len(pass1)}, Pass 2: {len(pass2_df)})")


if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Write 07c_subsample_reference.py (create subsampled reference h5ad for Module B)**

This creates the subsampled reference that Module B's Harmony mapping needs.

```python
#!/usr/bin/env python
"""Create subsampled reference h5ad for Harmony mapping (Module B).

Loads Allen WMB-10Xv3, subsamples to max_per_cluster cells, saves h5ad.

Usage:
    conda activate gencic
    python cge_subtype/scripts/07c_subsample_reference.py --max-per-cluster 200
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-per-cluster", type=int, default=200)
    parser.add_argument("--outdir", default=str(PROJECT_ROOT / "results" / "preprocessed"))
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load metadata
    meta = pd.read_csv(
        "/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/"
        "20230830/views/cell_metadata_with_cluster_annotation.csv"
    )
    log.info(f"Total cells in metadata: {len(meta)}")

    # Subsample per cluster
    rng = np.random.default_rng(42)
    sampled_cells = set()
    for clust in meta["cluster_alias"].unique():
        clust_ids = meta[meta["cluster_alias"] == clust]["cell_label"].values
        if len(clust_ids) > args.max_per_cluster:
            clust_ids = rng.choice(clust_ids, args.max_per_cluster, replace=False)
        sampled_cells.update(clust_ids)
    log.info(f"Sampled {len(sampled_cells)} cells from {meta['cluster_alias'].nunique()} clusters")

    # Load from h5ad files, keeping only sampled cells
    parts = []
    for part_num in [1, 2]:
        h5ad_path = (
            f"/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/"
            f"WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-{part_num}-raw.h5ad"
        )
        log.info(f"Loading {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path, backed="r")
        mask = np.isin(adata.obs_names.values, list(sampled_cells))
        if mask.sum() > 0:
            chunk = adata[mask].to_memory()
            parts.append(chunk)
            log.info(f"  Part {part_num}: kept {mask.sum()} cells")
        adata.file.close()

    # Concatenate
    ref = ad.concat(parts)

    # Add metadata
    meta_indexed = meta.set_index("cell_label")
    for col in ["cluster_alias", "subclass", "class"]:
        if col in meta_indexed.columns:
            ref.obs[col] = meta_indexed.loc[ref.obs_names, col].values

    # Save
    out_path = outdir / "reference_subsampled.h5ad"
    ref.write_h5ad(out_path)
    log.info(f"Saved subsampled reference: {ref.shape} to {out_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 8: Run subsampled reference creation (must complete before Module B sweep)**

Run:
```bash
conda run -n gencic python cge_subtype/scripts/07c_subsample_reference.py --max-per-cluster 200
```

Expected: `results/preprocessed/reference_subsampled.h5ad` created (~50K cells). Monitor RAM.

- [ ] **Step 9: Run MetaNeighbor Pass 2**

Run:
```bash
conda run -n gencic python cge_subtype/scripts/07b_metaneighbor_pass2.py
```

Expected: Output in `results/cluster_bridge/`:
- `metaneighbor_auroc_matrix.csv` — AUROC scores for unresolved clusters
- `pass2_metaneighbor_pairs.csv` — best matches from Pass 2
- `all_cluster_correspondences.csv` — merged Pass 1 + Pass 2

- [ ] **Step 10: Commit Module A scripts**

```bash
git add cge_subtype/scripts/06_compute_pseudobulk.py cge_subtype/scripts/07_cluster_bridge.py \
    cge_subtype/scripts/07b_metaneighbor_pass2.py cge_subtype/scripts/07c_subsample_reference.py
git commit -m "Add Module A scripts: pseudobulk, cluster bridge (Pass 1 + Pass 2), subsampled reference"
```

---

## Chunk 2: Module B — Improved Mouse Patch-seq Mapping

### Task 3: Core Library — harmony_mapping.py

**Files:**
- Create: `cge_subtype/src/harmony_mapping.py`

- [ ] **Step 1: Implement harmony_mapping.py**

```python
# cge_subtype/src/harmony_mapping.py
"""Harmony-based mapping of patch-seq cells to reference atlas."""

import logging
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def subsample_reference(adata, cluster_col, max_per_cluster=200, seed=42):
    """Subsample reference atlas to max_per_cluster cells per cluster.

    Parameters
    ----------
    adata : ad.AnnData
        Full reference atlas.
    cluster_col : str
        Column in adata.obs for cluster identity.
    max_per_cluster : int
        Maximum cells per cluster.
    seed : int
        Random seed.

    Returns
    -------
    ad.AnnData
        Subsampled AnnData.
    """
    rng = np.random.default_rng(seed)
    keep_idx = []
    for clust in adata.obs[cluster_col].unique():
        mask = adata.obs[cluster_col] == clust
        idx = np.where(mask)[0]
        if len(idx) > max_per_cluster:
            idx = rng.choice(idx, max_per_cluster, replace=False)
        keep_idx.extend(idx)
    keep_idx = sorted(keep_idx)
    log.info(f"Subsampled {len(keep_idx)} cells from {adata.n_obs} ({adata.obs[cluster_col].nunique()} clusters)")
    return adata[keep_idx].copy()


def run_harmony_mapping(ref_adata, query_adata, cluster_col, n_pcs=50, n_hvgs=3000,
                        theta=2.0, n_neighbors=30, batch_key="technology"):
    """Map query cells to reference using Harmony + kNN label transfer.

    Parameters
    ----------
    ref_adata : ad.AnnData
        Subsampled reference atlas (raw counts in .X or .layers["counts"]).
    query_adata : ad.AnnData
        Query patch-seq cells (raw counts).
    cluster_col : str
        Cluster label column in ref_adata.obs.
    n_pcs : int
        Number of PCA components.
    n_hvgs : int
        Number of highly variable genes.
    theta : float
        Harmony diversity penalty.
    n_neighbors : int
        k for kNN label transfer.
    batch_key : str
        Column distinguishing reference vs query in combined AnnData.

    Returns
    -------
    pd.DataFrame
        Per-query-cell: predicted_cluster, confidence, n_neighbors_agreeing.
    """
    import harmonypy

    # Tag datasets
    ref_adata.obs[batch_key] = "reference"
    query_adata.obs[batch_key] = "query"

    # Combine
    combined = ad.concat([ref_adata, query_adata], join="inner")
    log.info(f"Combined: {combined.n_obs} cells, {combined.n_vars} genes")

    # Normalize
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    # HVG on reference only
    ref_mask = combined.obs[batch_key] == "reference"
    sc.pp.highly_variable_genes(combined, n_top_genes=n_hvgs, batch_key=batch_key, subset=False)
    combined = combined[:, combined.var["highly_variable"]].copy()
    log.info(f"After HVG selection: {combined.n_vars} genes")

    # PCA
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, n_comps=n_pcs)

    # Harmony
    ho = harmonypy.run_harmony(combined.obsm["X_pca"], combined.obs, batch_key, theta=theta)
    combined.obsm["X_harmony"] = ho.Z_corr.T

    # kNN label transfer
    from sklearn.neighbors import NearestNeighbors

    ref_idx = np.where(combined.obs[batch_key] == "reference")[0]
    query_idx = np.where(combined.obs[batch_key] == "query")[0]

    ref_emb = combined.obsm["X_harmony"][ref_idx]
    query_emb = combined.obsm["X_harmony"][query_idx]
    ref_labels = combined.obs[cluster_col].values[ref_idx]

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(ref_emb)
    distances, indices = nn.kneighbors(query_emb)

    results = []
    query_cell_ids = combined.obs_names[query_idx]
    for i, cell_id in enumerate(query_cell_ids):
        neighbor_labels = ref_labels[indices[i]]
        label_counts = pd.Series(neighbor_labels).value_counts()
        predicted = label_counts.index[0]
        confidence = label_counts.iloc[0] / n_neighbors
        results.append({
            "cell_id": cell_id,
            "predicted_cluster": predicted,
            "confidence": confidence,
            "n_neighbors_agreeing": label_counts.iloc[0],
        })

    return pd.DataFrame(results)


def evaluate_mapping(predictions, ground_truth_col, metadata):
    """Evaluate mapping accuracy against Cre-line ground truth.

    Parameters
    ----------
    predictions : pd.DataFrame
        Output of run_harmony_mapping.
    ground_truth_col : str
        Column in metadata containing ground truth labels.
    metadata : pd.DataFrame
        Metadata with cell_id index and ground truth column.

    Returns
    -------
    dict
        Accuracy metrics: overall_accuracy, subclass_accuracy, confusion_matrix,
        non_neuronal_rate.
    """
    merged = predictions.merge(
        metadata[[ground_truth_col]],
        left_on="cell_id",
        right_index=True,
        how="inner",
    )

    # Extract subclass from ground truth (first word of AIT alias)
    merged["gt_subclass"] = merged[ground_truth_col].str.split().str[0]

    # Extract subclass from predicted cluster (need mapping — dataset-specific)
    # For now, return raw confusion matrix
    from sklearn.metrics import accuracy_score, confusion_matrix

    overall_acc = accuracy_score(merged["gt_subclass"], merged["predicted_cluster"])

    # Non-neuronal rate
    non_neuronal_types = ["Astro", "Oligo", "OPC", "Micro", "Endo", "VLMC", "Peri"]
    nn_mask = merged["predicted_cluster"].isin(non_neuronal_types)
    nn_rate = nn_mask.mean()

    cm = pd.crosstab(merged["gt_subclass"], merged["predicted_cluster"])

    return {
        "overall_accuracy": overall_acc,
        "non_neuronal_rate": nn_rate,
        "n_cells": len(merged),
        "confusion_matrix": cm,
    }
```

- [ ] **Step 2: Commit**

```bash
git add cge_subtype/src/harmony_mapping.py
git commit -m "Add Harmony mapping library (Module B): subsample, integrate, kNN transfer, evaluate"
```

---

### Task 4: Scripts — Harmony Mapping + Hyperparameter Sweep

**Files:**
- Create: `cge_subtype/scripts/08_harmony_mapping.py`
- Create: `cge_subtype/scripts/09_harmony_sweep.py`

- [ ] **Step 1: Write 08_harmony_mapping.py (single-run mapper)**

```python
#!/usr/bin/env python
"""Map mouse patch-seq cells to Allen WMB-10Xv3 reference using Harmony.

Usage:
    conda activate gencic
    python cge_subtype/scripts/08_harmony_mapping.py \
        --dataset M1 --n-pcs 50 --n-hvgs 3000 --theta 2.0 --n-neighbors 30
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from cge_subtype.src.harmony_mapping import (
    evaluate_mapping,
    run_harmony_mapping,
    subsample_reference,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATASET_CONFIGS = {
    "M1": {
        "counts": "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_counts.csv.gz",
        "metadata": "/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv",
        "gt_col": "RNA family",
    },
    "V1": {
        "counts": "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_counts.csv",
        "metadata": "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_metadata.csv",
        "gt_col": "corresponding_AIT2.3.1_alias",
    },
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["M1", "V1"], required=True)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument("--theta", type=float, default=2.0)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--max-per-cluster", type=int, default=200)
    parser.add_argument("--neuronal-only", action="store_true")
    parser.add_argument("--outdir", default=str(PROJECT_ROOT / "results" / "harmony_mapping"))
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cfg = DATASET_CONFIGS[args.dataset]

    # Load reference (subsampled)
    log.info("Loading reference atlas (subsampled)")
    # Load preprocessed reference if available, otherwise subsample from raw
    ref_preprocessed = PROJECT_ROOT / "results" / "preprocessed" / "reference_subsampled.h5ad"
    if ref_preprocessed.exists():
        ref = ad.read_h5ad(ref_preprocessed)
    else:
        log.info("No preprocessed reference found. Loading and subsampling from raw...")
        # This would need the full atlas loading logic from 06_compute_pseudobulk.py
        # For now, raise an error — the sweep script handles this
        raise FileNotFoundError(
            f"Preprocessed reference not found at {ref_preprocessed}. "
            "Run the subsample step first."
        )

    if args.neuronal_only:
        # Allen WMB class names may use various formats — check both
        neuronal_classes = [
            "Glutamatergic", "GABAergic",  # simple format
        ]
        # Also match pattern-based names like "01 IT-ET Glut", "06 CTX-CGE GABA"
        class_col = "class"
        if class_col in ref.obs.columns:
            is_neuronal = (
                ref.obs[class_col].isin(neuronal_classes)
                | ref.obs[class_col].str.contains("Glut|GABA", case=False, na=False)
            )
            ref = ref[is_neuronal].copy()
        log.info(f"Neuronal-only: {ref.n_obs} cells")

    # Load query
    log.info(f"Loading {args.dataset} patch-seq data")
    if cfg["counts"].endswith(".gz"):
        counts = pd.read_csv(cfg["counts"], index_col=0, compression="gzip")
    else:
        counts = pd.read_csv(cfg["counts"], index_col=0)
    query = ad.AnnData(counts)
    log.info(f"Query: {query.n_obs} cells, {query.n_vars} genes")

    # Map
    log.info("Running Harmony mapping")
    predictions = run_harmony_mapping(
        ref, query,
        cluster_col="cluster_alias",
        n_pcs=args.n_pcs,
        n_hvgs=args.n_hvgs,
        theta=args.theta,
        n_neighbors=args.n_neighbors,
    )

    # Save predictions
    tag = f"{args.dataset}_pcs{args.n_pcs}_hvg{args.n_hvgs}_theta{args.theta}_k{args.n_neighbors}"
    if args.neuronal_only:
        tag += "_neuronalonly"
    predictions.to_csv(outdir / f"{tag}_predictions.csv", index=False)

    # Evaluate
    metadata = pd.read_csv(cfg["metadata"], sep="\t" if cfg["metadata"].endswith(".tsv") else ",", index_col=0)
    metrics = evaluate_mapping(predictions, cfg["gt_col"], metadata)
    log.info(f"Accuracy: {metrics['overall_accuracy']:.3f}, NN rate: {metrics['non_neuronal_rate']:.3f}")

    # Save metrics
    pd.Series({
        "dataset": args.dataset,
        "n_pcs": args.n_pcs,
        "n_hvgs": args.n_hvgs,
        "theta": args.theta,
        "n_neighbors": args.n_neighbors,
        "neuronal_only": args.neuronal_only,
        "accuracy": metrics["overall_accuracy"],
        "non_neuronal_rate": metrics["non_neuronal_rate"],
        "n_cells": metrics["n_cells"],
    }).to_csv(outdir / f"{tag}_metrics.csv")

    metrics["confusion_matrix"].to_csv(outdir / f"{tag}_confusion.csv")
    log.info(f"Results saved to {outdir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write 09_harmony_sweep.py (hyperparameter sweep)**

```python
#!/usr/bin/env python
"""Hyperparameter sweep for Harmony mapping.

Runs all combinations of hyperparameters and collects metrics.

Usage:
    conda activate gencic
    python cge_subtype/scripts/09_harmony_sweep.py --dataset M1 --n-processes 10
"""

import argparse
import itertools
import logging
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]

HARMONY_GRID = {
    "n_pcs": [30, 50, 75],
    "n_hvgs": [2000, 3000, 5000],
    "theta": [1.0, 2.0],
    "n_neighbors": [15, 30, 50],
}


def run_single_config(config_tuple):
    """Run a single Harmony mapping configuration."""
    dataset, n_pcs, n_hvgs, theta, n_neighbors, neuronal_only = config_tuple
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "08_harmony_mapping.py"),
        "--dataset", dataset,
        "--n-pcs", str(n_pcs),
        "--n-hvgs", str(n_hvgs),
        "--theta", str(theta),
        "--n-neighbors", str(n_neighbors),
    ]
    if neuronal_only:
        cmd.append("--neuronal-only")

    tag = f"{dataset}_pcs{n_pcs}_hvg{n_hvgs}_theta{theta}_k{n_neighbors}"
    if neuronal_only:
        tag += "_neuronalonly"

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            log.error(f"FAILED {tag}: {result.stderr[-200:]}")
            return {"tag": tag, "status": "failed", "error": result.stderr[-200:]}
        log.info(f"DONE {tag}")
        return {"tag": tag, "status": "success"}
    except subprocess.TimeoutExpired:
        log.error(f"TIMEOUT {tag}")
        return {"tag": tag, "status": "timeout"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["M1", "V1", "both"], default="both")
    parser.add_argument("--n-processes", type=int, default=10)
    parser.add_argument("--include-neuronal-only", action="store_true")
    args = parser.parse_args()

    datasets = ["M1", "V1"] if args.dataset == "both" else [args.dataset]

    # Build config grid
    configs = []
    for dataset in datasets:
        for n_pcs, n_hvgs, theta, n_neighbors in itertools.product(
            HARMONY_GRID["n_pcs"],
            HARMONY_GRID["n_hvgs"],
            HARMONY_GRID["theta"],
            HARMONY_GRID["n_neighbors"],
        ):
            configs.append((dataset, n_pcs, n_hvgs, theta, n_neighbors, False))
            if args.include_neuronal_only:
                configs.append((dataset, n_pcs, n_hvgs, theta, n_neighbors, True))

    log.info(f"Running {len(configs)} configurations with {args.n_processes} processes")

    with Pool(args.n_processes) as pool:
        results = pool.map(run_single_config, configs)

    # Collect all metrics into summary
    outdir = PROJECT_ROOT / "results" / "harmony_mapping"
    metric_files = list(outdir.glob("*_metrics.csv"))
    all_metrics = []
    for mf in metric_files:
        m = pd.read_csv(mf, index_col=0, header=None).squeeze()
        all_metrics.append(m)

    if all_metrics:
        summary = pd.DataFrame(all_metrics)
        summary = summary.sort_values("accuracy", ascending=False)
        summary.to_csv(outdir / "sweep_summary.csv", index=False)
        log.info(f"Sweep summary saved ({len(summary)} configs)")
        log.info(f"Best config:\n{summary.iloc[0]}")
    else:
        log.warning("No metric files found")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the sweep (overnight)**

Run:
```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic python cge_subtype/scripts/09_harmony_sweep.py --dataset M1 --n-processes 10
```

Expected: ~24 configurations run, results in `results/harmony_mapping/sweep_summary.csv`. Review best configuration by CGE accuracy.

- [ ] **Step 4: Commit Module B scripts**

```bash
git add cge_subtype/scripts/08_harmony_mapping.py cge_subtype/scripts/09_harmony_sweep.py
git commit -m "Add Module B scripts: Harmony mapping and hyperparameter sweep"
```

---

## Chunk 3: Modules C + D — Concordance and Ephys Convergence

### Task 5: Module C — Path Concordance

**Files:**
- Create: `cge_subtype/src/concordance.py`
- Create: `cge_subtype/tests/test_concordance.py`
- Create: `cge_subtype/notebooks/Path_Concordance.py`

- [ ] **Step 1: Write concordance tests**

```python
# cge_subtype/tests/test_concordance.py
import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.concordance import (
    build_indirect_path,
    compute_concordance,
    compute_cohens_kappa,
)


def test_build_indirect_path():
    """Chain mouse patch-seq -> mouse cluster -> human cluster."""
    patchseq_to_mouse = pd.DataFrame({
        "cell_id": ["c1", "c2", "c3"],
        "mouse_cluster": ["mA", "mB", "mA"],
    })
    mouse_to_human = pd.DataFrame({
        "mouse_cluster": ["mA", "mB"],
        "human_cluster": ["hX", "hY"],
        "correlation": [0.9, 0.85],
    })
    result = build_indirect_path(patchseq_to_mouse, mouse_to_human)
    assert len(result) == 3
    assert result[result["cell_id"] == "c1"]["indirect_human_cluster"].values[0] == "hX"
    assert result[result["cell_id"] == "c2"]["indirect_human_cluster"].values[0] == "hY"


def test_concordance_perfect():
    """Perfect agreement should give 1.0."""
    direct = pd.Series(["A", "B", "C", "A"])
    indirect = pd.Series(["A", "B", "C", "A"])
    assert compute_concordance(direct, indirect) == 1.0


def test_concordance_none():
    """No agreement should give 0.0."""
    direct = pd.Series(["A", "B", "C"])
    indirect = pd.Series(["X", "Y", "Z"])
    assert compute_concordance(direct, indirect) == 0.0


def test_cohens_kappa_perfect():
    """Perfect agreement kappa = 1.0."""
    a = pd.Series(["A", "B", "A", "B"])
    b = pd.Series(["A", "B", "A", "B"])
    assert compute_cohens_kappa(a, b) == pytest.approx(1.0)
```

- [ ] **Step 2: Run tests to verify failure**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_concordance.py -v`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement concordance.py**

```python
# cge_subtype/src/concordance.py
"""Path concordance validation: chaining, agreement metrics, Cohen's kappa."""

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score


def build_indirect_path(patchseq_to_mouse, mouse_to_human):
    """Chain mouse patch-seq -> mouse cluster -> human cluster.

    Parameters
    ----------
    patchseq_to_mouse : pd.DataFrame
        Columns: cell_id, mouse_cluster (+ optional confidence).
    mouse_to_human : pd.DataFrame
        Columns: mouse_cluster, human_cluster (+ optional correlation).

    Returns
    -------
    pd.DataFrame
        Original patchseq_to_mouse with added indirect_human_cluster column.
    """
    result = patchseq_to_mouse.merge(
        mouse_to_human[["mouse_cluster", "human_cluster"]],
        on="mouse_cluster",
        how="left",
    )
    result = result.rename(columns={"human_cluster": "indirect_human_cluster"})
    return result


def compute_concordance(direct_labels, indirect_labels):
    """Fraction of cells where direct and indirect paths agree.

    Parameters
    ----------
    direct_labels : pd.Series
        Direct path cluster assignments.
    indirect_labels : pd.Series
        Indirect path cluster assignments.

    Returns
    -------
    float
        Agreement rate (0 to 1).
    """
    valid = ~(direct_labels.isna() | indirect_labels.isna())
    if valid.sum() == 0:
        return np.nan
    return (direct_labels[valid] == indirect_labels[valid]).mean()


def compute_cohens_kappa(labels_a, labels_b):
    """Compute Cohen's kappa (chance-corrected agreement).

    Parameters
    ----------
    labels_a, labels_b : pd.Series
        Two sets of categorical labels.

    Returns
    -------
    float
        Cohen's kappa (-1 to 1).
    """
    valid = ~(labels_a.isna() | labels_b.isna())
    return cohen_kappa_score(labels_a[valid], labels_b[valid])


def concordance_at_levels(direct_clusters, indirect_clusters, cluster_to_subclass, cluster_to_supercluster):
    """Compute concordance at cluster, subclass, and supercluster levels.

    Parameters
    ----------
    direct_clusters : pd.Series
        Direct path cluster assignments.
    indirect_clusters : pd.Series
        Indirect path cluster assignments.
    cluster_to_subclass : dict
        Mapping from cluster ID to subclass label.
    cluster_to_supercluster : dict
        Mapping from cluster ID to supercluster label.

    Returns
    -------
    dict
        Concordance and kappa at each level.
    """
    results = {}

    # Cluster level
    results["cluster_concordance"] = compute_concordance(direct_clusters, indirect_clusters)
    results["cluster_kappa"] = compute_cohens_kappa(direct_clusters, indirect_clusters)

    # Subclass level
    direct_sub = direct_clusters.map(cluster_to_subclass)
    indirect_sub = indirect_clusters.map(cluster_to_subclass)
    results["subclass_concordance"] = compute_concordance(direct_sub, indirect_sub)
    results["subclass_kappa"] = compute_cohens_kappa(direct_sub, indirect_sub)

    # Supercluster level
    direct_super = direct_clusters.map(cluster_to_supercluster)
    indirect_super = indirect_clusters.map(cluster_to_supercluster)
    results["supercluster_concordance"] = compute_concordance(direct_super, indirect_super)
    results["supercluster_kappa"] = compute_cohens_kappa(direct_super, indirect_super)

    return results
```

- [ ] **Step 4: Run tests to verify pass**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_concordance.py -v`
Expected: all PASS

- [ ] **Step 5: Write Path_Concordance.py notebook**

Create `cge_subtype/notebooks/Path_Concordance.py` as a jupytext percent-format notebook:

```python
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Module C: Path Concordance Validation
#
# Compare direct (mouse patch-seq -> human SC) vs indirect
# (mouse patch-seq -> mouse SC -> human SC) mapping paths.

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from cge_subtype.src.concordance import (
    build_indirect_path,
    compute_concordance,
    concordance_at_levels,
)

PROJECT_ROOT = Path("..").resolve()

# %%
# Load Module A output: mouse cluster -> human cluster correspondence
rbh = pd.read_csv(PROJECT_ROOT / "results" / "cluster_bridge" / "rbh_resolved_pairs.csv")
print(f"RBH pairs: {len(rbh)}")

# %%
# Load Module B output: mouse patch-seq -> mouse cluster (best Harmony config)
# TODO: Update path to best sweep result
harmony_pred = pd.read_csv(
    sorted((PROJECT_ROOT / "results" / "harmony_mapping").glob("M1*_predictions.csv"))[-1]
)
print(f"Harmony predictions: {len(harmony_pred)}")

# %%
# Load existing direct mapping: mouse patch-seq -> human SC
direct_mapping = pd.read_csv(PROJECT_ROOT / "results" / "patchseq_mapping_results.csv")
print(f"Direct mapping: {len(direct_mapping)}")

# %%
# Build indirect path
indirect = build_indirect_path(
    harmony_pred[["cell_id", "predicted_cluster"]].rename(
        columns={"predicted_cluster": "mouse_cluster"}
    ),
    rbh[["mouse_cluster", "human_cluster"]],
)
print(f"Indirect path: {len(indirect)} cells, {indirect['indirect_human_cluster'].notna().sum()} with assignments")

# %%
# Merge direct and indirect
merged = direct_mapping.merge(indirect[["cell_id", "indirect_human_cluster"]], on="cell_id", how="inner")
print(f"Cells with both paths: {len(merged)}")

# %%
# Concordance at 3 levels
# TODO: Load cluster -> subclass and cluster -> supercluster mappings from Siletti metadata
# cluster_to_subclass = ...
# cluster_to_supercluster = ...
# results = concordance_at_levels(merged["direct_cluster"], merged["indirect_human_cluster"],
#                                  cluster_to_subclass, cluster_to_supercluster)

# For now, compute cluster-level concordance
cluster_conc = compute_concordance(
    merged["predicted_human_cluster"],  # direct
    merged["indirect_human_cluster"],   # indirect
)
print(f"Cluster-level concordance: {cluster_conc:.3f}")

# %%
# CCKBC-focused analysis
# TODO: Filter to CCKBC cells (Sncg subclass + Vip Sncg type)
# and check if both paths agree on human CGE clusters

# %%
# Visualization: Sankey diagram (direct vs indirect)
# TODO: Implement Sankey using matplotlib or plotly

# %%
# Save results
# TODO: Save concordance table and figures
```

- [ ] **Step 6: Sync notebook and commit**

```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic jupytext --to notebook cge_subtype/notebooks/Path_Concordance.py
git add cge_subtype/src/concordance.py cge_subtype/tests/test_concordance.py \
    cge_subtype/notebooks/Path_Concordance.py cge_subtype/notebooks/Path_Concordance.ipynb
git commit -m "Add Module C: path concordance validation (library, tests, notebook)"
```

---

### Task 6: Module D — Ephys Convergence

**Files:**
- Create: `cge_subtype/src/ephys_harmonization.py`
- Create: `cge_subtype/tests/test_ephys_harmonization.py`
- Create: `cge_subtype/scripts/10_aggregate_ephys.py`
- Create: `cge_subtype/notebooks/Ephys_Convergence.py`

- [ ] **Step 1: Write ephys harmonization tests**

```python
# cge_subtype/tests/test_ephys_harmonization.py
import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.ephys_harmonization import (
    aggregate_cell_features,
    select_representative_sweep,
    zscore_within_species,
    permutation_test_cluster_similarity,
)


def test_select_representative_sweep():
    """Select sweep with firing rate closest to target range."""
    sweeps = pd.DataFrame({
        "spike_frequency_Hz": [0, 0, 10, 25, 60, 100],
        "avg_peak_voltage_mV": [0, 0, 30, 35, 40, 45],
    })
    idx = select_representative_sweep(sweeps, target_min=5, target_max=60)
    assert idx in [2, 3, 4]  # any sweep in 5-60 Hz range


def test_zscore_within_species():
    """Z-scoring should produce mean=0, std=1 within each species."""
    df = pd.DataFrame({
        "feat1": [1, 2, 3, 10, 20, 30],
        "feat2": [4, 5, 6, 40, 50, 60],
    })
    species = pd.Series(["mouse", "mouse", "mouse", "human", "human", "human"])
    result = zscore_within_species(df, species)
    # Check mouse group
    mouse_vals = result.loc[species == "mouse", "feat1"]
    assert abs(mouse_vals.mean()) < 1e-10
    assert abs(mouse_vals.std(ddof=0) - 1.0) < 0.1 or len(mouse_vals) < 4


def test_permutation_test_smoke():
    """Permutation test should return a p-value between 0 and 1."""
    np.random.seed(42)
    features = pd.DataFrame(np.random.randn(20, 5))
    cluster_labels = pd.Series(["A"] * 5 + ["B"] * 5 + ["A"] * 5 + ["B"] * 5)
    species = pd.Series(["mouse"] * 10 + ["human"] * 10)
    p = permutation_test_cluster_similarity(features, cluster_labels, species, n_perm=100)
    assert 0 <= p <= 1
```

- [ ] **Step 2: Run tests to verify failure**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_ephys_harmonization.py -v`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement ephys_harmonization.py**

```python
# cge_subtype/src/ephys_harmonization.py
"""Ephys feature harmonization and cross-species comparison."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

log = logging.getLogger(__name__)

# Columns to keep from EphysSumStats analysis.csv (biologically meaningful)
BIO_FEATURES = [
    "avg_peak_voltage_mV",
    "avg_threshold_voltage_mV",
    "avg_threshold_to_peak_mV",
    "avg_fast_trough_at_hyperpolarization_mV",
    "avg_height_mV",
    "avg_half_height_mV",
    "avg_upstroke_to_peak_ms",
    "avg_peak_to_downstroke_ms",
    "avg_ap_width_ms",
    "avg_threshold_fast_trough_width_ms",
    "avg_upstroke_downstroke_ratio",
    "avg_upstroke_mVms",
    "avg_downstroke_mVms",
    "spike_frequency_Hz",
    "resting_vm_mean_mV",
    "filtered_resting_vm_mean_mV",
    "current_threshold_pA",
]

# ISI bins (shared: bins 0-11)
ISI_FEATURES = [f"bin_{i}_isi_ms" for i in range(12)]
CV_ISI_FEATURES = [f"bin_{i}_cv_isi" for i in range(12)]

# Features to log-transform before z-scoring
LOG_FEATURES = ["spike_frequency_Hz", "current_threshold_pA"] + ISI_FEATURES


def select_representative_sweep(sweeps_df, target_min=5, target_max=60):
    """Select a representative spiking sweep.

    Parameters
    ----------
    sweeps_df : pd.DataFrame
        Per-sweep features (one row per sweep).
    target_min, target_max : float
        Target firing rate range (Hz).

    Returns
    -------
    int
        Index of the selected sweep.
    """
    spiking = sweeps_df[sweeps_df["spike_frequency_Hz"] > 0]
    if len(spiking) == 0:
        return sweeps_df.index[0]

    in_range = spiking[
        (spiking["spike_frequency_Hz"] >= target_min)
        & (spiking["spike_frequency_Hz"] <= target_max)
    ]

    if len(in_range) > 0:
        # Pick sweep closest to 20 Hz (middle of range)
        return in_range.iloc[(in_range["spike_frequency_Hz"] - 20).abs().argmin()].name
    else:
        # Fallback: lowest-current spiking sweep
        return spiking.iloc[0].name


def aggregate_cell_features(processed_dir, species_name):
    """Aggregate per-cell features from EphysSumStats analysis.csv files.

    Parameters
    ----------
    processed_dir : str or Path
        Path to /mnt/data0/DANDI/Processed/{dataset_id}/
    species_name : str
        "mouse" or "human" (for tagging).

    Returns
    -------
    pd.DataFrame
        One row per cell, columns = bio features + ISI bins.
    """
    processed_dir = Path(processed_dir)
    all_features = []

    # Find all analysis.csv files
    analysis_files = list(processed_dir.rglob("analysis.csv"))
    log.info(f"Found {len(analysis_files)} analysis.csv files in {processed_dir}")

    feature_cols = BIO_FEATURES + ISI_FEATURES + CV_ISI_FEATURES
    available_cols = None

    for af in analysis_files:
        cell_id = af.parent.name
        try:
            sweeps = pd.read_csv(af)
        except Exception as e:
            log.warning(f"Failed to read {af}: {e}")
            continue

        # Filter to available columns
        if available_cols is None:
            available_cols = [c for c in feature_cols if c in sweeps.columns]
            log.info(f"Available features: {len(available_cols)}/{len(feature_cols)}")

        # Select representative sweep
        if "spike_frequency_Hz" in sweeps.columns:
            sweep_idx = select_representative_sweep(sweeps)
            row = sweeps.loc[sweep_idx, available_cols].to_dict()
        else:
            row = sweeps.iloc[0][available_cols].to_dict() if len(available_cols) > 0 else {}

        row["cell_id"] = cell_id
        row["species"] = species_name
        all_features.append(row)

    result = pd.DataFrame(all_features)
    log.info(f"Aggregated {len(result)} cells for {species_name}")
    return result


def zscore_within_species(features_df, species_labels):
    """Z-score features within each species.

    Parameters
    ----------
    features_df : pd.DataFrame
        Cells x features (numeric only).
    species_labels : pd.Series
        "mouse" or "human" per cell.

    Returns
    -------
    pd.DataFrame
        Z-scored features.
    """
    result = features_df.copy()
    for sp in species_labels.unique():
        mask = species_labels == sp
        sp_data = result.loc[mask]
        result.loc[mask] = (sp_data - sp_data.mean()) / (sp_data.std() + 1e-10)
    return result


def permutation_test_cluster_similarity(features_df, cluster_labels, species_labels,
                                         n_perm=1000, seed=42):
    """Test if same-cluster cross-species pairs are more similar than expected.

    For each cluster with cells from both species, compute mean Euclidean distance
    between species. Compare to null (shuffled cluster labels).

    Parameters
    ----------
    features_df : pd.DataFrame
        Z-scored features.
    cluster_labels : pd.Series
        Cluster assignment per cell.
    species_labels : pd.Series
        Species per cell.
    n_perm : int
        Number of permutations.
    seed : int
        Random seed.

    Returns
    -------
    float
        P-value.
    """
    rng = np.random.default_rng(seed)

    def mean_cross_species_distance(feat, clusters, species):
        distances = []
        for clust in clusters.unique():
            mask = clusters == clust
            mouse_cells = feat.loc[mask & (species == "mouse")]
            human_cells = feat.loc[mask & (species == "human")]
            if len(mouse_cells) == 0 or len(human_cells) == 0:
                continue
            mouse_mean = mouse_cells.mean().values
            human_mean = human_cells.mean().values
            dist = np.sqrt(((mouse_mean - human_mean) ** 2).sum())
            distances.append(dist)
        return np.mean(distances) if distances else np.nan

    observed = mean_cross_species_distance(features_df, cluster_labels, species_labels)
    if np.isnan(observed):
        return np.nan

    null_dists = []
    for _ in range(n_perm):
        shuffled = cluster_labels.copy()
        shuffled.iloc[:] = rng.permutation(shuffled.values)
        null_val = mean_cross_species_distance(features_df, shuffled, species_labels)
        if not np.isnan(null_val):
            null_dists.append(null_val)

    if len(null_dists) == 0:
        return np.nan

    # One-sided: observed should be SMALLER than null (more similar)
    p_value = np.mean([n <= observed for n in null_dists])
    return p_value
```

- [ ] **Step 4: Run tests to verify pass**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_ephys_harmonization.py -v`
Expected: all PASS

- [ ] **Step 5: Write 10_aggregate_ephys.py**

```python
#!/usr/bin/env python
"""Aggregate EphysSumStats features for mouse and human patch-seq cells.

Usage:
    conda activate gencic
    python cge_subtype/scripts/10_aggregate_ephys.py
"""

import logging
from pathlib import Path

from cge_subtype.src.ephys_harmonization import aggregate_cell_features, LOG_FEATURES

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main():
    outdir = PROJECT_ROOT / "results" / "ephys_convergence"
    outdir.mkdir(parents=True, exist_ok=True)

    # Aggregate mouse (DANDI 000008)
    log.info("Aggregating mouse M1 patch-seq ephys features")
    mouse_feat = aggregate_cell_features("/mnt/data0/DANDI/Processed/000008/", "mouse")
    mouse_feat.to_csv(outdir / "mouse_ephys_features.csv", index=False)
    log.info(f"Mouse: {len(mouse_feat)} cells")

    # Aggregate human (DANDI 000636)
    log.info("Aggregating human cortical GABA ephys features")
    human_feat = aggregate_cell_features("/mnt/data0/DANDI/Processed/000636/", "human")
    human_feat.to_csv(outdir / "human_ephys_features.csv", index=False)
    log.info(f"Human: {len(human_feat)} cells")

    # Combine
    combined = pd.concat([mouse_feat, human_feat], ignore_index=True)

    # Log-transform skewed features
    for col in LOG_FEATURES:
        if col in combined.columns:
            combined[col] = np.log1p(combined[col].clip(lower=0))

    combined.to_csv(outdir / "combined_ephys_features.csv", index=False)
    log.info(f"Combined: {len(combined)} cells, {combined.shape[1]} columns")


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Write Ephys_Convergence.py notebook**

Create `cge_subtype/notebooks/Ephys_Convergence.py` as a jupytext percent-format notebook:

```python
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Module D: Cross-Species Ephys Convergence
#
# Validate that mouse and human patch-seq cells mapping to the same
# human cluster share electrophysiological features.

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

from cge_subtype.src.ephys_harmonization import (
    zscore_within_species,
    permutation_test_cluster_similarity,
    BIO_FEATURES,
)

PROJECT_ROOT = Path("..").resolve()

# %%
# Load aggregated ephys features
combined = pd.read_csv(PROJECT_ROOT / "results" / "ephys_convergence" / "combined_ephys_features.csv")
print(f"Combined: {combined.shape}")
print(f"Species: {combined['species'].value_counts().to_dict()}")

# %%
# Load cluster assignments for both species
# Mouse: from existing direct mapping OR Module B harmony results
# Human: from existing scANVI mapping
mouse_mapping = pd.read_csv(PROJECT_ROOT / "results" / "patchseq_mapping_results.csv")
human_mapping = pd.read_csv(PROJECT_ROOT / "results" / "human_scvi_mapping_results.csv")

# %%
# Merge cluster assignments with ephys features
# TODO: Match cell IDs between ephys and mapping results

# %%
# Tier 1: Within-species z-scoring
feature_cols = [c for c in BIO_FEATURES if c in combined.columns]
z_features = zscore_within_species(combined[feature_cols], combined["species"])
print(f"Z-scored features: {z_features.shape}")

# %%
# Identify shared clusters (both mouse and human cells)
# TODO: After merging, find clusters with cells from both species

# %%
# Global permutation test
# TODO: Run after cluster assignments are merged
# p = permutation_test_cluster_similarity(z_features, cluster_labels, combined["species"])
# print(f"Global permutation test p-value: {p:.4f}")

# %%
# Per-cluster ephys comparison
# TODO: For each shared cluster, compare mouse vs human profiles

# %%
# CCKBC-specific analysis
# TODO: Compare ephys profiles of putative CCKBC vs ISI clusters

# %%
# Tier 2: ComBat sensitivity (optional)
# TODO: Install neuroCombat and run if Tier 1 results are interesting

# %%
# Save results and figures
```

- [ ] **Step 7: Sync notebooks and commit**

```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic jupytext --to notebook cge_subtype/notebooks/Ephys_Convergence.py
git add cge_subtype/src/ephys_harmonization.py cge_subtype/tests/test_ephys_harmonization.py \
    cge_subtype/scripts/10_aggregate_ephys.py \
    cge_subtype/notebooks/Ephys_Convergence.py cge_subtype/notebooks/Ephys_Convergence.ipynb
git commit -m "Add Modules C+D: concordance validation and ephys convergence (libraries, tests, scripts, notebooks)"
```

---

## Chunk 4: Module E — Updated 22q Bias + Final Integration

### Task 7: Module E — Multi-Modal Classification + Updated Bias

**Files:**
- Create: `cge_subtype/src/multimodal_classification.py`
- Create: `cge_subtype/tests/test_multimodal_classification.py`
- Create: `cge_subtype/notebooks/Updated_22q_Bias.py`

- [ ] **Step 1: Write multimodal classification tests**

```python
# cge_subtype/tests/test_multimodal_classification.py
import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.multimodal_classification import (
    compute_cckbc_confidence,
    classify_clusters,
)


def test_confidence_all_evidence():
    """Cluster with all 4 lines of evidence = high confidence."""
    evidence = pd.DataFrame({
        "cluster_id": [276],
        "module_a_is_sncg": [True],
        "module_c_direct_cckbc": [True],
        "module_c_indirect_cckbc": [True],
        "module_d_fast_spiking": [True],
        "marker_cck_positive": [True],
        "marker_sncg_positive": [True],
    })
    result = compute_cckbc_confidence(evidence)
    assert result.loc[0, "cckbc_confidence"] >= 3
    assert result.loc[0, "cckbc_tier"] == "high-confidence CCKBC"


def test_confidence_no_evidence():
    """Cluster with no evidence = ISI VIP."""
    evidence = pd.DataFrame({
        "cluster_id": [290],
        "module_a_is_sncg": [False],
        "module_c_direct_cckbc": [False],
        "module_c_indirect_cckbc": [False],
        "module_d_fast_spiking": [False],
        "marker_cck_positive": [False],
        "marker_sncg_positive": [False],
    })
    result = compute_cckbc_confidence(evidence)
    assert result.loc[0, "cckbc_confidence"] == 0
    assert result.loc[0, "cckbc_tier"] == "ISI VIP"


def test_classify_clusters():
    """Classification tiers: 3+ = high, 1-2 = tentative, 0 = ISI."""
    evidence = pd.DataFrame({
        "cluster_id": [276, 279, 290],
        "module_a_is_sncg": [True, True, False],
        "module_c_direct_cckbc": [True, False, False],
        "module_c_indirect_cckbc": [True, False, False],
        "module_d_fast_spiking": [True, True, False],
        "marker_cck_positive": [True, False, False],
        "marker_sncg_positive": [True, False, False],
    })
    result = classify_clusters(evidence)
    tiers = result.set_index("cluster_id")["cckbc_tier"]
    assert tiers[276] == "high-confidence CCKBC"
    assert tiers[279] == "tentative CCKBC"
    assert tiers[290] == "ISI VIP"
```

- [ ] **Step 2: Run tests to verify failure**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_multimodal_classification.py -v`
Expected: FAIL

- [ ] **Step 3: Implement multimodal_classification.py**

```python
# cge_subtype/src/multimodal_classification.py
"""Multi-modal CCKBC confidence scoring and classification."""

import numpy as np
import pandas as pd


# Evidence columns that count toward CCKBC confidence
EVIDENCE_COLUMNS = [
    "module_a_is_sncg",         # Module A: mouse cluster is Sncg/CCKBC subclass
    "module_c_direct_cckbc",    # Module C: direct path maps CCKBC cells here
    "module_c_indirect_cckbc",  # Module C: indirect path maps CCKBC cells here
    "module_d_fast_spiking",    # Module D: ephys resembles fast-spiking CCKBC
    "marker_cck_positive",      # Existing: CCK+ marker expression
    "marker_sncg_positive",     # Existing: SNCG+ marker expression
]


def compute_cckbc_confidence(evidence_df):
    """Compute CCKBC confidence score per cluster.

    Parameters
    ----------
    evidence_df : pd.DataFrame
        One row per cluster. Must contain EVIDENCE_COLUMNS as boolean columns.

    Returns
    -------
    pd.DataFrame
        Original df with added cckbc_confidence (int) and cckbc_tier (str) columns.
    """
    result = evidence_df.copy()
    available_cols = [c for c in EVIDENCE_COLUMNS if c in result.columns]
    result["cckbc_confidence"] = result[available_cols].sum(axis=1).astype(int)

    result["cckbc_tier"] = "ISI VIP"
    result.loc[result["cckbc_confidence"] >= 1, "cckbc_tier"] = "tentative CCKBC"
    result.loc[result["cckbc_confidence"] >= 3, "cckbc_tier"] = "high-confidence CCKBC"

    return result


def classify_clusters(evidence_df):
    """Classify CGE clusters into CCKBC tiers.

    Convenience wrapper around compute_cckbc_confidence.
    """
    return compute_cckbc_confidence(evidence_df)
```

- [ ] **Step 4: Run tests to verify pass**

Run: `cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP && conda run -n gencic python -m pytest cge_subtype/tests/test_multimodal_classification.py -v`
Expected: all PASS

- [ ] **Step 5: Write Updated_22q_Bias.py notebook**

Create `cge_subtype/notebooks/Updated_22q_Bias.py`:

```python
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Module E: Updated 22q Bias with Multi-Modal CCKBC Classification
#
# Recalculate 22q bias using CCKBC/ISI groupings supported by
# converging evidence from Modules A-D.

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

import sys
sys.path.insert(0, str(Path("../..").resolve()))
from src.CellType_PSY import HumanCT_AvgZ_Weighted, AddPvalue_optimized

from cge_subtype.src.multimodal_classification import classify_clusters

PROJECT_ROOT = Path("..").resolve()
MAIN_ROOT = PROJECT_ROOT.parent

# %%
# Load expression specificity matrix
spec_mat = pd.read_csv(
    MAIN_ROOT / "dat" / "ExpMats" / "HumanCT.spec.csv",
    index_col=0,
)
print(f"Specificity matrix: {spec_mat.shape}")

# %%
# Load 22q gene sets
gene_set_names = ["22q_del", "22q_mouse", "22q_HighExp", "22q_DEG_d75"]
gene_sets = {}
for name in gene_set_names:
    gw_path = MAIN_ROOT / "dat" / "GeneWeights" / f"{name}.gw.csv"
    if gw_path.exists():
        gene_sets[name] = pd.read_csv(gw_path, header=None, names=["Entrez_ID", "Weight"])
        print(f"{name}: {len(gene_sets[name])} genes")

# %%
# Assemble evidence table from all modules
# TODO: Load outputs from Modules A, C, D and existing marker data
# evidence = pd.DataFrame({
#     "cluster_id": range(276, 297),
#     "module_a_is_sncg": [...],
#     "module_c_direct_cckbc": [...],
#     "module_c_indirect_cckbc": [...],
#     "module_d_fast_spiking": [...],
#     "marker_cck_positive": [...],
#     "marker_sncg_positive": [...],
# })

# %%
# Classify clusters
# classified = classify_clusters(evidence)
# print(classified[["cluster_id", "cckbc_confidence", "cckbc_tier"]])

# %%
# Calculate 22q bias per cluster
# TODO: For each gene set, compute bias using HumanCT_AvgZ_Weighted()
# Focus on CGE clusters 276-296

# %%
# Group comparisons
# Primary: high-confidence CCKBC vs ISI VIP (Mann-Whitney U)
# Secondary: 3-way split by VIP expression
# Correlation: CCKBC confidence vs 22q bias (Spearman)

# %%
# Sensitivity analysis
# - Transcriptomic evidence only vs ephys only
# - Across all 4 gene sets
# - Bootstrap CIs

# %%
# Visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.patch.set_alpha(0)

# Panel A: Updated boxplot (CCKBC vs ISI)
# axes[0]: boxplot of 22q bias by tier

# Panel B: Scatter (confidence score vs bias)
# axes[1]: scatter with Spearman r annotation

# Panel C: Evidence summary panel
# axes[2]: heatmap of evidence per cluster

for ax in axes:
    ax.patch.set_alpha(0)

plt.tight_layout()
# plt.savefig(PROJECT_ROOT / "results" / "updated_22q_bias_multimodal.pdf", transparent=True)

# %%
# Save results
# TODO: Save updated bias table, classification, statistics
```

- [ ] **Step 6: Sync notebooks and commit**

```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
conda run -n gencic jupytext --to notebook cge_subtype/notebooks/Updated_22q_Bias.py
git add cge_subtype/src/multimodal_classification.py cge_subtype/tests/test_multimodal_classification.py \
    cge_subtype/notebooks/Updated_22q_Bias.py cge_subtype/notebooks/Updated_22q_Bias.ipynb
git commit -m "Add Module E: multi-modal CCKBC classification and updated 22q bias (library, tests, notebook)"
```

---

## Execution Notes

### Parallelization Strategy

Modules A and B are independent and should be run in parallel:
- **Module A** (Tasks 1-2): Pseudobulk + cluster bridge — ~2-4 hours
- **Module B** (Tasks 3-4): Harmony sweep — ~8-12 hours (overnight)
- **Module D data prep** (Task 6, steps 5-6): Ephys aggregation — ~2 hours (can overlap with A/B)

After A + B complete:
- **Module C** (Task 5): Concordance — ~1 hour
- **Module D analysis** (Task 6, notebook): Ephys convergence — ~1 hour
- **Module E** (Task 7): Updated bias — ~1 hour

### Data Exploration Checkpoints

Several steps require manual inspection of intermediate results:

1. **After Task 2 Step 2**: Verify atlas h5ad structure (column names for clusters, gene symbols)
2. **After Task 2 Step 5**: Inspect RBH results — how many pairs resolved? Are CGE clusters unresolved?
3. **After Task 4 Step 3**: Review sweep summary — best Harmony config by CGE accuracy
4. **After Task 6 Step 5 aggregation**: Check ephys feature distributions, missing values
5. **After Task 7**: Review final CCKBC classification and bias comparison

### Known Risks

1. **RAM during pseudobulk**: Monitor with `htop`. If backed mode still exceeds limits, reduce chunk_size or process one supercluster file at a time.
2. **MetaNeighbor Pass 2**: Requires subsampled cells from unresolved clusters. The script for this is not fully implemented — it needs cell-level loading from the atlas for specific clusters only.
3. **Ortholog column names**: The ortholog mapping CSV column names need verification — `load_ortholog_mapping()` assumes columns 0 and 1 are mouse and human genes respectively.
4. **Cell ID matching**: Merging ephys features with transcriptomic mappings requires matching cell IDs across different naming conventions (DANDI stem vs cell_specimen_id vs cell_label).
