"""Harmony-based mapping of patch-seq cells to reference atlas."""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import harmonypy
from sklearn.neighbors import NearestNeighbors


# ---------------------------------------------------------------------------
# Non-neuronal cell type patterns used by evaluate_mapping
# ---------------------------------------------------------------------------

_NON_NEURONAL_TYPES = {"Astro", "Oligo", "OPC", "Micro", "Endo", "VLMC", "Peri"}
_NON_NEURONAL_PREFIX = "NN "


def subsample_reference(
    adata: ad.AnnData,
    cluster_col: str,
    max_per_cluster: int = 200,
    seed: int = 42,
) -> ad.AnnData:
    """Subsample a reference atlas to at most *max_per_cluster* cells per cluster.

    Parameters
    ----------
    adata:
        Reference AnnData object.  Must have ``adata.obs[cluster_col]``.
    cluster_col:
        Column in ``adata.obs`` that contains cluster labels.
    max_per_cluster:
        Maximum number of cells to retain per cluster.
    seed:
        Random seed for reproducibility (uses ``np.random.default_rng``).

    Returns
    -------
    AnnData
        A *copy* of the subsampled data.
    """
    rng = np.random.default_rng(seed)
    clusters = adata.obs[cluster_col].values
    keep_indices: list[int] = []

    for cluster_label in np.unique(clusters):
        idx = np.where(clusters == cluster_label)[0]
        if len(idx) > max_per_cluster:
            chosen = rng.choice(idx, size=max_per_cluster, replace=False)
        else:
            chosen = idx
        keep_indices.extend(chosen.tolist())

    # Sort to preserve original ordering
    keep_indices = sorted(keep_indices)
    return adata[keep_indices].copy()


def run_harmony_mapping(
    ref_adata: ad.AnnData,
    query_adata: ad.AnnData,
    cluster_col: str,
    n_pcs: int = 50,
    n_hvgs: int = 3000,
    theta: float = 2.0,
    n_neighbors: int = 30,
    batch_key: str = "technology",
) -> pd.DataFrame:
    """Map query cells onto a reference atlas using Harmony integration + kNN transfer.

    Steps
    -----
    1. Tag each dataset with "reference" / "query" in *batch_key*.
    2. Concatenate ref + query into a single AnnData.
    3. Normalize (``normalize_total`` to 1e4 + ``log1p``).
    4. Select highly variable genes (``n_top_genes=n_hvgs``, ``batch_key``).
    5. Scale and run PCA (``n_comps=n_pcs``).
    6. Run Harmony on the PCA embedding to remove batch effects.
    7. kNN label transfer: for each query cell, find the *n_neighbors* nearest
       reference cells in Harmony space and assign the majority-vote cluster.

    Parameters
    ----------
    ref_adata:
        Annotated reference atlas.  Must have ``obs[cluster_col]``.
    query_adata:
        Query cells (e.g. patch-seq dataset) to be mapped.
    cluster_col:
        Column in ``ref_adata.obs`` holding the cluster labels to transfer.
    n_pcs:
        Number of principal components for PCA.
    n_hvgs:
        Number of highly variable genes to select.
    theta:
        Harmony diversity penalty (higher = stronger batch correction).
    n_neighbors:
        Number of nearest reference neighbors used for label transfer.
    batch_key:
        Column name used to tag ref vs. query batches.

    Returns
    -------
    DataFrame with columns:
        - ``cell_id``: query cell index.
        - ``predicted_cluster``: majority-vote label from reference neighbors.
        - ``confidence``: fraction of the *n_neighbors* agreeing on the label.
        - ``n_neighbors_agreeing``: absolute count of agreeing neighbors.
    """
    # ------------------------------------------------------------------
    # 1. Tag batches
    # ------------------------------------------------------------------
    ref = ref_adata.copy()
    query = query_adata.copy()
    ref.obs[batch_key] = "reference"
    query.obs[batch_key] = "query"

    # ------------------------------------------------------------------
    # 2. Concatenate
    # ------------------------------------------------------------------
    combined = ad.concat(
        [ref, query],
        label=batch_key,
        keys=["reference", "query"],
        merge="same",
    )
    # The concat label column may be redundant — ensure batch_key is intact
    combined.obs[batch_key] = (
        combined.obs[batch_key]
        if batch_key in combined.obs.columns
        else combined.obs["batch"]
    )

    # ------------------------------------------------------------------
    # 3. Normalize
    # ------------------------------------------------------------------
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    # ------------------------------------------------------------------
    # 4. HVG selection
    # ------------------------------------------------------------------
    sc.pp.highly_variable_genes(
        combined,
        n_top_genes=n_hvgs,
        batch_key=batch_key,
        subset=False,
    )
    combined = combined[:, combined.var["highly_variable"]].copy()

    # ------------------------------------------------------------------
    # 5. Scale + PCA
    # ------------------------------------------------------------------
    sc.pp.scale(combined)
    sc.tl.pca(combined, n_comps=n_pcs)

    # ------------------------------------------------------------------
    # 6. Harmony
    # ------------------------------------------------------------------
    pca_embedding = combined.obsm["X_pca"]
    meta = combined.obs[[batch_key]].copy()

    ho = harmonypy.run_harmony(
        pca_embedding,
        meta,
        batch_key,
        theta=theta,
    )
    combined.obsm["X_harmony"] = ho.Z_corr.T  # harmony returns (n_pcs, n_cells)

    # ------------------------------------------------------------------
    # 7. kNN label transfer
    # ------------------------------------------------------------------
    is_ref = combined.obs[batch_key] == "reference"
    is_query = combined.obs[batch_key] == "query"

    harmony_emb = combined.obsm["X_harmony"]
    ref_emb = harmony_emb[is_ref.values]
    query_emb = harmony_emb[is_query.values]

    # Cluster labels for reference cells (aligned with ref_emb rows)
    # Use the original cluster_col from ref_adata; the obs index in combined
    # for reference cells preserves ref_adata obs names.
    ref_obs = combined.obs.loc[is_ref]
    ref_labels = ref_obs[cluster_col].values

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean", n_jobs=1)
    nn.fit(ref_emb)
    _, indices = nn.kneighbors(query_emb)

    # Majority vote
    neighbor_labels = ref_labels[indices]  # shape (n_query, n_neighbors)

    predicted_clusters = []
    confidences = []
    n_agreeing_list = []

    for row in neighbor_labels:
        # Count occurrences of each label
        unique, counts = np.unique(row, return_counts=True)
        best_idx = np.argmax(counts)
        winner = unique[best_idx]
        n_agree = counts[best_idx]
        conf = n_agree / n_neighbors
        predicted_clusters.append(winner)
        confidences.append(conf)
        n_agreeing_list.append(int(n_agree))

    query_cell_ids = combined.obs_names[is_query.values].tolist()

    return pd.DataFrame(
        {
            "cell_id": query_cell_ids,
            "predicted_cluster": predicted_clusters,
            "confidence": confidences,
            "n_neighbors_agreeing": n_agreeing_list,
        }
    )


def evaluate_mapping(
    predictions: pd.DataFrame,
    ground_truth_col: str,
    metadata: pd.DataFrame,
) -> dict:
    """Evaluate label-transfer predictions against ground-truth annotations.

    Parameters
    ----------
    predictions:
        DataFrame from :func:`run_harmony_mapping` with at least ``cell_id``
        and ``predicted_cluster`` columns.
    ground_truth_col:
        Column in *metadata* containing ground-truth cluster labels.
    metadata:
        DataFrame indexed by cell ID (or with ``cell_id`` column) holding
        ground-truth annotations.

    Returns
    -------
    dict with keys:
        - ``overall_accuracy``: fraction of cells correctly predicted at the
          subclass level (first word of label string).
        - ``non_neuronal_rate``: fraction of predicted labels that are
          non-neuronal.
        - ``n_cells``: total number of cells evaluated.
        - ``confusion_matrix``: DataFrame (rows = ground truth subclass,
          cols = predicted subclass) of cell counts.
    """
    # ------------------------------------------------------------------
    # Merge predictions with metadata
    # ------------------------------------------------------------------
    if "cell_id" in metadata.columns:
        meta = metadata.set_index("cell_id")
    else:
        meta = metadata.copy()

    merged = predictions.set_index("cell_id").join(meta[[ground_truth_col]], how="inner")

    gt_labels = merged[ground_truth_col].astype(str)
    pred_labels = merged["predicted_cluster"].astype(str)

    # ------------------------------------------------------------------
    # Extract subclass (first word) from each label
    # ------------------------------------------------------------------
    def _subclass(label: str) -> str:
        return label.split()[0] if label else label

    gt_subclass = gt_labels.map(_subclass)
    pred_subclass = pred_labels.map(_subclass)

    n_cells = len(merged)

    # ------------------------------------------------------------------
    # Overall accuracy at subclass level
    # ------------------------------------------------------------------
    overall_accuracy = (gt_subclass == pred_subclass).sum() / n_cells if n_cells > 0 else 0.0

    # ------------------------------------------------------------------
    # Non-neuronal rate in predictions
    # ------------------------------------------------------------------
    def _is_non_neuronal(label: str) -> bool:
        sub = _subclass(label)
        if sub in _NON_NEURONAL_TYPES:
            return True
        if label.startswith(_NON_NEURONAL_PREFIX):
            return True
        return False

    non_neuronal_mask = pred_labels.map(_is_non_neuronal)
    non_neuronal_rate = non_neuronal_mask.sum() / n_cells if n_cells > 0 else 0.0

    # ------------------------------------------------------------------
    # Confusion matrix
    # ------------------------------------------------------------------
    all_labels = sorted(set(gt_subclass) | set(pred_subclass))
    confusion = pd.crosstab(
        gt_subclass.rename("ground_truth"),
        pred_subclass.rename("predicted"),
    )

    return {
        "overall_accuracy": float(overall_accuracy),
        "non_neuronal_rate": float(non_neuronal_rate),
        "n_cells": int(n_cells),
        "confusion_matrix": confusion,
    }
