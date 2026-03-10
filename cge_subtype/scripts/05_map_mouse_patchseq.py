#!/usr/bin/env python3
"""
Step 5: Map mouse Patch-seq cells via scArches into the cross-species latent space.

Uses scArches surgery to add mouse Patch-seq cells to the pre-trained
cross-species scVI model, then performs hierarchical assignment to
Siletti atlas clusters.

Output:
  - results/cross_species/mapped/patchseq_mapping_results.csv
  - results/cross_species/mapped/patchseq_query_mapped.h5ad
"""

import sys
import gc
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
import yaml
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import cdist

warnings.filterwarnings("ignore")

_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_MATCHING_DIR = SCRIPT_DIR.parent.parent
PROJECT_ROOT = ATLAS_MATCHING_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(ATLAS_MATCHING_DIR))


def load_config():
    config_path = ATLAS_MATCHING_DIR / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def hierarchical_assignment_cross_species(
    adata_ref,
    adata_query,
    latent_key="X_latent",
    supercluster_key="supercluster_term",
    cluster_key="cluster_id",
    k_candidates=30,
    n_neighbors_global=100,
    min_neighbors_candidates=10,
    conf_threshold_high=0.8,
    conf_threshold_medium=0.5,
    conf_threshold_low=0.3,
    distance_weighted=True,
):
    """
    Hierarchical assignment for cross-species mapping.

    Same algorithm as hierarchical_assignment_human (02h_train_and_map_human.py)
    but only assigns supercluster + cluster (no subcluster, since mouse cells
    don't have meaningful subcluster assignments).

    Only uses HUMAN reference cells for assignment (mouse reference cells are
    excluded since we want to map TO human clusters).
    """
    print("=" * 60)
    print("HIERARCHICAL ATLAS MAPPING (CROSS-SPECIES)")
    print("=" * 60)

    # Filter reference to human cells only for assignment
    human_mask = adata_ref.obs["species"] == "human"
    ref_latent = adata_ref.obsm[latent_key][human_mask]
    ref_obs = adata_ref.obs[human_mask]
    query_latent = adata_query.obsm[latent_key]
    n_query = len(query_latent)

    ref_superclusters = ref_obs[supercluster_key].values.astype(str)
    ref_clusters = ref_obs[cluster_key].values.astype(str)

    print(f"  Reference (human only): {len(ref_latent):,} cells")
    print(f"  Query (mouse Patch-seq): {n_query} cells")

    # Stage 1: Cluster centroids
    print(f"\n[Stage 1] Computing cluster centroids (k={k_candidates})...")
    clusters_unique = sorted(pd.Series(ref_clusters).dropna().unique())
    n_clusters = len(clusters_unique)
    centroids = np.zeros((n_clusters, ref_latent.shape[1]))
    cluster_to_idx = {}
    for idx, clust in enumerate(clusters_unique):
        mask = ref_clusters == clust
        centroids[idx] = ref_latent[mask].mean(axis=0)
        cluster_to_idx[clust] = idx
    idx_to_cluster = {v: k for k, v in cluster_to_idx.items()}
    print(f"  {n_clusters} cluster centroids computed")

    dists_to_centroids = cdist(query_latent, centroids, metric="euclidean")
    candidate_indices = np.argsort(dists_to_centroids, axis=1)[:, :k_candidates]
    candidate_clusters_per_cell = np.array(
        [[idx_to_cluster[idx] for idx in row] for row in candidate_indices]
    )

    # Stage 2: Global KNN
    k_actual = min(n_neighbors_global, len(ref_latent))
    print(f"\n[Stage 2] Fitting global KNN (k={k_actual})...")
    nn = NearestNeighbors(n_neighbors=k_actual, metric="euclidean", n_jobs=10)
    nn.fit(ref_latent)
    distances, indices = nn.kneighbors(query_latent)

    neighbor_superclusters = ref_superclusters[indices]
    neighbor_clusters = ref_clusters[indices]

    # Stage 3: Voting
    print(f"\n[Stage 3] Voting from candidate-filtered neighbors...")

    assigned_supercluster = np.empty(n_query, dtype=object)
    conf_supercluster = np.zeros(n_query)
    assigned_cluster = np.empty(n_query, dtype=object)
    conf_cluster = np.zeros(n_query)
    n_neighbors_used = np.zeros(n_query, dtype=int)
    used_fallback = np.zeros(n_query, dtype=bool)

    def _weighted_vote(labels, weights):
        labels = np.asarray(labels)
        weights = np.asarray(weights, dtype=float)
        unique = np.unique(labels)
        scores = {lbl: weights[labels == lbl].sum() for lbl in unique}
        winner = max(scores, key=scores.get)
        total = sum(scores.values())
        return winner, scores[winner] / total if total > 0 else 0.0

    for i in range(n_query):
        cand_set = set(candidate_clusters_per_cell[i])
        in_candidates = np.array([neighbor_clusters[i, j] in cand_set for j in range(k_actual)])
        n_in = in_candidates.sum()

        if n_in >= min_neighbors_candidates:
            use_mask = in_candidates
            used_fallback[i] = False
        else:
            use_mask = np.ones(k_actual, dtype=bool)
            used_fallback[i] = True

        n_neighbors_used[i] = use_mask.sum()
        dists_sel = distances[i, use_mask]

        if distance_weighted:
            weights = 1.0 / (dists_sel + 1e-6)
            weights = weights / weights.sum()
        else:
            weights = np.ones(use_mask.sum()) / use_mask.sum()

        assigned_supercluster[i], conf_supercluster[i] = _weighted_vote(
            neighbor_superclusters[i, use_mask], weights
        )
        assigned_cluster[i], conf_cluster[i] = _weighted_vote(
            neighbor_clusters[i, use_mask], weights
        )

    # Stage 4: Tiered confidence
    print(f"\n[Stage 4] Tiered confidence scoring...")
    mapping_tier = np.empty(n_query, dtype=object)

    for i in range(n_query):
        if conf_supercluster[i] >= conf_threshold_high and conf_cluster[i] >= conf_threshold_medium:
            mapping_tier[i] = "high_confidence"
        elif conf_supercluster[i] >= conf_threshold_medium and conf_cluster[i] >= conf_threshold_low:
            mapping_tier[i] = "medium_confidence"
        elif conf_supercluster[i] >= conf_threshold_low:
            mapping_tier[i] = "low_confidence"
        else:
            mapping_tier[i] = "rejected"

    # Build results
    results = pd.DataFrame(
        {
            "assigned_human_supercluster": assigned_supercluster,
            "conf_supercluster": conf_supercluster,
            "assigned_human_cluster": assigned_cluster,
            "conf_cluster": conf_cluster,
            "n_neighbors_used": n_neighbors_used,
            "used_fallback": used_fallback,
            "mapping_tier": mapping_tier,
        },
        index=adata_query.obs.index,
    )

    # Add query metadata
    for col in ["dataset", "is_cckbc", "species"]:
        if col in adata_query.obs.columns:
            results[col] = adata_query.obs[col].values

    # Add mouse subclass/supertype if available
    for col in ["subclass", "supertype", "corresponding_AIT2.3.1_alias",
                 "RNA type", "RNA family"]:
        if col in adata_query.obs.columns:
            results[f"mouse_{col}"] = adata_query.obs[col].values

    # Summary
    print(f"\n{'=' * 60}")
    print("MAPPING RESULTS SUMMARY")
    print(f"{'=' * 60}")

    tier_counts = pd.Series(mapping_tier).value_counts()
    for tier, count in tier_counts.items():
        print(f"  {tier}: {count} cells ({100*count/n_query:.1f}%)")

    print(f"\n  Mean supercluster confidence: {conf_supercluster.mean():.3f}")
    print(f"  Mean cluster confidence: {conf_cluster.mean():.3f}")

    print(f"\n  Supercluster distribution:")
    for lbl, cnt in pd.Series(assigned_supercluster).value_counts().items():
        print(f"    {lbl}: {cnt} ({100*cnt/n_query:.1f}%)")

    # CCKBC-specific summary
    if "is_cckbc" in adata_query.obs.columns:
        cckbc_mask = adata_query.obs["is_cckbc"].values.astype(bool)
        if cckbc_mask.any():
            print(f"\n  CCKBC cells ({cckbc_mask.sum()}):")
            cckbc_results = results[cckbc_mask]
            for lbl, cnt in cckbc_results["assigned_human_supercluster"].value_counts().items():
                print(f"    -> {lbl}: {cnt}")
            print(f"  CCKBC cluster assignments:")
            for clust, cnt in cckbc_results["assigned_human_cluster"].value_counts().head(10).items():
                print(f"    -> cluster {clust}: {cnt}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Map mouse Patch-seq via scArches")
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--model-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-gpu", action="store_true")
    args = parser.parse_args()

    config = load_config()
    preproc_dir = args.preprocessed_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["preprocessed"])
    model_dir = args.model_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["models"])
    output_dir = args.output_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["mapped"])
    output_dir.mkdir(parents=True, exist_ok=True)
    use_gpu = not args.no_gpu

    scvi_config = config["scvi"]
    scarches_config = config["scarches"]
    mapping_config = config["mapping"]

    model_name = f"scvi_cross_species_lat{scvi_config['n_latent']}_hid{scvi_config['n_hidden']}"
    model_path = model_dir / model_name

    print("=" * 70)
    print("STEP 5: MAP MOUSE PATCH-SEQ VIA SCARCHES")
    print("=" * 70)

    # ================================================================
    # Load data
    # ================================================================
    print(f"\n[1/4] Loading data...")
    adata_ref = sc.read_h5ad(preproc_dir / "reference_with_latent.h5ad")
    adata_query = sc.read_h5ad(preproc_dir / "mouse_patchseq_query.h5ad")
    print(f"  Reference: {adata_ref.shape}")
    print(f"  Query: {adata_query.shape}")

    # Ensure query has same genes as reference
    ref_genes = set(adata_ref.var_names)
    query_genes = set(adata_query.var_names)
    common_genes = sorted(ref_genes & query_genes)
    print(f"  Common genes: {len(common_genes)} (ref={len(ref_genes)}, query={len(query_genes)})")

    # Save obs metadata (may be lost during gene manipulation with newer anndata)
    query_obs_saved = adata_query.obs.copy()

    adata_query = adata_query[:, common_genes].copy()
    if "counts" not in adata_query.layers:
        adata_query.layers["counts"] = adata_query.X.copy()

    # Pad query with zeros for genes in reference but not in query
    missing_genes = sorted(ref_genes - query_genes)
    if missing_genes:
        print(f"  Adding {len(missing_genes)} zero-padded genes to query...")
        zero_matrix = sparse.csr_matrix((adata_query.n_obs, len(missing_genes)), dtype=np.float32)
        zero_adata = ad.AnnData(
            X=zero_matrix,
            obs=pd.DataFrame(index=adata_query.obs_names),
            var=pd.DataFrame(index=missing_genes),
        )
        zero_adata.layers["counts"] = zero_adata.X.copy()
        adata_query = ad.concat([adata_query, zero_adata], axis=1, join="outer")
        adata_query.layers["counts"] = adata_query.X.copy()

    # Reorder to match reference gene order
    adata_query = adata_query[:, adata_ref.var_names].copy()
    adata_query.layers["counts"] = adata_query.X.copy()

    # Restore obs metadata
    adata_query.obs = query_obs_saved.loc[adata_query.obs_names]
    print(f"  Query aligned: {adata_query.shape}")
    print(f"  Obs columns: {len(adata_query.obs.columns)} (sample_species: {'sample_species' in adata_query.obs.columns})")

    # ================================================================
    # scArches surgery
    # ================================================================
    print(f"\n[2/4] Performing scArches surgery...")
    from scvi.model import SCVI

    # Prepare query batch
    adata_query_copy = adata_query.copy()

    # Clean categorical columns
    for col in adata_query_copy.obs.columns:
        if hasattr(adata_query_copy.obs[col], "cat"):
            adata_query_copy.obs[col] = adata_query_copy.obs[col].astype(str).astype("category")
        elif adata_query_copy.obs[col].dtype == "object":
            adata_query_copy.obs[col] = adata_query_copy.obs[col].astype(str)

    # Set query batch
    adata_query_copy.obs["_scvi_batch"] = (
        "PatchSeq_" + adata_query_copy.obs["sample_species"].astype(str)
    )
    adata_query_copy.obs["_scvi_batch"] = adata_query_copy.obs["_scvi_batch"].astype("category")
    print(f"  Query batches: {adata_query_copy.obs['_scvi_batch'].unique().tolist()}")

    # scArches: load query data into pre-trained model
    print(f"  Loading model from {model_path}...")
    model_query = SCVI.load_query_data(
        adata_query_copy,
        str(model_path),
        freeze_dropout=scarches_config["freeze_dropout"],
    )

    # Train batch parameters
    train_kwargs = {
        "max_epochs": scarches_config["max_epochs"],
        "batch_size": scarches_config["batch_size"],
        "early_stopping": True,
        "early_stopping_patience": 10,
        "check_val_every_n_epoch": 5,
        "enable_progress_bar": True,
    }
    if use_gpu:
        train_kwargs["accelerator"] = "gpu"
        train_kwargs["devices"] = "auto"
    else:
        train_kwargs["accelerator"] = "cpu"

    model_query.train(**train_kwargs)
    print(f"  scArches surgery complete")

    # Get query latent
    adata_query_copy.obsm["X_latent"] = model_query.get_latent_representation(adata_query_copy)
    print(f"  Query latent: {adata_query_copy.obsm['X_latent'].shape}")

    # Copy latent back to original query adata
    adata_query.obsm["X_latent"] = adata_query_copy.obsm["X_latent"]
    del adata_query_copy
    gc.collect()

    # ================================================================
    # Hierarchical assignment
    # ================================================================
    print(f"\n[3/4] Hierarchical cell type assignment...")
    mapping_results = hierarchical_assignment_cross_species(
        adata_ref,
        adata_query,
        latent_key="X_latent",
        supercluster_key="supercluster_term",
        cluster_key="cluster_id",
        k_candidates=mapping_config["k_candidates"],
        n_neighbors_global=mapping_config["n_neighbors_global"],
        min_neighbors_candidates=mapping_config["min_neighbors_candidates"],
        conf_threshold_high=mapping_config["conf_threshold_high"],
        conf_threshold_medium=mapping_config["conf_threshold_medium"],
        conf_threshold_low=mapping_config["conf_threshold_low"],
        distance_weighted=mapping_config["distance_weighted"],
    )

    # ================================================================
    # Save results
    # ================================================================
    print(f"\n[4/4] Saving results...")

    # Mapping CSV
    csv_path = output_dir / "patchseq_mapping_results.csv"
    mapping_results.to_csv(csv_path)
    print(f"  Mapping CSV: {csv_path}")

    # Query with latent — drop duplicate columns from mapping_results
    new_cols = [c for c in mapping_results.columns if c not in adata_query.obs.columns]
    adata_query.obs = pd.concat([adata_query.obs, mapping_results[new_cols]], axis=1)
    h5ad_path = output_dir / "patchseq_query_mapped.h5ad"
    adata_query.write(h5ad_path)
    print(f"  Query h5ad: {h5ad_path}")

    # Save latent
    np.save(output_dir / "patchseq_query_latent.npy", adata_query.obsm["X_latent"])

    print(f"\n{'=' * 70}")
    print("MAPPING COMPLETE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
