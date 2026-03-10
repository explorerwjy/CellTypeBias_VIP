#!/usr/bin/env python3
"""
Step 2H-scANVI: Train scANVI on Human Reference Atlas + Map Patch-seq Query

This script uses scANVI (semi-supervised) instead of scVI (unsupervised).
Key difference: scANVI learns from expression + supercluster labels, enforcing
interneuron vs non-interneuron separation in the latent space.

This script:
1. Loads preprocessed reference (Siletti et al. 2023 interneurons + decoys) and query (Patch-seq)
2. Trains scANVI on reference with supercluster_term labels (semi-supervised)
3. Maps Patch-seq query via scArches surgery with new batch category + unlabeled labels
4. Performs hierarchical cell type assignment (supercluster → cluster → subcluster)
5. Saves trained model, latent embeddings, and mapping results

Config: atlas_matching/configs/human_interneuron_config.yaml
Linear: TRA-30 Phase 2 (scANVI variant)
"""

import sys
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch

warnings.filterwarnings("ignore")

# Fix: PyTorch 2.6+ defaults weights_only=True; scvi-tools 1.2.0 predates this
_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

# Add module to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "atlas_matching"))

from patch_seq_transcriptome_mapping import map_query_scarches


# =============================================================================
# Hierarchical assignment adapted for Siletti atlas hierarchy
# =============================================================================

def hierarchical_assignment_human(
    adata_ref: ad.AnnData,
    adata_query: ad.AnnData,
    latent_key: str = "X_latent",
    supercluster_key: str = "supercluster_term",
    cluster_key: str = "cluster_id",
    subcluster_key: str = "subcluster_id",
    k_candidates: int = 30,
    n_neighbors_global: int = 100,
    min_neighbors_candidates: int = 10,
    conf_threshold_high: float = 0.8,
    conf_threshold_medium: float = 0.5,
    conf_threshold_low: float = 0.3,
    distance_weighted: bool = True,
) -> pd.DataFrame:
    """
    Hierarchical assignment for the Siletti human brain atlas.

    Uses the 3-level hierarchy: Supercluster → Cluster → Subcluster.
    Single global KNN on all reference cells, then derive assignments
    at all three levels from the same neighbor set.

    Parameters
    ----------
    adata_ref : ad.AnnData
        Reference atlas with latent embeddings and labels
    adata_query : ad.AnnData
        Query Patch-seq data with latent embeddings
    latent_key : str
        Key in obsm containing latent embeddings
    supercluster_key : str
        Column name for supercluster labels (e.g., "CGE interneuron")
    cluster_key : str
        Column name for cluster IDs (e.g., 276)
    subcluster_key : str
        Column name for subcluster IDs (e.g., 276.1)
    k_candidates : int
        Number of candidate clusters in centroid stage
    n_neighbors_global : int
        Total neighbors from global KNN
    min_neighbors_candidates : int
        Minimum valid neighbors in candidate clusters before fallback
    conf_threshold_high : float
        Supercluster confidence threshold for high-confidence tier
    conf_threshold_medium : float
        Medium confidence threshold
    conf_threshold_low : float
        Low confidence threshold
    distance_weighted : bool
        Weight votes by inverse distance

    Returns
    -------
    pd.DataFrame
        Mapping results with supercluster, cluster, subcluster assignments
    """
    from sklearn.neighbors import NearestNeighbors
    from scipy.spatial.distance import cdist

    print("=" * 60)
    print("HIERARCHICAL ATLAS MAPPING (HUMAN INTERNEURON)")
    print("=" * 60)

    ref_latent = adata_ref.obsm[latent_key]
    query_latent = adata_query.obsm[latent_key]
    n_query = len(query_latent)

    ref_superclusters = adata_ref.obs[supercluster_key].values.astype(str)
    ref_clusters = adata_ref.obs[cluster_key].values.astype(str)
    ref_subclusters = adata_ref.obs[subcluster_key].values.astype(str)

    # ------------------------------------------------------------------
    # Stage 1: Cluster centroids and candidate selection
    # ------------------------------------------------------------------
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
    print(f"  Computed {n_clusters} cluster centroids")

    # Find candidates per query cell
    dists_to_centroids = cdist(query_latent, centroids, metric="euclidean")
    candidate_indices = np.argsort(dists_to_centroids, axis=1)[:, :k_candidates]
    candidate_clusters_per_cell = np.array(
        [[idx_to_cluster[idx] for idx in row] for row in candidate_indices]
    )

    # ------------------------------------------------------------------
    # Stage 2: Global KNN
    # ------------------------------------------------------------------
    k_actual = min(n_neighbors_global, len(ref_latent))
    print(f"\n[Stage 2] Fitting global KNN (k={k_actual}) on {len(ref_latent):,} reference cells...")
    nn = NearestNeighbors(n_neighbors=k_actual, metric="euclidean", n_jobs=-1)
    nn.fit(ref_latent)
    distances, indices = nn.kneighbors(query_latent)
    print(f"  Done. Neighbor matrix: {indices.shape}")

    # Precompute neighbor labels
    neighbor_superclusters = ref_superclusters[indices]
    neighbor_clusters = ref_clusters[indices]
    neighbor_subclusters = ref_subclusters[indices]

    # ------------------------------------------------------------------
    # Stage 3: Voting
    # ------------------------------------------------------------------
    print(f"\n[Stage 3] Voting from candidate-filtered neighbors...")

    assigned_supercluster = np.empty(n_query, dtype=object)
    conf_supercluster = np.zeros(n_query)
    assigned_cluster = np.empty(n_query, dtype=object)
    conf_cluster = np.zeros(n_query)
    assigned_subcluster = np.empty(n_query, dtype=object)
    conf_subcluster = np.zeros(n_query)
    n_neighbors_used = np.zeros(n_query, dtype=int)
    used_fallback = np.zeros(n_query, dtype=bool)

    def _weighted_vote(labels, weights):
        """Weighted majority vote."""
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
        assigned_subcluster[i], conf_subcluster[i] = _weighted_vote(
            neighbor_subclusters[i, use_mask], weights
        )

    n_fallback = used_fallback.sum()
    print(f"  Fallback: {n_fallback}/{n_query} cells ({100 * n_fallback / n_query:.1f}%)")

    # ------------------------------------------------------------------
    # Stage 4: Tiered confidence
    # ------------------------------------------------------------------
    print("\n[Stage 4] Tiered confidence scoring...")
    mapping_tier = np.empty(n_query, dtype=object)
    mapping_status = np.empty(n_query, dtype=object)

    for i in range(n_query):
        if conf_supercluster[i] >= conf_threshold_high and conf_cluster[i] >= conf_threshold_medium:
            mapping_tier[i] = "high_confidence"
            mapping_status[i] = "ok_cluster"
        elif conf_supercluster[i] >= conf_threshold_medium and conf_cluster[i] >= conf_threshold_low:
            mapping_tier[i] = "medium_confidence"
            mapping_status[i] = "ok_cluster"
        elif conf_supercluster[i] >= conf_threshold_low:
            mapping_tier[i] = "low_confidence"
            mapping_status[i] = "ok_supercluster_only"
        else:
            mapping_tier[i] = "rejected"
            mapping_status[i] = "rejected"

    # ------------------------------------------------------------------
    # Interneuron purity check
    # ------------------------------------------------------------------
    interneuron_superclusters = {"CGE interneuron", "MGE interneuron", "LAMP5-LHX6 and Chandelier"}
    is_interneuron = np.array([s in interneuron_superclusters for s in assigned_supercluster])
    n_interneuron = is_interneuron.sum()

    # ------------------------------------------------------------------
    # Build output
    # ------------------------------------------------------------------
    results = pd.DataFrame(
        {
            "assigned_supercluster": assigned_supercluster,
            "conf_supercluster": conf_supercluster,
            "assigned_cluster": assigned_cluster,
            "conf_cluster": conf_cluster,
            "assigned_subcluster": assigned_subcluster,
            "conf_subcluster": conf_subcluster,
            "n_neighbors_used": n_neighbors_used,
            "used_fallback": used_fallback,
            "mapping_tier": mapping_tier,
            "mapping_status": mapping_status,
            "is_interneuron": is_interneuron,
        },
        index=adata_query.obs.index,
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("MAPPING RESULTS SUMMARY")
    print("=" * 60)

    tier_counts = pd.Series(mapping_tier).value_counts()
    for tier, count in tier_counts.items():
        print(f"  {tier}: {count} cells ({100 * count / n_query:.1f}%)")

    print(f"\n  Mean supercluster confidence: {conf_supercluster.mean():.3f}")
    print(f"  Mean cluster confidence: {conf_cluster.mean():.3f}")
    print(f"  Mean subcluster confidence: {conf_subcluster.mean():.3f}")
    print(f"  Fallback rate: {100 * n_fallback / n_query:.1f}%")

    print(f"\n  Interneuron purity: {n_interneuron}/{n_query} ({100 * n_interneuron / n_query:.1f}%)")
    print(f"  Non-interneuron assignments: {n_query - n_interneuron}")

    print(f"\n  Supercluster assignment distribution:")
    for lbl, cnt in pd.Series(assigned_supercluster).value_counts().items():
        print(f"    {lbl}: {cnt} ({100 * cnt / n_query:.1f}%)")

    return results


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Step 2H-scANVI: Train scANVI + Map Human Patch-seq",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Directories
    parser.add_argument(
        "--preprocessed-dir",
        type=Path,
        default=PROJECT_ROOT / "atlas_matching" / "results" / "human_interneuron" / "preprocessed",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=PROJECT_ROOT / "atlas_matching" / "results" / "human_interneuron" / "models",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "atlas_matching" / "results" / "human_interneuron" / "scanvi_mapped",
    )

    # scANVI training
    parser.add_argument("--n-latent", type=int, default=50)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--n-hidden", type=int, default=256)
    parser.add_argument("--max-epochs", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--batch-key", type=str, default="sample_id",
                        help="Batch key for within-reference correction")
    parser.add_argument("--labels-key", type=str, default="supercluster_term",
                        help="Column with cell type labels for semi-supervised training")
    parser.add_argument("--no-gpu", action="store_true")
    parser.add_argument("--force-retrain", action="store_true")

    # scArches query mapping
    parser.add_argument("--max-epochs-surgery", type=int, default=200)

    # Hierarchical assignment
    parser.add_argument("--k-candidates", type=int, default=30)
    parser.add_argument("--n-neighbors-global", type=int, default=100)
    parser.add_argument("--min-neighbors-candidates", type=int, default=10)
    parser.add_argument("--conf-threshold-high", type=float, default=0.8)
    parser.add_argument("--conf-threshold-medium", type=float, default=0.5)
    parser.add_argument("--conf-threshold-low", type=float, default=0.3)

    args = parser.parse_args()
    use_gpu = not args.no_gpu

    # Setup directories
    args.model_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    model_name = f"scanvi_human_latent{args.n_latent}_hidden{args.n_hidden}"
    model_path = args.model_dir / model_name

    print("=" * 80)
    print("STEP 2H-scANVI: TRAIN SCANVI + MAP HUMAN PATCH-SEQ")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  scANVI: latent={args.n_latent}, layers={args.n_layers}, hidden={args.n_hidden}")
    print(f"  Labels key: {args.labels_key} (semi-supervised)")
    print(f"  Max epochs (train): {args.max_epochs}")
    print(f"  Max epochs (surgery): {args.max_epochs_surgery}")
    print(f"  Batch key: {args.batch_key}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  GPU: {use_gpu}")
    print(f"\nHierarchical assignment:")
    print(f"  K candidates: {args.k_candidates}")
    print(f"  N neighbors: {args.n_neighbors_global}")
    print(f"  Confidence thresholds: {args.conf_threshold_high}/{args.conf_threshold_medium}/{args.conf_threshold_low}")
    print(f"\nPaths:")
    print(f"  Preprocessed: {args.preprocessed_dir}")
    print(f"  Model: {model_path}")
    print(f"  Output: {args.output_dir}")
    print("=" * 80)

    # ================================================================
    # Load preprocessed data
    # ================================================================
    print("\n[LOAD] Loading preprocessed data...")
    ref_path = args.preprocessed_dir / "human_reference_preprocessed.h5ad"
    query_path = args.preprocessed_dir / "human_patchseq_query_preprocessed.h5ad"

    if not ref_path.exists():
        raise FileNotFoundError(f"Reference not found: {ref_path}\nRun 01h_preprocess_human_interneuron.py first")
    if not query_path.exists():
        raise FileNotFoundError(f"Query not found: {query_path}\nRun 01h_preprocess_human_interneuron.py first")

    adata_ref = sc.read_h5ad(ref_path)
    adata_query = sc.read_h5ad(query_path)

    print(f"  Reference: {adata_ref.shape}")
    print(f"  Query: {adata_query.shape}")

    # Ensure counts layers exist
    if "counts" not in adata_ref.layers:
        print("  Adding counts layer to reference (copying .X)")
        adata_ref.layers["counts"] = adata_ref.X.copy()
    if "counts" not in adata_query.layers:
        print("  Adding counts layer to query (copying .X)")
        adata_query.layers["counts"] = adata_query.X.copy()

    # Verify gene match
    if not all(adata_query.var_names == adata_ref.var_names):
        raise ValueError("Gene names don't match between reference and query!")
    print(f"  Verified: {adata_ref.n_vars} matching genes")

    # Print reference obs columns
    print(f"  Reference obs: {list(adata_ref.obs.columns)}")
    print(f"  Query obs: {list(adata_query.obs.columns)}")

    # Check batch key exists
    if args.batch_key not in adata_ref.obs.columns:
        available = list(adata_ref.obs.columns)
        raise ValueError(
            f"Batch key '{args.batch_key}' not found in reference obs. "
            f"Available: {available}"
        )
    n_batches = adata_ref.obs[args.batch_key].nunique()
    print(f"  Batch key '{args.batch_key}': {n_batches} unique values")

    # ================================================================
    # Phase 2a: Train scANVI on reference (semi-supervised)
    # ================================================================
    model_exists = model_path.exists() and (model_path / "model.pt").exists()

    if model_exists and not args.force_retrain:
        print(f"\n[TRAIN] Existing model found at {model_path}, skipping training")
        print(f"  Use --force-retrain to retrain")

        # Load pre-computed reference latent from saved model adata
        saved_adata_path = model_path / "adata.h5ad"
        if saved_adata_path.exists():
            saved_adata = ad.read_h5ad(saved_adata_path)
            if "X_latent" in saved_adata.obsm:
                adata_ref.obsm["X_latent"] = saved_adata.obsm["X_latent"]
                print(f"  Loaded reference latent: {adata_ref.obsm['X_latent'].shape}")
            del saved_adata
    else:
        print(f"\n[TRAIN] Training scANVI on reference atlas (semi-supervised)...")

        import scvi
        from scvi.model import SCANVI

        # Make working copy for training
        adata_train = adata_ref.copy()

        # Setup batch key — ensure string type for scvi-tools
        # CRITICAL: Pass ORIGINAL column names to setup_anndata, NOT _scvi_* columns.
        # scvi-tools creates _scvi_batch/_scvi_labels internally from the original columns.
        # Naming our columns _scvi_* causes a collision where codes overwrite the original data.
        batch_col = "_batch_str"
        adata_train.obs[batch_col] = adata_train.obs[args.batch_key].astype(str).values
        print(f"  Batch key: {args.batch_key} ({adata_train.obs[batch_col].nunique()} categories)")

        # Setup labels for semi-supervised training
        labels_col = "_labels_str"
        adata_train.obs[labels_col] = adata_train.obs[args.labels_key].astype(str).values
        n_labels = adata_train.obs[labels_col].nunique()
        print(f"  Labels key: {args.labels_key} ({n_labels} types)")
        for lbl, cnt in adata_train.obs[labels_col].value_counts().items():
            print(f"    {lbl}: {cnt:,}")

        # Register AnnData with scANVI
        SCANVI.setup_anndata(
            adata_train,
            batch_key=batch_col,
            labels_key=labels_col,
            unlabeled_category="Unknown",
            layer="counts",
        )

        # Create scANVI model
        model = SCANVI(
            adata_train,
            n_latent=args.n_latent,
            n_layers=args.n_layers,
            n_hidden=args.n_hidden,
            gene_likelihood="nb",
        )

        print(f"  Model: scANVI, latent={args.n_latent}, layers={args.n_layers}, hidden={args.n_hidden}")
        print(f"  Gene likelihood: Negative Binomial")

        # Train
        train_kwargs = {
            "max_epochs": args.max_epochs,
            "train_size": 0.9,
            "batch_size": args.batch_size,
            "early_stopping": True,
            "early_stopping_patience": 20,
            "check_val_every_n_epoch": 5,
            "enable_progress_bar": True,
        }
        if use_gpu:
            train_kwargs["accelerator"] = "gpu"
            train_kwargs["devices"] = "auto"
        else:
            train_kwargs["accelerator"] = "cpu"

        model.train(**train_kwargs)

        # Get latent representation
        adata_train.obsm["X_latent"] = model.get_latent_representation(adata_train)
        print(f"  Latent embedding: {adata_train.obsm['X_latent'].shape}")

        # Copy latent back to adata_ref
        adata_ref.obsm["X_latent"] = adata_train.obsm["X_latent"]

        # Save model
        print(f"\n[SAVE] Saving trained model to {model_path}...")
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model.save(model_path, overwrite=True, save_anndata=True)

        # Fix saved adata (scANVI encodes categoricals as int8)
        print("  Fixing saved adata categoricals...")
        saved_adata_path = model_path / "adata.h5ad"
        saved_adata = ad.read_h5ad(saved_adata_path)

        # Get original categories from the registered AnnData
        original_batch_cats = sorted(adata_train.obs[batch_col].unique())
        original_label_cats = sorted(adata_train.obs[labels_col].unique())

        # scvi-tools stores codes in _scvi_batch and _scvi_labels
        if "_scvi_batch" in saved_adata.obs.columns and saved_adata.obs["_scvi_batch"].dtype == "int8":
            saved_adata.obs["_scvi_batch"] = pd.Categorical.from_codes(
                saved_adata.obs["_scvi_batch"].values, categories=original_batch_cats
            )
            print(f"  Fixed _scvi_batch: int8 -> categorical ({len(original_batch_cats)} categories)")

        if "_scvi_labels" in saved_adata.obs.columns and saved_adata.obs["_scvi_labels"].dtype == "int8":
            # scANVI adds Unknown as last category
            label_cats_with_unknown = original_label_cats + ["Unknown"]
            saved_adata.obs["_scvi_labels"] = pd.Categorical.from_codes(
                saved_adata.obs["_scvi_labels"].values, categories=label_cats_with_unknown
            )
            print(f"  Fixed _scvi_labels: int8 -> categorical ({len(label_cats_with_unknown)} categories)")

        saved_adata.write(saved_adata_path)
        print("  Model saved successfully")

    # ================================================================
    # Phase 2b: Map query via scANVI scArches surgery (inline)
    # ================================================================
    # NOTE: We bypass map_query_scarches() and do surgery directly here because
    # the model was trained with custom column names (_batch_str, _labels_str)
    # to work around a scvi-tools 1.2.0 naming collision bug. load_query_data()
    # expects those exact columns in the query adata.
    print(f"\n[MAP] Mapping Patch-seq query via scANVI scArches surgery...")
    print(f"  Loading model from: {model_path}")

    from scvi.model import SCANVI

    adata_query_for_mapping = adata_query.copy()

    # Set up the SAME column names used during training
    # _batch_str: new Patch-seq batch category (distinct from all reference batches)
    adata_query_for_mapping.obs["_batch_str"] = "PatchSeq_Human_GABA"
    # _labels_str: "Unknown" so scANVI treats query cells as unlabeled
    adata_query_for_mapping.obs["_labels_str"] = "Unknown"

    print(f"  Query batch (_batch_str): PatchSeq_Human_GABA (new category)")
    print(f"  Query labels (_labels_str): Unknown (unlabeled for prediction)")

    # Ensure counts layer exists
    if "counts" not in adata_query_for_mapping.layers:
        raise ValueError("Query adata must have layers['counts'] with raw counts")

    # scArches surgery: load_query_data freezes reference weights,
    # adds batch-specific parameters for the new query batch
    print(f"  Performing scArches surgery (max_epochs={args.max_epochs_surgery})...")
    model_query = SCANVI.load_query_data(
        adata_query_for_mapping,
        str(model_path),
        freeze_dropout=True,
    )

    # Train only the batch-specific parameters
    train_kwargs = {
        "max_epochs": args.max_epochs_surgery,
        "batch_size": 128,
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

    # Get latent representations
    print(f"  Computing latent embeddings...")
    adata_query_for_mapping.obsm["X_latent"] = model_query.get_latent_representation(
        adata_query_for_mapping
    )

    # Use pre-computed reference latent from training phase.
    # Reference weights are frozen during surgery, so the latent space is unchanged.
    adata_ref_mapped = adata_ref.copy()
    if "X_latent" not in adata_ref_mapped.obsm:
        raise RuntimeError(
            "Reference latent not found. Either run training first "
            "(remove --force-retrain skip) or ensure the saved model has X_latent."
        )

    adata_query_mapped = adata_query_for_mapping
    print(f"  Reference latent: {adata_ref_mapped.obsm['X_latent'].shape}")
    print(f"  Query latent: {adata_query_mapped.obsm['X_latent'].shape}")

    # Get scANVI predicted labels for query cells
    try:
        predicted_labels = model_query.predict(adata_query_mapped)
        adata_query_mapped.obs["scanvi_predicted"] = predicted_labels
        print(f"\n  scANVI predicted labels:")
        for lbl, cnt in pd.Series(predicted_labels).value_counts().items():
            print(f"    {lbl}: {cnt}")
    except Exception as e:
        print(f"  Warning: Could not get scANVI predictions: {e}")
        print(f"  Falling back to KNN-based assignment only")

    # ================================================================
    # Phase 2c: Hierarchical assignment
    # ================================================================
    print(f"\n[ASSIGN] Hierarchical cell type assignment...")

    mapping_results = hierarchical_assignment_human(
        adata_ref_mapped,
        adata_query_mapped,
        latent_key="X_latent",
        supercluster_key="supercluster_term",
        cluster_key="cluster_id",
        subcluster_key="subcluster_id",
        k_candidates=args.k_candidates,
        n_neighbors_global=args.n_neighbors_global,
        min_neighbors_candidates=args.min_neighbors_candidates,
        conf_threshold_high=args.conf_threshold_high,
        conf_threshold_medium=args.conf_threshold_medium,
        conf_threshold_low=args.conf_threshold_low,
    )

    # Merge mapping results into query obs
    adata_query_mapped.obs = pd.concat([adata_query_mapped.obs, mapping_results], axis=1)

    # ================================================================
    # Save results
    # ================================================================
    print(f"\n[SAVE] Saving results to {args.output_dir}...")

    # Latent embeddings (small, save first)
    np.save(args.output_dir / "reference_latent.npy", adata_ref_mapped.obsm["X_latent"])
    np.save(args.output_dir / "query_latent.npy", adata_query_mapped.obsm["X_latent"])
    print(f"  Reference latent: {args.output_dir / 'reference_latent.npy'} "
          f"({adata_ref_mapped.obsm['X_latent'].shape})")
    print(f"  Query latent: {args.output_dir / 'query_latent.npy'} "
          f"({adata_query_mapped.obsm['X_latent'].shape})")

    # Mapping results CSV (most important output)
    mapping_csv = args.output_dir / "mapping_results.csv"
    mapping_results.to_csv(mapping_csv)
    print(f"  Mapping CSV: {mapping_csv}")

    # Query mapped h5ad (778 cells — small)
    query_out = args.output_dir / "query_mapped.h5ad"
    adata_query_mapped.write(query_out)
    print(f"  Query h5ad: {query_out}")

    # Reference obs with latent index (save as CSV to avoid huge h5ad write)
    ref_obs_out = args.output_dir / "reference_obs.csv"
    adata_ref_mapped.obs.to_csv(ref_obs_out)
    print(f"  Reference obs: {ref_obs_out}")

    # Quick concordance check with Patch-seq labels
    if "subclass_label" in adata_query_mapped.obs.columns:
        print("\n" + "=" * 60)
        print("CONCORDANCE CHECK: Patch-seq labels vs Atlas assignment")
        print("=" * 60)

        obs = adata_query_mapped.obs
        cross = pd.crosstab(
            obs["subclass_label"],
            obs["assigned_supercluster"],
            margins=True,
        )
        print(cross.to_string())

        # Expected mapping check
        expected = {
            "PVALB": "MGE interneuron",
            "SST": "MGE interneuron",
            "VIP": "CGE interneuron",
        }
        total_correct = 0
        total_checked = 0
        for ps_label, expected_sc in expected.items():
            mask = obs["subclass_label"] == ps_label
            n_total = mask.sum()
            if n_total == 0:
                continue
            n_correct = (obs.loc[mask, "assigned_supercluster"] == expected_sc).sum()
            pct = 100 * n_correct / n_total
            print(f"  {ps_label} → {expected_sc}: {n_correct}/{n_total} ({pct:.1f}%)")
            total_correct += n_correct
            total_checked += n_total
        if total_checked > 0:
            print(f"  Overall concordance (PVALB/SST/VIP): {total_correct}/{total_checked} "
                  f"({100 * total_correct / total_checked:.1f}%)")

    print(f"\n{'=' * 80}")
    print("PHASE 2 (scANVI) COMPLETE")
    print(f"{'=' * 80}")
    print(f"  Model: {model_path}")
    print(f"  Model type: scANVI (semi-supervised)")
    print(f"  Results: {args.output_dir}")
    print(f"\nNext step: evaluation (Phase 3/4)")


if __name__ == "__main__":
    main()
