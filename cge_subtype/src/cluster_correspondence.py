"""Cross-species cluster correspondence: pseudobulk, Spearman RBH, MetaNeighbor AUROC."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.metrics import roc_auc_score


# ---------------------------------------------------------------------------
# Round 1: Pseudobulk
# ---------------------------------------------------------------------------

def compute_pseudobulk(expr: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Compute pseudobulk centroids by averaging expression within each cluster.

    Parameters
    ----------
    expr : pd.DataFrame
        Cell x gene expression matrix. Index = cell barcodes, columns = gene names.
    labels : pd.Series
        Cluster label for each cell. Index must match ``expr.index``.

    Returns
    -------
    pd.DataFrame
        Cluster x gene DataFrame of mean expression per cluster.
        Index = cluster labels, columns = gene names.
    """
    return expr.groupby(labels).mean()


# ---------------------------------------------------------------------------
# Round 2: Spearman correlation + RBH
# ---------------------------------------------------------------------------

def compute_spearman_corr_matrix(
    mouse_centroids: pd.DataFrame,
    human_centroids: pd.DataFrame,
) -> pd.DataFrame:
    """Compute Spearman correlation between all mouse and human cluster centroids.

    Only genes present in both datasets (shared genes) are used.

    Parameters
    ----------
    mouse_centroids : pd.DataFrame
        Mouse cluster x gene pseudobulk matrix.
    human_centroids : pd.DataFrame
        Human cluster x gene pseudobulk matrix.

    Returns
    -------
    pd.DataFrame
        (n_mouse_clusters) x (n_human_clusters) Spearman correlation matrix.
        Index = mouse cluster labels, columns = human cluster labels.
    """
    shared_genes = mouse_centroids.columns.intersection(human_centroids.columns)
    mouse_vals = mouse_centroids[shared_genes].values  # shape: (n_mouse, n_genes)
    human_vals = human_centroids[shared_genes].values  # shape: (n_human, n_genes)

    n_mouse = mouse_vals.shape[0]
    n_human = human_vals.shape[0]

    # Rank-transform each row (cluster) across genes — Spearman = Pearson on ranks
    def _rank_rows(mat: np.ndarray) -> np.ndarray:
        ranked = np.apply_along_axis(rankdata, axis=1, arr=mat)
        return ranked

    mouse_ranked = _rank_rows(mouse_vals)
    human_ranked = _rank_rows(human_vals)

    # Mean-center for efficient correlation via dot product
    mouse_centered = mouse_ranked - mouse_ranked.mean(axis=1, keepdims=True)
    human_centered = human_ranked - human_ranked.mean(axis=1, keepdims=True)

    mouse_norm = np.linalg.norm(mouse_centered, axis=1, keepdims=True)
    human_norm = np.linalg.norm(human_centered, axis=1, keepdims=True)

    # Avoid division by zero
    mouse_norm[mouse_norm == 0] = 1.0
    human_norm[human_norm == 0] = 1.0

    mouse_unit = mouse_centered / mouse_norm
    human_unit = human_centered / human_norm

    corr_matrix = mouse_unit @ human_unit.T  # (n_mouse, n_human)

    return pd.DataFrame(
        corr_matrix,
        index=mouse_centroids.index,
        columns=human_centroids.index,
    )


def find_reciprocal_best_hits(
    corr_matrix: pd.DataFrame,
    threshold: float | None = None,
) -> pd.DataFrame:
    """Find reciprocal best hit (RBH) pairs from a Spearman correlation matrix.

    A pair (mouse_cluster, human_cluster) is an RBH when:
    - human_cluster is the best-matching human cluster for mouse_cluster, AND
    - mouse_cluster is the best-matching mouse cluster for human_cluster.

    Parameters
    ----------
    corr_matrix : pd.DataFrame
        (n_mouse) x (n_human) Spearman correlation matrix.
        Index = mouse cluster labels, columns = human cluster labels.
    threshold : float or None, optional
        Minimum correlation required for a pair to be marked ``is_rbh=True``.
        If None, no threshold is applied beyond the RBH criterion itself.

    Returns
    -------
    pd.DataFrame
        One row per RBH candidate with columns:
        ``mouse_cluster``, ``human_cluster``, ``correlation``,
        ``is_rbh`` (bool), ``method`` (str).
    """
    vals = corr_matrix.values
    mouse_labels = corr_matrix.index.tolist()
    human_labels = corr_matrix.columns.tolist()

    # Best human match for each mouse cluster
    best_human_for_mouse = np.argmax(vals, axis=1)  # shape: (n_mouse,)
    # Best mouse match for each human cluster
    best_mouse_for_human = np.argmax(vals, axis=0)  # shape: (n_human,)

    records = []
    for mi, mouse_cl in enumerate(mouse_labels):
        hi = best_human_for_mouse[mi]
        human_cl = human_labels[hi]
        corr_val = float(vals[mi, hi])

        is_rbh_topology = best_mouse_for_human[hi] == mi
        is_rbh_thresh = (threshold is None) or (corr_val >= threshold)
        is_rbh = bool(is_rbh_topology and is_rbh_thresh)

        records.append({
            "mouse_cluster": mouse_cl,
            "human_cluster": human_cl,
            "correlation": corr_val,
            "is_rbh": is_rbh,
            "method": "spearman_rbh",
        })

    return pd.DataFrame(records)


def determine_rbh_threshold(
    corr_matrix: pd.DataFrame,
    method: str = "permutation",
    n_perm: int = 1000,
    seed: int = 42,
) -> float:
    """Determine the correlation threshold for RBH significance via permutation.

    Shuffles cluster labels ``n_perm`` times, computes the best-hit correlation
    for each permuted matrix, and returns the 95th percentile of that null
    distribution as the threshold.

    Parameters
    ----------
    corr_matrix : pd.DataFrame
        (n_mouse) x (n_human) Spearman correlation matrix.
    method : str
        Currently only ``"permutation"`` is supported.
    n_perm : int
        Number of permutations to draw.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    float
        95th-percentile null RBH correlation (suggested threshold).
    """
    if method != "permutation":
        raise ValueError(f"Unsupported method: {method!r}. Use 'permutation'.")

    rng = np.random.default_rng(seed)
    vals = corr_matrix.values.copy()
    n_mouse, n_human = vals.shape

    null_rbh_corrs = []
    for _ in range(n_perm):
        # Shuffle rows (mouse clusters) of the correlation matrix
        perm_idx = rng.permutation(n_mouse)
        perm_vals = vals[perm_idx, :]

        # Collect best-hit correlation for each mouse cluster under permutation
        best_corrs = perm_vals.max(axis=1)
        null_rbh_corrs.extend(best_corrs.tolist())

    threshold = float(np.percentile(null_rbh_corrs, 95))
    # Clip to valid correlation range
    return float(np.clip(threshold, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Round 3: MetaNeighbor AUROC
# ---------------------------------------------------------------------------

def compute_metaneighbor_auroc(
    expr_df: pd.DataFrame,
    labels_mouse: pd.Series,
    labels_human: pd.Series,
    species: pd.Series,
) -> pd.DataFrame:
    """Compute MetaNeighbor AUROC between mouse and human cluster labels.

    Algorithm (Python-native MetaNeighbor):
    1. Rank-transform each cell's expression across genes (Spearman approximation).
    2. Compute cell-cell Spearman correlations as inner products of rank vectors.
    3. For each mouse cluster X and each human cluster Y, use correlation scores
       from mouse cluster X cells to all human cells as the ranking signal,
       and compute AUROC vs the binary label "is human cell in cluster Y".
    4. Average AUROC across cells within each mouse cluster.

    Parameters
    ----------
    expr_df : pd.DataFrame
        Cell x gene expression matrix covering all cells (mouse + human).
        Index = cell barcodes, columns = gene names.
    labels_mouse : pd.Series
        Cluster labels for mouse cells. Index = mouse cell barcodes.
    labels_human : pd.Series
        Cluster labels for human cells. Index = human cell barcodes.
    species : pd.Series
        Species assignment (``"mouse"`` or ``"human"``) for each cell.
        Index = all cell barcodes.

    Returns
    -------
    pd.DataFrame
        (n_mouse_clusters) x (n_human_clusters) AUROC matrix.
        Index = mouse cluster labels, columns = human cluster labels.
        AUROC > 0.5 indicates mouse cluster is more similar to the matching
        human cluster than to other human cells.
    """
    mouse_cells = labels_mouse.index.tolist()
    human_cells = labels_human.index.tolist()
    all_cells = mouse_cells + human_cells

    # Rank-transform each gene across cells (rank along axis=0, i.e., per column).
    # This yields cell-cell Spearman correlation when combined with a Pearson
    # inner product — matching the original MetaNeighbor formulation.
    expr_sub = expr_df.loc[all_cells].values.astype(float)  # (n_cells, n_genes)
    ranked = np.apply_along_axis(rankdata, axis=0, arr=expr_sub)  # rank cells per gene

    # Mean-center per gene so that inner product gives Pearson on ranks = Spearman
    centered = ranked - ranked.mean(axis=0, keepdims=True)

    # Unit-normalize per cell for cosine similarity
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit_ranked = centered / norms

    n_mouse = len(mouse_cells)
    n_human = len(human_cells)

    # Cell-cell correlation: only mouse→human block needed
    # mouse_unit: (n_mouse, n_genes), human_unit: (n_human, n_genes)
    mouse_unit = unit_ranked[:n_mouse, :]
    human_unit = unit_ranked[n_mouse:, :]
    corr_mouse_human = mouse_unit @ human_unit.T  # (n_mouse, n_human)

    mouse_cluster_labels = np.array([labels_mouse[c] for c in mouse_cells])
    human_cluster_labels = np.array([labels_human[c] for c in human_cells])

    unique_mouse = sorted(set(mouse_cluster_labels))
    unique_human = sorted(set(human_cluster_labels))

    auroc_matrix = np.zeros((len(unique_mouse), len(unique_human)))

    for mi, mc in enumerate(unique_mouse):
        mc_mask = mouse_cluster_labels == mc  # boolean (n_mouse,)
        mc_indices = np.where(mc_mask)[0]

        for hi, hc in enumerate(unique_human):
            hc_mask = (human_cluster_labels == hc).astype(int)  # binary label

            # For each cell in mouse cluster mc, compute AUROC using correlation
            # scores to human cells as ranking signal vs hc_mask as truth
            cell_aurocs = []
            for cell_idx in mc_indices:
                scores = corr_mouse_human[cell_idx, :]  # (n_human,)
                # Need at least one positive and one negative
                if hc_mask.sum() == 0 or hc_mask.sum() == n_human:
                    cell_aurocs.append(0.5)
                    continue
                try:
                    auc = roc_auc_score(hc_mask, scores)
                except ValueError:
                    auc = 0.5
                cell_aurocs.append(auc)

            auroc_matrix[mi, hi] = float(np.mean(cell_aurocs))

    return pd.DataFrame(
        auroc_matrix,
        index=unique_mouse,
        columns=unique_human,
    )
