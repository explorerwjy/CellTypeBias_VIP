"""Cross-species cluster correspondence: pseudobulk, Spearman RBH, MetaNeighbor AUROC.

The MetaNeighbor implementation in this module is a Python port of Bakken
et al. 2021's BICCN_M1_Evo R helpers (`metaneighbor.R::compute_best_hits`),
which is itself a re-implementation of MetaNeighbor (Crow et al. 2018).
We mirror Bakken 2021 because that is the paper we cite for cross-species
interneuron comparison.

Reference:
  https://github.com/AllenInstitute/BICCN_M1_Evo/blob/master/MetaNeighbor/metaneighbor.R
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata


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
        One row **per mouse cluster** (regardless of whether it is an RBH).
        Each row contains the best-matching human cluster for that mouse cluster,
        plus ``is_rbh=True`` if the pair passes both the reciprocal topology check
        and the optional correlation threshold, ``is_rbh=False`` otherwise.
        Columns: ``mouse_cluster``, ``human_cluster``, ``correlation``,
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

    # Collect the per-permutation maximum best-hit correlation across all mouse
    # clusters.  Shuffling columns independently per row destroys the real
    # structure while preserving the marginal distribution of each row.
    perm_maxima = []
    for _ in range(n_perm):
        perm_vals = vals.copy()
        for row in range(n_mouse):
            perm_vals[row] = perm_vals[row, rng.permutation(n_human)]

        # Best-hit correlation for each mouse cluster, then take the max across
        # clusters so that the null reflects the family-wise null (one threshold
        # controls the experiment as a whole rather than being inflated by pooling
        # n_mouse × n_perm values from a single flat list).
        best_corrs = perm_vals.max(axis=1)
        perm_maxima.append(best_corrs.max())

    threshold = float(np.percentile(perm_maxima, 95))
    # Clip to valid correlation range
    return float(np.clip(threshold, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Round 3: MetaNeighbor AUROC (canonical, after Bakken 2021)
# ---------------------------------------------------------------------------
#
# Python port of Bakken et al. 2021 BICCN_M1_Evo `metaneighbor.R` helpers.
# The reference R implementation is reproduced below for traceability:
#
#     normalize_cols = function(M, ranked = TRUE) {
#       M = as.matrix(M)
#       if (ranked) M = matrixStats::colRanks(M, ties.method = "average",
#                                             preserveShape = TRUE)
#       return(scale_cols(M))
#     }
#
#     compute_aurocs = function(predictors, label_matrix) {
#       n_positives = colSums(label_matrix)
#       n_negatives = nrow(label_matrix) - n_positives
#       ranks = matrixStats::colRanks(predictors, ties.method = "average",
#                                     preserveShape=TRUE)
#       sum_of_positive_ranks = crossprod(label_matrix, ranks)
#       result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) /
#                n_negatives
#       return(result)
#     }
#
#     compute_best_hits = function(dataset, labels, one_vs_all=TRUE) {
#       normalized_data = normalize_cols(assay(dataset))
#       colnames(normalized_data) = paste(dataset$study_id, labels, sep = "|")
#       voter_id = design_matrix(colnames(normalized_data))
#       voters = normalized_data %*% voter_id
#       result = c()
#       for (study in unique(dataset$study_id)) {
#         candidates = normalized_data[, dataset$study_id == study]
#         votes = crossprod(candidates, voters)
#         aurocs = compute_aurocs(votes, design_matrix(rownames(votes)))
#         result = rbind(result, aurocs)
#       }
#       result = result[, rownames(result)]
#       return(result)
#     }


def _normalize_cols_metaneighbor(M: np.ndarray) -> np.ndarray:
    """Per-cell rank-normalize expression (Bakken 2021 ``normalize_cols``).

    Parameters
    ----------
    M : np.ndarray
        ``(n_genes, n_cells)`` raw expression matrix. Each column is one cell.

    Returns
    -------
    np.ndarray
        Same shape; each column has been replaced by its gene ranks (average
        ties), then mean-centered and L2-normalized so that the inner product
        of two columns equals their Spearman correlation.
    """
    M = np.asarray(M, dtype=float)
    # Rank gene expression within each cell (rank down each column)
    ranked = np.apply_along_axis(rankdata, axis=0, arr=M)
    centered = ranked - ranked.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centered, axis=0, keepdims=True)
    norms[norms == 0] = 1.0
    return centered / norms


def _design_matrix_metaneighbor(labels: np.ndarray) -> tuple[np.ndarray, list]:
    """One-hot encoding (Bakken 2021 ``design_matrix``).

    Parameters
    ----------
    labels : np.ndarray
        ``(n_cells,)`` array of categorical labels (any hashable).

    Returns
    -------
    tuple
        ``(M, categories)`` where ``M`` is ``(n_cells, n_categories)`` binary
        and ``categories`` is the column order.
    """
    labels = np.asarray(labels)
    categories = sorted(set(labels.tolist()))
    n = len(labels)
    k = len(categories)
    M = np.zeros((n, k), dtype=float)
    cat_to_idx = {c: i for i, c in enumerate(categories)}
    for i, lbl in enumerate(labels):
        M[i, cat_to_idx[lbl]] = 1.0
    return M, categories


def _compute_aurocs_mannwhitney(
    predictors: np.ndarray, label_matrix: np.ndarray
) -> np.ndarray:
    """Vectorized Mann-Whitney AUROC (Bakken 2021 ``compute_aurocs``).

    Computes one AUROC for each (label class, predictor column) pair using
    the Mann-Whitney U statistic in closed form (no per-cell loops):

        auroc = (sum_of_positive_ranks / n_pos - (n_pos+1)/2) / n_neg

    Parameters
    ----------
    predictors : np.ndarray
        ``(n_obs, n_columns)`` matrix of scores. Each column is a separate
        ranking signal (e.g., similarity to one voter group).
    label_matrix : np.ndarray
        ``(n_obs, n_classes)`` binary matrix encoding which observations
        belong to each class.

    Returns
    -------
    np.ndarray
        ``(n_classes, n_columns)`` AUROC matrix. Entry ``[c, j]`` is the
        AUROC for ranking class-``c`` observations using column-``j`` scores.
    """
    predictors = np.asarray(predictors, dtype=float)
    label_matrix = np.asarray(label_matrix, dtype=float)
    n_obs = label_matrix.shape[0]
    n_positives = label_matrix.sum(axis=0)  # (n_classes,)
    n_negatives = n_obs - n_positives
    # Rank each column independently (with average ties)
    ranks = np.apply_along_axis(rankdata, axis=0, arr=predictors)  # (n_obs, n_columns)
    # Sum of ranks of positive observations, per (class, column)
    sum_of_positive_ranks = label_matrix.T @ ranks  # (n_classes, n_columns)
    n_pos_col = n_positives[:, None]
    n_neg_col = n_negatives[:, None]
    # Avoid division by zero — return 0.5 (chance) for degenerate cases
    safe_pos = np.where(n_pos_col == 0, 1.0, n_pos_col)
    safe_neg = np.where(n_neg_col == 0, 1.0, n_neg_col)
    auroc = (sum_of_positive_ranks / safe_pos - (n_pos_col + 1) / 2.0) / safe_neg
    auroc[(n_pos_col == 0).ravel(), :] = 0.5
    auroc[:, ...][np.broadcast_to(n_neg_col == 0, auroc.shape)] = 0.5
    return auroc


def compute_best_hits_metaneighbor(
    expr_genes_x_cells: np.ndarray,
    labels: np.ndarray,
    study_ids: np.ndarray,
    cell_names: list | None = None,
) -> pd.DataFrame:
    """Cross-study MetaNeighbor AUROC matrix (Bakken 2021 ``compute_best_hits``).

    Implements the leave-one-study-out cross-validation procedure: for each
    study, its cells are used as test queries and **all** cells (including
    other studies) provide voter centroids; the resulting per-test-cluster
    AUROC rows are concatenated across studies.

    Parameters
    ----------
    expr_genes_x_cells : np.ndarray
        ``(n_genes, n_cells)`` raw expression matrix.
    labels : np.ndarray
        ``(n_cells,)`` cluster label for each cell.
    study_ids : np.ndarray
        ``(n_cells,)`` study/dataset identifier (e.g., ``"mouse"``, ``"human"``).
    cell_names : list, optional
        Cell barcodes for traceability. Not used in computation.

    Returns
    -------
    pd.DataFrame
        Square AUROC matrix indexed by ``"<study>|<cluster>"``. Entry
        ``[A, B]`` is the AUROC of group ``B`` voters discriminating cells
        of group ``A`` from other cells in group ``A``'s test study.
        Following Bakken's R code, the result is column-reordered to match
        the row order so that the matrix is square.
    """
    expr_genes_x_cells = np.asarray(expr_genes_x_cells, dtype=float)
    labels = np.asarray(labels).astype(str)
    study_ids = np.asarray(study_ids).astype(str)

    # Step 1: rank-normalize all cells (matches normalize_cols)
    normalized = _normalize_cols_metaneighbor(expr_genes_x_cells)

    # Step 2: build (study|cluster) labels and voter centroids
    cell_groups = np.array([f"{s}|{c}" for s, c in zip(study_ids, labels)])
    voter_id, voter_categories = _design_matrix_metaneighbor(cell_groups)
    # voters[:, j] = sum of normalized cell vectors in voter group j  (n_genes, n_groups)
    voters = normalized @ voter_id

    # Step 3: leave-one-study-out — concatenate per-study AUROC blocks
    result_rows = []
    result_row_names: list[str] = []
    for test_study in sorted(set(study_ids.tolist())):
        test_mask = study_ids == test_study
        candidates = normalized[:, test_mask]  # (n_genes, n_test_cells)
        test_groups = cell_groups[test_mask]
        votes = candidates.T @ voters  # (n_test_cells, n_voter_groups)
        test_label_matrix, test_categories = _design_matrix_metaneighbor(test_groups)
        aurocs = _compute_aurocs_mannwhitney(votes, test_label_matrix)
        # aurocs shape: (n_test_clusters, n_voter_groups)
        result_rows.append(aurocs)
        result_row_names.extend(test_categories)

    full = np.vstack(result_rows)
    df = pd.DataFrame(full, index=result_row_names, columns=voter_categories)
    # Match Bakken's `result = result[, rownames(result)]` — reorder columns to
    # the row order so the matrix is square and aligned.
    df = df[result_row_names]
    return df


def compute_metaneighbor_auroc(
    expr_df: pd.DataFrame,
    labels_mouse: pd.Series,
    labels_human: pd.Series,
    species: pd.Series,
) -> pd.DataFrame:
    """MetaNeighbor AUROC matrix between mouse and human cluster labels.

    Wraps :func:`compute_best_hits_metaneighbor` to preserve the legacy
    function signature used elsewhere in the pipeline. Internally:

    1. Builds a single ``"<species>|<cluster>"`` label per cell.
    2. Calls :func:`compute_best_hits_metaneighbor` to obtain the full LOSO
       AUROC matrix.
    3. Extracts the cross-species block (rows = ``mouse|*``, columns =
       ``human|*``) and strips the ``"mouse|"`` / ``"human|"`` prefixes.

    Parameters
    ----------
    expr_df : pd.DataFrame
        ``(n_cells, n_genes)`` expression. Index = cell barcodes.
    labels_mouse : pd.Series
        Cluster labels for mouse cells; index aligns with ``expr_df``.
    labels_human : pd.Series
        Cluster labels for human cells; index aligns with ``expr_df``.
    species : pd.Series
        ``"mouse"`` or ``"human"`` per cell; index aligns with ``expr_df``.

    Returns
    -------
    pd.DataFrame
        ``(n_mouse_clusters, n_human_clusters)`` AUROC. Higher values mean
        the mouse cluster more strongly resembles the human cluster than
        other human cells in the same test study.
    """
    mouse_cells = species[species == "mouse"].index.intersection(labels_mouse.index).tolist()
    human_cells = species[species == "human"].index.intersection(labels_human.index).tolist()
    all_cells = mouse_cells + human_cells

    # expr_df is (n_cells, n_genes); transpose to (n_genes, n_cells) for the
    # Bakken-style routines.
    expr_sub = expr_df.loc[all_cells].values.astype(float).T

    cluster_labels = np.array(
        [str(labels_mouse[c]) for c in mouse_cells]
        + [str(labels_human[c]) for c in human_cells]
    )
    study_ids = np.array(["mouse"] * len(mouse_cells) + ["human"] * len(human_cells))

    full_auroc = compute_best_hits_metaneighbor(expr_sub, cluster_labels, study_ids)

    # Extract the (mouse rows × human columns) block and strip the prefixes.
    mouse_rows = [r for r in full_auroc.index if r.startswith("mouse|")]
    human_cols = [c for c in full_auroc.columns if c.startswith("human|")]
    block = full_auroc.loc[mouse_rows, human_cols].copy()
    block.index = [r.split("|", 1)[1] for r in mouse_rows]
    block.columns = [c.split("|", 1)[1] for c in human_cols]
    block = block.sort_index().sort_index(axis=1)
    return block


