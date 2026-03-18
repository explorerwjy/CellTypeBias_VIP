"""Path concordance validation: chaining, agreement metrics, Cohen's kappa."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score


# ---------------------------------------------------------------------------
# Path chaining
# ---------------------------------------------------------------------------

def build_indirect_path(
    patchseq_to_mouse: pd.DataFrame,
    mouse_to_human: pd.DataFrame,
) -> pd.DataFrame:
    """Chain mouse Patch-seq → mouse SC → human SC to build the indirect path.

    Parameters
    ----------
    patchseq_to_mouse:
        DataFrame with at least columns ``cell_id`` and ``mouse_cluster``.
        Represents direct Patch-seq → mouse single-cell cluster assignments.
    mouse_to_human:
        DataFrame with at least columns ``mouse_cluster`` and ``human_cluster``.
        Represents the cross-species mouse SC → human SC correspondence
        (e.g. from Module A RBH pairs or Module B Harmony transfer).

    Returns
    -------
    pd.DataFrame
        A copy of *patchseq_to_mouse* with an additional column
        ``indirect_human_cluster``.  Cells whose ``mouse_cluster`` has no match
        in *mouse_to_human* will have ``NaN`` in that column.
    """
    mapping = mouse_to_human[["mouse_cluster", "human_cluster"]].copy()
    mapping = mapping.rename(columns={"human_cluster": "indirect_human_cluster"})

    result = patchseq_to_mouse.copy()
    result = result.merge(mapping, on="mouse_cluster", how="left")
    return result


# ---------------------------------------------------------------------------
# Agreement metrics
# ---------------------------------------------------------------------------

def compute_concordance(
    direct_labels: pd.Series,
    indirect_labels: pd.Series,
) -> float:
    """Compute fraction of cells where direct and indirect labels agree.

    Only cells where *both* labels are non-NaN are considered.

    Parameters
    ----------
    direct_labels:
        Series of cluster labels from the direct mapping path.
    indirect_labels:
        Series of cluster labels from the indirect mapping path.
        Must have the same index as *direct_labels*.

    Returns
    -------
    float
        Fraction of valid-pair cells where the two labels match (0.0 – 1.0),
        or ``float('nan')`` if there are no valid (non-NaN) pairs.
    """
    valid_mask = direct_labels.notna() & indirect_labels.notna()
    n_valid = valid_mask.sum()
    if n_valid == 0:
        return float("nan")
    agree = (direct_labels[valid_mask] == indirect_labels[valid_mask]).sum()
    return float(agree / n_valid)


def compute_cohens_kappa(
    labels_a: pd.Series,
    labels_b: pd.Series,
) -> float:
    """Compute Cohen's kappa between two label series.

    Rows where either series is NaN are dropped before computing kappa.

    Parameters
    ----------
    labels_a:
        First set of categorical labels.
    labels_b:
        Second set of categorical labels.  Must share the index of *labels_a*.

    Returns
    -------
    float
        Cohen's kappa in the range [-1, 1].  Returns ``float('nan')`` if fewer
        than two valid pairs remain after filtering NaN values.
    """
    valid_mask = labels_a.notna() & labels_b.notna()
    a_clean = labels_a[valid_mask]
    b_clean = labels_b[valid_mask]
    if len(a_clean) < 2:
        return float("nan")
    return float(cohen_kappa_score(a_clean, b_clean))


# ---------------------------------------------------------------------------
# Multi-level concordance
# ---------------------------------------------------------------------------

def concordance_at_levels(
    direct_clusters: pd.Series,
    indirect_clusters: pd.Series,
    cluster_to_subclass: dict[str, str],
    cluster_to_supercluster: dict[str, str],
) -> dict[str, float]:
    """Compute concordance and Cohen's kappa at cluster, subclass, and supercluster levels.

    Parameters
    ----------
    direct_clusters:
        Per-cell cluster labels from the direct mapping path
        (patch-seq → human SC).
    indirect_clusters:
        Per-cell cluster labels from the indirect mapping path
        (patch-seq → mouse SC → human SC).  Must share the index of
        *direct_clusters*.
    cluster_to_subclass:
        Mapping from cluster ID to subclass label.
    cluster_to_supercluster:
        Mapping from cluster ID to supercluster label.

    Returns
    -------
    dict
        Keys:
        ``cluster_concordance``, ``cluster_kappa``,
        ``subclass_concordance``, ``subclass_kappa``,
        ``supercluster_concordance``, ``supercluster_kappa``.
        Values are floats (or NaN when insufficient data).
    """
    # --- Cluster level ---
    cluster_concordance = compute_concordance(direct_clusters, indirect_clusters)
    cluster_kappa = compute_cohens_kappa(direct_clusters, indirect_clusters)

    # --- Subclass level ---
    def _map_series(s: pd.Series, mapping: dict) -> pd.Series:
        return s.map(mapping)  # unmapped keys → NaN

    direct_subclass = _map_series(direct_clusters, cluster_to_subclass)
    indirect_subclass = _map_series(indirect_clusters, cluster_to_subclass)
    subclass_concordance = compute_concordance(direct_subclass, indirect_subclass)
    subclass_kappa = compute_cohens_kappa(direct_subclass, indirect_subclass)

    # --- Supercluster level ---
    direct_super = _map_series(direct_clusters, cluster_to_supercluster)
    indirect_super = _map_series(indirect_clusters, cluster_to_supercluster)
    supercluster_concordance = compute_concordance(direct_super, indirect_super)
    supercluster_kappa = compute_cohens_kappa(direct_super, indirect_super)

    return {
        "cluster_concordance": cluster_concordance,
        "cluster_kappa": cluster_kappa,
        "subclass_concordance": subclass_concordance,
        "subclass_kappa": subclass_kappa,
        "supercluster_concordance": supercluster_concordance,
        "supercluster_kappa": supercluster_kappa,
    }
