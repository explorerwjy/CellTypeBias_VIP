"""Ephys feature harmonization and cross-species comparison."""

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Feature constants
# ---------------------------------------------------------------------------

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

ISI_FEATURES = [f"bin_{i}_isi_ms" for i in range(12)]
CV_ISI_FEATURES = [f"bin_{i}_cv_isi" for i in range(12)]

LOG_FEATURES = ["spike_frequency_Hz", "current_threshold_pA"] + ISI_FEATURES

# All features extracted per cell
ALL_FEATURES = BIO_FEATURES + ISI_FEATURES + CV_ISI_FEATURES


# ---------------------------------------------------------------------------
# Sweep selection
# ---------------------------------------------------------------------------


def select_representative_sweep(sweeps_df, target_min=5, target_max=60):
    """Select the most representative sweep for a cell.

    Selects the sweep whose firing rate is in [target_min, target_max] Hz and
    closest to 20 Hz. Falls back to the lowest-current spiking sweep when no
    sweep is in range.

    Parameters
    ----------
    sweeps_df : pd.DataFrame
        DataFrame with one row per sweep. Must contain ``spike_frequency_Hz``
        and ``avg_injected_current_pA`` columns.
    target_min : float
        Lower bound of the target firing rate range (Hz).
    target_max : float
        Upper bound of the target firing rate range (Hz).

    Returns
    -------
    int
        Row index (iloc position) of the selected sweep.
    """
    TARGET_FREQ = 20.0

    freqs = sweeps_df["spike_frequency_Hz"].values
    currents = sweeps_df["avg_injected_current_pA"].values

    # Sweeps in the target firing-rate window
    in_range = (freqs >= target_min) & (freqs <= target_max)

    if in_range.any():
        candidates = np.where(in_range)[0]
        distances = np.abs(freqs[candidates] - TARGET_FREQ)
        return int(candidates[np.argmin(distances)])

    # Fallback: lowest-current spiking sweep (frequency > 0)
    spiking = np.where(freqs > 0)[0]
    if len(spiking) == 0:
        # No spiking sweep at all — return the first row
        return 0

    lowest_current_idx = spiking[np.argmin(currents[spiking])]
    return int(lowest_current_idx)


# ---------------------------------------------------------------------------
# Feature aggregation
# ---------------------------------------------------------------------------


def aggregate_cell_features(processed_dir, species_name):
    """Walk a processed DANDI directory and aggregate per-cell ephys features.

    For each cell directory containing an ``analysis.csv``, selects the
    representative sweep (via :func:`select_representative_sweep`) and
    extracts :data:`ALL_FEATURES` (BIO + ISI + CV_ISI).  Features missing from
    a given CSV are filled with NaN.

    Parameters
    ----------
    processed_dir : str or Path
        Root directory to walk.  Each subdirectory may contain an
        ``analysis.csv`` produced by the EphysSumStats pipeline.
    species_name : str
        Label added to every row as the ``species`` column.

    Returns
    -------
    pd.DataFrame
        One row per cell.  Columns: ALL_FEATURES + ``cell_id`` + ``species``.
    """
    processed_dir = Path(processed_dir)
    records = []

    for root, _dirs, files in os.walk(processed_dir):
        if "analysis.csv" not in files:
            continue

        csv_path = Path(root) / "analysis.csv"
        try:
            sweeps = pd.read_csv(csv_path)
        except Exception as exc:
            warnings.warn(f"Could not read {csv_path}: {exc}")
            continue

        if sweeps.empty or "spike_frequency_Hz" not in sweeps.columns:
            continue

        idx = select_representative_sweep(sweeps)
        row = sweeps.iloc[idx]

        record = {"cell_id": Path(root).name, "species": species_name}
        for feat in ALL_FEATURES:
            record[feat] = row[feat] if feat in sweeps.columns else np.nan

        records.append(record)

    if not records:
        return pd.DataFrame(columns=["cell_id", "species"] + ALL_FEATURES)

    df = pd.DataFrame(records)
    # Reorder columns: cell_id, species, then features
    feature_cols = [c for c in ALL_FEATURES if c in df.columns]
    return df[["cell_id", "species"] + feature_cols]


# ---------------------------------------------------------------------------
# Within-species z-scoring
# ---------------------------------------------------------------------------


def zscore_within_species(features_df, species_labels):
    """Z-score each feature within each species group.

    Parameters
    ----------
    features_df : pd.DataFrame
        Numeric feature matrix (cells × features).
    species_labels : array-like
        Species label for each cell (same length as ``features_df``).

    Returns
    -------
    pd.DataFrame
        Z-scored DataFrame with the same shape and index/columns as
        ``features_df``.  Features that are constant within a species are set
        to 0 rather than NaN (epsilon guard applied).
    """
    EPSILON = 1e-10

    species = np.asarray(species_labels)
    result = features_df.copy().astype(float)

    for sp in np.unique(species):
        mask = species == sp
        subset = result.loc[mask]
        mu = subset.mean(axis=0)
        sigma = subset.std(axis=0, ddof=0)
        # Guard against constant features
        sigma = sigma.where(sigma > EPSILON, other=EPSILON)
        result.loc[mask] = (subset - mu) / sigma

    return result


# ---------------------------------------------------------------------------
# Permutation test
# ---------------------------------------------------------------------------


def permutation_test_cluster_similarity(
    features_df, cluster_labels, species_labels, n_perm=1000, seed=42
):
    """Test whether matching cross-species clusters are more similar than chance.

    For each cluster that contains cells from both species, the Euclidean
    distance between the two species centroids is computed.  The observed
    statistic is the mean of these within-cluster cross-species distances
    (smaller = more similar).

    A null distribution is built by permuting cluster labels and recomputing
    the statistic ``n_perm`` times.  The one-sided p-value is the fraction of
    null values ≤ the observed statistic (i.e. the fraction of times the null
    is *at least as similar* as the observed, testing that true matching is
    more similar than shuffled).

    Parameters
    ----------
    features_df : pd.DataFrame
        Numeric feature matrix (cells × features).
    cluster_labels : array-like
        Cluster label for each cell.
    species_labels : array-like
        Species label for each cell.
    n_perm : int
        Number of permutations.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    float
        Empirical one-sided p-value in [0, 1].
    """
    rng = np.random.default_rng(seed)

    X = features_df.values.astype(float)
    clusters = np.asarray(cluster_labels)
    species = np.asarray(species_labels)

    def _mean_cross_species_distance(clust_arr):
        distances = []
        for cl in np.unique(clust_arr):
            mask_cl = clust_arr == cl
            unique_sp = np.unique(species[mask_cl])
            if len(unique_sp) < 2:
                continue
            sp0_cells = mask_cl & (species == unique_sp[0])
            sp1_cells = mask_cl & (species == unique_sp[1])
            # Ignore clusters with no valid (non-NaN) cells in either species
            if sp0_cells.sum() == 0 or sp1_cells.sum() == 0:
                continue
            centroid0 = np.nanmean(X[sp0_cells], axis=0)
            centroid1 = np.nanmean(X[sp1_cells], axis=0)
            dist = np.sqrt(np.nansum((centroid0 - centroid1) ** 2))
            distances.append(dist)
        if len(distances) == 0:
            return np.nan
        return float(np.mean(distances))

    observed = _mean_cross_species_distance(clusters)

    null_values = np.empty(n_perm)
    for i in range(n_perm):
        shuffled = rng.permutation(clusters)
        null_values[i] = _mean_cross_species_distance(shuffled)

    if np.isnan(observed):
        return np.nan

    # One-sided p-value: fraction of null values <= observed
    p_value = float(np.mean(null_values <= observed))
    return p_value
