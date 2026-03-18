"""Tests for ephys_harmonization module (Module D)."""

import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.ephys_harmonization import (
    select_representative_sweep,
    zscore_within_species,
    permutation_test_cluster_similarity,
)


# ---------------------------------------------------------------------------
# select_representative_sweep
# ---------------------------------------------------------------------------


class TestSelectRepresentativeSweep:
    """Tests for select_representative_sweep."""

    def _make_sweeps(self, freqs, currents=None):
        """Helper: build a sweeps DataFrame."""
        if currents is None:
            currents = list(range(len(freqs)))
        return pd.DataFrame(
            {
                "spike_frequency_Hz": freqs,
                "avg_injected_current_pA": currents,
            }
        )

    def test_select_representative_sweep(self):
        """Sweeps with rates [0, 0, 10, 25, 60, 100] Hz; selection must be in [5, 60]."""
        freqs = [0, 0, 10, 25, 60, 100]
        currents = [0, 10, 20, 30, 40, 50]
        df = self._make_sweeps(freqs, currents)

        idx = select_representative_sweep(df, target_min=5, target_max=60)

        selected_freq = freqs[idx]
        assert 5 <= selected_freq <= 60, (
            f"Selected sweep {idx} has frequency {selected_freq} Hz outside [5, 60]"
        )

    def test_select_representative_sweep_closest_to_20(self):
        """Among valid sweeps, picks the one closest to 20 Hz."""
        freqs = [0, 0, 10, 25, 60, 100]
        df = self._make_sweeps(freqs)

        idx = select_representative_sweep(df, target_min=5, target_max=60)

        # 25 Hz is closer to 20 Hz than 10 Hz (|25-20|=5 vs |10-20|=10)
        # 60 Hz edge case: |60-20|=40, excluded since it equals target_max (inclusive)
        assert freqs[idx] == 25, (
            f"Expected sweep at 25 Hz (closest to 20), got {freqs[idx]} Hz"
        )

    def test_select_representative_sweep_fallback(self):
        """No sweep in [5, 60] Hz range — picks lowest-current spiking sweep."""
        # Only sweeps outside the target range
        freqs = [0, 0, 0, 80, 120, 200]
        currents = [0, 10, 20, 100, 150, 200]
        df = self._make_sweeps(freqs, currents)

        idx = select_representative_sweep(df, target_min=5, target_max=60)

        # Fallback: lowest-current spiking sweep
        # Spiking sweeps: indices 3, 4, 5 with currents 100, 150, 200
        # Lowest current among spiking = 100 at index 3
        assert freqs[idx] > 0, "Fallback should select a spiking sweep"
        # Must be the spiking sweep with the lowest current
        spiking_currents = [(i, currents[i]) for i in range(len(freqs)) if freqs[i] > 0]
        min_current = min(c for _, c in spiking_currents)
        expected_indices = [i for i, c in spiking_currents if c == min_current]
        assert idx in expected_indices, (
            f"Fallback selected index {idx} (current={currents[idx]}), "
            f"expected minimum-current spiking sweep at {expected_indices}"
        )

    def test_select_no_spiking_sweeps(self):
        """If no sweeps have spikes, returns index 0."""
        freqs = [0, 0, 0, 0]
        df = self._make_sweeps(freqs)

        idx = select_representative_sweep(df, target_min=5, target_max=60)

        assert idx == 0, f"Expected index 0 for all-silent cell, got {idx}"


# ---------------------------------------------------------------------------
# zscore_within_species
# ---------------------------------------------------------------------------


class TestZscoreWithinSpecies:
    """Tests for zscore_within_species."""

    def test_zscore_within_species(self):
        """3 mouse + 3 human cells; each species group should have mean ~0."""
        data = pd.DataFrame(
            {
                "feature_A": [1.0, 2.0, 3.0, 10.0, 20.0, 30.0],
                "feature_B": [5.0, 6.0, 7.0, 50.0, 60.0, 70.0],
            }
        )
        species = ["mouse", "mouse", "mouse", "human", "human", "human"]

        result = zscore_within_species(data, species)

        # Mean within each species should be ~0
        for sp, indices in [("mouse", [0, 1, 2]), ("human", [3, 4, 5])]:
            sp_means = result.iloc[indices].mean()
            for feat, mean_val in sp_means.items():
                assert abs(mean_val) < 1e-9, (
                    f"Species '{sp}', feature '{feat}': expected mean ~0, got {mean_val}"
                )

    def test_zscore_standardizes_variance(self):
        """Each species group should have std ~1 (population std) after z-scoring."""
        data = pd.DataFrame(
            {
                "feature_A": [1.0, 2.0, 3.0, 100.0, 200.0, 300.0],
            }
        )
        species = ["mouse", "mouse", "mouse", "human", "human", "human"]

        result = zscore_within_species(data, species)

        for sp, indices in [("mouse", [0, 1, 2]), ("human", [3, 4, 5])]:
            sp_std = result.iloc[indices]["feature_A"].std(ddof=0)
            assert abs(sp_std - 1.0) < 1e-6, (
                f"Species '{sp}': expected std=1.0, got {sp_std}"
            )

    def test_zscore_handles_constant(self):
        """A feature constant within a species must not produce NaN (epsilon guard)."""
        data = pd.DataFrame(
            {
                "constant_feat": [5.0, 5.0, 5.0, 10.0, 20.0, 30.0],
                "normal_feat": [1.0, 2.0, 3.0, 1.0, 2.0, 3.0],
            }
        )
        species = ["mouse", "mouse", "mouse", "human", "human", "human"]

        result = zscore_within_species(data, species)

        # constant_feat for mouse (all 5.0) should be all 0.0, not NaN
        mouse_constant = result.iloc[[0, 1, 2]]["constant_feat"]
        assert not mouse_constant.isna().any(), (
            "Constant feature within species produced NaN values"
        )
        np.testing.assert_allclose(
            mouse_constant.values, 0.0, atol=1e-6,
            err_msg="Constant feature within species should z-score to 0",
        )

    def test_zscore_preserves_shape(self):
        """Output DataFrame has same shape, index, and columns as input."""
        rng = np.random.default_rng(0)
        data = pd.DataFrame(
            rng.standard_normal((6, 4)),
            columns=["f1", "f2", "f3", "f4"],
            index=[f"cell{i}" for i in range(6)],
        )
        species = ["mouse"] * 3 + ["human"] * 3

        result = zscore_within_species(data, species)

        assert result.shape == data.shape
        assert list(result.columns) == list(data.columns)
        assert list(result.index) == list(data.index)


# ---------------------------------------------------------------------------
# permutation_test_cluster_similarity
# ---------------------------------------------------------------------------


class TestPermutationTestClusterSimilarity:
    """Tests for permutation_test_cluster_similarity."""

    def test_permutation_test_returns_valid_pvalue(self):
        """Random data should produce a p-value strictly between 0 and 1."""
        rng = np.random.default_rng(1)
        n_cells = 40
        data = pd.DataFrame(
            rng.standard_normal((n_cells, 5)),
            columns=[f"f{i}" for i in range(5)],
        )
        clusters = np.array(["A"] * 20 + ["B"] * 20)
        species = np.array(["mouse"] * 10 + ["human"] * 10 + ["mouse"] * 10 + ["human"] * 10)

        p = permutation_test_cluster_similarity(
            data, clusters, species, n_perm=200, seed=42
        )

        assert 0.0 <= p <= 1.0, f"p-value {p} outside [0, 1]"

    def test_permutation_test_well_separated_clusters_low_pvalue(self):
        """Well-matched clusters (same centroid per cluster across species) yield low p.

        Uses 100 cells × 20 features with a large offset (+/-20) so the
        cross-species centroid distances are reliably small relative to the
        shuffled null.
        """
        rng = np.random.default_rng(7)
        n_per = 100
        n_feat = 20
        offset = 20.0
        # Cluster A: centroid at +offset in all features; Cluster B: -offset
        A_mouse = rng.standard_normal((n_per, n_feat)) + offset
        B_mouse = rng.standard_normal((n_per, n_feat)) - offset
        A_human = rng.standard_normal((n_per, n_feat)) + offset
        B_human = rng.standard_normal((n_per, n_feat)) - offset

        data = pd.DataFrame(
            np.vstack([A_mouse, B_mouse, A_human, B_human]),
            columns=[f"f{i}" for i in range(n_feat)],
        )
        clusters = np.array(["A"] * n_per + ["B"] * n_per + ["A"] * n_per + ["B"] * n_per)
        species = np.array(["mouse"] * (2 * n_per) + ["human"] * (2 * n_per))

        p = permutation_test_cluster_similarity(
            data, clusters, species, n_perm=500, seed=42
        )

        # Well-matched clusters should be more similar than shuffled → low p-value
        assert p < 0.05, f"Expected p < 0.05 for well-matched clusters, got {p:.3f}"

    def test_permutation_test_reproducible(self):
        """Same seed gives the same p-value."""
        rng = np.random.default_rng(99)
        data = pd.DataFrame(rng.standard_normal((20, 4)), columns=list("ABCD"))
        clusters = np.array(["X"] * 10 + ["Y"] * 10)
        species = np.array(["mouse"] * 5 + ["human"] * 5 + ["mouse"] * 5 + ["human"] * 5)

        p1 = permutation_test_cluster_similarity(data, clusters, species, n_perm=100, seed=0)
        p2 = permutation_test_cluster_similarity(data, clusters, species, n_perm=100, seed=0)

        assert p1 == pytest.approx(p2), "Same seed should produce the same p-value"
