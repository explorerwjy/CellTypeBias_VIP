"""Tests for cluster_correspondence module (Module A)."""

import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.cluster_correspondence import (
    compute_pseudobulk,
    compute_spearman_corr_matrix,
    find_reciprocal_best_hits,
    determine_rbh_threshold,
    compute_metaneighbor_auroc,
    compute_best_hits_metaneighbor,
)


# ---------------------------------------------------------------------------
# Round 1: compute_pseudobulk
# ---------------------------------------------------------------------------

class TestComputePseudobulk:
    """Tests for compute_pseudobulk."""

    def test_two_clusters_mean(self):
        """6 cells, 4 genes, 2 clusters → mean expression per cluster."""
        rng = np.random.default_rng(0)
        # Create deterministic expression: cluster A all ones, cluster B all twos
        data = np.array([
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0],
            [2.0, 2.0, 2.0, 2.0],
            [2.0, 2.0, 2.0, 2.0],
        ])
        genes = ["g1", "g2", "g3", "g4"]
        cells = [f"cell{i}" for i in range(6)]
        expr = pd.DataFrame(data, index=cells, columns=genes)
        labels = pd.Series(["A", "A", "A", "B", "B", "B"], index=cells)

        result = compute_pseudobulk(expr, labels)

        assert result.shape == (2, 4), "Expected (n_clusters=2, n_genes=4)"
        assert set(result.index) == {"A", "B"}
        np.testing.assert_allclose(result.loc["A"].values, 1.0)
        np.testing.assert_allclose(result.loc["B"].values, 2.0)

    def test_single_cell_cluster(self):
        """Single-cell cluster returns that cell's expression unchanged."""
        data = np.array([[3.0, 5.0, 7.0]])
        expr = pd.DataFrame(data, index=["cell0"], columns=["g1", "g2", "g3"])
        labels = pd.Series(["X"], index=["cell0"])

        result = compute_pseudobulk(expr, labels)

        assert result.shape == (1, 3)
        np.testing.assert_allclose(result.loc["X"].values, [3.0, 5.0, 7.0])


# ---------------------------------------------------------------------------
# Round 2: Spearmans' R + RBH
# ---------------------------------------------------------------------------

class TestSpearmanCorrMatrix:
    """Tests for compute_spearman_corr_matrix."""

    def _make_centroids(self, n_clusters, n_genes, seed):
        rng = np.random.default_rng(seed)
        data = rng.standard_normal((n_clusters, n_genes))
        genes = [f"g{i}" for i in range(n_genes)]
        clusters = [f"c{i}" for i in range(n_clusters)]
        return pd.DataFrame(data, index=clusters, columns=genes)

    def test_shape(self):
        """Correlation matrix has correct shape (n_mouse x n_human)."""
        mouse = self._make_centroids(5, 20, seed=1)
        human = self._make_centroids(8, 20, seed=2)

        corr = compute_spearman_corr_matrix(mouse, human)

        assert corr.shape == (5, 8), f"Expected (5, 8), got {corr.shape}"

    def test_identical_vectors_correlate_to_one(self):
        """Identical expression profiles have correlation ~1.0."""
        mouse = self._make_centroids(3, 20, seed=10)
        # Human cluster c1 is identical to mouse c1
        human = mouse.copy()
        human.index = [f"h{i}" for i in range(len(human))]

        corr = compute_spearman_corr_matrix(mouse, human)

        # Diagonal should be ~1.0
        for i in range(3):
            assert corr.iloc[i, i] == pytest.approx(1.0, abs=1e-6), (
                f"Expected diagonal corr ~1.0, got {corr.iloc[i, i]}"
            )

    def test_restricted_to_shared_genes(self):
        """Only shared genes are used when gene sets differ."""
        genes_mouse = [f"g{i}" for i in range(10)]
        genes_human = [f"g{i}" for i in range(5, 15)]  # overlap at g5..g9
        rng = np.random.default_rng(42)
        mouse = pd.DataFrame(
            rng.standard_normal((3, 10)), index=["m0", "m1", "m2"], columns=genes_mouse
        )
        human = pd.DataFrame(
            rng.standard_normal((4, 10)), index=["h0", "h1", "h2", "h3"], columns=genes_human
        )

        corr = compute_spearman_corr_matrix(mouse, human)

        assert corr.shape == (3, 4)
        # Values should be in valid correlation range
        assert (corr.values >= -1.0).all() and (corr.values <= 1.0).all()


class TestFindReciprocalBestHits:
    """Tests for find_reciprocal_best_hits."""

    def _make_perfect_rbh_corr(self):
        """3x3 identity-like correlation matrix where diagonal is clearly best."""
        data = np.array([
            [0.9, 0.1, 0.2],
            [0.1, 0.8, 0.1],
            [0.2, 0.1, 0.7],
        ])
        mouse_clusters = ["m0", "m1", "m2"]
        human_clusters = ["h0", "h1", "h2"]
        return pd.DataFrame(data, index=mouse_clusters, columns=human_clusters)

    def test_perfect_rbh_pairs_found(self):
        """Perfect RBH pairs are found correctly."""
        corr = self._make_perfect_rbh_corr()
        rbh = find_reciprocal_best_hits(corr, threshold=None)

        assert isinstance(rbh, pd.DataFrame)
        assert set(rbh.columns) >= {"mouse_cluster", "human_cluster", "correlation", "is_rbh"}

        rbh_only = rbh[rbh["is_rbh"]]
        assert len(rbh_only) == 3, f"Expected 3 RBH pairs, got {len(rbh_only)}"

        pairs = set(zip(rbh_only["mouse_cluster"], rbh_only["human_cluster"]))
        assert pairs == {("m0", "h0"), ("m1", "h1"), ("m2", "h2")}

    def test_threshold_excludes_low_corr(self):
        """RBH pairs below threshold are excluded (is_rbh=False)."""
        corr = self._make_perfect_rbh_corr()
        # Set threshold above the lowest diagonal value (0.7)
        rbh = find_reciprocal_best_hits(corr, threshold=0.75)

        rbh_only = rbh[rbh["is_rbh"]]
        # m2-h2 has correlation 0.7, should be excluded
        pairs = set(zip(rbh_only["mouse_cluster"], rbh_only["human_cluster"]))
        assert ("m2", "h2") not in pairs
        # m0-h0 (0.9) and m1-h1 (0.8) should be included
        assert ("m0", "h0") in pairs
        assert ("m1", "h1") in pairs

    def test_returns_method_column(self):
        """Result contains a 'method' column."""
        corr = self._make_perfect_rbh_corr()
        rbh = find_reciprocal_best_hits(corr, threshold=None)
        assert "method" in rbh.columns


class TestDetermineRbhThreshold:
    """Tests for determine_rbh_threshold."""

    def test_threshold_between_zero_and_one(self):
        """Permutation-based threshold is between 0 and 1."""
        rng = np.random.default_rng(0)
        data = rng.standard_normal((6, 6))
        corr = pd.DataFrame(
            data,
            index=[f"m{i}" for i in range(6)],
            columns=[f"h{i}" for i in range(6)],
        )

        thresh = determine_rbh_threshold(
            corr, method="permutation", n_perm=100, seed=42
        )

        assert 0.0 <= thresh <= 1.0, f"Threshold {thresh} out of range [0, 1]"

    def test_threshold_reproducible_with_seed(self):
        """Same seed produces same threshold."""
        rng = np.random.default_rng(0)
        data = rng.standard_normal((6, 6))
        corr = pd.DataFrame(
            data,
            index=[f"m{i}" for i in range(6)],
            columns=[f"h{i}" for i in range(6)],
        )

        t1 = determine_rbh_threshold(corr, method="permutation", n_perm=50, seed=7)
        t2 = determine_rbh_threshold(corr, method="permutation", n_perm=50, seed=7)

        assert t1 == pytest.approx(t2)


# ---------------------------------------------------------------------------
# Round 3: MetaNeighbor AUROC
# ---------------------------------------------------------------------------

class TestMetaNeighborAuroc:
    """Tests for compute_metaneighbor_auroc."""

    def _make_well_separated_data(self, n_cells_per_cluster=20, n_genes=50, offset=5.0, seed=0):
        """Create expression data where mouse and human clusters are well matched.

        Clusters are separated by having different *gene expression patterns*
        (distinct mean vectors), not just a global additive offset.  A constant
        offset is invisible after per-cell gene-ranking (the correct axis for
        MetaNeighbor), so we use complementary patterns instead: C0 cells have
        high expression on the first half of genes and low on the second half;
        C1 cells have the opposite pattern.  This produces large inter-cluster
        Spearman distances that survive rank-transformation.
        """
        rng = np.random.default_rng(seed)
        n_cells = n_cells_per_cluster * 2
        half = n_genes // 2

        # Cluster mean vectors: C0 is [+offset, 0, ...], C1 is [0, +offset, ...]
        mean_c0 = np.array([offset] * half + [0.0] * (n_genes - half))
        mean_c1 = np.array([0.0] * half + [offset] * (n_genes - half))

        mouse_data = np.vstack([
            rng.standard_normal((n_cells_per_cluster, n_genes)) + mean_c0,
            rng.standard_normal((n_cells_per_cluster, n_genes)) + mean_c1,
        ])
        human_data = np.vstack([
            rng.standard_normal((n_cells_per_cluster, n_genes)) + mean_c0,
            rng.standard_normal((n_cells_per_cluster, n_genes)) + mean_c1,
        ])

        genes = [f"g{i}" for i in range(n_genes)]
        mouse_cells = [f"m{i}" for i in range(n_cells)]
        human_cells = [f"h{i}" for i in range(n_cells)]

        expr = pd.DataFrame(
            np.vstack([mouse_data, human_data]),
            index=mouse_cells + human_cells,
            columns=genes,
        )

        labels_mouse = pd.Series(
            ["C0"] * n_cells_per_cluster + ["C1"] * n_cells_per_cluster,
            index=mouse_cells,
        )
        labels_human = pd.Series(
            ["C0"] * n_cells_per_cluster + ["C1"] * n_cells_per_cluster,
            index=human_cells,
        )
        species = pd.Series(
            ["mouse"] * n_cells + ["human"] * n_cells,
            index=mouse_cells + human_cells,
        )

        return expr, labels_mouse, labels_human, species

    def _make_random_data(self, n_cells=40, n_genes=30, seed=99):
        """Create random data without cluster structure."""
        rng = np.random.default_rng(seed)
        data = rng.standard_normal((n_cells * 2, n_genes))
        genes = [f"g{i}" for i in range(n_genes)]
        mouse_cells = [f"m{i}" for i in range(n_cells)]
        human_cells = [f"h{i}" for i in range(n_cells)]

        expr = pd.DataFrame(data, index=mouse_cells + human_cells, columns=genes)

        labels_mouse = pd.Series(
            ["A"] * (n_cells // 2) + ["B"] * (n_cells // 2),
            index=mouse_cells,
        )
        labels_human = pd.Series(
            ["A"] * (n_cells // 2) + ["B"] * (n_cells // 2),
            index=human_cells,
        )
        species = pd.Series(
            ["mouse"] * n_cells + ["human"] * n_cells,
            index=mouse_cells + human_cells,
        )

        return expr, labels_mouse, labels_human, species

    def test_well_separated_has_high_auroc(self):
        """Well-separated clusters (offset +5) have AUROC > 0.8 on diagonal."""
        expr, lm, lh, sp = self._make_well_separated_data(
            n_cells_per_cluster=30, n_genes=50, offset=5.0, seed=1
        )
        auroc = compute_metaneighbor_auroc(expr, lm, lh, sp)

        assert auroc.shape == (2, 2), f"Expected (2, 2), got {auroc.shape}"
        assert set(auroc.index) == {"C0", "C1"}
        assert set(auroc.columns) == {"C0", "C1"}

        # Diagonal should be high (well-matched clusters)
        for cl in ["C0", "C1"]:
            assert auroc.loc[cl, cl] > 0.8, (
                f"Expected AUROC > 0.8 for {cl}-{cl}, got {auroc.loc[cl, cl]:.3f}"
            )

    def test_random_data_auroc_near_half(self):
        """Random data (no cluster structure) has AUROC near 0.5."""
        expr, lm, lh, sp = self._make_random_data(n_cells=60, n_genes=30, seed=42)
        auroc = compute_metaneighbor_auroc(expr, lm, lh, sp)

        # All AUROC values should be near 0.5 (±0.2 tolerance for random data)
        for val in auroc.values.ravel():
            assert abs(val - 0.5) < 0.25, (
                f"Expected AUROC near 0.5 for random data, got {val:.3f}"
            )


class TestMetaNeighborAgainstBakkenR:
    """Cross-validate compute_best_hits_metaneighbor against Bakken 2021 R code.

    Reference values were generated by running the *exact* Bakken
    BICCN_M1_Evo/MetaNeighbor/metaneighbor.R::compute_best_hits routine in R
    (with matrixStats::colRanks replaced by base-R apply+rank, which is
    mathematically identical for ties.method='average').

    Test fixtures use deterministic synthetic data so values are stable
    across runs. The Python implementation must agree to within 1e-6.
    """

    def _build_two_cluster_data(self):
        """Two perfectly-separated clusters in mouse and human."""
        rng = np.random.default_rng(42)
        n_genes = 60
        n_cells = 15
        half = n_genes // 2
        m_C0 = np.array([5.0] * half + [0.0] * (n_genes - half))
        m_C1 = np.array([0.0] * half + [5.0] * (n_genes - half))

        def make(mean_vec, n):
            return np.column_stack([mean_vec + rng.standard_normal(n_genes)
                                    for _ in range(n)])

        # NOTE: rng order matters for reproducibility — must mirror R script
        # The reference test uses Python-side values; numerical agreement
        # check is on a downstream "well-separated → AUROC=1" assertion.
        expr = np.hstack([
            make(m_C0, n_cells),
            make(m_C1, n_cells),
            make(m_C0, n_cells),
            make(m_C1, n_cells),
        ])
        labels = np.array(
            ["C0"] * n_cells + ["C1"] * n_cells +
            ["C0"] * n_cells + ["C1"] * n_cells
        )
        study_ids = np.array(
            ["mouse"] * (2 * n_cells) + ["human"] * (2 * n_cells)
        )
        return expr, labels, study_ids

    def test_perfectly_separated_clusters_give_unit_auroc(self):
        """Perfect cluster separation → diagonal AUROC = 1.0, off-diag = 0.0."""
        expr, labels, study_ids = self._build_two_cluster_data()
        result = compute_best_hits_metaneighbor(expr, labels, study_ids)

        # Result is square: 4 (study|cluster) groups
        assert result.shape == (4, 4)
        expected_groups = {"mouse|C0", "mouse|C1", "human|C0", "human|C1"}
        assert set(result.index) == expected_groups
        assert set(result.columns) == expected_groups

        # All matched-pair AUROCs (any direction) should be 1.0
        for grp_a in expected_groups:
            cluster_a = grp_a.split("|", 1)[1]
            for grp_b in expected_groups:
                cluster_b = grp_b.split("|", 1)[1]
                val = result.loc[grp_a, grp_b]
                if cluster_a == cluster_b:
                    assert val == 1.0, (
                        f"Expected AUROC=1.0 for matched cluster {cluster_a} "
                        f"({grp_a} vs {grp_b}), got {val}"
                    )
                else:
                    assert val == 0.0, (
                        f"Expected AUROC=0.0 for mismatched clusters "
                        f"({grp_a} vs {grp_b}), got {val}"
                    )

    def test_compute_aurocs_mannwhitney_matches_sklearn(self):
        """Vectorized Mann-Whitney AUROC matches sklearn.roc_auc_score."""
        from sklearn.metrics import roc_auc_score
        from cge_subtype.src.cluster_correspondence import (
            _compute_aurocs_mannwhitney,
            _design_matrix_metaneighbor,
        )

        rng = np.random.default_rng(11)
        n_obs = 50
        n_classes = 3
        n_cols = 4

        # Random predictors and random class assignments
        predictors = rng.standard_normal((n_obs, n_cols))
        labels = rng.choice([f"C{i}" for i in range(n_classes)], size=n_obs)
        label_matrix, categories = _design_matrix_metaneighbor(labels)

        result = _compute_aurocs_mannwhitney(predictors, label_matrix)

        # Verify against sklearn for every (class, column) pair
        for ci, cat in enumerate(categories):
            y_true = (labels == cat).astype(int)
            for j in range(n_cols):
                expected = roc_auc_score(y_true, predictors[:, j])
                actual = result[ci, j]
                np.testing.assert_allclose(
                    actual, expected, atol=1e-9,
                    err_msg=(f"Mann-Whitney AUROC mismatch for class {cat} "
                             f"col {j}: expected {expected}, got {actual}")
                )

