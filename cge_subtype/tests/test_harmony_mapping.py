"""Tests for harmony_mapping module (Module B)."""

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse

from cge_subtype.src.harmony_mapping import (
    subsample_reference,
    evaluate_mapping,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_adata(n_cells: int, n_genes: int, cluster_labels, seed: int = 0) -> ad.AnnData:
    """Create a minimal AnnData with random dense expression and cluster labels."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame(
        {"cluster": cluster_labels},
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"gene{j}" for j in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# Tests for subsample_reference
# ---------------------------------------------------------------------------

class TestSubsampleReference:
    """Tests for subsample_reference."""

    def test_subsample_caps_cluster_size(self):
        """Large cluster (500 cells) is capped to max_per_cluster."""
        # Cluster A: 500 cells (large), Cluster B: 50 cells (small)
        labels = ["A"] * 500 + ["B"] * 50
        adata = _make_adata(550, 20, labels, seed=1)

        result = subsample_reference(adata, cluster_col="cluster", max_per_cluster=200, seed=42)

        counts = result.obs["cluster"].value_counts()
        assert counts["A"] == 200, f"Expected 200 cells for cluster A, got {counts['A']}"

    def test_subsample_preserves_small_clusters(self):
        """Small clusters (< max_per_cluster) are retained in full."""
        labels = ["A"] * 500 + ["B"] * 50
        adata = _make_adata(550, 20, labels, seed=2)

        result = subsample_reference(adata, cluster_col="cluster", max_per_cluster=200, seed=42)

        counts = result.obs["cluster"].value_counts()
        assert counts["B"] == 50, f"Expected all 50 cells for cluster B, got {counts['B']}"

    def test_subsample_total_size(self):
        """Total output size equals sum of capped cluster sizes."""
        labels = ["A"] * 500 + ["B"] * 50 + ["C"] * 300
        adata = _make_adata(850, 20, labels, seed=3)

        result = subsample_reference(adata, cluster_col="cluster", max_per_cluster=100, seed=42)

        # A: 100, B: 50 (unchanged), C: 100 → total 250
        assert len(result) == 250, f"Expected 250 total cells, got {len(result)}"

    def test_subsample_returns_copy(self):
        """Result is a copy; modifying it does not affect original."""
        labels = ["A"] * 50
        adata = _make_adata(50, 10, labels, seed=4)

        result = subsample_reference(adata, cluster_col="cluster", max_per_cluster=200, seed=42)

        # Modify result obs
        result.obs["cluster"] = "Z"
        # Original should be unchanged
        assert (adata.obs["cluster"] == "A").all()

    def test_subsample_reproducible_with_seed(self):
        """Two calls with the same seed return the same cells."""
        labels = ["A"] * 500 + ["B"] * 500
        adata = _make_adata(1000, 10, labels, seed=5)

        r1 = subsample_reference(adata, cluster_col="cluster", max_per_cluster=100, seed=99)
        r2 = subsample_reference(adata, cluster_col="cluster", max_per_cluster=100, seed=99)

        assert list(r1.obs_names) == list(r2.obs_names)

    def test_subsample_different_seeds_differ(self):
        """Different seeds produce different subsamples for large clusters."""
        labels = ["A"] * 500
        adata = _make_adata(500, 10, labels, seed=6)

        r1 = subsample_reference(adata, cluster_col="cluster", max_per_cluster=50, seed=1)
        r2 = subsample_reference(adata, cluster_col="cluster", max_per_cluster=50, seed=2)

        # It's astronomically unlikely that two random seeds pick the same 50 cells
        assert set(r1.obs_names) != set(r2.obs_names)


# ---------------------------------------------------------------------------
# Tests for evaluate_mapping
# ---------------------------------------------------------------------------

class TestEvaluateMapping:
    """Tests for evaluate_mapping."""

    def _make_perfect_predictions(self):
        """All predictions match ground truth."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1", "c2", "c3", "c4"],
            "predicted_cluster": ["VIP Htr3a", "VIP Htr3a", "SST Chodl", "SST Chodl", "Pvalb Tac1"],
        })
        metadata = pd.DataFrame(
            {
                "ground_truth": ["VIP Htr3a", "VIP Htr3a", "SST Chodl", "SST Chodl", "Pvalb Tac1"],
            },
            index=["c0", "c1", "c2", "c3", "c4"],
        )
        return predictions, metadata

    def _make_partial_predictions(self):
        """3 of 5 predictions correct at subclass level."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1", "c2", "c3", "c4"],
            # c0: VIP vs VIP ✓, c1: VIP vs SST ✗, c2: SST vs SST ✓,
            # c3: Pvalb vs Pvalb ✓, c4: VIP vs Astro ✗
            "predicted_cluster": ["VIP Htr3a", "SST Chodl", "SST Chodl", "Pvalb Tac1", "Astro Fgfr3"],
        })
        metadata = pd.DataFrame(
            {
                "ground_truth": ["VIP Htr3a", "VIP Htr3a", "SST Chodl", "Pvalb Tac1", "VIP Htr3a"],
            },
            index=["c0", "c1", "c2", "c3", "c4"],
        )
        return predictions, metadata

    def test_evaluate_perfect_accuracy(self):
        """Perfect predictions yield overall_accuracy == 1.0."""
        preds, meta = self._make_perfect_predictions()
        result = evaluate_mapping(preds, ground_truth_col="ground_truth", metadata=meta)

        assert result["overall_accuracy"] == pytest.approx(1.0)
        assert result["n_cells"] == 5

    def test_evaluate_known_accuracy(self):
        """Provide known predictions and verify accuracy == 3/5 == 0.6."""
        preds, meta = self._make_partial_predictions()
        result = evaluate_mapping(preds, ground_truth_col="ground_truth", metadata=meta)

        # c0: VIP == VIP ✓
        # c1: VIP != SST ✗
        # c2: SST == SST ✓
        # c3: Pvalb == Pvalb ✓
        # c4: VIP != Astro ✗
        assert result["overall_accuracy"] == pytest.approx(3 / 5)
        assert result["n_cells"] == 5

    def test_evaluate_non_neuronal_rate(self):
        """Non-neuronal predictions are counted correctly."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1", "c2", "c3", "c4"],
            # Astro, Oligo are non-neuronal; VIP, SST, Pvalb are neuronal
            "predicted_cluster": ["Astro Fgfr3", "Oligo Opalin", "VIP Htr3a", "SST Chodl", "Pvalb Tac1"],
        })
        metadata = pd.DataFrame(
            {"ground_truth": ["Astro", "Oligo", "VIP", "SST", "Pvalb"]},
            index=["c0", "c1", "c2", "c3", "c4"],
        )
        result = evaluate_mapping(predictions, ground_truth_col="ground_truth", metadata=metadata)

        # 2 out of 5 predicted labels are non-neuronal
        assert result["non_neuronal_rate"] == pytest.approx(2 / 5)

    def test_evaluate_nn_prefix(self):
        """Labels starting with 'NN ' are counted as non-neuronal."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1"],
            "predicted_cluster": ["NN Astro", "VIP Htr3a"],
        })
        metadata = pd.DataFrame(
            {"ground_truth": ["NN Astro", "VIP Htr3a"]},
            index=["c0", "c1"],
        )
        result = evaluate_mapping(predictions, ground_truth_col="ground_truth", metadata=metadata)

        assert result["non_neuronal_rate"] == pytest.approx(0.5)

    def test_evaluate_returns_confusion_matrix(self):
        """Return value includes a confusion matrix DataFrame."""
        preds, meta = self._make_partial_predictions()
        result = evaluate_mapping(preds, ground_truth_col="ground_truth", metadata=meta)

        cm = result["confusion_matrix"]
        assert isinstance(cm, pd.DataFrame)
        # Rows/cols should be subclass labels (first words)
        assert "VIP" in cm.index or "VIP" in cm.columns

    def test_evaluate_zero_accuracy(self):
        """All wrong predictions yield overall_accuracy == 0.0."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1"],
            "predicted_cluster": ["SST Chodl", "Pvalb Tac1"],
        })
        metadata = pd.DataFrame(
            {"ground_truth": ["VIP Htr3a", "VIP Htr3a"]},
            index=["c0", "c1"],
        )
        result = evaluate_mapping(predictions, ground_truth_col="ground_truth", metadata=metadata)

        assert result["overall_accuracy"] == pytest.approx(0.0)

    def test_evaluate_metadata_with_cell_id_column(self):
        """metadata with a cell_id column (not index) is handled correctly."""
        predictions = pd.DataFrame({
            "cell_id": ["c0", "c1"],
            "predicted_cluster": ["VIP Htr3a", "SST Chodl"],
        })
        # metadata has cell_id as a regular column
        metadata = pd.DataFrame({
            "cell_id": ["c0", "c1"],
            "ground_truth": ["VIP Htr3a", "SST Chodl"],
        })
        result = evaluate_mapping(predictions, ground_truth_col="ground_truth", metadata=metadata)

        assert result["overall_accuracy"] == pytest.approx(1.0)
        assert result["n_cells"] == 2
