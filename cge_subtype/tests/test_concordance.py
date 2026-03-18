"""Tests for concordance module (Module C)."""

import math

import numpy as np
import pandas as pd
import pytest

from cge_subtype.src.concordance import (
    build_indirect_path,
    compute_concordance,
    compute_cohens_kappa,
    concordance_at_levels,
)


# ---------------------------------------------------------------------------
# build_indirect_path
# ---------------------------------------------------------------------------

class TestBuildIndirectPath:
    """Tests for build_indirect_path."""

    def test_basic_chaining(self):
        """3 cells, 2 mouse clusters mapping to 2 human clusters."""
        patchseq_to_mouse = pd.DataFrame({
            "cell_id": ["c0", "c1", "c2"],
            "mouse_cluster": ["Sncg_1", "Sncg_2", "Sncg_1"],
        })
        mouse_to_human = pd.DataFrame({
            "mouse_cluster": ["Sncg_1", "Sncg_2"],
            "human_cluster": ["CGE_A", "CGE_B"],
        })

        result = build_indirect_path(patchseq_to_mouse, mouse_to_human)

        assert "indirect_human_cluster" in result.columns
        assert result.loc[result["cell_id"] == "c0", "indirect_human_cluster"].iloc[0] == "CGE_A"
        assert result.loc[result["cell_id"] == "c1", "indirect_human_cluster"].iloc[0] == "CGE_B"
        assert result.loc[result["cell_id"] == "c2", "indirect_human_cluster"].iloc[0] == "CGE_A"

    def test_unmapped_cluster_gives_nan(self):
        """A mouse cluster with no human match results in NaN."""
        patchseq_to_mouse = pd.DataFrame({
            "cell_id": ["c0", "c1"],
            "mouse_cluster": ["Sncg_1", "Sncg_unknown"],
        })
        mouse_to_human = pd.DataFrame({
            "mouse_cluster": ["Sncg_1"],
            "human_cluster": ["CGE_A"],
        })

        result = build_indirect_path(patchseq_to_mouse, mouse_to_human)

        assert result.loc[result["cell_id"] == "c0", "indirect_human_cluster"].iloc[0] == "CGE_A"
        assert pd.isna(result.loc[result["cell_id"] == "c1", "indirect_human_cluster"].iloc[0])

    def test_preserves_original_columns(self):
        """Output DataFrame retains all original columns."""
        patchseq_to_mouse = pd.DataFrame({
            "cell_id": ["c0"],
            "mouse_cluster": ["Sncg_1"],
            "extra_col": [42.0],
        })
        mouse_to_human = pd.DataFrame({
            "mouse_cluster": ["Sncg_1"],
            "human_cluster": ["CGE_A"],
        })

        result = build_indirect_path(patchseq_to_mouse, mouse_to_human)

        assert "extra_col" in result.columns
        assert result.loc[0, "extra_col"] == 42.0

    def test_does_not_modify_input(self):
        """Input DataFrames are not modified in place."""
        patchseq_to_mouse = pd.DataFrame({
            "cell_id": ["c0"],
            "mouse_cluster": ["Sncg_1"],
        })
        mouse_to_human = pd.DataFrame({
            "mouse_cluster": ["Sncg_1"],
            "human_cluster": ["CGE_A"],
        })
        original_cols_ps = list(patchseq_to_mouse.columns)
        original_cols_m2h = list(mouse_to_human.columns)

        build_indirect_path(patchseq_to_mouse, mouse_to_human)

        assert list(patchseq_to_mouse.columns) == original_cols_ps
        assert list(mouse_to_human.columns) == original_cols_m2h


# ---------------------------------------------------------------------------
# compute_concordance
# ---------------------------------------------------------------------------

class TestComputeConcordance:
    """Tests for compute_concordance."""

    def test_concordance_perfect(self):
        """Identical labels → concordance == 1.0."""
        labels = pd.Series(["A", "B", "C", "A"])
        assert compute_concordance(labels, labels.copy()) == pytest.approx(1.0)

    def test_concordance_none(self):
        """Completely different labels → concordance == 0.0."""
        direct = pd.Series(["A", "A", "A", "A"])
        indirect = pd.Series(["B", "B", "B", "B"])
        assert compute_concordance(direct, indirect) == pytest.approx(0.0)

    def test_concordance_partial(self):
        """2 of 4 cells agree → concordance == 0.5."""
        direct = pd.Series(["A", "B", "C", "D"])
        indirect = pd.Series(["A", "B", "X", "Y"])
        assert compute_concordance(direct, indirect) == pytest.approx(0.5)

    def test_concordance_nan_excluded(self):
        """NaN entries are excluded from the count."""
        direct = pd.Series(["A", "B", np.nan, "C"])
        indirect = pd.Series(["A", "X", "C", np.nan])
        # Valid pairs: index 0 (A==A agree) and index 1 (B!=X disagree) → 1/2 = 0.5
        # Indices 2 and 3 are excluded because one side is NaN
        assert compute_concordance(direct, indirect) == pytest.approx(0.5)

    def test_concordance_all_nan_returns_nan(self):
        """All-NaN input returns float NaN."""
        direct = pd.Series([np.nan, np.nan])
        indirect = pd.Series([np.nan, np.nan])
        result = compute_concordance(direct, indirect)
        assert math.isnan(result)


# ---------------------------------------------------------------------------
# compute_cohens_kappa
# ---------------------------------------------------------------------------

class TestComputeCohensKappa:
    """Tests for compute_cohens_kappa."""

    def test_kappa_perfect(self):
        """Identical labels → kappa == 1.0."""
        labels = pd.Series(["A", "B", "A", "B", "C"])
        assert compute_cohens_kappa(labels, labels.copy()) == pytest.approx(1.0)

    def test_kappa_nan_filtered(self):
        """NaN entries are removed before kappa computation."""
        a = pd.Series(["A", "B", np.nan, "A"])
        b = pd.Series(["A", "B", "C", np.nan])
        # Only cells 0 and 1 are valid; both agree → kappa == 1.0
        assert compute_cohens_kappa(a, b) == pytest.approx(1.0)

    def test_kappa_returns_float(self):
        """Return value is a Python float."""
        a = pd.Series(["X", "Y", "X", "Y"])
        b = pd.Series(["X", "Y", "Y", "X"])
        result = compute_cohens_kappa(a, b)
        assert isinstance(result, float)

    def test_kappa_range(self):
        """Kappa is in valid range [-1, 1]."""
        a = pd.Series(["A", "B", "A", "B", "C", "C"])
        b = pd.Series(["B", "A", "B", "A", "C", "C"])
        result = compute_cohens_kappa(a, b)
        assert -1.0 <= result <= 1.0

    def test_kappa_insufficient_data_returns_nan(self):
        """Fewer than 2 valid pairs → returns NaN."""
        a = pd.Series(["A", np.nan])
        b = pd.Series([np.nan, "B"])
        result = compute_cohens_kappa(a, b)
        assert math.isnan(result)


# ---------------------------------------------------------------------------
# concordance_at_levels
# ---------------------------------------------------------------------------

class TestConcordanceAtLevels:
    """Tests for concordance_at_levels."""

    def _make_data(self):
        """6 cells with known cluster, subclass, supercluster structure."""
        # direct path: all cells labelled with cluster IDs
        direct = pd.Series(["c1", "c1", "c2", "c2", "c3", "c3"], name="direct")
        # indirect path: agrees on clusters 1 & 2, disagrees on 3 (c3 vs c4)
        indirect = pd.Series(["c1", "c1", "c2", "c4", "c4", "c4"], name="indirect")

        # c1,c2 → subclass A | c3,c4 → subclass B
        cluster_to_subclass = {
            "c1": "SubA",
            "c2": "SubA",
            "c3": "SubB",
            "c4": "SubB",
        }
        # All clusters → same supercluster
        cluster_to_supercluster = {
            "c1": "Super1",
            "c2": "Super1",
            "c3": "Super1",
            "c4": "Super1",
        }
        return direct, indirect, cluster_to_subclass, cluster_to_supercluster

    def test_returns_six_keys(self):
        """Result dict has exactly the 6 expected keys."""
        d, i, c2s, c2sc = self._make_data()
        result = concordance_at_levels(d, i, c2s, c2sc)

        expected_keys = {
            "cluster_concordance",
            "cluster_kappa",
            "subclass_concordance",
            "subclass_kappa",
            "supercluster_concordance",
            "supercluster_kappa",
        }
        assert set(result.keys()) == expected_keys

    def test_cluster_concordance_value(self):
        """Cluster-level concordance: 3 of 6 cells agree → 0.5."""
        d, i, c2s, c2sc = self._make_data()
        # direct:   c1 c1 c2 c2 c3 c3
        # indirect: c1 c1 c2 c4 c4 c4
        # agree:    Y  Y  Y  N  N  N  → 3/6 = 0.5
        result = concordance_at_levels(d, i, c2s, c2sc)
        assert result["cluster_concordance"] == pytest.approx(0.5)

    def test_subclass_concordance_value(self):
        """Subclass-level concordance: 4 of 6 cells agree."""
        d, i, c2s, c2sc = self._make_data()
        # direct_sub:   SubA SubA SubA SubA SubB SubB
        # indirect_sub: SubA SubA SubA SubB SubB SubB
        # agree:         Y    Y    Y    N    Y    Y  → 5/6
        result = concordance_at_levels(d, i, c2s, c2sc)
        assert result["subclass_concordance"] == pytest.approx(5 / 6)

    def test_supercluster_concordance_perfect(self):
        """All clusters map to the same supercluster → concordance == 1.0."""
        d, i, c2s, c2sc = self._make_data()
        result = concordance_at_levels(d, i, c2s, c2sc)
        assert result["supercluster_concordance"] == pytest.approx(1.0)

    def test_all_values_are_floats(self):
        """All returned values are Python floats (or nan)."""
        d, i, c2s, c2sc = self._make_data()
        result = concordance_at_levels(d, i, c2s, c2sc)
        for key, val in result.items():
            assert isinstance(val, float), f"{key} is not float: {type(val)}"

    def test_kappa_keys_present(self):
        """Kappa values are present and within [-1, 1] (or NaN)."""
        d, i, c2s, c2sc = self._make_data()
        result = concordance_at_levels(d, i, c2s, c2sc)
        for key in ["cluster_kappa", "subclass_kappa", "supercluster_kappa"]:
            v = result[key]
            if not math.isnan(v):
                assert -1.0 <= v <= 1.0, f"{key}={v} out of range"
