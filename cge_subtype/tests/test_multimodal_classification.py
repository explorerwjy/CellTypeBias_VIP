"""Tests for multimodal_classification module (Module E)."""

import pandas as pd
import pytest

from cge_subtype.src.multimodal_classification import (
    EVIDENCE_COLUMNS,
    compute_cckbc_confidence,
    classify_clusters,
)


def _make_row(values: list[bool], cluster_id: int = 1) -> pd.DataFrame:
    """Build a single-row evidence DataFrame from a list of booleans."""
    assert len(values) == len(EVIDENCE_COLUMNS)
    row = dict(zip(EVIDENCE_COLUMNS, values))
    row["cluster_id"] = cluster_id
    return pd.DataFrame([row])


# ---------------------------------------------------------------------------
# test_confidence_all_evidence
# ---------------------------------------------------------------------------

class TestConfidenceAllEvidence:
    """All six evidence flags True -> confidence 6, tier high-confidence CCKBC."""

    def test_score_is_six(self):
        df = _make_row([True] * 6)
        result = compute_cckbc_confidence(df)
        assert result["cckbc_confidence"].iloc[0] == 6

    def test_tier_is_high_confidence(self):
        df = _make_row([True] * 6)
        result = compute_cckbc_confidence(df)
        assert result["cckbc_tier"].iloc[0] == "high-confidence CCKBC"


# ---------------------------------------------------------------------------
# test_confidence_no_evidence
# ---------------------------------------------------------------------------

class TestConfidenceNoEvidence:
    """All six evidence flags False -> confidence 0, tier ISI VIP."""

    def test_score_is_zero(self):
        df = _make_row([False] * 6)
        result = compute_cckbc_confidence(df)
        assert result["cckbc_confidence"].iloc[0] == 0

    def test_tier_is_isi_vip(self):
        df = _make_row([False] * 6)
        result = compute_cckbc_confidence(df)
        assert result["cckbc_tier"].iloc[0] == "ISI VIP"


# ---------------------------------------------------------------------------
# test_confidence_partial
# ---------------------------------------------------------------------------

class TestConfidencePartial:
    """Exactly 2 True flags -> confidence 2, tier tentative CCKBC."""

    def test_score_is_two(self):
        df = _make_row([True, True, False, False, False, False])
        result = compute_cckbc_confidence(df)
        assert result["cckbc_confidence"].iloc[0] == 2

    def test_tier_is_tentative(self):
        df = _make_row([True, True, False, False, False, False])
        result = compute_cckbc_confidence(df)
        assert result["cckbc_tier"].iloc[0] == "tentative CCKBC"


# ---------------------------------------------------------------------------
# test_confidence_threshold_boundary
# ---------------------------------------------------------------------------

class TestConfidenceThresholdBoundary:
    """Exactly 3 True flags is the lower boundary for high-confidence."""

    def test_score_is_three(self):
        df = _make_row([True, True, True, False, False, False])
        result = compute_cckbc_confidence(df)
        assert result["cckbc_confidence"].iloc[0] == 3

    def test_tier_is_high_confidence(self):
        df = _make_row([True, True, True, False, False, False])
        result = compute_cckbc_confidence(df)
        assert result["cckbc_tier"].iloc[0] == "high-confidence CCKBC"

    def test_score_two_is_not_high_confidence(self):
        """Score of 2 must NOT reach high-confidence tier."""
        df = _make_row([True, True, False, False, False, False])
        result = compute_cckbc_confidence(df)
        assert result["cckbc_tier"].iloc[0] != "high-confidence CCKBC"


# ---------------------------------------------------------------------------
# test_classify_multiple_clusters
# ---------------------------------------------------------------------------

class TestClassifyMultipleClusters:
    """Three clusters with varying evidence receive correct tiers."""

    def _make_multi(self):
        rows = [
            # cluster 1: 5 True -> high-confidence
            dict(zip(EVIDENCE_COLUMNS, [True, True, True, True, True, False]),
                 cluster_id=1),
            # cluster 2: 1 True -> tentative
            dict(zip(EVIDENCE_COLUMNS, [False, False, False, False, True, False]),
                 cluster_id=2),
            # cluster 3: 0 True -> ISI VIP
            dict(zip(EVIDENCE_COLUMNS, [False] * 6),
                 cluster_id=3),
        ]
        return pd.DataFrame(rows)

    def test_correct_tiers(self):
        df = self._make_multi()
        result = classify_clusters(df)
        tiers = result.set_index("cluster_id")["cckbc_tier"]
        assert tiers[1] == "high-confidence CCKBC"
        assert tiers[2] == "tentative CCKBC"
        assert tiers[3] == "ISI VIP"

    def test_output_has_three_rows(self):
        df = self._make_multi()
        result = classify_clusters(df)
        assert len(result) == 3

    def test_confidence_values(self):
        df = self._make_multi()
        result = classify_clusters(df)
        confs = result.set_index("cluster_id")["cckbc_confidence"]
        assert confs[1] == 5
        assert confs[2] == 1
        assert confs[3] == 0


# ---------------------------------------------------------------------------
# test_handles_missing_columns
# ---------------------------------------------------------------------------

class TestHandlesMissingColumns:
    """Only 3 of 6 evidence columns present -> function still works."""

    def _make_partial_df(self):
        """DataFrame with only the first three EVIDENCE_COLUMNS."""
        partial_cols = EVIDENCE_COLUMNS[:3]
        rows = [
            dict(zip(partial_cols, [True, True, True]), cluster_id=10),
            dict(zip(partial_cols, [True, False, False]), cluster_id=11),
            dict(zip(partial_cols, [False, False, False]), cluster_id=12),
        ]
        return pd.DataFrame(rows)

    def test_no_exception(self):
        """Should not raise even with only 3 columns available."""
        df = self._make_partial_df()
        result = compute_cckbc_confidence(df)
        assert "cckbc_confidence" in result.columns
        assert "cckbc_tier" in result.columns

    def test_scores_sum_available_columns_only(self):
        df = self._make_partial_df()
        result = compute_cckbc_confidence(df).set_index("cluster_id")
        # cluster 10: all 3 available True -> confidence 3 -> high-confidence
        assert result.loc[10, "cckbc_confidence"] == 3
        assert result.loc[10, "cckbc_tier"] == "high-confidence CCKBC"
        # cluster 11: 1 True -> tentative
        assert result.loc[11, "cckbc_confidence"] == 1
        assert result.loc[11, "cckbc_tier"] == "tentative CCKBC"
        # cluster 12: 0 -> ISI VIP
        assert result.loc[12, "cckbc_confidence"] == 0
        assert result.loc[12, "cckbc_tier"] == "ISI VIP"

    def test_empty_evidence_df_all_zero(self):
        """DataFrame with no evidence columns at all -> confidence 0 for all."""
        df = pd.DataFrame({"cluster_id": [1, 2]})
        result = compute_cckbc_confidence(df)
        assert (result["cckbc_confidence"] == 0).all()
        assert (result["cckbc_tier"] == "ISI VIP").all()
