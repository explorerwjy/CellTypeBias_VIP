"""Multi-modal CCKBC confidence scoring and classification."""

import pandas as pd

EVIDENCE_COLUMNS = [
    "module_a_is_sncg",         # Module A: mouse cluster is Sncg/CCKBC subclass
    "module_c_direct_cckbc",    # Module C: direct path maps CCKBC cells here
    "module_c_indirect_cckbc",  # Module C: indirect path maps CCKBC cells here
    "module_d_fast_spiking",    # Module D: ephys resembles fast-spiking CCKBC
    "marker_cck_positive",      # Existing: CCK+ marker expression
    "marker_sncg_positive",     # Existing: SNCG+ marker expression
]

_TIER_HIGH = "high-confidence CCKBC"
_TIER_TENTATIVE = "tentative CCKBC"
_TIER_ISI = "ISI VIP"


def compute_cckbc_confidence(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Compute CCKBC confidence score and tier for each cluster.

    Parameters
    ----------
    evidence_df : pd.DataFrame
        DataFrame with a cluster_id column (or index) and one or more boolean
        columns from EVIDENCE_COLUMNS.  Missing evidence columns are silently
        ignored; only the columns that are present are summed.

    Returns
    -------
    pd.DataFrame
        Copy of *evidence_df* with two additional columns:

        ``cckbc_confidence`` : int
            Number of True evidence flags present for each cluster.
        ``cckbc_tier`` : str
            "high-confidence CCKBC"  (confidence >= 3)
            "tentative CCKBC"        (confidence in [1, 2])
            "ISI VIP"                (confidence == 0)
    """
    df = evidence_df.copy()

    # Only sum columns that are actually present in the DataFrame
    available = [c for c in EVIDENCE_COLUMNS if c in df.columns]

    if available:
        df["cckbc_confidence"] = df[available].astype(bool).sum(axis=1).astype(int)
    else:
        df["cckbc_confidence"] = 0

    def _assign_tier(score: int) -> str:
        if score >= 3:
            return _TIER_HIGH
        elif score >= 1:
            return _TIER_TENTATIVE
        else:
            return _TIER_ISI

    df["cckbc_tier"] = df["cckbc_confidence"].map(_assign_tier)

    return df


def classify_clusters(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Convenience wrapper around compute_cckbc_confidence.

    Parameters
    ----------
    evidence_df : pd.DataFrame
        See :func:`compute_cckbc_confidence`.

    Returns
    -------
    pd.DataFrame
        Same as :func:`compute_cckbc_confidence`.
    """
    return compute_cckbc_confidence(evidence_df)
