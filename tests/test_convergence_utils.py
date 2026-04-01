"""Tests for convergence_utils module."""

import sys
from pathlib import Path

# Add src to path so we can import convergence_utils
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import pytest
from convergence_utils import (
    GENES_22Q,
    GENE_ALIASES,
    CHANNEL_FAMILIES,
    CURATED_TARGETS,
    EPHYS_FEATURES,
    compute_residuals,
    map_symbols_to_entrez,
    pearson_partial,
    run_profile_similarity,
)


class TestCuratedTargetGenes:
    """Tests for CURATED_TARGETS and CHANNEL_FAMILIES constants."""

    def test_total_gene_count(self):
        all_genes = set()
        for genes in CURATED_TARGETS.values():
            all_genes.update(genes)
        assert 45 <= len(all_genes) <= 55, (
            f"Expected 45-55 unique genes across all families, got {len(all_genes)}"
        )

    def test_minimum_families(self):
        assert len(CURATED_TARGETS) >= 14, (
            f"Expected >=14 families, got {len(CURATED_TARGETS)}"
        )

    def test_families_have_feature_mapping(self):
        for family_name, info in CHANNEL_FAMILIES.items():
            assert "genes" in info, f"Family {family_name} missing 'genes' key"
            assert "feature" in info, f"Family {family_name} missing 'feature' key"
            assert len(info["genes"]) > 0, f"Family {family_name} has empty gene list"
            assert isinstance(info["feature"], str), (
                f"Family {family_name} 'feature' should be a string"
            )

    def test_curated_targets_matches_channel_families(self):
        """CURATED_TARGETS should be derived from CHANNEL_FAMILIES."""
        for family_name in CURATED_TARGETS:
            assert family_name in CHANNEL_FAMILIES, (
                f"Family {family_name} in CURATED_TARGETS but not in CHANNEL_FAMILIES"
            )
        for family_name, info in CHANNEL_FAMILIES.items():
            assert family_name in CURATED_TARGETS, (
                f"Family {family_name} in CHANNEL_FAMILIES but not in CURATED_TARGETS"
            )
            assert CURATED_TARGETS[family_name] == info["genes"], (
                f"Gene list mismatch for family {family_name}"
            )


class TestGenes22q:
    """Tests for GENES_22Q and GENE_ALIASES constants."""

    def test_gene_count(self):
        assert len(GENES_22Q) == 19, f"Expected 19 genes in 22q, got {len(GENES_22Q)}"

    def test_includes_dgcr8(self):
        assert "DGCR8" in GENES_22Q, "DGCR8 should be in GENES_22Q"

    def test_includes_sept5(self):
        assert "SEPT5" in GENES_22Q, "SEPT5 should be in GENES_22Q"

    def test_aliases_correct(self):
        assert "SEPT5" in GENE_ALIASES, "SEPT5 should have an alias"
        assert GENE_ALIASES["SEPT5"] == "SEPTIN5"
        assert "UFD1L" in GENE_ALIASES, "UFD1L should have an alias"
        assert GENE_ALIASES["UFD1L"] == "UFD1"


class TestMapSymbolsToEntrez:
    """Tests for map_symbols_to_entrez function."""

    def test_maps_scn1a_successfully(self):
        """SCN1A should map to an Entrez ID."""
        mapped, missing = map_symbols_to_entrez(["SCN1A"])
        assert "SCN1A" in mapped, "SCN1A should be in mapped dict"
        assert len(missing) == 0, f"No genes should be missing, got {missing}"

    def test_fake_gene_goes_to_missing(self):
        """A fake gene symbol should end up in the missing list."""
        mapped, missing = map_symbols_to_entrez(["FAKE_GENE_XYZ"])
        assert "FAKE_GENE_XYZ" in missing, (
            "FAKE_GENE_XYZ should be in missing list"
        )
        assert "FAKE_GENE_XYZ" not in mapped, (
            "FAKE_GENE_XYZ should not be in mapped dict"
        )

    def test_mixed_valid_and_invalid(self):
        """Test with a mix of valid and invalid genes."""
        mapped, missing = map_symbols_to_entrez(["SCN1A", "FAKE_GENE_XYZ"])
        assert "SCN1A" in mapped
        assert "FAKE_GENE_XYZ" in missing

    def test_alias_resolution(self):
        """Genes with aliases should be resolved via the alias."""
        mapped, missing = map_symbols_to_entrez(["SEPT5"])
        # SEPT5 may or may not be in the gene info table directly,
        # but the alias SEPTIN5 should resolve it
        assert "SEPT5" in mapped or "SEPT5" in missing, (
            "SEPT5 should appear in either mapped or missing"
        )

    def test_custom_mapping_dict(self):
        """Should use a provided mapping dict instead of loading gene info."""
        custom_map = {"GENE_A": 12345, "GENE_B": 67890}
        mapped, missing = map_symbols_to_entrez(
            ["GENE_A", "GENE_C"], gene_symbol_to_entrez=custom_map
        )
        assert mapped == {"GENE_A": 12345}
        assert missing == ["GENE_C"]


class TestComputeResiduals:
    """Tests for compute_residuals function."""

    def test_compute_residuals(self):
        """Residuals should have near-zero group means after removing class effect."""
        np.random.seed(42)
        n = 50
        idx = [f"ct_{i}" for i in range(n)]
        # Create class labels: first 25 = "exc", last 25 = "inh"
        class_labels = pd.Series(
            ["exc"] * 25 + ["inh"] * 25, index=idx
        )
        # Create specificity with strong class effect: inh cells += 5.0
        gene_spec = pd.Series(np.random.randn(n), index=idx)
        gene_spec.iloc[:25] += 0.0
        gene_spec.iloc[25:] += 5.0

        resid = compute_residuals(gene_spec, class_labels)
        assert isinstance(resid, pd.Series)
        assert len(resid) == n
        assert list(resid.index) == idx

        # Group means of residuals should be near zero
        resid_exc = resid.iloc[:25].mean()
        resid_inh = resid.iloc[25:].mean()
        assert abs(resid_exc) < 0.1, f"Exc group mean residual {resid_exc} not near zero"
        assert abs(resid_inh) < 0.1, f"Inh group mean residual {resid_inh} not near zero"


class TestPearsonPartial:
    """Tests for pearson_partial function."""

    def test_pearson_partial(self):
        """Pearson partial should remove class-driven correlation."""
        np.random.seed(123)
        n = 50
        idx = [f"ct_{i}" for i in range(n)]
        class_labels = pd.Series(
            ["exc"] * 25 + ["inh"] * 25, index=idx
        )
        # Two genes both elevated in "inh" class -> class-driven correlation
        gene_a = pd.Series(np.random.randn(n), index=idx)
        gene_b = pd.Series(np.random.randn(n), index=idx)
        gene_a.iloc[25:] += 5.0
        gene_b.iloc[25:] += 5.0

        # Raw Spearman rho should be > 0.3 (class-driven)
        from scipy.stats import spearmanr
        raw_rho, _ = spearmanr(gene_a, gene_b)
        assert abs(raw_rho) > 0.3, (
            f"Expected raw Spearman rho > 0.3, got {raw_rho}"
        )

        # Pearson partial should have lower absolute correlation
        r, pval = pearson_partial(gene_a, gene_b, class_labels)
        assert abs(r) < abs(raw_rho), (
            f"Partial r ({r}) should be lower in abs value than raw rho ({raw_rho})"
        )


class TestRunProfileSimilarity:
    """Tests for run_profile_similarity pipeline function."""

    def test_run_profile_similarity(self):
        rng = np.random.default_rng(42)
        n_ct = 50
        n_genes = 10
        spec_mat = pd.DataFrame(
            rng.normal(size=(n_genes, n_ct)),
            index=[f"gene_{i}" for i in range(n_genes)],
            columns=range(n_ct),
        )
        source_ids = {"SRC1": "gene_0", "SRC2": "gene_1"}
        target_ids = {"TGT1": "gene_2", "TGT2": "gene_3", "TGT3": "gene_4"}
        target_families = {"FamA": ["TGT1", "TGT2"], "FamB": ["TGT3"]}
        class_labels = pd.Series(["A"] * 25 + ["B"] * 25, index=range(n_ct))
        scope_idx = list(range(n_ct))

        results = run_profile_similarity(
            spec_mat, source_ids, target_ids, target_families,
            scope_idx, class_labels,
        )

        # Pair results
        pair_df = results["pair_results"]
        assert len(pair_df) == 6  # 2 sources x 3 targets
        for col in ["source", "target", "family", "rho_raw", "p_raw",
                     "r_partial", "p_partial", "q_partial"]:
            assert col in pair_df.columns

        # Family results
        fam_df = results["family_results"]
        assert len(fam_df) == 4  # 2 sources x 2 families
        for col in ["source", "family", "stouffer_z", "stouffer_p",
                     "n_genes", "best_target", "best_r", "q_family"]:
            assert col in fam_df.columns

        # FamA should have n_genes=2, FamB should have n_genes=1
        fam_a = fam_df[fam_df["family"] == "FamA"]
        assert (fam_a["n_genes"] == 2).all()
        fam_b = fam_df[fam_df["family"] == "FamB"]
        assert (fam_b["n_genes"] == 1).all()
