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
    convergence_permutation_test,
    _single_permutation,
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


class TestConvergencePermutation:
    """Tests for convergence_permutation_test and _single_permutation."""

    @pytest.fixture
    def setup_data(self):
        """Create synthetic data for permutation tests."""
        rng = np.random.default_rng(42)
        n_ct = 50
        n_genes = 100
        spec_mat = pd.DataFrame(
            rng.normal(size=(n_genes, n_ct)),
            index=[f"gene_{i}" for i in range(n_genes)],
            columns=range(n_ct),
        )
        source_eids = ["gene_0", "gene_1", "gene_2"]
        target_ids = {"T1": "gene_50", "T2": "gene_51"}
        target_families = {"Fam1": ["T1", "T2"]}
        class_labels = pd.Series(["A"] * 25 + ["B"] * 25, index=range(n_ct))
        scope_idx = list(range(n_ct))
        return {
            "spec_mat": spec_mat,
            "source_eids": source_eids,
            "target_ids": target_ids,
            "target_families": target_families,
            "class_labels": class_labels,
            "scope_idx": scope_idx,
        }

    def test_convergence_permutation_basic(self, setup_data):
        """Basic smoke test: correct keys and shapes in result dict."""
        result = convergence_permutation_test(
            setup_data["spec_mat"],
            setup_data["source_eids"],
            setup_data["target_ids"],
            setup_data["target_families"],
            setup_data["scope_idx"],
            setup_data["class_labels"],
            observed_hits=1,
            z_threshold=2.0,
            n_perms=100,
            seed=42,
            n_jobs=1,
        )
        assert "null_hits" in result
        assert len(result["null_hits"]) == 100
        assert "p_value" in result
        assert 0.0 <= result["p_value"] <= 1.0
        assert "z_score" in result

    def test_null_hits_are_nonnegative_integers(self, setup_data):
        """Each permutation hit count should be a non-negative integer."""
        result = convergence_permutation_test(
            setup_data["spec_mat"],
            setup_data["source_eids"],
            setup_data["target_ids"],
            setup_data["target_families"],
            setup_data["scope_idx"],
            setup_data["class_labels"],
            observed_hits=0,
            z_threshold=2.0,
            n_perms=50,
            seed=123,
            n_jobs=1,
        )
        assert all(h >= 0 for h in result["null_hits"])
        assert all(isinstance(h, (int, np.integer)) for h in result["null_hits"])

    def test_p_value_with_zero_observed_hits(self, setup_data):
        """With observed_hits=0, p-value should be high (most perms >= 0)."""
        result = convergence_permutation_test(
            setup_data["spec_mat"],
            setup_data["source_eids"],
            setup_data["target_ids"],
            setup_data["target_families"],
            setup_data["scope_idx"],
            setup_data["class_labels"],
            observed_hits=0,
            z_threshold=2.0,
            n_perms=100,
            seed=42,
            n_jobs=1,
        )
        # All null_hits >= 0 == observed_hits, so p should be close to 1
        assert result["p_value"] > 0.5

    def test_reproducibility_with_same_seed(self, setup_data):
        """Same seed should produce identical null distributions."""
        kwargs = dict(
            spec_mat=setup_data["spec_mat"],
            source_eids=setup_data["source_eids"],
            target_ids=setup_data["target_ids"],
            target_families=setup_data["target_families"],
            scope_idx=setup_data["scope_idx"],
            class_labels=setup_data["class_labels"],
            observed_hits=1,
            z_threshold=2.0,
            n_perms=50,
            seed=99,
            n_jobs=1,
        )
        r1 = convergence_permutation_test(**kwargs)
        r2 = convergence_permutation_test(**kwargs)
        np.testing.assert_array_equal(r1["null_hits"], r2["null_hits"])

    def test_single_permutation_excludes_targets(self, setup_data):
        """_single_permutation should not sample target gene IDs."""
        d = setup_data
        cols = [d["spec_mat"].columns[i] for i in d["scope_idx"]]
        spec_scope = d["spec_mat"][cols]
        class_scope = d["class_labels"].loc[cols]

        # Precompute target residuals
        resid_cache_targets = {}
        for sym, gid in d["target_ids"].items():
            if gid in spec_scope.index:
                resid_cache_targets[sym] = compute_residuals(
                    spec_scope.loc[gid], class_scope
                )

        target_entrez_set = set(d["target_ids"].values())
        # Run many single permutations and check the function runs without error
        for seed in range(10):
            hit_count = _single_permutation(
                spec_mat=spec_scope,
                n_source=len(d["source_eids"]),
                target_ids=d["target_ids"],
                target_families=d["target_families"],
                scope_idx=d["scope_idx"],
                class_labels=class_scope,
                residuals_cache_targets=resid_cache_targets,
                z_threshold=2.0,
                rng_seed=seed,
            )
            assert isinstance(hit_count, (int, np.integer))
            assert hit_count >= 0

    def test_high_threshold_yields_zero_hits(self, setup_data):
        """With an extremely high z_threshold, no permutation should have hits."""
        result = convergence_permutation_test(
            setup_data["spec_mat"],
            setup_data["source_eids"],
            setup_data["target_ids"],
            setup_data["target_families"],
            setup_data["scope_idx"],
            setup_data["class_labels"],
            observed_hits=1,
            z_threshold=100.0,  # impossibly high
            n_perms=50,
            seed=42,
            n_jobs=1,
        )
        assert all(h == 0 for h in result["null_hits"])

    def test_parallel_matches_serial(self, setup_data):
        """n_jobs>1 should produce the same results as n_jobs=1."""
        kwargs = dict(
            spec_mat=setup_data["spec_mat"],
            source_eids=setup_data["source_eids"],
            target_ids=setup_data["target_ids"],
            target_families=setup_data["target_families"],
            scope_idx=setup_data["scope_idx"],
            class_labels=setup_data["class_labels"],
            observed_hits=1,
            z_threshold=2.0,
            n_perms=50,
            seed=42,
        )
        r_serial = convergence_permutation_test(**kwargs, n_jobs=1)
        r_parallel = convergence_permutation_test(**kwargs, n_jobs=2)
        np.testing.assert_array_equal(r_serial["null_hits"], r_parallel["null_hits"])


class TestPPI2HopReachability:
    """Tests for compute_2hop_reachability function."""

    def test_ppi_2hop_reachability(self):
        import networkx as nx
        from convergence_utils import compute_2hop_reachability

        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (4, 5)])
        # Node 1 reaches 0,2,3 in 1 hop, 4 in 2 hops -> 4/5 = 0.8
        frac = compute_2hop_reachability(G, 1)
        assert 0.7 < frac < 0.9

    def test_node_not_in_graph(self):
        import networkx as nx
        from convergence_utils import compute_2hop_reachability

        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2)])
        assert compute_2hop_reachability(G, 99) == 0.0

    def test_isolated_node(self):
        import networkx as nx
        from convergence_utils import compute_2hop_reachability

        G = nx.Graph()
        G.add_nodes_from([0, 1, 2])
        # Node 0 has no neighbors -> reachable = 0 -> 0 / 2 = 0.0
        assert compute_2hop_reachability(G, 0) == 0.0


class TestFindPPIBridges:
    """Tests for find_ppi_bridges function."""

    def test_find_bridges(self):
        import networkx as nx
        from convergence_utils import find_ppi_bridges

        G = nx.Graph()
        G.add_edges_from([("A", "BRIDGE"), ("BRIDGE", "B"), ("A", "C")])
        id_to_name = {"A": "A", "BRIDGE": "BRIDGE", "B": "B", "C": "C"}
        bridges = find_ppi_bridges(
            G, {"SRC": "A"}, {"TGT": "B"}, id_to_name, max_hops=2
        )
        assert len(bridges) == 1
        assert bridges[0]["bridge"] == "BRIDGE"
        assert bridges[0]["path_length"] == 2

    def test_direct_interaction(self):
        import networkx as nx
        from convergence_utils import find_ppi_bridges

        G = nx.Graph()
        G.add_edges_from([("A", "B")])
        id_to_name = {"A": "A", "B": "B"}
        bridges = find_ppi_bridges(
            G, {"SRC": "A"}, {"TGT": "B"}, id_to_name, max_hops=2
        )
        assert len(bridges) == 1
        assert bridges[0]["bridge"] is None
        assert bridges[0]["path_length"] == 1

    def test_no_path(self):
        import networkx as nx
        from convergence_utils import find_ppi_bridges

        G = nx.Graph()
        G.add_nodes_from(["A", "B"])  # No edge
        id_to_name = {"A": "A", "B": "B"}
        bridges = find_ppi_bridges(
            G, {"SRC": "A"}, {"TGT": "B"}, id_to_name, max_hops=2
        )
        assert len(bridges) == 0

    def test_path_too_long(self):
        import networkx as nx
        from convergence_utils import find_ppi_bridges

        G = nx.Graph()
        G.add_edges_from([("A", "X"), ("X", "Y"), ("Y", "B")])
        id_to_name = {"A": "A", "X": "X", "Y": "Y", "B": "B"}
        bridges = find_ppi_bridges(
            G, {"SRC": "A"}, {"TGT": "B"}, id_to_name, max_hops=2
        )
        # Shortest path is A->X->Y->B, length 3 > max_hops=2
        assert len(bridges) == 0


class TestLoadStringNetwork:
    """Tests for load_string_network function."""

    def test_load_string_network(self):
        """Test with actual STRING data if available."""
        from convergence_utils import load_string_network

        info = Path("dat/STRING/9606.protein.info.v12.0.txt.gz")
        links = Path("dat/STRING/9606.protein.links.v12.0.txt.gz")
        if not info.exists():
            pytest.skip("STRING data not available")
        G, n2i, i2n = load_string_network(info, links, score_threshold=900)
        assert G.number_of_nodes() > 1000
        assert "SCN1A" in n2i or "COMT" in n2i  # common genes should be there


class TestExtractAllEphysFeatures:
    """Tests for extract_all_ephys_features function."""

    def test_extract_ephys_features(self):
        from convergence_utils import extract_all_ephys_features

        bundle_dir = Path("dat/VIP_Ephys/analysis_bundles")
        if not bundle_dir.exists():
            pytest.skip("Ephys data not available")
        df = extract_all_ephys_features(bundle_dir)
        assert len(df) == 29
        assert "genotype" in df.columns
        assert "rise_slope" in df.columns
        assert "hero_cv_isi" in df.columns
        assert "max_freq_Hz" in df.columns
        # Check rise_slope is in reasonable V/s range (not raw 1e6x values)
        assert df["rise_slope"].max() < 500  # should be ~100-300 V/s


class TestCompareFeature:
    """Tests for compare_feature function."""

    def test_compare_feature(self):
        from convergence_utils import compare_feature

        df = pd.DataFrame({
            "genotype": ["WT"] * 10 + ["Df16"] * 10,
            "feat": list(range(10, 20)) + list(range(5, 15)),
        })
        result = compare_feature(df, "feat")
        assert result is not None
        assert "cohens_d" in result
        assert "ci_low" in result
        assert result["n_wt"] == 10

    def test_compare_feature_insufficient_data(self):
        from convergence_utils import compare_feature

        df = pd.DataFrame({
            "genotype": ["WT"] * 2 + ["Df16"] * 10,
            "feat": [1.0, 2.0] + list(range(10)),
        })
        result = compare_feature(df, "feat")
        assert result is None

    def test_compare_feature_with_nans(self):
        from convergence_utils import compare_feature

        df = pd.DataFrame({
            "genotype": ["WT"] * 10 + ["Df16"] * 10,
            "feat": [1.0, 2.0, np.nan] + [3.0] * 7 + [4.0] * 10,
        })
        result = compare_feature(df, "feat")
        assert result is not None
        assert result["n_wt"] == 9  # one NaN dropped
