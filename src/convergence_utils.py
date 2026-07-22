"""
Convergence analysis utilities: gene set constants, mapping functions,
and partial-correlation helpers.

Defines curated ion channel / scaffolding gene families, 22q11.2 deletion
genes, ephys feature metadata, and statistical functions for computing
Pearson partial correlations controlling for cell-class membership.
"""

import logging
import os
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr, spearmanr, norm, ttest_ind
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests
from statsmodels.tools import add_constant

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 22q11.2 deletion region genes (19 genes)
# ---------------------------------------------------------------------------
GENES_22Q = [
    "DGCR2", "DGCR8", "DGCR14", "TBX1", "COMT", "PRODH",
    "GP1BB", "SEPT5", "UFD1L", "CDC45", "CLDN5", "HIRA",
    "MRPL40", "SLC25A1", "TANGO2", "TXNRD2", "ZDHHC8",
    "RTN4R", "RANBP1",
]

# Old HGNC symbols -> current approved symbols
GENE_ALIASES = {
    "SEPT5": "SEPTIN5",
    "UFD1L": "UFD1",
}

# ---------------------------------------------------------------------------
# Curated ion-channel / scaffolding gene families
# ---------------------------------------------------------------------------
CHANNEL_FAMILIES = {
    # --- Rise slope / Peak voltage ---
    "Nav-alpha": {
        "genes": ["SCN1A", "SCN2A", "SCN3A", "SCN8A"],
        "feature": "Rise slope/Peak voltage",
    },
    "Nav-beta": {
        "genes": ["SCN1B", "SCN2B", "SCN3B", "SCN4B"],
        "feature": "Rise slope/Peak voltage",
    },
    # --- AP width / Decay slope ---
    "Kv3": {
        "genes": ["KCNC1", "KCNC2", "KCNC3", "KCNC4"],
        "feature": "AP width/Decay slope",
    },
    "BK": {
        "genes": ["KCNMA1"],
        "feature": "AP width/Decay slope",
    },
    "Kv1": {
        "genes": ["KCNA1", "KCNA2", "KCNA3"],
        "feature": "AP width/Decay slope",
    },
    # --- RMP / Input resistance ---
    "K2P": {
        "genes": ["KCNK2", "KCNK3"],
        "feature": "RMP/Input resistance",
    },
    "HCN": {
        "genes": ["HCN1", "HCN2"],
        "feature": "RMP/Input resistance",
    },
    "Kir": {
        "genes": ["KCNJ2", "KCNJ4"],
        "feature": "RMP/Input resistance",
    },
    # --- Scaffolding ---
    "Nav-scaffold": {
        "genes": ["ANK3", "FGF14", "NFASC", "RANGRF"],
        "feature": "Scaffolding",
    },
    "Kv-scaffold": {
        "genes": ["DLG1", "DLG4", "CNTNAP2", "AKAP5"],
        "feature": "Scaffolding",
    },
    "Leak-scaffold": {
        "genes": ["PEX5L", "SHANK1", "SHANK2", "SHANK3"],
        "feature": "Scaffolding",
    },
    # --- Reverse prediction families ---
    "SK": {
        "genes": ["KCNN1", "KCNN2", "KCNN3"],
        "feature": "Reverse prediction",
    },
    "SK-aux": {
        "genes": ["CALM1", "CALM2", "CALM3"],
        "feature": "Reverse prediction",
    },
    "Kv2": {
        "genes": ["KCNB1", "KCNB2"],
        "feature": "Reverse prediction",
    },
    "Kv2-aux": {
        "genes": ["VAPA", "VAPB"],
        "feature": "Reverse prediction",
    },
    "Kv4-aux": {
        "genes": ["KCNIP1", "KCNIP2", "KCNIP3", "KCNIP4", "DPP6", "DPP10"],
        "feature": "Reverse prediction",
    },
}

# Flat dict: {family_name: [gene_symbols]} derived from CHANNEL_FAMILIES
CURATED_TARGETS = {
    family: info["genes"] for family, info in CHANNEL_FAMILIES.items()
}

# ---------------------------------------------------------------------------
# Ephys feature metadata
# ---------------------------------------------------------------------------
# 6 significant/marginal features + 3 reverse-prediction features
EPHYS_FEATURES = {
    "AP_rise_slope": {
        "label": "AP rise slope",
        "direction": "decreased",
        "p": 0.007,
        "d": -0.45,
    },
    "AP_peak_voltage": {
        "label": "AP peak voltage",
        "direction": "decreased",
        "p": 0.012,
        "d": -0.42,
    },
    "AP_width": {
        "label": "AP width",
        "direction": "increased",
        "p": 0.003,
        "d": 0.50,
    },
    "AP_decay_slope": {
        "label": "AP decay slope",
        "direction": "decreased",
        "p": 0.015,
        "d": -0.40,
    },
    "RMP": {
        "label": "Resting membrane potential",
        "direction": "depolarized",
        "p": 0.045,
        "d": 0.33,
    },
    "input_resistance": {
        "label": "Input resistance",
        "direction": "increased",
        "p": 0.062,
        "d": 0.31,
    },
    # Reverse predictions
    "AHP_amplitude": {
        "label": "AHP amplitude",
        "direction": "predicted_decreased",
        "p": None,
        "d": None,
    },
    "firing_rate": {
        "label": "Firing rate",
        "direction": "predicted_decreased",
        "p": None,
        "d": None,
    },
    "ISI_CV": {
        "label": "ISI coefficient of variation",
        "direction": "predicted_increased",
        "p": None,
        "d": None,
    },
}


# ---------------------------------------------------------------------------
# Mapping function
# ---------------------------------------------------------------------------
def map_symbols_to_entrez(symbols, gene_symbol_to_entrez=None):
    """Map gene symbols to Entrez IDs.

    Parameters
    ----------
    symbols : list of str
        Gene symbols to map.
    gene_symbol_to_entrez : dict or None
        Pre-built {symbol: entrez_id} mapping. If None, loads via
        ``CellType_PSY.LoadGeneINFO()``.

    Returns
    -------
    mapped : dict
        ``{symbol: entrez_id}`` for successfully mapped genes.
    missing : list
        Symbols that could not be mapped.
    """
    if gene_symbol_to_entrez is None:
        from CellType_PSY import LoadGeneINFO
        _hgnc, _ens2ent, gene_symbol_to_entrez, _ent2sym = LoadGeneINFO()

    mapped = {}
    missing = []
    for sym in symbols:
        # Try direct lookup first
        if sym in gene_symbol_to_entrez:
            mapped[sym] = gene_symbol_to_entrez[sym]
        # Try alias
        elif sym in GENE_ALIASES and GENE_ALIASES[sym] in gene_symbol_to_entrez:
            mapped[sym] = gene_symbol_to_entrez[GENE_ALIASES[sym]]
        else:
            missing.append(sym)
    return mapped, missing


# ---------------------------------------------------------------------------
# Partial correlation helpers
# ---------------------------------------------------------------------------
def compute_residuals(gene_spec, class_labels):
    """Regress out cell-class membership from gene specificity values.

    Fits an OLS model: gene_spec ~ one-hot(class_labels) + constant,
    and returns the residuals.

    Parameters
    ----------
    gene_spec : pd.Series
        Gene specificity values indexed by cell-type ID.
    class_labels : pd.Series
        Cell-class labels (e.g. "exc", "inh", "glia") with the same index
        as *gene_spec*.

    Returns
    -------
    pd.Series
        OLS residuals, same index as *gene_spec*.
    """
    dummies = pd.get_dummies(class_labels, drop_first=True, dtype=float)
    X = add_constant(dummies)
    model = OLS(gene_spec.values, X).fit()
    return pd.Series(model.resid, index=gene_spec.index)


def pearson_partial(gene_a, gene_b, class_labels):
    """Pearson correlation between two genes after regressing out cell class.

    Computes OLS residuals for each gene (removing the linear effect of
    cell-class dummies) and then returns Pearson *r* and *p*-value on the
    residuals.

    **Methodological note**: OLS residualization is a linear operation that
    destroys rank structure. Pearson correlation on residuals is a
    well-defined partial correlation. Spearman on OLS residuals would be
    methodologically unsound because rank ordering is not preserved under
    linear projection.

    Parameters
    ----------
    gene_a, gene_b : pd.Series
        Gene specificity values indexed by cell-type ID.
    class_labels : pd.Series
        Cell-class labels with the same index.

    Returns
    -------
    r : float
        Pearson correlation coefficient on residuals.
    pvalue : float
        Two-sided p-value for the correlation.
    """
    resid_a = compute_residuals(gene_a, class_labels)
    resid_b = compute_residuals(gene_b, class_labels)
    r, pvalue = pearsonr(resid_a, resid_b)
    return r, pvalue


# ---------------------------------------------------------------------------
# Profile similarity pipeline
# ---------------------------------------------------------------------------
def run_profile_similarity(spec_mat, source_ids, target_ids, target_families,
                           scope_idx, class_labels):
    """Run full profile similarity analysis (Tiers 1-2) with family collapsing.

    For every source x target gene pair, computes:
      - Tier 1: Spearmans' R on raw specificity profiles
      - Tier 2: Pearson r on OLS residuals (controlling for cell class)

    Results are FDR-corrected and collapsed to families via Stouffer's
    signed z method.

    Parameters
    ----------
    spec_mat : pd.DataFrame
        Genes (rows) x cell types (columns). Index = gene IDs (e.g. Entrez).
    source_ids : dict
        {symbol: gene_id} for source genes (22q).
    target_ids : dict
        {symbol: gene_id} for target genes (ion channels, scaffolds).
    target_families : dict
        {family_name: [symbol, ...]} grouping targets into families.
    scope_idx : list of int
        Column indices (cell types) to restrict analysis to.
    class_labels : pd.Series
        Cell-class labels indexed by cell-type ID, matching the columns
        selected by *scope_idx*.

    Returns
    -------
    dict with keys:
        ``"pair_results"`` : pd.DataFrame
            Columns: source, target, family, rho_raw, p_raw, r_partial,
            p_partial, q_partial
        ``"family_results"`` : pd.DataFrame
            Columns: source, family, stouffer_z, stouffer_p, n_genes,
            best_target, best_r, q_family
    """
    # --- 1. Build symbol -> family lookup ---
    symbol_to_family = {}
    for fam_name, symbols in target_families.items():
        for sym in symbols:
            symbol_to_family[sym] = fam_name

    # --- 2. Subset spec_mat to scope and precompute residuals ---
    cols = [spec_mat.columns[i] for i in scope_idx]
    spec_scope = spec_mat[cols]
    class_scope = class_labels.loc[cols]

    # Precompute OLS residuals for ALL genes in scope (cache per gene_id)
    all_gene_ids = set()
    for gid in source_ids.values():
        all_gene_ids.add(gid)
    for gid in target_ids.values():
        all_gene_ids.add(gid)

    resid_cache = {}
    for gid in all_gene_ids:
        if gid in spec_scope.index:
            gene_spec = spec_scope.loc[gid]
            resid_cache[gid] = compute_residuals(gene_spec, class_scope)

    # --- 3. Pairwise analysis ---
    pair_records = []
    for src_sym, src_id in source_ids.items():
        if src_id not in spec_scope.index:
            continue
        src_raw = spec_scope.loc[src_id]
        src_resid = resid_cache[src_id]

        for tgt_sym, tgt_id in target_ids.items():
            if tgt_id not in spec_scope.index:
                continue
            tgt_raw = spec_scope.loc[tgt_id]
            tgt_resid = resid_cache[tgt_id]

            # Tier 1: Spearmans' R on raw specificity
            rho_raw, p_raw = spearmanr(src_raw.values, tgt_raw.values)

            # Tier 2: Pearson on cached residuals
            r_partial, p_partial = pearsonr(src_resid.values, tgt_resid.values)

            pair_records.append({
                "source": src_sym,
                "target": tgt_sym,
                "family": symbol_to_family.get(tgt_sym, "unknown"),
                "rho_raw": rho_raw,
                "p_raw": p_raw,
                "r_partial": r_partial,
                "p_partial": p_partial,
            })

    pair_df = pd.DataFrame(pair_records)

    # --- 4. BH FDR on pair-level p_partial ---
    if len(pair_df) > 0:
        _, q_vals, _, _ = multipletests(pair_df["p_partial"], method="fdr_bh")
        pair_df["q_partial"] = q_vals
    else:
        pair_df["q_partial"] = []

    # --- 5. Family-level collapsing via Stouffer's signed z ---
    family_records = []
    if len(pair_df) > 0:
        for (src_sym, fam_name), group in pair_df.groupby(["source", "family"]):
            n = len(group)
            # Convert each pair's two-sided p_partial to signed z
            # CRITICAL: pearsonr returns two-sided p-values.
            # Use p/2 to get one-tail probability, then norm.ppf(1 - p/2)
            # gives the z-magnitude. Multiply by sign(r) for directionality.
            signed_zs = np.sign(group["r_partial"].values) * norm.ppf(
                1 - group["p_partial"].values / 2
            )
            # Stouffer's combined z
            stouffer_z = np.sum(signed_zs) / np.sqrt(n)
            # Two-sided p-value for the combined z
            stouffer_p = 2 * (1 - norm.cdf(np.abs(stouffer_z)))

            # Best individual target (highest |r_partial|)
            best_idx = group["r_partial"].abs().idxmax()
            best_target = group.loc[best_idx, "target"]
            best_r = group.loc[best_idx, "r_partial"]

            family_records.append({
                "source": src_sym,
                "family": fam_name,
                "stouffer_z": stouffer_z,
                "stouffer_p": stouffer_p,
                "n_genes": n,
                "best_target": best_target,
                "best_r": best_r,
            })

    fam_df = pd.DataFrame(family_records)

    # BH FDR on family-level p-values
    if len(fam_df) > 0:
        _, q_fam, _, _ = multipletests(fam_df["stouffer_p"], method="fdr_bh")
        fam_df["q_family"] = q_fam
    else:
        fam_df["q_family"] = []

    return {
        "pair_results": pair_df,
        "family_results": fam_df,
    }


# ---------------------------------------------------------------------------
# Convergence permutation test
# ---------------------------------------------------------------------------
def _single_permutation(spec_mat, n_source, target_ids, target_families,
                        scope_idx, class_labels, residuals_cache_targets,
                        z_threshold, rng_seed):
    """Run one permutation iteration and return the number of family-level hits.

    Samples ``n_source`` random genes from *spec_mat* (excluding target gene
    Entrez IDs), computes OLS residuals for each, then for every sampled-gene x
    target-family pair computes Stouffer's signed z.  Returns the count of
    source x family pairs where ``|Z| >= z_threshold``.

    Parameters
    ----------
    spec_mat : pd.DataFrame
        Expression specificity matrix (genes x cell types), already scoped
        to the cell types of interest.
    n_source : int
        Number of source genes to sample (matches the real source gene set size).
    target_ids : dict
        ``{symbol: gene_id}`` for target genes.
    target_families : dict
        ``{family_name: [symbol, ...]}`` grouping targets into families.
    scope_idx : list of int
        Column indices (unused here; kept for API consistency).
    class_labels : pd.Series
        Cell-class labels indexed by cell-type ID, matching *spec_mat* columns.
    residuals_cache_targets : dict
        ``{symbol: pd.Series}`` of precomputed OLS residuals for each target gene.
    z_threshold : float
        Minimum ``|Stouffer Z|`` to count a source x family pair as a hit.
    rng_seed : int
        Seed for the random number generator (reproducibility).

    Returns
    -------
    int
        Total number of source x family hits in this permutation.
    """
    rng = np.random.default_rng(rng_seed)

    # Build set of target Entrez IDs to exclude from sampling
    target_entrez_set = set(target_ids.values())

    # Pool of eligible genes (everything in spec_mat except targets)
    eligible = [gid for gid in spec_mat.index if gid not in target_entrez_set]

    # Sample n_source genes
    sampled_idx = rng.choice(len(eligible), size=n_source, replace=False)
    sampled_genes = [eligible[i] for i in sampled_idx]

    # Compute OLS residuals for each sampled gene
    sampled_resids = {}
    for gid in sampled_genes:
        gene_spec = spec_mat.loc[gid]
        sampled_resids[gid] = compute_residuals(gene_spec, class_labels)

    # For each sampled gene x family, compute Stouffer's signed z
    hit_count = 0
    for gid in sampled_genes:
        src_resid = sampled_resids[gid]

        for fam_name, symbols in target_families.items():
            signed_zs = []
            for sym in symbols:
                if sym not in residuals_cache_targets:
                    continue
                tgt_resid = residuals_cache_targets[sym]

                r, p = pearsonr(src_resid.values, tgt_resid.values)

                # Convert two-sided p to signed z
                # Clamp p away from 0 to avoid inf
                p_clamped = max(p, 1e-300)
                z_i = np.sign(r) * norm.ppf(1 - p_clamped / 2)
                signed_zs.append(z_i)

            if len(signed_zs) == 0:
                continue

            # Stouffer's combined z
            n_genes = len(signed_zs)
            stouffer_z = np.sum(signed_zs) / np.sqrt(n_genes)

            if np.abs(stouffer_z) >= z_threshold:
                hit_count += 1

    return int(hit_count)


def convergence_permutation_test(spec_mat, source_eids, target_ids,
                                 target_families, scope_idx, class_labels,
                                 observed_hits, z_threshold,
                                 n_perms=10000, seed=42, n_jobs=10):
    """Permutation test for convergence significance using Stouffer's z.

    Assesses whether the observed number of significant source x family hits
    (defined by ``|Stouffer Z| >= z_threshold``) exceeds what is expected by
    chance, by repeatedly sampling random gene sets of the same size as the
    source set and recomputing the test statistic.

    Parameters
    ----------
    spec_mat : pd.DataFrame
        Genes (rows) x cell types (columns).  Index = gene IDs.
    source_eids : list
        Gene IDs (Entrez or equivalent) for the source gene set (e.g. 22q genes).
    target_ids : dict
        ``{symbol: gene_id}`` for target genes.
    target_families : dict
        ``{family_name: [symbol, ...]}`` grouping targets into families.
    scope_idx : list of int
        Column indices (cell types) to restrict analysis to.
    class_labels : pd.Series
        Cell-class labels indexed by cell-type ID.
    observed_hits : int
        Number of source x family pairs with ``|Stouffer Z| >= z_threshold``
        in the real data.
    z_threshold : float
        Minimum ``|Stouffer Z|`` for a hit (calibrated from observed data,
        typically the min |stouffer_z| among family-level hits with
        q_family < 0.05).
    n_perms : int
        Number of permutations (default 10,000).
    seed : int
        Master random seed for reproducibility.
    n_jobs : int
        Number of parallel workers (default 10).

    Returns
    -------
    dict
        ``"null_hits"`` : np.ndarray of int, length *n_perms*
            Hit count from each permutation.
        ``"p_value"`` : float
            Empirical p-value: ``(sum(null >= observed) + 1) / (n_perms + 1)``.
        ``"z_score"`` : float
            Z-score of observed vs null distribution:
            ``(observed - mean(null)) / std(null)``.
    """
    # --- 1. Scope the data ---
    cols = [spec_mat.columns[i] for i in scope_idx]
    spec_scope = spec_mat[cols]
    class_scope = class_labels.loc[cols]

    # --- 2. Precompute target residuals (invariant across permutations) ---
    residuals_cache_targets = {}
    for sym, gid in target_ids.items():
        if gid in spec_scope.index:
            residuals_cache_targets[sym] = compute_residuals(
                spec_scope.loc[gid], class_scope
            )

    n_source = len(source_eids)

    # --- 3. Generate reproducible per-permutation seeds ---
    master_rng = np.random.default_rng(seed)
    perm_seeds = master_rng.integers(0, 2**31, size=n_perms).tolist()

    # --- 4. Run permutations in parallel ---
    null_hits_list = Parallel(n_jobs=n_jobs)(
        delayed(_single_permutation)(
            spec_mat=spec_scope,
            n_source=n_source,
            target_ids=target_ids,
            target_families=target_families,
            scope_idx=scope_idx,
            class_labels=class_scope,
            residuals_cache_targets=residuals_cache_targets,
            z_threshold=z_threshold,
            rng_seed=s,
        )
        for s in perm_seeds
    )
    null_hits = np.array(null_hits_list, dtype=int)

    # --- 5. Compute empirical p-value (with pseudocount) ---
    p_value = (np.sum(null_hits >= observed_hits) + 1) / (n_perms + 1)

    # --- 6. Compute z-score ---
    null_mean = np.mean(null_hits)
    null_std = np.std(null_hits)
    if null_std > 0:
        z_score = (observed_hits - null_mean) / null_std
    else:
        # All null hits identical; use sign-based z
        z_score = np.inf if observed_hits > null_mean else 0.0

    return {
        "null_hits": null_hits,
        "p_value": float(p_value),
        "z_score": float(z_score),
    }


# ---------------------------------------------------------------------------
# GO-derived ion channel gene set
# ---------------------------------------------------------------------------

# Target GO terms for ion channel genes
_GO_ION_CHANNEL_TERMS = [
    "GO:0005216",  # ion channel activity
    "GO:0034765",  # regulation of ion transmembrane transport
]

_OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
_GAF_URL = "http://geneontology.org/gene-associations/goa_human.gaf.gz"


def _find_or_download(cache_dir, prefix, suffix, url):
    """Find a date-pinned file or download it.

    Looks for ``<prefix>.<YYYY-MM><suffix>`` in *cache_dir*.  If a file
    matching the current month exists, reuse it.  If only an older-dated
    file exists, warn and reuse it.  Otherwise, download from *url*.

    Returns the path to the (possibly just-downloaded) file.
    """
    import urllib.request

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    current_month = datetime.now().strftime("%Y-%m")
    current_name = f"{prefix}.{current_month}{suffix}"
    current_path = cache_dir / current_name

    if current_path.exists():
        logger.info("Reusing cached %s", current_path)
        return current_path

    # Check for any older-dated file
    pattern = f"{prefix}.*{suffix}"
    import glob as _glob
    existing = sorted(_glob.glob(str(cache_dir / f"{prefix}.*{suffix}")))
    if existing:
        old_path = existing[-1]  # most recent by name sort
        warnings.warn(
            f"Found older GO file {os.path.basename(old_path)} "
            f"(current month is {current_month}). Reusing it. "
            f"Delete it and re-run to force a fresh download.",
            stacklevel=2,
        )
        return Path(old_path)

    # Download (use a proper User-Agent to avoid 403 errors)
    logger.info("Downloading %s -> %s", url, current_path)
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "goatools-downloader/1.0"},
    )
    with urllib.request.urlopen(req) as response:
        with open(str(current_path), "wb") as out_f:
            import shutil as _shutil
            _shutil.copyfileobj(response, out_f)
    return current_path


def get_go_ion_channel_genes(cache_dir=None):
    """Get human gene symbols annotated with ion channel GO terms.

    Uses ``goatools`` to parse the Gene Ontology OBO file and the human
    GAF annotation file.  Collects all descendant terms of:

    - GO:0005216 (ion channel activity)
    - GO:0034765 (regulation of ion transmembrane transport)

    and returns gene symbols annotated to any of these terms.

    Parameters
    ----------
    cache_dir : str or Path or None
        Directory for caching downloaded GO files.  Defaults to ``dat/GO/``
        relative to the project root.

    Returns
    -------
    set of str
        Gene symbols annotated with ion channel GO terms.
    """
    from goatools.obo_parser import GODag
    from goatools.anno.gaf_reader import GafReader

    if cache_dir is None:
        # Default: dat/GO/ relative to this source file's project root
        project_root = Path(__file__).resolve().parents[1]
        cache_dir = project_root / "dat" / "GO"

    cache_dir = Path(cache_dir)

    # 1. Download / find OBO and GAF files
    obo_path = _find_or_download(cache_dir, "go-basic", ".obo", _OBO_URL)
    gaf_path = _find_or_download(
        cache_dir, "goa_human", ".gaf.gz", _GAF_URL
    )

    # 2. Parse OBO
    godag = GODag(str(obo_path), prt=None)

    # 3. Collect all descendant GO terms (children of children, recursively)
    target_go_ids = set()
    for root_term in _GO_ION_CHANNEL_TERMS:
        if root_term not in godag:
            warnings.warn(f"GO term {root_term} not found in OBO file")
            continue
        # BFS to collect all descendants (terms that have this as ancestor)
        target_go_ids.add(root_term)
        _collect_descendants(godag, root_term, target_go_ids)

    logger.info(
        "Collected %d GO terms (from %d root terms)",
        len(target_go_ids),
        len(_GO_ION_CHANNEL_TERMS),
    )

    # 4. Parse GAF — need to decompress .gaf.gz first for goatools
    #    GafReader expects a plain text file, not gzipped
    import gzip as _gzip
    import shutil
    import tempfile

    # Check if we need to decompress
    gaf_plain = cache_dir / gaf_path.name.replace(".gz", "")
    if not gaf_plain.exists():
        logger.info("Decompressing %s", gaf_path)
        with _gzip.open(str(gaf_path), "rb") as f_in:
            with open(str(gaf_plain), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    # Read all namespaces (MF for ion channel activity, BP for regulation)
    gaf_reader = GafReader(str(gaf_plain), godag=godag, prt=None,
                           namespaces={"MF", "BP", "CC"})
    associations = gaf_reader.get_associations()

    # 5. Find gene symbols annotated to any target GO term
    gene_symbols = set()
    for nt in associations:
        if nt.GO_ID in target_go_ids:
            gene_symbols.add(nt.DB_Symbol)

    logger.info(
        "Found %d gene symbols annotated with ion channel GO terms",
        len(gene_symbols),
    )

    return gene_symbols


def _collect_descendants(godag, go_id, collected):
    """Recursively collect all descendant GO terms via 'is_a' and 'part_of'.

    Traverses the GO DAG downward from *go_id*, adding all child terms
    to *collected*.

    Parameters
    ----------
    godag : GODag
        Parsed Gene Ontology DAG.
    go_id : str
        Root GO term ID (e.g. "GO:0005216").
    collected : set
        Accumulator set; modified in place.
    """
    term = godag[go_id]
    for child in term.children:
        if child.id not in collected:
            collected.add(child.id)
            _collect_descendants(godag, child.id, collected)


# ---------------------------------------------------------------------------
# PPI bridge utilities (STRING network)
# ---------------------------------------------------------------------------
import gzip
import networkx as nx


def load_string_network(info_path, links_path, score_threshold=400):
    """Load STRING PPI network as a networkx graph.

    Parameters
    ----------
    info_path : str or Path
        Path to gzipped STRING protein info file
        (tab-separated: string_protein_id, preferred_name, ...).
    links_path : str or Path
        Path to gzipped STRING protein links file
        (space-separated: protein1, protein2, combined_score).
    score_threshold : int
        Minimum combined_score to include an edge (default 400).

    Returns
    -------
    G : nx.Graph
        Undirected graph with STRING protein IDs as nodes.
    name_to_id : dict
        ``{gene_symbol: string_protein_id}``
    id_to_name : dict
        ``{string_protein_id: gene_symbol}``
    """
    # --- 1. Parse info file to build name <-> id mappings ---
    name_to_id = {}
    id_to_name = {}
    with gzip.open(info_path, "rt") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            string_id = parts[0]
            preferred_name = parts[1]
            name_to_id[preferred_name] = string_id
            id_to_name[string_id] = preferred_name

    # --- 2. Parse links file and build graph ---
    G = nx.Graph()
    # Add all protein nodes
    for string_id in id_to_name:
        G.add_node(string_id)

    with gzip.open(links_path, "rt") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            p1 = parts[0]
            p2 = parts[1]
            score = int(parts[2])
            if score >= score_threshold:
                G.add_edge(p1, p2)

    return G, name_to_id, id_to_name


def compute_2hop_reachability(G, node):
    """Fraction of graph nodes reachable from *node* within 2 hops.

    Parameters
    ----------
    G : nx.Graph
        Undirected graph.
    node
        Node identifier to start from.

    Returns
    -------
    float
        ``len(reachable) / (G.number_of_nodes() - 1)``, or 0.0 if *node*
        is not in the graph.
    """
    if node not in G:
        return 0.0

    total = G.number_of_nodes() - 1
    if total <= 0:
        return 0.0

    # 1-hop neighbors
    neighbors_1 = set(G.neighbors(node))
    # 2-hop: neighbors of neighbors
    neighbors_2 = set()
    for n1 in neighbors_1:
        neighbors_2.update(G.neighbors(n1))
    # Union of 1-hop and 2-hop, excluding the node itself
    reachable = (neighbors_1 | neighbors_2) - {node}

    return len(reachable) / total


def find_ppi_bridges(G, source_ids, target_ids, id_to_name, max_hops=2):
    """Find bridge proteins on shortest paths between source and target genes.

    Parameters
    ----------
    G : nx.Graph
        PPI network graph.
    source_ids : dict
        ``{symbol: node_id}`` for source genes.
    target_ids : dict
        ``{symbol: node_id}`` for target genes.
    id_to_name : dict
        ``{node_id: gene_symbol}`` for translating node IDs to names.
    max_hops : int
        Maximum path length (number of edges) to consider (default 2).

    Returns
    -------
    list of dict
        Each dict has keys: ``source``, ``target``, ``bridge``,
        ``path_length``, ``path_names``.

        - ``bridge`` is ``None`` for direct interactions (path_length == 1)
        - ``bridge`` is the gene symbol of the intermediate node for
          path_length == 2
    """
    results = []
    for src_sym, src_id in source_ids.items():
        if src_id not in G:
            continue
        for tgt_sym, tgt_id in target_ids.items():
            if tgt_id not in G:
                continue
            try:
                path = nx.shortest_path(G, src_id, tgt_id)
            except nx.NetworkXNoPath:
                continue

            path_length = len(path) - 1  # number of edges
            if path_length > max_hops:
                continue

            path_names = [id_to_name.get(n, n) for n in path]

            if path_length == 1:
                bridge = None
            elif path_length == 2:
                bridge = id_to_name.get(path[1], path[1])
            else:
                bridge = None  # shouldn't happen with max_hops=2

            results.append({
                "source": src_sym,
                "target": tgt_sym,
                "bridge": bridge,
                "path_length": path_length,
                "path_names": path_names,
            })

    return results


# ---------------------------------------------------------------------------
# Ephys feature extraction and comparison
# ---------------------------------------------------------------------------

# EphysSumStats pipeline has a 1e6x factor error in dV/dt:
# spike_detection_new.py line 391 uses *1000 instead of /1000
DVDT_CORRECTION = 1e6


def extract_all_ephys_features(bundle_dir):
    """Extract per-cell electrophysiology summary from all analysis bundles.

    Iterates over subdirectories in *bundle_dir*, reads each cell's
    ``analysis.csv``, identifies rheobase and hero sweeps, and extracts
    a standard feature set.

    Parameters
    ----------
    bundle_dir : str or Path
        Path to directory containing one subdirectory per cell, each with
        an ``analysis.csv`` file.

    Returns
    -------
    pd.DataFrame
        One row per cell with columns: cell_id, genotype, and all extracted
        electrophysiology features.
    """
    from pathlib import Path
    bundle_dir = Path(bundle_dir)

    records = []
    for cell_dir in sorted(bundle_dir.iterdir()):
        if not cell_dir.is_dir():
            continue
        analysis_file = cell_dir / "analysis.csv"
        if not analysis_file.exists():
            continue

        cell_id = cell_dir.name
        genotype = "WT" if cell_id.startswith("WT") else "Df16"

        df = pd.read_csv(analysis_file)

        # Filter to spiking sweeps
        spiking = df[df["spike_frequency_Hz"] > 0].copy()
        if len(spiking) == 0:
            # No spiking sweeps -- skip this cell
            continue

        # Rheobase sweep = first spiking sweep (lowest injected current)
        rheo_idx = spiking["avg_injected_current_pA"].idxmin()
        rheo = spiking.loc[rheo_idx]

        # Hero sweep = sweep with max spike_frequency_Hz
        hero_idx = spiking["spike_frequency_Hz"].idxmax()
        hero = spiking.loc[hero_idx]

        # Max current among all spiking sweeps
        max_current = spiking["avg_injected_current_pA"].max()
        rheo_current = rheo["avg_injected_current_pA"]

        rec = {
            "cell_id": cell_id,
            "genotype": genotype,
            "rheobase_pA": rheo_current,
            "RMP_mV": df["resting_vm_mean_mV"].mean(),
            "rise_slope": rheo["avg_upstroke_mVms"] / DVDT_CORRECTION,
            "decay_slope": rheo["avg_downstroke_mVms"] / DVDT_CORRECTION,
            "ap_width": rheo["avg_ap_width_ms"],
            "peak_voltage": rheo["avg_peak_voltage_mV"],
            "threshold_mV": rheo["avg_threshold_voltage_mV"],
            "amplitude_mV": rheo["avg_threshold_to_peak_mV"],
            "ud_ratio": rheo["avg_upstroke_downstroke_ratio"],
            "hero_freq_Hz": hero["spike_frequency_Hz"],
            "hero_mean_isi_ms": hero["mean_isi_ms"],
            "hero_cv_isi": hero["cv_isi"],
            "max_freq_Hz": spiking["spike_frequency_Hz"].max(),
            "adapt_width_ratio": hero["last_first_AP_width_ratio"],
            "adapt_peak_thresh_ratio": hero["last_first_AP_peak_threshold_ratio"],
            "adapt_trough_ratio": hero["last_first_AP_fast_trough_ratio"],
        }

        # F-I slope: max_freq / (max_current - rheobase_current)
        denom = max_current - rheo_current
        rec["fi_slope"] = rec["max_freq_Hz"] / denom if denom > 0 else np.nan

        records.append(rec)

    return pd.DataFrame(records)


def compare_feature(ephys_df, feature_col):
    """Compare an electrophysiology feature between WT and Df16 genotypes.

    Computes a two-sample t-test, Cohen's d, and a 95% confidence interval
    on the mean difference.

    Parameters
    ----------
    ephys_df : pd.DataFrame
        Must contain a ``"genotype"`` column (values ``"WT"`` or ``"Df16"``)
        and the column specified by *feature_col*.
    feature_col : str
        Name of the numeric feature column to compare.

    Returns
    -------
    dict or None
        Dictionary with keys: wt_mean, wt_sd, n_wt, df16_mean, df16_sd,
        n_df16, t, p, cohens_d, ci_low, ci_high.
        Returns ``None`` if either group has fewer than 3 non-null values.
    """
    wt = ephys_df.loc[ephys_df["genotype"] == "WT", feature_col].dropna()
    df16 = ephys_df.loc[ephys_df["genotype"] == "Df16", feature_col].dropna()

    if len(wt) < 3 or len(df16) < 3:
        return None

    wt_mean = wt.mean()
    df16_mean = df16.mean()
    wt_sd = wt.std(ddof=1)
    df16_sd = df16.std(ddof=1)
    n_wt = len(wt)
    n_df16 = len(df16)

    # Two-sample t-test (Welch's by default)
    t_stat, p_val = ttest_ind(wt, df16, equal_var=False)

    # Cohen's d with pooled SD
    pooled_sd = np.sqrt(
        ((n_wt - 1) * wt_sd**2 + (n_df16 - 1) * df16_sd**2)
        / (n_wt + n_df16 - 2)
    )
    cohens_d = (wt_mean - df16_mean) / pooled_sd if pooled_sd > 0 else np.nan

    # 95% CI on the difference
    diff = wt_mean - df16_mean
    se_diff = np.sqrt(wt_sd**2 / n_wt + df16_sd**2 / n_df16)
    ci_low = diff - 1.96 * se_diff
    ci_high = diff + 1.96 * se_diff

    return {
        "wt_mean": wt_mean,
        "wt_sd": wt_sd,
        "n_wt": n_wt,
        "df16_mean": df16_mean,
        "df16_sd": df16_sd,
        "n_df16": n_df16,
        "t": t_stat,
        "p": p_val,
        "cohens_d": cohens_d,
        "ci_low": ci_low,
        "ci_high": ci_high,
    }
