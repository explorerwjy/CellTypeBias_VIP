"""
Convergence analysis utilities: gene set constants, mapping functions,
and partial-correlation helpers.

Defines curated ion channel / scaffolding gene families, 22q11.2 deletion
genes, ephys feature metadata, and statistical functions for computing
Pearson partial correlations controlling for cell-class membership.
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, norm
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests
from statsmodels.tools import add_constant

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
      - Tier 1: Spearman rho on raw specificity profiles
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

            # Tier 1: Spearman on raw specificity
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
