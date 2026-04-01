"""
Convergence analysis utilities: gene set constants, mapping functions,
and partial-correlation helpers.

Defines curated ion channel / scaffolding gene families, 22q11.2 deletion
genes, ephys feature metadata, and statistical functions for computing
Pearson partial correlations controlling for cell-class membership.
"""

import pandas as pd
from scipy.stats import pearsonr
from statsmodels.regression.linear_model import OLS
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
