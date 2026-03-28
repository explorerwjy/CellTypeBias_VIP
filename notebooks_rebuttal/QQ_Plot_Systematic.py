# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %% [markdown]
# # QQ Plots: Observed vs Expected P-values Across Null Methods
#
# Supplementary figure comparing p-value distributions under Random vs
# Best-of-N (matched) null models across all disorders. A well-calibrated
# null should produce p-values that follow the diagonal; enrichment of
# small p-values above the diagonal indicates real biological signal that
# survives confounder matching.

# %% [markdown]
# ## Setup

# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import OrderedDict

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

FIG_DIR = PROJ_DIR / "results" / "figures" / "systematic_comparison"
FIG_DIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['font.size'] = 14

# %% [markdown]
# ## Configuration
#
# Toggle `SHOW_METHODS` to control which null methods appear on the plot.

# %%
# All available methods and their result directories
ALL_METHODS = OrderedDict([
    ("Random", "main_results/random"),
    ("Best-of-N", "main_results/matched_WB_mean_phastCons_n_CDS_bases_Best1000"),
])

# --- TOGGLE HERE: which methods to show ---
# Default: Random + Best-of-N for the main supplementary figure
SHOW_METHODS = ["Random", "Best-of-N"]
# Uncomment below to show all methods:
# SHOW_METHODS = list(ALL_METHODS.keys())

# Gene sets (disorders) to plot
GENE_SETS = OrderedDict([
    ("SCZ", "SCZ"),
    ("ASD (all)", "ASD_All"),
    ("ASD w/o ID", "ASD_HIQ"),
    ("ASD with ID", "ASD_LIQ"),
    ("DDD", "DDD_61"),
    ("22q11.2 del", "22q_del"),
    ("UKBB VNR+", "UKBB_VNR_Pos"),
    ("UKBB VNR-", "UKBB_VNR_Neg"),
])

# Negative controls (non-brain traits, expect no cell-type enrichment)
NEG_CTRL_SETS = OrderedDict([
    ("HDL cholesterol", "NegCtrl_HDL"),
    ("Alanine", "NegCtrl_Alanine"),
    ("Red Blood Cell", "NegCtrl_RBC"),
    ("IBD", "NegCtrl_IBD"),
])

ANALYSIS = "Centering"

# Style
METHOD_STYLE = {
    "Random":       {"color": "#aaaaaa", "marker": "o", "zorder": 2},
    "Gene-by-gene": {"color": "#e67e22", "marker": "s", "zorder": 3},
    "Rejection":    {"color": "#f1c40f", "marker": "D", "zorder": 3},
    "Best-of-N":    {"color": "#2563eb", "marker": "o", "zorder": 4},
    "PropWeight":   {"color": "#16a34a", "marker": "^", "zorder": 3},
    "SIS":          {"color": "#9b59b6", "marker": "v", "zorder": 3},
}

# %% [markdown]
# ## Load P-values

# %%
ALL_GENE_SETS = OrderedDict(list(GENE_SETS.items()) + list(NEG_CTRL_SETS.items()))

pvals = {}       # {method: {gene_set_label: array of p-values}}
is_neuron = {}   # {method: {gene_set_label: bool array, True = neuronal}}
for method in SHOW_METHODS:
    result_dir = ALL_METHODS[method]
    pvals[method] = {}
    is_neuron[method] = {}
    for gs_label, gs_key in ALL_GENE_SETS.items():
        fpath = PROJ_DIR / "results" / result_dir / ANALYSIS / f"{gs_key}_bias_addP.csv"
        if fpath.exists():
            df = pd.read_csv(fpath, index_col=0)
            pvals[method][gs_label] = df["P-value"].values
            is_neuron[method][gs_label] = (df["Class"] == "NEUR").values
        else:
            print(f"Missing: {method} / {gs_key}")

n_loaded = sum(len(v) for v in pvals.values())
print(f"Loaded {n_loaded} p-value vectors "
      f"({len(SHOW_METHODS)} methods × {len(ALL_GENE_SETS)} gene sets)")

# %% [markdown]
# ## QQ Plot

# %%
def qq_plot(pvals_dict, gene_sets, method_style, figsize=None, title=None,
            neg_ctrl_labels=None, ncols=None, neuron_dict=None):
    """
    Multi-panel QQ plot of observed vs expected -log10(p).

    Parameters
    ----------
    pvals_dict : dict
        {method_name: {gene_set_label: np.array of p-values}}
    gene_sets : OrderedDict
        {display_label: config_key} — determines panel order
    method_style : dict
        {method_name: {"color": ..., "marker": ..., "zorder": ...}}
    neg_ctrl_labels : set or None
        Labels to mark as negative controls (italic title, gray background)
    ncols : int or None
        Number of columns per row (default: all in one row)
    neuron_dict : dict or None
        {method_name: {gene_set_label: bool array}} — True = neuronal.
        When provided, neuronal cells are plotted as circles and
        non-neuronal cells as diamonds.
    """
    n_panels = len(gene_sets)
    if neg_ctrl_labels is None:
        neg_ctrl_labels = set()

    if ncols is None:
        ncols = n_panels
    nrows = int(np.ceil(n_panels / ncols))

    if figsize is None:
        figsize = (5 * ncols, 5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    for idx, gs_label in enumerate(gene_sets):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        global_max = 0
        for method, gs_pvals in pvals_dict.items():
            if gs_label not in gs_pvals:
                continue
            p = gs_pvals[gs_label]
            sort_idx = np.argsort(p)
            obs_p = np.maximum(p[sort_idx], 1e-10)
            n = len(obs_p)
            expected = (np.arange(1, n + 1)) / (n + 1)

            neg_log_exp = -np.log10(expected)
            neg_log_obs = -np.log10(obs_p)

            global_max = max(global_max, neg_log_exp.max(), neg_log_obs.max())

            style = method_style.get(method, {})
            base_color = style.get("color", "grey")
            base_zorder = style.get("zorder", 2)

            if neuron_dict is not None and method in neuron_dict and gs_label in neuron_dict[method]:
                neur_mask = neuron_dict[method][gs_label][sort_idx]

                # Neuronal cells: circles
                ax.scatter(
                    neg_log_exp[neur_mask], neg_log_obs[neur_mask],
                    s=18, alpha=0.7, linewidths=0.3, edgecolors="white",
                    color=base_color, marker="o", zorder=base_zorder,
                    label=f"{method} (neuron)",
                )
                # Non-neuronal cells: diamonds
                ax.scatter(
                    neg_log_exp[~neur_mask], neg_log_obs[~neur_mask],
                    s=28, alpha=0.85, linewidths=0.4, edgecolors="white",
                    color=base_color, marker="D", zorder=base_zorder + 1,
                    label=f"{method} (non-neuron)",
                )
            else:
                ax.scatter(
                    neg_log_exp, neg_log_obs,
                    s=18, alpha=0.7, linewidths=0.3, edgecolors="white",
                    color=base_color,
                    marker=style.get("marker", "o"),
                    zorder=base_zorder,
                    label=method,
                )

        # Diagonal and axis limits
        upper = global_max * 1.08
        ax.plot([0, upper], [0, upper], ls="--", lw=1, color="black", zorder=1)
        ax.set_xlim(0, upper)
        ax.set_ylim(0, upper)
        ax.set_aspect("equal")

        ax.set_xlabel("Expected $-\\log_{10}(p)$")
        if col == 0:
            ax.set_ylabel("Observed $-\\log_{10}(p)$")

        # Style negative controls differently
        if gs_label in neg_ctrl_labels:
            ax.set_title(gs_label, fontweight="bold", fontstyle="italic",
                         color="#666666", pad=8)
            ax.set_facecolor("#f5f5f5")
        else:
            ax.set_title(gs_label, fontweight="bold", pad=8)

        if idx == n_panels - 1:
            ax.legend(fontsize=10, framealpha=0.8, loc="lower right")

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].set_visible(False)

    if title:
        fig.suptitle(title, fontsize=16, fontweight="bold", y=1.02)

    fig.subplots_adjust(hspace=0.35, wspace=0.30)
    return fig


# Main figure: disorders only
fig = qq_plot(pvals, GENE_SETS, METHOD_STYLE, neuron_dict=is_neuron,
              title="Cell-Type Bias P-value QQ Plots")

fig.savefig(FIG_DIR / "FigS_QQ_plot.pdf", dpi=300, transparent=True,
            bbox_inches="tight")
fig.savefig(FIG_DIR / "FigS_QQ_plot.png", dpi=300, transparent=True,
            bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR / 'FigS_QQ_plot.pdf'}")

# %% [markdown]
# ## QQ Plot — With Negative Controls
#
# Disorders (5 panels) + negative controls (2 panels, italic titles, gray
# background). Negative controls should hug the diagonal, confirming that
# non-brain traits produce no cell-type enrichment.

# %%
fig_nc = qq_plot(pvals, ALL_GENE_SETS, METHOD_STYLE,
                 neg_ctrl_labels=set(NEG_CTRL_SETS.keys()),
                 ncols=4, neuron_dict=is_neuron,
                 title="Cell-Type Bias QQ Plots (with negative controls)")

fig_nc.savefig(FIG_DIR / "FigS_QQ_plot_with_negctrl.pdf", dpi=300,
               transparent=True, bbox_inches="tight")
fig_nc.savefig(FIG_DIR / "FigS_QQ_plot_with_negctrl.png", dpi=300,
               transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR / 'FigS_QQ_plot_with_negctrl.pdf'}")

# %% [markdown]
# ## QQ Plot — All Methods (optional)
#
# Uncomment and run this cell to generate the version with all 6 methods.

# %%
# Reload with all methods and all gene sets (including neg controls)
pvals_all = {}
is_neuron_all = {}
for method in ALL_METHODS:
    result_dir = ALL_METHODS[method]
    pvals_all[method] = {}
    is_neuron_all[method] = {}
    for gs_label, gs_key in ALL_GENE_SETS.items():
        fpath = PROJ_DIR / "results" / result_dir / ANALYSIS / f"{gs_key}_bias_addP.csv"
        if fpath.exists():
            df = pd.read_csv(fpath, index_col=0)
            pvals_all[method][gs_label] = df["P-value"].values
            is_neuron_all[method][gs_label] = (df["Class"] == "NEUR").values

fig_all = qq_plot(pvals_all, ALL_GENE_SETS, METHOD_STYLE,
                  neg_ctrl_labels=set(NEG_CTRL_SETS.keys()),
                  ncols=4, neuron_dict=is_neuron_all,
                  title="Cell-Type Bias QQ Plots — All Null Methods")

fig_all.savefig(FIG_DIR / "FigS_QQ_plot_all_methods.pdf", dpi=300,
                transparent=True, bbox_inches="tight")
fig_all.savefig(FIG_DIR / "FigS_QQ_plot_all_methods.png", dpi=300,
                transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR / 'FigS_QQ_plot_all_methods.pdf'}")

# %% [markdown]
# ## Genomic Inflation Factor ($\lambda$)

# %%
from scipy import stats

def genomic_inflation_factor(pvalues):
    """Compute genomic inflation factor lambda."""
    chi2 = stats.chi2.ppf(1 - pvalues, df=1)
    return np.median(chi2) / stats.chi2.ppf(0.5, df=1)

rows = []
for method in ALL_METHODS:
    result_dir = ALL_METHODS[method]
    for gs_label, gs_key in ALL_GENE_SETS.items():
        fpath = PROJ_DIR / "results" / result_dir / ANALYSIS / f"{gs_key}_bias_addP.csv"
        if fpath.exists():
            df = pd.read_csv(fpath, index_col=0)
            p = df["P-value"].values
            p = p[(p > 0) & (p < 1)]  # exclude exact 0 or 1
            lam = genomic_inflation_factor(p)
            rows.append({"Method": method, "Gene Set": gs_label,
                          "lambda": round(lam, 3)})

lambda_df = pd.DataFrame(rows).pivot(index="Gene Set", columns="Method", values="lambda")
# Reorder columns
lambda_df = lambda_df[[m for m in ALL_METHODS if m in lambda_df.columns]]
print("Genomic inflation factor (lambda):\n")
print(lambda_df.to_string())
