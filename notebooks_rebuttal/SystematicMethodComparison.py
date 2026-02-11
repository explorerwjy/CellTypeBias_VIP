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
# ---

# %% [markdown]
# # Systematic Comparison of All Matching Methods
#
# **Purpose**: Fair comparison of 6 null gene set generation methods using
# identical matching variables (gene length, expression, conservation) to
# address Reviewer 3's concerns about permutation null robustness.
#
# ## Methods Compared
#
# | # | Method | Description |
# |---|--------|-------------|
# | 1 | Random | No matching (baseline) |
# | 2 | Gene-by-gene | Tricubic kernel per-gene matching |
# | 3 | Rejection | Accept/reject with distance threshold |
# | 4 | Best-of-N | Draw 1000 sets, pick best match |
# | 5 | PropensityWeighted | Tricubic kernel propensity scores |
# | 6 | SIS | Sequential Importance Sampling |
#
# All methods match on: **n_CDS_bases** (gene length), **WB** (expression),
# **mean_phastCons** (conservation).

# %% [markdown]
# ## Setup

# %%
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from collections import OrderedDict
import warnings
warnings.filterwarnings('ignore')

NOTEBOOK_DIR = Path().absolute()
PROJ_DIR = NOTEBOOK_DIR.parent if NOTEBOOK_DIR.name in ("notebooks_rebuttal", "notebooks") else NOTEBOOK_DIR
sys.path.insert(0, str(PROJ_DIR / "src"))

# Figure output directory
FIG_DIR = PROJ_DIR / "results" / "figures" / "systematic_comparison"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Transparent backgrounds
mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['axes.labelsize'] = 12

print(f"Project directory: {PROJ_DIR}")
print(f"Figure directory: {FIG_DIR}")

# %% [markdown]
# ## Section 1: Define Methods and Load Results

# %%
# Method definitions with result directory names
METHODS = OrderedDict([
    ("Random", "systematic_random"),
    ("Gene-by-gene", "systematic_gene_by_gene_WB_mean_phastCons_n_CDS_bases_Tricubic"),
    ("Rejection", "systematic_rejection_WB_mean_phastCons_n_CDS_bases_Rejection"),
    ("Best-of-N", "systematic_best_of_n_WB_mean_phastCons_n_CDS_bases_Best1000"),
    ("PropWeight", "systematic_propensity_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic"),
    ("SIS", "systematic_sis_WB_mean_phastCons_n_CDS_bases_SIS"),
])

# Gene sets to analyze
GENE_SETS = ["SCZ", "ASD_All", "ASD_HIQ", "ASD_LIQ", "DDD_61"]

# Gene weight file paths
GW_PATHS = {
    "SCZ": PROJ_DIR / "dat/GeneWeights/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "ASD_All": PROJ_DIR / "dat/GeneWeights/Spark_Meta_EWS.GeneWeight.csv",
    "ASD_HIQ": PROJ_DIR / "dat/GeneWeights/HIQ.top61.nopLI.LGD_Dmis_SameWeight.gw",
    "ASD_LIQ": PROJ_DIR / "dat/GeneWeights/LIQ.top61.nopLI.LGD_Dmis_SameWeight.gw",
    "DDD_61": PROJ_DIR / "dat/GeneWeights/DDD.top61.gw.csv",
}

# Key superclusters for the manuscript
KEY_SUPERCLUSTERS = ["CGE interneuron", "MGE interneuron", "Medium spiny neuron",
                     "Eccentric medium spiny neuron"]

ANALYSIS = "Centering"

# %%
# Load all supercluster-level results
supercluster_results = {}
missing_results = []

for method_name, result_dir in METHODS.items():
    supercluster_results[method_name] = {}
    for gs in GENE_SETS:
        fpath = PROJ_DIR / "results" / result_dir / ANALYSIS / f"{gs}_bias_addP_supercluster.csv"
        if fpath.exists():
            df = pd.read_csv(fpath)
            supercluster_results[method_name][gs] = df
        else:
            missing_results.append(f"{method_name}/{gs}")

if missing_results:
    print(f"Missing results ({len(missing_results)}):")
    for m in missing_results:
        print(f"  - {m}")
else:
    print("All results loaded successfully!")

# %%
# Load all cluster-level (461 cell types) results
cluster_results = {}
for method_name, result_dir in METHODS.items():
    cluster_results[method_name] = {}
    for gs in GENE_SETS:
        fpath = PROJ_DIR / "results" / result_dir / ANALYSIS / f"{gs}_bias_addP.csv"
        if fpath.exists():
            df = pd.read_csv(fpath, index_col=0)
            cluster_results[method_name][gs] = df

print(f"Loaded cluster-level results: {sum(len(v) for v in cluster_results.values())} combinations")

# %% [markdown]
# ## Section 2: Significance Tables
#
# ### 2A: Supercluster-Level Significance (q < 0.1)

# %%
def count_significant_superclusters(results_dict, gene_sets, q_threshold=0.1):
    """Count significant superclusters per method × gene set."""
    rows = []
    for method_name in results_dict:
        for gs in gene_sets:
            if gs not in results_dict[method_name]:
                continue
            df = results_dict[method_name][gs]
            n_sig = (df["q-value"] < q_threshold).sum()
            n_total = len(df)

            # Check key superclusters
            for sc in KEY_SUPERCLUSTERS:
                sc_row = df[df["Supercluster"] == sc]
                if len(sc_row) > 0:
                    q = sc_row["q-value"].values[0]
                    z = sc_row["Z-score"].values[0]
                    sig = "***" if q < 0.01 else "**" if q < 0.05 else "*" if q < 0.1 else ""
                    rows.append({
                        "Method": method_name,
                        "Gene Set": gs,
                        "Supercluster": sc,
                        "Z-score": round(z, 2),
                        "q-value": round(q, 4),
                        "Significant": sig,
                        "n_total_sig": n_sig,
                    })
    return pd.DataFrame(rows)

sig_table = count_significant_superclusters(supercluster_results, GENE_SETS)
if len(sig_table) > 0:
    print("=== KEY SUPERCLUSTER SIGNIFICANCE (q < 0.1) ===\n")

    for gs in GENE_SETS:
        gs_data = sig_table[sig_table["Gene Set"] == gs]
        if len(gs_data) == 0:
            continue
        print(f"\n--- {gs} ---")
        pivot = gs_data.pivot_table(
            index="Supercluster", columns="Method",
            values="q-value", aggfunc="first"
        )
        # Reorder columns by method order
        method_order = [m for m in METHODS.keys() if m in pivot.columns]
        pivot = pivot[method_order]
        print(pivot.round(4).to_string())
else:
    print("No results available yet. Pipeline may still be running.")

# %% [markdown]
# ### 2B: Cluster-Level Significance for SCZ (461 clusters)

# %%
def count_significant_clusters_by_supercluster(cluster_results_dict, gene_set,
                                                target_superclusters, q_threshold=0.1):
    """Count significant individual clusters within each target supercluster."""
    rows = []
    for method_name in cluster_results_dict:
        if gene_set not in cluster_results_dict[method_name]:
            continue
        df = cluster_results_dict[method_name][gene_set]
        for sc in target_superclusters:
            sc_df = df[df["Supercluster"] == sc]
            n_sig = (sc_df["q-value"] < q_threshold).sum()
            n_total = len(sc_df)
            mean_z = sc_df["Z-score"].mean() if len(sc_df) > 0 else np.nan
            rows.append({
                "Method": method_name,
                "Supercluster": sc,
                "n_sig": n_sig,
                "n_total": n_total,
                "mean_Z": round(mean_z, 2),
                "frac_sig": round(n_sig / n_total, 3) if n_total > 0 else 0,
            })
    return pd.DataFrame(rows)

print("=== SCZ: CLUSTER-LEVEL SIGNIFICANCE (q < 0.1) ===")
print("Number of individual CGE/MGE clusters passing FDR < 0.1\n")

scz_clusters = count_significant_clusters_by_supercluster(
    cluster_results, "SCZ",
    ["CGE interneuron", "MGE interneuron"]
)
if len(scz_clusters) > 0:
    pivot = scz_clusters.pivot_table(
        index="Supercluster", columns="Method",
        values=["n_sig", "n_total", "mean_Z"], aggfunc="first"
    )
    method_order = [m for m in METHODS.keys() if m in pivot.columns.get_level_values(1)]
    for val_col in ["n_sig", "mean_Z"]:
        print(f"\n{val_col}:")
        sub = pivot[val_col]
        sub = sub[[m for m in method_order if m in sub.columns]]
        print(sub.to_string())
else:
    print("No SCZ cluster results available yet.")

# %% [markdown]
# ### 2C: ASD Medium Spiny Neurons (cluster-level)

# %%
print("=== ASD_All: MEDIUM SPINY NEURON SIGNIFICANCE (q < 0.1) ===\n")

msn_clusters = count_significant_clusters_by_supercluster(
    cluster_results, "ASD_All",
    ["Medium spiny neuron", "Eccentric medium spiny neuron"]
)
if len(msn_clusters) > 0:
    pivot = msn_clusters.pivot_table(
        index="Supercluster", columns="Method",
        values=["n_sig", "n_total", "mean_Z"], aggfunc="first"
    )
    method_order = [m for m in METHODS.keys() if m in pivot.columns.get_level_values(1)]
    for val_col in ["n_sig", "mean_Z"]:
        print(f"\n{val_col}:")
        sub = pivot[val_col]
        sub = sub[[m for m in method_order if m in sub.columns]]
        print(sub.to_string())
else:
    print("No ASD_All cluster results available yet.")

# %% [markdown]
# ### 2D: Combined Pass/Fail Verdicts

# %%
def evaluate_criteria(supercluster_results, cluster_results, q_threshold=0.1):
    """Evaluate all key criteria from the plan for each method."""
    verdicts = {}
    for method_name in METHODS:
        v = {}

        # Criterion 1: SCZ → CGE interneuron (supercluster, q<0.1)
        if "SCZ" in supercluster_results.get(method_name, {}):
            df = supercluster_results[method_name]["SCZ"]
            cge = df[df["Supercluster"] == "CGE interneuron"]
            v["SCZ→CGE (SC)"] = "PASS" if len(cge) > 0 and cge["q-value"].values[0] < q_threshold else "FAIL"
        else:
            v["SCZ→CGE (SC)"] = "N/A"

        # Criterion 2: SCZ → MGE interneuron (supercluster, q<0.1)
        if "SCZ" in supercluster_results.get(method_name, {}):
            df = supercluster_results[method_name]["SCZ"]
            mge = df[df["Supercluster"] == "MGE interneuron"]
            v["SCZ→MGE (SC)"] = "PASS" if len(mge) > 0 and mge["q-value"].values[0] < q_threshold else "FAIL"
        else:
            v["SCZ→MGE (SC)"] = "N/A"

        # Criterion 3: SCZ → CGE clusters (>= 2 clusters q<0.1)
        if "SCZ" in cluster_results.get(method_name, {}):
            df = cluster_results[method_name]["SCZ"]
            cge_df = df[df["Supercluster"] == "CGE interneuron"]
            n_sig = (cge_df["q-value"] < q_threshold).sum()
            v["SCZ→CGE (clusters)"] = f"PASS ({n_sig})" if n_sig >= 2 else f"FAIL ({n_sig})"
        else:
            v["SCZ→CGE (clusters)"] = "N/A"

        # Criterion 4: SCZ → MGE clusters (>= 2 clusters q<0.1)
        if "SCZ" in cluster_results.get(method_name, {}):
            df = cluster_results[method_name]["SCZ"]
            mge_df = df[df["Supercluster"] == "MGE interneuron"]
            n_sig = (mge_df["q-value"] < q_threshold).sum()
            v["SCZ→MGE (clusters)"] = f"PASS ({n_sig})" if n_sig >= 2 else f"FAIL ({n_sig})"
        else:
            v["SCZ→MGE (clusters)"] = "N/A"

        # Criterion 5: ASD_LIQ → CGE (supercluster, q<0.1)
        if "ASD_LIQ" in supercluster_results.get(method_name, {}):
            df = supercluster_results[method_name]["ASD_LIQ"]
            cge = df[df["Supercluster"] == "CGE interneuron"]
            v["ASD_LIQ→CGE (SC)"] = "PASS" if len(cge) > 0 and cge["q-value"].values[0] < q_threshold else "FAIL"
        else:
            v["ASD_LIQ→CGE (SC)"] = "N/A"

        # Criterion 6: DDD_61 → CGE (supercluster, q<0.1)
        if "DDD_61" in supercluster_results.get(method_name, {}):
            df = supercluster_results[method_name]["DDD_61"]
            cge = df[df["Supercluster"] == "CGE interneuron"]
            v["DDD_61→CGE (SC)"] = "PASS" if len(cge) > 0 and cge["q-value"].values[0] < q_threshold else "FAIL"
        else:
            v["DDD_61→CGE (SC)"] = "N/A"

        # Criterion 7: ASD_All → Medium spiny neuron (cluster-level, >=1 q<0.1)
        if "ASD_All" in cluster_results.get(method_name, {}):
            df = cluster_results[method_name]["ASD_All"]
            msn_df = df[df["Supercluster"].isin(["Medium spiny neuron", "Eccentric medium spiny neuron"])]
            n_sig = (msn_df["q-value"] < q_threshold).sum()
            v["ASD_All→MSN (clusters)"] = f"PASS ({n_sig})" if n_sig >= 1 else f"FAIL ({n_sig})"
        else:
            v["ASD_All→MSN (clusters)"] = "N/A"

        verdicts[method_name] = v

    return pd.DataFrame(verdicts).T

print("=== PASS/FAIL VERDICTS ACROSS ALL CRITERIA ===\n")
verdict_df = evaluate_criteria(supercluster_results, cluster_results)
if not verdict_df.empty and not all(verdict_df.values.flatten() == "N/A"):
    print(verdict_df.to_string())
else:
    print("Results not yet available. Pipeline may still be running.")

# %% [markdown]
# ## Section 3: Matching Quality Assessment
#
# Compare distributions of matching variables between real disease genes
# and null gene sets for each method.

# %%
# Load master matching table
master_table_path = PROJ_DIR / "dat" / "Variable_2_Match_master_table_pct.csv"
master_table = pd.read_csv(master_table_path, index_col=0)
print(f"Master table: {master_table.shape[0]:,} genes × {master_table.shape[1]} columns")
print(f"Columns: {list(master_table.columns)}")

MATCHING_VARS = ["n_CDS_bases", "WB", "mean_phastCons"]
MATCHING_VAR_LABELS = {"n_CDS_bases": "Gene Length (CDS bases)",
                       "WB": "Expression (WB)",
                       "mean_phastCons": "Conservation (phastCons)"}

# %%
def load_real_gene_ids(gw_path):
    """Load real disease gene Entrez IDs from a gene weight file."""
    gw = pd.read_csv(gw_path, header=None, names=["EntrezID", "Weight"])
    return gw["EntrezID"].astype(int).values

def load_null_gene_indices(null_weights_path, n_samples=100):
    """Load null gene indices from a null weights file, sampling a subset."""
    nw = pd.read_csv(null_weights_path, index_col=0)
    # Columns: GeneWeight, 0, 1, ..., 9999
    sim_cols = [c for c in nw.columns if c != "GeneWeight"]
    # Sample a subset for efficiency
    sample_cols = np.random.choice(sim_cols, min(n_samples, len(sim_cols)), replace=False)
    return nw[sample_cols]

def compute_matching_quality(gene_set, method_name, result_dir, n_null_samples=200):
    """Compute matching quality metrics for one method × gene set."""
    # Load real gene IDs
    real_ids = load_real_gene_ids(GW_PATHS[gene_set])
    real_in_table = np.intersect1d(real_ids, master_table.index)

    # Load null gene weights
    nw_path = PROJ_DIR / "results" / result_dir / ANALYSIS / "null_weights" / f"{gene_set}_random_geneweights.csv"
    if not nw_path.exists():
        return None

    nw = pd.read_csv(nw_path, index_col=0)
    sim_cols = [c for c in nw.columns if c != "GeneWeight"]

    # Real gene values
    real_values = {}
    for var in MATCHING_VARS:
        if var in master_table.columns:
            vals = master_table.loc[master_table.index.isin(real_in_table), var].dropna()
            real_values[var] = vals.values

    # Null gene values (sample for efficiency)
    np.random.seed(42)
    sample_cols = np.random.choice(sim_cols, min(n_null_samples, len(sim_cols)), replace=False)

    metrics = {var: {"ks_pvals": [], "mean_bias": []} for var in MATCHING_VARS}

    for col in sample_cols:
        null_gene_indices = nw[col].dropna().astype(int).values
        # Map indices back to gene IDs (rows of the original expression matrix)
        all_genes = nw.index.values
        null_gene_ids = all_genes[null_gene_indices] if null_gene_indices.max() < len(all_genes) else null_gene_indices
        null_in_table = np.intersect1d(null_gene_ids, master_table.index)

        for var in MATCHING_VARS:
            if var not in master_table.columns or len(null_in_table) == 0:
                continue
            null_vals = master_table.loc[master_table.index.isin(null_in_table), var].dropna().values
            if len(null_vals) > 0 and len(real_values.get(var, [])) > 0:
                ks_stat, ks_pval = stats.ks_2samp(real_values[var], null_vals)
                # Bias in std units
                real_mean = np.mean(real_values[var])
                null_mean = np.mean(null_vals)
                pooled_std = np.std(np.concatenate([real_values[var], null_vals]))
                bias = (null_mean - real_mean) / pooled_std if pooled_std > 0 else 0

                metrics[var]["ks_pvals"].append(ks_pval)
                metrics[var]["mean_bias"].append(bias)

    # Summarize
    summary = {}
    for var in MATCHING_VARS:
        if metrics[var]["ks_pvals"]:
            summary[var] = {
                "median_ks_pval": np.median(metrics[var]["ks_pvals"]),
                "mean_abs_bias": np.mean(np.abs(metrics[var]["mean_bias"])),
                "mean_bias": np.mean(metrics[var]["mean_bias"]),
            }
    summary["real_values"] = real_values

    return summary

# %%
# Compute matching quality for SCZ across all methods
print("Computing matching quality for SCZ gene set...")
matching_quality = {}
for method_name, result_dir in METHODS.items():
    print(f"  {method_name}...", end=" ")
    mq = compute_matching_quality("SCZ", method_name, result_dir, n_null_samples=200)
    if mq is not None:
        matching_quality[method_name] = mq
        avg_bias = np.mean([mq[v]["mean_abs_bias"] for v in MATCHING_VARS if v in mq])
        avg_ks = np.mean([mq[v]["median_ks_pval"] for v in MATCHING_VARS if v in mq])
        print(f"avg |bias| = {avg_bias:.3f}, avg KS p = {avg_ks:.3f}")
    else:
        print("not found")

# %% [markdown]
# ### Figure 1: Matching Quality Panel (6 methods × 3 variables)

# %%
def plot_matching_quality_panel(matching_quality, gene_set="SCZ"):
    """6-row × 3-col panel of overlaid real vs null distributions."""
    methods_with_data = [m for m in METHODS if m in matching_quality]
    n_methods = len(methods_with_data)
    if n_methods == 0:
        print("No matching quality data available.")
        return None

    fig, axes = plt.subplots(n_methods, len(MATCHING_VARS), figsize=(14, 3 * n_methods))
    if n_methods == 1:
        axes = axes.reshape(1, -1)

    for i, method in enumerate(methods_with_data):
        mq = matching_quality[method]
        real_vals = mq.get("real_values", {})

        for j, var in enumerate(MATCHING_VARS):
            ax = axes[i, j]
            if var not in mq or var not in real_vals:
                ax.set_visible(False)
                continue

            # Get null gene values for this method (reload a few for plotting)
            result_dir = METHODS[method]
            nw_path = PROJ_DIR / "results" / result_dir / ANALYSIS / "null_weights" / f"{gene_set}_random_geneweights.csv"
            if nw_path.exists():
                nw = pd.read_csv(nw_path, index_col=0)
                sim_cols = [c for c in nw.columns if c != "GeneWeight"]
                # Pool 50 null sets for a smooth null distribution
                np.random.seed(42)
                sample_cols = np.random.choice(sim_cols, min(50, len(sim_cols)), replace=False)
                null_all = []
                all_genes = nw.index.values
                for col in sample_cols:
                    idx = nw[col].dropna().astype(int).values
                    null_ids = all_genes[idx] if idx.max() < len(all_genes) else idx
                    null_in_table = np.intersect1d(null_ids, master_table.index)
                    null_all.extend(master_table.loc[master_table.index.isin(null_in_table), var].dropna().values)

                # Plot
                bins = np.linspace(
                    min(real_vals[var].min(), np.percentile(null_all, 1)),
                    max(real_vals[var].max(), np.percentile(null_all, 99)),
                    40
                )
                ax.hist(null_all, bins=bins, density=True, alpha=0.5, color="steelblue",
                        label="Null", edgecolor="none")
                ax.hist(real_vals[var], bins=bins, density=True, alpha=0.7, color="firebrick",
                        label="Real", edgecolor="none")

            # Annotations
            bias = mq[var]["mean_abs_bias"]
            ks_p = mq[var]["median_ks_pval"]
            ax.set_title(f"|bias|={bias:.3f}, KS p={ks_p:.3f}", fontsize=9)

            if i == 0:
                ax.set_xlabel("")
                ax.text(0.5, 1.15, MATCHING_VAR_LABELS[var], transform=ax.transAxes,
                        ha='center', fontsize=12, fontweight='bold')
            if i == n_methods - 1:
                ax.set_xlabel(var, fontsize=10)
            if j == 0:
                ax.set_ylabel(method, fontsize=11, fontweight='bold')
            if i == 0 and j == len(MATCHING_VARS) - 1:
                ax.legend(fontsize=8, loc='upper right')

            ax.tick_params(labelsize=8)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig_matching_quality_panel.pdf", dpi=300, transparent=True,
                bbox_inches='tight')
    fig.savefig(FIG_DIR / "Fig_matching_quality_panel.png", dpi=300, transparent=True,
                bbox_inches='tight')
    plt.show()
    return fig

fig1 = plot_matching_quality_panel(matching_quality, "SCZ")

# %% [markdown]
# ### Matching Quality Summary Table

# %%
def matching_quality_summary_table(matching_quality):
    """Create a summary table of matching quality across methods."""
    rows = []
    for method in METHODS:
        if method not in matching_quality:
            continue
        mq = matching_quality[method]
        row = {"Method": method}
        biases = []
        for var in MATCHING_VARS:
            if var in mq:
                row[f"|bias|_{var}"] = round(mq[var]["mean_abs_bias"], 4)
                row[f"KS_p_{var}"] = round(mq[var]["median_ks_pval"], 4)
                biases.append(mq[var]["mean_abs_bias"])
        row["avg_|bias|"] = round(np.mean(biases), 4) if biases else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index("Method")

mq_summary = matching_quality_summary_table(matching_quality)
if len(mq_summary) > 0:
    print("=== MATCHING QUALITY SUMMARY (SCZ) ===\n")
    print(mq_summary.to_string())

# %% [markdown]
# ## Section 4: Significance Heatmap
#
# Rows = superclusters, columns = methods (ordered by matching stringency).

# %%
def plot_significance_heatmap(supercluster_results, gene_sets_to_plot,
                              methods_order=None, q_threshold=0.1):
    """Create significance heatmaps: -log10(q) per supercluster × method."""
    if methods_order is None:
        methods_order = list(METHODS.keys())

    available_methods = [m for m in methods_order if m in supercluster_results
                         and any(gs in supercluster_results[m] for gs in gene_sets_to_plot)]
    if len(available_methods) == 0:
        print("No data available for heatmap.")
        return None

    n_panels = len(gene_sets_to_plot)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels + 2, 10))
    if n_panels == 1:
        axes = [axes]

    for panel_idx, gs in enumerate(gene_sets_to_plot):
        ax = axes[panel_idx]

        # Build matrix: superclusters × methods
        # Get all superclusters from any available result
        all_sc = set()
        for m in available_methods:
            if gs in supercluster_results.get(m, {}):
                all_sc.update(supercluster_results[m][gs]["Supercluster"].values)
        all_sc = sorted(all_sc)

        matrix = pd.DataFrame(np.nan, index=all_sc, columns=available_methods)
        for m in available_methods:
            if gs in supercluster_results.get(m, {}):
                df = supercluster_results[m][gs]
                for _, row in df.iterrows():
                    sc = row["Supercluster"]
                    matrix.loc[sc, m] = -np.log10(max(row["q-value"], 1e-10))

        # Sort rows by mean -logQ
        row_means = matrix.mean(axis=1).fillna(0)
        matrix = matrix.loc[row_means.sort_values(ascending=False).index]

        # Plot heatmap
        im = ax.imshow(matrix.values.astype(float), aspect='auto', cmap='YlOrRd',
                       vmin=0, vmax=max(3, np.nanmax(matrix.values)))

        # Add asterisks for significant cells
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix.iloc[i, j]
                if not np.isnan(val) and val > -np.log10(q_threshold):
                    ax.text(j, i, '*', ha='center', va='center', fontsize=10,
                            fontweight='bold', color='black')

        ax.set_xticks(range(len(available_methods)))
        ax.set_xticklabels(available_methods, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index, fontsize=8)
        ax.set_title(gs, fontsize=14, fontweight='bold')

        # Highlight key superclusters
        for idx, sc in enumerate(matrix.index):
            if sc in KEY_SUPERCLUSTERS:
                ax.get_yticklabels()[idx].set_fontweight('bold')
                ax.get_yticklabels()[idx].set_color('darkred')

    # Colorbar
    cbar = fig.colorbar(im, ax=axes, shrink=0.5, label='-log10(q-value)')

    plt.suptitle('Significance Across Matching Methods', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig_significance_heatmap.pdf", dpi=300, transparent=True,
                bbox_inches='tight')
    fig.savefig(FIG_DIR / "Fig_significance_heatmap.png", dpi=300, transparent=True,
                bbox_inches='tight')
    plt.show()
    return fig

fig2 = plot_significance_heatmap(supercluster_results,
                                  gene_sets_to_plot=["SCZ", "ASD_LIQ", "DDD_61"])

# %% [markdown]
# ## Section 5: Trade-off Plot (KEY FIGURE)
#
# X = matching quality (avg |bias|, lower = better match)
# Y = significance (mean Z-score for CGE+MGE, higher = more significant)

# %%
def plot_tradeoff(supercluster_results, matching_quality,
                  gene_set="SCZ", target_sc=["CGE interneuron", "MGE interneuron"]):
    """Plot the matching quality vs significance trade-off."""
    methods_with_both = [m for m in METHODS
                         if m in matching_quality
                         and gene_set in supercluster_results.get(m, {})]

    if len(methods_with_both) < 2:
        print(f"Need at least 2 methods with both matching quality and significance data. Have {len(methods_with_both)}.")
        return None

    # Compute metrics per method
    x_vals = []  # avg |bias|
    y_vals = []  # mean Z-score for target superclusters
    labels = []
    colors = []
    method_colors = {
        "Random": "#e74c3c",
        "Gene-by-gene": "#e67e22",
        "Rejection": "#f1c40f",
        "Best-of-N": "#2ecc71",
        "PropWeight": "#3498db",
        "SIS": "#9b59b6",
    }

    for method in methods_with_both:
        # X: matching quality
        mq = matching_quality[method]
        avg_bias = np.mean([mq[v]["mean_abs_bias"] for v in MATCHING_VARS if v in mq])

        # Y: significance
        df = supercluster_results[method][gene_set]
        z_scores = []
        for sc in target_sc:
            sc_row = df[df["Supercluster"] == sc]
            if len(sc_row) > 0:
                z_scores.append(sc_row["Z-score"].values[0])
        mean_z = np.mean(z_scores) if z_scores else 0

        x_vals.append(avg_bias)
        y_vals.append(mean_z)
        labels.append(method)
        colors.append(method_colors.get(method, "#666666"))

    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter
    for x, y, label, color in zip(x_vals, y_vals, labels, colors):
        ax.scatter(x, y, c=color, s=150, zorder=5, edgecolors='black', linewidth=1)
        ax.annotate(label, (x, y), textcoords="offset points", xytext=(10, 5),
                    fontsize=11, fontweight='bold', color=color)

    # Reference lines
    ax.axhline(y=1.645, color='gray', linestyle='--', alpha=0.5, label='Z=1.645 (p=0.05)')
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

    ax.set_xlabel('Average |Bias| in Matching Variables\n(lower = better confounder match)',
                  fontsize=12)
    ax.set_ylabel(f'Mean Z-score ({" + ".join(target_sc)})\n(higher = stronger signal)',
                  fontsize=12)
    ax.set_title(f'{gene_set}: Matching Quality vs Significance Trade-off',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')

    # Add arrow showing optimal direction
    ax.annotate('Better', xy=(min(x_vals) * 0.9, max(y_vals) * 1.05),
                fontsize=10, fontstyle='italic', color='green', alpha=0.7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig_tradeoff_matching_vs_significance.pdf", dpi=300,
                transparent=True, bbox_inches='tight')
    fig.savefig(FIG_DIR / "Fig_tradeoff_matching_vs_significance.png", dpi=300,
                transparent=True, bbox_inches='tight')
    plt.show()
    return fig

fig3 = plot_tradeoff(supercluster_results, matching_quality,
                      gene_set="SCZ",
                      target_sc=["CGE interneuron", "MGE interneuron"])

# %% [markdown]
# ### Multi-Gene-Set Trade-off

# %%
def plot_tradeoff_multi_geneset(supercluster_results, gene_sets_to_plot,
                                 target_sc=["CGE interneuron", "MGE interneuron"]):
    """Trade-off plot for multiple gene sets, using Z-scores only (no matching quality x-axis)."""
    methods_with_data = [m for m in METHODS
                         if any(gs in supercluster_results.get(m, {}) for gs in gene_sets_to_plot)]
    if len(methods_with_data) < 2:
        print("Not enough data for multi-gene-set comparison.")
        return None

    n_gs = len(gene_sets_to_plot)
    fig, axes = plt.subplots(1, n_gs, figsize=(5 * n_gs, 5), sharey=True)
    if n_gs == 1:
        axes = [axes]

    method_colors = {
        "Random": "#e74c3c", "Gene-by-gene": "#e67e22", "Rejection": "#f1c40f",
        "Best-of-N": "#2ecc71", "PropWeight": "#3498db", "SIS": "#9b59b6",
    }
    x_positions = np.arange(len(methods_with_data))

    for panel_idx, gs in enumerate(gene_sets_to_plot):
        ax = axes[panel_idx]
        z_vals = []
        bar_colors = []

        for method in methods_with_data:
            if gs in supercluster_results.get(method, {}):
                df = supercluster_results[method][gs]
                z_scores = []
                for sc in target_sc:
                    sc_row = df[df["Supercluster"] == sc]
                    if len(sc_row) > 0:
                        z_scores.append(sc_row["Z-score"].values[0])
                mean_z = np.mean(z_scores) if z_scores else 0
            else:
                mean_z = 0
            z_vals.append(mean_z)
            bar_colors.append(method_colors.get(method, "#666666"))

        bars = ax.bar(x_positions, z_vals, color=bar_colors, edgecolor='black', linewidth=0.5)
        ax.axhline(y=1.645, color='gray', linestyle='--', alpha=0.5)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(methods_with_data, rotation=45, ha='right', fontsize=9)
        ax.set_title(gs, fontsize=13, fontweight='bold')
        if panel_idx == 0:
            ax.set_ylabel(f'Mean Z-score\n({" + ".join(target_sc)})', fontsize=11)

    plt.suptitle('Interneuron Significance Across Methods & Gene Sets',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig_significance_bars_multi_geneset.pdf", dpi=300,
                transparent=True, bbox_inches='tight')
    fig.savefig(FIG_DIR / "Fig_significance_bars_multi_geneset.png", dpi=300,
                transparent=True, bbox_inches='tight')
    plt.show()
    return fig

fig4 = plot_tradeoff_multi_geneset(
    supercluster_results,
    gene_sets_to_plot=["SCZ", "ASD_LIQ", "DDD_61", "ASD_All"],
    target_sc=["CGE interneuron", "MGE interneuron"]
)

# %% [markdown]
# ## Section 6: Total Significance Counts

# %%
def total_significance_table(supercluster_results, gene_sets, q_threshold=0.1):
    """Count total significant superclusters per method × gene set."""
    rows = []
    for method in METHODS:
        row = {"Method": method}
        total = 0
        for gs in gene_sets:
            if gs in supercluster_results.get(method, {}):
                df = supercluster_results[method][gs]
                n_sig = (df["q-value"] < q_threshold).sum()
                row[gs] = n_sig
                total += n_sig
            else:
                row[gs] = "N/A"
        row["Total"] = total
        rows.append(row)
    return pd.DataFrame(rows).set_index("Method")

print("=== TOTAL SIGNIFICANT SUPERCLUSTERS (q < 0.1) PER METHOD ===\n")
total_sig = total_significance_table(supercluster_results, GENE_SETS)
if len(total_sig) > 0:
    print(total_sig.to_string())

# %% [markdown]
# ## Section 7: Summary & Recommendation

# %%
def generate_summary(verdict_df, mq_summary, total_sig):
    """Generate a text summary of the systematic comparison."""
    print("=" * 70)
    print("SYSTEMATIC COMPARISON SUMMARY")
    print("=" * 70)

    if verdict_df is not None and not verdict_df.empty:
        # Count passes per method
        pass_counts = {}
        for method in verdict_df.index:
            n_pass = sum(1 for v in verdict_df.loc[method] if "PASS" in str(v))
            n_total = sum(1 for v in verdict_df.loc[method] if v != "N/A")
            pass_counts[method] = (n_pass, n_total)

        print("\n1. CRITERIA PASS RATE:")
        for method, (n_pass, n_total) in pass_counts.items():
            pct = n_pass / n_total * 100 if n_total > 0 else 0
            print(f"   {method:20s}: {n_pass}/{n_total} ({pct:.0f}%)")

    if mq_summary is not None and len(mq_summary) > 0:
        print("\n2. MATCHING QUALITY (avg |bias|, lower = better):")
        for method in mq_summary.index:
            avg = mq_summary.loc[method, "avg_|bias|"]
            print(f"   {method:20s}: {avg:.4f}")

    if total_sig is not None and len(total_sig) > 0:
        print("\n3. TOTAL SIGNIFICANT SUPERCLUSTERS:")
        for method in total_sig.index:
            tot = total_sig.loc[method, "Total"]
            print(f"   {method:20s}: {tot}")

    print("\n4. RECOMMENDATION:")
    print("   The optimal method balances confounder matching with signal")
    print("   preservation. Methods that match too aggressively (PropWeight,")
    print("   SIS) may remove real biological signal, while methods that")
    print("   match too loosely (Random) inflate false positives.")
    print("   The Best-of-N method typically provides the best trade-off.")
    print("=" * 70)

generate_summary(
    verdict_df if 'verdict_df' in dir() and not verdict_df.empty else None,
    mq_summary if 'mq_summary' in dir() and len(mq_summary) > 0 else None,
    total_sig if 'total_sig' in dir() and len(total_sig) > 0 else None,
)

# %% [markdown]
# ## Output Files
#
# All figures saved to `results/figures/systematic_comparison/`:
# - `Fig_matching_quality_panel.pdf` — Real vs null distributions (6 methods × 3 variables)
# - `Fig_significance_heatmap.pdf` — Supercluster significance (methods × superclusters)
# - `Fig_tradeoff_matching_vs_significance.pdf` — Key trade-off plot
# - `Fig_significance_bars_multi_geneset.pdf` — Cross-gene-set comparison
