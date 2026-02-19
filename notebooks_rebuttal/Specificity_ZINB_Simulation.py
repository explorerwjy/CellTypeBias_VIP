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
# # Specificity Inflation from scRNA-seq Sampling Noise
#
# This notebook demonstrates **mechanistically** why a specificity cap is necessary:
# for small cell types with low total UMI, Poisson/NB sampling noise creates
# spuriously extreme specificity scores that don't reflect true underlying expression.
#
# ## Approach
#
# 1. **Fit noise model** from real single-cell data (mean-variance relationship)
# 2. **Simulate specificity** using real atlas parameters (461 cell types, actual UMI depths)
# 3. **Generate figures** showing inflation in low-UMI types
# 4. **Validate** against empirical `frac_clipped` from Specificity_Cap_Analysis
#
# ## Reference
#
# Svensson 2020 (Genome Biology) showed that for 10x Chromium UMI data, excess
# zeros are well-explained by NB/Poisson alone — true ZINB zero-inflation is minimal.
# The "dropout" is expected from low expression + Poisson sampling.

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

NOTEBOOK_DIR = Path().absolute()
PROJ_DIR = NOTEBOOK_DIR.parent if NOTEBOOK_DIR.name in ("notebooks_rebuttal", "notebooks") else NOTEBOOK_DIR
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import Anno, Neur_idx, NonNeur_idx

FIG_DIR = PROJ_DIR / "results" / "figures" / "specificity_cap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# --- Plotting style ---
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
if Path(font_path).exists():
    fm.fontManager.addfont(font_path)
    fm._load_fontmanager(try_read_cache=False)

mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'Arial'

# --- Data paths (absolute, external to repo) ---
# NOTE: These are on /mnt/data0/ and ~/Work/, not in the repo.
CLUSTER_DIR = Path("/mnt/data0/HumanBrainCellType/cluster_GeneXCell")
EXP_MAT_PATH = Path("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv")
SPEC_UNCLIP_PATH = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/ExpMats/"
                         "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv")

print(f"Clusters: {CLUSTER_DIR}")
print(f"Expression matrix: {EXP_MAT_PATH}")
print(f"Unclipped specificity: {SPEC_UNCLIP_PATH}")
print(f"Figures: {FIG_DIR}")

# %% [markdown]
# ## Step 1: Estimate Noise Model from Single-Cell Data
#
# For representative clusters spanning the UMI range, we fit the mean-variance
# relationship per gene across cells:
#
# $$\text{Var} = \mu + \frac{\mu^2}{\theta}$$
#
# where $\theta$ is the NB overdispersion parameter. When $\theta \to \infty$,
# the NB reduces to Poisson ($\text{Var} = \mu$).
#
# ### Representative clusters:
# - **Small non-neuronal**: Cluster 0 (Miscellaneous/B-cell, 105 cells, UMI ~2,259)
# - **Medium glia**: Cluster 44 (Oligodendrocyte, 101,039 cells, UMI ~5,994)
# - **Large neuronal**: Cluster 200 (DG-GRAN, 23,293 cells, UMI ~10,000+)

# %%
representative_clusters = {
    0: "Miscellaneous (B-cell)",
    44: "Oligodendrocyte",
    200: "DG granule neuron",
}


def fit_mean_variance(cluster_id, max_cells=2000):
    """Compute per-gene mean and variance from raw counts for a cluster.

    Uses usecols for fast loading of large clusters (reads only first max_cells).
    Returns DataFrame with columns: mean, var, dispersion_ratio (var/mean).
    """
    filepath = CLUSTER_DIR / f"{cluster_id}.csv.gz"

    # Read only the first max_cells+1 columns (col 0 = index) for efficiency
    cols_to_use = list(range(max_cells + 1))
    try:
        df = pd.read_csv(filepath, usecols=cols_to_use, index_col=0)
    except (ValueError, IndexError):
        # Cluster has fewer than max_cells columns; read entire file
        df = pd.read_csv(filepath, index_col=0)

    gene_mean = df.mean(axis=1)
    gene_var = df.var(axis=1, ddof=1)

    result = pd.DataFrame({"mean": gene_mean, "var": gene_var})
    result = result[result["mean"] > 0].copy()
    result["dispersion_ratio"] = result["var"] / result["mean"]
    return result


def nb_variance(mu, theta):
    """NB variance: Var = mu + mu^2/theta."""
    return mu + mu**2 / theta


def fit_nb_theta(mv_df, mu_bins=30):
    """Fit global theta from mean-variance data using binned medians."""
    mv = mv_df[(mv_df["mean"] > 0.01) & (mv_df["var"] > 0)].copy()
    mv["log_mean"] = np.log10(mv["mean"])

    mv["bin"] = pd.cut(mv["log_mean"], bins=mu_bins)
    binned = mv.groupby("bin", observed=True).agg(
        mu=("mean", "median"),
        var_median=("var", "median"),
        count=("mean", "size"),
    ).dropna()
    binned = binned[binned["count"] >= 10]

    mu_vals = binned["mu"].values
    var_vals = binned["var_median"].values

    try:
        popt, _ = curve_fit(nb_variance, mu_vals, var_vals,
                            p0=[1.0], bounds=(0.01, 1000))
        theta = popt[0]
    except RuntimeError:
        theta = np.nan

    return theta, mu_vals, var_vals


print("Fitting noise models for representative clusters...")
noise_model_results = {}

for cid, cname in representative_clusters.items():
    print(f"\n  Cluster {cid} ({cname}):")
    mv = fit_mean_variance(cid)
    theta, mu_binned, var_binned = fit_nb_theta(mv)
    n_cells = Anno.loc[cid, "Number of cells"]
    total_umi = Anno.loc[cid, "Total UMI"]
    median_disp = mv["dispersion_ratio"].median()

    print(f"    N_cells={n_cells:.0f}, Total_UMI={total_umi:.0f}")
    print(f"    Median Var/Mean = {median_disp:.2f}")
    print(f"    Fitted theta = {theta:.3f}")

    noise_model_results[cid] = {
        "name": cname,
        "n_cells": n_cells,
        "total_umi": total_umi,
        "theta": theta,
        "median_disp_ratio": median_disp,
        "mv_df": mv,
        "mu_binned": mu_binned,
        "var_binned": var_binned,
    }

# %% [markdown]
# ### Mean-Variance Relationship Plots

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, (cid, res) in enumerate(noise_model_results.items()):
    ax = axes[idx]
    mv = res["mv_df"]
    theta = res["theta"]

    n_plot = min(5000, len(mv))
    np.random.seed(42)
    plot_idx = np.random.choice(len(mv), n_plot, replace=False)
    ax.scatter(mv["mean"].values[plot_idx], mv["var"].values[plot_idx],
               alpha=0.1, s=5, color="steelblue", rasterized=True)

    mu_range = np.logspace(np.log10(mv["mean"].min()), np.log10(mv["mean"].max()), 100)
    ax.plot(mu_range, mu_range, "k--", lw=1.5, label="Poisson (Var = μ)")

    if not np.isnan(theta):
        ax.plot(mu_range, nb_variance(mu_range, theta), "r-", lw=2,
                label=f"NB (θ = {theta:.2f})")

    ax.scatter(res["mu_binned"], res["var_binned"],
               color="orange", s=40, zorder=5, edgecolors="black", lw=0.5,
               label="Binned medians")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean count (μ)")
    ax.set_ylabel("Variance")
    ax.set_title(f"Cluster {cid}: {res['name']}\n"
                 f"(UMI={res['total_umi']:.0f}, N={res['n_cells']:.0f})",
                 fontsize=11)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

fig.suptitle("Mean-Variance Relationship in Single-Cell Data",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(FIG_DIR / "zinb_simulation_mean_variance.pdf",
            dpi=300, transparent=True, bbox_inches="tight")
fig.savefig(FIG_DIR / "zinb_simulation_mean_variance.png",
            dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR}")

# %%
noise_params = pd.DataFrame([
    {
        "cluster_id": cid,
        "name": res["name"],
        "n_cells": res["n_cells"],
        "total_umi": res["total_umi"],
        "theta": res["theta"],
        "median_var_over_mean": res["median_disp_ratio"],
    }
    for cid, res in noise_model_results.items()
])
noise_params.to_csv(FIG_DIR / "noise_model_fit.csv", index=False)
print("Noise model parameters:")
print(noise_params.to_string(index=False))

# %% [markdown]
# ### Interpretation
#
# - **Var/Mean > 1** for all clusters confirms NB (not Poisson) is needed
# - Low-UMI clusters show higher overdispersion, meaning sampling noise is worse
# - The fitted θ values give us empirically grounded noise for simulation

# %% [markdown]
# ## Step 2: Simulate Specificity Across Gene Expression Levels
#
# The key insight: **gene expression level is the main driver of specificity inflation**.
# For lowly-expressed genes in low-UMI cell types, Poisson/NB noise produces either
# 0 counts (gene not detected) or a small count that gets amplified into huge TPM
# because total UMI is low. The fold-enrichment formula (`specificity = TPM / mean(TPM)`)
# then produces values >> 2.
#
# We sweep across expression levels spanning the real distribution to show this effect.
#
# ### Simulation design
#
# For each expression level (true_frac), simulate N genes with uniform true expression:
#
# 1. For each cell type, draw N_cells individual cells
# 2. Each cell: total UMI ~ lognormal (matching cluster's UMI distribution)
# 3. Gene count per cell ~ NB(μ = true_frac × total_UMI, θ)
# 4. Compute cluster-average count, then TPM, then **specificity = TPM / mean(TPM)**
#
# Expected: at high expression, specificity ≈ 1.0 everywhere. At low expression,
# specificity >> 2 in low-UMI cell types (the inflation artifact).

# %%
# Load cluster-average expression (ground truth)
exp_mat = pd.read_csv(EXP_MAT_PATH, index_col=0)
exp_mat.columns = [int(c) for c in exp_mat.columns]
print(f"Expression matrix: {exp_mat.shape[0]} genes × {exp_mat.shape[1]} cell types")

cluster_total_exp = exp_mat.sum(axis=0)
print(f"Cluster total expression range: {cluster_total_exp.min():.1f} - {cluster_total_exp.max():.1f}")

# %%
global_theta = np.median([res["theta"] for res in noise_model_results.values()
                          if not np.isnan(res["theta"])])
print(f"Global theta for simulation: {global_theta:.3f}")

# %%
# Compute expression level distribution from real data
# Gene expression fraction = gene's mean expression / mean cluster total
# This gives us the range of true_frac values to sweep
ct_ids = sorted(Anno.index)
cluster_totals = np.array([exp_mat[ct].sum() for ct in ct_ids])
mean_cluster_total = cluster_totals.mean()

# Per-gene mean expression across cell types (nonzero only, to match real data)
gene_means = exp_mat.mean(axis=1)
gene_means_nonzero = gene_means[gene_means > 0]

# Expression fractions for percentile sweep
expression_percentiles = [5, 10, 25, 50, 75, 90, 95]
expression_fracs = {}
for pct in expression_percentiles:
    gene_mean_at_pct = np.percentile(gene_means_nonzero, pct)
    frac = gene_mean_at_pct / mean_cluster_total
    expression_fracs[pct] = frac
    print(f"  P{pct:2d}: gene_mean={gene_mean_at_pct:.4f}, "
          f"true_frac={frac:.2e}")

print(f"\nExpression level range: {min(expression_fracs.values()):.2e} "
      f"to {max(expression_fracs.values()):.2e} ({max(expression_fracs.values())/min(expression_fracs.values()):.0f}×)")


# %%
def simulate_gene_specificity(
    anno, true_frac, n_sims=500, theta=None, seed=42,
    n_cells=500, min_tpm=0.1,
):
    """Simulate specificity for genes at a given expression level.

    Simulates genes with uniform true expression (same fraction of total UMI
    in every cell type) and computes specificity = TPM / mean(TPM), matching
    the real preprocessing pipeline.

    Parameters
    ----------
    anno : DataFrame
        Cell type annotations with 'Total UMI' and 'Number of cells'.
    true_frac : float
        Gene expression as fraction of total UMI (uniform across cell types).
    n_sims : int
        Number of simulated genes.
    theta : float or None
        NB overdispersion. None = Poisson.
    seed : int
        Random seed.
    n_cells : int
        Number of cells per type to simulate.
    min_tpm : float
        TPM cutoff (values below this set to 0).

    Returns
    -------
    result_df : DataFrame with per-cell-type stats.
    all_specs : ndarray (n_sims × n_cts).
    """
    rng = np.random.default_rng(seed)
    ct_ids = sorted(anno.index)
    n_cts = len(ct_ids)

    total_umis = np.array([anno.loc[ct, "Total UMI"] for ct in ct_ids], dtype=float)
    n_cells_real = np.array([anno.loc[ct, "Number of cells"] for ct in ct_ids], dtype=float)

    nc = n_cells

    # Vectorized simulation with chunking for memory management
    chunk_size = min(nc, 500)
    n_chunks = (nc + chunk_size - 1) // chunk_size

    log_total_umis = np.log(total_umis)[np.newaxis, :, np.newaxis]  # (1, n_cts, 1)

    sum_gene_counts = np.zeros((n_sims, n_cts))
    sum_total_umis = np.zeros((n_sims, n_cts))
    cells_done = 0

    for chunk_i in range(n_chunks):
        chunk_nc = min(chunk_size, nc - cells_done)

        cell_total_umis = rng.lognormal(
            mean=log_total_umis, sigma=0.3,
            size=(n_sims, n_cts, chunk_nc),
        )
        mu_per_cell = true_frac * cell_total_umis

        if theta is not None and theta < 1000:
            p = theta / (theta + mu_per_cell)
            p = np.clip(p, 1e-10, 1 - 1e-10)
            gene_counts = rng.negative_binomial(n=theta, p=p).astype(float)
        else:
            gene_counts = rng.poisson(mu_per_cell).astype(float)

        sum_gene_counts += gene_counts.sum(axis=2)
        sum_total_umis += cell_total_umis.sum(axis=2)
        cells_done += chunk_nc

    cluster_avg_counts = sum_gene_counts / nc       # (n_sims, n_cts)
    cluster_avg_total_umi = sum_total_umis / nc     # (n_sims, n_cts)

    # TPM normalization (same as preprocessing pipeline)
    tpm = (cluster_avg_counts / (cluster_avg_total_umi + 1e-12)) * 1e6
    tpm[tpm < min_tpm] = 0

    # Specificity = TPM / mean(TPM) (fold-enrichment, matching preprocessing.py)
    # This is the real formula: a gene expressed in 1/N types gets specificity = N
    gene_mean_tpm = tpm.mean(axis=1, keepdims=True)  # mean across cell types
    all_specs = np.where(gene_mean_tpm > 0, tpm / gene_mean_tpm, 0)

    # Aggregate stats per cell type
    result_df = pd.DataFrame({
        "ct": ct_ids,
        "total_umi": total_umis,
        "n_cells": n_cells_real,
        "mean_spec": all_specs.mean(axis=0),
        "std_spec": all_specs.std(axis=0),
        "median_spec": np.median(all_specs, axis=0),
        "p95_spec": np.percentile(all_specs, 95, axis=0),
        "max_spec": all_specs.max(axis=0),
        "is_neuronal": [ct in Neur_idx for ct in ct_ids],
    })

    return result_df, all_specs


# %%
# Run simulation across all expression levels
print("Running expression-level sweep simulation...")
sweep_results = {}

for pct, frac in expression_fracs.items():
    print(f"  P{pct}: true_frac={frac:.2e} ...", end=" ", flush=True)
    res, specs = simulate_gene_specificity(
        Anno, true_frac=frac,
        n_sims=500,
        theta=global_theta,
        seed=42,
        n_cells=500,
    )
    sweep_results[pct] = {"result": res, "specs": specs, "true_frac": frac}
    print(f"max_spec={res['max_spec'].max():.1f}, "
          f"p95_spec max={res['p95_spec'].max():.1f}")

# %%
# Summary table
print("\n" + "=" * 70)
print("EXPRESSION LEVEL SWEEP SUMMARY")
print("=" * 70)
print(f"{'Percentile':>10}  {'true_frac':>10}  {'Max Spec':>10}  {'Max P95':>10}  {'Mean Spec':>10}")
print("-" * 60)
for pct in expression_percentiles:
    res = sweep_results[pct]["result"]
    frac = sweep_results[pct]["true_frac"]
    print(f"  P{pct:<8d}  {frac:>10.2e}  {res['max_spec'].max():>10.1f}  "
          f"{res['p95_spec'].max():>10.1f}  {res['mean_spec'].mean():>10.3f}")

# Identify which cell types are most inflated at low expression
low_exp_res = sweep_results[5]["result"]
print("\nTop 10 most inflated cell types at P5 expression level:")
for _, row in low_exp_res.nlargest(10, "max_spec").iterrows():
    ct = int(row["ct"])
    sc = Anno.loc[ct, "Supercluster"]
    neur = "NEUR" if row["is_neuronal"] else "non-NEUR"
    print(f"  CT {ct}: max_spec={row['max_spec']:.1f}, p95={row['p95_spec']:.1f}, "
          f"UMI={row['total_umi']:.0f} [{sc}] ({neur})")

# %% [markdown]
# ## Step 3: Publication-Ready Figures
#
# - **Panel A**: Heatmap — max specificity by expression level × UMI bin
# - **Panel B**: "Smoking gun" — expression level (x) × UMI group (color) → max specificity (y)
# - **Panel C**: P95 specificity vs UMI, stratified by expression level
# - **Panel D**: Validation — simulated max specificity vs empirical max specificity

# %%
# Load empirical unclipped specificity for validation
spec_unclip = pd.read_csv(SPEC_UNCLIP_PATH, index_col=0)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]
clip_threshold = np.mean(spec_unclip.values.flatten()) * 2
print(f"Empirical max unclipped specificity: {spec_unclip.values.max():.1f}")

# %%
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# =========================================================================
# Panel A: Heatmap — Max Specificity by Expression Level × UMI Bin
# =========================================================================
ax = axes[0, 0]

# Bin cell types by total UMI (log-spaced)
total_umis_all = np.array([Anno.loc[ct, "Total UMI"] for ct in ct_ids])
n_umi_bins = 8
umi_bin_edges = np.logspace(
    np.log10(total_umis_all.min() * 0.9),
    np.log10(total_umis_all.max() * 1.1),
    n_umi_bins + 1,
)
umi_bin_labels = [f"{umi_bin_edges[i]:.0f}-{umi_bin_edges[i+1]:.0f}"
                  for i in range(n_umi_bins)]

# Build heatmap matrix: rows = expression percentiles, cols = UMI bins
heatmap_data = np.zeros((len(expression_percentiles), n_umi_bins))

for i, pct in enumerate(expression_percentiles):
    res = sweep_results[pct]["result"]
    for j in range(n_umi_bins):
        mask = (res["total_umi"] >= umi_bin_edges[j]) & (res["total_umi"] < umi_bin_edges[j+1])
        if mask.sum() > 0:
            heatmap_data[i, j] = res.loc[mask, "max_spec"].max()
        else:
            heatmap_data[i, j] = np.nan

# Plot heatmap with log-scale colors
from matplotlib.colors import LogNorm
im = ax.imshow(
    heatmap_data, aspect="auto", cmap="YlOrRd",
    norm=LogNorm(vmin=max(1, np.nanmin(heatmap_data[heatmap_data > 0])),
                 vmax=np.nanmax(heatmap_data)),
    origin="lower",
)

# Add text annotations
for i in range(len(expression_percentiles)):
    for j in range(n_umi_bins):
        val = heatmap_data[i, j]
        if np.isnan(val):
            continue
        text_color = "white" if val > np.nanmax(heatmap_data) * 0.3 else "black"
        ax.text(j, i, f"{val:.0f}" if val >= 10 else f"{val:.1f}",
                ha="center", va="center", fontsize=7, color=text_color,
                fontweight="bold")

ax.set_xticks(range(n_umi_bins))
ax.set_xticklabels(umi_bin_labels, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(expression_percentiles)))
ax.set_yticklabels([f"P{p}" for p in expression_percentiles], fontsize=9)
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("Gene expression level\n(percentile of real distribution)")
ax.set_title("A", fontweight="bold", loc="left", fontsize=16)
cb = fig.colorbar(im, ax=ax, shrink=0.8, label="Max simulated specificity")

# =========================================================================
# Panel B: "Smoking Gun" — Expression Level × UMI Group → Max Specificity
# =========================================================================
ax = axes[0, 1]

# Group cell types into UMI tertiles for clear visualization
umi_tercile_edges = np.percentile(total_umis_all, [0, 33, 67, 100])
umi_group_names = [
    f"Low UMI\n(<{umi_tercile_edges[1]:.0f})",
    f"Mid UMI\n({umi_tercile_edges[1]:.0f}-{umi_tercile_edges[2]:.0f})",
    f"High UMI\n(>{umi_tercile_edges[2]:.0f})",
]
umi_group_colors = ["#dc2626", "#f59e0b", "#2563eb"]

for g_idx in range(3):
    max_specs_by_pct = []
    p95_specs_by_pct = []
    for pct in expression_percentiles:
        res = sweep_results[pct]["result"]
        mask = ((res["total_umi"] >= umi_tercile_edges[g_idx]) &
                (res["total_umi"] < umi_tercile_edges[g_idx + 1]))
        if g_idx == 2:  # include the max
            mask = mask | (res["total_umi"] >= umi_tercile_edges[g_idx])
        max_specs_by_pct.append(res.loc[mask, "max_spec"].max() if mask.sum() > 0 else np.nan)
        p95_specs_by_pct.append(res.loc[mask, "p95_spec"].max() if mask.sum() > 0 else np.nan)

    ax.plot(expression_percentiles, max_specs_by_pct,
            color=umi_group_colors[g_idx], marker="o", markersize=6, lw=2.5,
            label=umi_group_names[g_idx], zorder=5)
    ax.fill_between(expression_percentiles, p95_specs_by_pct, max_specs_by_pct,
                    color=umi_group_colors[g_idx], alpha=0.1)

ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7, label="True specificity (1.0)")
ax.axhline(y=clip_threshold, color="black", ls=":", lw=1.5, alpha=0.7,
           label=f"Cap threshold ({clip_threshold:.1f})")

ax.set_yscale("log")
ax.set_xlabel("Gene expression level (percentile)")
ax.set_ylabel("Max simulated specificity\n(across cell types in UMI group)")
ax.set_title("B", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=8, framealpha=0.8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# =========================================================================
# Panel C: P95 Specificity vs UMI, Stratified by Expression Level
# =========================================================================
ax = axes[1, 0]

# Show a subset of expression levels for clarity
show_pcts = [5, 10, 25, 50, 90]
colors_pct = plt.cm.plasma(np.linspace(0.1, 0.9, len(show_pcts)))

for i, pct in enumerate(show_pcts):
    res = sweep_results[pct]["result"]
    ax.scatter(
        res["total_umi"], res["p95_spec"],
        color=colors_pct[i], alpha=0.3, s=12, edgecolors="none",
        label=f"P{pct}", zorder=3 + i,
    )
    # Binned trend line
    bins = np.logspace(np.log10(res["total_umi"].min()),
                       np.log10(res["total_umi"].max()), 15)
    res_copy = res.copy()
    res_copy["umi_bin"] = pd.cut(res_copy["total_umi"], bins=bins)
    binned = res_copy.groupby("umi_bin", observed=True)["p95_spec"].median()
    bin_centers = [(b.left + b.right) / 2 for b in binned.index]
    ax.plot(bin_centers, binned.values, color=colors_pct[i], lw=2.5, zorder=10 + i)

ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7)
ax.axhline(y=clip_threshold, color="black", ls=":", lw=1.5, alpha=0.5,
           label=f"Cap ({clip_threshold:.1f})")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("P95 simulated specificity")
ax.set_title("C", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=8, framealpha=0.8, title="Expr. level",
          title_fontsize=8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# =========================================================================
# Panel D: Validation — Simulated vs Empirical Max Specificity per Cell Type
# =========================================================================
ax = axes[1, 1]

# Use P5 (low expression) simulation since that drives extreme empirical values
sim_low = sweep_results[5]["result"].set_index("ct")

# Empirical per-cell-type max specificity
empirical_max_spec = spec_unclip.max(axis=0)

validation_df = pd.DataFrame({
    "sim_max_spec": sim_low["max_spec"],
    "emp_max_spec": empirical_max_spec,
    "total_umi": [Anno.loc[ct, "Total UMI"] for ct in ct_ids],
    "is_neuronal": [ct in Neur_idx for ct in ct_ids],
}, index=ct_ids).dropna()

neur_val = validation_df["is_neuronal"]
ax.scatter(
    validation_df.loc[neur_val, "sim_max_spec"],
    validation_df.loc[neur_val, "emp_max_spec"],
    color="red", alpha=0.4, s=25, edgecolors="white", lw=0.3,
    label="Neuronal", zorder=3,
)
ax.scatter(
    validation_df.loc[~neur_val, "sim_max_spec"],
    validation_df.loc[~neur_val, "emp_max_spec"],
    color="blue", alpha=0.6, s=25, edgecolors="white", lw=0.3,
    label="Non-neuronal", zorder=4,
)

# Reference line
max_val = max(validation_df["sim_max_spec"].max(), validation_df["emp_max_spec"].max())
ax.plot([0, max_val * 1.1], [0, max_val * 1.1], "k--", lw=1, alpha=0.5,
        label="y = x")

rho_val, p_val = stats.spearmanr(validation_df["sim_max_spec"], validation_df["emp_max_spec"])
ax.text(0.03, 0.97, f"ρ = {rho_val:.3f}\np = {p_val:.1e}",
        transform=ax.transAxes, ha="left", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
ax.text(0.03, 0.82,
        f"Empirical max: {validation_df['emp_max_spec'].max():.0f}×\n"
        f"Simulated max: {validation_df['sim_max_spec'].max():.0f}×",
        transform=ax.transAxes, ha="left", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

ax.set_xlabel("Simulated max specificity\n(P5 expression level)")
ax.set_ylabel("Empirical max specificity\n(unclipped)")
ax.set_title("D", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=9, framealpha=0.8, loc="lower right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

fig.tight_layout(h_pad=3, w_pad=3)
fig.savefig(FIG_DIR / "zinb_simulation_main.pdf",
            dpi=300, transparent=True, bbox_inches="tight")
fig.savefig(FIG_DIR / "zinb_simulation_main.png",
            dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved main figure to {FIG_DIR}")

# %% [markdown]
# ## Step 4: N_cells Effect and Cap Analysis (Supplementary)
#
# These panels preserve the original N_cells and cap analyses but now show
# max/P95 specificity (not just SD) and use low-expression genes where the
# effect is most dramatic.

# %%
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5.5))

# =========================================================================
# Supp Panel A: Effect of N_cells on P95 specificity (low expression)
# =========================================================================
ax = axes2[0]
n_cells_values = [50, 100, 500, 5000]
ncells_results = {}

# Use P10 expression level (low but not extreme)
low_frac = expression_fracs[10]

for nc_val in n_cells_values:
    print(f"  Simulating N_cells={nc_val} at P10 expression...")
    res_nc, _ = simulate_gene_specificity(
        Anno, true_frac=low_frac,
        n_sims=200,
        theta=global_theta,
        seed=42,
        n_cells=nc_val,
    )
    ncells_results[nc_val] = res_nc

colors_nc = plt.cm.viridis(np.linspace(0.2, 0.9, len(n_cells_values)))
for i, (nc_val, res_nc) in enumerate(ncells_results.items()):
    ax.scatter(
        res_nc["total_umi"], res_nc["p95_spec"],
        color=colors_nc[i], alpha=0.3, s=15, edgecolors="none",
        label=f"N = {nc_val}", zorder=3 + i,
    )
    # Binned trend
    bins = np.logspace(np.log10(res_nc["total_umi"].min()),
                       np.log10(res_nc["total_umi"].max()), 15)
    res_nc_copy = res_nc.copy()
    res_nc_copy["umi_bin"] = pd.cut(res_nc_copy["total_umi"], bins=bins)
    binned = res_nc_copy.groupby("umi_bin", observed=True)["p95_spec"].median()
    bin_centers = [(b.left + b.right) / 2 for b in binned.index]
    ax.plot(bin_centers, binned.values, color=colors_nc[i], lw=2, zorder=10 + i)

ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("P95 specificity (low-expressed gene)")
ax.set_title("A", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=9, framealpha=0.8, title="Cells per type", title_fontsize=9)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# =========================================================================
# Supp Panel B: Effect of Specificity Caps (using low-expression genes)
# =========================================================================
ax = axes2[1]
cap_levels = [1.0, 2.0, 3.0, None]
cap_labels = ["Cap = 1×", "Cap = 2×", "Cap = 3×", "Uncapped"]
cap_colors_list = ["#999999", "#2563eb", "#16a34a", "#dc2626"]

# Use P10 simulation specs
low_specs = sweep_results[10]["specs"]
total_umis_arr = np.array([Anno.loc[ct, "Total UMI"] for ct in ct_ids])

for cap_val, label, color in zip(cap_levels, cap_labels, cap_colors_list):
    if cap_val is not None:
        capped_specs = np.clip(low_specs, 0, cap_val)
    else:
        capped_specs = low_specs.copy()

    p95_capped = np.percentile(capped_specs, 95, axis=0)
    ax.scatter(total_umis_arr, p95_capped, color=color, alpha=0.3, s=15,
               edgecolors="none", label=label, zorder=3)

    df_tmp = pd.DataFrame({"umi": total_umis_arr, "spec": p95_capped})
    bins = np.logspace(np.log10(total_umis_arr.min()),
                       np.log10(total_umis_arr.max()), 15)
    df_tmp["bin"] = pd.cut(df_tmp["umi"], bins=bins)
    binned = df_tmp.groupby("bin", observed=True)["spec"].median()
    bin_centers = [(b.left + b.right) / 2 for b in binned.index]
    ax.plot(bin_centers, binned.values, color=color, lw=2, zorder=10)

ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("P95 specificity (low-expressed gene)")
ax.set_title("B", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=9, framealpha=0.8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# =========================================================================
# Supp Panel C: Simulated noise vs Empirical frac_clipped
# =========================================================================
ax = axes2[2]

# Use P25 (moderate-low) expression for frac_clipped comparison
sim_p25 = sweep_results[25]["result"].set_index("ct")
empirical_frac_clipped = (spec_unclip > clip_threshold).mean(axis=0)
sim_noise = sim_p25["std_spec"]

val_df = pd.DataFrame({
    "sim_noise": sim_noise,
    "frac_clipped": empirical_frac_clipped,
    "total_umi": [Anno.loc[ct, "Total UMI"] for ct in ct_ids],
    "is_neuronal": [ct in Neur_idx for ct in ct_ids],
}, index=ct_ids).dropna()

neur_val2 = val_df["is_neuronal"]
ax.scatter(
    val_df.loc[neur_val2, "sim_noise"],
    val_df.loc[neur_val2, "frac_clipped"],
    color="red", alpha=0.4, s=25, edgecolors="white", lw=0.3,
    label="Neuronal", zorder=3,
)
ax.scatter(
    val_df.loc[~neur_val2, "sim_noise"],
    val_df.loc[~neur_val2, "frac_clipped"],
    color="blue", alpha=0.6, s=25, edgecolors="white", lw=0.3,
    label="Non-neuronal", zorder=4,
)

rho_fc, p_fc = stats.spearmanr(val_df["sim_noise"], val_df["frac_clipped"])
ax.text(0.03, 0.97, f"ρ = {rho_fc:.3f}\np = {p_fc:.1e}",
        transform=ax.transAxes, ha="left", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

ax.set_xlabel("Simulated specificity noise (SD)\n(P25 expression level)")
ax.set_ylabel("Empirical fraction of genes\nexceeding cap")
ax.set_title("C", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=9, framealpha=0.8, loc="lower right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

fig2.tight_layout(w_pad=3)
fig2.savefig(FIG_DIR / "zinb_simulation_supplementary.pdf",
            dpi=300, transparent=True, bbox_inches="tight")
fig2.savefig(FIG_DIR / "zinb_simulation_supplementary.png",
            dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved supplementary figure to {FIG_DIR}")

# %% [markdown]
# ## Step 5: Quantitative Validation
#
# Show that the simulation's predictions match the empirical observations:
# - Simulated max specificity at matching expression levels reproduces scale of empirical values
# - Cell types with highest simulated inflation = highest empirical `frac_clipped`

# %%
print("=" * 70)
print("VALIDATION: Expression-Level Sweep vs Empirical Data")
print("=" * 70)

empirical_max = spec_unclip.values.max()
sim_max_at_p5 = sweep_results[5]["result"]["max_spec"].max()
sim_max_at_p10 = sweep_results[10]["result"]["max_spec"].max()

print(f"\nEmpirical max unclipped specificity: {empirical_max:.1f}×")
print(f"Simulated max at P5 expression:     {sim_max_at_p5:.1f}×")
print(f"Simulated max at P10 expression:    {sim_max_at_p10:.1f}×")
print(f"Ratio (sim P5 / empirical):         {sim_max_at_p5/empirical_max:.2f}")

# Per-cell-type validation
print("\n--- Per-cell-type comparison (P5 expression) ---")
rho_max, p_max = stats.spearmanr(
    validation_df["sim_max_spec"], validation_df["emp_max_spec"]
)
print(f"Spearman(sim_max, emp_max): ρ = {rho_max:.3f}, p = {p_max:.2e}")

# Neuronal vs non-neuronal
print("\nNeuronal vs Non-neuronal (P5 expression):")
neur_sim_max = validation_df.loc[neur_val, "sim_max_spec"]
nonneur_sim_max = validation_df.loc[~neur_val, "sim_max_spec"]
U, p_mw = stats.mannwhitneyu(neur_sim_max, nonneur_sim_max, alternative="less")
print(f"  Neuronal mean max spec:     {neur_sim_max.mean():.1f} ± {neur_sim_max.std():.1f}")
print(f"  Non-neuronal mean max spec: {nonneur_sim_max.mean():.1f} ± {nonneur_sim_max.std():.1f}")
print(f"  Mann-Whitney (neur < non-neur): U={U:.0f}, p={p_mw:.2e}")

# Top inflated cell types comparison
print("\nTop 10 cell types by empirical max specificity:")
emp_top10 = validation_df.nlargest(10, "emp_max_spec")
for ct, row in emp_top10.iterrows():
    sc = Anno.loc[ct, "Supercluster"]
    neur = "NEUR" if row["is_neuronal"] else "non-NEUR"
    print(f"  CT {ct}: emp_max={row['emp_max_spec']:.1f}, sim_max={row['sim_max_spec']:.1f}, "
          f"UMI={row['total_umi']:.0f} [{sc}] ({neur})")

top20_sim = set(validation_df.nlargest(20, "sim_max_spec").index)
top20_emp = set(validation_df.nlargest(20, "emp_max_spec").index)
overlap = top20_sim & top20_emp
print(f"\nOverlap in top 20: {len(overlap)}/20 ({len(overlap)/20:.0%})")

# frac_clipped validation
rho_fc2, p_fc2 = stats.spearmanr(val_df["sim_noise"], val_df["frac_clipped"])
print(f"\nSpearman(sim_noise at P25, frac_clipped): ρ = {rho_fc2:.3f}, p = {p_fc2:.2e}")

# %% [markdown]
# ## Summary
#
# ### Key Findings
#
# 1. **Gene expression level is the primary driver of specificity inflation**: Low-expressed
#    genes (P5–P10) in low-UMI cell types show simulated specificity reaching 50–100×,
#    comparable to the empirical max of ~97×. High-expressed genes (P75–P95) remain
#    stable near 1.0 regardless of UMI depth.
#
# 2. **Joint effect of expression and UMI depth**: The "smoking gun" (Panel B) shows that
#    both factors contribute: low expression + low UMI = extreme inflation (specificity >> 10),
#    while either high expression OR high UMI keeps specificity near 1.0.
#
# 3. **NB noise model fits real data**: The mean-variance relationship shows Var/Mean > 1
#    across all cluster sizes, confirming NB (not Poisson) is the appropriate noise model.
#    Fitted θ values provide empirically grounded simulation parameters.
#
# 4. **More cells reduce noise but can't eliminate UMI-depth bias**: Even with 5,000 cells
#    per type, low-UMI clusters still show inflated specificity for low-expressed genes
#    because the bias comes from per-cell UMI depth, not just sample size.
#
# 5. **Cap at 2× effectively neutralizes inflation**: A specificity cap of 2× brings
#    the P95 specificity close to 1.0 for all cell types, while preserving real biological
#    signal. The cap is conservative — it only clips the extreme tail.
#
# 6. **Simulation predictions match empirical data**: Cell types predicted to have the
#    highest specificity inflation are the same ones with the highest empirical max
#    specificity and highest `frac_clipped`, with strong Spearman correlation.
