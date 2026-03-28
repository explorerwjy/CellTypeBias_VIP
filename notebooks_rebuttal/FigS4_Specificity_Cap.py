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
# # Figure S4: Specificity Cap Validation
#
# This notebook produces a consolidated 6-panel supplementary figure demonstrating
# why the specificity cap at 2× mean is necessary, robust, and biologically appropriate.
#
# **Panel layout (3 rows × 2 columns):**
# - **A**: Empirical specificity inflation vs UMI depth (the problem)
# - **B**: NB simulation — max specificity by expression level × UMI group (the mechanism)
# - **C**: Spearman ρ across cap levels for 3 disorders (robustness)
# - **D**: Supercluster rank across cap levels (neuronal vs non-neuronal divergence)
# - **E**: Neuronal fraction in top-N: capped vs TDEP uncapped (external validation)
# - **F**: Leave-one-out stability: capped vs uncapped (single-gene dominance)
#
# **Source notebooks** (each section derived from):
# - Panel A: `Specificity_Cap_Analysis.py`
# - Panel B: `Specificity_ZINB_Simulation.py`
# - Panels C-D: `Cap_Sensitivity_Figure.py`
# - Panels E-F: `TDEP_Specificity_Comparison.py`

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import (
    HumanCT_AvgZ_Weighted, AnnotateCTDat, Fil2Dict,
    Anno, Neur_idx, NonNeur_idx, LoadGeneINFO,
)

FIG_DIR = PROJ_DIR / "results" / "figures" / "FigS4"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Font setup
font_path = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
if Path(font_path).exists():
    fm.fontManager.addfont(font_path)
    fm._load_fontmanager(try_read_cache=False)

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "font.size": 12,
    "font.family": "Arial",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

# %% [markdown]
# ---
# ## Section 1: Load Shared Data

# %%
# Unclipped specificity matrix (clip100 is effectively unclipped; max observed ~97)
spec_unclip = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv",
    index_col=0,
)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]

mean_spec = np.mean(spec_unclip.values.flatten())
clip_threshold = mean_spec * 2
n_genes, n_cts = spec_unclip.shape
print(f"Matrix: {n_genes} genes × {n_cts} cell types")
print(f"Global mean specificity: {mean_spec:.4f}")
print(f"Clip threshold (2× mean): {clip_threshold:.4f}")
print(f"Max unclipped specificity: {spec_unclip.values.max():.1f}×")

# %%
# Gene weights for three disorders
gw_files = {
    "SCZ": PROJ_DIR / "dat" / "GeneWeights" / "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "ASD (HIQ)": PROJ_DIR / "dat" / "GeneWeights" / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "DDD": PROJ_DIR / "dat" / "GeneWeights" / "DDD.top61.gw.bgmr.csv",
}
gw_dicts = {}
for label, fpath in gw_files.items():
    gw_dicts[label] = Fil2Dict(str(fpath))
    print(f"{label}: {len(gw_dicts[label])} genes")

# %%
# Annotation
anno = pd.read_excel(PROJ_DIR / "dat" / "annotation.xlsx", index_col=0)
anno = anno[anno.index.notna()]
anno.index = anno.index.astype(int)

# Cell type IDs
ct_ids = sorted(Anno.index)

# %% [markdown]
# ---
# ## Section 2: Panel A — Empirical Specificity Inflation
#
# Shows that low-UMI cell types have a higher fraction of genes exceeding
# the specificity cap, confirming that extreme specificity is a technical artifact.

# %%
ct_stats = pd.DataFrame(index=Anno.index)
for ct in ct_stats.index:
    vals = spec_unclip[ct].values
    ct_stats.loc[ct, "frac_clipped"] = np.mean(vals > clip_threshold)
    ct_stats.loc[ct, "max_spec"] = np.max(vals)

ct_stats["Total_UMI"] = Anno["Total UMI"]
ct_stats["N_cells"] = Anno["Number of cells"]
ct_stats["Supercluster"] = Anno["Supercluster"]
ct_stats["is_neuronal"] = ct_stats.index.isin(Neur_idx)

for col in ["frac_clipped", "max_spec", "Total_UMI", "N_cells"]:
    ct_stats[col] = pd.to_numeric(ct_stats[col])

neur_mask = ct_stats["is_neuronal"]
rho_umi, p_umi = stats.spearmanr(ct_stats["Total_UMI"], ct_stats["frac_clipped"])
print(f"Spearman(Total UMI, frac_clipped): ρ = {rho_umi:.3f}, p = {p_umi:.2e}")
print(f"Neuronal: {neur_mask.sum()}, Non-neuronal: {(~neur_mask).sum()}")

# %% [markdown]
# ---
# ## Section 3: Panel B — NB Simulation
#
# Demonstrates mechanistically that sampling noise in low-UMI cell types
# inflates specificity for lowly expressed genes.

# %%
# Try to load pre-computed ZINB data; if not available, recompute
PLOT_DATA_DIR = PROJ_DIR / "results" / "figures" / "plot_data"
zinb_pkl = PLOT_DATA_DIR / "zinb_data.pkl"

if zinb_pkl.exists():
    import pickle
    with open(zinb_pkl, "rb") as f:
        zinb_data = pickle.load(f)
    sweep_results_df = zinb_data["sweep_results"]
    expression_percentiles = zinb_data["expression_percentiles"]
    expression_fracs = zinb_data["expression_fracs"]
    zinb_clip_threshold = zinb_data["clip_threshold"]
    global_theta = zinb_data["global_theta"]
    print(f"Loaded pre-computed ZINB data from {zinb_pkl}")
    print(f"  Expression percentiles: {expression_percentiles}")
    print(f"  Global theta: {global_theta:.3f}")
    print(f"  Clip threshold: {zinb_clip_threshold:.4f}")
    USE_PRECOMPUTED = True
else:
    print("Pre-computed ZINB data not found. Running simulation...")
    USE_PRECOMPUTED = False

# %%
if not USE_PRECOMPUTED:
    from scipy.optimize import curve_fit

    # External data paths
    CLUSTER_DIR = Path("/mnt/data0/HumanBrainCellType/cluster_GeneXCell")
    EXP_MAT_PATH = Path("/home/jw3514/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv")

    # Fit noise model from representative clusters
    representative_clusters = {0: "Miscellaneous (B-cell)", 44: "Oligodendrocyte", 200: "DG granule neuron"}

    def fit_mean_variance(cluster_id, max_cells=2000):
        filepath = CLUSTER_DIR / f"{cluster_id}.csv.gz"
        cols_to_use = list(range(max_cells + 1))
        try:
            df = pd.read_csv(filepath, usecols=cols_to_use, index_col=0)
        except (ValueError, IndexError):
            df = pd.read_csv(filepath, index_col=0)
        gene_mean = df.mean(axis=1)
        gene_var = df.var(axis=1, ddof=1)
        result = pd.DataFrame({"mean": gene_mean, "var": gene_var})
        result = result[result["mean"] > 0].copy()
        result["dispersion_ratio"] = result["var"] / result["mean"]
        return result

    def nb_variance(mu, theta):
        return mu + mu**2 / theta

    def fit_nb_theta(mv_df, mu_bins=30):
        mv = mv_df[(mv_df["mean"] > 0.01) & (mv_df["var"] > 0)].copy()
        mv["log_mean"] = np.log10(mv["mean"])
        mv["bin"] = pd.cut(mv["log_mean"], bins=mu_bins)
        binned = mv.groupby("bin", observed=True).agg(
            mu=("mean", "median"), var_median=("var", "median"), count=("mean", "size"),
        ).dropna()
        binned = binned[binned["count"] >= 10]
        try:
            popt, _ = curve_fit(nb_variance, binned["mu"].values, binned["var_median"].values,
                                p0=[1.0], bounds=(0.01, 1000))
            return popt[0]
        except RuntimeError:
            return np.nan

    thetas = []
    for cid in representative_clusters:
        mv = fit_mean_variance(cid)
        theta = fit_nb_theta(mv)
        thetas.append(theta)
        print(f"  Cluster {cid}: theta = {theta:.3f}")

    global_theta = np.median([t for t in thetas if not np.isnan(t)])
    print(f"Global theta: {global_theta:.3f}")

    # Expression level sweep
    exp_mat = pd.read_csv(EXP_MAT_PATH, index_col=0)
    exp_mat.columns = [int(c) for c in exp_mat.columns]
    cluster_totals = np.array([exp_mat[ct].sum() for ct in ct_ids])
    mean_cluster_total = cluster_totals.mean()
    gene_means = exp_mat.mean(axis=1)
    gene_means_nonzero = gene_means[gene_means > 0]

    expression_percentiles = [10, 25, 50, 75, 90, 95]
    expression_fracs = {}
    for pct in expression_percentiles:
        expression_fracs[pct] = np.percentile(gene_means_nonzero, pct) / mean_cluster_total

    def simulate_gene_specificity(anno_df, true_frac, n_sims=500, theta=None, seed=42, n_cells=500, min_tpm=0.1):
        rng = np.random.default_rng(seed)
        ct_ids_sim = sorted(anno_df.index)
        n_cts_sim = len(ct_ids_sim)
        total_umis = np.array([anno_df.loc[ct, "Total UMI"] for ct in ct_ids_sim], dtype=float)
        nc = n_cells
        chunk_size = min(nc, 500)
        n_chunks = (nc + chunk_size - 1) // chunk_size
        log_total_umis = np.log(total_umis)[np.newaxis, :, np.newaxis]
        sum_gene_counts = np.zeros((n_sims, n_cts_sim))
        sum_total_umis_arr = np.zeros((n_sims, n_cts_sim))
        cells_done = 0
        for _ in range(n_chunks):
            chunk_nc = min(chunk_size, nc - cells_done)
            cell_total_umis = rng.lognormal(mean=log_total_umis, sigma=0.3, size=(n_sims, n_cts_sim, chunk_nc))
            mu_per_cell = true_frac * cell_total_umis
            if theta is not None and theta < 1000:
                p = theta / (theta + mu_per_cell)
                p = np.clip(p, 1e-10, 1 - 1e-10)
                gene_counts = rng.negative_binomial(n=theta, p=p).astype(float)
            else:
                gene_counts = rng.poisson(mu_per_cell).astype(float)
            sum_gene_counts += gene_counts.sum(axis=2)
            sum_total_umis_arr += cell_total_umis.sum(axis=2)
            cells_done += chunk_nc
        cluster_avg_counts = sum_gene_counts / nc
        cluster_avg_total_umi = sum_total_umis_arr / nc
        tpm = (cluster_avg_counts / (cluster_avg_total_umi + 1e-12)) * 1e6
        tpm[tpm < min_tpm] = 0
        gene_mean_tpm = tpm.mean(axis=1, keepdims=True)
        all_specs = np.where(gene_mean_tpm > 0, tpm / gene_mean_tpm, 0)
        n_cells_real = np.array([anno_df.loc[ct, "Number of cells"] for ct in ct_ids_sim], dtype=float)
        result_df = pd.DataFrame({
            "ct": ct_ids_sim, "total_umi": total_umis, "n_cells": n_cells_real,
            "mean_spec": all_specs.mean(axis=0), "max_spec": all_specs.max(axis=0),
            "p95_spec": np.percentile(all_specs, 95, axis=0),
            "is_neuronal": [ct in Neur_idx for ct in ct_ids_sim],
        })
        return result_df, all_specs

    sweep_results_df = {}
    for pct, frac in expression_fracs.items():
        print(f"  P{pct}: true_frac={frac:.2e} ...", end=" ", flush=True)
        res, _ = simulate_gene_specificity(Anno, true_frac=frac, n_sims=500, theta=global_theta, seed=42, n_cells=500)
        sweep_results_df[pct] = res
        print(f"max_spec={res['max_spec'].max():.1f}")

    zinb_clip_threshold = clip_threshold

# %%
# Prepare Panel B data: max specificity by UMI group × expression level
# Use round UMI thresholds that align with Panel A's log-scale x-axis
total_umis_all = np.array([Anno.loc[ct, "Total UMI"] for ct in ct_ids])
umi_group_edges = [total_umis_all.min(), 10_000, 30_000, total_umis_all.max() + 1]
umi_group_names = [
    "Low UMI (<10k)",
    "Mid UMI (10k–30k)",
    "High UMI (>30k)",
]

n_per_group = []
for g_idx in range(3):
    mask = (total_umis_all >= umi_group_edges[g_idx]) & (total_umis_all < umi_group_edges[g_idx + 1])
    n_per_group.append(mask.sum())
print(f"UMI groups: {list(zip(umi_group_names, n_per_group))}")

panelB_data = {}  # {group_idx: [max_spec_at_each_pct]}
for g_idx in range(3):
    max_specs = []
    for pct in expression_percentiles:
        res = sweep_results_df[pct]
        mask = ((res["total_umi"] >= umi_group_edges[g_idx]) &
                (res["total_umi"] < umi_group_edges[g_idx + 1]))
        max_specs.append(res.loc[mask, "max_spec"].max() if mask.sum() > 0 else np.nan)
    panelB_data[g_idx] = max_specs

print("Panel B data prepared (max specificity by UMI group × expression level)")

# %% [markdown]
# ---
# ## Section 4: Panels C & D — Cap Sensitivity Sweep
#
# Recomputes mutation bias at 7 cap levels for 3 disorders.

# %%
CAP_LEVELS = [1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 100.0]
ref_cap = 2.0

bias_all = {}
for disorder, gw_dict in gw_dicts.items():
    bias_all[disorder] = {}
    for cap in CAP_LEVELS:
        threshold = mean_spec * cap
        spec_clipped = spec_unclip.clip(lower=0, upper=threshold)
        # No mean-centering: use raw clipped specificity so cap sensitivity
        # reflects the full effect of extreme values without partial correction
        bias_df = HumanCT_AvgZ_Weighted(spec_clipped, gw_dict)
        bias_df = AnnotateCTDat(bias_df, Anno)
        bias_all[disorder][cap] = bias_df
    print(f"{disorder}: computed bias at {len(CAP_LEVELS)} cap levels")

# %%
# Panel C: Spearman correlation vs cap=2x
# NOTE: AnnotateCTDat sorts by EFFECT descending, so indices differ across cap levels.
# Must align by cell type ID before computing Spearman correlation.
corr_results = {}
for disorder in gw_dicts:
    ref_df = bias_all[disorder][ref_cap].set_index(bias_all[disorder][ref_cap].index)
    corr_results[disorder] = {}
    for cap in CAP_LEVELS:
        other_df = bias_all[disorder][cap].set_index(bias_all[disorder][cap].index)
        common = ref_df.index.intersection(other_df.index)
        rho, _ = stats.spearmanr(ref_df.loc[common, "EFFECT"], other_df.loc[common, "EFFECT"])
        corr_results[disorder][cap] = rho

corr_df = pd.DataFrame(corr_results)
print("Spearman ρ vs cap=2×:")
print(corr_df.round(4))

# %%
# Panel D: Supercluster rank across cap levels (SCZ)
# Select key neuronal + non-neuronal superclusters
panel_scs = ["CGE interneuron", "MGE interneuron", "Medium spiny neuron",
             "Vascular", "Microglia", "Astrocyte", "Ependymal"]

show_disorder = "SCZ"
sc_rank_data = {}
for cap in CAP_LEVELS:
    df = bias_all[show_disorder][cap]
    sc_means = df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
    for rank_i, (sc_name, _) in enumerate(sc_means.items(), 1):
        if sc_name not in sc_rank_data:
            sc_rank_data[sc_name] = {}
        sc_rank_data[sc_name][cap] = rank_i

print(f"\nSupercluster ranks at cap=2× vs 10× ({show_disorder}):")
for sc in panel_scs:
    r2 = sc_rank_data[sc][2.0]
    r10 = sc_rank_data[sc][10.0]
    arrow = "↑" if r10 < r2 else ("↓" if r10 > r2 else "=")
    print(f"  {sc:30s}: rank {r2:2d} (2×) → {r10:2d} (10×)  {arrow}")

# %% [markdown]
# ---
# ## Section 5: Panels E & F — TDEP Comparison and LOO Stability
#
# Panel E: Compares neuronal fraction in top-N between capped and TDEP uncapped specificity.
# Panel F: Leave-one-out stability (capped vs uncapped).

# %%
# Load TDEP specificity
tdep_spec = pd.read_csv(
    "/home/jw3514/Work/NeuralP/TDEP-sLDSC/data/cluster.specificity_matrix_entrez.csv",
    index_col=0,
)
tdep_spec.columns = [int(c) - 1 for c in tdep_spec.columns]  # 1-based → 0-based
print(f"TDEP specificity: {tdep_spec.shape[0]} genes × {tdep_spec.shape[1]} cell types")

# Scale TDEP proportions to fold-enrichment and mean-center
n_cts_tdep = tdep_spec.shape[1]
tdep_fold = tdep_spec * n_cts_tdep
tdep_mc = tdep_fold.sub(tdep_fold.mean(axis=1), axis=0)
print(f"TDEP max fold-enrichment: {tdep_fold.values.max():.1f}×")

# %%
# Our capped, mean-centered specificity
our_spec = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    index_col=0,
)
our_spec.columns = [int(c) for c in our_spec.columns]

# Gene weights for TDEP comparison
# Use full ASD set (not HIQ subset) — HIQ shows low neuronal enrichment even with capping
gw_asd_full = Fil2Dict(str(PROJ_DIR / "dat" / "GeneWeights" / "Spark_Meta_EWS.GeneWeight.bgmr.csv"))
print(f"ASD (full): {len(gw_asd_full)} genes")

gw_tdep = {
    "ASD": gw_asd_full,
    "SCZ": gw_dicts["SCZ"],
    "DDD": gw_dicts["DDD"],
}

# Compute bias with both metrics
results_our = {}
results_tdep = {}
for label, gw in gw_tdep.items():
    bias_our = HumanCT_AvgZ_Weighted(our_spec, gw)
    bias_our = AnnotateCTDat(bias_our, anno)
    bias_our["is_neuronal"] = bias_our.index.isin(Neur_idx)
    results_our[label] = bias_our

    bias_tdep = HumanCT_AvgZ_Weighted(tdep_mc, gw)
    bias_tdep = AnnotateCTDat(bias_tdep, anno)
    bias_tdep["is_neuronal"] = bias_tdep.index.isin(Neur_idx)
    results_tdep[label] = bias_tdep

    top20_our = bias_our.head(20)["is_neuronal"].sum()
    top20_tdep = bias_tdep.head(20)["is_neuronal"].sum()
    print(f"{label}: Top-20 neuronal — Ours: {top20_our}/20, TDEP: {top20_tdep}/20")

# %%
# Panel E: Neuronal fraction vs top-N cutoff
cutoffs = [5, 10, 15, 20, 30, 50]
panelE_data = {}
for disorder in gw_tdep:
    panelE_data[disorder] = {"capped": [], "tdep": []}
    for k in cutoffs:
        n_o = results_our[disorder].head(k)["is_neuronal"].sum()
        n_t = results_tdep[disorder].head(k)["is_neuronal"].sum()
        panelE_data[disorder]["capped"].append(n_o / k)
        panelE_data[disorder]["tdep"].append(n_t / k)

# %%
# Panel F: LOO stability — capped vs uncapped for SCZ and ASD
# Uncapped, mean-centered specificity
spec_uncap = spec_unclip.subtract(spec_unclip.mean(axis=1), axis=0)

loo_gw = {"SCZ": gw_tdep["SCZ"], "ASD": gw_tdep["ASD"]}


def run_loo_analysis(spec_mat, gw_dict, label=""):
    """Run leave-one-out analysis: remove each gene, recompute bias, measure rank stability.

    Returns dict with 'rhos', 'genes', 'full_bias', and 'loo_biases' (per-gene LOO bias).
    """
    full_bias = HumanCT_AvgZ_Weighted(spec_mat, gw_dict)
    genes = list(gw_dict.keys())
    rhos = []
    loo_biases = []
    for g in genes:
        gw_loo = {k: v for k, v in gw_dict.items() if k != g}
        loo_bias = HumanCT_AvgZ_Weighted(spec_mat, gw_loo)
        common_cts = full_bias.index.intersection(loo_bias.index)
        rho, _ = spearmanr(full_bias.loc[common_cts, "EFFECT"].values,
                           loo_bias.loc[common_cts, "EFFECT"].values)
        rhos.append(rho)
        loo_biases.append(loo_bias)
    if label:
        print(f"  {label}: {len(genes)} LOO iterations")
    return {
        "rhos": np.array(rhos),
        "genes": genes,
        "full_bias": full_bias,
        "loo_biases": loo_biases,
    }


loo_results = {}
for disorder, gw in loo_gw.items():
    loo_results[disorder] = {}
    for spec_label, spec_mat in [("Capped (2×)", our_spec), ("Uncapped", spec_uncap)]:
        key = f"{disorder} | {spec_label}"
        print(f"Running LOO: {key}")
        loo_results[disorder][spec_label] = run_loo_analysis(spec_mat, gw, label=key)

# %%
# LOO summary
print("\nLOO Stability Summary (Spearman ρ vs full gene set)")
print("=" * 70)
for disorder in loo_gw:
    for spec_label in ["Capped (2×)", "Uncapped"]:
        rhos = loo_results[disorder][spec_label]["rhos"]
        print(f"  {disorder:12s} | {spec_label:12s} | "
              f"min={rhos.min():.4f}  median={np.median(rhos):.4f}  "
              f"mean={rhos.mean():.4f}  max={rhos.max():.4f}")

# %% [markdown]
# ---
# ## Assemble Figure S4

# %%
fig, axes = plt.subplots(3, 2, figsize=(14, 16), dpi=150)

# ============================================================
# Panel A: Empirical specificity inflation vs UMI depth
# ============================================================
ax = axes[0, 0]

ax.scatter(ct_stats.loc[neur_mask, "Total_UMI"],
           ct_stats.loc[neur_mask, "frac_clipped"],
           color="#D04040", alpha=0.5, s=20, edgecolors="white", lw=0.3,
           label="Neuronal", zorder=3)
ax.scatter(ct_stats.loc[~neur_mask, "Total_UMI"],
           ct_stats.loc[~neur_mask, "frac_clipped"],
           color="#3070B0", alpha=0.6, s=20, edgecolors="white", lw=0.3,
           label="Non-neuronal", zorder=4)

ax.axvline(x=10_000, color="gray", ls=":", lw=1, alpha=0.5, zorder=1)
ax.axvline(x=30_000, color="gray", ls=":", lw=1, alpha=0.5, zorder=1)

ax.set_xscale("log")
ax.set_xlabel("Library size (total UMI)")
ax.set_ylabel("Fraction of genes with specificity > 2")
ax.legend(fontsize=9, framealpha=0.8, loc="upper right")
ax.text(0.03, 0.97, f"ρ = {rho_umi:.3f}\np = {p_umi:.1e}",
        transform=ax.transAxes, ha="left", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
ax.set_title("A", fontweight="bold", loc="left", fontsize=16, pad=8)

# ============================================================
# Panel B: NB simulation — max specificity by expr × UMI group
# ============================================================
ax = axes[0, 1]

umi_group_colors = ["#dc2626", "#f59e0b", "#2563eb"]
for g_idx in range(3):
    ax.plot(expression_percentiles, panelB_data[g_idx],
            color=umi_group_colors[g_idx], marker="o", markersize=5, lw=2.2,
            label=umi_group_names[g_idx], zorder=5)

ax.axhline(y=1.0, color="gray", ls="--", lw=1.2, alpha=0.7, label="True specificity (1.0)")
ax.axhline(y=clip_threshold, color="black", ls=":", lw=1.2, alpha=0.7,
           label=f"Cap threshold ({clip_threshold:.1f})")

ax.set_yscale("log")
ax.set_xlabel("Gene expression level (percentile)")
ax.set_ylabel("Max simulated specificity")
ax.legend(fontsize=8.5, framealpha=0.8, loc="upper right")
ax.set_title("B", fontweight="bold", loc="left", fontsize=16, pad=8)

# ============================================================
# Panel C: Spearman ρ across cap levels
# ============================================================
ax = axes[1, 0]

disorder_colors = {"SCZ": "#2563eb", "ASD (HIQ)": "#dc2626", "DDD": "#16a34a"}
disorder_markers = {"SCZ": "o", "ASD (HIQ)": "s", "DDD": "D"}

# Shade optimal range
ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)

for disorder in corr_df.columns:
    rho_vals = corr_df[disorder].values
    ax.plot(CAP_LEVELS, rho_vals,
            color=disorder_colors[disorder],
            marker=disorder_markers[disorder],
            markersize=6, lw=2, label=disorder, zorder=3)

ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
cap_tick_labels = [f"{c:.0f}×" if c == int(c) else f"{c}×" for c in CAP_LEVELS]
cap_tick_labels[-1] = "No cap"  # 100× is effectively uncapped
ax.set_xticks(CAP_LEVELS)
ax.set_xticklabels(cap_tick_labels, fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("Specificity cap level (× mean)")
ax.set_ylabel("Spearman ρ vs. cap = 2×")
ax.set_ylim(min(corr_df.values.min() - 0.05, 0.35), 1.02)
ax.legend(fontsize=9, frameon=True, framealpha=0.9, edgecolor="none",
          loc="lower left", handlelength=1.5)
ax.set_title("C", fontweight="bold", loc="left", fontsize=16, pad=8)

# ============================================================
# Panel D: Supercluster rank across cap levels (SCZ)
# ============================================================
ax = axes[1, 1]

sc_style = {
    "CGE interneuron":      {"color": "#e11d48", "marker": "o", "ls": "-",  "lw": 2.2},
    "MGE interneuron":      {"color": "#7c3aed", "marker": "s", "ls": "-",  "lw": 2.2},
    "Medium spiny neuron":  {"color": "#0891b2", "marker": "D", "ls": "-",  "lw": 2.2},
    "Vascular":             {"color": "#64748b", "marker": "v", "ls": "--", "lw": 1.5},
    "Microglia":            {"color": "#92400e", "marker": "^", "ls": "--", "lw": 1.5},
    "Astrocyte":            {"color": "#d97706", "marker": "P", "ls": "--", "lw": 1.5},
    "Ependymal":            {"color": "#6b7280", "marker": "X", "ls": "--", "lw": 1.5},
}
sc_short = {
    "CGE interneuron": "CGE IN", "MGE interneuron": "MGE IN",
    "Medium spiny neuron": "MSN", "Vascular": "Vascular",
    "Microglia": "Microglia", "Astrocyte": "Astrocyte", "Ependymal": "Ependymal",
}

# Shade optimal range
ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)

neuronal_scs = [sc for sc in panel_scs if sc_style[sc]["ls"] == "-"]
nonneuronal_scs = [sc for sc in panel_scs if sc_style[sc]["ls"] == "--"]

for sc in neuronal_scs + nonneuronal_scs:
    sty = sc_style[sc]
    ranks = [sc_rank_data[sc][cap] for cap in CAP_LEVELS]
    ax.plot(CAP_LEVELS, ranks,
            color=sty["color"], marker=sty["marker"], markersize=5,
            lw=sty["lw"], ls=sty["ls"], label=sc_short[sc], zorder=3)

ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
ax.set_xticks(CAP_LEVELS)
ax.set_xticklabels(cap_tick_labels, fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("Specificity cap level (× mean)")
ax.set_ylabel("Supercluster rank (1 = highest bias)")
ax.invert_yaxis()
ax.set_title("D", fontweight="bold", loc="left", fontsize=16, pad=8)

leg = ax.legend(fontsize=8, frameon=True, framealpha=0.9, edgecolor="none",
                loc="lower right", ncol=2, columnspacing=0.8, handlelength=1.5,
                title="Neuronal          Non-neuronal", title_fontsize=8)
leg._legend_box.align = "left"

# ============================================================
# Panel E: Neuronal fraction — capped vs TDEP uncapped
# ============================================================
ax = axes[2, 0]

DISORDERS_E = ["ASD", "SCZ", "DDD"]
disorder_colors_e = {"ASD": "#2166AC", "SCZ": "#B2182B", "DDD": "#1B7837"}

# Small jitter for overlapping lines at 1.0
_jitter_map = {"ASD": 0.0, "SCZ": -0.012, "DDD": 0.012}
for disorder in DISORDERS_E:
    c = disorder_colors_e[disorder]
    jit = _jitter_map[disorder]
    fracs_cap = panelE_data[disorder]["capped"]
    fracs_tdep = panelE_data[disorder]["tdep"]
    fracs_cap_jit = [f + jit if f >= 0.99 else f for f in fracs_cap]
    ax.plot(cutoffs, fracs_cap_jit, "o-", color=c, linewidth=2, markersize=5,
            label=f"{disorder} (capped)", zorder=5)
    ax.plot(cutoffs, fracs_tdep, "s--", color=c, linewidth=1.5, markersize=4,
            alpha=0.6, label=f"{disorder} (TDEP)", zorder=3)

expected_frac = len(Neur_idx) / (len(Neur_idx) + len(NonNeur_idx))
ax.axhline(expected_frac, color="gray", linestyle=":", linewidth=1.2, zorder=1,
           label=f"Expected ({expected_frac:.0%})")

ax.set_xlabel("Top-N cutoff")
ax.set_ylabel("Fraction neuronal")
ax.set_ylim(-0.05, 1.08)
ax.set_xticks(cutoffs)
ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9, ncol=2)
ax.set_title("E", fontweight="bold", loc="left", fontsize=16, pad=8)

# ============================================================
# Panel F: LOO stability violins
# ============================================================
ax = axes[2, 1]

cap_color = "#2166AC"
uncap_color = "#B2182B"

loo_disorders = ["SCZ", "ASD"]
positions = []
all_rhos = []
all_colors = []
xtick_labels = []
pos = 0
for disorder in loo_disorders:
    for spec_label, color in [("Capped (2×)", cap_color), ("Uncapped", uncap_color)]:
        rhos = loo_results[disorder][spec_label]["rhos"]
        positions.append(pos)
        all_rhos.append(rhos)
        all_colors.append(color)
        short_d = "SCZ" if disorder == "SCZ" else "ASD"
        short_s = "Cap" if "Cap" in spec_label else "Uncap"
        xtick_labels.append(f"{short_d}\n{short_s}")
        pos += 1
    pos += 0.5

parts = ax.violinplot(all_rhos, positions=positions, showmeans=False,
                      showmedians=True, widths=0.7)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(all_colors[i])
    body.set_alpha(0.3)
for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
    parts[partname].set_edgecolor("gray")
    parts[partname].set_linewidth(1)
parts["cmedians"].set_edgecolor("black")
parts["cmedians"].set_linewidth(1.5)

np.random.seed(42)
for i, (pos_i, rhos, color) in enumerate(zip(positions, all_rhos, all_colors)):
    jitter = np.random.normal(0, 0.05, len(rhos))
    ax.scatter(pos_i + jitter, rhos, color=color, alpha=0.5, s=8, zorder=3)
    ax.text(pos_i, rhos.min() - 0.003, f"{rhos.min():.3f}",
            ha="center", va="top", fontsize=7, color=color, fontweight="bold")

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
ax.set_xticks(positions)
ax.set_xticklabels(xtick_labels, fontsize=9)
ax.set_ylabel("Spearman ρ (LOO vs full)")

all_rhos_flat = np.concatenate(all_rhos)
ymin = min(all_rhos_flat.min() - 0.015, 0.90)
ax.set_ylim(ymin, 1.005)

legend_elems = [
    Patch(facecolor=cap_color, alpha=0.5, label="Capped (2×)"),
    Patch(facecolor=uncap_color, alpha=0.5, label="Uncapped"),
]
ax.legend(handles=legend_elems, fontsize=9, loc="lower left", framealpha=0.8)
ax.set_title("F", fontweight="bold", loc="left", fontsize=16, pad=8)

# ============================================================
# Save
# ============================================================
fig.tight_layout(h_pad=2.5, w_pad=2.0)

fig.savefig(FIG_DIR / "FigS4_specificity_cap.pdf",
            dpi=300, transparent=True, bbox_inches="tight")
fig.savefig(FIG_DIR / "FigS4_specificity_cap.png",
            dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"\nFigure S4 saved to {FIG_DIR}")

# %% [markdown]
# ---
# ## Summary Statistics for Supplementary Note

# %%
print("=" * 70)
print("SUMMARY STATISTICS FOR SUPPLEMENTARY NOTE 1")
print("=" * 70)

print("\n--- Section 1.1: Empirical inflation ---")
print(f"Max unclipped specificity: {spec_unclip.values.max():.1f}×")
print(f"Clip threshold (2× mean): {clip_threshold:.4f}")
total_vals = spec_unclip.values.size
total_clipped = np.sum(spec_unclip.values > clip_threshold)
print(f"Values exceeding cap: {total_clipped:,} / {total_vals:,} ({100*total_clipped/total_vals:.2f}%)")
print(f"Spearman(UMI depth, frac_clipped): ρ = {rho_umi:.3f}, p = {p_umi:.1e}")

neur_frac_mean = ct_stats.loc[neur_mask, "frac_clipped"].mean()
nonneur_frac_mean = ct_stats.loc[~neur_mask, "frac_clipped"].mean()
print(f"Mean frac clipped — Neuronal: {neur_frac_mean:.2%}, Non-neuronal: {nonneur_frac_mean:.2%}")

print("\n--- Section 1.2: NB simulation ---")
print(f"Global theta: {global_theta:.3f}")
for pct in expression_percentiles:
    res = sweep_results_df[pct]
    print(f"  P{pct}: max_spec = {res['max_spec'].max():.1f}×")

print("\n--- Section 1.3: Cap sensitivity ---")
for disorder in corr_df.columns:
    rhos_range = [corr_results[disorder][c] for c in [1.5, 2.5, 3.0]]
    print(f"{disorder}: ρ at 1.5–3× = {min(rhos_range):.4f}–{max(rhos_range):.4f}, "
          f"at 1× = {corr_results[disorder][1.0]:.4f}, "
          f"at 10× = {corr_results[disorder][10.0]:.4f}")

print("\n--- Section 1.4: LOO stability ---")
for disorder in loo_disorders:
    for spec_label in ["Capped (2×)", "Uncapped"]:
        rhos = loo_results[disorder][spec_label]["rhos"]
        print(f"  {disorder} {spec_label}: min ρ = {rhos.min():.4f}, "
              f"median = {np.median(rhos):.4f}")

# %% [markdown]
# ---
# ## Bonus: Worst-Case Gene Rank Scatter (not in Figure S4)
#
# Shows that removing the single most influential gene from the uncapped specificity
# causes large rank rearrangements across cell types. Each dot is one cell type;
# deviation from the diagonal = rank change when that gene is removed.
# This illustrates how one gene with extreme uncapped specificity can dominate
# the entire bias ranking.

# %%
HGNC, _, _, Entrez2Symbol = LoadGeneINFO()

NEUR_COLOR = "#D04040"
NONNEUR_COLOR = "#3070B0"

fig_bonus, axes_bonus = plt.subplots(1, 2, figsize=(14, 6), dpi=150)

for ax_idx, disorder_b in enumerate(loo_disorders):
    ax = axes_bonus[ax_idx]

    loo_uncap = loo_results[disorder_b]["Uncapped"]
    rhos_uncap = loo_uncap["rhos"]
    worst_idx = np.argmin(rhos_uncap)
    worst_gene = loo_uncap["genes"][worst_idx]
    worst_rho = rhos_uncap[worst_idx]
    worst_loo_bias = loo_uncap["loo_biases"][worst_idx]
    full_bias_uncap = loo_uncap["full_bias"]

    worst_symbol = Entrez2Symbol.get(worst_gene, Entrez2Symbol.get(int(worst_gene), str(worst_gene)))

    # Check uncapped specificity of this gene
    gene_id = int(worst_gene) if int(worst_gene) in spec_unclip.index else worst_gene
    if gene_id in spec_unclip.index:
        max_spec_gene = spec_unclip.loc[gene_id].max()
        max_ct = spec_unclip.loc[gene_id].idxmax()
        max_ct_name = Anno.loc[int(max_ct), "Supercluster"] if int(max_ct) in Anno.index else "?"
    else:
        max_spec_gene = np.nan
        max_ct_name = "?"

    # Compute ranks (align by cell type index)
    common_cts = full_bias_uncap.index.intersection(worst_loo_bias.index)
    rank_full = full_bias_uncap.loc[common_cts, "EFFECT"].rank(ascending=False)
    rank_loo = worst_loo_bias.loc[common_cts, "EFFECT"].rank(ascending=False)

    neur_mask_b = np.array([int(ct) in Neur_idx for ct in common_cts])

    ax.scatter(rank_full.values[neur_mask_b], rank_loo.values[neur_mask_b],
               color=NEUR_COLOR, alpha=0.4, s=18, label="Neuronal", zorder=3)
    ax.scatter(rank_full.values[~neur_mask_b], rank_loo.values[~neur_mask_b],
               color=NONNEUR_COLOR, alpha=0.5, s=18, label="Non-neuronal", zorder=4)

    max_rank = max(rank_full.max(), rank_loo.max())
    ax.plot([1, max_rank], [1, max_rank], "k--", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("Rank (all genes)")
    ax.set_ylabel(f"Rank (remove {worst_symbol})")
    ax.set_title(
        f"{disorder_b} uncapped: remove {worst_symbol}\n"
        f"(max spec = {max_spec_gene:.1f}× in {max_ct_name}, "
        f"LOO ρ = {worst_rho:.3f})",
        fontsize=10,
    )
    ax.legend(fontsize=9, loc="upper left", framealpha=0.8)

    print(f"{disorder_b}: worst gene = {worst_symbol} (Entrez {worst_gene}), "
          f"max uncapped spec = {max_spec_gene:.1f}×, LOO ρ = {worst_rho:.4f}")

fig_bonus.tight_layout(w_pad=3)
fig_bonus.savefig(FIG_DIR / "worst_gene_rank_scatter.pdf",
                  dpi=300, transparent=True, bbox_inches="tight")
fig_bonus.savefig(FIG_DIR / "worst_gene_rank_scatter.png",
                  dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR / 'worst_gene_rank_scatter.pdf'}")
