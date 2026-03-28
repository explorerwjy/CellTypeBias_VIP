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
# # TDEP-sLDSC Specificity vs Capped Specificity: Bias Comparison
#
# The TDEP-sLDSC pipeline (Bryois et al., Nature Communications 2020) defines
# cell-type specificity as a **proportion**: $S_{g,ct} = \text{TPM}_{g,ct} / \sum_{ct'} \text{TPM}_{g,ct'}$.
# Values sum to 1.0 per gene and are used only in a **rank-based** top-10% framework,
# never as continuous multipliers.
#
# Our pipeline uses **fold-enrichment** specificity capped at 2x the mean,
# then mean-centered. Values are used as **continuous weights** in a weighted sum:
# $\text{Bias}(ct) = \sum_g w_g \cdot S_{g,ct} / \sum_g w_g$.
#
# This notebook shows that importing TDEP's uncapped proportion-based
# specificity into our continuous weighting framework produces results
# dominated by non-neuronal cell types -- because rare, low-UMI cell types
# accumulate inflated specificity for a handful of detected genes,
# and continuous weighting amplifies those extremes.

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## Setup

# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import (
    HumanCT_AvgZ_Weighted,
    AnnotateCTDat,
    Fil2Dict,
    Neur_idx,
    NonNeur_idx,
    Anno,
)

FIG_DIR = PROJ_DIR / "results" / "figures" / "specificity_cap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

font_path = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
if Path(font_path).exists():
    fm.fontManager.addfont(font_path)
    fm._load_fontmanager(try_read_cache=False)

mpl.rcParams["figure.facecolor"] = "none"
mpl.rcParams["axes.facecolor"] = "none"
mpl.rcParams["savefig.facecolor"] = "none"
mpl.rcParams["font.size"] = 11
mpl.rcParams["font.family"] = "Arial"

# %% [markdown]
# ## Load Data

# %%
# --- TDEP specificity (proportion-based, uncapped) ---
tdep_spec = pd.read_csv(
    "/home/jw3514/Work/NeuralP/TDEP-sLDSC/data/cluster.specificity_matrix_entrez.csv",
    index_col=0,
)
# Fix column IDs: TDEP uses 1-based, we use 0-based
tdep_spec.columns = [int(c) - 1 for c in tdep_spec.columns]
print(f"TDEP specificity: {tdep_spec.shape[0]} genes x {tdep_spec.shape[1]} cell types")
print(f"  Row sums (should be ~1.0): {tdep_spec.iloc[:3].sum(axis=1).values}")
print(f"  Value range: [{tdep_spec.values.min():.6f}, {tdep_spec.values.max():.6f}]")

# %%
# --- Our capped, mean-centered specificity ---
our_spec = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    index_col=0,
)
our_spec.columns = [int(c) for c in our_spec.columns]
print(f"Our specificity: {our_spec.shape[0]} genes x {our_spec.shape[1]} cell types")
print(f"  Value range: [{our_spec.values.min():.4f}, {our_spec.values.max():.4f}]")

# %%
# --- Load annotation ---
anno = pd.read_excel(PROJ_DIR / "dat" / "annotation.xlsx", index_col=0)
anno = anno[anno.index.notna()]
anno.index = anno.index.astype(int)

# %% [markdown]
# ## Process TDEP Specificity for Our Framework
#
# The TDEP proportions sum to 1.0 per gene. To make them comparable to
# fold-enrichment (where the mean across cell types = 1.0), we multiply
# by the number of cell types (461). Then we mean-center per gene,
# exactly as our pipeline does -- **but without the 2x cap**.

# %%
n_cts = tdep_spec.shape[1]
print(f"Number of cell types: {n_cts}")

# Scale proportions to fold-enrichment scale (mean = 1.0 per gene)
tdep_fold = tdep_spec * n_cts

# Mean-center per gene (subtract row mean from each value)
tdep_mc = tdep_fold.sub(tdep_fold.mean(axis=1), axis=0)

print(f"After scaling to fold-enrichment:")
print(f"  Value range: [{tdep_fold.values.min():.2f}, {tdep_fold.values.max():.2f}]")
print(f"  Max fold-enrichment: {tdep_fold.values.max():.1f}x (ours caps at 2x)")
print(f"After mean-centering:")
print(f"  Value range: [{tdep_mc.values.min():.4f}, {tdep_mc.values.max():.4f}]")

# %% [markdown]
# ## Load Gene Weights

# %%
gw_files = {
    "ASD w/o ID": PROJ_DIR / "dat" / "GeneWeights" / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw",
    "SCZ": PROJ_DIR / "dat" / "GeneWeights" / "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    "DDD": PROJ_DIR / "dat" / "GeneWeights" / "DDD.top61.gw.bgmr.csv",
}

gene_weights = {}
for label, path in gw_files.items():
    gw = Fil2Dict(str(path))
    gene_weights[label] = gw
    print(f"{label}: {len(gw)} genes")

# %% [markdown]
# ## Compute Bias Scores

# %%
results_our = {}
results_tdep = {}

for label, gw in gene_weights.items():
    # Our capped + mean-centered specificity
    bias_our = HumanCT_AvgZ_Weighted(our_spec, gw)
    bias_our = AnnotateCTDat(bias_our, anno)
    bias_our["is_neuronal"] = bias_our.index.isin(Neur_idx)
    results_our[label] = bias_our

    # TDEP uncapped + mean-centered specificity
    bias_tdep = HumanCT_AvgZ_Weighted(tdep_mc, gw)
    bias_tdep = AnnotateCTDat(bias_tdep, anno)
    bias_tdep["is_neuronal"] = bias_tdep.index.isin(Neur_idx)
    results_tdep[label] = bias_tdep

    # Quick summary
    top20_our_neur = bias_our.head(20)["is_neuronal"].sum()
    top20_tdep_neur = bias_tdep.head(20)["is_neuronal"].sum()
    print(f"{label}: Top-20 neuronal — Ours: {top20_our_neur}/20, TDEP: {top20_tdep_neur}/20")

# %% [markdown]
# ## Summary Statistics

# %%
print("=" * 72)
print("SUMMARY: Neuronal fraction in top-ranked cell types")
print("=" * 72)
n_neur_total = len(Neur_idx)
n_nonneur_total = len(NonNeur_idx)
expected_frac = n_neur_total / (n_neur_total + n_nonneur_total)
print(f"Expected by chance: {n_neur_total}/{n_neur_total + n_nonneur_total} = {expected_frac:.2%}")
print()

for cutoff in [5, 10, 15, 20, 30, 50]:
    print(f"--- Top {cutoff} ---")
    for label in gene_weights:
        n_our = results_our[label].head(cutoff)["is_neuronal"].sum()
        n_tdep = results_tdep[label].head(cutoff)["is_neuronal"].sum()
        print(f"  {label:12s}  Ours: {n_our:2d}/{cutoff} ({n_our/cutoff:.0%})  "
              f"TDEP: {n_tdep:2d}/{cutoff} ({n_tdep/cutoff:.0%})")
    print()

# %%
# Show the non-neuronal interlopers in TDEP top 20 for ASD w/o ID
print("TDEP top-20 for ASD w/o ID — non-neuronal cell types highlighted:")
top20 = results_tdep["ASD w/o ID"].head(20)
for i, (ct_idx, row) in enumerate(top20.iterrows()):
    marker = "  <<<" if not row["is_neuronal"] else ""
    sc_short = str(row["Supercluster"])[:25]
    cls = str(row["Class"]) if pd.notna(row["Class"]) else "N/A"
    print(f"  Rank {i+1:2d}: ct={ct_idx:3d}  {sc_short:25s}  "
          f"Class={cls:5s}  EFFECT={row['EFFECT']:.4f}{marker}")

# %% [markdown]
# ## Combined Figure: TDEP Comparison + LOO Stability (MC2d)
#
# Three panels:
# - **A**: Rank scatter (Ours vs TDEP) for ASD w/o ID — non-neuronal types
#   jump to top ranks under TDEP uncapped specificity
# - **B**: Neuronal fraction vs top-N cutoff across all 3 disorders
# - **C**: LOO stability violins (SCZ + ASD, capped vs uncapped)

# %%
# (Combined figure is drawn after LOO computation below)

# %% [markdown]
# ## Conclusion (TDEP comparison)
#
# Using the TDEP-sLDSC proportion-based specificity (uncapped) in our
# continuous-weighting framework produces dramatically different results:
#
# - **ASD w/o ID**: Non-neuronal cell types (B-cells, T-cells, microglia)
#   dominate the top rankings, displacing the CGE interneurons found by
#   our capped metric.
# - **SCZ and DDD**: Less affected because their gene sets happen to have
#   less overlap with the inflated-specificity genes, but still show
#   more non-neuronal contamination than our capped metric.
#
# The key insight: TDEP's specificity metric was designed for a
# **rank-based** top-10% thresholding approach, where extreme values
# are irrelevant once a gene is classified as "specific" or "not specific."
# Our method uses specificity as **continuous multipliers**, so extreme
# values in rare, low-UMI cell types can dominate the weighted sum.
# The 2x cap prevents this artifact.

# %% [markdown]
# ---
# # Leave-One-Out Gene Stability Test: Capped vs Uncapped
#
# **Key insight:** With capped specificity (cap = 2x mean), removing any single gene
# from the risk gene set barely changes the bias ranking. Without capping, a single
# gene with extreme specificity (e.g., 21x in Purkinje neurons) can dominate
# the entire result, making the bias ranking unstable.
#
# The bias is a weighted mean of specificity:
# `bias(ct) = sum(weight * specificity) / sum(weight)`
#
# Without capping, a gene with specificity 100x in one cell type contributes 50x more
# than a gene at 2x. Removing that ONE gene can completely change which cell types rank
# highest. With capping at 2x, no single gene can contribute more than 2x, so removing
# any gene has a small effect.

# %% [markdown]
# ## Load Uncapped Specificity

# %%
from scipy.stats import spearmanr
import seaborn as sns
from CellType_PSY import LoadGeneINFO

# Load unclipped specificity (clip100 is effectively unclipped)
spec_unclip = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv",
    index_col=0,
)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]

# Mean-center per gene (same preprocessing as our capped pipeline, but without cap)
spec_uncap = spec_unclip.subtract(spec_unclip.mean(axis=1), axis=0)
print(f"Uncapped specificity: {spec_uncap.shape[0]} genes x {spec_uncap.shape[1]} cell types")
print(f"  Max unclipped value: {spec_unclip.values.max():.1f}x")

# Gene symbol mapping for identifying dominant genes
HGNC, _, _, Entrez2Symbol = LoadGeneINFO()

# %% [markdown]
# ## LOO Computation
#
# For each disorder and each specificity matrix (capped vs uncapped):
# 1. Compute full bias with all genes
# 2. For each gene, remove it and recompute bias
# 3. Compute Spearman rho between full and LOO bias rankings

# %%
# Use SCZ and ASD w/o ID for LOO (same gene weights as above)
loo_gw_dicts = {k: v for k, v in gene_weights.items() if k in ["SCZ", "ASD w/o ID"]}


def run_loo_analysis(spec_mat, gw_dict, label=""):
    """Run leave-one-out analysis for a single disorder x specificity matrix."""
    full_bias = HumanCT_AvgZ_Weighted(spec_mat, gw_dict)
    loo_results = []
    genes = list(gw_dict.keys())
    for g in genes:
        gw_loo = {k: v for k, v in gw_dict.items() if k != g}
        loo_bias = HumanCT_AvgZ_Weighted(spec_mat, gw_loo)
        common_cts = full_bias.index.intersection(loo_bias.index)
        rho, _ = spearmanr(
            full_bias.loc[common_cts, "EFFECT"].values,
            loo_bias.loc[common_cts, "EFFECT"].values,
        )
        loo_results.append({"gene_id": g, "rho": rho, "loo_bias": loo_bias})
    if label:
        print(f"  {label}: {len(genes)} LOO iterations completed")
    return full_bias, loo_results


# %%
# Run LOO for all combinations
# our_spec = capped (already loaded above), spec_uncap = uncapped
loo_results = {}
for disorder, gw_dict in loo_gw_dicts.items():
    loo_results[disorder] = {}
    for spec_label, spec_mat in [("Capped (2x)", our_spec), ("Uncapped", spec_uncap)]:
        key = f"{disorder} | {spec_label}"
        print(f"Running LOO: {key}")
        full_bias, loo = run_loo_analysis(spec_mat, gw_dict, label=key)
        loo_results[disorder][spec_label] = {
            "full_bias": full_bias,
            "loo": loo,
            "rhos": np.array([r["rho"] for r in loo]),
        }

# %% [markdown]
# ## LOO Summary Statistics

# %%
print("LOO Stability Summary (Spearman rho vs full gene set)")
print("=" * 70)
for disorder in loo_gw_dicts:
    for spec_label in ["Capped (2x)", "Uncapped"]:
        rhos = loo_results[disorder][spec_label]["rhos"]
        print(f"  {disorder:12s} | {spec_label:12s} | "
              f"min={rhos.min():.4f}  median={np.median(rhos):.4f}  "
              f"max={rhos.max():.4f}  mean={rhos.mean():.4f}")
print("=" * 70)

# %% [markdown]
# ## Complete Combined Figure (Panel C: LOO violins)

# %%
# Combined 1x3 figure: A=Neuronal fraction, B=LOO violins, C=Worst-case gene rank scatter
from matplotlib.patches import Patch
from scipy.stats import spearmanr as _spearmanr

NEUR_COLOR = "#D04040"
NONNEUR_COLOR = "#3070B0"
cap_color = "#2166AC"
uncap_color = "#B2182B"
DISORDERS = ["ASD w/o ID", "SCZ", "DDD"]
disorder_colors = {"ASD w/o ID": "#2166AC", "SCZ": "#B2182B", "DDD": "#1B7837"}

fig, axes = plt.subplots(1, 3, figsize=(21, 6), dpi=150)

# ── Panel A: Neuronal fraction vs top-N ──
ax_a = axes[0]
cutoffs = [5, 10, 15, 20, 30, 50]

# Print the raw fractions to debug SCZ capped visibility
for disorder in DISORDERS:
    fracs_ours = []
    fracs_tdep = []
    for k in cutoffs:
        n_o = results_our[disorder].head(k)["is_neuronal"].sum()
        n_t = results_tdep[disorder].head(k)["is_neuronal"].sum()
        fracs_ours.append(n_o / k)
        fracs_tdep.append(n_t / k)
    print(f"{disorder} capped:  {fracs_ours}")
    print(f"{disorder} TDEP:    {fracs_tdep}")

# Plot with slight y-offset for overlapping capped lines at 1.0
_jitter_map = {"ASD w/o ID": 0.0, "SCZ": -0.012, "DDD": 0.012}
for disorder in DISORDERS:
    fracs_ours = []
    fracs_tdep = []
    for k in cutoffs:
        n_o = results_our[disorder].head(k)["is_neuronal"].sum()
        n_t = results_tdep[disorder].head(k)["is_neuronal"].sum()
        fracs_ours.append(n_o / k)
        fracs_tdep.append(n_t / k)

    c = disorder_colors[disorder]
    jit = _jitter_map[disorder]
    fracs_ours_jit = [f + jit if f >= 0.99 else f for f in fracs_ours]
    ax_a.plot(cutoffs, fracs_ours_jit, "o-", color=c, linewidth=2.2, markersize=7,
              label=f"{disorder} (capped)", zorder=5)
    ax_a.plot(cutoffs, fracs_tdep, "s--", color=c, linewidth=2, markersize=6,
              alpha=0.6, label=f"{disorder} (TDEP)", zorder=3)

expected_frac = len(Neur_idx) / (len(Neur_idx) + len(NonNeur_idx))
ax_a.axhline(expected_frac, color="gray", linestyle=":", linewidth=1.5, zorder=1,
             label=f"Expected ({expected_frac:.0%})")

ax_a.set_xlabel("Top-N cutoff", fontsize=11)
ax_a.set_ylabel("Fraction neuronal", fontsize=11)
ax_a.set_ylim(-0.05, 1.05)
ax_a.set_xticks(cutoffs)
ax_a.legend(fontsize=7.5, loc="lower right", framealpha=0.9, ncol=2)
ax_a.set_title("A", fontweight="bold", loc="left", fontsize=14)
ax_a.set_title("Neuronal enrichment in top-ranked", fontsize=11, loc="center")
for sp in ["top", "right"]:
    ax_a.spines[sp].set_visible(False)

# ── Panel B: LOO stability violins ──
ax_b = axes[1]

loo_disorders = ["SCZ", "ASD w/o ID"]
positions = []
all_rhos = []
all_colors = []
xtick_labels = []
pos = 0
for d_i, disorder in enumerate(loo_disorders):
    for s_i, (spec_label, color) in enumerate([("Capped (2x)", cap_color), ("Uncapped", uncap_color)]):
        rhos = loo_results[disorder][spec_label]["rhos"]
        positions.append(pos)
        all_rhos.append(rhos)
        all_colors.append(color)
        short_d = "SCZ" if disorder == "SCZ" else "ASD"
        short_s = "Cap" if "Cap" in spec_label else "Uncap"
        xtick_labels.append(f"{short_d}\n{short_s}")
        pos += 1
    pos += 0.5

parts = ax_b.violinplot(all_rhos, positions=positions, showmeans=False,
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
    ax_b.scatter(pos_i + jitter, rhos, color=color, alpha=0.5, s=12, zorder=3)
    ax_b.text(pos_i, rhos.min() - 0.003, f"{rhos.min():.3f}",
              ha="center", va="top", fontsize=7.5, color=color, fontweight="bold")

ax_b.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
ax_b.set_xticks(positions)
ax_b.set_xticklabels(xtick_labels, fontsize=9)
ax_b.set_ylabel("Spearman $\\rho$ (LOO vs full)", fontsize=11)

all_rhos_flat = np.concatenate(all_rhos)
ymin = min(all_rhos_flat.min() - 0.015, 0.90)
ax_b.set_ylim(ymin, 1.005)

legend_elems = [
    Patch(facecolor=cap_color, alpha=0.5, label="Capped (2×)"),
    Patch(facecolor=uncap_color, alpha=0.5, label="Uncapped"),
]
ax_b.legend(handles=legend_elems, fontsize=9, loc="lower left", framealpha=0.8)
ax_b.set_title("B", fontweight="bold", loc="left", fontsize=14)
ax_b.set_title("Leave-one-out stability", fontsize=11, loc="center")
for sp in ["top", "right"]:
    ax_b.spines[sp].set_visible(False)

# ── Panel C: Worst-case gene rank scatter (SCZ uncapped, CACNA1G) ──
ax_c = axes[2]
disorder_c = "SCZ"

loo_uncap_scz = loo_results[disorder_c]["Uncapped"]["loo"]
rhos_uncap_scz = loo_results[disorder_c]["Uncapped"]["rhos"]
worst_idx = np.argmin(rhos_uncap_scz)
worst_gene = loo_uncap_scz[worst_idx]["gene_id"]
worst_rho = rhos_uncap_scz[worst_idx]
worst_loo_bias = loo_uncap_scz[worst_idx]["loo_bias"]
full_bias_uncap = loo_results[disorder_c]["Uncapped"]["full_bias"]

worst_symbol = Entrez2Symbol.get(worst_gene, Entrez2Symbol.get(int(worst_gene), str(worst_gene)))

common_cts_c = full_bias_uncap.index.intersection(worst_loo_bias.index)
rank_full = full_bias_uncap.loc[common_cts_c, "EFFECT"].rank(ascending=False)
rank_loo = worst_loo_bias.loc[common_cts_c, "EFFECT"].rank(ascending=False)

neur_mask_c = np.array([int(ct) in Neur_idx for ct in common_cts_c])

ax_c.scatter(rank_full.values[neur_mask_c], rank_loo.values[neur_mask_c],
             color=NEUR_COLOR, alpha=0.4, s=14, label="Neuronal", zorder=3)
ax_c.scatter(rank_full.values[~neur_mask_c], rank_loo.values[~neur_mask_c],
             color=NONNEUR_COLOR, alpha=0.4, s=14, label="Non-neuronal", zorder=3)

max_rank_c = max(rank_full.max(), rank_loo.max())
ax_c.plot([1, max_rank_c], [1, max_rank_c], "k--", linewidth=0.8, alpha=0.5)

ax_c.set_xlabel("Rank (all genes)", fontsize=11)
ax_c.set_ylabel(f"Rank (remove *{worst_symbol}*)", fontsize=11)
ax_c.set_title("C", fontweight="bold", loc="left", fontsize=14)
ax_c.set_title(
    f"SCZ uncapped: remove {worst_symbol} ($\\rho$ = {worst_rho:.3f})",
    fontsize=10, loc="center",
)
ax_c.legend(fontsize=9, loc="upper left", framealpha=0.8)
for sp in ["top", "right"]:
    ax_c.spines[sp].set_visible(False)

plt.tight_layout(pad=2.0)

fig.savefig(FIG_DIR / "tdep_loo_combined.pdf", bbox_inches="tight", transparent=True, dpi=300)
fig.savefig(FIG_DIR / "tdep_loo_combined.png", bbox_inches="tight", transparent=True, dpi=300)
print(f"Combined figure saved to {FIG_DIR / 'tdep_loo_combined.pdf'}")
plt.show()

# %% [markdown]
# ## Identify Dominant Genes
#
# For each disorder, find the genes whose removal causes the biggest rank change
# in the uncapped case.

# %%
for disorder in loo_gw_dicts:
    print(f"\n{'=' * 70}")
    print(f"DOMINANT GENES — {disorder} (Uncapped)")
    print(f"{'=' * 70}")

    loo_uncap_d = loo_results[disorder]["Uncapped"]["loo"]
    rhos_d = loo_results[disorder]["Uncapped"]["rhos"]

    sorted_idx = np.argsort(rhos_d)

    print(f"\n{'Gene':>10s}  {'Symbol':>10s}  {'LOO rho':>8s}  "
          f"{'Max Spec':>9s}  {'Max Spec CT':>12s}  {'CT Name'}")
    print("-" * 90)

    for rank, idx in enumerate(sorted_idx[:10]):
        gene_id = loo_uncap_d[idx]["gene_id"]
        rho = rhos_d[idx]
        symbol = Entrez2Symbol.get(gene_id, Entrez2Symbol.get(int(gene_id), "?"))

        if gene_id in spec_unclip.index:
            gene_spec = spec_unclip.loc[gene_id]
        elif int(gene_id) in spec_unclip.index:
            gene_spec = spec_unclip.loc[int(gene_id)]
        else:
            gene_spec = pd.Series(dtype=float)

        if len(gene_spec) > 0:
            max_spec = gene_spec.max()
            max_ct = gene_spec.idxmax()
            ct_name = anno.loc[int(max_ct), "Subtype auto-annotation"] if int(max_ct) in anno.index else "?"
        else:
            max_spec = np.nan
            max_ct = "?"
            ct_name = "?"

        print(f"{gene_id:>10}  {symbol:>10s}  {rho:>8.4f}  "
              f"{max_spec:>9.2f}  {max_ct:>12}  {ct_name}")

# %%
# Compare capped vs uncapped for the worst-case gene
for disorder in loo_gw_dicts:
    print(f"\n--- {disorder} ---")

    rhos_uncap_d = loo_results[disorder]["Uncapped"]["rhos"]
    loo_uncap_d = loo_results[disorder]["Uncapped"]["loo"]
    worst_idx = np.argmin(rhos_uncap_d)
    worst_gene = loo_uncap_d[worst_idx]["gene_id"]
    worst_symbol = Entrez2Symbol.get(worst_gene, Entrez2Symbol.get(int(worst_gene), str(worst_gene)))

    loo_cap_d = loo_results[disorder]["Capped (2x)"]["loo"]
    rhos_cap_d = loo_results[disorder]["Capped (2x)"]["rhos"]
    cap_gene_ids = [r["gene_id"] for r in loo_cap_d]

    if worst_gene in cap_gene_ids:
        cap_idx = cap_gene_ids.index(worst_gene)
        cap_rho = rhos_cap_d[cap_idx]
    else:
        cap_rho = np.nan

    print(f"  Worst gene (uncapped): {worst_symbol} (Entrez {worst_gene})")
    print(f"  Uncapped LOO rho:  {rhos_uncap_d[worst_idx]:.4f}")
    print(f"  Capped LOO rho:    {cap_rho:.4f}")
    print(f"  --> Capping protects: rho improves by {cap_rho - rhos_uncap_d[worst_idx]:.4f}")

# %% [markdown]
# ## Overall Conclusion
#
# 1. **TDEP comparison**: Using the TDEP-sLDSC uncapped specificity in our
#    continuous-weighting framework produces results dominated by non-neuronal
#    cell types, especially for ASD w/o ID.
#
# 2. **LOO stability**: Capping at 2x ensures that no single gene can dominate
#    the bias ranking. The minimum LOO Spearman rho is consistently high
#    (>0.95) with capping, while uncapped specificity shows lower stability
#    (min rho ~0.91) due to genes like *CACNA1G* (21x in Purkinje) and
#    *TBR1* (11.7x in L6 IT) that can single-handedly shift rankings.
