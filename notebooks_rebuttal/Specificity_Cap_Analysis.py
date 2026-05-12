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
# # Specificity Cap Analysis: Why Clipping at 2× Mean Is Necessary
#
# Expression specificity is computed as fold-enrichment: $S_{g,ct} = TPM_{g,ct} \,/\, \overline{TPM}_g$.
# Without clipping, small cell types with low sequencing depth develop inflated
# specificity scores (up to ~97×) because sparse expression + TPM normalization
# concentrates signal in the few detected genes.
#
# This notebook demonstrates that the cap at 2× mean specificity:
# 1. Primarily affects low-UMI non-neuronal cell types
# 2. Has minimal impact on neuronal populations that drive all biological conclusions
# 3. Is a conservative guard against technical artifacts, not a data manipulation

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
import seaborn as sns
from scipy import stats

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import Anno, Neur_idx, NonNeur_idx

FIG_DIR = PROJ_DIR / "results" / "figures" / "specificity_cap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
if Path(font_path).exists():
    fm.fontManager.addfont(font_path)
    fm._load_fontmanager(try_read_cache=False)

mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'Arial'

# %% [markdown]
# ## Load Data

# %%
# Use clip100.0 as effectively unclipped (max observed ~97, well below 100)
spec_unclip = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv",
    index_col=0
)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]

clip_threshold = np.mean(spec_unclip.values.flatten()) * 2
n_genes, n_cts = spec_unclip.shape
print(f"Matrix: {n_genes} genes × {n_cts} cell types")
print(f"Clip threshold (2× mean specificity): {clip_threshold:.4f}")
print(f"Max unclipped specificity: {spec_unclip.values.max():.1f}")

# %% [markdown]
# ## Per-Cell-Type Statistics

# %%
ct_stats = pd.DataFrame(index=Anno.index)

for ct in ct_stats.index:
    vals = spec_unclip[ct].values
    ct_stats.loc[ct, "frac_clipped"] = np.mean(vals > clip_threshold)
    ct_stats.loc[ct, "n_clipped"] = np.sum(vals > clip_threshold)
    ct_stats.loc[ct, "max_spec"] = np.max(vals)
    ct_stats.loc[ct, "p95_spec"] = np.percentile(vals, 95)

ct_stats["Total_UMI"] = Anno["Total UMI"]
ct_stats["N_cells"] = Anno["Number of cells"]
ct_stats["Class"] = Anno["Class auto-annotation"]
ct_stats["Supercluster"] = Anno["Supercluster"]
ct_stats["is_neuronal"] = ct_stats.index.isin(Neur_idx)

# Ensure numeric
for col in ["frac_clipped", "n_clipped", "max_spec", "p95_spec", "Total_UMI", "N_cells"]:
    ct_stats[col] = pd.to_numeric(ct_stats[col])

neur_mask = ct_stats["is_neuronal"]
print(f"Cell types: {neur_mask.sum()} neuronal, {(~neur_mask).sum()} non-neuronal")

# %% [markdown]
# ## Summary Statistics

# %%
total_vals = spec_unclip.values.size
total_clipped = np.sum(spec_unclip.values > clip_threshold)
print(f"Overall: {total_clipped:,} / {total_vals:,} values exceed cap ({100*total_clipped/total_vals:.2f}%)")
print()

neur_frac = ct_stats.loc[neur_mask, "frac_clipped"]
nonneur_frac = ct_stats.loc[~neur_mask, "frac_clipped"]
print(f"Mean fraction clipped — Neuronal: {neur_frac.mean():.2%}, Non-neuronal: {nonneur_frac.mean():.2%}")
print(f"Ratio: {nonneur_frac.mean() / neur_frac.mean():.1f}× higher in non-neuronal")
print()

U, p_mw = stats.mannwhitneyu(neur_frac, nonneur_frac, alternative="less")
print(f"Mann-Whitney U (neuronal < non-neuronal): U={U:.0f}, p={p_mw:.2e}")

rho_umi, p_umi = stats.spearmanr(ct_stats["Total_UMI"], ct_stats["frac_clipped"])
print(f"Spearmans' R (Total UMI, frac_clipped): ρ={rho_umi:.3f}, p={p_umi:.2e}")

rho_cells, p_cells = stats.spearmanr(ct_stats["N_cells"], ct_stats["frac_clipped"])
print(f"Spearmans' R (N cells, frac_clipped): ρ={rho_cells:.3f}, p={p_cells:.2e}")

# %% [markdown]
# ## Main Figure (1×2)

# %%
from scipy.stats import gaussian_kde

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# --- Panel A: Total UMI vs Fraction Clipped ---
ax = axes[0]
ax.scatter(ct_stats.loc[neur_mask, "Total_UMI"], ct_stats.loc[neur_mask, "frac_clipped"],
           color="red", alpha=0.5, s=30, edgecolors="white", lw=0.3, label="Neuronal", zorder=3)
ax.scatter(ct_stats.loc[~neur_mask, "Total_UMI"], ct_stats.loc[~neur_mask, "frac_clipped"],
           color="blue", alpha=0.6, s=30, edgecolors="white", lw=0.3, label="Non-neuronal", zorder=4)
ax.set_xscale("log")
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("Fraction of genes exceeding cap")
ax.set_title("A", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=10, framealpha=0.8)
ax.text(0.97, 0.97, f"ρ = {rho_umi:.3f}\np = {p_umi:.1e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# --- Panel B: Specificity Distribution — Example Cell Types ---
ax = axes[1]

# Non-neuronal: Vascular and Fibroblast (lowest UMI, most extreme tails)
vasc_cts = ct_stats.loc[(ct_stats["Supercluster"] == "Vascular") & ~neur_mask]
vasc_ct = vasc_cts.sort_values("Total_UMI").index[0]

fibro_cts = ct_stats.loc[(ct_stats["Supercluster"] == "Fibroblast") & ~neur_mask]
fibro_ct = fibro_cts.sort_values("Total_UMI").index[0]

# Glia: Astrocyte (median UMI)
astro_cts = ct_stats.loc[ct_stats["Supercluster"] == "Astrocyte"]
astro_ct = astro_cts.sort_values("Total_UMI").iloc[len(astro_cts)//2].name

# Cerebellum neuron: Cerebellar inhibitory (median UMI — includes Purkinje)
cereb_cts = ct_stats.loc[(ct_stats["Supercluster"] == "Cerebellar inhibitory") & neur_mask]
cereb_ct = cereb_cts.sort_values("Total_UMI").iloc[len(cereb_cts)//2].name

# CGE interneuron (median UMI) — our key cell type
cge_cts = ct_stats.loc[(ct_stats["Supercluster"] == "CGE interneuron") & neur_mask]
cge_ct = cge_cts.sort_values("Total_UMI").iloc[len(cge_cts)//2].name

# Largest neuronal — cleanest distribution
large_neur_ct = ct_stats.loc[neur_mask].sort_values("Total_UMI").index[-1]

examples = [
    (vasc_ct, "darkblue", f"Vascular (UMI={ct_stats.loc[vasc_ct, 'Total_UMI']:.0f})"),
    (fibro_ct, "royalblue", f"Fibroblast (UMI={ct_stats.loc[fibro_ct, 'Total_UMI']:.0f})"),
    (astro_ct, "#7B68EE", f"Astrocyte (UMI={ct_stats.loc[astro_ct, 'Total_UMI']:.0f})"),
    (cereb_ct, "#FF8C00", f"Cereb. IN (UMI={ct_stats.loc[cereb_ct, 'Total_UMI']:.0f})"),
    (cge_ct, "salmon", f"CGE IN (UMI={ct_stats.loc[cge_ct, 'Total_UMI']:.0f})"),
    (large_neur_ct, "red", f"Neuron (UMI={ct_stats.loc[large_neur_ct, 'Total_UMI']:.0f})"),
]

for ct, color, label in examples:
    vals = spec_unclip[ct].values
    vals_pos = vals[vals > 0.01]  # skip near-zero for cleaner KDE
    log_vals = np.log10(vals_pos)
    kde = gaussian_kde(log_vals, bw_method=0.15)
    x_grid = np.linspace(np.log10(0.01), np.log10(max(vals_pos.max(), 100) * 1.1), 500)
    ax.plot(10**x_grid, kde(x_grid), color=color, lw=2, label=label)
    ax.fill_between(10**x_grid, kde(x_grid), alpha=0.12, color=color)

ax.axvline(x=clip_threshold, color="black", ls="--", lw=1.5, alpha=0.7, label=f"Cap = {clip_threshold:.1f}")
ax.set_xscale("log")
ax.set_xlabel("Specificity (fold-enrichment)")
ax.set_ylabel("Density")
ax.set_title("B", fontweight="bold", loc="left", fontsize=16)
ax.legend(fontsize=8, framealpha=0.8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

fig.tight_layout(w_pad=3)

fig.savefig(FIG_DIR / "specificity_cap_analysis.pdf", dpi=300, transparent=True, bbox_inches="tight")
fig.savefig(FIG_DIR / "specificity_cap_analysis.png", dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR}")

# %% [markdown]
# ## Specificity Distribution: Heavy Tail Dominates Bias
#
# Although only ~5.6% of values exceed the cap, these extreme values carry
# disproportionate weight in the bias calculation (weighted mean specificity).
# Without capping, a single gene with specificity 50× would outweigh 25 genes
# at specificity 2×.

# %%
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))

# --- Left: Overall specificity distribution ---
ax = axes2[0]
all_vals = spec_unclip.values.flatten()
all_pos = all_vals[all_vals > 0.01]

ax.hist(all_pos, bins=np.logspace(np.log10(0.01), np.log10(100), 150),
        color="#888888", alpha=0.7, density=True)
ax.axvline(x=clip_threshold, color="red", ls="--", lw=2, label=f"Cap = {clip_threshold:.1f}")

# Shade the tail
tail_vals = all_pos[all_pos > clip_threshold]
ax.hist(tail_vals, bins=np.logspace(np.log10(clip_threshold), np.log10(100), 50),
        color="red", alpha=0.4, density=True, label=f"Tail: {len(tail_vals)/len(all_pos):.1%} of values")

ax.set_xscale("log")
ax.set_xlabel("Specificity (fold-enrichment)")
ax.set_ylabel("Density")
ax.set_title("Overall Specificity Distribution", fontweight="bold", fontsize=13)
ax.legend(fontsize=10, framealpha=0.8)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# --- Right: Contribution to weighted sum ---
# For each cell type, compute what fraction of the total specificity sum
# comes from genes exceeding the cap (unclipped values)
ax = axes2[1]

frac_sum_from_tail = []
for ct in Anno.index:
    vals = spec_unclip[ct].values
    total_sum = vals.sum()
    tail_sum = vals[vals > clip_threshold].sum()
    frac_sum_from_tail.append(tail_sum / total_sum if total_sum > 0 else 0)

ct_stats["frac_sum_tail"] = frac_sum_from_tail

ax.scatter(ct_stats.loc[neur_mask, "Total_UMI"], ct_stats.loc[neur_mask, "frac_sum_tail"],
           color="red", alpha=0.5, s=30, edgecolors="white", lw=0.3, label="Neuronal", zorder=3)
ax.scatter(ct_stats.loc[~neur_mask, "Total_UMI"], ct_stats.loc[~neur_mask, "frac_sum_tail"],
           color="blue", alpha=0.6, s=30, edgecolors="white", lw=0.3, label="Non-neuronal", zorder=4)
ax.set_xscale("log")
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("Fraction of total specificity sum\nfrom genes exceeding cap")
ax.set_title("Tail Contribution to Bias Calculation", fontweight="bold", fontsize=13)
ax.legend(fontsize=10, framealpha=0.8)

rho_tail, p_tail = stats.spearmanr(ct_stats["Total_UMI"], ct_stats["frac_sum_tail"])
ax.text(0.97, 0.97, f"ρ = {rho_tail:.3f}\np = {p_tail:.1e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# Print summary
neur_tail = ct_stats.loc[neur_mask, "frac_sum_tail"]
nonneur_tail = ct_stats.loc[~neur_mask, "frac_sum_tail"]
print(f"Fraction of total specificity sum from tail (genes > cap):")
print(f"  Neuronal:     {neur_tail.mean():.1%} (median {neur_tail.median():.1%})")
print(f"  Non-neuronal: {nonneur_tail.mean():.1%} (median {nonneur_tail.median():.1%})")
print(f"  Ratio: {nonneur_tail.mean() / neur_tail.mean():.1f}×")

fig2.tight_layout(w_pad=3)
fig2.savefig(FIG_DIR / "specificity_tail_contribution.pdf", dpi=300, transparent=True, bbox_inches="tight")
fig2.savefig(FIG_DIR / "specificity_tail_contribution.png", dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR}")

# %% [markdown]
# ## SCZ Bias Across Cap Values
#
# Using SCZ as an example, we show how the specificity cap affects mutation
# bias rankings:
# - **Cap = 1×**: Over-clips — nearly all specificity variation removed, signal destroyed
# - **Cap = 2×** (default): Preserves biological signal while controlling inflation
# - **Cap = 3×**: Similar to 2×, confirming robustness
# - **No cap**: Non-neuronal cell types with inflated specificity dominate the top ranks

# %%
from CellType_PSY import HumanCT_AvgZ_Weighted, AnnotateCTDat

# Load SCZ gene weights
gw_df = pd.read_csv(
    PROJ_DIR / "dat" / "GeneWeights" / "SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
    header=None, names=["Entrez", "Weight"]
)
scz_weights = dict(zip(gw_df["Entrez"].values, gw_df["Weight"].values))
print(f"SCZ: {len(scz_weights)} genes")

# Prepare cap levels: compute clipped matrices on the fly from unclipped
mean_spec = np.mean(spec_unclip.values.flatten())
cap_configs = {
    "Cap = 1×": 1.0,
    "Cap = 2× (default)": 2.0,
    "Cap = 3×": 3.0,
    "No cap": 100.0,
}

bias_results = {}
for label, multiplier in cap_configs.items():
    threshold = mean_spec * multiplier
    spec_clipped = spec_unclip.clip(lower=0, upper=threshold)
    # Mean-center per cell type (same as preprocessing pipeline)
    spec_centered = spec_clipped.subtract(spec_clipped.mean(axis=0), axis=1)
    # Compute bias
    bias = HumanCT_AvgZ_Weighted(spec_centered, scz_weights)
    bias = AnnotateCTDat(bias, Anno)
    bias_results[label] = bias
    # Print top 5
    top5 = bias.head(5)
    print(f"\n{label} (threshold={threshold:.2f}) — Top 5:")
    for rank, (ct, row) in enumerate(top5.iterrows(), 1):
        sc = row.get("Supercluster", "")
        cls = "NEUR" if int(ct) in Neur_idx else "non-NEUR"
        print(f"  {rank}. CT {ct}: EFFECT={row['EFFECT']:.4f}  [{sc}] ({cls})")

# %%
fig3, axes3 = plt.subplots(1, 4, figsize=(22, 6), sharey=False)

cap_colors = {
    "Cap = 1×": "#999999",
    "Cap = 2× (default)": "#2563eb",
    "Cap = 3×": "#16a34a",
    "No cap": "#dc2626",
}

for idx, (label, bias) in enumerate(bias_results.items()):
    ax = axes3[idx]
    color = cap_colors[label]

    # Get top 20 cell types
    top20 = bias.head(20).copy()
    top20["ct_id"] = top20.index.astype(int)
    top20["is_neur"] = top20["ct_id"].isin(Neur_idx)

    # Determine short labels
    names = []
    for ct in top20.index:
        ct_int = int(ct)
        sc = Anno.loc[ct_int, "Supercluster"] if ct_int in Anno.index else ""
        # Abbreviate supercluster
        short = sc[:20] + "..." if len(sc) > 20 else sc
        names.append(f"{ct_int} ({short})")
    top20["name"] = names

    # Horizontal bar plot
    y_pos = np.arange(len(top20))
    bar_colors = ["red" if n else "blue" for n in top20["is_neur"]]
    ax.barh(y_pos, top20["EFFECT"].values, color=bar_colors, alpha=0.7, edgecolor="white", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top20["name"].values, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Mutation Bias (EFFECT)")
    ax.set_title(label, fontweight="bold", fontsize=13, color=color)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Count neuronal in top 20
    n_neur_top20 = top20["is_neur"].sum()
    ax.text(0.97, 0.97, f"Neuronal: {n_neur_top20}/20",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# Add shared legend
from matplotlib.patches import Patch
axes3[0].legend(handles=[Patch(facecolor="red", alpha=0.7, label="Neuronal"),
                         Patch(facecolor="blue", alpha=0.7, label="Non-neuronal")],
                fontsize=9, loc="lower right")

fig3.suptitle("SCZ Mutation Bias — Top 20 Cell Types at Different Specificity Caps",
              fontsize=15, fontweight="bold", y=1.02)
fig3.tight_layout(w_pad=2)
fig3.savefig(FIG_DIR / "scz_bias_cap_comparison.pdf", dpi=300, transparent=True, bbox_inches="tight")
fig3.savefig(FIG_DIR / "scz_bias_cap_comparison.png", dpi=300, transparent=True, bbox_inches="tight")
plt.show()
print(f"Saved to {FIG_DIR}")

# %% [markdown]
# ## Summary Table

# %%
summary = pd.DataFrame({
    "Neuronal (n=378)": [
        f"{ct_stats.loc[neur_mask, 'Total_UMI'].mean():.0f}",
        f"{ct_stats.loc[neur_mask, 'N_cells'].median():.0f}",
        f"{ct_stats.loc[neur_mask, 'frac_clipped'].mean():.2%}",
        f"{ct_stats.loc[neur_mask, 'max_spec'].median():.1f}",
        f"{ct_stats.loc[neur_mask, 'max_spec'].min():.1f} – {ct_stats.loc[neur_mask, 'max_spec'].max():.1f}",
    ],
    "Non-neuronal (n=83)": [
        f"{ct_stats.loc[~neur_mask, 'Total_UMI'].mean():.0f}",
        f"{ct_stats.loc[~neur_mask, 'N_cells'].median():.0f}",
        f"{ct_stats.loc[~neur_mask, 'frac_clipped'].mean():.2%}",
        f"{ct_stats.loc[~neur_mask, 'max_spec'].median():.1f}",
        f"{ct_stats.loc[~neur_mask, 'max_spec'].min():.1f} – {ct_stats.loc[~neur_mask, 'max_spec'].max():.1f}",
    ],
}, index=["Mean Total UMI", "Median cell count", "Mean fraction clipped", "Median max specificity", "Max specificity range"])

print(summary.to_string())
