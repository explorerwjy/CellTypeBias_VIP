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
# # Dopamine Receptor Expression Across Neuronal Superclusters
#
# Reviewer 2 Q5 asked whether VIP interneurons express dopamine receptors.
# This notebook shows DRD1–DRD5 expression specificity across all neuronal
# superclusters in the Siletti et al. human brain atlas.

# %% [markdown]
# ## Setup

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

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

from CellType_PSY import Anno, Neur_idx

FIG_DIR = PROJ_DIR / "results" / "figures" / "drd_expression"
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
# ## Configuration
#
# Set `MODE` to switch between expression measures:
# - `"mean_centered"` — mean-centered specificity (default for bias analysis)
# - `"tpm"` — raw TPM expression levels

# %%
# ---- Toggle this to switch modes ----
MODE = "mean_centered"  # "mean_centered" or "tpm"
# --------------------------------------

DATA_FILES = {
    "mean_centered": "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    "tpm":           "HumanCT.TPM.0.1.Filt.csv",
}
AXIS_LABELS = {
    "mean_centered": "Mean-Centered Specificity",
    "tpm":           "TPM Expression",
}
# Reference line: 0 for mean-centered, None for TPM
REF_LINE = {"mean_centered": 0.0, "tpm": None}
# Heatmap colormap: diverging for mean-centered, sequential for TPM
HEATMAP_CMAP = {"mean_centered": "RdBu_r", "tpm": "YlOrRd"}

assert MODE in DATA_FILES, f"Unknown MODE '{MODE}', use 'mean_centered' or 'tpm'"
XLABEL = AXIS_LABELS[MODE]
print(f"MODE = {MODE}")

# %% [markdown]
# ## Load Expression Data

# %%
spec = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / DATA_FILES[MODE],
    index_col=0
)
spec.columns = [int(c) for c in spec.columns]
print(f"Expression matrix ({MODE}): {spec.shape[0]} genes × {spec.shape[1]} cell types")
print(f"Value range: [{spec.values.min():.3f}, {spec.values.max():.3f}], mean={spec.values.mean():.4f}")

# DRD genes: Entrez IDs
DRD_GENES = {
    "DRD1": 1812,
    "DRD2": 1813,
    "DRD3": 1814,
    "DRD4": 1815,
    "DRD5": 1816,
}

# Verify all present
for name, eid in DRD_GENES.items():
    assert eid in spec.index, f"{name} (Entrez {eid}) not in specificity matrix"
print("All DRD1-DRD5 found in specificity matrix")

# %% [markdown]
# ## Compute Per-Supercluster DRD Specificity

# %%
# Build cluster → supercluster mapping for neuronal clusters only
neur_anno = Anno.loc[[c for c in Neur_idx if c in Anno.index]]
cluster_to_sc = neur_anno["Supercluster"].to_dict()

# Extract DRD specificity for neuronal clusters
drd_entrez = list(DRD_GENES.values())
drd_names = list(DRD_GENES.keys())
neur_cols = [c for c in spec.columns if c in cluster_to_sc]

drd_spec = spec.loc[drd_entrez, neur_cols].copy()
drd_spec.index = drd_names

# Map columns to superclusters
col_sc = pd.Series({c: cluster_to_sc[c] for c in neur_cols})

# For each supercluster: average specificity per DRD gene across clusters,
# then average across the 5 DRD genes → mean ± SEM
records = []
for sc in sorted(col_sc.unique()):
    sc_cols = col_sc[col_sc == sc].index.tolist()
    # Per-gene mean across clusters in this supercluster
    gene_means = drd_spec[sc_cols].mean(axis=1)  # 5 values (one per DRD)
    records.append({
        "Supercluster": sc,
        "mean_spec": gene_means.mean(),
        "sem_spec": gene_means.sem(),
        "n_clusters": len(sc_cols),
    })

sc_df = pd.DataFrame(records).sort_values("mean_spec", ascending=True).reset_index(drop=True)
print(sc_df.to_string(index=False))

# %% [markdown]
# ## Bar Plot: Mean DRD1–DRD5 Specificity by Supercluster

# %%
# Highlight CGE interneuron
colors = ["#d62728" if sc == "CGE interneuron" else "#1f77b4" for sc in sc_df["Supercluster"]]

fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(
    sc_df["Supercluster"], sc_df["mean_spec"],
    xerr=sc_df["sem_spec"], capsize=3,
    color=colors, edgecolor="none", alpha=0.85
)
ax.set_xlabel(f"{XLABEL} (DRD1–DRD5 ± SEM)")
ax.set_title("Dopamine Receptor Expression Across Neuronal Superclusters")
if REF_LINE[MODE] is not None:
    ax.axvline(REF_LINE[MODE], ls="--", color="grey", alpha=0.5,
               label=f"Genome-wide mean ({REF_LINE[MODE]})")
    ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG_DIR / "DRD_mean_specificity_by_supercluster.pdf", transparent=True, dpi=300)
fig.savefig(FIG_DIR / "DRD_mean_specificity_by_supercluster.png", transparent=True, dpi=300)
plt.show()
print(f"Saved to {FIG_DIR}")

# %% [markdown]
# ## Heatmap: Individual DRD Gene Specificity by Supercluster

# %%
# Build supercluster × DRD gene matrix
heat_data = pd.DataFrame(index=sc_df["Supercluster"], columns=drd_names, dtype=float)
for sc in heat_data.index:
    sc_cols = col_sc[col_sc == sc].index.tolist()
    heat_data.loc[sc] = drd_spec[sc_cols].mean(axis=1).values

fig, ax = plt.subplots(figsize=(6, 8))
heat_vals = heat_data.values.astype(float)
cmap = HEATMAP_CMAP[MODE]
if MODE == "mean_centered":
    vabs = max(abs(heat_vals.min()), abs(heat_vals.max()))
    im = ax.imshow(heat_vals, aspect="auto", cmap=cmap, vmin=-vabs, vmax=vabs)
else:
    im = ax.imshow(heat_vals, aspect="auto", cmap=cmap)
    vabs = heat_vals.max()
ax.set_xticks(range(len(drd_names)))
ax.set_xticklabels(drd_names, rotation=45, ha="right")
ax.set_yticks(range(len(heat_data)))
ax.set_yticklabels(heat_data.index)
plt.colorbar(im, ax=ax, label=XLABEL, shrink=0.6)
ax.set_title(f"DRD1–DRD5 {XLABEL} by Supercluster")

# Annotate values
fmt = ".2f" if MODE == "mean_centered" else ".1f"
for i in range(heat_data.shape[0]):
    for j in range(heat_data.shape[1]):
        val = heat_data.iloc[i, j]
        if MODE == "mean_centered":
            color = "white" if abs(val) > vabs * 0.65 else "black"
        else:
            color = "white" if val > vabs * 0.7 else "black"
        ax.text(j, i, f"{val:{fmt}}", ha="center", va="center", fontsize=7, color=color)

fig.tight_layout()
fig.savefig(FIG_DIR / "DRD_heatmap_by_supercluster.pdf", transparent=True, dpi=300)
fig.savefig(FIG_DIR / "DRD_heatmap_by_supercluster.png", transparent=True, dpi=300)
plt.show()

# %% [markdown]
# ## Detailed View: Per-Cluster Specificity for Each DRD Gene
#
# For each DRD gene, show specificity of every individual cell type cluster,
# grouped by supercluster. This gives a more detailed picture than the
# supercluster-level summary above.

# %%
def plot_gene_spec_by_supercluster(gene_spec, gene_name, anno, neur_idx,
                                   highlight_sc="CGE interneuron",
                                   figsize=(10, 7), ax=None,
                                   sc_order=None, ref_line=None,
                                   xlabel="Expression"):
    """
    Horizontal strip plot of per-cluster expression specificity for one gene,
    grouped and sorted by supercluster median.

    Each dot = one cell type cluster. Superclusters sorted by median specificity.
    Box shows IQR + median line; individual clusters shown as jittered points.

    If sc_order is provided, use that ordering; otherwise sort by median.
    ref_line: x-value for vertical reference line, or None to skip.
    """
    # Build long-form DataFrame
    neur_clusters = [c for c in neur_idx if c in anno.index and c in gene_spec.index]
    rows = []
    for c in neur_clusters:
        rows.append({
            "cluster": c,
            "Supercluster": anno.loc[c, "Supercluster"],
            "spec": gene_spec[c],
        })
    df = pd.DataFrame(rows)

    # Sort superclusters by median spec, or use provided order
    if sc_order is None:
        sc_order = (df.groupby("Supercluster")["spec"]
                    .median()
                    .sort_values(ascending=True)
                    .index.tolist())

    n_sc = len(sc_order)

    show_plot = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Box plot — thin grey boxes
    bp = ax.boxplot(
        [df.loc[df["Supercluster"] == sc, "spec"].values for sc in sc_order],
        positions=range(n_sc), vert=False, patch_artist=False, widths=0.55,
        boxprops=dict(color="#999999", linewidth=0.8),
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="#aaaaaa", linewidth=0.7),
        capprops=dict(color="#aaaaaa", linewidth=0.7),
        flierprops=dict(marker="none"),
    )

    # Overlay individual clusters as jittered dots
    # Grey for everything, red for highlighted supercluster
    for i, sc in enumerate(sc_order):
        vals = df.loc[df["Supercluster"] == sc, "spec"].values
        jitter = np.random.default_rng(42).normal(0, 0.12, size=len(vals))
        y = i + jitter
        is_hl = (sc == highlight_sc)
        ax.scatter(
            vals, y,
            s=35 if is_hl else 16,
            color="#d62728" if is_hl else "#7f7f7f",
            alpha=0.9 if is_hl else 0.45,
            edgecolors="white", linewidth=0.3,
            zorder=5 if is_hl else 3,
        )

    ax.set_yticks(range(n_sc))
    ax.set_yticklabels(sc_order, fontsize=10)
    if ref_line is not None:
        ax.axvline(ref_line, ls="--", color="grey", alpha=0.4, lw=1)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_title(gene_name, fontsize=14, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)

    if show_plot:
        fig.tight_layout()
    return ax


def style_supercluster_labels(ax, highlight_sc="CGE interneuron"):
    """Bold + red the highlighted supercluster y-label. Call after fig.tight_layout()."""
    ax.figure.canvas.draw()
    for label in ax.get_yticklabels():
        if label.get_text() == highlight_sc:
            label.set_fontweight("bold")
            label.set_color("#d62728")
        else:
            label.set_fontweight("normal")
            label.set_color("black")


# %%
# Plot all 5 DRD genes as a 1×5 layout with shared y-axis order
# Use raw TPM expression for the x-axis in this figure.
tpm = pd.read_csv(
    PROJ_DIR / "dat" / "ExpMats" / DATA_FILES["tpm"],
    index_col=0,
)
tpm.columns = [int(c) for c in tpm.columns]
TPM_XLABEL = AXIS_LABELS["tpm"]

# Compute shared sort order: mean of per-gene median TPM across all 5 DRDs
neur_clusters = [c for c in Neur_idx if c in Anno.index and c in tpm.columns]
_tmp_df = pd.DataFrame({
    "cluster": neur_clusters,
    "Supercluster": [Anno.loc[c, "Supercluster"] for c in neur_clusters],
})
# For each DRD gene, get median TPM per supercluster; then average across genes
_medians = {}
for gname, eid in DRD_GENES.items():
    _s = tpm.loc[eid, neur_clusters]
    _tmp_df[gname] = _s.values
    _medians[gname] = _tmp_df.groupby("Supercluster")[gname].median()
_median_df = pd.DataFrame(_medians)
shared_sc_order = _median_df.mean(axis=1).sort_values(ascending=True).index.tolist()

# Increase base font size by 50%
LABEL_FONTSIZE = 18   # xlabels
TITLE_FONTSIZE = 21   # main plot title
SUPTITLE_FONTSIZE = 24
TICK_FONTSIZE = 15
YLABEL_FONTSIZE = 15

fig, axes = plt.subplots(1, 5, figsize=(36, 8), sharey=False)

for idx, (gene_name, entrez_id) in enumerate(DRD_GENES.items()):
    ax = axes[idx]
    gene_tpm = tpm.loc[entrez_id]
    plot_gene_spec_by_supercluster(
        gene_tpm, gene_name, Anno, Neur_idx,
        highlight_sc="CGE interneuron", ax=ax,
        sc_order=shared_sc_order,
        ref_line=REF_LINE["tpm"], xlabel=TPM_XLABEL,
    )
    # Set x-axis to show the TPM range with padding
    neur_vals = gene_tpm[[c for c in tpm.columns if c in cluster_to_sc]].values
    xmin_val, xmax_val = np.percentile(neur_vals, [1, 99])
    pad = (xmax_val - xmin_val) * 0.15
    ax.set_xlim(max(0, xmin_val - pad), xmax_val + pad)
    # Only show y-labels on the leftmost panel
    if idx > 0:
        ax.set_yticklabels([])
        ax.set_ylabel("")
    # Adjust x/y label font size and title
    ax.set_xlabel(TPM_XLABEL, fontsize=LABEL_FONTSIZE)
    ax.set_title(gene_name, fontsize=TITLE_FONTSIZE, fontweight="bold")
    ax.tick_params(axis="x", labelsize=TICK_FONTSIZE)
    ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    for label in ax.get_xticklabels():
        label.set_fontsize(TICK_FONTSIZE)
    if idx == 0:
        try:
            ax.set_ylabel(ax.get_ylabel(), fontsize=YLABEL_FONTSIZE)
            for label in ax.get_yticklabels():
                label.set_fontsize(TICK_FONTSIZE)
        except Exception:
            pass

# fig.suptitle(
#     "Dopamine Receptor TPM Expression — Individual Cell Type Clusters",
#     fontsize=SUPTITLE_FONTSIZE, fontweight="bold", y=1.01
# )
# fig.tight_layout()

# Style labels on leftmost panel (only one with visible labels)
style_supercluster_labels(axes[0], "CGE interneuron")

fig.savefig(FIG_DIR / "DRD_per_cluster_all5.pdf", transparent=True, dpi=300, bbox_inches="tight")
fig.savefig(FIG_DIR / "DRD_per_cluster_all5.png", transparent=True, dpi=300, bbox_inches="tight")
plt.show()

# %%
# Also plot each DRD gene individually (larger, more readable)
for gene_name, entrez_id in DRD_GENES.items():
    gene_spec = spec.loc[entrez_id]
    fig, ax = plt.subplots(figsize=(10, 7))
    plot_gene_spec_by_supercluster(
        gene_spec, gene_name, Anno, Neur_idx,
        highlight_sc="CGE interneuron", ax=ax,
        ref_line=REF_LINE[MODE], xlabel=XLABEL,
    )
    fig.tight_layout()
    style_supercluster_labels(ax, "CGE interneuron")
    fig.savefig(FIG_DIR / f"{gene_name}_per_cluster.pdf", transparent=True, dpi=300)
    fig.savefig(FIG_DIR / f"{gene_name}_per_cluster.png", transparent=True, dpi=300)
    plt.show()

# %% [markdown]
# ## Summary Statistics

# %%
cge_row = sc_df[sc_df["Supercluster"] == "CGE interneuron"].iloc[0]
print(f"MODE: {MODE}")
print(f"CGE interneuron DRD mean {XLABEL}: {cge_row['mean_spec']:.3f} ± {cge_row['sem_spec']:.3f}")
print(f"CGE interneuron rank: {sc_df[sc_df['Supercluster'] == 'CGE interneuron'].index[0] + 1} / {len(sc_df)}")
print()

# Per-gene details for CGE
cge_cols = col_sc[col_sc == "CGE interneuron"].index.tolist()
for gene in drd_names:
    val = drd_spec.loc[gene, cge_cols].mean()
    print(f"  {gene}: {val:.3f}")
