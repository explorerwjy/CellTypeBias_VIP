# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (gencic)
#     language: python
#     name: gencic
# ---

# %% [markdown]
# # Supplementary Figures and Tables
#
# Figures S1–S11 ordered by first reference in the revised manuscript.
# S1, S3 are pre-made external figures (not generated in this notebook).
#
# | Fig | Content | Source |
# |-----|---------|--------|
# | **S1** | Full-scale IQ distribution (SPARK & ASC) | External |
# | **S2** | Temporal expression patterns across brain development | Inline |
# | **S3** | Workflow for computing cell type mutation biases | External |
# | **S4** | Specificity capping validation & noise from low expression | `Specificity_Cap_Analysis`, `ZINB_Simulation`, `Cap_Sensitivity` |
# | **S5** | Impact of gene set size + real vs random expansion | `Number_Gene_Effect` |
# | **S6** | Impact of genetic archeture / downsample mutations | `Number_Gene_Effect` |
# | **S7** | Comprehensive mutation biases across brain cell types | Inline |
# | **S8** | Negative controls (non-brain traits) & SCZ protective genes | `NegativeControl_BiasPlot`, `SCZ_Protective_BiasPlot` |
# | **S9** | Impact of gene expression levels on ASD-SCZ bias correlation | `Similarity_ASD_SCZ.spec` |
# | **S10** | ASD and SCZ mutation bias towards different cell type superclusters | Inline |
# | **S11** | Mutation bias comparison: MGE & LAMP5-LHX6/Chandelier | Inline |
# | **S12** | Comprehensive analysis of mutation biases in 22q11.2 deletion | Inline |
#
# **Prerequisites:**
# - Run `Number_Gene_Effect.ipynb` for gene sweep figures (S5)
# - Run `Similarity_ASD_SCZ.spec.ipynb` for BrainSpan gene removal figure (S8)
# - Run rebuttal notebooks for: specificity cap (S4), negative controls & SCZ protective (S7)
#
# Tables S2–S11 follow the figures at the end.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os
import io
import shutil
import subprocess
from contextlib import contextmanager

from pathlib import Path
import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *
from matplotlib.image import imread
from PIL import Image as PILImage

import matplotlib.font_manager as fm
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)
fm._load_fontmanager(try_read_cache=False)
plt.style.use('seaborn-v0_8-whitegrid')

FIG_DIR = str(PROJ_DIR / "results/figures/supp/") + "/"
os.makedirs(FIG_DIR, exist_ok=True)

REBUTTAL_FIG_DIR = str(PROJ_DIR / "results/figures/") + "/"


# %%
# -- Utility functions --

@contextmanager
def save_panel(filepath, dpi=300):
    """Intercept plt.show() to save the figure before displaying inline."""
    _orig = plt.show
    plt.show = lambda *a, **kw: None
    try:
        yield
    finally:
        plt.show = _orig
        fig = plt.gcf()
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', transparent=True,
                    facecolor='none')
        plt.close(fig)
        from IPython.display import Image as IPImage, display as ipdisplay
        ipdisplay(IPImage(filename=filepath))


def assemble_panels(panel_paths, panel_labels, layout, figsize, out_path, dpi=300):
    """Assemble saved panel PNGs into a multi-panel composite figure."""
    nrows = max(r for r, _ in layout) + 1
    ncols = max(c.stop for _, c in layout)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='white')
    gs = fig.add_gridspec(nrows, ncols, hspace=0.08, wspace=0.08)

    for path, label, (row, col) in zip(panel_paths, panel_labels, layout):
        ax = fig.add_subplot(gs[row, col])
        img = imread(path)
        ax.imshow(img)
        ax.axis('off')
        ax.text(-0.02, 1.05, label, transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='bottom', ha='right')

    fig.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.show()
    print(f"Saved composite: {out_path}")


def pdf_to_png(pdf_path, dpi=300):
    """Convert PDF to PNG using pdftoppm. Returns PNG path."""
    png_path = pdf_path.replace('.pdf', '.png')
    if os.path.exists(png_path):
        return png_path
    prefix = png_path.replace('.png', '')
    subprocess.run(['pdftoppm', '-png', '-r', str(dpi), '-singlefile',
                    pdf_path, prefix], check=True)
    print(f"Converted: {pdf_path} → {png_path}")
    return png_path


def copy_figure(src_path, fig_name, convert_pdf=True):
    """Copy a pre-generated figure to FIG_DIR with proper naming."""
    ext = os.path.splitext(src_path)[1]
    if ext == '.png':
        shutil.copy2(src_path, FIG_DIR + fig_name + ".png")
        img = PILImage.open(FIG_DIR + fig_name + ".png")
        img.save(FIG_DIR + fig_name + ".pdf", "PDF", resolution=300)
    elif ext == '.pdf':
        shutil.copy2(src_path, FIG_DIR + fig_name + ".pdf")
        if convert_pdf:
            pdf_to_png(FIG_DIR + fig_name + ".pdf")
    print(f"Saved: {FIG_DIR}{fig_name}.pdf")
    from IPython.display import Image as IPImage, display as ipdisplay
    png = FIG_DIR + fig_name + ".png"
    if os.path.exists(png):
        ipdisplay(IPImage(filename=png))


# %% [markdown]
# ## Load data

# %%
# Bias DataFrames
Bias_Save_Dir = str(PROJ_DIR / "results/main_results/random/Centering/") + "/"

ASD_All_Bias = pd.read_csv(Bias_Save_Dir + "ASD_All_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)
VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)

# %%
# Pre-computed contrasts
CONTRAST_DIR = str(PROJ_DIR / "results/main_results/contrasts/") + "/"

ASD_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_SCZ_contrast.csv", index_col=0)
ASD_wID_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_wID_vs_SCZ_contrast.csv", index_col=0)
HIQ_LIQ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_ASD_wID_contrast.csv", index_col=0)
ASD_DDD_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_DDD_contrast.csv", index_col=0)
SCZ_ASD_wID_Contrast = pd.read_csv(CONTRAST_DIR + "SCZ_vs_ASD_wID_contrast.csv", index_col=0)
VNR_Contrast = pd.read_csv(CONTRAST_DIR + "VNR_neg_vs_pos_contrast.csv", index_col=0)
EDU_Contrast = pd.read_csv(CONTRAST_DIR + "EDU_neg_vs_pos_contrast.csv", index_col=0)
all_contrasts_df = pd.read_csv(CONTRAST_DIR + "all_contrasts_fdr.csv")
brainspan_df = pd.read_csv(CONTRAST_DIR + "brainspan_expression.csv")

# Neuron-filtered contrasts
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
HIQ_LIQ_Contrast_Neurons = HIQ_LIQ_Contrast[HIQ_LIQ_Contrast.index.isin(Neurons)]
ASD_DDD_Contrast_Neurons = ASD_DDD_Contrast[ASD_DDD_Contrast.index.isin(Neurons)]
VNR_Contrast_Neurons = VNR_Contrast[VNR_Contrast.index.isin(Neurons)]
EDU_Contrast_Neurons = EDU_Contrast[EDU_Contrast.index.isin(Neurons)]

EffLabel = "EFFECT"

# %% [markdown]
# ---
# ## Figure S2 — Temporal Expression Patterns Across Brain Development

# %%
palette = sns.color_palette("Set2", 5)

Time = ['mean_2A', 'mean_2B', 'mean_3A', 'mean_3B', 'mean_4', 'mean_5',
        'mean_6', 'mean_7', 'mean_8', 'mean_9', 'mean_10', 'mean_11']

labels_time = [
    "Embryonic", "Early fetal", "Early mid-fetal", "Late mid-fetal", "Late fetal",
    "Early infancy", "Late infancy", "Early Childhood", "Late childhood", "Adolescence",
    "Young adulthood", "Adulthood"
]

gene_set_styles = {
    "ASD w/o ID": dict(color=palette[0], marker='o', label="ASD genes"),
    "ASD with ID": dict(color=palette[1], marker='^', label="ASD with ID genes"),
    "SCZ": dict(color=palette[2], marker='s', label="SCZ genes"),
    "ALL": dict(color=palette[4], marker='*', label="ALL Genes", ls="--", markersize=12, alpha=0.8, zorder=2),
}

fig, ax = plt.subplots(dpi=300, figsize=(12, 8), constrained_layout=True)

xvals = np.arange(len(Time))
for gene_set, style in gene_set_styles.items():
    subset = brainspan_df[brainspan_df["gene_set"] == gene_set].sort_values("stage", key=lambda s: s.map({t: i for i, t in enumerate(Time)}))
    shift = {"ASD w/o ID": 0.05, "SCZ": -0.05}.get(gene_set, 0.0)
    ax.errorbar(
        x=xvals + shift,
        y=subset["mean"].values,
        yerr=subset["sem"].values,
        linewidth=2, markersize=style.get("markersize", 8),
        alpha=style.get("alpha", 0.95), zorder=style.get("zorder", 3),
        capsize=5, elinewidth=2, capthick=2,
        color=style["color"], marker=style["marker"], label=style["label"],
        ls=style.get("ls", "-"),
    )

ax.axvline(x=5, color='grey', ls="--", linewidth=1.5, zorder=1)

ax.set_xticks(xvals)
ax.set_xticklabels(labels_time, rotation=30, ha='right', fontsize=13, weight='bold')
ax.set_ylabel("Brain expression level (log10(RPKM))", fontsize=18, weight='bold')
ax.set_xlabel("Developmental stage", fontsize=16, weight='bold', labelpad=10)

sns.despine(ax=ax)
ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.5, zorder=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.7, 1), fancybox=True, shadow=True,
          ncol=2, fontsize=14, frameon=True, borderaxespad=0.5, title="Gene sets", title_fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=13, width=1.5)

fig.savefig(FIG_DIR + "FigureS2.pdf", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()
print(f"Saved: {FIG_DIR}FigureS2.pdf")

# %% [markdown]
# ---
# ## Figure S4 — Specificity Capping Validation & Noise from Low Expression
#
# | Panel | Content | Source |
# |-------|---------|--------|
# | A | Empirical specificity inflation by UMI depth | `Specificity_Cap_Analysis.ipynb` |
# | B | NB simulation: sampling noise inflates specificity | `Specificity_ZINB_Simulation.ipynb` |
# | C | Cap sensitivity sweep (1x–10x) | `Cap_Sensitivity_Figure.ipynb` |

# %%
# --- S4 Data Loading ---
# Load unclipped specificity matrix (clip100 is effectively unclipped)
from scipy.stats import gaussian_kde
spec_unclip = pd.read_csv(
    str(PROJ_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv"), index_col=0)
spec_unclip.columns = [int(c) for c in spec_unclip.columns]
clip_threshold = np.mean(spec_unclip.values.flatten()) * 2
mean_spec = np.mean(spec_unclip.values.flatten())
print(f"Clip threshold (2x mean): {clip_threshold:.4f}")

# Compute per-cell-type clipping statistics
ct_stats = pd.DataFrame(index=Anno.index)
for ct in ct_stats.index:
    vals = spec_unclip[ct].values
    ct_stats.loc[ct, "frac_clipped"] = np.mean(vals > clip_threshold)
    ct_stats.loc[ct, "max_spec"] = np.max(vals)
ct_stats["Total_UMI"] = Anno["Total UMI"]
ct_stats["Supercluster"] = Anno["Supercluster"]
ct_stats["is_neuronal"] = ct_stats.index.isin(Neur_idx)
for col in ["frac_clipped", "max_spec", "Total_UMI"]:
    ct_stats[col] = pd.to_numeric(ct_stats[col])
neur_mask_cap = ct_stats["is_neuronal"]

from scipy import stats as sp_stats
rho_umi, p_umi = sp_stats.spearmanr(ct_stats["Total_UMI"], ct_stats["frac_clipped"])

# %% [markdown]
# ### Panel A — Total UMI vs Fraction of Genes Exceeding Cap

# %%
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(ct_stats.loc[neur_mask_cap, "Total_UMI"], ct_stats.loc[neur_mask_cap, "frac_clipped"],
           color="red", alpha=0.5, s=30, edgecolors="white", lw=0.3, label="Neuronal", zorder=3)
ax.scatter(ct_stats.loc[~neur_mask_cap, "Total_UMI"], ct_stats.loc[~neur_mask_cap, "frac_clipped"],
           color="blue", alpha=0.6, s=30, edgecolors="white", lw=0.3, label="Non-neuronal", zorder=4)
ax.set_xscale("log")
ax.set_xlabel("Total UMI per cell type", fontsize=12)
ax.set_ylabel("Fraction of genes exceeding cap", fontsize=12)
ax.legend(fontsize=10, framealpha=0.8)
ax.text(0.97, 0.97, f"ρ = {rho_umi:.3f}\np = {p_umi:.1e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_A.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel B — Specificity Distribution: Example Cell Types

# %%
# Pick example cell types spanning UMI range
vasc_ct = ct_stats.loc[(ct_stats["Supercluster"] == "Vascular") & ~neur_mask_cap].sort_values("Total_UMI").index[0]
astro_cts = ct_stats.loc[ct_stats["Supercluster"] == "Astrocyte"]
astro_ct = astro_cts.sort_values("Total_UMI").iloc[len(astro_cts)//2].name
cge_cts = ct_stats.loc[(ct_stats["Supercluster"] == "CGE interneuron") & neur_mask_cap]
cge_ct = cge_cts.sort_values("Total_UMI").iloc[len(cge_cts)//2].name
large_neur_ct = ct_stats.loc[neur_mask_cap].sort_values("Total_UMI").index[-1]

examples = [
    (vasc_ct, "darkblue", f"Vascular (UMI={ct_stats.loc[vasc_ct, 'Total_UMI']:.0f})"),
    (astro_ct, "#7B68EE", f"Astrocyte (UMI={ct_stats.loc[astro_ct, 'Total_UMI']:.0f})"),
    (cge_ct, "salmon", f"CGE IN (UMI={ct_stats.loc[cge_ct, 'Total_UMI']:.0f})"),
    (large_neur_ct, "red", f"Neuron (UMI={ct_stats.loc[large_neur_ct, 'Total_UMI']:.0f})"),
]

fig, ax = plt.subplots(figsize=(7, 6))
for ct, color, label in examples:
    vals = spec_unclip[ct].values
    vals_pos = vals[vals > 0.01]
    log_vals = np.log10(vals_pos)
    kde = gaussian_kde(log_vals, bw_method=0.15)
    x_grid = np.linspace(np.log10(0.01), np.log10(max(vals_pos.max(), 100) * 1.1), 500)
    ax.plot(10**x_grid, kde(x_grid), color=color, lw=2, label=label)
    ax.fill_between(10**x_grid, kde(x_grid), alpha=0.12, color=color)
ax.axvline(x=clip_threshold, color="black", ls="--", lw=1.5, alpha=0.7, label=f"Cap = {clip_threshold:.1f}")
ax.set_xscale("log")
ax.set_xlabel("Specificity (fold-enrichment)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.legend(fontsize=8, framealpha=0.8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_B.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel C — ZINB Simulation: Max Specificity by Expression × UMI
#
# Pre-computed by `Specificity_ZINB_Simulation.ipynb`; loaded from saved data.

# %%
import pickle
PLOT_DATA_DIR = str(PROJ_DIR / "results/figures/plot_data/") + "/"
with open(PLOT_DATA_DIR + "zinb_data.pkl", "rb") as f:
    zinb_data = pickle.load(f)

zinb_sweep = zinb_data["sweep_results"]
zinb_pcts = zinb_data["expression_percentiles"]
zinb_ct_ids = zinb_data["ct_ids"]
zinb_clip = zinb_data["clip_threshold"]
total_umis_zinb = np.array([Anno.loc[ct, "Total UMI"] for ct in zinb_ct_ids])

# Heatmap: max specificity by expression level x UMI bin
n_umi_bins = 8
umi_bin_edges = np.logspace(np.log10(total_umis_zinb.min() * 0.9),
                            np.log10(total_umis_zinb.max() * 1.1), n_umi_bins + 1)
umi_bin_labels = [f"{umi_bin_edges[i]:.0f}-{umi_bin_edges[i+1]:.0f}" for i in range(n_umi_bins)]

heatmap_data = np.zeros((len(zinb_pcts), n_umi_bins))
for i, pct in enumerate(zinb_pcts):
    res = zinb_sweep[pct]
    for j in range(n_umi_bins):
        mask = (res["total_umi"] >= umi_bin_edges[j]) & (res["total_umi"] < umi_bin_edges[j+1])
        heatmap_data[i, j] = res.loc[mask, "max_spec"].median() if mask.sum() > 0 else np.nan

from matplotlib.colors import LogNorm
fig, ax = plt.subplots(figsize=(8, 5))
im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd",
    norm=LogNorm(vmin=max(1, np.nanmin(heatmap_data[heatmap_data > 0])),
                 vmax=np.nanmax(heatmap_data)), origin="lower")
for i in range(len(zinb_pcts)):
    for j in range(n_umi_bins):
        val = heatmap_data[i, j]
        if not np.isnan(val):
            tc = "white" if val > np.nanmax(heatmap_data) * 0.3 else "black"
            ax.text(j, i, f"{val:.0f}" if val >= 10 else f"{val:.1f}",
                    ha="center", va="center", fontsize=7, color=tc, fontweight="bold")
ax.set_xticks(range(n_umi_bins))
ax.set_xticklabels(umi_bin_labels, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(zinb_pcts)))
ax.set_yticklabels([f"P{p}" for p in zinb_pcts], fontsize=9)
ax.set_xlabel("Total UMI per cell type")
ax.set_ylabel("Gene expression level (percentile)")
fig.colorbar(im, ax=ax, shrink=0.8, label="Median of max simulated specificity")
fig.savefig(FIG_DIR + "FigS4_C.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel D — ZINB: Expression Level × UMI Group → Max Specificity

# %%
umi_tercile_edges = np.percentile(total_umis_zinb, [0, 33, 67, 100])
umi_group_names = [f"Low UMI (<{umi_tercile_edges[1]:.0f})",
                   f"Mid UMI ({umi_tercile_edges[1]:.0f}-{umi_tercile_edges[2]:.0f})",
                   f"High UMI (>{umi_tercile_edges[2]:.0f})"]
umi_group_colors = ["#dc2626", "#f59e0b", "#2563eb"]

fig, ax = plt.subplots(figsize=(7, 5))
for g_idx in range(3):
    max_specs, p95_specs = [], []
    for pct in zinb_pcts:
        res = zinb_sweep[pct]
        mask = ((res["total_umi"] >= umi_tercile_edges[g_idx]) &
                (res["total_umi"] < umi_tercile_edges[g_idx + 1]))
        if g_idx == 2:
            mask = mask | (res["total_umi"] >= umi_tercile_edges[g_idx])
        max_specs.append(res.loc[mask, "max_spec"].max() if mask.sum() > 0 else np.nan)
        p95_specs.append(res.loc[mask, "p95_spec"].max() if mask.sum() > 0 else np.nan)
    ax.plot(zinb_pcts, max_specs, color=umi_group_colors[g_idx], marker="o",
            markersize=6, lw=2.5, label=umi_group_names[g_idx], zorder=5)
    ax.fill_between(zinb_pcts, p95_specs, max_specs, color=umi_group_colors[g_idx], alpha=0.1)
ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7, label="True specificity (1.0)")
ax.axhline(y=zinb_clip, color="black", ls=":", lw=1.5, alpha=0.7, label=f"Cap ({zinb_clip:.1f})")
ax.set_yscale("log")
ax.set_xlabel("Gene expression level (percentile)", fontsize=12)
ax.set_ylabel("Max simulated specificity", fontsize=12)
ax.legend(fontsize=8, framealpha=0.8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_D.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel E — ZINB: P95 Specificity vs UMI by Expression Level

# %%
from statsmodels.nonparametric.smoothers_lowess import lowess

show_pcts = [10, 25, 50, 75, 90, 95]
colors_pct = plt.cm.plasma(np.linspace(0.1, 0.9, len(show_pcts)))

fig, ax = plt.subplots(figsize=(7, 5))
for i, pct in enumerate(show_pcts):
    res = zinb_sweep[pct]
    ax.scatter(res["total_umi"], res["p95_spec"], color=colors_pct[i], alpha=0.3, s=12,
               edgecolors="none", label=f"P{pct}", zorder=3 + i)
    pos_mask = res["p95_spec"].values > 0.01
    if pos_mask.sum() > 10:
        log_umi = np.log10(res["total_umi"].values[pos_mask])
        log_spec = np.log10(res["p95_spec"].values[pos_mask])
        smooth = lowess(log_spec, log_umi, frac=0.4, return_sorted=True)
        ax.plot(10**smooth[:, 0], 10**smooth[:, 1], color=colors_pct[i], lw=2.5, zorder=10 + i)
ax.axhline(y=1.0, color="gray", ls="--", lw=1.5, alpha=0.7)
ax.axhline(y=zinb_clip, color="black", ls=":", lw=1.5, alpha=0.5, label=f"Cap ({zinb_clip:.1f})")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total UMI per cell type", fontsize=12)
ax.set_ylabel("P95 simulated specificity", fontsize=12)
ax.legend(fontsize=8, framealpha=0.8, title="Expr. level", title_fontsize=8, loc="upper right")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_E.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel F — Cap Sensitivity: Spearmans' R vs Cap=2x

# %%
# Load gene weights for cap sensitivity analysis
gw_files_cap = {
    "SCZ": str(PROJ_DIR / "dat/GeneWeights/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw"),
    "ASD (HIQ)": str(PROJ_DIR / "dat/GeneWeights/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    "DDD": str(PROJ_DIR / "dat/GeneWeights/DDD.top61.gw.bgmr.csv"),
}
gw_dicts_cap = {label: Fil2Dict(fpath) for label, fpath in gw_files_cap.items()}

CAP_LEVELS = [1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0]
ref_cap = 2.0

# Compute bias at each cap level
bias_all_cap = {}
for disorder, gw_dict in gw_dicts_cap.items():
    bias_all_cap[disorder] = {}
    for cap in CAP_LEVELS:
        threshold = mean_spec * cap
        spec_clipped = spec_unclip.clip(lower=0, upper=threshold)
        bias_df = HumanCT_AvgZ_Weighted(spec_clipped, gw_dict)
        bias_df = AnnotateCTDat(bias_df, Anno)
        bias_all_cap[disorder][cap] = bias_df

# Spearmans' R vs cap=2x
corr_results_cap = {}
for disorder in gw_dicts_cap.keys():
    ref_bias = bias_all_cap[disorder][ref_cap]["EFFECT"]
    corr_results_cap[disorder] = {}
    for cap in CAP_LEVELS:
        other_bias = bias_all_cap[disorder][cap]["EFFECT"]
        common = ref_bias.index.intersection(other_bias.index)
        rho, _ = sp_stats.spearmanr(ref_bias.loc[common], other_bias.loc[common])
        corr_results_cap[disorder][cap] = rho
corr_df_cap = pd.DataFrame(corr_results_cap)

# %%
disorder_colors = {"SCZ": "#2563eb", "ASD (HIQ)": "#dc2626", "DDD": "#16a34a"}
disorder_markers = {"SCZ": "o", "ASD (HIQ)": "s", "DDD": "D"}

fig, ax = plt.subplots(figsize=(7, 5))
ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)
for disorder in corr_df_cap.columns:
    ax.plot(CAP_LEVELS, corr_df_cap[disorder].values, color=disorder_colors[disorder],
            marker=disorder_markers[disorder], markersize=7, lw=2.2, label=disorder, zorder=3)
ax.set_xticks(CAP_LEVELS)
ax.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=9)
ax.set_xlabel("Specificity cap level (× mean)", fontsize=11)
ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
ax.set_ylabel("Spearmans' R vs. cap = 2×", fontsize=11)
ax.set_ylim(min(corr_df_cap.values.min() - 0.03, 0.35), 1.01)
ax.legend(fontsize=9, frameon=True, framealpha=0.9, loc="lower left")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_F.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel G — Cap Sensitivity: Supercluster Bias Profiles Across Cap Levels

# %%
panel_scs = ["Medium spiny neuron", "CGE interneuron", "MGE interneuron",
             "Vascular", "Microglia", "Astrocyte", "Ependymal"]
sc_style = {
    "CGE interneuron":      {"color": "#e11d48", "marker": "o", "ls": "-",  "lw": 2.5},
    "MGE interneuron":      {"color": "#7c3aed", "marker": "s", "ls": "-",  "lw": 2.5},
    "Medium spiny neuron":  {"color": "#0891b2", "marker": "D", "ls": "-",  "lw": 2.5},
    "Vascular":             {"color": "#64748b", "marker": "v", "ls": "--", "lw": 1.8},
    "Microglia":            {"color": "#92400e", "marker": "^", "ls": "--", "lw": 1.8},
    "Astrocyte":            {"color": "#d97706", "marker": "P", "ls": "--", "lw": 1.8},
    "Ependymal":            {"color": "#6b7280", "marker": "X", "ls": "--", "lw": 1.8},
}
sc_short = {"CGE interneuron": "CGE IN", "MGE interneuron": "MGE IN",
            "Medium spiny neuron": "MSN", "Vascular": "Vascular",
            "Microglia": "Microglia", "Astrocyte": "Astrocyte", "Ependymal": "Ependymal"}

# Compute supercluster mean bias across cap levels (SCZ)
sc_disorder_bias = {}
for sc in panel_scs:
    for disorder in ["SCZ"]:
        key = f"{sc} | {disorder}"
        sc_disorder_bias[key] = [bias_all_cap[disorder][cap].loc[
            bias_all_cap[disorder][cap]["Supercluster"] == sc, "EFFECT"].mean() for cap in CAP_LEVELS]
sc_disorder_df = pd.DataFrame(sc_disorder_bias, index=CAP_LEVELS)

fig, ax = plt.subplots(figsize=(7, 5))
ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)
for sc in panel_scs:
    key = f"{sc} | SCZ"
    sty = sc_style[sc]
    ax.plot(CAP_LEVELS, sc_disorder_df[key].values, color=sty["color"], marker=sty["marker"],
            markersize=5, lw=sty["lw"], ls=sty["ls"], label=sc_short[sc], zorder=3)
ax.set_xticks(CAP_LEVELS)
ax.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=9)
ax.set_xlabel("Specificity cap level (× mean)", fontsize=11)
ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
ax.set_ylabel("Mean supercluster bias", fontsize=11)
ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, loc="lower right", ncol=2)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_G.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel H — Cap Sensitivity: Supercluster Rank Stability

# %%
sc_rank_data = {}
for disorder in ["SCZ"]:
    sc_rank_data[disorder] = {}
    for cap in CAP_LEVELS:
        df = bias_all_cap[disorder][cap]
        sc_means = df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
        for rank_i, (sc_name, _) in enumerate(sc_means.items(), 1):
            if sc_name not in sc_rank_data[disorder]:
                sc_rank_data[disorder][sc_name] = {}
            sc_rank_data[disorder][sc_name][cap] = rank_i

fig, ax = plt.subplots(figsize=(7, 5))
ax.axvspan(1.5, 3.0, color="#e0f2fe", alpha=0.35, zorder=0)
for sc in panel_scs:
    sty = sc_style[sc]
    ranks = [sc_rank_data["SCZ"][sc][cap] for cap in CAP_LEVELS]
    ax.plot(CAP_LEVELS, ranks, color=sty["color"], marker=sty["marker"], markersize=5,
            lw=sty["lw"], ls=sty["ls"], label=sc_short[sc], zorder=3)
ax.set_xticks(CAP_LEVELS)
ax.set_xticklabels([f"{c}×" for c in CAP_LEVELS], fontsize=9)
ax.set_xlabel("Specificity cap level (× mean)", fontsize=11)
ax.axvline(x=ref_cap, color="gray", ls="--", lw=1, alpha=0.4, zorder=1)
ax.invert_yaxis()
ax.set_ylabel("Supercluster rank (1 = highest bias)", fontsize=11)
ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, loc="lower right", ncol=2)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.savefig(FIG_DIR + "FigS4_H.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ## Figure S5 — Impact of Number of Genes on Mutation Bias Analysis
#
# | Panel | Content | Source |
# |-------|---------|--------|
# | A | SCZ gene set size sweep | `Number_Gene_Effect.ipynb` |
# | B | ASD w/o ID gene set size sweep | `Number_Gene_Effect.ipynb` |
# | C | ASD with ID gene set size sweep | `Number_Gene_Effect.ipynb` |
# | D | DD/ID gene set size sweep | `Number_Gene_Effect.ipynb` |
# | E | Real vs random gene additions (SCZ, ASD w/ID, DDD) | `Number_Gene_Effect.ipynb` |

# %%
# Load pre-computed gene sweep data
import pickle
import matplotlib.ticker as mticker

with open(PLOT_DATA_DIR + "gene_sweep_data.pkl", "rb") as f:
    gene_sweep = pickle.load(f)
with open(PLOT_DATA_DIR + "gene_expansion_data.pkl", "rb") as f:
    gene_expansion = pickle.load(f)


def plot_gene_set_correlation_inline(GeneIdx, Corr, GeneSig, ylabel_corr, ax1=None,
                                     Corr_unweighted=None, color_corr='#1f77b4',
                                     color_pval='#d62728', color_unweighted='#2ca02c'):
    """Plot gene set correlation and -log10(p-value) vs number of genes."""
    if ax1 is None:
        fig, ax1 = plt.subplots(figsize=(6, 4.5), dpi=120)
    else:
        fig = ax1.figure
    ax1.plot(GeneIdx, Corr, color=color_corr, linewidth=2, marker='o', markersize=4, label="Bias Correlation")
    # if Corr_unweighted is not None:
    #     ax1.plot(GeneIdx, Corr_unweighted, color=color_unweighted, linewidth=2, marker='^',
    #              markersize=4, label="Bias Correlation (Unweighted)", linestyle='--')
    ax1.set_xlabel("Number of Genes", fontsize=12)
    ax1.set_ylabel(ylabel_corr, color=color_corr, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_corr)
    ax1.set_ylim(0, 1.01)
    ax1.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax1.grid(True, which='major', axis='y', linestyle='--', alpha=0.5, zorder=1)
    ax2 = ax1.twinx()
    ax2.plot(GeneIdx, GeneSig, color=color_pval, linewidth=2, marker='s', markersize=4,
             label=r"$-\log_{10}$(max p-value)", zorder=3)
    ax2.set_ylabel(r"$-\log_{10}$(max p-value)", color=color_pval, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_pval)
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower left', fontsize=10, frameon=False)
    fig.tight_layout(pad=2)
    return fig, ax1, ax2


# %% [markdown]
# ### Panel A–D — Gene Set Size Sweeps (2×2)

# %%
panel_configs = [
    ('SCZ', gene_sweep['SCZ'], "Bias Correlation with Main SCZ Set", "A"),
    ('ASD_woID', gene_sweep['ASD_woID'], "Bias Correlation with Main ASD w/o ID Set", "B"),
    ('ASD_wID', gene_sweep['ASD_wID'], "Bias Correlation with Main ASD with ID Set", "C"),
    ('DDD', gene_sweep['DDD'], "Bias Correlation with Main DD/ID Set", "D"),
]

fig, axes_grid = plt.subplots(2, 2, figsize=(12, 8), dpi=120)
axes_flat = axes_grid.flatten()

for ax1, (key, d, ylabel, panel_label) in zip(axes_flat, panel_configs):
    color_corr, color_pval, color_uw = '#1f77b4', '#d62728', '#2ca02c'
    GeneIdx = gene_sweep['GeneIdx']

    ax1.plot(GeneIdx, d['Corr'], color=color_corr, linewidth=2, marker='o', markersize=3, label="Bias Correlation")
    SHOW_UNWEIGHTED = False
    if SHOW_UNWEIGHTED and d.get('Corr_unweighted'):
        ax1.plot(GeneIdx, d['Corr_unweighted'], color=color_uw, linewidth=2, marker='^',
                 markersize=3, label="Bias Correlation (Unweighted)", linestyle='--')
    ax1.set_xlabel("Number of Genes", fontsize=14)
    ax1.set_ylabel(ylabel, color=color_corr, fontsize=13)
    ax1.tick_params(axis='y', labelcolor=color_corr)
    ax1.set_ylim(0, 1.01)
    ax1.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax1.grid(True, which='major', axis='y', linestyle='--', alpha=0.5, zorder=1)
    ax1.tick_params(axis='both', labelsize=12)

    ax2 = ax1.twinx()
    ax2.plot(GeneIdx, d['GeneSig'], color=color_pval, linewidth=2, marker='s', markersize=3,
             label=r"$-\log_{10}$(max p-value)", zorder=3)
    ax2.set_ylabel(r"$-\log_{10}$(max p-value)", color=color_pval, fontsize=13)
    ax2.tick_params(axis='y', labelcolor=color_pval, labelsize=12)

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower left', fontsize=10, frameon=False)

    ax1.text(-0.02, 1.05, panel_label, transform=ax1.transAxes,
             fontsize=22, fontweight='bold', va='bottom', ha='right')

fig.tight_layout(pad=1.0, h_pad=1.0, w_pad=2.0)
fig.patch.set_alpha(0)
fig.savefig(FIG_DIR + "FigS5_ABCD.png", dpi=300, bbox_inches="tight", transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel E–G — Real vs Random Gene Additions (SCZ, ASD with ID, DDD)

# %%
total_genes = gene_expansion['total_genes']
disorders_exp = gene_expansion['disorders']
disorder_order = [("SCZ", "#ff7f0e", "E"), ("ASD with ID", "#1f77b4", "F"), ("DDD", "#2ca02c", "G")]

fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=120, sharey=True)
for ax, (disorder_name, color, panel_label) in zip(axes, disorder_order):
    res = disorders_exp[disorder_name]
    ax.fill_between(total_genes, res["rand_lo"], res["rand_hi"],
                    color="#999999", alpha=0.25, label="Random genes (95% CI)")
    ax.plot(total_genes, res["rand_mean"], color="#999999", lw=2, ls="--",
            marker="s", markersize=4, label="Random genes (mean)")
    ax.plot(total_genes, res["real"], color=color, lw=2.5,
            marker="o", markersize=5, label=f"Ranked {disorder_name} genes", zorder=5)
    ax.axvline(61, color="gray", ls=":", lw=1.5, alpha=0.5)
    ax.set_xlabel("Total number of genes", fontsize=16)
    ax.set_title(disorder_name, fontweight="bold", fontsize=18)
    ax.text(-0.02, 1.05, panel_label, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='bottom', ha='right')
    ax.legend(fontsize=12, framealpha=0.8, loc="lower left")
    ax.set_ylim(0.3, 1.02)
    ax.tick_params(axis='both', labelsize=14)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
axes[0].set_ylabel("Spearmans' R with top-61 bias profile", fontsize=16)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(FIG_DIR + "FigS5_EFG.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ---
# ## Figure S6 — Impact of Genetic Architecture / Downsample Mutations
#
# **TODO:** Add panels from `Number_Gene_Effect.ipynb` (e.g., split-half bias comparison,
# sliding window correlation decay, downsampling stability).

# %%
# Placeholder — add S6 panels here
print("Figure S6: TODO — add genetic architecture / downsampling panels")

# %% [markdown]
# ---
# ## Figure S7 — Comprehensive Mutation Biases Across Brain Cell Types in Psychiatric Disorders
#
# | Panel | Content |
# |-------|---------|
# | A | ASD (all) |
# | B | ASD w/o ID |
# | C | ASD with ID |
# | D | SCZ |
# | E | DD/ID |

# %% [markdown]
# ### Panel A — ASD

# %%
with save_panel(FIG_DIR + "FigS7_A.png"):
    SuperClusterBias_BoxPlot(ASD_All_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.05)

# %% [markdown]
# ### Panel B — ASD w/o ID

# %%
with save_panel(FIG_DIR + "FigS7_B.png"):
    SuperClusterBias_BoxPlot(HighIQ_ASD_Bias, "ASD w/o ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.1)

# %% [markdown]
# ### Panel C — ASD with ID

# %%
with save_panel(FIG_DIR + "FigS7_C.png"):
    SuperClusterBias_BoxPlot(LowIQ_ASD_Bias, "ASD with ID", NeuroOnly=False, sortby="mean", EffectCol="-logP")

# %% [markdown]
# ### Panel D — SCZ

# %%
with save_panel(FIG_DIR + "FigS7_D.png"):
    SuperClusterBias_BoxPlot(SCZ_Bias, "SCZ", NeuroOnly=False, sortby="mean", EffectCol="-logP")

# %% [markdown]
# ### Panel E — DD/ID

# %%
with save_panel(FIG_DIR + "FigS7_E.png"):
    SuperClusterBias_BoxPlot(DDD_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP")


# %% [markdown]
# ---
# ## Figure S8 — Negative Controls (Non-Brain Traits) & SCZ Protective Genes
#
# | Panel | Content | Source |
# |-------|---------|--------|
# | A | Non-brain trait bias (HDL, IBD, ALT) | `NegativeControl_BiasPlot.ipynb` |
# | B | SCZ protective-direction (OR < 1) bias | `SCZ_Protective_BiasPlot.ipynb` |
#
# CGE signal is specific to psychiatric risk genes; absent in non-brain traits
# and inverted for protective-direction SCZ genes.

# %%
# Load negative control and SCZ protective bias+pval from main pipeline results
_matched_dir = str(PROJ_DIR / "results/main_results/matched_WB_mean_phastCons_n_CDS_bases_Best1000/Centering/") + "/"

HDL_Pval = pd.read_csv(_matched_dir + "NegCtrl_HDL_bias_addP.csv", index_col=0)
IBD_Pval = pd.read_csv(_matched_dir + "NegCtrl_IBD_bias_addP.csv", index_col=0)
ALT_Pval = pd.read_csv(_matched_dir + "NegCtrl_Alanine_bias_addP.csv", index_col=0)
SCZ_Protect_Pval = pd.read_csv(_matched_dir + "SCZ_protect_bias_addP.csv", index_col=0)

# Compute bias inline from gene weights (consistent with pipeline)
_gw_dir = str(PROJ_DIR / "dat/GeneWeights/") + "/"
_exp_mat = pd.read_csv(str(PROJ_DIR / _cfg['analysis_types']['Centering']), index_col=0)
_exp_mat.columns = _exp_mat.columns.astype(int)

HDL_Bias = AnnotateCTDat(HumanCT_AvgZ_Weighted(_exp_mat, Fil2Dict(_gw_dir + "NegCtrl_HDL.gw.csv")), Anno)
IBD_Bias = AnnotateCTDat(HumanCT_AvgZ_Weighted(_exp_mat, Fil2Dict(_gw_dir + "NegCtrl_IBD.gw.csv")), Anno)
ALT_Bias = AnnotateCTDat(HumanCT_AvgZ_Weighted(_exp_mat, Fil2Dict(_gw_dir + "NegCtrl_Alanine.gw.csv")), Anno)
SCZ_Protect_Bias = AnnotateCTDat(HumanCT_AvgZ_Weighted(_exp_mat, Fil2Dict(_gw_dir + "SCZ.top61.protect.gw")), Anno)

# %% [markdown]
# ### Panel A — Non-brain trait bias (EFFECT)

# %%
fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor="none")
for ax, (bias_df, name) in zip(axes, [(HDL_Bias, "HDL Cholesterol"), (IBD_Bias, "IBD"), (ALT_Bias, "Alanine AT")]):
    ax.patch.set_alpha(0)
    SuperClusterBias_BoxPlot(bias_df, name, ax=ax)
fig.patch.set_alpha(0)
plt.tight_layout()
fig.savefig(FIG_DIR + "FigS8_A_bias.png", dpi=300, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ### Panel B — Non-brain trait significance (-logP)

# %%
fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor="none")
for ax, (pval_df, name) in zip(axes, [(HDL_Pval, "HDL Cholesterol"), (IBD_Pval, "IBD"), (ALT_Pval, "Alanine AT")]):
    ax.patch.set_alpha(0)
    SuperClusterBias_BoxPlot(pval_df, name, EffectCol="-logP", fdr_cut=0.1, ax=ax)
fig.patch.set_alpha(0)
plt.tight_layout()
fig.savefig(FIG_DIR + "FigS8_B_pval.png", dpi=300, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ### Panel C — SCZ protective genes (EFFECT + significance)

# %%
fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="none")
axes[0].patch.set_alpha(0)
SuperClusterBias_BoxPlot(SCZ_Protect_Bias, "SCZ Protective (OR < 1)", ax=axes[0])
axes[1].patch.set_alpha(0)
SuperClusterBias_BoxPlot(SCZ_Protect_Pval, "SCZ Protective (OR < 1)", EffectCol="-logP", fdr_cut=0.1, ax=axes[1])
fig.patch.set_alpha(0)
plt.tight_layout()
fig.savefig(FIG_DIR + "FigS8_C_protect.png", dpi=300, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ---
# ## Figure S9 — Impact of Gene Expression Levels on ASD-SCZ Bias Correlation
#
# Pre-generated by Similarity_ASD_SCZ.spec notebook.

# %%
from multiprocessing import Pool

with open(str(PROJ_DIR / "config/config.yaml")) as file:
    config = yaml.safe_load(file)

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# Load expression matrix
expression_matrix = config['analysis_types']['Centering']
HCT_Z2_MAT = pd.read_csv(str(PROJ_DIR / expression_matrix), index_col=0)
HCT_Z2_MAT.columns = HCT_Z2_MAT.columns.astype(int)

# Gene weights
ASD_GW = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"))
SCZ_GW = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw"))

# BrainSpan expression
BrainSpan = pd.read_csv("/home/jw3514/Work/CellType_Psy/dat2/ExpMatch/BrainSpan.MatchDF.csv", index_col=0)

# gnomAD constraint (v4 primary, v2 fallback)
gnomad4 = pd.read_csv("/home/jw3514/Work/data/gnomad/gnomad.v4.0.constraint_metrics.tsv", sep="\t")
gnomad4 = gnomad4[(gnomad4["transcript"].str.contains('ENST')) & (gnomad4["mane_select"] == True)]
gnomad4["Entrez"] = gnomad4["gene"].map(GeneSymbol2Entrez).fillna(0).astype(int)
gnomad4 = gnomad4[gnomad4["Entrez"] != 0][["Entrez", "gene", "lof.pLI", "lof.z_score", "lof.oe_ci.upper"]].copy()

gnomad2 = pd.read_csv("/home/jw3514/Work/data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t")
gnomad2["Entrez"] = gnomad2["gene"].map(GeneSymbol2Entrez).fillna(0).astype(int)
gnomad2 = gnomad2[["Entrez", "gene", "pLI", "lof_z", "oe_lof_upper"]].copy()
gnomad2.columns = ["Entrez", "gene", "lof.pLI", "lof.z_score", "lof.oe_ci.upper"]
missing_in_v4 = set(gnomad2["Entrez"]) - set(gnomad4["Entrez"])
gnomad = pd.concat([gnomad4, gnomad2[gnomad2["Entrez"].isin(missing_in_v4)]], ignore_index=True)
gnomad = gnomad.drop_duplicates(subset="Entrez", keep="first").sort_values("lof.z_score", ascending=False)

# Annotate genes
def _annotate_gw(gw_dict, gnomad_df, brainspan_df):
    df = gnomad_df[gnomad_df["Entrez"].isin(gw_dict.keys())].copy()
    df["GW"] = df["Entrez"].map(gw_dict)
    df["BrainSpan"] = df["Entrez"].map(lambda x: brainspan_df.loc[x, "WB"] if x in brainspan_df.index else np.nan)
    return df

ASD_Genes = _annotate_gw(ASD_GW, gnomad, BrainSpan)
SCZ_Genes = _annotate_gw(SCZ_GW, gnomad, BrainSpan)

# %%
# Pre-compute numpy structures for fast bias correlation
expr_np = HCT_Z2_MAT.values
expr_gene_set = set(HCT_Z2_MAT.index)
expr_gene_to_row = {g: i for i, g in enumerate(HCT_Z2_MAT.index)}
ct_cols = HCT_Z2_MAT.columns.values
neur_col_mask = np.array([int(c) in Neur_idx for c in ct_cols])

def fast_bias_corr(asd_entrez, asd_weights, scz_entrez, scz_weights):
    asd_rows = np.array([expr_gene_to_row[g] for g in asd_entrez])
    asd_bias = np.average(expr_np[asd_rows], axis=0, weights=asd_weights)
    scz_rows = np.array([expr_gene_to_row[g] for g in scz_entrez])
    scz_bias = np.average(expr_np[scz_rows], axis=0, weights=scz_weights)
    from scipy.stats import spearmanr
    r, _ = spearmanr(asd_bias[neur_col_mask], scz_bias[neur_col_mask])
    return r

def prepare_gene_arrays(genes_df):
    mask = genes_df["Entrez"].isin(expr_gene_set)
    return genes_df.loc[mask, "Entrez"].values, genes_df.loc[mask, "GW"].values

N_REMOVAL_STEPS = 31

def compute_removal_curve(asd_sorted_df, scz_sorted_df):
    asd_entrez, asd_weights = prepare_gene_arrays(asd_sorted_df)
    scz_entrez, scz_weights = prepare_gene_arrays(scz_sorted_df)
    return [fast_bias_corr(asd_entrez[i:], asd_weights[i:], scz_entrez[i:], scz_weights[i:]) for i in range(N_REMOVAL_STEPS)]

# %%
# BrainSpan removal curves
ASD_by_BS = ASD_Genes.dropna(subset=["BrainSpan"]).sort_values("BrainSpan", ascending=False)
SCZ_by_BS = SCZ_Genes.dropna(subset=["BrainSpan"]).sort_values("BrainSpan", ascending=False)
Y_BS_high_first = compute_removal_curve(ASD_by_BS, SCZ_by_BS)

ASD_by_BS_rev = ASD_Genes.dropna(subset=["BrainSpan"]).sort_values("BrainSpan", ascending=True)
SCZ_by_BS_rev = SCZ_Genes.dropna(subset=["BrainSpan"]).sort_values("BrainSpan", ascending=True)
Y_BS_low_first = compute_removal_curve(ASD_by_BS_rev, SCZ_by_BS_rev)

# Load cached random null
CACHE_FILE = PROJ_DIR / "dat/Other/ASD_SCZ_RandomGeneRemoval_Null_v2.npy"
RandNull = np.load(CACHE_FILE)
rand_mean = RandNull.mean(axis=0)
rand_std = RandNull.std(axis=0)

X = list(range(N_REMOVAL_STEPS))

# %%
fig, ax = plt.subplots(dpi=150, figsize=(9.5, 6), facecolor='none')
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

ax.plot(X, Y_BS_low_first, label="Remove lowest expressed genes first",
        color="red", linestyle='-', marker='o', markersize=8,
        markeredgecolor='black', markeredgewidth=1, zorder=10)
ax.plot(X, Y_BS_high_first, label="Remove highest expressed genes first",
        color="blue", linestyle='--', marker='s', markersize=8,
        markeredgecolor='black', markeredgewidth=1, zorder=10)
ax.errorbar(X, rand_mean, yerr=rand_std, fmt='-', color="grey", ecolor='grey',
            elinewidth=2, capsize=4, capthick=2, label="Random removal", zorder=5)

ax.set_xlabel("Number of Genes Removed", fontsize=25)
ax.set_ylabel("Mutation Bias Correlation", fontsize=25)
ax.legend(fontsize=18, loc='best', frameon=False)
ax.tick_params(axis='both', labelsize=15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
fig.savefig(FIG_DIR + "FigureS9.png", dpi=300, transparent=True, bbox_inches='tight')
fig.savefig(FIG_DIR + "FigureS9.pdf", dpi=300, transparent=True, bbox_inches='tight')
plt.show()
print(f"Saved: {FIG_DIR}FigureS9.pdf")

# %% [markdown]
# ---
# ## Figure S10 — ASD and SCZ Mutation Bias Towards Different Cell Type Superclusters
#
# Individual cell type comparisons across disorders.
#
# | Panel | Content |
# |-------|---------|
# | A | Hippo CA1-3 (ASD w/o ID vs SCZ) |
# | B | Upper IT (ASD w/o ID vs SCZ) |
# | C | Deep IT (ASD w/o ID vs SCZ) |
# | D | Deep CT6b (ASD w/o ID vs SCZ) |
# | E | ASD w/o ID vs ASD with ID bar plot |
# | F | CGE: ASD with ID vs SCZ |
# | F2 | MGE: ASD with ID vs SCZ |
# | G | CGE: VNR+ vs VNR- |
# | H | CGE: ASD w/o ID vs DD/ID |

# %% [markdown]
# ### Panel A — Hippocampal CA1-3

# %%
with save_panel(FIG_DIR + "FigS10_A.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Hippocampal CA1-3", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.15, 0.12))

# %% [markdown]
# ### Panel B — Upper-layer IT

# %%
with save_panel(FIG_DIR + "FigS10_B.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Upper-layer intratelencephalic", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, pval="Mann_Whitney_FDR", loc=(0.15, 0.12))

# %% [markdown]
# ### Panel C — Deep-layer IT

# %%
with save_panel(FIG_DIR + "FigS10_C.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer intratelencephalic", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.3), pval="Mann_Whitney_FDR")

# %% [markdown]
# ### Panel D — Deep-layer CT6b

# %%
with save_panel(FIG_DIR + "FigS10_D.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer corticothalamic and 6b", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.1, 0.05))

# %% [markdown]
# ### Panel E — ASD w/o ID vs ASD with ID bar plot

# %%
SCZ_ASD_wID_Contrast_Neurons = SCZ_ASD_wID_Contrast[SCZ_ASD_wID_Contrast.index.isin(Neurons)]

# %%
with save_panel(FIG_DIR + "FigS10_E.png"):
    plot_bias_comparison(SCZ_ASD_wID_Contrast_Neurons, "SCZ", "ASD with ID",
                         p_test="Mann_Whitney_FDR", legend_anchor=(0.9, 1.0))

# %% [markdown]
# ### Panel F — CGE: ASD with ID vs SCZ

# %%
with save_panel(FIG_DIR + "FigS10_F.png"):
    CompareSingleCT(LowIQ_ASD_Bias, SCZ_Bias, "CGE interneuron", ASD_wID_SCZ_Contrast,
                    "ASD with ID Mutation Bias", "SCZ Mutation Bias", loc=(0.08, 0.21))

# %% [markdown]
# ### Panel F2 — MGE: ASD with ID vs SCZ

# %%
with save_panel(FIG_DIR + "FigS10_F2.png"):
    CompareSingleCT(LowIQ_ASD_Bias, SCZ_Bias, "MGE interneuron", ASD_wID_SCZ_Contrast,
                    "ASD with ID Mutation Bias", "SCZ Mutation Bias", loc=(0.08, 0.21))

# %% [markdown]
# ### Panel G — CGE: VNR+ vs VNR-

# %%
with save_panel(FIG_DIR + "FigS10_G.png"):
    CompareSingleCT(VNR_Pos_Bias, VNR_Neg_Bias, "CGE interneuron", VNR_Contrast,
                    "VNR + Mutation Bias", "VNR - Mutation Bias", loc=(0.0, -0.05))

# %% [markdown]
# ### Panel H — CGE: ASD w/o ID vs DD/ID

# %%
with save_panel(FIG_DIR + "FigS10_H.png"):
    CompareSingleCT(HighIQ_ASD_Bias, DDD_Bias, "CGE interneuron", ASD_DDD_Contrast,
                    "ASD w/o ID Mutation Bias", "DD/ID Mutation Bias", loc=(0.1, 0.05))


# %% [markdown]
# ---
# ## Figure S11 — Mutation Bias Comparison Across Psychiatric Disorders for MGE & LAMP5-LHX6/Chandelier
#
# | Panel | Content |
# |-------|---------|
# | A | MGE interneuron |
# | B | LAMP5-LHX6 and Chandelier |

# %%
datasets = {
    'ASD w/o ID': HighIQ_ASD_Bias,
    'ASD with ID': LowIQ_ASD_Bias,
    'VNR+': VNR_Pos_Bias,
    'VNR-': VNR_Neg_Bias,
    'DD/ID': DDD_Bias,
    'SCZ': SCZ_Bias
}
TestPairs = [("VNR+", "VNR-"), ("ASD w/o ID", "SCZ"), ("ASD w/o ID", "ASD with ID"),
             ("SCZ", "ASD with ID"), ("ASD w/o ID", "DD/ID")]

# %% [markdown]
# ### Panel A — MGE interneuron

# %%
with save_panel(FIG_DIR + "FigS11_A.png"):
    plot_mutation_bias_comparison_V2("MGE interneuron", datasets, Anno, all_contrasts_df, TestPairs=TestPairs)

# %% [markdown]
# ### Panel B — LAMP5-LHX6 and Chandelier

# %%
with save_panel(FIG_DIR + "FigS11_B.png"):
    plot_mutation_bias_comparison_V2("LAMP5-LHX6 and Chandelier", datasets, Anno, all_contrasts_df, TestPairs=TestPairs)


# %% [markdown]
# ---
# ## Figure S12 — Comprehensive Analysis of Mutation Biases Across Brain Cell Types in 22q11.2 Deletion

# %%
with save_panel(FIG_DIR + "FigureS12.png"):
    SuperClusterBias_BoxPlot(X22q_Bias, "22q11.2", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.1)

img = PILImage.open(FIG_DIR + "FigureS12.png")
img.save(FIG_DIR + "FigureS12.pdf", "PDF", resolution=300)
print(f"Saved: {FIG_DIR}FigureS12.pdf")

# %% [markdown]
# ---
# # Supplementary Tables

# %% [markdown]
# ---
# ## Supplementary Tables S2-S7: Cluster-level biases

# %%
import openpyxl

def prepare_bias_table_v1(df):
    """Prepare a bias DataFrame for supplementary table output."""
    df_out = df.copy(deep=True)
    df_out.index.name = "Cluster"
    columns_to_keep = [
        "EFFECT", "P-value", "q-value", "Class", "Supercluster", "Subtype", "Neurotransmitter",
        "Top three regions", "Top three dissections", "Number of cells"]
    df_out = df_out[columns_to_keep]
    df_out.rename(columns={"EFFECT": "Bias"}, inplace=True)
    df_out = df_out.sort_values(by=["P-value", "Bias"], ascending=[True, False])
    return df_out

def prepare_bias_table_v2(BiasPos, BiasNeg, Name1, Name2):
    df_out = pd.DataFrame()
    df_out["Bias {}".format(Name1)] = BiasPos["EFFECT"]
    df_out["P-value {}".format(Name1)] = BiasPos["P-value"]
    df_out["q-value {}".format(Name1)] = BiasPos["q-value"]
    df_out["Bias {}".format(Name2)] = BiasNeg["EFFECT"]
    df_out["P-value {}".format(Name2)] = BiasNeg["P-value"]
    df_out["q-value {}".format(Name2)] = BiasNeg["q-value"]
    df_out["Class"] = BiasPos["Class"]
    df_out["Supercluster"] = BiasPos["Supercluster"]
    df_out["Subtype"] = BiasPos["Subtype"]
    df_out["Neurotransmitter"] = BiasPos["Neurotransmitter"]
    df_out["Top three regions"] = BiasPos["Top three regions"]
    df_out["Top three dissections"] = BiasPos["Top three dissections"]
    df_out["Number of cells"] = BiasPos["Number of cells"]
    df_out.index.name = "Cluster"
    df_out = df_out.sort_values(by=["P-value {}".format(Name2), "Bias {}".format(Name2)], ascending=[True, False])
    return df_out

SCZ_Bias_toST = prepare_bias_table_v1(SCZ_Bias)
HighIQ_ASD_Bias_toST = prepare_bias_table_v1(HighIQ_ASD_Bias)
LowIQ_ASD_Bias_toST = prepare_bias_table_v1(LowIQ_ASD_Bias)
X22q_Bias_toST = prepare_bias_table_v1(X22q_Bias)
DDD_Bias_toST = prepare_bias_table_v1(DDD_Bias)
UKBB_VNR_Bias_toST = prepare_bias_table_v2(VNR_Pos_Bias, VNR_Neg_Bias, "VNR+", "VNR-")

# %%
SuppTabOutDir = str(PROJ_DIR / "dat/suppl.data/") + "/"
excel_path = SuppTabOutDir + "SupTab_VIP.xlsx"

# Remove existing data sheets (keep Table_of_contents and experiment sheets)
wb = openpyxl.load_workbook(excel_path)
sheets_to_remove = [s for s in wb.sheetnames if s not in [
    "Table_of_contents",
    "Mouse Experiment Statistical detail  Statistical detail summary.",
    "Experiement Statistical detail "
]]
for s in sheets_to_remove:
    wb.remove(wb[s])
wb.save(excel_path)
wb.close()

# %%
with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    Name_Dict = {
        "SCZ_Bias": "Table_S2_Cluster_Bias_SCZ",
        "HighIQ_ASD_Bias": "Table_S3_Cluster_Bias_ASD_woID",
        "LowIQ_ASD_Bias": "Table_S4_Cluster_Bias_ASD_ID",
        "X22q_Bias": "Table_S5_Cluster_Bias_22q11.2",
        "DDD_Bias": "Table_S6_Cluster_Bias_DD",
        "VNR_Bias": "Table_S7_Cluster_Bias_VNR",
    }
    DF_list = [
        ("SCZ_Bias", SCZ_Bias_toST),
        ("HighIQ_ASD_Bias", HighIQ_ASD_Bias_toST),
        ("LowIQ_ASD_Bias", LowIQ_ASD_Bias_toST),
        ("X22q_Bias", X22q_Bias_toST),
        ("DDD_Bias", DDD_Bias_toST),
        ("VNR_Bias", UKBB_VNR_Bias_toST),
    ]
    for df_name, DF in DF_list:
        sheet_name = Name_Dict.get(df_name, df_name)
        DF.to_excel(writer, sheet_name=sheet_name)

print("Tables S2-S7 written.")

# %% [markdown]
# ---
# ## Supplementary Tables S8-S11: SuperCluster bias contrasts

# %%
def process_BiasContrast_df(df):
    columns_to_drop = ["Wilcoxon_P", "Wilcoxon_FDR", "Bonferroni_P"]
    df = df.drop(columns=[col for col in columns_to_drop if col in df.columns])
    disorder_name_1 = df.columns[0].replace("Bias_", "")
    disorder_name_2 = df.columns[1].replace("Bias_", "")
    df.rename(columns={"Bias_Diff": "Bias_Diff_{}_{}".format(disorder_name_1, disorder_name_2)}, inplace=True)
    return df

ASD_SCZ_Contrast_toST = process_BiasContrast_df(ASD_SCZ_Contrast_Neurons)
HIQ_LIQ_Contrast_toST = process_BiasContrast_df(HIQ_LIQ_Contrast_Neurons)
DDD_ASD_Contrast_toST = process_BiasContrast_df(ASD_DDD_Contrast_Neurons)
VNR_Contrast_toST = process_BiasContrast_df(VNR_Contrast_Neurons)

# %%
with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    Name_Dict = {
        "BiasContrast_ASD_SCZ": "Table_S8_BiasContrast_ASD_SCZ",
        "BiasContrast_HIQ_LIQ": "Table_S9_BiasContrast_HIQ_LIQ",
        "BiasContrast_HIQ_DDD": "Table_S10_BiasContrast_HIQ_DDD",
        "BiasContrast_VNR": "Table_S11_BiasContrast_VNR",
    }
    DF_list = [
        ("BiasContrast_ASD_SCZ", ASD_SCZ_Contrast_toST),
        ("BiasContrast_HIQ_LIQ", HIQ_LIQ_Contrast_toST),
        ("BiasContrast_HIQ_DDD", DDD_ASD_Contrast_toST),
        ("BiasContrast_VNR", VNR_Contrast_toST),
    ]
    for df_name, DF in DF_list:
        sheet_name = Name_Dict.get(df_name, df_name)
        DF.to_excel(writer, sheet_name=sheet_name)

print("Tables S8-S11 written.")
