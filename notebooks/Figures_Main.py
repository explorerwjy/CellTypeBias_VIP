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
# # Main Manuscript Figures
#
# Figures 2, 3, 4.
#
# **This notebook only loads pre-computed data and plots.**
# Run `Compute_BiasContrasts.ipynb` first to generate the required CSV files.
#
# Each figure section produces individual panels, then assembles them into a
# multi-panel composite PDF at the end.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os
import io
from contextlib import contextmanager

from pathlib import Path
import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *
from PIL import Image as PILImage
from matplotlib.image import imread

import matplotlib.font_manager as fm
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)
fm._load_fontmanager(try_read_cache=False)
plt.style.use('seaborn-v0_8-whitegrid')

FIG_DIR = str(PROJ_DIR / "results/figures/main/") + "/"
os.makedirs(FIG_DIR, exist_ok=True)


# %%
# -- Panel capture utility --
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
        # Display the saved image inline
        from IPython.display import Image as IPImage, display as ipdisplay
        ipdisplay(IPImage(filename=filepath))


def assemble_panels(panel_paths, panel_labels, layout, figsize, out_path, dpi=300):
    """Assemble saved panel PNGs into a multi-panel composite figure.

    Parameters
    ----------
    panel_paths : list of str
        Paths to individual panel PNG files.
    panel_labels : list of str
        Labels for each panel (e.g., "A", "B", ...).
    layout : list of (row, col_slice)
        GridSpec placement for each panel.  Each entry is
        ``(row_index, col_slice)`` where *col_slice* is a ``slice`` object.
    figsize : tuple
        Figure size in inches.
    out_path : str
        Output PDF/PNG path.
    dpi : int
        Resolution.
    """
    # Determine grid dimensions from layout
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
X22q_mousemodel = pd.read_csv(Bias_Save_Dir + "22q_small_del_bias_addP.csv", index_col=0)

# SCZ protective-direction genes (OR < 1)
SCZ_Protect_Bias = pd.read_csv(
    str(Bias_Save_Dir + "/SCZ_protect_bias_addP.csv"),
    index_col=0
)

# %%
# Pre-computed contrasts
CONTRAST_DIR = str(PROJ_DIR / "results/main_results/contrasts/") + "/"

ASD_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_SCZ_contrast.csv", index_col=0)
ASD_wID_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_wID_vs_SCZ_contrast.csv", index_col=0)
HIQ_LIQ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_ASD_wID_contrast.csv", index_col=0)
ASD_DDD_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_DDD_contrast.csv", index_col=0)
SCZ_ASD_wID_Contrast = pd.read_csv(CONTRAST_DIR + "SCZ_vs_ASD_wID_contrast.csv", index_col=0)
VNR_Contrast = pd.read_csv(CONTRAST_DIR + "VNR_neg_vs_pos_contrast.csv", index_col=0)
all_contrasts_df = pd.read_csv(CONTRAST_DIR + "all_contrasts_fdr.csv")
cluster_biasdiff = pd.read_csv(CONTRAST_DIR + "cluster_biasdiff.csv")
cluster_biasdiff_tests = pd.read_csv(CONTRAST_DIR + "cluster_biasdiff_tests.csv")
CGE_anno_22q = pd.read_csv(CONTRAST_DIR + "CGE_VIP_annotation_22q.csv", index_col=0)

# Neuron-filtered contrasts
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
HIQ_LIQ_Contrast_Neurons = HIQ_LIQ_Contrast[HIQ_LIQ_Contrast.index.isin(Neurons)]

# Compute SCZ vs SCZ Protective contrast and append to all_contrasts_df
SCZ_Protect_Contrast = compare_biases(SCZ_Bias, SCZ_Protect_Bias, name1="SCZ", name2="SCZ Protective", neurons=Neurons)
_protect_rows = []
for supercluster, row in SCZ_Protect_Contrast.iterrows():
    _protect_rows.append({
        "Pair": "SCZ - SCZ Protective",
        "SuperCluster": supercluster,
        "Bias1": row["Bias_SCZ"],
        "Bias2": row["Bias_SCZ Protective"],
        "BiasDiff": row["Bias_Diff"],
        "Pval": row["Mann_Whitney_P"],
        "MWU_FDR": row["Mann_Whitney_FDR"],
    })
_protect_df = pd.DataFrame(_protect_rows)
_protect_df["ALL_FDR"] = _protect_df["MWU_FDR"]  # Will be approximate; not pooled with other tests
all_contrasts_df = pd.concat([all_contrasts_df, _protect_df], ignore_index=True)

EffLabel = "EFFECT"

# %% [markdown]
# ---
# ## Figure 2
#
# | Panel | Content |
# |-------|---------|
# | A | ASD w/o ID vs SCZ — bias correlation scatter |
# | B | LOEUF gene removal — bias correlation drop |
# | C | ASD w/o ID vs SCZ — neuron bias bar plot |
# | D | Medium spiny neuron scatter |
# | E | MGE interneuron scatter |
# | F | CGE interneuron scatter (ASD w/o ID vs SCZ) |
# | G | CGE interneuron scatter (ASD w/o ID vs ASD with ID) |

# %% [markdown]
# ### Panel A — ASD w/o ID vs SCZ Bias Correlation

# %%
HIQ_ASD_SCZ, _ = CompareCT(HighIQ_ASD_Bias, SCZ_Bias, "ASD w/o ID", "SCZ",
                            effectlabel="EFFECT", SuperClusters=ALL_CTs)

with save_panel(FIG_DIR + "Fig2_A.png"):
    PlotBiasContrast(HIQ_ASD_SCZ,
                     label1="EFFECT_ASD w/o ID", label2="EFFECT_SCZ",
                     name1="ASD w/o ID Mutation Bias", name2="SCZ Mutation Bias",
                     title="", neur_only=True)

# %%
LIQ_ASD_SCZ, _ = CompareCT(LowIQ_ASD_Bias, SCZ_Bias, "ASD w/o ID", "SCZ",
                            effectlabel="EFFECT", SuperClusters=ALL_CTs)

#with save_panel(FIG_DIR + "FigSX.png"):
PlotBiasContrast(LIQ_ASD_SCZ,
                    label1="EFFECT_ASD w/o ID", label2="EFFECT_SCZ",
                    name1="ASD w/o ID Mutation Bias", name2="SCZ Mutation Bias",
                    title="", neur_only=True)

# %% [markdown]
# ### Panel B — LOEUF Gene Removal Correlation
#
# Pre-generated by Similarity_ASD_SCZ.spec notebook. Copy into figure panel dir.

# %%
import shutil
LOEUF_SRC = str(PROJ_DIR / "results/figures/LOEUF_gene_removal_ASD_SCZ_correlation.png")
shutil.copy2(LOEUF_SRC, FIG_DIR + "Fig2_B.png")

# %% [markdown]
# ### Panel C — ASD w/o ID vs SCZ Neuron Bias Comparison

# %%
with save_panel(FIG_DIR + "Fig2_C.png"):
    plot_bias_comparison(ASD_SCZ_Contrast_Neurons, "ASD w/o ID", "SCZ",
                         p_test="Mann_Whitney_FDR", legend_anchor=(0.01, 1.1))

# %%
with save_panel(FIG_DIR + "Fig2_D_bias_comparison.png"):
    plot_bias_comparison(HIQ_LIQ_Contrast_Neurons, "ASD w/o ID", "ASD with ID",
                         p_test="Mann_Whitney_FDR", legend_anchor=(0.65, 1.1))

# %% [markdown]
# ### Panel D — Medium Spiny Neuron

# %%
with save_panel(FIG_DIR + "Fig2_D.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Medium spiny neuron", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias",
                    efflabel=EffLabel, pval="Mann_Whitney_FDR", loc=(0.05, 0.23))

# %% [markdown]
# ### Panel E — MGE Interneuron

# %%
with save_panel(FIG_DIR + "Fig2_E.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "MGE interneuron", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.05))

# %% [markdown]
# ### Panel F — CGE Interneuron (ASD w/o ID vs SCZ)

# %%
with save_panel(FIG_DIR + "Fig2_F.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "CGE interneuron", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.05))

# %% [markdown]
# ### Panel G — CGE Interneuron (ASD w/o ID vs ASD with ID)

# %%
with save_panel(FIG_DIR + "Fig2_G.png"):
    CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "CGE interneuron",
                    HIQ_LIQ_Contrast_Neurons,
                    "ASD w/o ID Mutation Bias", "ASD with ID\nMutation Bias",
                    loc=(0.12, 0.05), efflabel=EffLabel)

# %% [markdown]
# ### Panel D — ASD with ID vs ASD w/o ID Neuron Bias Comparison

# %%

# %% [markdown]
# ---
# ## Figure 3
#
# | Panel | Content |
# |-------|---------|
# | A | Multi-disorder CGE comparison with FDR |
# | B | Cluster-level bias differences across interneuron superclusters |

# %% [markdown]
# ### Panel A — Multi-disorder CGE comparison

# %%
datasets = {
    'ASD w/o ID': HighIQ_ASD_Bias,
    'ASD with ID': LowIQ_ASD_Bias,
    'VNR+': VNR_Pos_Bias,
    'VNR-': VNR_Neg_Bias,
    'DD/ID': DDD_Bias,
    'SCZ': SCZ_Bias,
    'SCZ Protective': SCZ_Protect_Bias,
}
TestPairs = [("VNR+", "VNR-"), ("ASD w/o ID", "SCZ"), ("ASD w/o ID", "ASD with ID"),
             ("SCZ", "ASD with ID"), ("ASD w/o ID", "DD/ID"), ("SCZ", "SCZ Protective")]

with save_panel(FIG_DIR + "Fig3_A.png"):
    plot_mutation_bias_comparison_V2("CGE interneuron", datasets, Anno,
                                     all_contrasts_df, TestPairs=TestPairs)

# %% [markdown]
# ### Panel B — Cluster-level bias differences

# %%
from scipy.stats import mannwhitneyu
import matplotlib.patches as mpatches

plot_df = cluster_biasdiff
test_df = cluster_biasdiff_tests

superclusters = ["CGE interneuron", "MGE interneuron", "LAMP5-LHX6 and Chandelier"]
disorder_pairs = [
    ("ASD with ID", "SCZ"),
    ("SCZ", "ASD w/o ID"),
    ("ASD with ID", "ASD w/o ID"),
    ("DD/ID", "ASD w/o ID"),
    ("VNR-", "VNR+"),
]

def disorderpair_label(d1, d2):
    return f"{d1} to {d2}"

disorderpair_order = [disorderpair_label(d1, d2) for d1, d2 in disorder_pairs]
supercluster_order = superclusters

def align_disorderpair_direction(df, disorder_pairs, value_cols=()):
    """Return rows in requested pair directions, flipping signed values as needed."""
    aligned_rows = []
    available_pairs = set(df["DisorderPair"])
    for d1, d2 in disorder_pairs:
        label = disorderpair_label(d1, d2)
        reverse_label = disorderpair_label(d2, d1)
        if label in available_pairs:
            pair_rows = df[df["DisorderPair"] == label].copy()
        elif reverse_label in available_pairs:
            pair_rows = df[df["DisorderPair"] == reverse_label].copy()
            pair_rows["DisorderPair"] = label
            for col in value_cols:
                if col in pair_rows:
                    pair_rows[col] = -pair_rows[col]
        else:
            raise ValueError(f"Missing disorder pair '{label}' and reverse pair '{reverse_label}'")
        aligned_rows.append(pair_rows)
    return pd.concat(aligned_rows, ignore_index=True)

plot_df = align_disorderpair_direction(plot_df, disorder_pairs, value_cols=["BiasDiff"])
test_df = align_disorderpair_direction(
    test_df,
    disorder_pairs,
    value_cols=["CGE_median", "MGE_median", "LAMP5_median"],
)

fig, ax = plt.subplots(figsize=(8, 6), dpi=300, facecolor='none')
fig.patch.set_alpha(0.0)
ax.set_facecolor('none')

palette = sns.color_palette("Set2", len(supercluster_order))
color_dict = dict(zip(supercluster_order, palette))

sns.boxplot(
    data=plot_df, x="DisorderPair", y="BiasDiff", hue="Supercluster",
    order=disorderpair_order, hue_order=supercluster_order, dodge=True,
    showfliers=False, width=0.6,
    boxprops=dict(alpha=0.7), medianprops={"color": "black"},
    whiskerprops={"color": "black"}, ax=ax, palette=color_dict
)

# Calculate box positions for scatter overlay
n_pairs = len(disorderpair_order)
n_super = len(supercluster_order)
xticks = np.arange(n_pairs)
width = 0.6
total_width = width
each_width = total_width / n_super
offsets = np.linspace(-total_width/2 + each_width/2, total_width/2 - each_width/2, n_super)
pos_dict = {}
for i, dp in enumerate(disorderpair_order):
    for j, sc in enumerate(supercluster_order):
        pos_dict[(dp, sc)] = xticks[i] + offsets[j]

for (dp, sc), group in plot_df.groupby(["DisorderPair", "Supercluster"]):
    if (dp, sc) not in pos_dict:
        continue
    x = pos_dict[(dp, sc)]
    y = group["BiasDiff"].values
    jitter = np.random.uniform(-each_width/4, each_width/4, size=len(y))
    ax.scatter(np.full_like(y, x) + jitter, y, color=color_dict[sc],
               edgecolor="black", alpha=0.8, linewidth=0.7, s=60, label=None)

ax.set_ylabel("Bias Difference Group1 - Group2", fontsize=15)
ax.set_xlabel("")
ax.axhline(0, color='gray', linestyle='--', linewidth=1)
ax.set_xticks(xticks)
ax.set_xticklabels(disorderpair_order, rotation=30, ha='right', rotation_mode='anchor', fontsize=15)
ax.tick_params(axis='x', which='major', bottom=True, length=4, width=1,
               direction='out', color='black', pad=8)

ax.spines['left'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Annotate p-values from pre-computed test results
cap_height = 0.01
for i, dp in enumerate(disorderpair_order):
    cge_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "CGE interneuron")]["BiasDiff"].values
    mge_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "MGE interneuron")]["BiasDiff"].values
    lamp_diffs = plot_df[(plot_df["DisorderPair"] == dp) & (plot_df["Supercluster"] == "LAMP5-LHX6 and Chandelier")]["BiasDiff"].values

    row_cge_mge = test_df[(test_df["DisorderPair"] == dp) & (test_df["Test"] == "CGE != MGE")]
    if not row_cge_mge.empty and len(cge_diffs) > 0 and len(mge_diffs) > 0:
        pval = row_cge_mge["MWU_pval"].values[0]
        if pval < 0.05:
            x1 = pos_dict[(dp, "CGE interneuron")]
            x2 = pos_dict[(dp, "MGE interneuron")]
            y_max = max(np.max(cge_diffs), np.max(mge_diffs))
            y = y_max + 0.03
            ax.plot([x1, x2], [y, y], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.plot([x1, x1], [y, y - cap_height], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.plot([x2, x2], [y, y - cap_height], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.text((x1+x2)/2, y+0.01, f"{format_pval_scientific(pval)}", ha='center', va='bottom', fontsize=10, zorder=11)

    row_cge_lamp = test_df[(test_df["DisorderPair"] == dp) & (test_df["Test"] == "CGE != LAMP5")]
    if not row_cge_lamp.empty and len(cge_diffs) > 0 and len(lamp_diffs) > 0:
        pval = row_cge_lamp["MWU_pval"].values[0]
        if pval < 0.05:
            x1 = pos_dict[(dp, "CGE interneuron")]
            x2 = pos_dict[(dp, "LAMP5-LHX6 and Chandelier")]
            y_max = max(np.max(cge_diffs), np.max(lamp_diffs))
            y = y_max + 0.07
            ax.plot([x1, x2], [y, y], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.plot([x1, x1], [y, y-cap_height], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.plot([x2, x2], [y, y-cap_height], color='k', lw=1.2, zorder=10, ls='--', alpha=0.7)
            ax.text((x1+x2)/2, y+0.01, f"{format_pval_scientific(pval)}", ha='center', va='bottom', fontsize=10, zorder=11)

plt.tight_layout()
handles, labels = ax.get_legend_handles_labels()
seen = set()
new_handles, new_labels = [], []
for h, l in zip(handles, labels):
    if l not in seen and l in supercluster_order:
        new_handles.append(h)
        new_labels.append(l)
        seen.add(l)
ax.legend(new_handles, new_labels, loc='upper left', bbox_to_anchor=(0.005, 1.08), borderaxespad=0., fontsize=15)
plt.ylim(-0.1, 0.35)

fig.savefig(FIG_DIR + "Fig3_B.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Composite Figure 3

# %%
assemble_panels(
    panel_paths=[FIG_DIR + "Fig3_A.png", FIG_DIR + "Fig3_B.png"],
    panel_labels=["A", "B"],
    layout=[(0, slice(0, 1)), (0, slice(1, 2))],
    figsize=(18, 8),
    out_path=FIG_DIR + "Figure3.pdf",
)

# %% [markdown]
# ---
# ## Figure 4
#
# | Panel | Content |
# |-------|---------|
# | A | 22q11.2: CGE vs LAMP5 vs MGE |
# | B | 22q11.2: VIP+ vs VIP- CGE interneurons |

# %% [markdown]
# ### Panel A — 22q11.2: CGE vs LAMP5 vs MGE

# %%
from scipy.stats import mannwhitneyu

X22q_Bias_CGE = X22q_Bias[X22q_Bias["Supercluster"] == "CGE interneuron"]
X22q_Bias_MGE = X22q_Bias[X22q_Bias["Supercluster"] == "MGE interneuron"]
X22q_Bias_LAMP5 = X22q_Bias[X22q_Bias["Supercluster"] == "LAMP5-LHX6 and Chandelier"]

labels_4a = ['CGE Interneuron', 'LAMP5/LHX6 and \nChandelier Interneuron', 'MGE Interneuron']
dat = [X22q_Bias_CGE["EFFECT"], X22q_Bias_LAMP5["EFFECT"], X22q_Bias_MGE["EFFECT"]]

stat_cge_lamp5, p_cge_lamp5 = mannwhitneyu(dat[0], dat[1])
stat_lamp5_mge, p_lamp5_mge = mannwhitneyu(dat[1], dat[2])
stat_cge_mge, p_cge_mge = mannwhitneyu(dat[0], dat[2])

fig, ax = plt.subplots(figsize=(5, 5), dpi=240, facecolor='none')
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

pos_boxes = ax.boxplot(dat, positions=np.arange(len(labels_4a)), widths=0.5,
                       patch_artist=False,
                       boxprops=dict(color="red", linewidth=1.5),
                       medianprops=dict(color="darkred", linewidth=2),
                       whiskerprops=dict(color="red", linewidth=1.5),
                       capprops=dict(color="red", linewidth=1.5))

for i in range(len(labels_4a)):
    x_pos = np.random.normal(i, 0.04, size=len(dat[i]))
    ax.scatter(x_pos, dat[i], color='red', alpha=0.6, s=20, edgecolor='k', linewidth=0.3)

y_maxes = [max(d) for d in dat]
h = 0.02
y_lamp5_mge = max(y_maxes[1], y_maxes[2]) + 0.05
y_cge_lamp5 = max(y_maxes[0], y_maxes[1]) + 0.06
y_cge_mge = max(y_maxes[0], y_maxes[2]) + 0.10

ax.plot([1, 1, 2, 2], [y_lamp5_mge, y_lamp5_mge+h, y_lamp5_mge+h, y_lamp5_mge], linewidth=1.0, ls='--', alpha=0.7, color='black')
ax.text(1.5, y_lamp5_mge + h + 0.01, format_pval_scientific(p_lamp5_mge), ha='center', va='bottom', fontsize=10)

ax.plot([0, 0, 1, 1], [y_cge_lamp5, y_cge_lamp5+h, y_cge_lamp5+h, y_cge_lamp5], linewidth=1.0, ls='--', alpha=0.7, color='black')
ax.text(0.5, y_cge_lamp5 + h + 0.01, format_pval_scientific(p_cge_lamp5), ha='center', va='bottom', fontsize=10)

ax.plot([0, 0, 2, 2], [y_cge_mge, y_cge_mge+h, y_cge_mge+h, y_cge_mge], linewidth=1.0, ls='--', alpha=0.7, color='black')
ax.text(1, y_cge_mge + h + 0.01, format_pval_scientific(p_cge_mge), ha='center', va='bottom', fontsize=10)

ax.set_xticks(range(len(labels_4a)))
ax.set_xticklabels(labels_4a, fontsize=10, fontweight='normal', rotation=30)
ax.tick_params(axis='y', labelsize=12)
ax.set_ylabel("Cell Type Bias", fontsize=14, fontweight='normal')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('black')
ax.spines['left'].set_alpha(0.9)
ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_alpha(0.9)

plt.tight_layout()
fig.savefig(FIG_DIR + "Fig4_A.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Panel B — 22q11.2: VIP+ vs VIP- CGE interneurons

# %%
from scipy.stats import mannwhitneyu

cutoff = 1.0
vip_pos = CGE_anno_22q[CGE_anno_22q["VIP"] >= cutoff]
vip_neg = CGE_anno_22q[CGE_anno_22q["VIP"] < cutoff]

VIP1 = [vip_pos["EFFECT_22q.11"], vip_pos["EFFECT_22q.11_mouse_gene"]]
VIP2 = [vip_neg["EFFECT_22q.11"], vip_neg["EFFECT_22q.11_mouse_gene"]]
ticks = ['Human 22q11.2\nCNV genes', 'Mouse 16qA13 \nCNV orthologs genes']

fig, ax = plt.subplots(figsize=(6, 5), dpi=240, facecolor='none')
fig.patch.set_alpha(0.0)
ax.patch.set_alpha(0.0)

summer_rain_plot = ax.boxplot(
    VIP1, positions=np.array(np.arange(len(VIP1)))*2.0-0.35, widths=0.5,
    patch_artist=False,
    boxprops=dict(color="red", linewidth=1.5),
    medianprops=dict(color="darkred", linewidth=2),
    whiskerprops=dict(color="red", linewidth=1.5),
    capprops=dict(color="red", linewidth=1.5)
)

winter_rain_plot = ax.boxplot(
    VIP2, positions=np.array(np.arange(len(VIP2)))*2.0+0.35, widths=0.5,
    patch_artist=False,
    boxprops=dict(color="blue", linewidth=1.5),
    medianprops=dict(color="darkblue", linewidth=2),
    whiskerprops=dict(color="blue", linewidth=1.5),
    capprops=dict(color="blue", linewidth=1.5)
)

ax.plot([], c='#D7191C', label='VIP+')
ax.plot([], c='#2C7BB6', label='VIP-')
ax.legend(frameon=False, fontsize=15, loc='center left', bbox_to_anchor=(1.02, 0.8), borderaxespad=0)

for i in range(len(VIP1)):
    x1 = np.random.normal(i*2.0-0.35, 0.04, size=len(VIP1[i]))
    x2 = np.random.normal(i*2.0+0.35, 0.04, size=len(VIP2[i]))
    ax.scatter(x1, VIP1[i], color='#D7191C', alpha=0.6, s=20, edgecolor='k', linewidth=0.3)
    ax.scatter(x2, VIP2[i], color='#2C7BB6', alpha=0.6, s=20, edgecolor='k', linewidth=0.3)

    stat, p = mannwhitneyu(VIP1[i], VIP2[i])
    x1_box, x2_box = i*2.0-0.35, i*2.0+0.35
    y1_max = summer_rain_plot['whiskers'][2*i+1].get_ydata()[1]
    y2_max = winter_rain_plot['whiskers'][2*i+1].get_ydata()[1]
    y = max(y1_max, y2_max) + 0.005
    h = 0.006
    ax.plot([x1_box, x1_box, x2_box, x2_box], [y, y+h, y+h, y], lw=1.0, ls='dashed', alpha=0.8, color='black')
    ax.text((x1_box + x2_box)*.5, y+h, format_pval_scientific(p), ha='center', va='bottom', color='k', fontsize=12, fontweight='normal')

ax.set_xticks(np.arange(0, len(ticks) * 2, 2))
ax.set_xticklabels(ticks, fontsize=12, fontweight='normal')
ax.tick_params(axis='y', labelsize=12)
ax.set_ylabel("CGE Interneuron Bias", fontsize=14, fontweight='normal')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('black')
ax.spines['left'].set_alpha(0.9)
ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_alpha(0.9)

fig.tight_layout()
fig.savefig(FIG_DIR + "Fig4_B.png", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()

# %% [markdown]
# ### Composite Figure 4

# %%
assemble_panels(
    panel_paths=[FIG_DIR + "Fig4_A.png", FIG_DIR + "Fig4_B.png"],
    panel_labels=["A", "B"],
    layout=[(0, slice(0, 1)), (0, slice(1, 2))],
    figsize=(14, 7),
    out_path=FIG_DIR + "Figure4.pdf",
)
