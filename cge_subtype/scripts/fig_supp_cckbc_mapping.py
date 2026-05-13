#!/usr/bin/env python3
"""
Build supplementary figure: cross-species CCKBC mapping evidence.

4-panel figure:
  A. Harmony + scVI CCKBC fraction per human CGE cluster
  B. MetaNeighbor: best mouse subclass match per human CGE cluster (AUROC)
  C. Cross-species ephys feature comparison (Mouse CCKBC vs Mouse VIP-other
     vs Human VIP+ ISI / VIP+ CCKBC / VIP- CCKBC)
  D. 22q11.2 mutation bias by group (VIP-, VIP+ CCKBC, VIP+ ISI)

Output:
  cge_subtype/results/fig_supp_cckbc_mapping.pdf
  cge_subtype/results/fig_supp_cckbc_mapping.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
RESULTS_DIR = PROJECT_DIR / "results"
M1_META = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_metadata.csv")
M1_EPHYS = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_ephys_features.csv")
HUMAN_EPHYS = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")
HUMAN_MAPPING = RESULTS_DIR / "harmony_human_patchseq_validation.csv"

OUT_PDF = RESULTS_DIR / "fig_supp_cckbc_mapping.pdf"
OUT_PNG = RESULTS_DIR / "fig_supp_cckbc_mapping.png"

# ---------------------------------------------------------------------------
# Cluster classification
# ---------------------------------------------------------------------------
HUMAN_VIP_NEG = [277, 278]
HUMAN_VIP_POS_CCKBC = [279, 280, 281]
HUMAN_VIP_POS_ISI = [276, 282, 283, 284, 285, 286, 287, 288, 289, 290,
                     291, 292, 293, 294, 295, 296]

# Mouse CCKBC / VIP-other types (from cross_species_ephys_alignment.py)
MOUSE_CCKBC_TYPES = [
    "Vip Sncg", "Vip Serpinf1_1", "Vip Serpinf1_2",
    "Sncg Npy2r", "Sncg Col14a1", "Sncg Calb1_1", "Sncg Calb1_2",
]
MOUSE_VIP_OTHER_TYPES = [
    "Vip Mybpc1_1", "Vip Mybpc1_2", "Vip Mybpc1_3",
    "Vip Gpc3", "Vip Chat_1", "Vip C1ql1", "Vip Htr1f", "Vip Col15a1",
]

# Colors
COL_HARMONY = "#1f77b4"
COL_SCVI = "#ff7f0e"
COL_SNCG = "#d62728"
COL_VIP = "#9467bd"
COL_LAMP5 = "#2ca02c"
COL_OTHER = "#7f7f7f"
COL_VIP_NEG = "#2ca02c"
COL_VIP_POS_CCKBC = "#d62728"
COL_VIP_POS_ISI = "#1f77b4"

# Plot styling
matplotlib.rcParams["savefig.transparent"] = True
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.size"] = 8


def style_ax(ax):
    ax.patch.set_alpha(0)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(direction="out", length=3)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data...")
summary = pd.read_csv(RESULTS_DIR / "cckbc_convergent_bias_summary.csv")
classification = pd.read_csv(RESULTS_DIR / "updated_22q_bias" / "multimodal_classification.csv")

# Merge cluster_id from summary's "Unnamed: 0" column
summary = summary.rename(columns={"Unnamed: 0": "cluster_id"})
summary["cluster_id"] = summary["cluster_id"].astype(int)
classification["cluster_id"] = classification["cluster_id"].astype(int)

# Reindex to all 21 CGE clusters (276..296)
all_clusters = list(range(276, 297))
ordered = pd.DataFrame({"cluster_id": all_clusters})
summary = ordered.merge(summary, on="cluster_id", how="left")
# Then add AUROC + best_mouse_subclass from classification (always available)
summary = summary.merge(
    classification[["cluster_id", "best_mouse_subclass", "best_auroc"]],
    on="cluster_id",
    how="left",
)

# Fill NaN cckbc fractions with 0
summary["cckbc_frac_harmony"] = summary["cckbc_frac_harmony"].fillna(0)
summary["cckbc_frac_scvi"] = summary["cckbc_frac_scvi"].fillna(0)

# Categorize each cluster
def cluster_category(cid):
    if cid in HUMAN_VIP_NEG:
        return "VIP- (FP)"
    elif cid in HUMAN_VIP_POS_CCKBC:
        return "VIP+ CCKBC"
    else:
        return "VIP+ ISI"
summary["category"] = summary["cluster_id"].apply(cluster_category)

# Best mouse subclass simplification
def simplify_mouse_subclass(s):
    if pd.isna(s):
        return "Other"
    s = str(s)
    if "Sncg" in s:
        return "Sncg"
    elif "Vip" in s:
        return "Vip"
    elif "Lamp5" in s:
        return "Lamp5"
    else:
        return "Other"
summary["mouse_subclass_simple"] = summary["best_mouse_subclass"].apply(simplify_mouse_subclass)

print(summary[["cluster_id", "category", "cckbc_frac_harmony", "cckbc_frac_scvi",
                "mouse_subclass_simple", "best_auroc", "bias_22q_del"]].to_string())


# ---------------------------------------------------------------------------
# Load ephys data
# ---------------------------------------------------------------------------
print("\nLoading ephys data...")
mouse_meta = pd.read_csv(M1_META, sep="\t")
mouse_meta = mouse_meta.set_index("Cell")
mouse_ephys = pd.read_csv(M1_EPHYS).set_index("cell id")

# Mouse features
MOUSE_FEATS = {
    "AP width (ms)": "AP width",
    "ISI CV": "ISI CV",
    "R_input (MΩ)": "R_input",
    "tau (ms)": "tau",
}

mouse_common = mouse_meta.index.intersection(mouse_ephys.index)
mouse_meta = mouse_meta.loc[mouse_common]
mouse_ephys = mouse_ephys.loc[mouse_common]

mouse_groups = pd.Series("Other", index=mouse_common)
mouse_groups[mouse_meta["RNA type"].isin(MOUSE_CCKBC_TYPES)] = "Mouse CCKBC"
mouse_groups[mouse_meta["RNA type"].isin(MOUSE_VIP_OTHER_TYPES)] = "Mouse VIP-other"
print(f"Mouse: CCKBC={(mouse_groups=='Mouse CCKBC').sum()}, "
      f"VIP-other={(mouse_groups=='Mouse VIP-other').sum()}")

# Human ephys
human_ephys = pd.read_csv(HUMAN_EPHYS, index_col=0)
print(f"Human ephys columns: {list(human_ephys.columns)[:20]}")

# Map human cells to clusters via harmony validation
human_map = pd.read_csv(HUMAN_MAPPING, index_col=0)
print(f"Human mapping columns: {list(human_map.columns)}")

# Find cluster column
h_cluster_col = None
for c in human_map.columns:
    if "cluster" in c.lower() and "super" not in c.lower() and "conf" not in c.lower():
        h_cluster_col = c
        break
print(f"Using human cluster column: {h_cluster_col}")

human_common = human_ephys.index.intersection(human_map.index)
human_ephys = human_ephys.loc[human_common]
human_clusters = human_map.loc[human_common, h_cluster_col].astype(float).astype(int)

human_groups = pd.Series("Other", index=human_common)
human_groups[human_clusters.isin(HUMAN_VIP_NEG)] = "Human VIP- CCKBC"
human_groups[human_clusters.isin(HUMAN_VIP_POS_CCKBC)] = "Human VIP+ CCKBC"
human_groups[human_clusters.isin(HUMAN_VIP_POS_ISI)] = "Human VIP+ ISI"
print(f"Human groups: {human_groups.value_counts().to_dict()}")

# Map mouse → human ephys feature names
HUMAN_FEAT_MAP = {
    "AP width (ms)": "width_ramp",
    "ISI CV": "isi_cv_hero",
    "R_input (MΩ)": "input_resistance",
    "tau (ms)": "tau",
}

# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
print("\nBuilding figure...")
fig = plt.figure(figsize=(11, 13), constrained_layout=False)
fig.patch.set_alpha(0)
gs = gridspec.GridSpec(4, 2, figure=fig, height_ratios=[1, 1, 1.4, 1.0],
                       hspace=0.85, wspace=0.32,
                       left=0.08, right=0.97, top=0.96, bottom=0.07)

# Panel A: Harmony + scVI CCKBC fraction
axA = fig.add_subplot(gs[0, :])
style_ax(axA)

x = np.arange(len(summary))
width = 0.38
axA.bar(x - width/2, summary["cckbc_frac_harmony"], width,
        label="Harmony", color=COL_HARMONY, edgecolor="white", linewidth=0.5)
axA.bar(x + width/2, summary["cckbc_frac_scvi"], width,
        label="scVI", color=COL_SCVI, edgecolor="white", linewidth=0.5)

# Highlight category bands
for i, row in summary.iterrows():
    if row["category"] == "VIP+ CCKBC":
        axA.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_POS_CCKBC, zorder=0)
    elif row["category"] == "VIP- (FP)":
        axA.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_NEG, zorder=0)

axA.set_xticks(x)
axA.set_xticklabels(summary["cluster_id"].astype(int), rotation=45, ha="right", fontsize=7)
axA.set_xlabel("Human CGE cluster ID")
axA.set_ylabel("Mouse CCKBC fraction")
axA.set_title("A. Cross-species transcriptomic mapping (M1 patch-seq → human Siletti CGE)",
              loc="left", fontweight="bold", fontsize=10)
axA.set_ylim(0, 1.05)
axA.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
axA.legend(loc="upper right", fontsize=8, frameon=False)

# Panel B: MetaNeighbor best mouse subclass AUROC
axB = fig.add_subplot(gs[1, :])
style_ax(axB)

subclass_colors = {"Sncg": COL_SNCG, "Vip": COL_VIP, "Lamp5": COL_LAMP5, "Other": COL_OTHER}
bar_colors = [subclass_colors[s] for s in summary["mouse_subclass_simple"]]

axB.bar(x, summary["best_auroc"].fillna(0), color=bar_colors,
        edgecolor="white", linewidth=0.5)

# Same category bands
for i, row in summary.iterrows():
    if row["category"] == "VIP+ CCKBC":
        axB.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_POS_CCKBC, zorder=0)
    elif row["category"] == "VIP- (FP)":
        axB.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_NEG, zorder=0)

axB.set_xticks(x)
axB.set_xticklabels(summary["cluster_id"].astype(int), rotation=45, ha="right", fontsize=7)
axB.set_xlabel("Human CGE cluster ID")
axB.set_ylabel("Best mouse cluster\nAUROC")
axB.set_title("B. MetaNeighbor cross-species cluster correspondence (best mouse subclass)",
              loc="left", fontweight="bold", fontsize=10)
axB.set_ylim(0, 1.0)
axB.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

# Legend for subclass colors
from matplotlib.patches import Patch
legend_elems = [Patch(facecolor=c, edgecolor="white", label=l)
                for l, c in subclass_colors.items()]
legend_elems.append(Patch(facecolor=COL_VIP_POS_CCKBC, alpha=0.3,
                          label="Confirmed CCKBC (279-281)"))
legend_elems.append(Patch(facecolor=COL_VIP_NEG, alpha=0.3,
                          label="Harmony FP (277-278)"))
axB.legend(handles=legend_elems, loc="upper right", fontsize=7, frameon=False, ncol=2)

# Panel C: Ephys feature comparison (4 features in a row)
groups_to_plot = ["Mouse\nVIP-other", "Mouse\nCCKBC", "Human\nVIP+ ISI",
                  "Human\nVIP- CCKBC", "Human\nVIP+ CCKBC"]
group_colors = ["#aec7e8", "#1f77b4", "#bdbdbd", "#2ca02c", "#d62728"]

ephys_features_to_plot = [
    ("AP width (ms)", "AP width", "width_ramp", 1000),    # human in seconds -> ms
    ("ISI CV", "ISI CV", "isi_cv_hero", 1),
    ("R_input (MΩ)", "R_input", "input_resistance", 1),
    ("tau (ms)", "tau", "tau", 1000),                      # human in seconds -> ms
]

# Sub-gridspec for Panel C: 1 row x 4 cols spanning both columns of the main grid
gsC = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[2, :], wspace=0.55)

def add_stat(ax, i1, i2, vals1, vals2, y_axes_frac=0.92):
    """Add a stat bracket using axes-fraction y position."""
    if len(vals1) >= 3 and len(vals2) >= 3:
        _, p = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        # Convert axes fraction to data coordinates
        ylim = ax.get_ylim()
        y = ylim[0] + y_axes_frac * (ylim[1] - ylim[0])
        ax.plot([i1+1, i2+1], [y, y], color="black", linewidth=0.8)
        ax.text((i1+i2)/2 + 1, y, sig, ha="center", va="bottom", fontsize=7)

for fi, (feat_name, mouse_feat, human_feat, scale) in enumerate(ephys_features_to_plot):
    axC = fig.add_subplot(gsC[0, fi])
    style_ax(axC)

    # Mouse: feature is in seconds for some, mouse data is already in correct unit
    if "width" in mouse_feat.lower() or "tau" in mouse_feat.lower():
        # Mouse AP width and tau in M1 ephys are in ms already
        m_scale = 1
    else:
        m_scale = 1

    data_groups = []
    data_groups.append(mouse_ephys.loc[mouse_groups == "Mouse VIP-other", mouse_feat].dropna().values * m_scale)
    data_groups.append(mouse_ephys.loc[mouse_groups == "Mouse CCKBC", mouse_feat].dropna().values * m_scale)
    data_groups.append(human_ephys.loc[human_groups == "Human VIP+ ISI", human_feat].dropna().values * scale)
    data_groups.append(human_ephys.loc[human_groups == "Human VIP- CCKBC", human_feat].dropna().values * scale)
    data_groups.append(human_ephys.loc[human_groups == "Human VIP+ CCKBC", human_feat].dropna().values * scale)

    bp = axC.boxplot(data_groups, tick_labels=groups_to_plot, patch_artist=True,
                     widths=0.6, showfliers=False,
                     medianprops=dict(color="black", linewidth=1.0),
                     boxprops=dict(linewidth=0.4),
                     whiskerprops=dict(linewidth=0.4),
                     capprops=dict(linewidth=0.4))
    for patch, c in zip(bp["boxes"], group_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)

    # Strip plot
    for i, vals in enumerate(data_groups):
        if len(vals) > 0:
            jitter = np.random.RandomState(i + fi).uniform(-0.16, 0.16, size=len(vals))
            axC.scatter(np.full(len(vals), i+1) + jitter, vals,
                        s=3, alpha=0.4, color=group_colors[i], edgecolors="none")

    axC.set_ylabel(feat_name, fontsize=8)
    plt.setp(axC.get_xticklabels(), rotation=45, ha="right", fontsize=6)

    # Set ylim and add stats
    all_vals = np.concatenate([v for v in data_groups if len(v) > 0])
    if len(all_vals) > 0:
        v_lo, v_hi = float(np.min(all_vals)), float(np.max(all_vals))
        span = v_hi - v_lo if v_hi > v_lo else max(abs(v_hi), 1.0) * 0.05
        axC.set_ylim(v_lo - 0.1 * span, v_hi + 0.18 * span)

    add_stat(axC, 0, 1, data_groups[0], data_groups[1], y_axes_frac=0.85)
    add_stat(axC, 2, 4, data_groups[2], data_groups[4], y_axes_frac=0.95)

    if fi == 0:
        axC.set_title("C. Cross-species electrophysiology features",
                      loc="left", fontweight="bold", fontsize=10, pad=12)

# Panel D: 22q bias by group
axD = fig.add_subplot(gs[3, :])
style_ax(axD)

bias_col = "bias_22q_del"
group_data = {
    "VIP-\n(277, 278)": summary[summary["category"] == "VIP- (FP)"][bias_col].dropna().values,
    "VIP+ CCKBC\n(279-281)": summary[summary["category"] == "VIP+ CCKBC"][bias_col].dropna().values,
    "VIP+ ISI\n(others)": summary[summary["category"] == "VIP+ ISI"][bias_col].dropna().values,
}

box_colors_d = [COL_VIP_NEG, COL_VIP_POS_CCKBC, COL_VIP_POS_ISI]
positions = [1, 2, 3]
bp = axD.boxplot(list(group_data.values()), positions=positions,
                 tick_labels=list(group_data.keys()), patch_artist=True, widths=0.55,
                 showfliers=False,
                 medianprops=dict(color="black", linewidth=1.2),
                 boxprops=dict(linewidth=0.5),
                 whiskerprops=dict(linewidth=0.5),
                 capprops=dict(linewidth=0.5))
for patch, c in zip(bp["boxes"], box_colors_d):
    patch.set_facecolor(c)
    patch.set_alpha(0.6)

# Individual cluster points
for i, (name, vals) in enumerate(group_data.items()):
    jitter = np.random.RandomState(i).uniform(-0.12, 0.12, size=len(vals))
    axD.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                s=20, color=box_colors_d[i], edgecolors="white", linewidth=0.5,
                zorder=3)

axD.set_ylabel("22q11.2 deletion bias (EFFECT)")
axD.set_title("D. 22q11.2 mutation bias by group", loc="left",
              fontweight="bold", fontsize=10)
axD.tick_params(axis="x", labelsize=7)

# Stat tests
vals_neg = group_data["VIP-\n(277, 278)"]
vals_cckbc = group_data["VIP+ CCKBC\n(279-281)"]
vals_isi = group_data["VIP+ ISI\n(others)"]

# CCKBC vs ISI
if len(vals_cckbc) >= 2 and len(vals_isi) >= 2:
    _, p_cckbc_isi = stats.mannwhitneyu(vals_cckbc, vals_isi, alternative="two-sided")
else:
    p_cckbc_isi = np.nan

# VIP- vs ISI
if len(vals_neg) >= 2 and len(vals_isi) >= 2:
    _, p_neg_isi = stats.mannwhitneyu(vals_neg, vals_isi, alternative="two-sided")
else:
    p_neg_isi = np.nan

ymax = max(np.max(np.concatenate([vals_neg, vals_cckbc, vals_isi])),
           0.2)
axD.set_ylim(-0.02, ymax * 1.35)

# CCKBC vs ISI bracket
y1 = ymax * 1.05
axD.plot([2, 3], [y1, y1], color="black", linewidth=0.8)
sig_text = (f"P={p_cckbc_isi:.2f}" if not np.isnan(p_cckbc_isi)
            else "N/A")
axD.text(2.5, y1, f" n.s. ({sig_text})", ha="center", va="bottom", fontsize=7)

# VIP- vs ISI bracket
y2 = ymax * 1.20
axD.plot([1, 3], [y2, y2], color="black", linewidth=0.8)
sig_text2 = (f"P={p_neg_isi:.3f}" if not np.isnan(p_neg_isi)
             else "N/A")
star = "*" if not np.isnan(p_neg_isi) and p_neg_isi < 0.05 else ""
axD.text(2, y2, f" {star}{sig_text2}", ha="center", va="bottom", fontsize=7)

# Save
plt.savefig(OUT_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"\nSaved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")
