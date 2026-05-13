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
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Supplementary Figure: Cross-Species CCKBC Mapping
#
# **4-panel figure:**
# - **A.** Harmony + scVI CCKBC fraction per human CGE cluster
# - **B.** MetaNeighbor: best mouse subclass match per cluster (AUROC)
# - **C.** Cross-species electrophysiology comparison (4 features)
# - **D.** 22q11.2 mutation bias by group (VIP-, VIP+ CCKBC, VIP+ ISI)
#
# **Key result:** VIP+ CCKBC clusters (279-281) show 22q bias comparable to
# VIP+ ISI clusters (P=0.29, n.s.). The primary split is VIP- vs VIP+
# (P=0.026), not CCKBC vs ISI identity.

# %%
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy import stats

# %% [markdown]
# ## Paths and constants

# %%
PROJECT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
RESULTS_DIR = PROJECT_DIR / "results"
M1_META = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_metadata.csv")
# Use the PUBLISHED Scala et al. 2021 ephys features from the mini-atlas repo.
# This is the canonical feature set from the paper (1328 cells, 29 features with
# proper names+units), NOT a re-extraction. Column names follow the published
# convention (e.g., "AP width (ms)", "Input resistance (MOhm)").
M1_EPHYS = Path("/home/jw3514/Work/NeurSim/mini-atlas/data/m1_patchseq_ephys_features.csv")
HUMAN_EPHYS = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")
HUMAN_MAPPING = RESULTS_DIR / "harmony_human_patchseq_validation.csv"

OUT_PDF = RESULTS_DIR / "fig_supp_cckbc_mapping.pdf"
OUT_PNG = RESULTS_DIR / "fig_supp_cckbc_mapping.png"

# Cluster groupings (refined per marker-gene evidence in Supp Note 2.3)
#   277: False positive — low CCK/CNR1, high RELN/M2R (neurogliaform-like)
#   278: VIP- CCKBC — very high CCK (rank 4/461), high CNR1/CXCL14/CRH (marker-confirmed)
#   279-281: VIP+ CCKBC — transcriptomic + MetaNeighbor Sncg homologs
HUMAN_FP = [277]                         # Harmony false positive, NOT a CCKBC
HUMAN_VIP_NEG_CCKBC = [278]              # VIP-, marker-confirmed CCKBC
HUMAN_VIP_POS_CCKBC = [279, 280, 281]    # VIP+, Sncg-confirmed CCKBC
HUMAN_VIP_POS_ISI = [276, 282, 283, 284, 285, 286, 287, 288, 289, 290,
                     291, 292, 293, 294, 295, 296]
# Backward compat: VIP- includes 277+278 for some analyses
HUMAN_VIP_NEG = HUMAN_FP + HUMAN_VIP_NEG_CCKBC

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
COL_FP = "#7f7f7f"             # gray - false positive (277)
COL_VIP_NEG_CCKBC = "#ff7f0e"  # orange - VIP- CCKBC (278, marker-confirmed)
COL_VIP_POS_CCKBC = "#d62728"  # red - VIP+ CCKBC (279-281)
COL_VIP_POS_ISI = "#1f77b4"    # blue - VIP+ ISI
COL_VIP_NEG = COL_VIP_NEG_CCKBC  # back-compat

matplotlib.rcParams["savefig.transparent"] = True
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.size"] = 8


def style_ax(ax):
    ax.patch.set_alpha(0)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(direction="out", length=3)


# %% [markdown]
# ## Load and prepare cluster summary
#
# Compute per-cluster CCKBC fraction in three flavors:
# - **All** (V1+M1, CGE families): the original analysis
# - **M1 only** (M1, CGE families): cleaner but smaller sample
# - **V1 only** (V1, CGE families): larger sample but only V1 dataset

# %%
classification = pd.read_csv(RESULTS_DIR / "updated_22q_bias" / "multimodal_classification.csv")
classification["cluster_id"] = classification["cluster_id"].astype(int)

# Per-cell mapping data
harmony_map = pd.read_csv(RESULTS_DIR / "harmony_patchseq_mapping_results.csv", index_col=0)
scvi_map = pd.read_csv(RESULTS_DIR / "patchseq_mapping_results.csv", index_col=0)
# Harmony CSV doesn't have dataset; pull from scVI mapping (same cell IDs)
harmony_map["dataset"] = scvi_map.loc[harmony_map.index, "dataset"].values
harmony_map["mouse_RNA family"] = scvi_map.loc[harmony_map.index, "mouse_RNA family"].values

CGE_FAMILIES = ["Sncg", "Vip", "Lamp5", "Serpinf1"]  # Serpinf1 only appears in V1 alias names


def get_cge_mask(df):
    """CGE membership: M1 uses RNA family; V1 uses first word of AIT alias."""
    # M1: RNA family in CGE families
    fam = df["mouse_RNA family"].astype(str)
    m1_cge = fam.isin(CGE_FAMILIES)
    # V1: alias starts with one of the CGE family prefixes
    alias = df["mouse_corresponding_AIT2.3.1_alias"].astype(str)
    v1_cge = pd.Series(False, index=df.index)
    for fam_prefix in CGE_FAMILIES:
        v1_cge |= alias.str.startswith(fam_prefix + " ", na=False)
    return m1_cge | v1_cge


def filter_subset(df, dataset_filter, cluster_col):
    """Filter to dataset (or 'all') + CGE cells."""
    cge_mask = get_cge_mask(df)
    if dataset_filter == "all":
        mask = cge_mask
    else:
        mask = (df["dataset"] == dataset_filter) & cge_mask
    return df[mask]


def cckbc_frac_per_cluster(df, cluster_col):
    out = {}
    for cid, grp in df.groupby(cluster_col):
        try:
            cid_int = int(float(cid))
        except (ValueError, TypeError):
            continue
        out[cid_int] = (grp["is_cckbc"].sum() / len(grp), len(grp), int(grp["is_cckbc"].sum()))
    return out


def per_cluster_table(harmony_df, scvi_df):
    h_frac = cckbc_frac_per_cluster(harmony_df, "assigned_cluster")
    s_frac = cckbc_frac_per_cluster(scvi_df, "assigned_human_cluster")
    rows = []
    for cid in range(276, 297):
        h = h_frac.get(cid, (0, 0, 0))
        s = s_frac.get(cid, (0, 0, 0))
        rows.append({
            "cluster_id": cid,
            "cckbc_frac_harmony": h[0], "n_harmony": h[1], "n_cckbc_harmony": h[2],
            "cckbc_frac_scvi": s[0], "n_scvi": s[1], "n_cckbc_scvi": s[2],
        })
    return pd.DataFrame(rows)


# Three subsets
h_all = filter_subset(harmony_map, "all", "assigned_cluster")
s_all = filter_subset(scvi_map, "all", "assigned_human_cluster")
h_m1 = filter_subset(harmony_map, "M1", "assigned_cluster")
s_m1 = filter_subset(scvi_map, "M1", "assigned_human_cluster")
h_v1 = filter_subset(harmony_map, "V1", "assigned_cluster")
s_v1 = filter_subset(scvi_map, "V1", "assigned_human_cluster")

print(f"Harmony — All CGE: {len(h_all)} cells ({h_all['is_cckbc'].sum()} CCKBC)")
print(f"Harmony — M1 CGE:  {len(h_m1)} cells ({h_m1['is_cckbc'].sum()} CCKBC)")
print(f"Harmony — V1 CGE:  {len(h_v1)} cells ({h_v1['is_cckbc'].sum()} CCKBC)")

summary_all = per_cluster_table(h_all, s_all)
summary_m1 = per_cluster_table(h_m1, s_m1)
summary_v1 = per_cluster_table(h_v1, s_v1)

# Add bias + MetaNeighbor (cluster-level, dataset-independent) for downstream panels
old_summary = pd.read_csv(RESULTS_DIR / "cckbc_convergent_bias_summary.csv")
old_summary = old_summary.rename(columns={"Unnamed: 0": "cluster_id"})
old_summary["cluster_id"] = old_summary["cluster_id"].astype(int)


def annotate(summary):
    summary = summary.merge(
        old_summary[["cluster_id", "bias_22q_del", "bias_22q_mouse"]],
        on="cluster_id", how="left",
    )
    summary = summary.merge(
        classification[["cluster_id", "best_mouse_subclass", "best_auroc"]],
        on="cluster_id", how="left",
    )
    return summary


summary_all = annotate(summary_all)
summary_m1 = annotate(summary_m1)
summary_v1 = annotate(summary_v1)

# Use M1 as the "main" summary for downstream (Panel B/D use cluster-level data only)
summary = summary_m1


def cluster_category(cid):
    if cid in HUMAN_FP:
        return "FP (277)"
    elif cid in HUMAN_VIP_NEG_CCKBC:
        return "VIP- CCKBC (278)"
    elif cid in HUMAN_VIP_POS_CCKBC:
        return "VIP+ CCKBC (279-281)"
    else:
        return "VIP+ ISI"


def simplify_mouse_subclass(s):
    """Extract the Allen WMB subclass id-label (e.g. '047 Sncg Gaba') from
    the multimodal_classification's full string '<subclass> (<cluster>)'."""
    if pd.isna(s):
        return "Unknown"
    s = str(s).strip()
    # Format: "047 Sncg Gaba (0681 Sncg Gaba_5)" → take the part before " ("
    if " (" in s:
        return s.split(" (", 1)[0]
    return s


for s in [summary_all, summary_m1, summary_v1]:
    s["category"] = s["cluster_id"].apply(cluster_category)
    s["mouse_subclass_simple"] = s["best_mouse_subclass"].apply(simplify_mouse_subclass)

summary = summary_m1
summary[["cluster_id", "category", "cckbc_frac_harmony", "cckbc_frac_scvi",
         "mouse_subclass_simple", "best_auroc", "bias_22q_del"]]

# %% [markdown]
# ## Load ephys data

# %%
mouse_meta = pd.read_csv(M1_META, sep="\t").set_index("Cell")
mouse_ephys = pd.read_csv(M1_EPHYS, index_col=0)

mouse_common = mouse_meta.index.intersection(mouse_ephys.index)
mouse_meta = mouse_meta.loc[mouse_common]
mouse_ephys = mouse_ephys.loc[mouse_common]

mouse_groups = pd.Series("Other", index=mouse_common)
mouse_groups[mouse_meta["RNA type"].isin(MOUSE_CCKBC_TYPES)] = "Mouse CCKBC"
mouse_groups[mouse_meta["RNA type"].isin(MOUSE_VIP_OTHER_TYPES)] = "Mouse VIP-other"

human_ephys = pd.read_csv(HUMAN_EPHYS, index_col=0)

# Human cell classification uses the GROUND-TRUTH Lee/Dalley transcriptomic
# labels, NOT the Harmony-mapped Siletti cluster. This is the honest
# classification: "Human CCKBC" = cells whose published Lee/Dalley
# Transcriptomic_type matches 'Inh L2-5 VIP SERPINF1' (the only CCKBC-like
# type in the human MTG patch-seq dataset, n=3 — small but unambiguous).
# All other VIP cells are labeled "Human VIP-other".
LEEDALLEY_META = Path(
    "/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_manuscript_metadata.csv"
)
HUMAN_CCKBC_TRANSCRIPTOMIC_TYPES = ["Inh L2-5 VIP SERPINF1"]

leedalley = pd.read_csv(LEEDALLEY_META)
# The transcriptomic_type column has a trailing space in the header
_ttype_col = [c for c in leedalley.columns if "Transcriptomic_type" in c][0]
leedalley = leedalley.set_index("specimen_id")

# Intersect with ephys cells
human_common = human_ephys.index.intersection(leedalley.index)
human_ephys = human_ephys.loc[human_common]
lee_sub = leedalley.loc[human_common]

# STRICT filter: cell must be both VIP subclass AND have "VIP" in its
# Transcriptomic_type name. This excludes 8 cells that are classified as
# VIP subclass but have non-VIP transcriptomic types (4 ADARB2 MC4R, 3 SST
# NMBR, 1 SST BAGE2) — leaving a cleaner "VIP-other" pool (n=100 with ephys).
is_vip_subclass = lee_sub["Revised_subclass_label"] == "VIP"
ttype = lee_sub[_ttype_col].astype(str)
is_vip_typed = ttype.str.contains(" VIP ", case=False, regex=False)
is_strict_vip = is_vip_subclass & is_vip_typed

human_groups = pd.Series("Other", index=human_common)
human_groups[is_strict_vip & ttype.isin(HUMAN_CCKBC_TRANSCRIPTOMIC_TYPES)] = "Human CCKBC"
human_groups[is_strict_vip & ~ttype.isin(HUMAN_CCKBC_TRANSCRIPTOMIC_TYPES)] = "Human VIP-other"

print(f"Mouse: CCKBC={(mouse_groups=='Mouse CCKBC').sum()}, "
      f"VIP-other={(mouse_groups=='Mouse VIP-other').sum()}")
print(f"Human (ground-truth Lee/Dalley labels): "
      f"{human_groups.value_counts().to_dict()}")

# %% [markdown]
# ## Panel A: Cross-species transcriptomic mapping — split by dataset
#
# Two sub-panels showing CCKBC fraction per human CGE cluster:
# - Top: M1 only (51 CCKBC, 257 cells)
# - Bottom: V1 only (282 CCKBC, ~5 000 cells)
#
# Splitting reveals whether results are consistent across datasets.

# %%
def plot_panel_A(ax, summary_df, title, n_cells, n_cckbc):
    style_ax(ax)
    x = np.arange(len(summary_df))
    width = 0.38
    ax.bar(x - width/2, summary_df["cckbc_frac_harmony"], width,
           label="Harmony", color=COL_HARMONY, edgecolor="white", linewidth=0.5)
    ax.bar(x + width/2, summary_df["cckbc_frac_scvi"], width,
           label="scVI", color=COL_SCVI, edgecolor="white", linewidth=0.5)
    for i, row in summary_df.iterrows():
        if row["category"] == "VIP+ CCKBC (279-281)":
            ax.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_POS_CCKBC, zorder=0)
        elif row["category"] == "VIP- CCKBC (278)":
            ax.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_VIP_NEG_CCKBC, zorder=0)
        elif row["category"] == "FP (277)":
            ax.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_FP, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df["cluster_id"].astype(int), rotation=45, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_xlabel("Human CGE cluster ID", fontsize=12)
    ax.set_ylabel("Mouse CCKBC fraction", fontsize=12)
    ax.set_title(f"{title}  (n={n_cells} cells, {n_cckbc} CCKBC)",
                 loc="left", fontweight="bold", fontsize=13)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.legend(loc="upper right", fontsize=11, frameon=False)


fig_A = plt.figure(figsize=(11, 6.5))
fig_A.patch.set_alpha(0)
gs_A = gridspec.GridSpec(2, 1, figure=fig_A, hspace=0.55,
                          left=0.08, right=0.97, top=0.96, bottom=0.08)

ax_A1 = fig_A.add_subplot(gs_A[0, 0])
plot_panel_A(ax_A1, summary_m1,
             "M1 only (Scala et al. 2021)",
             len(h_m1), int(h_m1["is_cckbc"].sum()))

ax_A2 = fig_A.add_subplot(gs_A[1, 0])
plot_panel_A(ax_A2, summary_v1,
             "V1 only (Gouwens et al. 2020)",
             len(h_v1), int(h_v1["is_cckbc"].sum()))

OUT_A_PDF = RESULTS_DIR / "fig_supp_panelA_mapping.pdf"
OUT_A_PNG = RESULTS_DIR / "fig_supp_panelA_mapping.png"
plt.savefig(OUT_A_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_A_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_A_PDF}")
plt.show()

# %% [markdown]
# ## Panel B — Within-mouse validation UMAP: Sncg subclass = CCKBC
#
# Mouse M1 patch-seq cells projected onto the WMB-10Xv3 cortical GABAergic
# atlas via Harmony + k-NN. The atlas is shown as transparent gray dots
# with each subclass outlined by a smoothed 2D-KDE contour (one contour
# per DBSCAN-detected island, so disconnected components stay separate).
# M1 patch-seq query cells are overlaid as colored dots, with Sncg/CCKBC
# cells highlighted in amber. CCKBC dots fall directly on the WMB Sncg
# Gaba region (red contour), independently confirming that the operational
# CCKBC definition (Sncg subclass + Vip Sncg + Vip Serpinf1 supertypes;
# Gouwens 2020) is biologically correct.

# %%
import scanpy as sc  # noqa: E402
from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.patheffects as patheffects
from sklearn.cluster import DBSCAN

# Load Sncg validation UMAP h5ad (within-mouse mapping)
val_umap_h5ad = RESULTS_DIR / "sncg_validation_umap.h5ad"
adata_val_umap = sc.read_h5ad(val_umap_h5ad)
umap_v = adata_val_umap.obsm["X_umap"]
obs_v = adata_val_umap.obs.copy()
ref_mask_v = (obs_v["source"] == "atlas").values
query_mask_v = (obs_v["source"] == "patchseq").values

# Subclass color scheme — chosen to be maximally distinct on a single panel.
# Vip is hot magenta (was purple, too close to Pvalb blue), Sst is dark
# brown (was gray, washed out).
SUBCLASS_COLORS_VAL = {
    "047 Sncg Gaba":            "#d62728",  # red
    "046 Vip Gaba":             "#e91e63",  # magenta/hot pink
    "049 Lamp5 Gaba":           "#2ca02c",  # green
    "050 Lamp5 Lhx6 Gaba":      "#bcbd22",  # olive
    "052 Pvalb Gaba":           "#1f77b4",  # blue
    "053 Sst Gaba":             "#8c564b",  # dark brown
    "048 RHP-COA Ndnf Gaba":    "#17becf",  # cyan
    "051 Pvalb chandelier Gaba": "#9467bd",  # purple
}

# Query (M1 patch-seq) coloring matches contour colors
QUERY_COLORS = {
    "Sncg":   "#d62728",
    "Vip":    "#e91e63",
    "Lamp5":  "#2ca02c",
    "Pvalb":  "#1f77b4",
    "Sst":    "#8c564b",
}

fig_single = plt.figure(figsize=(8.5, 6.0))
fig_single.patch.set_alpha(0)
ax_s = fig_single.add_subplot(111)
style_ax(ax_s)

ref_subclass = obs_v.loc[ref_mask_v, "subclass"].astype(str).values
ref_coords = umap_v[ref_mask_v]

# Order: draw subclasses in a stable order so legend order is consistent
SUBCLASS_DRAW_ORDER = [
    "047 Sncg Gaba",
    "046 Vip Gaba",
    "049 Lamp5 Gaba",
    "050 Lamp5 Lhx6 Gaba",
    "052 Pvalb Gaba",
    "053 Sst Gaba",
    "048 RHP-COA Ndnf Gaba",
    "051 Pvalb chandelier Gaba",
]

# Filter out very small subclasses that are not informative for the CCKBC
# question (e.g., 048 RHP-COA Ndnf has only ~60 cells).
SUBCLASS_MIN_CELLS = 200
shown_subclasses = [
    s for s in SUBCLASS_DRAW_ORDER
    if (ref_subclass == s).sum() >= SUBCLASS_MIN_CELLS
]
shown_mask = np.isin(ref_subclass, shown_subclasses)

# DBSCAN params (used twice: once for the main-island axis-limit pass,
# again for per-island contours below)
DBSCAN_EPS = 0.8
DBSCAN_MIN = 10

# Compute axis limits using only the LARGEST DBSCAN island per subclass —
# this drops the small "outlier islands" that pull the bounding box wide
# (e.g., a few stray Sst cells in the Pvalb region).
def main_island_coords(coords_arr):
    """Return coordinates of the largest DBSCAN-detected island."""
    if len(coords_arr) < DBSCAN_MIN:
        return coords_arr
    db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN).fit(coords_arr)
    labels = db.labels_
    valid = labels[labels >= 0]
    if len(valid) == 0:
        return coords_arr
    sizes = {cid: (labels == cid).sum() for cid in np.unique(valid)}
    biggest = max(sizes, key=sizes.get)
    return coords_arr[labels == biggest]


relevant_xs = []
relevant_ys = []
for s in shown_subclasses:
    sub_coords = ref_coords[ref_subclass == s]
    main = main_island_coords(sub_coords)
    relevant_xs.append(main[:, 0])
    relevant_ys.append(main[:, 1])
# Add query cells too — they should always be inside the visible region
relevant_xs.append(umap_v[query_mask_v, 0])
relevant_ys.append(umap_v[query_mask_v, 1])
relevant_x = np.concatenate(relevant_xs)
relevant_y = np.concatenate(relevant_ys)
xmin = float(relevant_x.min()) - 1.0
xmax = float(relevant_x.max()) + 1.0
ymin = float(relevant_y.min()) - 1.0
ymax = float(relevant_y.max()) + 1.0

# 1) Atlas cells as transparent gray background — only the shown subclasses,
#    AND only those falling inside the cropped axis range
ref_coords_shown = ref_coords[shown_mask]
in_view = (
    (ref_coords_shown[:, 0] >= xmin) & (ref_coords_shown[:, 0] <= xmax)
    & (ref_coords_shown[:, 1] >= ymin) & (ref_coords_shown[:, 1] <= ymax)
)
ax_s.scatter(ref_coords_shown[in_view, 0], ref_coords_shown[in_view, 1],
             s=2.5, c="#bdbdbd", alpha=0.20, rasterized=True,
             edgecolors="none", zorder=1)

# 2) Smoothed 2D KDE contour outline per subclass — only shown subclasses
xx, yy = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
grid_pos = np.vstack([xx.ravel(), yy.ravel()])

# Track centroid for each drawn subclass (for label placement) — use the
# centroid of the LARGEST DBSCAN island so labels never fall in empty space
# between disconnected components.
subclass_centroids = {}

ISLAND_MIN_POINTS = 15  # don't draw a contour for islands smaller than this

for sub_label in shown_subclasses:
    mask = ref_subclass == sub_label
    n = mask.sum()
    if n < 30:
        continue
    color = SUBCLASS_COLORS_VAL.get(sub_label, "#cccccc")
    coords = ref_coords[mask]

    # 1) Find UMAP islands via DBSCAN (so we draw separate contours for
    #    physically disconnected components of the same subclass)
    db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN).fit(coords)
    island_labels = db.labels_

    largest_island_size = 0
    largest_island_centroid = (coords[:, 0].mean(), coords[:, 1].mean())

    # 2) For each island, fit its own KDE and draw a contour
    for cluster_id in np.unique(island_labels[island_labels >= 0]):
        island_coords = coords[island_labels == cluster_id]
        if len(island_coords) < ISLAND_MIN_POINTS:
            continue

        try:
            kde = gaussian_kde(island_coords.T, bw_method=0.30)
            z = kde(grid_pos).reshape(xx.shape)
        except Exception:
            continue

        # Contour level capturing ~80% of this island's cell density
        z_at_cells = kde(island_coords.T)
        level = float(np.quantile(z_at_cells, 0.20))
        z_top = float(z.max())
        # Skip degenerate islands where the KDE level is too close to the
        # max (happens for very small/dense islands)
        if z_top - level < 1e-9 or not np.isfinite(level) or not np.isfinite(z_top):
            continue

        # Filled (very transparent) + outlined contour
        try:
            ax_s.contourf(xx, yy, z, levels=[level, z_top * 1.01 + 1e-9],
                          colors=[color], alpha=0.10, zorder=2)
            ax_s.contour(xx, yy, z, levels=[level],
                         colors=[color], linewidths=2.0, zorder=3)
        except ValueError:
            continue

        # Track largest island for label placement
        if len(island_coords) > largest_island_size:
            largest_island_size = len(island_coords)
            largest_island_centroid = (island_coords[:, 0].mean(),
                                       island_coords[:, 1].mean())

    if largest_island_size > 0:
        subclass_centroids[sub_label] = largest_island_centroid

# Add subclass text labels at centroids (small, outlined for legibility)
for sub_label, (cx, cy) in subclass_centroids.items():
    color = SUBCLASS_COLORS_VAL.get(sub_label, "#cccccc")
    short_label = sub_label.replace(" Gaba", "")  # "047 Sncg Gaba" → "047 Sncg"
    txt = ax_s.text(cx, cy, short_label, ha="center", va="center",
                    fontsize=8, fontweight="bold", color=color, zorder=10)
    txt.set_path_effects([
        patheffects.Stroke(linewidth=2, foreground="white"),
        patheffects.Normal(),
    ])

# 3) M1 patch-seq query cells overlaid (CCKBCs highlighted in amber)
q_family = obs_v.loc[query_mask_v, "RNA_family"].astype(str).values
q_iscckbc = obs_v.loc[query_mask_v, "is_cckbc"].astype(str).values
q_coords = umap_v[query_mask_v]

# Non-CCKBC by RNA family
q_legend_handles = []
non_cckbc_idx = q_iscckbc != "1.0"
for fam in ["Sst", "Pvalb", "Lamp5", "Vip"]:
    sel = non_cckbc_idx & (q_family == fam)
    if sel.sum() == 0:
        continue
    color = QUERY_COLORS[fam]
    ax_s.scatter(q_coords[sel, 0], q_coords[sel, 1],
                 s=14, c=color, alpha=0.85,
                 edgecolors="white", linewidth=0.5, zorder=4)
    q_legend_handles.append(
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color,
               markeredgecolor="white", markersize=6,
               label=f"M1 {fam} (n={int(sel.sum())})", linestyle="none")
    )

# CCKBC overlay (highlighted)
cckbc_idx = q_iscckbc == "1.0"
ax_s.scatter(q_coords[cckbc_idx, 0], q_coords[cckbc_idx, 1],
             s=40, c="#fbb03b", alpha=1.0,
             edgecolors="black", linewidth=0.8, zorder=5)
q_legend_handles.append(
    Line2D([0], [0], marker="o", color="black", markerfacecolor="#fbb03b",
           markeredgewidth=0.8, markersize=8,
           label=f"M1 CCKBC (n={int(cckbc_idx.sum())})", linestyle="none")
)

ax_s.set_xlabel("UMAP 1", fontsize=13)
ax_s.set_ylabel("UMAP 2", fontsize=13)
ax_s.tick_params(axis="both", labelsize=11)
ax_s.set_xlim(xmin, xmax)
ax_s.set_ylim(ymin, ymax)

# Place legend outside the plot to the right so it never overlaps cells
ax_s.legend(handles=q_legend_handles,
            loc="center left", bbox_to_anchor=(1.02, 0.5),
            fontsize=9, frameon=False, title="M1 patch-seq",
            title_fontsize=9)

OUT_SINGLE_PDF = RESULTS_DIR / "fig_supp_sncg_validation_umap_single.pdf"
OUT_SINGLE_PNG = RESULTS_DIR / "fig_supp_sncg_validation_umap_single.png"
plt.savefig(OUT_SINGLE_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_SINGLE_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_SINGLE_PDF}")
plt.show()

# %% [markdown]
# ## Panel C: MetaNeighbor cross-species cluster correspondence
#
# Standalone figure: best mouse subclass match per human cluster (AUROC).
# Dataset-independent (uses cluster centroids, not individual cells).

# %%
fig_B = plt.figure(figsize=(11.5, 4.8))
fig_B.patch.set_alpha(0)
ax_B = fig_B.add_subplot(111)
style_ax(ax_B)

# Use M1 summary just for category/cluster ordering — these are dataset-independent
summary_for_B = summary_m1.copy()
x = np.arange(len(summary_for_B))

# Map exact Allen WMB subclass labels to colors. Match by substring on the
# subclass family so we still get a consistent color even if a different
# Sncg/Vip/Lamp5/Sst variant appears.
SUBCLASS_FAMILY_COLORS = [
    ("Sncg", COL_SNCG),
    ("Vip", COL_VIP),
    ("Lamp5", COL_LAMP5),
    ("Sst", "#7f7f7f"),  # gray (was COL_OTHER) — only one cluster (276) hits Sst
]


def subclass_to_color(label):
    for family, color in SUBCLASS_FAMILY_COLORS:
        if family in label:
            return color
    return "#cccccc"


bar_colors = [subclass_to_color(s) for s in summary_for_B["mouse_subclass_simple"]]

ax_B.bar(x, summary_for_B["best_auroc"].fillna(0), color=bar_colors,
         edgecolor="white", linewidth=0.5, zorder=2)

# Single highlight band over the Harmony/scVI-identified CCKBC candidates
# (clusters 277-281). The bar colors below (Sncg/Vip/Lamp5/Other) tell the
# real story about which of these are genuine Sncg homologs.
HARMONY_SCVI_CANDIDATES = HUMAN_FP + HUMAN_VIP_NEG_CCKBC + HUMAN_VIP_POS_CCKBC
COL_CANDIDATE_BAND = "#fbb03b"  # neutral amber
for i, row in summary_for_B.iterrows():
    if int(row["cluster_id"]) in HARMONY_SCVI_CANDIDATES:
        ax_B.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=COL_CANDIDATE_BAND, zorder=0)

ax_B.set_xticks(x)
ax_B.set_xticklabels(summary_for_B["cluster_id"].astype(int), rotation=45, ha="right", fontsize=11)
ax_B.tick_params(axis="y", labelsize=11)
ax_B.set_xlabel("Human CGE cluster ID", fontsize=12)
ax_B.set_ylabel("Best mouse cluster AUROC", fontsize=12)
ax_B.set_title(
    "MetaNeighbor cross-species cluster correspondence (best mouse subclass)",
    loc="left", fontweight="bold", fontsize=13,
)
ax_B.set_ylim(0, 1.0)
ax_B.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, zorder=1)

# Build legend dynamically from subclasses actually present in the bars,
# preserving the order Sncg → Vip → Lamp5 → Sst → other
present_subclasses = []
seen = set()
order_keys = ["Sncg", "Vip", "Lamp5", "Sst"]
for fam in order_keys:
    for s in summary_for_B["mouse_subclass_simple"]:
        if fam in s and s not in seen:
            present_subclasses.append(s)
            seen.add(s)
# Catch anything that doesn't match the canonical families
for s in summary_for_B["mouse_subclass_simple"]:
    if s not in seen:
        present_subclasses.append(s)
        seen.add(s)

legend_elems = [
    Patch(facecolor=subclass_to_color(s), edgecolor="white", label=s)
    for s in present_subclasses
]
legend_elems.append(Patch(facecolor=COL_CANDIDATE_BAND, alpha=0.3,
                          label="Harmony/scVI candidates (277–281)"))
fig_B.subplots_adjust(left=0.10, right=0.74, bottom=0.18, top=0.94)
fig_B.legend(
    handles=legend_elems,
    loc="upper right",
    bbox_to_anchor=(1.0, 1.0),
    bbox_transform=fig_B.transFigure,
    fontsize=10,
    frameon=False,
    ncol=2,
    columnspacing=1.0,
    labelspacing=0.75,
    title="Bar color: mouse subclass",
    title_fontsize=11,
)

OUT_B_PDF = RESULTS_DIR / "fig_supp_panelB_metaneighbor.pdf"
OUT_B_PNG = RESULTS_DIR / "fig_supp_panelB_metaneighbor.png"
plt.savefig(OUT_B_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_B_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_B_PDF}")
plt.show()

# %% [markdown]
# ## Panel D — CCKBC marker gene expression heatmap
#
# Pseudobulk TPM expression of canonical CCKBC markers (CCK, CNR1/CB1R,
# CXCL14, CRH, SNCG) and discriminating markers (VIP, CHRM2/M2R, RELN)
# across all 21 human Siletti CGE clusters. The 5 Harmony/scVI candidate
# CCKBC clusters (277-281) are highlighted with an amber band, matching
# the convention used in Panel C (MetaNeighbor). The bar values let the
# reader judge each candidate independently:
# - **279, 280, 281**: high CCK + CNR1, VIP+ → confirmed VIP+ CCKBC
# - **278**: very high CCK (rank 4/461), high CNR1/CXCL14/CRH, VIP-negative → marker-confirmed VIP- CCKBC
# - **277**: low CCK, low CNR1, very high CHRM2/M2R, very high RELN → neurogliaform-like, NOT a CCKBC

# %%
TPM_FILE = PROJECT_DIR.parent / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.csv"

# Markers (Entrez ID, role)
MARKERS_HM = [
    ("CCK",     885,    "Neuropeptide (CCKBC marker)"),
    ("CNR1",    1268,   "CB1R (CCKBC marker)"),
    ("CXCL14",  9547,   "Sncg-enriched"),
    ("CRH",     1392,   "Sncg-associated"),
    ("SNCG",    6623,   "Subclass marker"),
    ("VIP",     7432,   "VIP+ marker"),
    ("CHRM2",   1129,   "ISI3 marker"),
    ("RELN",    5649,   "Non-VIP CGE / neurogliaform"),
]

print(f"Loading {TPM_FILE.name} ...")
tpm_df = pd.read_csv(TPM_FILE, index_col=0)
print(f"  TPM shape: {tpm_df.shape}")

# Build markers × CGE clusters matrix
CGE_LIST = list(range(276, 297))
marker_mat = pd.DataFrame(
    index=[m[0] for m in MARKERS_HM],
    columns=CGE_LIST,
    dtype=float,
)
for sym, eid, _ in MARKERS_HM:
    if eid not in tpm_df.index:
        continue
    for c in CGE_LIST:
        marker_mat.loc[sym, c] = float(tpm_df.loc[eid, str(c)])

# Build figure
fig_hm = plt.figure(figsize=(11, 5.5))
fig_hm.patch.set_alpha(0)
gs_hm = gridspec.GridSpec(1, 1, figure=fig_hm,
                           left=0.10, right=0.97, top=0.86, bottom=0.18)
ax_hm = fig_hm.add_subplot(gs_hm[0, 0])
style_ax(ax_hm)

log_data_hm = np.log10(marker_mat.values + 1.0)
im = ax_hm.imshow(log_data_hm, aspect="auto", cmap="YlOrRd", interpolation="nearest", zorder=1)

ax_hm.set_xticks(np.arange(len(CGE_LIST)))
ax_hm.set_xticklabels(CGE_LIST, rotation=45, ha="right", fontsize=8)
ax_hm.set_yticks(np.arange(len(MARKERS_HM)))
ax_hm.set_yticklabels([m[0] for m in MARKERS_HM], fontsize=9)
ax_hm.set_xlabel("Human CGE cluster ID", fontsize=10)

# Harmony/scVI candidate columns — same amber as MetaNeighbor; rectangles above imshow
HARMONY_SCVI_CANDIDATES = HUMAN_FP + HUMAN_VIP_NEG_CCKBC + HUMAN_VIP_POS_CCKBC
COL_CANDIDATE_BAND = "#fbb03b"
n_markers = len(MARKERS_HM)
for cid in HARMONY_SCVI_CANDIDATES:
    j = CGE_LIST.index(cid)
    ax_hm.add_patch(
        plt.Rectangle(
            (j - 0.5, -0.5), 1, n_markers,
            facecolor=COL_CANDIDATE_BAND, edgecolor="none",
            alpha=0.22, zorder=5, clip_on=False,
        )
    )
    ax_hm.add_patch(
        plt.Rectangle(
            (j - 0.5, -0.5), 1, n_markers,
            facecolor="none", edgecolor=COL_CANDIDATE_BAND, linewidth=2.25,
            zorder=6, clip_on=False,
        )
    )

cbar = plt.colorbar(im, ax=ax_hm, fraction=0.025, pad=0.02)
cbar.set_label("log₁₀(TPM + 1)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

ax_hm.set_title(
    "CCKBC marker gene expression across human Siletti CGE clusters (pseudobulk TPM)",
    loc="left", fontweight="bold", fontsize=11, pad=24,
)

hm_legend_elems = [
    Patch(facecolor=COL_CANDIDATE_BAND, alpha=0.3,
          label="Harmony/scVI candidates (277–281)"),
]
ax_hm.legend(handles=hm_legend_elems, loc="upper right",
             bbox_to_anchor=(1.0, 1.20), fontsize=7.5, frameon=False)

OUT_HM_PDF = RESULTS_DIR / "fig_supp_cckbc_markers.pdf"
OUT_HM_PNG = RESULTS_DIR / "fig_supp_cckbc_markers.png"
plt.savefig(OUT_HM_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_HM_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_HM_PDF}")
plt.show()

# Print key values for verification
print()
print("Key marker values (TPM) for clusters 277-281:")
print(marker_mat[[277, 278, 279, 280, 281]].round(0).to_string())

# %% [markdown]
# ## Panel E: Cross-species electrophysiology comparison (standalone)
#
# Four groups (2 per species), all from ground-truth transcriptomic labels:
# - **Mouse VIP-other**: Gouwens RNA type in Vip Mybpc1/Gpc3/Chat/C1ql1/Htr1f/Col15a1
# - **Mouse CCKBC**: Gouwens RNA type in Vip Sncg/Vip Serpinf1_{1,2} + Sncg subclass
# - **Human VIP-other**: Lee/Dalley VIP subclass, any transcriptomic type except VIP SERPINF1
# - **Human CCKBC**: Lee/Dalley `Inh L2-5 VIP SERPINF1` (the only CCKBC-like human VIP type)
#
# Human n=3 for CCKBC (small but honest — the Lee/Dalley MTG dataset has few
# SERPINF1 cells).

# %%
groups_to_plot = ["Mouse\nVIP-other", "Mouse\nCCKBC",
                  "Human\nVIP-other", "Human\nCCKBC"]
group_colors = ["#aec7e8", "#1f77b4", "#f4a582", "#d62728"]

# Features selected from the published Scala et al. 2021 ephys feature set
# (mini-atlas M1 patch-seq). ALL 4 significantly discriminate mouse CCKBCs
# from mouse VIP-other cells (Mann-Whitney U, P < 0.001 in published
# values), and each has a well-defined counterpart in the Lee/Dalley human
# MTG ephys feature set. None of the 4 features show a significant
# within-species difference on the human side, i.e., they fail to transfer.
ephys_features_to_plot = [
    # (display name,                mouse col,                        human col,                scale)
    ("ISI CV",                       "ISI coefficient of variation",   "isi_cv_hero",            1.0),
    ("Membrane tau (ms)",            "Membrane time constant (ms)",    "tau",                    1000.0),  # human tau in s → ms
    ("AHP (mV)",                     "Afterhyperpolarization (mV)",    "trough_deltav_hero",     1.0),
    ("Max # APs / firing rate (Hz)", "Max number of APs",              "avg_rate_hero",          1.0),
]


def add_stat(ax, i1, i2, vals1, vals2, y_axes_frac=0.92, fontsize=10, require_min_n=3):
    if require_min_n is not None and (
        len(vals1) < require_min_n or len(vals2) < require_min_n
    ):
        return
    ylim = ax.get_ylim()
    y = ylim[0] + y_axes_frac * (ylim[1] - ylim[0])
    if len(vals1) == 0 or len(vals2) == 0:
        sig = "n.s."
    else:
        try:
            _, p = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        except ValueError:
            sig = "n.s."
    ax.plot([i1+1, i2+1], [y, y], color="black", linewidth=1.0)
    ax.text((i1+i2)/2 + 1, y, sig, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold")


fig_C = plt.figure(figsize=(12, 5))
fig_C.patch.set_alpha(0)
gs_C = gridspec.GridSpec(1, 4, figure=fig_C, wspace=0.50,
                          left=0.06, right=0.98, top=0.86, bottom=0.22)

for fi, (feat_name, mouse_feat, human_feat, scale) in enumerate(ephys_features_to_plot):
    axC = fig_C.add_subplot(gs_C[0, fi])
    style_ax(axC)

    data_groups = []
    data_groups.append(mouse_ephys.loc[mouse_groups == "Mouse VIP-other", mouse_feat].dropna().values)
    data_groups.append(mouse_ephys.loc[mouse_groups == "Mouse CCKBC", mouse_feat].dropna().values)
    data_groups.append(human_ephys.loc[human_groups == "Human VIP-other", human_feat].dropna().values * scale)
    data_groups.append(human_ephys.loc[human_groups == "Human CCKBC", human_feat].dropna().values * scale)

    bp = axC.boxplot(data_groups, tick_labels=groups_to_plot, patch_artist=True,
                     widths=0.6, showfliers=False,
                     medianprops=dict(color="black", linewidth=1.0),
                     boxprops=dict(linewidth=0.4),
                     whiskerprops=dict(linewidth=0.4),
                     capprops=dict(linewidth=0.4))
    for patch, c in zip(bp["boxes"], group_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    for i, vals in enumerate(data_groups):
        if len(vals) > 0:
            jitter = np.random.RandomState(i + fi).uniform(-0.16, 0.16, size=len(vals))
            axC.scatter(np.full(len(vals), i+1) + jitter, vals,
                        s=6, alpha=0.5, color=group_colors[i],
                        edgecolors="white", linewidth=0.3)
    axC.set_ylabel(feat_name, fontsize=9)
    plt.setp(axC.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    all_vals = np.concatenate([v for v in data_groups if len(v) > 0])
    if len(all_vals) > 0:
        # Full min/max so highs (e.g. Human VIP-other) are not clipped vs percentile caps.
        v_lo, v_hi = float(np.min(all_vals)), float(np.max(all_vals))
        span = v_hi - v_lo if v_hi > v_lo else max(abs(v_hi), 1.0) * 0.05
        axC.set_ylim(v_lo - 0.1 * span, v_hi + 0.18 * span)
    # Stat comparisons: mouse CCKBC vs mouse VIP-other; human CCKBC vs human VIP-other
    add_stat(axC, 0, 1, data_groups[0], data_groups[1], y_axes_frac=0.85)
    add_stat(axC, 2, 3, data_groups[2], data_groups[3],
             y_axes_frac=0.95, require_min_n=None)

fig_C.suptitle("E. Cross-species electrophysiology features", fontweight="bold",
               fontsize=11, x=0.06, y=0.96, ha="left")

OUT_C_PDF = RESULTS_DIR / "fig_supp_panelC_ephys.pdf"
OUT_C_PNG = RESULTS_DIR / "fig_supp_panelC_ephys.png"
plt.savefig(OUT_C_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_C_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_C_PDF}")
plt.show()

# %% [markdown]
# ## Panel F: 22q11.2 mutation bias by group (standalone)
#
# Three groups after dropping cluster 277 as a confirmed false positive
# (failed Harmony+scVI agreement, MetaNeighbor mapped to Lamp5, marker
# profile is neurogliaform-like; see Sections 2.2--2.6 of Supp Note 2):
# - **VIP- CCKBC** = cluster 278 (marker-confirmed CCKBC, VIP negative)
# - **VIP+ CCKBC** = clusters 279-281 (Sncg-confirmed CCKBC, VIP positive)
# - **VIP+ ISI** = remaining 16 VIP+ clusters

# %%
fig_D = plt.figure(figsize=(7, 5.4))
fig_D.patch.set_alpha(0)
ax_D = fig_D.add_subplot(111)
style_ax(ax_D)

bias_col = "bias_22q_del"
group_data = {
    "VIP- CCKBC\n(278)": summary[summary["category"] == "VIP- CCKBC (278)"][bias_col].dropna().values,
    "VIP+ CCKBC\n(279-281)": summary[summary["category"] == "VIP+ CCKBC (279-281)"][bias_col].dropna().values,
    "VIP+ ISI\n(others)": summary[summary["category"] == "VIP+ ISI"][bias_col].dropna().values,
}

box_colors_d = [COL_VIP_NEG_CCKBC, COL_VIP_POS_CCKBC, COL_VIP_POS_ISI]
positions = [1, 2, 3]
bp = ax_D.boxplot(list(group_data.values()), positions=positions,
                  tick_labels=list(group_data.keys()), patch_artist=True, widths=0.58,
                  showfliers=False,
                  medianprops=dict(color="black", linewidth=1.4),
                  boxprops=dict(linewidth=0.6),
                  whiskerprops=dict(linewidth=0.6),
                  capprops=dict(linewidth=0.6))
for patch, c in zip(bp["boxes"], box_colors_d):
    patch.set_facecolor(c)
    patch.set_alpha(0.62)

for i, (name, vals) in enumerate(group_data.items()):
    jitter = np.random.RandomState(i).uniform(-0.12, 0.12, size=len(vals))
    ax_D.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                 s=34, color=box_colors_d[i], edgecolors="white", linewidth=0.55,
                 zorder=3)

ax_D.set_ylabel("22q11.2 deletion bias (EFFECT)", fontsize=12)
ax_D.tick_params(axis="both", labelsize=11)
plt.setp(ax_D.get_xticklabels(), fontsize=10, ma="center")

# Gather vals
vals_neg_cc = group_data["VIP- CCKBC\n(278)"]
vals_cckbc = group_data["VIP+ CCKBC\n(279-281)"]
vals_isi = group_data["VIP+ ISI\n(others)"]

# Key comparison: all CCKBC (278+279-281) vs VIP+ ISI
vals_all_cckbc = np.concatenate([vals_neg_cc, vals_cckbc])
if len(vals_all_cckbc) >= 2 and len(vals_isi) >= 2:
    _, p_cckbc_isi = stats.mannwhitneyu(vals_all_cckbc, vals_isi, alternative="two-sided")
else:
    p_cckbc_isi = np.nan

if len(vals_neg_cc) >= 1 and len(vals_isi) >= 2:
    _, p_neg_isi = stats.mannwhitneyu(vals_neg_cc, vals_isi, alternative="two-sided")
else:
    p_neg_isi = np.nan

# Set limits with room for one comparison bracket
all_vals = np.concatenate([vals_neg_cc, vals_cckbc, vals_isi])
ymax = max(np.max(all_vals), 0.2) if len(all_vals) else 0.2
ax_D.set_ylim(-0.02, ymax * 1.34)

# Sub-bracket: 278 + 279–281; main bracket: pooled CCKBC midpoint vs VIP+ ISI
x_left_cckbc = float(positions[0])
x_right_cckbc = float(positions[1])
x_cc_mid = (x_left_cckbc + x_right_cckbc) / 2
x_isi = float(positions[2])
x_mid = (x_cc_mid + x_isi) / 2
ylo, yhi = ax_D.get_ylim()
y_span = yhi - ylo
data_top = float(np.max(all_vals)) if len(all_vals) else ymax

y_pool = data_top + 0.034 * y_span
cap_p = 0.013 * y_span
ax_D.plot([x_left_cckbc, x_left_cckbc], [y_pool - cap_p, y_pool], color="0.42", linewidth=1.1)
ax_D.plot([x_right_cckbc, x_right_cckbc], [y_pool - cap_p, y_pool], color="0.42", linewidth=1.1)
ax_D.plot([x_left_cckbc, x_right_cckbc], [y_pool, y_pool], color="0.42", linewidth=1.1)

y_h = y_pool + 0.058 * y_span
cap = 0.020 * y_span
ax_D.plot([x_cc_mid, x_cc_mid], [y_h - cap, y_h], color="black", linewidth=1.0)
ax_D.plot([x_isi, x_isi], [y_h - cap, y_h], color="black", linewidth=1.0)
ax_D.plot([x_cc_mid, x_isi], [y_h, y_h], color="black", linewidth=1.0)
y_dash_top = y_h + 0.048 * y_span
ax_D.plot([x_mid, x_mid], [y_h, y_dash_top], color="black", linewidth=0.9,
          linestyle="--", zorder=5)
hdr = "Pooled CCKBC (278 + 279–281) vs VIP+ ISI"
if np.isnan(p_cckbc_isi):
    sub = "— n.a."
elif p_cckbc_isi >= 0.05:
    sub = f"— n.s. (P = {p_cckbc_isi:.2f})"
else:
    sub = f"— P = {p_cckbc_isi:.3f}"
sig = f"{hdr}\n{sub}"
y_lbl = y_dash_top + 0.018 * y_span
ax_D.text(x_mid, y_lbl, sig, ha="center", va="bottom", fontsize=9)
new_top = max(yhi, y_lbl + 0.055 * y_span)
ax_D.set_ylim(ylo, new_top)

fig_D.subplots_adjust(left=0.14, right=0.97, top=0.94, bottom=0.16)

OUT_D_PDF = RESULTS_DIR / "fig_supp_panelD_22qbias.pdf"
OUT_D_PNG = RESULTS_DIR / "fig_supp_panelD_22qbias.png"
plt.savefig(OUT_D_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_D_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_D_PDF}")
plt.show()

# %% [markdown]
# ## Supporting figures
#
# Below: stacked-bar quantitative summary of the Panel B validation, plus
# the older 2-panel UMAP version (kept for reference).

# %% [markdown]
# ### Supporting: stacked bar of mouse M1 → WMB-10Xv3 mapping
#
# Quantitative version of the validation in Panel B. Shows the fraction of
# each mouse M1 patch-seq cell type mapping to each WMB GABAergic subclass.
# Sncg-family cells → 92% Sncg Gaba, non-CCKBC VIP → 100% Vip Gaba, etc.

# %%
val_csv = RESULTS_DIR / "cckbc_validation" / "cckbc_atlas_mapping_validation.csv"
val_df = pd.read_csv(val_csv, index_col=0)
print(f"Validation data: {val_df.shape}")

# CCKBC type definitions (Gouwens 2020 Table 1)
CCKBC_RNA_TYPES_VAL = {
    "Sncg Col14a1", "Sncg Slc17a8", "Sncg Calb1_1", "Sncg Calb1_2",
    "Sncg Npy2r",
    "Vip Sncg",
    "Vip Serpinf1_1", "Vip Serpinf1_2", "Vip Serpinf1_3",
}
val_df["is_cckbc_type"] = val_df["RNA_type"].astype(str).isin(CCKBC_RNA_TYPES_VAL) | (
    val_df["RNA_family"].astype(str) == "Sncg"
)

# Input groups: CCKBC subtypes | non-CCKBC CGE | MGE controls
INPUT_GROUPS_VAL = [
    ("Sncg family", val_df[val_df["RNA_family"].astype(str) == "Sncg"]),
    ("Vip Sncg", val_df[val_df["RNA_type"].astype(str) == "Vip Sncg"]),
    ("Vip Serpinf1", val_df[val_df["RNA_type"].astype(str).str.startswith("Vip Serpinf1")]),
    ("Vip non-CCKBC", val_df[(val_df["RNA_family"].astype(str) == "Vip")
                              & (~val_df["is_cckbc_type"])]),
    ("Lamp5", val_df[val_df["RNA_family"].astype(str) == "Lamp5"]),
    ("Pvalb", val_df[val_df["RNA_family"].astype(str) == "Pvalb"]),
    ("Sst", val_df[val_df["RNA_family"].astype(str) == "Sst"]),
]

# WMB subclass color scheme
WMB_COLORS_VAL = {
    "047 Sncg Gaba": "#d62728",
    "046 Vip Gaba": "#9467bd",
    "049 Lamp5 Gaba": "#2ca02c",
    "050 Lamp5 Lhx6 Gaba": "#bcbd22",
    "052 Pvalb Gaba": "#1f77b4",
    "053 Sst Gaba": "#7f7f7f",
    "051 Pvalb chandelier Gaba": "#17becf",
}


def color_for_subclass(label):
    return WMB_COLORS_VAL.get(str(label), "#cccccc")


# Build per-group fraction dict
val_bar_data = {}
for label, sub in INPUT_GROUPS_VAL:
    if len(sub) == 0:
        continue
    val_bar_data[label] = sub["predicted_subclass"].value_counts(normalize=True)

subclass_order_val = [
    "047 Sncg Gaba", "046 Vip Gaba", "049 Lamp5 Gaba", "050 Lamp5 Lhx6 Gaba",
    "052 Pvalb Gaba", "053 Sst Gaba", "051 Pvalb chandelier Gaba",
]
all_predicted_val = sorted(set(val_df["predicted_subclass"].dropna().unique()))
subclass_order_val += [s for s in all_predicted_val if s not in subclass_order_val]

fig_val = plt.figure(figsize=(11, 6.5))
fig_val.patch.set_alpha(0)
gs_val = gridspec.GridSpec(1, 1, figure=fig_val,
                            left=0.08, right=0.74, top=0.84, bottom=0.20)
ax_val = fig_val.add_subplot(gs_val[0, 0])
style_ax(ax_val)

x_labels_val = list(val_bar_data.keys())
x_val = np.arange(len(x_labels_val))
bottom_val = np.zeros(len(x_labels_val))
present_subclasses_val = set()
for sub_label in subclass_order_val:
    vals_arr = np.array([val_bar_data[lbl].get(sub_label, 0.0) for lbl in x_labels_val])
    if vals_arr.sum() == 0:
        continue
    present_subclasses_val.add(sub_label)
    ax_val.bar(x_val, vals_arr, bottom=bottom_val,
               color=color_for_subclass(sub_label),
               edgecolor="white", linewidth=0.5, label=sub_label)
    bottom_val += vals_arr

# Annotate cells with percent
cell_count_per_group_val = {lbl: int(len(sub)) for lbl, sub in INPUT_GROUPS_VAL if len(sub) > 0}
for xi, lbl in enumerate(x_labels_val):
    bot = 0.0
    for sub_label in subclass_order_val:
        v = val_bar_data[lbl].get(sub_label, 0.0)
        if v >= 0.10:
            ax_val.text(xi, bot + v / 2, f"{int(round(v*100))}%",
                        ha="center", va="center", fontsize=7,
                        color="white" if v > 0.18 else "black", weight="bold")
        bot += v

xtick_labels_val = [f"{lbl}\n(n={cell_count_per_group_val[lbl]})" for lbl in x_labels_val]
ax_val.set_xticks(x_val)
ax_val.set_xticklabels(xtick_labels_val, rotation=30, ha="right", fontsize=9)
ax_val.set_ylim(0, 1.18)
ax_val.set_ylabel("Fraction mapped to WMB subclass", fontsize=10)
ax_val.set_title(
    "Validation: mouse M1 patch-seq cell types → WMB-10Xv3 subclass\n"
    "(Harmony integration + k-NN, k=30; reference = WMB-10Xv3 cortex GABAergic)",
    loc="left", fontweight="bold", fontsize=11, pad=20,
)

# Group separators + labels
ax_val.axvline(2.5, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
ax_val.axvline(4.5, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
ax_val.text(1.0, 1.10, "CCKBC subtypes\n(Gouwens 2020)", ha="center",
            fontsize=8.5, fontweight="bold", color="#d62728")
ax_val.text(3.5, 1.10, "Non-CCKBC CGE", ha="center", fontsize=8.5,
            fontweight="bold", color="#2ca02c")
ax_val.text(5.5, 1.10, "MGE controls", ha="center", fontsize=8.5,
            fontweight="bold", color="#1f77b4")

# Legend
val_legend_elems = [
    Patch(facecolor=color_for_subclass(s), edgecolor="white", label=s)
    for s in subclass_order_val if s in present_subclasses_val
]
ax_val.legend(handles=val_legend_elems, loc="center left",
              bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False,
              title="Predicted WMB subclass", title_fontsize=8)

OUT_VAL_PDF = RESULTS_DIR / "fig_supp_sncg_validation.pdf"
OUT_VAL_PNG = RESULTS_DIR / "fig_supp_sncg_validation.png"
plt.savefig(OUT_VAL_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_VAL_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_VAL_PDF}")
plt.show()

# %% [markdown]
# ### Validation UMAP overlay
#
# UMAP of the same Harmony-corrected within-mouse mapping. WMB-10Xv3 cortical
# GABAergic atlas cells (background) colored by subclass; M1 patch-seq query
# cells (foreground) colored by Gouwens RNA family/CCKBC label. Mouse Sncg/CCKBC
# query cells should land directly on the WMB Sncg Gaba atlas region.
#
# Pre-computed by `cge_subtype/scripts/compute_validation_umap.py`.

# %%
import scanpy as sc  # noqa: E402

val_umap_h5ad = RESULTS_DIR / "sncg_validation_umap.h5ad"
adata_val_umap = sc.read_h5ad(val_umap_h5ad)
print(f"Validation UMAP data: {adata_val_umap.shape}")
print(f"Source: {adata_val_umap.obs['source'].value_counts().to_dict()}")

umap_v = adata_val_umap.obsm["X_umap"]
obs_v = adata_val_umap.obs.copy()

ref_mask_v = (obs_v["source"] == "atlas").values
query_mask_v = (obs_v["source"] == "patchseq").values

# Subclass color scheme — chosen to be maximally distinct on a single panel.
# Key changes from the original tab10 palette: Vip is hot magenta (was purple,
# too close to Pvalb blue), Sst is dark brown (was gray, washed out).
SUBCLASS_COLORS_VAL = {
    "047 Sncg Gaba":            "#d62728",  # red
    "046 Vip Gaba":             "#e91e63",  # magenta/hot pink
    "049 Lamp5 Gaba":           "#2ca02c",  # green
    "050 Lamp5 Lhx6 Gaba":      "#bcbd22",  # olive
    "052 Pvalb Gaba":           "#1f77b4",  # blue
    "053 Sst Gaba":             "#8c564b",  # dark brown
    "048 RHP-COA Ndnf Gaba":    "#17becf",  # cyan
    "051 Pvalb chandelier Gaba": "#9467bd",  # purple
}

# Query (M1 patch-seq) coloring matches contour colors so M1 dots and atlas
# contours are visually linked
QUERY_COLORS = {
    "Sncg":   "#d62728",   # red — CCKBC
    "Vip":    "#e91e63",   # magenta — Vip
    "Lamp5":  "#2ca02c",   # green
    "Pvalb":  "#1f77b4",   # blue
    "Sst":    "#8c564b",   # brown
}


def query_color(family, is_cckbc):
    if str(is_cckbc) == "1.0":
        return "#fbb03b"  # amber — highlight CCKBCs
    return QUERY_COLORS.get(str(family), "#cccccc")


# Build figure: 2 panels — left = atlas reference colored by subclass,
# right = same UMAP with M1 patch-seq overlaid
fig_val_umap = plt.figure(figsize=(13, 6))
fig_val_umap.patch.set_alpha(0)
gs_vu = gridspec.GridSpec(1, 2, figure=fig_val_umap, wspace=0.25,
                           left=0.05, right=0.97, top=0.92, bottom=0.08)

# --- Panel left: atlas reference colored by subclass ---
ax_l = fig_val_umap.add_subplot(gs_vu[0, 0])
style_ax(ax_l)

# Plot atlas cells by subclass (in order to handle layering)
ref_subclass = obs_v.loc[ref_mask_v, "subclass"].astype(str).values
for sc_name in sorted(set(ref_subclass)):
    sel = ref_mask_v.copy()
    sel[ref_mask_v] = (ref_subclass == sc_name)
    color = SUBCLASS_COLORS_VAL.get(sc_name, "#cccccc")
    ax_l.scatter(umap_v[sel, 0], umap_v[sel, 1],
                 s=2, c=color, alpha=0.5, rasterized=True,
                 edgecolors="none", label=sc_name if sc_name in SUBCLASS_COLORS_VAL else None)

ax_l.set_xlabel("UMAP 1")
ax_l.set_ylabel("UMAP 2")
ax_l.set_title("WMB-10Xv3 cortical GABAergic atlas\n(colored by subclass)",
               fontweight="bold", fontsize=10, loc="left")

# Subclass legend
from matplotlib.lines import Line2D
present_subs = [s for s in SUBCLASS_COLORS_VAL if s in ref_subclass]
legend_l = [Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=SUBCLASS_COLORS_VAL[s], markersize=6,
                   label=s, linestyle="none")
            for s in present_subs]
ax_l.legend(handles=legend_l, loc="upper right", fontsize=6.5,
            frameon=False)

# --- Panel right: atlas (light) + M1 patch-seq overlay ---
ax_r = fig_val_umap.add_subplot(gs_vu[0, 1])
style_ax(ax_r)

# Atlas cells faded as background
ax_r.scatter(umap_v[ref_mask_v, 0], umap_v[ref_mask_v, 1],
             s=1.5, c="#e0e0e0", alpha=0.4, rasterized=True, edgecolors="none")

# Highlight Sncg Gaba atlas cells as a faint anchor (so reader can see "this is the Sncg region")
sncg_atlas_mask = ref_mask_v.copy()
sncg_atlas_mask[ref_mask_v] = (ref_subclass == "047 Sncg Gaba")
ax_r.scatter(umap_v[sncg_atlas_mask, 0], umap_v[sncg_atlas_mask, 1],
             s=2, c="#d62728", alpha=0.45, rasterized=True, edgecolors="none",
             label="WMB 047 Sncg Gaba (atlas)")

# M1 patch-seq query cells, colored by RNA family / CCKBC status
q_family = obs_v.loc[query_mask_v, "RNA_family"].astype(str).values
q_iscckbc = obs_v.loc[query_mask_v, "is_cckbc"].astype(str).values

# Plot non-CCKBC query cells first (background), then CCKBC on top
non_cckbc_idx = q_iscckbc != "1.0"
cckbc_idx = q_iscckbc == "1.0"

for fam in ["Sst", "Pvalb", "Lamp5", "Vip"]:
    sel = non_cckbc_idx & (q_family == fam)
    if sel.sum() == 0:
        continue
    coords = umap_v[query_mask_v][sel]
    ax_r.scatter(coords[:, 0], coords[:, 1],
                 s=8, c=QUERY_COLORS[fam], alpha=0.7,
                 edgecolors="white", linewidth=0.3,
                 label=f"M1 {fam} (n={int(sel.sum())})")

# CCKBC overlay (Sncg + Vip Sncg/Serpinf1 — labeled True)
coords_cckbc = umap_v[query_mask_v][cckbc_idx]
ax_r.scatter(coords_cckbc[:, 0], coords_cckbc[:, 1],
             s=22, c="#fbb03b", alpha=0.95,
             edgecolors="black", linewidth=0.6, zorder=5,
             label=f"M1 CCKBC (n={int(cckbc_idx.sum())})")

ax_r.set_xlabel("UMAP 1")
ax_r.set_ylabel("UMAP 2")
ax_r.set_title("M1 patch-seq cells overlaid on WMB atlas\n(CCKBC cells highlighted in amber)",
               fontweight="bold", fontsize=10, loc="left")
ax_r.legend(loc="upper right", fontsize=6.5, frameon=False)

OUT_VAL_UMAP_PDF = RESULTS_DIR / "fig_supp_sncg_validation_umap.pdf"
OUT_VAL_UMAP_PNG = RESULTS_DIR / "fig_supp_sncg_validation_umap.png"
plt.savefig(OUT_VAL_UMAP_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_VAL_UMAP_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_VAL_UMAP_PDF}")
plt.show()

# %% [markdown]
# ## UMAP: Mouse CCKBC cells in cross-species CGE space
#
# **Filter:** Mouse query restricted to **M1 dataset, CGE families only**
# (Sncg + Vip + Lamp5). This gives 257 cells: 51 CCKBC + 206 non-CCKBC
# (115 Vip-other + 91 Lamp5).
#
# Two cross-species integration methods compared side-by-side:
# - **scVI UMAP** (top row): scVI latent space (50D), pre-trained on joint
#   mouse+human reference
# - **Harmony UMAP** (bottom row): Harmony-corrected PCA space with theta=4
#
# **Caveat: cross-species integration is incomplete.** Even after batch
# correction, mouse and human cells live in partially separate manifolds
# in latent space. Distance ratio: mouse CCKBC → nearest human is ~1.8×
# farther than mouse CCKBC → nearest other mouse cell. This is consistent
# with the known poor cross-species conservation of Sncg/CGE interneurons
# (Bakken et al. 2021). The k-NN-based mapping in Panel A picks the
# *least-far* human clusters, which are clusters 280, 281, 279 — matching
# the MetaNeighbor result in Panel B.
#
# Pre-computed by:
# - `cge_subtype/scripts/compute_cge_harmony_umap.py` (scVI)
# - `cge_subtype/scripts/compute_cge_harmony_umap_v2.py` (Harmony, theta=4)

# %%
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

UMAP_SCVI_H5AD = RESULTS_DIR / "harmony_umap_cge.h5ad"           # scVI latent
UMAP_HARMONY_H5AD = RESULTS_DIR / "harmony_umap_cge_harmony.h5ad"  # Harmony

adata_scvi = sc.read_h5ad(UMAP_SCVI_H5AD)
adata_harmony = sc.read_h5ad(UMAP_HARMONY_H5AD)
print(f"scVI UMAP:    {adata_scvi.shape}")
print(f"Harmony UMAP: {adata_harmony.shape}")
# Use scvi for backwards-compat with downstream code
adata_umap = adata_scvi
print(f"Source breakdown: {adata_umap.obs['source'].value_counts().to_dict()}")
print(f"Species breakdown: {adata_umap.obs['species'].value_counts().to_dict()}")

# %% [markdown]
# **Note on UMAP vs latent-space mapping discrepancy:** Panel A's "Mouse CCKBC
# fraction" is computed via k-NN voting in the **full 50-dim scVI latent space**.
# The UMAP below is a 2D projection of that same latent space, but UMAP can
# distort cross-species local neighborhoods.
#
# Direct verification (k-NN in scVI latent, k=30): top 3 clusters mouse CCKBCs
# map to are **280 (n=81), 281 (n=46), 279 (n=35)** — the Sncg-confirmed CCKBC
# clusters from Panel A. 50.5% of CCKBCs land in clusters 277-281. So the
# Panel A result is correct; the 2D UMAP just visually compresses cross-species
# distances.

# %%
# Build the UMAP comparison figure: scVI (top row) + Harmony (bottom row)
from matplotlib.lines import Line2D

cat_colors = {
    "Other": "#e0e0e0",
    "VIP+ ISI": COL_VIP_POS_ISI,
    "FP (277)": COL_FP,
    "VIP- CCKBC (278)": COL_VIP_NEG_CCKBC,
    "VIP+ CCKBC (279-281)": COL_VIP_POS_CCKBC,
}


def cluster_to_category(cid_str):
    try:
        cid = int(float(cid_str))
    except (ValueError, TypeError):
        return "Other"
    if cid in HUMAN_FP:
        return "FP (277)"
    elif cid in HUMAN_VIP_NEG_CCKBC:
        return "VIP- CCKBC (278)"
    elif cid in HUMAN_VIP_POS_CCKBC:
        return "VIP+ CCKBC (279-281)"
    elif cid in HUMAN_VIP_POS_ISI:
        return "VIP+ ISI"
    else:
        return "Other"


def plot_umap_panels(adata, ax_left, ax_right, method_label):
    """Plot two UMAP panels: clusters with CCKBC overlay (left) and CCKBC vs non (right)."""
    umap = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["cge_category"] = obs["cluster_id"].apply(cluster_to_category)

    # Convert is_cckbc to bool consistently (sometimes string, sometimes bool)
    is_cckbc = obs["is_cckbc"].astype(str).isin(["True", "true", "1"]).values

    ref_human_mask = ((obs["source"] == "reference") & (obs["species"] == "human")).values
    ref_mouse_mask = ((obs["source"] == "reference") & (obs["species"] != "human")).values
    query_mask = (obs["source"] == "query").values
    query_cckbc_mask = query_mask & is_cckbc
    query_noncckbc_mask = query_mask & (~is_cckbc)

    # --- Left subpanel: clusters + CCKBC overlay ---
    style_ax(ax_left)
    for cat in ["Other", "VIP+ ISI", "FP (277)", "VIP- CCKBC (278)", "VIP+ CCKBC (279-281)"]:
        sel = ref_human_mask & (obs["cge_category"] == cat).values
        if sel.sum() > 0:
            ax_left.scatter(umap[sel, 0], umap[sel, 1],
                            s=2, c=cat_colors[cat], alpha=0.5, rasterized=True,
                            edgecolors="none")
    ax_left.scatter(umap[ref_mouse_mask, 0], umap[ref_mouse_mask, 1],
                    s=2, c="#9e9e9e", alpha=0.20, rasterized=True, edgecolors="none")
    ax_left.scatter(umap[query_cckbc_mask, 0], umap[query_cckbc_mask, 1],
                    s=18, c="#fbb03b", alpha=0.95,
                    edgecolors="black", linewidth=0.4, zorder=5)
    ax_left.set_xlabel("UMAP 1", fontsize=8)
    ax_left.set_ylabel("UMAP 2", fontsize=8)
    ax_left.set_title(f"{method_label}: clusters + CCKBC overlay",
                      fontweight="bold", fontsize=9, loc="left")

    # --- Right subpanel: mouse query CCKBC vs non-CCKBC ---
    style_ax(ax_right)
    ref_all_mask = ref_human_mask | ref_mouse_mask
    ax_right.scatter(umap[ref_all_mask, 0], umap[ref_all_mask, 1],
                     s=1.5, c="#e0e0e0", alpha=0.4, rasterized=True, edgecolors="none")
    ax_right.scatter(umap[query_noncckbc_mask, 0], umap[query_noncckbc_mask, 1],
                     s=6, c="#1f77b4", alpha=0.5, edgecolors="none", rasterized=True)
    ax_right.scatter(umap[query_cckbc_mask, 0], umap[query_cckbc_mask, 1],
                     s=18, c="#fbb03b", alpha=0.95,
                     edgecolors="black", linewidth=0.4, zorder=5)
    ax_right.set_xlabel("UMAP 1", fontsize=8)
    ax_right.set_ylabel("UMAP 2", fontsize=8)
    ax_right.set_title(f"{method_label}: mouse query CCKBC vs non-CCKBC",
                       fontweight="bold", fontsize=9, loc="left")
    return query_cckbc_mask.sum(), query_noncckbc_mask.sum()


fig_umap = plt.figure(figsize=(11, 9), constrained_layout=False)
fig_umap.patch.set_alpha(0)
gs_u = gridspec.GridSpec(2, 2, figure=fig_umap, wspace=0.28, hspace=0.35,
                          left=0.07, right=0.97, top=0.94, bottom=0.07)

ax_scvi_l = fig_umap.add_subplot(gs_u[0, 0])
ax_scvi_r = fig_umap.add_subplot(gs_u[0, 1])
ax_harm_l = fig_umap.add_subplot(gs_u[1, 0])
ax_harm_r = fig_umap.add_subplot(gs_u[1, 1])

n_cckbc, n_noncckbc = plot_umap_panels(adata_scvi, ax_scvi_l, ax_scvi_r, "scVI")
plot_umap_panels(adata_harmony, ax_harm_l, ax_harm_r, "Harmony")

# Shared legend at top of figure
legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=cat_colors["VIP+ ISI"],
           markersize=6, label="Human VIP+ ISI", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=cat_colors["FP (277)"],
           markersize=6, label="Human FP (277)", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=cat_colors["VIP- CCKBC (278)"],
           markersize=6, label="Human VIP- CCKBC (278)", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=cat_colors["VIP+ CCKBC (279-281)"],
           markersize=6, label="Human VIP+ CCKBC (279-281)", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#e0e0e0",
           markersize=6, label="Human other CGE", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#9e9e9e",
           markersize=6, label="Mouse reference", linestyle="none"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f77b4",
           markersize=6, label=f"Mouse query non-CCKBC (n={n_noncckbc})", linestyle="none"),
    Line2D([0], [0], marker="o", color="black", markerfacecolor="#fbb03b",
           markersize=8, markeredgewidth=0.5,
           label=f"Mouse CCKBC query (n={n_cckbc})", linestyle="none"),
]
fig_umap.legend(handles=legend_handles, loc="upper center",
                bbox_to_anchor=(0.5, 1.005), ncol=4, fontsize=7,
                frameon=False)

# Save
UMAP_PDF = RESULTS_DIR / "fig_supp_cckbc_umap.pdf"
UMAP_PNG = RESULTS_DIR / "fig_supp_cckbc_umap.png"
plt.savefig(UMAP_PDF, transparent=True, bbox_inches="tight")
plt.savefig(UMAP_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {UMAP_PDF}")
print(f"Saved: {UMAP_PNG}")
plt.show()

# %%
# Verify CCKBC convergence: which clusters do CCKBCs land near?
# (using harmony_patchseq_mapping_results to get the assigned cluster)
mapping = pd.read_csv(RESULTS_DIR / "harmony_patchseq_mapping_results.csv", index_col=0)
cckbc_assign = mapping[mapping["is_cckbc"] == True]
print(f"Mouse CCKBC assigned clusters (top 10):")
print(cckbc_assign["assigned_cluster"].value_counts().head(10))
print()
in_target = cckbc_assign["assigned_cluster"].isin(
    HUMAN_VIP_NEG + HUMAN_VIP_POS_CCKBC
).sum()
print(f"CCKBCs landing in target clusters (277-281): {in_target}/{len(cckbc_assign)} "
      f"({100*in_target/len(cckbc_assign):.1f}%)")

# %% [markdown]
# ## Summary statistics

# %%
print(f"Mouse CCKBC (M1 patch-seq):                       n = {(mouse_groups=='Mouse CCKBC').sum()}")
print(f"Mouse VIP-other:                                  n = {(mouse_groups=='Mouse VIP-other').sum()}")
print(f"Human CCKBC (Lee/Dalley VIP SERPINF1):            n = {(human_groups=='Human CCKBC').sum()}")
print(f"Human VIP-other (Lee/Dalley VIP, non-SERPINF1):   n = {(human_groups=='Human VIP-other').sum()}")
print()
print("22q deletion bias (Mann-Whitney, after dropping 277 as FP):")
print(f"  All CCKBC (278+279-281) vs VIP+ ISI:  P = {p_cckbc_isi:.3f}")
print(f"  VIP- CCKBC (278) vs VIP+ ISI:         P = {p_neg_isi:.3f}")
print()
print("Confirmed Sncg homologs (best mouse subclass):")
for _, row in summary[summary["mouse_subclass_simple"] == "Sncg"].iterrows():
    print(f"  Cluster {row['cluster_id']}: AUROC={row['best_auroc']:.3f}")
