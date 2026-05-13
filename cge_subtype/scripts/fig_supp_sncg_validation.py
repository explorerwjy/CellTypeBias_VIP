#!/usr/bin/env python3
"""
Build supplementary figure: validation that Sncg subclass = CCKBC.

Maps mouse M1 patch-seq cells (Gouwens 2020 / Scala 2021 with known
RNA family/type labels) onto the Allen WMB-10Xv3 cortical GABAergic atlas
via Harmony integration + k-NN. Shows that:
  1. Sncg-family cells overwhelmingly land in WMB Sncg Gaba subclass
  2. Vip Sncg / Vip Serpinf1 supertypes (small/deep CCKBCs) split between
     Sncg Gaba and Vip Gaba — consistent with Vip-Sncg boundary
  3. Non-CCKBC VIP cells exclusively map to Vip Gaba (clean separation)
  4. Pvalb / Sst / Lamp5 controls map to their cognate WMB subclasses

This independently validates the Sncg-subclass = CCKBC mapping definition
used throughout Section 2 of the supplementary note.

Output:
  cge_subtype/results/fig_supp_sncg_validation.{pdf,png}
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

PROJECT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype")
RESULTS_DIR = PROJECT_DIR / "results"
VAL_CSV = RESULTS_DIR / "cckbc_validation" / "cckbc_atlas_mapping_validation.csv"

OUT_PDF = RESULTS_DIR / "fig_supp_sncg_validation.pdf"
OUT_PNG = RESULTS_DIR / "fig_supp_sncg_validation.png"

# CCKBC type definitions (from Gouwens 2020 Supplementary Table 1)
CCKBC_RNA_TYPES = {
    "Sncg Col14a1", "Sncg Slc17a8", "Sncg Calb1_1", "Sncg Calb1_2",
    "Sncg Npy2r",
    "Vip Sncg",
    "Vip Serpinf1_1", "Vip Serpinf1_2", "Vip Serpinf1_3",
}

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
# Load and process
# ---------------------------------------------------------------------------
print(f"Loading {VAL_CSV} ...")
df = pd.read_csv(VAL_CSV, index_col=0)
print(f"  rows: {len(df)}")

# Filter to CGE-related families (Sncg, Vip, Lamp5) — same as our analysis filter
CGE_FAMILIES = ["Sncg", "Vip", "Lamp5"]
df["is_cge"] = df["RNA_family"].astype(str).isin(CGE_FAMILIES)

# Add MGE controls (Pvalb, Sst) for the right panel
MGE_FAMILIES = ["Pvalb", "Sst"]
df["is_mge"] = df["RNA_family"].astype(str).isin(MGE_FAMILIES)

# CCKBC tag from RNA_type
df["is_cckbc"] = df["RNA_type"].astype(str).isin(CCKBC_RNA_TYPES) | (
    df["RNA_family"].astype(str) == "Sncg"
)

print(f"  CGE cells: {df['is_cge'].sum()}, of which CCKBC: {df.loc[df['is_cge'], 'is_cckbc'].sum()}")

# ---------------------------------------------------------------------------
# Build the per-cell-type confusion matrix
# ---------------------------------------------------------------------------
# Group input cells by RNA family/type, then count predicted WMB subclasses
def simplify_predicted(s):
    """Return a short subclass label from '047 Sncg Gaba' etc."""
    return str(s)


# Define the input groups we want to display (in order)
INPUT_GROUPS = [
    # CCKBC subtypes (from Gouwens definition)
    ("Sncg family", df[df["RNA_family"].astype(str) == "Sncg"]),
    ("Vip Sncg", df[df["RNA_type"].astype(str) == "Vip Sncg"]),
    ("Vip Serpinf1", df[df["RNA_type"].astype(str).str.startswith("Vip Serpinf1")]),
    # Non-CCKBC CGE controls
    ("Vip non-CCKBC", df[(df["RNA_family"].astype(str) == "Vip")
                          & (~df["is_cckbc"])]),
    ("Lamp5", df[df["RNA_family"].astype(str) == "Lamp5"]),
    # MGE controls
    ("Pvalb", df[df["RNA_family"].astype(str) == "Pvalb"]),
    ("Sst", df[df["RNA_family"].astype(str) == "Sst"]),
]

# WMB subclass color scheme (matches Panel B convention)
WMB_COLORS = {
    "047 Sncg Gaba": "#d62728",       # red — CCKBC homolog
    "046 Vip Gaba": "#9467bd",         # purple — VIP
    "049 Lamp5 Gaba": "#2ca02c",       # green — Lamp5
    "050 Lamp5 Lhx6 Gaba": "#bcbd22",  # olive — Lamp5 LHX6
    "052 Pvalb Gaba": "#1f77b4",       # blue — Pvalb
    "053 Sst Gaba": "#7f7f7f",         # gray — Sst
    "048 RHP-COA Ndnf Gaba": "#e377c2",
    "051 Pvalb chandelier Gaba": "#17becf",
}


def get_color(subclass_label):
    return WMB_COLORS.get(str(subclass_label), "#cccccc")


# Build the stacked-bar matrix: rows = input groups, columns = WMB subclasses
all_predicted = sorted(set(df["predicted_subclass"].dropna().unique()))
print(f"  Predicted WMB subclasses: {all_predicted}")

bar_data = {}
for label, sub in INPUT_GROUPS:
    if len(sub) == 0:
        continue
    counts = sub["predicted_subclass"].value_counts(normalize=True)
    bar_data[label] = counts


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 6.5), constrained_layout=False)
fig.patch.set_alpha(0)
gs = gridspec.GridSpec(1, 1, figure=fig,
                        left=0.08, right=0.74, top=0.84, bottom=0.20)
ax = fig.add_subplot(gs[0, 0])
style_ax(ax)

# Order subclasses: CCKBC-relevant first
subclass_order = [
    "047 Sncg Gaba",
    "046 Vip Gaba",
    "049 Lamp5 Gaba",
    "050 Lamp5 Lhx6 Gaba",
    "052 Pvalb Gaba",
    "053 Sst Gaba",
]
subclass_order += [s for s in all_predicted if s not in subclass_order]

x_labels = list(bar_data.keys())
x = np.arange(len(x_labels))
bottom = np.zeros(len(x_labels))

# Track which subclasses actually appear (for legend)
present_subclasses = set()
for sc in subclass_order:
    vals = np.array([bar_data[lbl].get(sc, 0.0) for lbl in x_labels])
    if vals.sum() == 0:
        continue
    present_subclasses.add(sc)
    ax.bar(x, vals, bottom=bottom, color=get_color(sc),
           edgecolor="white", linewidth=0.5, label=sc)
    # Cell-count annotations on segments > 0.10
    for xi, v in enumerate(vals):
        if v >= 0.10:
            n_cells = int(round(v * sum(bar_data[x_labels[xi]].sum() for _ in [0])
                                   * 0))  # placeholder; recompute below
    bottom += vals

# Re-annotate using true cell counts (not normalized)
cell_count_per_group = {lbl: int(len(sub)) for lbl, sub in INPUT_GROUPS if len(sub) > 0}
for xi, lbl in enumerate(x_labels):
    n_total = cell_count_per_group[lbl]
    counts_abs = (df[df["RNA_family"].astype(str) == "Sncg"]  # placeholder
                    if False else None)
    bottom_acc = 0.0
    for sc in subclass_order:
        v = bar_data[lbl].get(sc, 0.0)
        if v >= 0.10:
            n_seg = int(round(v * n_total))
            ax.text(xi, bottom_acc + v / 2, f"{int(round(v*100))}%",
                    ha="center", va="center", fontsize=7,
                    color="white" if v > 0.18 else "black", weight="bold")
        bottom_acc += v

# X-axis: input group names with cell counts
xtick_labels = [f"{lbl}\n(n={cell_count_per_group[lbl]})" for lbl in x_labels]
ax.set_xticks(x)
ax.set_xticklabels(xtick_labels, rotation=30, ha="right", fontsize=9)
ax.set_ylim(0, 1.18)  # leave headroom for group labels
ax.set_ylabel("Fraction mapped to WMB subclass", fontsize=10)
ax.set_title(
    "Validation: mouse M1 patch-seq cell types → WMB-10Xv3 subclass\n"
    "(Harmony integration + k-NN, k=30; reference = WMB-10Xv3 cortex GABAergic)",
    loc="left", fontweight="bold", fontsize=11, pad=20,
)

# Vertical separators between conceptual groups (CCKBC | non-CCKBC CGE | MGE)
ax.axvline(2.5, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
ax.axvline(4.5, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
ax.text(1.0, 1.10, "CCKBC subtypes\n(Gouwens 2020)", ha="center", fontsize=8.5,
        fontweight="bold", color="#d62728")
ax.text(3.5, 1.10, "Non-CCKBC CGE", ha="center", fontsize=8.5,
        fontweight="bold", color="#2ca02c")
ax.text(5.5, 1.10, "MGE controls", ha="center", fontsize=8.5,
        fontweight="bold", color="#1f77b4")

# Legend (only present subclasses, in display order)
from matplotlib.patches import Patch
legend_elems = [
    Patch(facecolor=get_color(sc), edgecolor="white", label=sc)
    for sc in subclass_order if sc in present_subclasses
]
ax.legend(handles=legend_elems, loc="center left",
          bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False,
          title="Predicted WMB subclass", title_fontsize=8)

plt.savefig(OUT_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# Print summary
print()
print("=" * 60)
print("Summary: M1 patch-seq → WMB subclass")
print("=" * 60)
for label, sub in INPUT_GROUPS:
    if len(sub) == 0:
        continue
    top = sub["predicted_subclass"].value_counts().head(2)
    top_str = "; ".join([f"{name}: {n}" for name, n in top.items()])
    print(f"  {label} (n={len(sub)}): {top_str}")
