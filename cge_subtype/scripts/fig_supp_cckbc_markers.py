#!/usr/bin/env python3
"""
Build supplementary figure: CCKBC marker gene expression across human CGE clusters.

Shows pseudobulk TPM expression of canonical CCKBC and discriminating markers
across all 21 human Siletti CGE clusters (276-296), to support the marker-level
classification:
  - 279, 280, 281: high CCK + CNR1 + low VIP/CHRM2 → confident CCKBC
  - 278: high CCK + CNR1 + high CXCL14/CRH + VIP-negative → marker-defined CCKBC
  - 277: low CCK, low CNR1, very high CHRM2/RELN → neurogliaform-like, NOT CCKBC

Output:
  cge_subtype/results/fig_supp_cckbc_markers.{pdf,png}
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

PROJECT_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
TPM_FILE = PROJECT_DIR / "dat/ExpMats/HumanCT.TPM.0.1.Filt.csv"
OUT_PDF = PROJECT_DIR / "cge_subtype/results/fig_supp_cckbc_markers.pdf"
OUT_PNG = PROJECT_DIR / "cge_subtype/results/fig_supp_cckbc_markers.png"

# Markers and their Entrez IDs
MARKERS = [
    ("CCK",     885,    "Neuropeptide (CCKBC marker)"),
    ("CNR1",    1268,   "CB1R (CCKBC marker)"),
    ("CXCL14",  9547,   "Sncg-enriched"),
    ("CRH",     1392,   "Sncg-associated"),
    ("SNCG",    6623,   "Subclass marker"),
    ("VIP",     7432,   "VIP+ marker"),
    ("CHRM2",   1129,   "ISI3 marker"),
    ("RELN",    5649,   "Non-VIP CGE / neurogliaform"),
]

CGE_CLUSTERS = list(range(276, 297))

# Cluster categories (matching the figure notebook)
HUMAN_FP = [277]
HUMAN_VIP_NEG_CCKBC = [278]
HUMAN_VIP_POS_CCKBC = [279, 280, 281]

COL_FP = "#7f7f7f"
COL_VIP_NEG_CCKBC = "#ff7f0e"
COL_VIP_POS_CCKBC = "#d62728"
COL_VIP_POS_ISI = "#1f77b4"
HARMONY_SCVI_CANDIDATES = HUMAN_FP + HUMAN_VIP_NEG_CCKBC + HUMAN_VIP_POS_CCKBC
COL_CANDIDATE_BAND = "#fbb03b"  # match MetaNeighbor panel (Harmony/scVI candidates)

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
# Load TPM and extract marker values for the 21 CGE clusters
# ---------------------------------------------------------------------------
print(f"Loading {TPM_FILE} ...")
tpm = pd.read_csv(TPM_FILE, index_col=0)
print(f"  TPM shape: {tpm.shape}")

# Build a (markers × CGE clusters) DataFrame
marker_matrix = pd.DataFrame(
    index=[m[0] for m in MARKERS],
    columns=CGE_CLUSTERS,
    dtype=float,
)
for sym, eid, _ in MARKERS:
    if eid not in tpm.index:
        print(f"  WARNING: {sym} ({eid}) not in TPM file")
        continue
    for c in CGE_CLUSTERS:
        col = str(c)
        marker_matrix.loc[sym, c] = float(tpm.loc[eid, col])

print(f"\nMarker matrix shape: {marker_matrix.shape}")

# ---------------------------------------------------------------------------
# Build figure: heatmap (markers × clusters) with cluster category bands
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 5.5))
fig.patch.set_alpha(0)

gs = gridspec.GridSpec(1, 1, figure=fig,
                        left=0.10, right=0.97, top=0.86, bottom=0.18)
ax = fig.add_subplot(gs[0, 0])
style_ax(ax)

# Use log10(TPM + 1) for visualization (TPM ranges 0 to >5000)
log_data = np.log10(marker_matrix.values + 1.0)

# Heatmap
im = ax.imshow(log_data, aspect="auto", cmap="YlOrRd", interpolation="nearest", zorder=1)

# Tick labels
ax.set_xticks(np.arange(len(CGE_CLUSTERS)))
ax.set_xticklabels(CGE_CLUSTERS, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(len(MARKERS)))
ax.set_yticklabels([m[0] for m in MARKERS], fontsize=9)
ax.set_xlabel("Human CGE cluster ID", fontsize=10)

# Harmony/scVI candidate columns — rectangles above imshow for reliable highlighting
n_markers_hm = len(MARKERS)
for cid in HARMONY_SCVI_CANDIDATES:
    j = CGE_CLUSTERS.index(cid)
    ax.add_patch(
        plt.Rectangle(
            (j - 0.5, -0.5), 1, n_markers_hm,
            facecolor=COL_CANDIDATE_BAND, edgecolor="none",
            alpha=0.22, zorder=5, clip_on=False,
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (j - 0.5, -0.5), 1, n_markers_hm,
            facecolor="none", edgecolor=COL_CANDIDATE_BAND, linewidth=2.25,
            zorder=6, clip_on=False,
        )
    )

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("log₁₀(TPM + 1)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# Title and legend
ax.set_title(
    "CCKBC marker gene expression across human Siletti CGE clusters (pseudobulk TPM)",
    loc="left", fontweight="bold", fontsize=11, pad=24,
)

legend_elems = [
    Patch(facecolor=COL_CANDIDATE_BAND, alpha=0.3,
          label="Harmony/scVI candidates (277–281)"),
]
ax.legend(handles=legend_elems, loc="upper right",
          bbox_to_anchor=(1.0, 1.20), fontsize=7.5, frameon=False)

plt.savefig(OUT_PDF, transparent=True, bbox_inches="tight")
plt.savefig(OUT_PNG, transparent=True, bbox_inches="tight", dpi=150)
print(f"\nSaved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# Print key values for verification
print()
print("=" * 70)
print("Key marker values (TPM) for clusters 277-281:")
print("=" * 70)
print(marker_matrix[[277, 278, 279, 280, 281]].round(0).to_string())
