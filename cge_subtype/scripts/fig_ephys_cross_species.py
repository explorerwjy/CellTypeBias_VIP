# ---
# jupyter:
#   jupytext:
#     formats: cge_subtype/notebooks//ipynb,cge_subtype/scripts//py:percent
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

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Cross-Species Ephys Comparison: Mouse CCKBC vs VIP → Human Transfer
#
# **Goal:** Show that electrophysiological features discriminating mouse CCKBC
# from mouse VIP-other do NOT transfer to human VIP subtypes.
#
# **Panel A (left):** Mouse CCKBC vs VIP-other — features that discriminate (p<0.01)
# **Panel A (right):** Human VIP-SERPINF1 (CCKBC-like) vs other VIP — same features, no discrimination

# %%
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

# %%
# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")

# Mouse M1 patchseq: ephys features + transcriptomic labels
M1_EPHYS = Path("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_ephys_features.csv")
HARMONY_MAP = REPO / "cge_subtype/results/harmony_patchseq_mapping_results.csv"

# Human patchseq: ephys features + ground truth labels
HUMAN_EPHYS = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")
HUMAN_META = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_manuscript_metadata.csv")

# Raw traces
MOUSE_DANDI = Path("/mnt/data0/DANDI/Processed/000008")
HUMAN_DANDI = Path("/mnt/data0/DANDI/Processed/000636")

OUTDIR = REPO / "results/figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

# %%
# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
COLOR_CCKBC = "#D55E00"      # vermilion (mouse CCKBC / human SERPINF1)
COLOR_VIP_OTHER = "#0072B2"  # blue (mouse VIP-other / human VIP-other)

# %% [markdown]
# ## Step 1: Load mouse data with labels

# %%
# Mouse M1 patchseq ephys (1,033 cells, custom feature extraction)
m1 = pd.read_csv(M1_EPHYS)
harmony = pd.read_csv(HARMONY_MAP)
harmony_labels = harmony[["Unnamed: 0", "is_cckbc", "mouse_RNA family"]].rename(
    columns={"Unnamed: 0": "cell id"}
)
m1 = m1.merge(harmony_labels, on="cell id", how="left")

mouse_cckbc = m1[m1["is_cckbc"].fillna(False)].copy()
mouse_vip = m1[
    (m1["mouse_RNA family"] == "Vip") & (~m1["is_cckbc"].fillna(False))
].copy()

print(f"Mouse CCKBC: {len(mouse_cckbc)}")
print(f"Mouse VIP-other: {len(mouse_vip)}")

# %% [markdown]
# ## Step 2: Load human data with ground truth labels

# %%
human_ephys = pd.read_csv(HUMAN_EPHYS)
human_meta = pd.read_csv(HUMAN_META)

# Merge ephys with metadata
human = human_ephys.merge(
    human_meta[["specimen_id", "Revised_subclass_label", "Transcriptomic_type "]],
    on="specimen_id",
    how="left",
)

# Select VIP cells only
human_vip_all = human[human["Revised_subclass_label"] == "VIP"].copy()

# Split: SERPINF1 types (CCKBC-like) vs other VIP
human_vip_all["is_serpinf1"] = human_vip_all["Transcriptomic_type "].str.contains(
    "SERPINF1", na=False
)
human_serpinf1 = human_vip_all[human_vip_all["is_serpinf1"]].copy()
human_vip_other = human_vip_all[~human_vip_all["is_serpinf1"]].copy()

print(f"Human VIP total: {len(human_vip_all)}")
print(f"Human VIP-SERPINF1 (CCKBC-like): {len(human_serpinf1)}")
print(f"Human VIP-other: {len(human_vip_other)}")
print(f"\nHuman VIP transcriptomic types:")
print(human_vip_all["Transcriptomic_type "].value_counts().to_string())

# %% [markdown]
# ## Step 3: Define matched feature pairs
#
# Features that discriminate mouse CCKBC from VIP-other (p<0.01),
# mapped to equivalent human IPFX features.

# %%
# (mouse_col, human_col, display_name)
FEATURE_PAIRS = [
    ("AP count",            "avg_rate_hero",           "Firing rate"),
    ("AP count 2nd half",   "adapt_hero",              "Sustained firing"),
    ("AP average amp adapt","upstroke_adapt_ratio",     "Waveform adaptation"),
    ("ISI CV",              "isi_cv_hero",             "ISI CV"),
    ("latency",             "latency_hero",            "First spike latency"),
    ("$V_{m}$ mean",        "v_baseline",              "Membrane potential"),
]

# Verify all columns exist
for m_col, h_col, name in FEATURE_PAIRS:
    assert m_col in m1.columns, f"Mouse column '{m_col}' not found"
    assert h_col in human_ephys.columns, f"Human column '{h_col}' not found"
print("All feature pairs verified.")

# %% [markdown]
# ## Step 4: Test discrimination in mouse and human

# %%
print("=" * 70)
print("MOUSE: CCKBC vs VIP-other")
print(f"{'Feature':25s} {'P-value':>12} {'CCKBC mean':>12} {'VIP mean':>12}")
print("-" * 70)
for m_col, h_col, name in FEATURE_PAIRS:
    cv = mouse_cckbc[m_col].dropna()
    vv = mouse_vip[m_col].dropna()
    if len(cv) >= 3 and len(vv) >= 3:
        u, p = mannwhitneyu(cv, vv, alternative="two-sided")
        print(f"  {name:25s} {p:12.2e} {cv.mean():12.3f} {vv.mean():12.3f}")

print()
print("=" * 70)
print("HUMAN: VIP-SERPINF1 (CCKBC-like) vs VIP-other")
print(f"{'Feature':25s} {'P-value':>12} {'SERP mean':>12} {'VIP mean':>12} {'n_serp':>7} {'n_vip':>7}")
print("-" * 70)
for m_col, h_col, name in FEATURE_PAIRS:
    sv = human_serpinf1[h_col].dropna()
    vv = human_vip_other[h_col].dropna()
    if len(sv) >= 2 and len(vv) >= 3:
        u, p = mannwhitneyu(sv, vv, alternative="two-sided")
        print(f"  {name:25s} {p:12.2e} {sv.mean():12.3f} {vv.mean():12.3f} {len(sv):>7} {len(vv):>7}")
    else:
        print(f"  {name:25s} {'too few':>12} {'-':>12} {'-':>12} {len(sv):>7} {len(vv):>7}")

# %% [markdown]
# ## Step 5: Build figure
#
# Each feature gets a subplot with 4 groups:
# Mouse CCKBC | Mouse VIP-other | Human VIP-other | Human SERPINF1
#
# Z-score within species so features are on comparable scales.
# Show p-values for CCKBC vs VIP-other (mouse) and SERPINF1 vs VIP-other (human).

# %%
def zscore_series(s):
    """Z-score a series, handling NaN."""
    return (s - s.mean()) / s.std()


def draw_violin(ax, data_list, colors, positions, labels):
    """Draw violin + strip + box for multiple groups."""
    rng = np.random.default_rng(42)

    # Filter out empty groups for violin
    valid_data = [(d, c, p) for d, c, p in zip(data_list, colors, positions) if len(d) >= 3]

    if valid_data:
        vdata, vcols, vpos = zip(*valid_data)
        parts = ax.violinplot(
            list(vdata), positions=list(vpos),
            showmeans=False, showmedians=False, showextrema=False
        )
        for body, col in zip(parts["bodies"], vcols):
            body.set_facecolor(col)
            body.set_alpha(0.35)
            body.set_edgecolor("none")

    # Box overlay + strip for all groups
    for vals, col, pos in zip(data_list, colors, positions):
        if len(vals) == 0:
            continue
        # Box (IQR)
        if len(vals) >= 3:
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.vlines(pos, q1, q3, color="k", lw=1.5, zorder=3)
            ax.scatter([pos], [med], color="white", edgecolor="k",
                       s=30, zorder=4, linewidths=0.8)
        # Strip
        if len(vals) < 10:
            # Few points: show as larger markers, no jitter
            ax.scatter(
                [pos] * len(vals), vals,
                color=col, alpha=0.8, s=40, edgecolors="k",
                linewidths=0.5, zorder=5, marker="D"
            )
        else:
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(
                pos + jitter, vals, color=col,
                alpha=0.15, s=8, edgecolors="none", zorder=2
            )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8)


def add_bracket(ax, x1, x2, y, p, color="k", lw=1.0):
    """Add significance bracket with p-value."""
    h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color=color, lw=lw)
    if p < 0.001:
        txt = f"p={p:.1e}"
    elif p < 0.05:
        txt = f"p={p:.3f}"
    else:
        txt = f"p={p:.2f}\nn.s."
    ax.text((x1 + x2) / 2, y + h, txt, ha="center", va="bottom",
            fontsize=6.5, color=color)


# %%
n_features = len(FEATURE_PAIRS)
fig, axes = plt.subplots(2, n_features, figsize=(3.0 * n_features, 8), dpi=150)

for col_idx, (m_col, h_col, name) in enumerate(FEATURE_PAIRS):
    # --- Top row: Mouse ---
    ax_m = axes[0, col_idx]

    mc = zscore_series(mouse_cckbc[m_col].dropna()).values
    mv = zscore_series(mouse_vip[m_col].dropna()).values
    # Z-score together (within mouse)
    all_mouse = pd.concat([mouse_cckbc[m_col], mouse_vip[m_col]]).dropna()
    mu, sd = all_mouse.mean(), all_mouse.std()
    mc_z = ((mouse_cckbc[m_col].dropna() - mu) / sd).values
    mv_z = ((mouse_vip[m_col].dropna() - mu) / sd).values

    draw_violin(
        ax_m,
        [mc_z, mv_z],
        [COLOR_CCKBC, COLOR_VIP_OTHER],
        [0, 1],
        ["CCKBC", "VIP-other"],
    )

    # P-value
    u, p_mouse = mannwhitneyu(mc_z, mv_z, alternative="two-sided")
    ylim = ax_m.get_ylim()
    bracket_y = ylim[1] + (ylim[1] - ylim[0]) * 0.05
    add_bracket(ax_m, 0, 1, bracket_y, p_mouse, color="k", lw=1.2)
    ax_m.set_ylim(ylim[0], bracket_y + (ylim[1] - ylim[0]) * 0.25)

    ax_m.set_title(name, fontsize=10, fontweight="bold")
    if col_idx == 0:
        ax_m.set_ylabel("z-scored value\n(within mouse)", fontsize=9)
    ax_m.spines["top"].set_visible(False)
    ax_m.spines["right"].set_visible(False)

    # --- Bottom row: Human ---
    ax_h = axes[1, col_idx]

    hs = human_serpinf1[h_col].dropna()
    hv = human_vip_other[h_col].dropna()
    # Z-score together (within human VIP)
    all_human = pd.concat([hs, hv])
    mu_h, sd_h = all_human.mean(), all_human.std()
    if sd_h > 0:
        hs_z = ((hs - mu_h) / sd_h).values
        hv_z = ((hv - mu_h) / sd_h).values
    else:
        hs_z = np.zeros(len(hs))
        hv_z = np.zeros(len(hv))

    draw_violin(
        ax_h,
        [hs_z, hv_z],
        [COLOR_CCKBC, COLOR_VIP_OTHER],
        [0, 1],
        [f"SERPINF1\n(n={len(hs_z)})", f"VIP-other\n(n={len(hv_z)})"],
    )

    # P-value
    if len(hs_z) >= 2 and len(hv_z) >= 3:
        u, p_human = mannwhitneyu(hs_z, hv_z, alternative="two-sided")
    else:
        p_human = 1.0
    ylim_h = ax_h.get_ylim()
    bracket_y_h = ylim_h[1] + (ylim_h[1] - ylim_h[0]) * 0.05
    add_bracket(ax_h, 0, 1, bracket_y_h, p_human, color="k", lw=1.2)
    ax_h.set_ylim(ylim_h[0], bracket_y_h + (ylim_h[1] - ylim_h[0]) * 0.25)

    if col_idx == 0:
        ax_h.set_ylabel("z-scored value\n(within human VIP)", fontsize=9)
    ax_h.spines["top"].set_visible(False)
    ax_h.spines["right"].set_visible(False)

# Row labels
axes[0, 0].annotate(
    "Mouse M1", xy=(-0.4, 0.5), xycoords="axes fraction",
    fontsize=12, fontweight="bold", rotation=90, va="center", ha="center",
)
axes[1, 0].annotate(
    "Human MTG", xy=(-0.4, 0.5), xycoords="axes fraction",
    fontsize=12, fontweight="bold", rotation=90, va="center", ha="center",
)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLOR_CCKBC, alpha=0.5, label=f"CCKBC / SERPINF1"),
    Patch(facecolor=COLOR_VIP_OTHER, alpha=0.5, label=f"VIP-other"),
]
fig.legend(
    handles=legend_elements, loc="upper center", ncol=2,
    fontsize=10, frameon=False, bbox_to_anchor=(0.5, 1.0)
)

fig.suptitle(
    "Mouse CCKBC-discriminating features do not transfer to human VIP subtypes",
    fontsize=12, fontweight="bold", y=1.03,
)
fig.tight_layout()
fig.patch.set_alpha(0)
for ax in axes.flat:
    ax.patch.set_alpha(0)

# Save
outpath = OUTDIR / "FigSX_ephys_cross_species_comparison.png"
fig.savefig(outpath, dpi=150, transparent=True, bbox_inches="tight")
outpath_pdf = OUTDIR / "FigSX_ephys_cross_species_comparison.pdf"
fig.savefig(outpath_pdf, transparent=True, bbox_inches="tight")
print(f"Saved: {outpath}")
print(f"Saved: {outpath_pdf}")

plt.show()
print("Done.")
