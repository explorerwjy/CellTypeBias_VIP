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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Allen IPFX ephys features — SVG export
#
# Plot every Allen-IPFX-extracted spike feature for the **first AP** at:
#
# - **R1 — rheobase** (smallest current that fires ≥ 1 AP)
# - **R2 — rheobase + 20 pA** (first AP of the next-up step)
#
# Source table: `dat/VIP_Ephys/ipfx_3way.csv` (built by `scripts/run_ipfx_3way.py`).
# Output: per-feature SVGs + combined grids in
# `notebooks_rebuttal/figures_allen_rheo_r20_svg/`.

# %%
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind

PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
DATA = PROJ / "dat" / "VIP_Ephys"
OUT_DIR = PROJ / "notebooks_rebuttal" / "figures_allen_rheo_r20_svg"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Match manuscript Figure 4 / S13 style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "svg.fonttype": "none",  # keep text editable in SVG
})

# Manuscript palette (Figure 4, Figure S13)
COL = {"WT": "#4a8a82", "Df16": "#b2557e"}
GENO_LABEL = {"WT": r"$\mathit{Wildtype}$",
              "Df16": r"$\mathit{Df(16)A^{+/-}}$"}

# %%
df = pd.read_csv(DATA / "ipfx_3way.csv")
df = df[~df.cell_id.str.contains("lostcell", case=False)].copy()
print(f"Cells loaded: {len(df)}  ({df.genotype.value_counts().to_dict()})")
print(f"R2 (rheobase + 20 pA) available: {df['R2_upstroke'].notna().sum()} / {len(df)}")

# %% [markdown]
# ## Feature definitions
#
# For each feature we apply small unit/sign transforms so plots use
# physiologist-friendly conventions: ms for AP width, positive
# downstroke magnitude, mV for AHP and amplitude.

# %%
# (col_in_csv, pretty_label, unit, transform)
FEATURES = [
    ("threshold_v",               "Threshold voltage",       "mV",     lambda x: x),
    ("peak_v",                    "Peak voltage",            "mV",     lambda x: x),
    ("upstroke",                  "Upstroke (peak dV/dt)",   "V/s",    lambda x: x),
    ("downstroke",                "|Downstroke| (peak)",     "V/s",    lambda x: -x),
    ("width",                     "AP width (half-amp.)",    "ms",     lambda x: x),
    ("upstroke_downstroke_ratio", "Upstroke / downstroke",   "ratio",  lambda x: x),
    ("trough_v",                  "Trough voltage",          "mV",     lambda x: x),
    ("fast_trough_v",             "Fast-trough voltage",     "mV",     lambda x: x),
    ("ahp_mV",                    "AHP depth",               "mV",     lambda x: x),
    ("amplitude_mV",              "AP amplitude",            "mV",     lambda x: x),
]

SWEEPS = [("R1", "Rheobase"),
          ("R2", "Rheobase + 20 pA")]

# %%
def stat_pair(wt, df16):
    wt = np.asarray(wt, dtype=float); wt = wt[~np.isnan(wt)]
    df16 = np.asarray(df16, dtype=float); df16 = df16[~np.isnan(df16)]
    t_p = mwu_p = np.nan
    if len(wt) >= 2 and len(df16) >= 2:
        try: t_p = ttest_ind(wt, df16, equal_var=False).pvalue
        except Exception: pass
        try: mwu_p = mannwhitneyu(wt, df16, alternative="two-sided").pvalue
        except Exception: pass
    return wt, df16, t_p, mwu_p


def p_to_stars(p):
    if pd.isna(p): return "n.s."
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."


# %% [markdown]
# ## Per-feature SVGs
#
# Each panel: WT vs Df16 box + jittered strip. Title shows MWU p-value
# (advisor's preferred test) with t-test p in parentheses. n's per group
# annotated under the x-axis.

# %%
def panel(ax, wt, df16, ylabel, title, t_p, mwu_p):
    """Manuscript-style box + jittered strip (matches Figure 4 / S13)."""
    rng = np.random.default_rng(0)
    box_data = [wt, df16]
    bp = ax.boxplot(box_data, positions=[0, 1], widths=0.55, showfliers=False,
                    patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.2),
                    whiskerprops=dict(color="black", linewidth=0.9),
                    capprops=dict(color="black", linewidth=0.9),
                    boxprops=dict(linewidth=0.9, edgecolor="black"))
    for patch, geno in zip(bp["boxes"], ("WT", "Df16")):
        patch.set_facecolor(COL[geno]); patch.set_alpha(0.55)
    for x, vals, geno in zip([0, 1], box_data, ("WT", "Df16")):
        jx = x + (rng.random(len(vals)) - 0.5) * 0.16
        ax.scatter(jx, vals, s=14, color=COL[geno],
                   edgecolor="black", linewidth=0.4, alpha=1.0, zorder=3)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([GENO_LABEL["WT"], GENO_LABEL["Df16"]], fontsize=9)
    ax.set_ylabel(ylabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=11)

    # Significance bracket above the boxes
    all_vals = np.concatenate([wt, df16])
    if all_vals.size > 0:
        ymax = float(np.nanmax(all_vals))
        ymin = float(np.nanmin(all_vals))
        yrange = ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)
        ytop = ymax + yrange * 0.10
        tick_h = yrange * 0.03
        sig = p_to_stars(mwu_p)
        ax.plot([0, 0, 1, 1],
                [ytop - tick_h, ytop, ytop, ytop - tick_h],
                color="black", linewidth=0.8, clip_on=False)
        ax.text(0.5, ytop + yrange * 0.02, sig, ha="center", va="bottom",
                fontsize=10 if sig != "n.s." else 9,
                fontweight="bold" if sig != "n.s." else "normal")
        ax.set_ylim(top=ytop + yrange * 0.18)

    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9, direction="out", length=3)
    ax.set_xlim(-0.6, 1.6)


# %%
# Per-feature: one SVG with R1 and R2 side-by-side
results_rows = []
for col, label, unit, tfm in FEATURES:
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6))
    for ax, (prefix, sweep_label) in zip(axes, SWEEPS):
        full_col = f"{prefix}_{col}"
        if full_col not in df.columns:
            ax.axis("off"); continue
        wt_raw = df[df.genotype == "WT"][full_col].values
        d16_raw = df[df.genotype == "Df16"][full_col].values
        wt = tfm(np.asarray(wt_raw, dtype=float))
        d16 = tfm(np.asarray(d16_raw, dtype=float))
        wt_v, d16_v, tp, mp = stat_pair(wt, d16)
        # Y-axis: just the unit (manuscript style); feature name goes in suptitle
        panel(ax, wt_v, d16_v, unit, sweep_label, tp, mp)
        results_rows.append({"feature": label, "sweep": sweep_label,
                             "wt_mean": wt_v.mean() if len(wt_v) else np.nan,
                             "wt_sd": wt_v.std(ddof=1) if len(wt_v) > 1 else np.nan,
                             "n_wt": len(wt_v),
                             "df16_mean": d16_v.mean() if len(d16_v) else np.nan,
                             "df16_sd": d16_v.std(ddof=1) if len(d16_v) > 1 else np.nan,
                             "n_df16": len(d16_v),
                             "t_p": tp, "mwu_p": mp})
    fig.suptitle(label, fontsize=12, fontweight="bold", y=1.02)
    fig.patch.set_alpha(0)
    for ax in axes: ax.patch.set_alpha(0)
    fig.tight_layout()
    safe = col.replace("/", "_")
    out_path = OUT_DIR / f"{safe}.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  saved {out_path.relative_to(PROJ)}")

results = pd.DataFrame(results_rows)
results.to_csv(OUT_DIR / "summary_stats.csv", index=False)
print(f"\nSummary table → {(OUT_DIR / 'summary_stats.csv').relative_to(PROJ)}")
results.style.format({"wt_mean": "{:.3f}", "wt_sd": "{:.3f}",
                       "df16_mean": "{:.3f}", "df16_sd": "{:.3f}",
                       "t_p": "{:.4f}", "mwu_p": "{:.4f}"})

# %% [markdown]
# ## Combined grids
#
# Two single-page SVGs — one for each sweep — laying out all features in
# a 2×5 grid for at-a-glance review.

# %%
def grid_svg(prefix, sweep_label, out_name):
    fig, axes = plt.subplots(2, 5, figsize=(16, 7))
    for ax, (col, label, unit, tfm) in zip(axes.flat, FEATURES):
        full_col = f"{prefix}_{col}"
        if full_col not in df.columns:
            ax.axis("off"); continue
        wt = tfm(df[df.genotype == "WT"][full_col].values.astype(float))
        d16 = tfm(df[df.genotype == "Df16"][full_col].values.astype(float))
        wt_v, d16_v, tp, mp = stat_pair(wt, d16)
        panel(ax, wt_v, d16_v, unit, "", tp, mp)
        ax.set_title(label, fontsize=11)
    fig.suptitle(f"Allen IPFX — {sweep_label} (1st AP)", fontsize=14,
                 fontweight="bold", y=1.00)
    fig.patch.set_alpha(0)
    for ax in axes.flat: ax.patch.set_alpha(0)
    fig.tight_layout()
    out_path = OUT_DIR / out_name
    fig.savefig(out_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  saved {out_path.relative_to(PROJ)}")


grid_svg("R1", "Rheobase",          "_grid_R1_rheobase.svg")
grid_svg("R2", "Rheobase + 20 pA",  "_grid_R2_rheo_plus20.svg")

# %% [markdown]
# ## Master raw-value table
#
# Three layouts in one xlsx so Steph can re-plot in whichever shape fits her
# workflow. Values use the same unit/sign transforms as the SVGs above
# (downstroke positive; AP width in ms; everything else native).
#
# - **Sheet `wide_cells_x_features`** — one row per cell, columns
#   `R1_<feature>` and `R2_<feature>`. Easiest to merge with other per-cell
#   tables.
# - **Sheet `long_tidy`** — `cell_id, genotype, feature, sweep, value`.
#   For ggplot/seaborn or group-by stats.
# - **One sheet per feature** — Prism-style: columns `WT_R1, Df16_R1,
#   WT_R2, Df16_R2`. Paste straight into Prism (column = group, row =
#   replicate).
#
# CSV mirror of the wide table is also written for non-Excel users.

# %%
def transformed_series(prefix, col, tfm):
    full = f"{prefix}_{col}"
    if full not in df.columns:
        return pd.Series(index=df.index, dtype=float)
    return tfm(df[full].astype(float))


# Wide: cells × (R1_<feat>, R2_<feat>)
wide_cols = {"cell_id": df["cell_id"].values,
             "genotype": df["genotype"].values,
             "rheobase_pA": df["rheobase_pA"].values,
             "R2_stim_pA": df["R2_stim_pA"].values}
for col, label, unit, tfm in FEATURES:
    for prefix, _ in SWEEPS:
        wide_cols[f"{prefix}_{col}"] = transformed_series(prefix, col, tfm).values
wide = pd.DataFrame(wide_cols)
wide.to_csv(OUT_DIR / "raw_values_wide.csv", index=False)
print(f"  wrote raw_values_wide.csv  ({wide.shape[0]} cells × {wide.shape[1]} cols)")

# Long: tidy
long_rows = []
for col, label, unit, tfm in FEATURES:
    for prefix, sweep_label in SWEEPS:
        s = transformed_series(prefix, col, tfm)
        for cid, geno, v in zip(df["cell_id"], df["genotype"], s):
            long_rows.append({"cell_id": cid, "genotype": geno,
                              "feature": label, "feature_col": col,
                              "unit": unit, "sweep": sweep_label,
                              "sweep_code": prefix, "value": v})
long_df = pd.DataFrame(long_rows)

# Prism-style per feature
def prism_block(col, tfm):
    block = {}
    for prefix, _ in SWEEPS:
        s = transformed_series(prefix, col, tfm)
        wt = s[df["genotype"] == "WT"].values
        d16 = s[df["genotype"] == "Df16"].values
        n = max(len(wt), len(d16))
        wt_pad = np.concatenate([wt, [np.nan] * (n - len(wt))])
        d16_pad = np.concatenate([d16, [np.nan] * (n - len(d16))])
        block[f"WT_{prefix}"] = wt_pad
        block[f"Df16_{prefix}"] = d16_pad
    return pd.DataFrame(block)


xlsx_path = OUT_DIR / "raw_values_master.xlsx"
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
    wide.to_excel(xw, sheet_name="wide_cells_x_features", index=False)
    long_df.to_excel(xw, sheet_name="long_tidy", index=False)
    for col, label, unit, tfm in FEATURES:
        # Excel sheet names ≤ 31 chars, no /, \, *, ?, :, [, ]
        sheet = col.replace("/", "_")[:31]
        prism_block(col, tfm).to_excel(xw, sheet_name=sheet, index=False)
print(f"  wrote {xlsx_path.relative_to(PROJ)}")
print(f"    sheets: wide_cells_x_features, long_tidy, + 1 per feature ({len(FEATURES)} total)")
wide.head()

# %% [markdown]
# ## Done
#
# **Figures**
# - Per-feature SVGs (R1 + R2 panels): `figures_allen_rheo_r20_svg/<feature>.svg`
# - Combined grids: `_grid_R1_rheobase.svg`, `_grid_R2_rheo_plus20.svg`
#
# **Tables**
# - Stats: `summary_stats.csv`
# - Raw values (3 layouts): `raw_values_master.xlsx`
# - Raw values wide (CSV mirror): `raw_values_wide.csv`
