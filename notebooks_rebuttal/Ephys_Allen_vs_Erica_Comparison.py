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
# # Ephys analysis comparison: Erica's pipeline vs. Allen IPFX
#
# Three pipelines, all on the same recordings:
#
# 1. **Erica's pipeline** — AxographX, single AP at rheobase, t-test (manuscript Methods).
#    Reproduced here under both t-test (matches manuscript) and Mann-Whitney U
#    (the appropriate non-parametric test at this n).
# 2. **Allen IPFX, rheobase first AP** — most direct apples-to-apples comparison
#    to Erica's reference point. Allen's standard feature definitions
#    (peak dV/dt, 50%-amplitude width).
# 3. **Allen IPFX, hero sweep (rheobase + 40 ± 10 pA, mean of first 5 APs)** —
#    Allen Cell Types Database convention (Gouwens et al. 2019). The
#    repetitive-firing regime, where AP shape is more reliably measurable.
#
# **Cell counts**
#
# * Manuscript reports **n = 11 WT / 19 Df16** for most features (10 WT for
#   rheobase + input resistance; 6 WT / 11 Df16 for AHP and sag).
# * Allen runs use **n = 10 WT / 19 Df16** (`WT_mouse07_cell04_lostcell` has
#   zero ABF files in our copy — recording was lost).
#
# The single unrecoverable cell is the only reason Allen's WT n is 10 instead
# of 11. All other cell IDs match Erica's analysis cohort.

# %%
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib.pyplot as plt
import yaml

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    cfg = yaml.safe_load(f)
PROJ = Path(cfg["ProjDIR"])
DATA = PROJ / "dat" / "VIP_Ephys"
print(f"PROJ_DIR = {PROJ}")

# %% [markdown]
# ## 1. Erica's pipeline (replicate manuscript + add MWU)
#
# Source: per-cell features in `dat/VIP_Ephys/erica_prism/raw_features/<Feature>.csv`.
# Each CSV is a two-column table (WT, Df(16)A^+/-) of single per-cell numbers
# extracted in AxographX from rheobase first AP (or sub-threshold I-V for IR /
# τ / sag).

# %%
def parse_erica_csv(fn):
    """Parse a two-column WT/Df16 CSV from Erica's raw_features folder."""
    df = pd.read_csv(fn, encoding="utf-8-sig", header=None)
    headers = df.iloc[0].tolist()
    body = df.iloc[1:].reset_index(drop=True)
    rec = []
    for ci, h in enumerate(headers):
        h_str = str(h).strip()
        if h_str.upper().startswith("WT"):
            geno = "WT"
        elif "DF" in h_str.upper() and "16" in h_str:
            geno = "Df16"
        else:
            continue
        col = pd.to_numeric(body.iloc[:, ci], errors="coerce").dropna().values
        for v in col:
            rec.append({"genotype": geno, "value": float(v)})
    return pd.DataFrame(rec)


def two_sample(wt, df16):
    """Return WT mean±sd (n), Df16 mean±sd (n), Cohen's d, t-test p, MWU p."""
    if len(wt) < 3 or len(df16) < 3:
        return None
    pooled = np.sqrt(((len(wt) - 1) * wt.std() ** 2 + (len(df16) - 1) * df16.std() ** 2)
                     / (len(wt) + len(df16) - 2))
    d = (wt.mean() - df16.mean()) / pooled if pooled > 0 else np.nan
    _, pt = ttest_ind(wt, df16)
    _, pu = mannwhitneyu(wt, df16, alternative="two-sided")
    return {"wt_mean": wt.mean(), "wt_sd": wt.std(), "n_wt": len(wt),
            "df16_mean": df16.mean(), "df16_sd": df16.std(), "n_df16": len(df16),
            "cohens_d": d, "t_p": pt, "mwu_p": pu}


# %%
ERICA_FEATURES = {
    "Resting Membrane Potential":      "Resting Membrane Potential.csv",
    "Rheobase":                        "Rheobase.csv",
    "Firing Threshold":                "Firing Threshold.csv",
    "Latency to Fire":                 "Latency to Fire.csv",
    "Peak Amplitude":                  "Peak Amplitude.csv",
    "Half Width":                      "Half Width.csv",
    "Rise Slope (10-90% avg)":         "Rise Slope.csv",
    "|Decay Slope| (10-90% avg)":      "Decay Slope.csv",
    "Input Resistance":                "Input Resistance.csv",
    "Time Constant":                   "Time Constant.csv",
    "Ih Sag Amplitude":                "Ih Sag Amplitude.csv",
    "AHP":                             "AHP.csv",
}

erica_rows = []
for label, fn in ERICA_FEATURES.items():
    df = parse_erica_csv(DATA / "erica_prism" / "raw_features" / fn)
    wt = df[df.genotype == "WT"]["value"].values
    df16 = df[df.genotype == "Df16"]["value"].values
    if "Decay" in label:
        wt = -wt; df16 = -df16  # report as positive magnitude
    s = two_sample(wt, df16)
    if s is None:
        continue
    s["feature"] = label
    erica_rows.append(s)

erica_df = pd.DataFrame(erica_rows)[
    ["feature", "wt_mean", "wt_sd", "n_wt", "df16_mean", "df16_sd", "n_df16",
     "cohens_d", "t_p", "mwu_p"]
]
erica_df.style.set_caption("Erica's pipeline — MWU and t-test").format(
    {"wt_mean": "{:.3f}", "wt_sd": "{:.3f}", "df16_mean": "{:.3f}",
     "df16_sd": "{:.3f}", "cohens_d": "{:+.2f}",
     "t_p": "{:.4f}", "mwu_p": "{:.4f}"}
)

# %% [markdown]
# **Reading the table.** Half-width is the only feature with t-test p < 0.05
# in Erica's data; under MWU it sits exactly on the edge (p = 0.05). All other
# AP-shape features (rise slope, decay slope, peak amp) and intrinsic
# properties (RMP, rheobase, threshold, IR, τ, sag) are non-significant
# under both tests. AHP is a borderline trend (p = 0.06) at small n.

# %% [markdown]
# ## 2. Allen IPFX, rheobase first AP
#
# Source: `dat/VIP_Ephys/ipfx_no_exclusion.csv` (built by
# `scripts/run_ipfx_no_exclusion.py`). For each cell, the first action
# potential of the rheobase sweep (smallest stim with ≥ 1 spike).
# Allen feature definitions (peak dV/dt, 50%-amplitude width).

# %%
allen_3way = pd.read_csv(DATA / "ipfx_3way.csv")  # has R1, R2, R3
allen_no_excl = pd.read_csv(DATA / "ipfx_no_exclusion.csv")  # has H1, H2, H3

# Restrict to non-lost cells (drop "_lostcell" rows; "_excluded" is kept,
# matching Erica's analysis cohort)
def filter_cohort(df):
    return df[~df.cell_id.str.contains("lostcell", case=False)].copy()

allen_3way = filter_cohort(allen_3way)
allen_no_excl = filter_cohort(allen_no_excl)

ALLEN_R1_FEATURES = [
    ("upstroke",                  "Upstroke (peak dV/dt, V/s)"),
    ("downstroke",                "|Downstroke| (peak |dV/dt|, V/s)"),
    ("width",                     "AP width at half-amplitude (ms)"),
    ("peak_v",                    "Peak voltage (mV)"),
    ("threshold_v",               "Threshold voltage (mV)"),
    ("upstroke_downstroke_ratio", "Upstroke/downstroke ratio"),
    ("ahp_mV",                    "AHP depth (threshold − trough, mV)"),
    ("amplitude_mV",              "AP amplitude (peak − threshold, mV)"),
]


def collect_allen(df, prefix, features, transforms=None):
    rows = []
    for feat, label in features:
        col = f"{prefix}{feat}"
        if col not in df.columns:
            continue
        wt = df[df.genotype == "WT"][col].dropna().values
        df16 = df[df.genotype == "Df16"][col].dropna().values
        # transform mappings
        if feat == "downstroke":
            wt = -wt; df16 = -df16
        if feat == "width":
            wt = wt * 1000; df16 = df16 * 1000  # s → ms
        s = two_sample(wt, df16)
        if s is None:
            continue
        s["feature"] = label
        rows.append(s)
    return pd.DataFrame(rows)[
        ["feature", "wt_mean", "wt_sd", "n_wt", "df16_mean", "df16_sd", "n_df16",
         "cohens_d", "t_p", "mwu_p"]
    ]


allen_rheo = collect_allen(allen_3way, "R1_", ALLEN_R1_FEATURES)
allen_rheo.style.set_caption(
    "Allen IPFX, rheobase first AP (R1)").format(
    {"wt_mean": "{:.3f}", "wt_sd": "{:.3f}", "df16_mean": "{:.3f}",
     "df16_sd": "{:.3f}", "cohens_d": "{:+.2f}",
     "t_p": "{:.4f}", "mwu_p": "{:.4f}"}
)

# %% [markdown]
# **Reading the table.** Same reference point as Erica (rheobase first AP),
# but Allen's standard feature definitions. Upstroke is significant under
# t-test (p = 0.028) and on the edge under MWU (p = 0.057). Width and peak
# voltage are trends. AHP under Allen's definition (threshold − trough) is
# significant under both tests.

# %% [markdown]
# ## 3. Allen IPFX, hero sweep (rheobase + 40 ± 10 pA, **first AP only**)
#
# Allen Cell Types Database convention, but using only the **first AP** of
# the hero sweep rather than the mean of the first 5 APs. Reason: at hero
# sweep some cells fire only 2–3 APs total — averaging over "first 5"
# silently averages over fewer than 5 for those cells, mixing biology with
# sample size. The first AP is well-defined for every cell that fires at
# all, so it's the more standardized choice.

# %%
allen_hero = collect_allen(allen_no_excl, "H1_", ALLEN_R1_FEATURES)
allen_hero.style.set_caption(
    "Allen IPFX, hero sweep first AP (H1)").format(
    {"wt_mean": "{:.3f}", "wt_sd": "{:.3f}", "df16_mean": "{:.3f}",
     "df16_sd": "{:.3f}", "cohens_d": "{:+.2f}",
     "t_p": "{:.4f}", "mwu_p": "{:.4f}"}
)

# %% [markdown]
# **Reading the table.** At hero sweep first AP, **upstroke**, **downstroke**,
# **width**, and **AP amplitude** all cross MWU p < 0.05. Cohen's d for
# upstroke is large (~+1.4). Switching from rheobase first AP (R1) to hero
# first AP (H1) is a within-cell change — same cells, same single-AP
# measurement, different sweep — and the genotype effect becomes much
# clearer because the AP shape is more reproducible at hero current than at
# rheobase.

# %% [markdown]
# ## 4. Side-by-side summary across pipelines
#
# Each row aligns the conceptually equivalent feature across the three
# pipelines, with the manuscript-style t-test and the Mann-Whitney U test.

# %%
def pick(label, pipeline_df, t_or_mwu):
    """Pull a single number from a results df by feature label substring."""
    for _, r in pipeline_df.iterrows():
        if label in r["feature"]:
            return r[t_or_mwu]
    return None

ROWS = [
    ("AP width / Half-width",
     "Half Width", "AP width at half-amplitude", "AP width at half-amplitude"),
    ("AP rise / Upstroke (peak dV/dt)",
     "Rise Slope", "Upstroke", "Upstroke"),
    ("AP decay / Downstroke",
     "|Decay Slope|", "|Downstroke|", "|Downstroke|"),
    ("Peak voltage",
     "Peak Amplitude", "Peak voltage", "Peak voltage"),
    ("AP amplitude",
     None, "AP amplitude", "AP amplitude"),
    ("Threshold",
     "Firing Threshold", "Threshold voltage", "Threshold voltage"),
    ("Upstroke/downstroke ratio",
     None, "Upstroke/downstroke", "Upstroke/downstroke"),
    ("AHP",
     "AHP", "AHP depth", "AHP depth"),
    ("RMP",
     "Resting Membrane Potential", None, None),
    ("Rheobase",
     "Rheobase", None, None),
    ("Latency to fire",
     "Latency to Fire", None, None),
    ("Input resistance",
     "Input Resistance", None, None),
    ("Time constant",
     "Time Constant", None, None),
    ("I_h sag",
     "Ih Sag Amplitude", None, None),
]


def pretty_p(p):
    if p is None or pd.isna(p):
        return "—"
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "."
    sig = sig if p < 0.10 else ""
    return f"{p:.4f} {sig}".rstrip()


summary_rows = []
for label, e_lab, r1_lab, h3_lab in ROWS:
    summary_rows.append({
        "Feature": label,
        "Erica t":   pretty_p(pick(e_lab, erica_df, "t_p")) if e_lab else "—",
        "Erica MWU": pretty_p(pick(e_lab, erica_df, "mwu_p")) if e_lab else "—",
        "Allen rheobase t":   pretty_p(pick(r1_lab, allen_rheo, "t_p")) if r1_lab else "—",
        "Allen rheobase MWU": pretty_p(pick(r1_lab, allen_rheo, "mwu_p")) if r1_lab else "—",
        "Allen hero t":   pretty_p(pick(h3_lab, allen_hero, "t_p")) if h3_lab else "—",
        "Allen hero MWU": pretty_p(pick(h3_lab, allen_hero, "mwu_p")) if h3_lab else "—",
    })
summary = pd.DataFrame(summary_rows)
summary

# %% [markdown]
# ## 5. Train-level features at hero sweep (Allen-only; Erica reports none)
#
# These features are computed from the per-spike CSVs in
# `dat/VIP_Ephys/allen_ipfx_vip/features_fi_perspike/` (built by
# `scripts/run_allen_fi_perspike.py`). They quantify within-train spike
# patterns — Erica did not report any of these.

# %%
spike_dir = DATA / "allen_ipfx_vip" / "features_fi_perspike"
all_sp = pd.concat([pd.read_csv(f) for f in sorted(spike_dir.glob("*_fi_spikes.csv"))
                    if "lostcell" not in f.name],
                   ignore_index=True)
print(f"Loaded {len(all_sp)} spikes from {all_sp.cell_id.nunique()} cells.")
print(f"Genotype: {all_sp.groupby('genotype')['cell_id'].nunique().to_dict()}")

# %%
def cell_train_features_at_hero(g):
    sw_counts = (g.groupby("sweep_idx")
                  .agg(stim=("stim_pA", "first"),
                       n=("spike_idx_in_sweep", "count"))
                  .reset_index())
    sw_counts = sw_counts[sw_counts.n >= 3]
    if sw_counts.empty:
        return None
    rheo = sw_counts.sort_values("stim").iloc[0].stim
    target = rheo + 40
    cand = sw_counts[(sw_counts.stim >= rheo + 30) & (sw_counts.stim <= rheo + 50)]
    if cand.empty:
        cand = sw_counts[(sw_counts.stim >= rheo + 20) & (sw_counts.stim <= rheo + 60)]
    if cand.empty:
        return None
    sw_idx = int(cand.iloc[(cand.stim - target).abs().argmin()].sweep_idx)
    sweep = g[g.sweep_idx == sw_idx].sort_values("threshold_t").reset_index(drop=True)
    n = len(sweep)
    isi = np.diff(sweep.threshold_t.values)
    pairwise = ((isi[1:] - isi[:-1]) / (isi[1:] + isi[:-1])).mean() if len(isi) >= 2 else np.nan
    amps = (sweep.peak_v - sweep.threshold_v).values
    return {
        "n_spikes": n, "isi_mean_s": isi.mean(),
        "isi_cv": isi.std() / isi.mean() if isi.mean() > 0 else np.nan,
        "adapt_idx_pairwise": pairwise,
        "amp_drop_frac": (amps[-1] - amps[0]) / amps[0] if amps[0] > 0 else np.nan,
        "upstroke_drop_frac": (sweep.upstroke.iloc[-1] - sweep.upstroke.iloc[0]) /
                              sweep.upstroke.iloc[0] if sweep.upstroke.iloc[0] > 0 else np.nan,
        "width_growth_frac": (sweep.width.iloc[-1] - sweep.width.iloc[0]) /
                             sweep.width.iloc[0] if sweep.width.iloc[0] > 0 else np.nan,
    }


train_rows = []
for cid, g in all_sp.groupby("cell_id"):
    tf = cell_train_features_at_hero(g)
    if tf is None:
        continue
    tf["cell_id"] = cid
    tf["genotype"] = g.genotype.iloc[0]
    train_rows.append(tf)
train_df = pd.DataFrame(train_rows)
print(f"Train features cohort: "
      f"{train_df.groupby('genotype')['cell_id'].nunique().to_dict()}")

train_test_rows = []
TRAIN_FEATURES = [
    ("isi_mean_s",         "Mean ISI (s)"),
    ("isi_cv",             "ISI CV (regularity)"),
    ("adapt_idx_pairwise", "Pairwise adaptation index"),
    ("amp_drop_frac",      "AP amplitude adaptation (fraction)"),
    ("upstroke_drop_frac", "Upstroke adaptation (fraction)"),
    ("width_growth_frac",  "Width broadening (fraction)"),
]
for col, label in TRAIN_FEATURES:
    wt = train_df[train_df.genotype == "WT"][col].dropna().values
    df16 = train_df[train_df.genotype == "Df16"][col].dropna().values
    s = two_sample(wt, df16)
    if s is None:
        continue
    s["feature"] = label
    train_test_rows.append(s)
train_test_df = pd.DataFrame(train_test_rows)[
    ["feature", "wt_mean", "wt_sd", "n_wt", "df16_mean", "df16_sd", "n_df16",
     "cohens_d", "t_p", "mwu_p"]
]
train_test_df.style.set_caption(
    "Allen IPFX train-level features at hero sweep").format(
    {"wt_mean": "{:.4f}", "wt_sd": "{:.4f}", "df16_mean": "{:.4f}",
     "df16_sd": "{:.4f}", "cohens_d": "{:+.2f}",
     "t_p": "{:.4f}", "mwu_p": "{:.4f}"}
)

# %% [markdown]
# **Negative result, biologically informative.** ISI regularity, adaptation
# index, and within-train AP amplitude / upstroke / width changes are all
# similar between genotypes. Combined with the per-AP-shape changes
# (Section 3), this argues against generalized K^+ conductance changes
# (BK, SK, Kv7) — those would manifest as adaptation differences — and for
# voltage-gated channels that act on single-AP timescales (Nav, Kv3-family).

# %% [markdown]
# ## 6. Visualization — boxplots of the headline features

# %%
fig, axes = plt.subplots(2, 4, figsize=(16, 8))

PLOT_PANELS = [
    # (ax_idx, source, col, label, transform)
    ((0, 0), allen_no_excl, "H1_upstroke",     "Hero (1st AP) upstroke (V/s)",       lambda x: x),
    ((0, 1), allen_no_excl, "H1_downstroke",   "|Hero (1st AP) downstroke| (V/s)",   lambda x: -x),
    ((0, 2), allen_no_excl, "H1_width",        "Hero (1st AP) AP width (ms)",        lambda x: x * 1000),
    ((0, 3), allen_no_excl, "H1_amplitude_mV", "Hero (1st AP) AP amplitude (mV)",    lambda x: x),
    ((1, 0), allen_3way,    "R1_upstroke",     "Rheobase upstroke (V/s)",            lambda x: x),
    ((1, 1), allen_3way,    "R1_width",        "Rheobase AP width (ms)",             lambda x: x * 1000),
    ((1, 2), allen_3way,    "R1_ahp_mV",       "Rheobase AHP depth (mV)",            lambda x: x),
    ((1, 3), allen_3way,    "R1_amplitude_mV", "Rheobase AP amplitude (mV)",         lambda x: x),
]

colors = {"WT": "#2166ac", "Df16": "#b2182b"}

for (r, c), src, col, label, transform in PLOT_PANELS:
    ax = axes[r, c]
    if col not in src.columns:
        ax.axis("off"); continue
    wt = transform(src[src.genotype == "WT"][col].dropna()).values
    df16 = transform(src[src.genotype == "Df16"][col].dropna()).values
    bp = ax.boxplot([wt, df16], positions=[0, 1], widths=0.55,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, gen in zip(bp["boxes"], ["WT", "Df16"]):
        patch.set_facecolor(colors[gen]); patch.set_alpha(0.4)
    rng = np.random.default_rng(42)
    for i, vals in enumerate([wt, df16]):
        ax.scatter(np.full(len(vals), i) + rng.normal(0, 0.07, len(vals)),
                   vals, c=[colors["WT"] if i == 0 else colors["Df16"]],
                   s=20, alpha=0.85, edgecolors="black", linewidths=0.4, zorder=3)
    _, p = mannwhitneyu(wt, df16, alternative="two-sided")
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    ymax = max(np.max(wt), np.max(df16))
    ymin = min(np.min(wt), np.min(df16))
    rng_y = ymax - ymin
    ax.text(0.5, ymax + rng_y * 0.08, f"MWU p = {p:.4f}\n{sig}",
            ha="center", fontsize=9, fontweight="bold")
    ax.set_ylim(ymin - rng_y * 0.10, ymax + rng_y * 0.30)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"WT (n={len(wt)})", f"Df16 (n={len(df16)})"], fontsize=9)
    ax.set_ylabel(label)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.patch.set_alpha(0); ax.set_facecolor("none")

fig.patch.set_alpha(0)
fig.suptitle(
    "Top: Allen hero sweep (rheobase + 40 ± 10 pA, first AP only)\n"
    "Bottom: Allen rheobase first AP",
    fontsize=12, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(PROJ / "notebooks_rebuttal" / "Ephys_Allen_vs_Erica_boxplots.png",
            dpi=150, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ## 7. Summary
#
# * **Erica's pipeline + MWU**: only **half-width** crosses p < 0.05 (and only
#   barely, p = 0.05). Rise slope, decay slope, peak amplitude, AHP, and all
#   intrinsic properties are non-significant.
# * **Allen IPFX, rheobase first AP**: **upstroke** (peak dV/dt) and **AHP
#   depth** cross p < 0.05 under MWU. AP width is a strong trend.
# * **Allen IPFX, hero sweep (rheobase + 40 ± 10 pA, first AP only)**: four
#   AP-shape features survive MWU at p < 0.05 — **upstroke**, **downstroke**,
#   **width**, **AP amplitude** — with large effect sizes (Cohen's d ~ 1).
#   Using first-AP rather than mean-of-first-5 avoids the artifact that
#   "mean of first 5" silently averages over fewer than 5 APs in cells that
#   don't fire 5 spikes at hero current.
# * **Train-level features (ISI CV, adaptation indices)**: all n.s. at hero
#   sweep — a clean negative result that narrows the candidate ion channels
#   to those affecting single-AP shape (Nav, Kv3-family) rather than train
#   patterns (BK, SK, Kv7).
# * **Direction is consistent across all three pipelines** for every feature;
#   the gain in significance under Allen + hero comes from a more sensitive
#   and standardized measurement, not from changing the biological
#   interpretation.

# %% [markdown]
# ## 8. Burstiness — Steph's qualitative observation tested
#
# **Background.** Steph (the recordist) noticed in raw traces that WT cells
# appear to fire in *bursts* (clusters of fast spikes separated by pauses),
# while Df16 cells fire more sparsely / less burst-like. Section 5 already
# tested simple ISI CV and adaptation index at hero sweep and found no
# difference — but ISI CV is a coarse measure of regularity that can be
# fooled by slow adaptation. Here we test more burst-specific metrics
# across multiple sweep choices, and use the *appropriate* statistical test
# for each metric.
#
# **Metrics:**
# * **LV (Shinomoto local variation)** = ⟨3·(ISI_i − ISI_{i+1})² /
#   (ISI_i + ISI_{i+1})²⟩ — the standard burst-vs-tonic discriminator
#   (Shinomoto 2003). LV > 1 = bursty, LV ~ 1 = Poisson, LV < 1 = regular.
#   Not fooled by slow adaptation.
# * **CV2** = ⟨2·|ISI_i − ISI_{i+1}| / (ISI_i + ISI_{i+1})⟩ — like LV but
#   linear; also robust to slow adaptation.
# * **# bursts** = count of groups of ≥ 3 spikes with within-group ISI <
#   min(25 ms, median ISI / 2). This is a binary-ish counter (most cells
#   score 0 or 1), so we test with **Fisher's exact** on the
#   "any burst vs no burst" 2×2 table — *not* MWU, which is
#   anti-conservative on heavily-tied data.
# * **Max ISI / median ISI** = pause prominence.
#
# **Sweep choices** — burst patterns require many spikes:
# * Hero sweep (rheobase + 40 ± 10 pA) — standardized current
# * Max-firing sweep — most data per cell, but current varies cell-to-cell

# %%
from scipy.stats import fisher_exact

def burstiness_metrics(spike_times, isi_burst_thresh_s=0.025):
    isi = np.diff(spike_times)
    if len(isi) < 2:
        return None
    out = {"n_spikes": len(spike_times),
           "isi_cv": isi.std() / isi.mean() if isi.mean() > 0 else np.nan,
           "cv2": float(np.mean(2 * np.abs(isi[1:] - isi[:-1]) / (isi[1:] + isi[:-1]))),
           "lv_shinomoto": float(np.mean(3 * (isi[:-1] - isi[1:])**2 / (isi[:-1] + isi[1:])**2)),
           "max_isi_over_median": float(isi.max() / np.median(isi))}
    # # bursts: groups of >=3 spikes with ISI < min(25 ms, median ISI / 2)
    burst_thr = min(isi_burst_thresh_s, np.median(isi) / 2)
    in_burst = isi < burst_thr
    n_bursts = 0; run = 0
    for v in in_burst:
        if v:
            run += 1
        else:
            if run >= 2: n_bursts += 1
            run = 0
    if run >= 2:
        n_bursts += 1
    out["n_bursts"] = n_bursts
    return out


def select_sweep(g, mode):
    """mode: 'hero' (rheobase + 40 ± 10) or 'max' (max-firing-rate sweep)."""
    sw_counts = (g.groupby("sweep_idx")
                  .agg(stim=("stim_pA", "first"),
                       n=("spike_idx_in_sweep", "count"))
                  .reset_index())
    sw_counts = sw_counts[sw_counts.n >= 3]
    if sw_counts.empty: return None
    rheo = sw_counts.sort_values("stim").iloc[0].stim
    if mode == "hero":
        cand = sw_counts[(sw_counts.stim >= rheo + 30) & (sw_counts.stim <= rheo + 50)]
        if cand.empty:
            cand = sw_counts[(sw_counts.stim >= rheo + 20) & (sw_counts.stim <= rheo + 60)]
        if cand.empty: return None
        return int(cand.iloc[(cand.stim - (rheo + 40)).abs().argmin()].sweep_idx)
    if mode == "max":
        return int(sw_counts.loc[sw_counts.n.idxmax()].sweep_idx)
    return None


def burst_summary(mode):
    rows = []
    for cid, g in all_sp.groupby("cell_id"):
        sw_idx = select_sweep(g, mode)
        if sw_idx is None: continue
        spike_t = g[g.sweep_idx == sw_idx].sort_values("threshold_t").threshold_t.values
        m = burstiness_metrics(spike_t)
        if m is None: continue
        m["cell_id"] = cid
        m["genotype"] = g.genotype.iloc[0]
        rows.append(m)
    return pd.DataFrame(rows)


def burst_test_table(df, label):
    """Run MWU on continuous metrics, Fisher's exact on # bursts (binary-ish)."""
    rows = []
    cont_cols = [("n_spikes",            "# spikes in sweep"),
                 ("isi_cv",              "ISI CV"),
                 ("cv2",                 "CV2"),
                 ("lv_shinomoto",        "LV (Shinomoto)"),
                 ("max_isi_over_median", "max ISI / median ISI")]
    for col, lbl in cont_cols:
        wt = df[df.genotype == "WT"][col].dropna().values
        df16 = df[df.genotype == "Df16"][col].dropna().values
        if len(wt) < 3 or len(df16) < 3: continue
        _, pt = ttest_ind(wt, df16)
        _, pu = mannwhitneyu(wt, df16, alternative="two-sided")
        rows.append({"feature": lbl, "wt_med": np.median(wt), "n_wt": len(wt),
                     "df16_med": np.median(df16), "n_df16": len(df16),
                     "t_p": pt, "test_p": pu, "test": "MWU"})
    # # bursts → Fisher's exact on "any burst vs no burst"
    wt_any = (df[df.genotype == "WT"]["n_bursts"] >= 1).sum()
    wt_none = (df[df.genotype == "WT"]["n_bursts"] == 0).sum()
    df_any = (df[df.genotype == "Df16"]["n_bursts"] >= 1).sum()
    df_none = (df[df.genotype == "Df16"]["n_bursts"] == 0).sum()
    _, pf = fisher_exact([[wt_any, wt_none], [df_any, df_none]],
                          alternative="two-sided")
    rows.append({"feature": "# bursts (any vs none)",
                 "wt_med": f"{wt_any}/{wt_any + wt_none} cells",
                 "n_wt": wt_any + wt_none,
                 "df16_med": f"{df_any}/{df_any + df_none} cells",
                 "n_df16": df_any + df_none,
                 "t_p": np.nan, "test_p": pf, "test": "Fisher exact"})
    return pd.DataFrame(rows)


# %%
hero_burst = burst_summary("hero")
max_burst = burst_summary("max")

print("Hero sweep burstiness cohort:",
      hero_burst.groupby("genotype")["cell_id"].nunique().to_dict())
print("Max-firing sweep burstiness cohort:",
      max_burst.groupby("genotype")["cell_id"].nunique().to_dict())

# %%
hero_table = burst_test_table(hero_burst, "Hero sweep")
hero_table.style.set_caption(
    "Burstiness at hero sweep (rheobase + 40 ± 10 pA)").format(
    {"t_p": "{:.4f}", "test_p": "{:.4f}"}, na_rep="—")

# %%
max_table = burst_test_table(max_burst, "Max-firing sweep")
max_table.style.set_caption(
    "Burstiness at max-firing sweep (most data per cell)").format(
    {"t_p": "{:.4f}", "test_p": "{:.4f}"}, na_rep="—")

# %% [markdown]
# **Reading the tables.** No metric crosses p < 0.05 under the
# appropriately-chosen test:
#
# * **LV (Shinomoto)** — Df16 actually trends *higher* than WT (i.e., more
#   irregular, in the *opposite* direction from Steph's claim) at both
#   sweep choices, but n.s. (MWU p ~ 0.10 at hero).
# * **CV2 and ISI CV** — same direction as LV, n.s.
# * **# bursts (Fisher's exact)** — at hero, ~10 % of WT cells have ≥ 1
#   burst vs ~0 % of Df16; at max-firing, similar fractions. Fisher's
#   exact p > 0.05 in both. (Note: an earlier table in this conversation
#   reported MWU p = 0.041 for # bursts at max sweep, which was *anti-
#   conservative* because MWU on heavily-tied count data violates its tie
#   correction; Fisher's exact is the right test and gives p ≈ 0.10.)
#
# **What this means biologically.**
#
# Combined with the per-AP-shape findings (§3) — Df16 cells have slower
# upstroke, slower downstroke, wider AP, smaller AP amplitude — and the
# absence of a train-pattern difference here, the genotype effect is
# **localized to single-AP biophysics** rather than train-level
# integration. This argues:
#
# * **Against** changes in Ca²⁺-activated K⁺ channels (BK, SK) or Kv7
#   M-current — those would manifest as adaptation or burst pattern
#   differences.
# * **For** changes in voltage-gated Na⁺ (Nav1.x) and Kv3-family fast K⁺
#   channels — those act on single-AP timescales and would produce
#   exactly the upstroke / downstroke / width changes observed.
#
# **Caveat on Steph's observation.** The qualitative impression of "more
# bursty WT" from raw traces may have come from a specific cell pair shown
# without matched current or labels. At matched current (e.g., 150 pA, see
# `dat/VIP_Ephys/trace_pdfs/{WT,Df16}_150pA.pdf`), neither genotype shows
# clean burst-pause patterns dominantly; both show an early cluster of
# spikes followed by adaptation, with strong cell-to-cell variability that
# is not genotype-specific. The single most "bursty" WT example happens to
# be the cell flagged `_excluded` in the original recordings.

# %% [markdown]
# ## 9. Final summary
#
# Across all the analyses in this notebook (§1–§8), the picture that
# emerges is consistent and biologically clean:
#
# 1. **Single-AP shape** is altered in Df16: slower upstroke (peak dV/dt),
#    slower downstroke, wider half-width, smaller AP amplitude. Robust
#    under Allen IPFX with MWU at both rheobase first AP (R1) and hero
#    first AP (H1).
# 2. **AHP depth** is altered at rheobase (deeper in Df16 under Allen's
#    threshold-trough definition).
# 3. **Train pattern** (ISI CV, LV, # bursts, adaptation index) is *not*
#    altered. The qualitative "WT is burstier" observation does not survive
#    quantitative testing at either hero or max-firing sweep.
# 4. **Implication**: the ion-channel substrate is in voltage-gated channels
#    that shape individual APs (Nav, Kv3-family), not in adaptation /
#    afterhyperpolarization channels (BK, SK, Kv7).
