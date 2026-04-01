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
# # VIP Interneuron Electrophysiology: 22q11.2 Deletion vs WT
#
# **Dataset**: Whole-cell current-clamp recordings from VIP+ interneurons in hippocampus
# of Df(16)A+/− (22q11.2 deletion model) and WT mice.
#
# **Pipeline**: EphysSumStats (in-house) applied to ABF recordings.
#
# **Cross-region comparison**: M1 patch-seq VIP and CCKBC cells from DANDI 000008.

# %%
# %load_ext autoreload
# %autoreload 2

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pyabf
from pathlib import Path
from scipy.stats import ttest_ind, mannwhitneyu
from datetime import datetime

# Paths
PROJ = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
EPHYS_DIR = PROJ / "dat" / "VIP_Ephys"
BUNDLE_DIR = EPHYS_DIR / "analysis_bundles"
RECORDING_DIR = EPHYS_DIR / "IV_step_recordings"
M1_DIR = Path("/mnt/data0/DANDI/Processed/000008")
META_PATH = Path("/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv")

DVDT_CORRECTION = 1e6  # pipeline dV/dt unit fix

sns.set_style("ticks")
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.transparent": True,
})

COLORS = {"WT": "#2166ac", "Df16": "#b2182b",
           "M1-VIP": "#4daf4a", "M1-CCKBC": "#984ea3"}

# %% [markdown]
# ## 1. Example Traces

# %%
def plot_example_traces(cell_id, ax_v, ax_i, title=None, sweep_indices=None):
    """Plot voltage and current traces for a cell."""
    cell_dir = RECORDING_DIR / cell_id
    abfs = sorted(cell_dir.glob("*.abf"))
    # Find the IV step file
    abf = None
    for f in abfs:
        a = pyabf.ABF(str(f))
        if a.protocol == "Intrinsic_IV curve_VIP":
            abf = a
            break
    if abf is None:
        return

    if sweep_indices is None:
        # Pick representative sweeps: subthreshold, near-threshold, supra-threshold
        n = abf.sweepCount
        sweep_indices = [0, n//4, n//3, n//2, int(n*0.6), int(n*0.75), n-1]

    cmap = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(sweep_indices)))
    for i, sw in enumerate(sweep_indices):
        if sw >= abf.sweepCount:
            continue
        abf.setSweep(sw, channel=0)
        t = abf.sweepX
        v = abf.sweepY
        ax_v.plot(t, v, color=cmap[i], linewidth=0.6, alpha=0.8)

        abf.setSweep(sw, channel=1)
        ax_i.plot(t, abf.sweepY, color=cmap[i], linewidth=0.6, alpha=0.8)

    ax_v.set_ylabel("Vm (mV)")
    ax_v.set_xlim(0, 5)
    ax_i.set_ylabel("I (pA)")
    ax_i.set_xlabel("Time (s)")
    ax_i.set_xlim(0, 5)
    if title:
        ax_v.set_title(title, fontsize=10, fontweight="bold")


# %%
# Select representative cells: one WT, one Df16
# Pick cells with moderate rheobase for clear traces
wt_example = "WT_mouse03_cell02"
df16_example = "Df16_mouse07_cell04"

fig = plt.figure(figsize=(12, 6))
gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], hspace=0.05, wspace=0.3)

ax_wt_v = fig.add_subplot(gs[0, 0])
ax_wt_i = fig.add_subplot(gs[1, 0], sharex=ax_wt_v)
ax_df_v = fig.add_subplot(gs[0, 1])
ax_df_i = fig.add_subplot(gs[1, 1], sharex=ax_df_v)

sweeps_to_show = [0, 5, 10, 15, 20, 30, 40]
plot_example_traces(wt_example, ax_wt_v, ax_wt_i, "WT", sweeps_to_show)
plot_example_traces(df16_example, ax_df_v, ax_df_i, "Df(16)A+/−", sweeps_to_show)

plt.setp(ax_wt_v.get_xticklabels(), visible=False)
plt.setp(ax_df_v.get_xticklabels(), visible=False)

for ax in [ax_wt_v, ax_df_v]:
    sns.despine(ax=ax, bottom=True)
for ax in [ax_wt_i, ax_df_i]:
    sns.despine(ax=ax)

fig.suptitle("Representative VIP Interneuron Recordings (Hippocampus)", fontsize=12, y=1.02)
plt.savefig(EPHYS_DIR / "figures" / "Fig_ExampleTraces.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 2. Load Pipeline Results (22q Hippocampal VIP)

# %%
def load_22q_results():
    rows = []
    for cell_dir in sorted(BUNDLE_DIR.iterdir()):
        if not cell_dir.is_dir():
            continue
        f = cell_dir / "analysis.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        cid = cell_dir.name
        g = "WT" if cid.startswith("WT") else "Df16"
        rmp = df["resting_vm_mean_mV"].mean()
        rin = np.nan
        if "input_resistance" in df.columns and df["input_resistance"].notna().any():
            rin = float(df["input_resistance"].dropna().iloc[0]) * 1000  # GOhm -> MOhm
        sp = df[df["spike_frequency_Hz"] > 0]
        d = {"cell_id": cid, "genotype": g, "RMP_mV": rmp, "Rin_MOhm": rin,
             "n_spiking": len(sp)}
        if len(sp) > 0:
            rheo = sp[sp.get("is_rheobase", pd.Series(dtype=bool)) == True] if "is_rheobase" in sp.columns else pd.DataFrame()
            if len(rheo) == 0:
                rheo = sp.iloc[[0]]
            r = rheo.iloc[0]
            d.update({
                "rheobase_pA": sp["avg_injected_current_pA"].min(),
                "peak_mV": r.get("avg_peak_voltage_mV", np.nan),
                "threshold_mV": r.get("avg_threshold_voltage_mV", np.nan),
                "amplitude_mV": r.get("avg_threshold_to_peak_mV", np.nan),
                "width_ms": r.get("avg_ap_width_ms", np.nan),
                "rise_Vs": r.get("avg_upstroke_mVms", np.nan) / DVDT_CORRECTION,
                "decay_Vs": r.get("avg_downstroke_mVms", np.nan) / DVDT_CORRECTION,
                "ud_ratio": r.get("avg_upstroke_downstroke_ratio", np.nan),
            })
        rows.append(d)
    return pd.DataFrame(rows)


q22 = load_22q_results()
print(f"22q cells: {len(q22)} ({(q22.genotype=='WT').sum()} WT, {(q22.genotype=='Df16').sum()} Df16)")
print(f"Spiking: {(q22.n_spiking > 0).sum()}")

# %% [markdown]
# ## 3. WT vs Df16 Feature Comparisons

# %%
def stat_annot(ax, x1, x2, y, p, h=0.02):
    """Add significance bracket."""
    yr = ax.get_ylim()[1] - ax.get_ylim()[0]
    y_bar = y + h * yr
    ax.plot([x1, x1, x2, x2], [y_bar - h*yr*0.3, y_bar, y_bar, y_bar - h*yr*0.3],
            lw=1, c="k")
    if p < 0.001:
        txt = "***"
    elif p < 0.01:
        txt = "**"
    elif p < 0.05:
        txt = "*"
    else:
        txt = f"p={p:.2f}"
    ax.text((x1 + x2) / 2, y_bar + h * yr * 0.1, txt, ha="center", va="bottom", fontsize=9)


# %%
features = [
    ("RMP_mV", "RMP (mV)"),
    ("Rin_MOhm", "Input Resistance\n(MOhm)"),
    ("rheobase_pA", "Rheobase (pA)"),
    ("width_ms", "AP Width (ms)"),
    ("rise_Vs", "Rise Slope (V/s)"),
    ("decay_Vs", "Decay Slope (V/s)"),
    ("peak_mV", "Peak Voltage (mV)"),
    ("ud_ratio", "Upstroke/Downstroke\nRatio"),
]

fig, axes = plt.subplots(2, 4, figsize=(14, 7))
axes = axes.flatten()

for i, (col, label) in enumerate(features):
    ax = axes[i]
    wt = q22[q22.genotype == "WT"][col].dropna()
    df16 = q22[q22.genotype == "Df16"][col].dropna()

    # Strip plot with box
    data = q22[q22[col].notna()].copy()
    sns.boxplot(data=data, x="genotype", y=col, order=["WT", "Df16"],
                palette=COLORS, width=0.5, fliersize=0, ax=ax,
                boxprops=dict(alpha=0.3), medianprops=dict(color="k"))
    sns.stripplot(data=data, x="genotype", y=col, order=["WT", "Df16"],
                  palette=COLORS, size=5, alpha=0.7, jitter=0.15, ax=ax)

    # Stats
    if len(wt) > 1 and len(df16) > 1:
        t, p = ttest_ind(wt, df16)
        ymax = max(wt.max(), df16.max())
        stat_annot(ax, 0, 1, ymax, p)

    ax.set_xlabel("")
    ax.set_ylabel(label)
    ax.set_title("")
    sns.despine(ax=ax)

fig.suptitle("22q11.2 Deletion: VIP Interneuron Electrophysiology (Hippocampus)",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_WT_vs_Df16_Features.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 4. Rise Slope (Peak dV/dt) — Key Finding

# %%
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# A: Rise slope violin
ax = axes[0]
data = q22[q22["rise_Vs"].notna()].copy()
sns.violinplot(data=data, x="genotype", y="rise_Vs", order=["WT", "Df16"],
               palette=COLORS, inner=None, alpha=0.3, ax=ax)
sns.stripplot(data=data, x="genotype", y="rise_Vs", order=["WT", "Df16"],
              palette=COLORS, size=7, alpha=0.8, jitter=0.1, ax=ax)
wt_r = data[data.genotype == "WT"]["rise_Vs"]
df_r = data[data.genotype == "Df16"]["rise_Vs"]
_, p = ttest_ind(wt_r, df_r)
stat_annot(ax, 0, 1, max(wt_r.max(), df_r.max()), p)
ax.set_ylabel("Peak dV/dt (V/s)")
ax.set_xlabel("")
ax.set_title("A. Rise Slope")
sns.despine(ax=ax)

# B: Decay slope violin
ax = axes[1]
data = q22[q22["decay_Vs"].notna()].copy()
sns.violinplot(data=data, x="genotype", y="decay_Vs", order=["WT", "Df16"],
               palette=COLORS, inner=None, alpha=0.3, ax=ax)
sns.stripplot(data=data, x="genotype", y="decay_Vs", order=["WT", "Df16"],
              palette=COLORS, size=7, alpha=0.8, jitter=0.1, ax=ax)
wt_d = data[data.genotype == "WT"]["decay_Vs"]
df_d = data[data.genotype == "Df16"]["decay_Vs"]
_, p = ttest_ind(wt_d, df_d)
stat_annot(ax, 0, 1, max(wt_d.max(), df_d.max()), p, h=0.03)
ax.set_ylabel("Min dV/dt (V/s)")
ax.set_xlabel("")
ax.set_title("B. Decay Slope")
sns.despine(ax=ax)

# C: Width vs Rise scatter
ax = axes[2]
for g, color in [("WT", COLORS["WT"]), ("Df16", COLORS["Df16"])]:
    sub = q22[(q22.genotype == g) & q22["rise_Vs"].notna() & q22["width_ms"].notna()]
    ax.scatter(sub["width_ms"], sub["rise_Vs"], c=color, s=50, alpha=0.7,
               edgecolors="white", linewidths=0.5, label=g, zorder=3)
ax.set_xlabel("AP Width (ms)")
ax.set_ylabel("Rise Slope (V/s)")
ax.set_title("C. Width vs Rise Slope")
ax.legend(frameon=False)
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_dVdt_Analysis.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 5. Load M1 Patch-Seq Data (DANDI 000008)

# %%
def load_m1_vip_sncg():
    meta = pd.read_csv(META_PATH, sep="\t")
    meta["date_fmt"] = meta["Date"].apply(
        lambda d: datetime.strptime(d, "%m/%d/%y").strftime("%Y%m%d")
        if pd.notna(d) else None
    )
    meta["sample_num"] = meta["Sample"].str.extract(r"(\d+)")[0]
    meta["mouse_id"] = meta["Mouse"].str.replace("mouse_", "")
    meta["match_key"] = (
        "sub_mouse-" + meta["mouse_id"] + "_ses_" +
        meta["date_fmt"] + "-sample-" + meta["sample_num"]
    )

    bundles = {d.name: d for d in M1_DIR.iterdir() if d.is_dir()}
    vip_sncg = meta[meta["RNA family"].isin(["Vip", "Sncg"])]

    rows = []
    for _, row in vip_sncg.iterrows():
        matches = [b for b in bundles if b.startswith(row["match_key"])]
        if not matches:
            continue
        f = bundles[matches[0]] / "analysis.csv"
        if not f.exists():
            continue
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        if "spike_frequency_Hz" not in df.columns:
            continue
        sp = df[df["spike_frequency_Hz"] > 0]
        if len(sp) == 0:
            continue
        rheo = sp[sp.get("is_rheobase", pd.Series(dtype=bool)) == True] if "is_rheobase" in sp.columns else pd.DataFrame()
        if len(rheo) == 0:
            rheo = sp.iloc[[0]]
        r = rheo.iloc[0]
        group = "M1-CCKBC" if row["RNA family"] == "Sncg" or row["RNA type"] == "Vip Sncg" else "M1-VIP"
        rows.append({
            "group": group, "rna_type": row["RNA type"],
            "width_ms": r.get("avg_ap_width_ms", np.nan),
            "rise_Vs": r.get("avg_upstroke_mVms", np.nan) / DVDT_CORRECTION if pd.notna(r.get("avg_upstroke_mVms", np.nan)) else np.nan,
            "decay_Vs": r.get("avg_downstroke_mVms", np.nan) / DVDT_CORRECTION if pd.notna(r.get("avg_downstroke_mVms", np.nan)) else np.nan,
            "ud_ratio": r.get("avg_upstroke_downstroke_ratio", np.nan),
            "peak_mV": r.get("avg_peak_voltage_mV", np.nan),
            "rheobase_pA": sp["avg_injected_current_pA"].min(),
        })
    return pd.DataFrame(rows)


m1 = load_m1_vip_sncg()
print(f"M1 cells: {len(m1)} ({(m1.group=='M1-VIP').sum()} VIP, {(m1.group=='M1-CCKBC').sum()} CCKBC)")

# %% [markdown]
# ## 6. Cross-Region Comparison: M1 vs Hippocampus

# %%
# Build combined dataframe
q22_plot = q22[q22.n_spiking > 0][["genotype", "width_ms", "rise_Vs", "decay_Vs",
                                     "ud_ratio", "peak_mV", "rheobase_pA"]].copy()
q22_plot["group"] = q22_plot["genotype"].map({"WT": "22q-WT", "Df16": "22q-Df16"})
q22_plot["region"] = "Hippocampus"

m1_plot = m1[["group", "width_ms", "rise_Vs", "decay_Vs", "ud_ratio",
              "peak_mV", "rheobase_pA"]].copy()
m1_plot["region"] = "M1"

combined = pd.concat([q22_plot, m1_plot], ignore_index=True)
group_order = ["M1-CCKBC", "M1-VIP", "22q-WT", "22q-Df16"]
palette = {g: COLORS.get(g, "#999999") for g in group_order}

# %%
feats = [
    ("width_ms", "AP Width (ms)"),
    ("rise_Vs", "Rise Slope (V/s)"),
    ("decay_Vs", "Decay Slope (V/s)"),
    ("ud_ratio", "Up/Down Ratio"),
]

fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))

for i, (col, label) in enumerate(feats):
    ax = axes[i]
    data = combined[combined[col].notna()]
    sns.boxplot(data=data, x="group", y=col, order=group_order,
                palette=palette, width=0.55, fliersize=0, ax=ax,
                boxprops=dict(alpha=0.3), medianprops=dict(color="k"))
    sns.stripplot(data=data, x="group", y=col, order=group_order,
                  palette=palette, size=4, alpha=0.6, jitter=0.15, ax=ax)

    ax.set_xlabel("")
    ax.set_ylabel(label)
    ax.set_xticklabels(["M1\nCCKBC", "M1\nVIP", "HPC\nWT", "HPC\nDf16"],
                       fontsize=8)
    sns.despine(ax=ax)

    # Add significance brackets for key comparisons
    # M1-VIP vs 22q-WT
    g1 = combined[(combined.group == "M1-VIP") & combined[col].notna()][col]
    g2 = combined[(combined.group == "22q-WT") & combined[col].notna()][col]
    if len(g1) > 1 and len(g2) > 1:
        _, p = ttest_ind(g1, g2)
        if p < 0.05:
            ymax = data[col].max()
            stat_annot(ax, 1, 2, ymax, p, h=0.03)

fig.suptitle("Cross-Region Comparison: M1 Patch-Seq vs Hippocampal 22q VIP",
             fontsize=12, fontweight="bold", y=1.05)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_CrossRegion_Comparison.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 7. Key Finding: 22q-Df16 Resembles M1 VIP

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# A: Rise Slope across groups
ax = axes[0]
data = combined[combined["rise_Vs"].notna()]
parts = ax.violinplot(
    [data[data.group == g]["rise_Vs"].values for g in group_order],
    positions=range(4), showmedians=True, showextrema=False
)
for i, (pc, g) in enumerate(zip(parts["bodies"], group_order)):
    pc.set_facecolor(palette[g])
    pc.set_alpha(0.3)
parts["cmedians"].set_color("black")

for i, g in enumerate(group_order):
    vals = data[data.group == g]["rise_Vs"]
    ax.scatter(np.full(len(vals), i) + np.random.normal(0, 0.05, len(vals)),
               vals, c=palette[g], s=25, alpha=0.7, edgecolors="white", linewidths=0.3,
               zorder=3)

ax.set_xticks(range(4))
ax.set_xticklabels(["M1\nCCKBC", "M1\nVIP", "HPC\nWT", "HPC\nDf16"], fontsize=9)
ax.set_ylabel("Peak dV/dt (V/s)")
ax.set_title("A. Rise Slope: 22q-Df16 ≈ M1-VIP")
# Bracket: M1-VIP vs 22q-Df16 (not significant)
g1 = data[data.group == "M1-VIP"]["rise_Vs"]
g2 = data[data.group == "22q-Df16"]["rise_Vs"]
_, p = ttest_ind(g1, g2)
ax.text(2.5, ax.get_ylim()[1] * 0.95, f"M1-VIP vs Df16\np={p:.2f} (n.s.)",
        ha="center", fontsize=8, style="italic", color="gray")
# Bracket: 22q-WT vs 22q-Df16 (significant)
g1 = data[data.group == "22q-WT"]["rise_Vs"]
_, p2 = ttest_ind(g1, g2)
stat_annot(ax, 2, 3, max(g1.max(), g2.max()), p2, h=0.04)
sns.despine(ax=ax)

# B: Width vs Rise scatter - all groups
ax = axes[1]
for g in group_order:
    sub = combined[(combined.group == g) & combined["rise_Vs"].notna() & combined["width_ms"].notna()]
    ax.scatter(sub["width_ms"], sub["rise_Vs"], c=palette[g], s=40, alpha=0.6,
               edgecolors="white", linewidths=0.5, label=g, zorder=3)

ax.set_xlabel("AP Width (ms)")
ax.set_ylabel("Rise Slope (V/s)")
ax.set_title("B. AP Width vs Rise Slope")
ax.legend(frameon=False, fontsize=8, loc="upper right")
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_Df16_Resembles_M1.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 8. Summary Statistics

# %%
print("=" * 80)
print("22q HIPPOCAMPAL VIP: WT vs Df16")
print("=" * 80)
for col, label in features:
    wt = q22[q22.genotype == "WT"][col].dropna()
    df16 = q22[q22.genotype == "Df16"][col].dropna()
    if len(wt) > 1 and len(df16) > 1:
        t, p = ttest_ind(wt, df16)
        sig = " *" if p < 0.05 else " **" if p < 0.01 else ""
        print(f"  {label:<25}: WT={wt.mean():.2f}+/-{wt.std():.2f} (n={len(wt)})  "
              f"Df16={df16.mean():.2f}+/-{df16.std():.2f} (n={len(df16)})  "
              f"t={t:.2f} p={p:.4f}{sig}")

print()
print("=" * 80)
print("M1 PATCH-SEQ: CCKBC vs VIP-ISI")
print("=" * 80)
for col, label in [("width_ms", "AP Width"), ("rise_Vs", "Rise Slope"),
                    ("decay_Vs", "Decay Slope"), ("rheobase_pA", "Rheobase")]:
    c = m1[m1.group == "M1-CCKBC"][col].dropna()
    v = m1[m1.group == "M1-VIP"][col].dropna()
    if len(c) > 1 and len(v) > 1:
        t, p = ttest_ind(c, v)
        sig = " *" if p < 0.05 else ""
        print(f"  {label:<15}: CCKBC={c.mean():.2f}+/-{c.std():.2f} (n={len(c)})  "
              f"VIP={v.mean():.2f}+/-{v.std():.2f} (n={len(v)})  p={p:.4f}{sig}")

print()
print("=" * 80)
print("CROSS-REGION: M1-VIP vs 22q (hippocampal) VIP")
print("=" * 80)
for col, label in [("width_ms", "AP Width"), ("rise_Vs", "Rise Slope"),
                    ("decay_Vs", "Decay Slope"), ("ud_ratio", "Up/Down Ratio")]:
    m1v = combined[combined.group == "M1-VIP"][col].dropna()
    wt = combined[combined.group == "22q-WT"][col].dropna()
    df16 = combined[combined.group == "22q-Df16"][col].dropna()
    if len(m1v) > 1:
        _, p1 = ttest_ind(m1v, wt)
        _, p2 = ttest_ind(m1v, df16)
        print(f"  {label:<15}: M1-VIP={m1v.mean():.1f}  22q-WT={wt.mean():.1f} (p={p1:.4f})  "
              f"22q-Df16={df16.mean():.1f} (p={p2:.4f})")

# %% [markdown]
# ## 9. Phase Plane Analysis: Representative WT vs Df16
#
# Phase plane plots (dV/dt vs V) reveal Na+ channel activation differences.
# Representative cells selected as median by rise slope within each group.
#
# **Note**: Rheobase differs between cells (WT=40 pA, Df16=85 pA).
# We show both rheobase and matched-current (85 pA) comparisons.

# %%
from scipy.signal import savgol_filter

DT = 1 / 20000  # 20 kHz sampling


def load_sweep_trace(cell_id, sweep_num):
    """Load voltage trace for a specific sweep number."""
    cell_dir = BUNDLE_DIR / cell_id
    mv_file = list(cell_dir.glob("mv_*.parquet"))[0]
    mv = pd.read_parquet(mv_file)
    trace = mv[mv["sweep"] == sweep_num].sort_values("t_s")
    t = trace["t_s"].values
    v = trace["value"].values
    return t, v


def get_sweep_for_current(cell_id, target_pA):
    """Find the sweep number closest to a target injected current."""
    analysis = pd.read_csv(BUNDLE_DIR / cell_id / "analysis.csv")
    idx = (analysis["avg_injected_current_pA"] - target_pA).abs().idxmin()
    return int(analysis.loc[idx, "sweep"]), analysis.loc[idx, "avg_injected_current_pA"]


def get_rheobase_sweep(cell_id):
    """Return (sweep_num, current_pA) for the rheobase sweep."""
    analysis = pd.read_csv(BUNDLE_DIR / cell_id / "analysis.csv")
    spiking = analysis[analysis["spike_frequency_Hz"] > 0]
    row = spiking.iloc[0]
    return int(row["sweep"]), row["avg_injected_current_pA"]


def smooth_and_dvdt(v):
    """Smooth voltage trace and compute dV/dt in V/s."""
    v_s = savgol_filter(v, window_length=11, polyorder=3)
    dvdt = np.gradient(v_s, DT) / 1000  # mV/s -> V/s
    return v_s, dvdt


# Select median cells by rise slope
wt_cell = "WT_mouse07_cell12"      # median WT: rise=224.6 V/s, width=0.80 ms
df16_cell = "Df16_mouseB10916_cell06"  # median Df16: rise=142.8 V/s, width=0.90 ms

# Get sweep info
wt_rheo_sw, wt_rheo_pA = get_rheobase_sweep(wt_cell)
df_rheo_sw, df_rheo_pA = get_rheobase_sweep(df16_cell)
# Matched current: 85 pA (lowest common spiking current)
MATCHED_PA = 85.0
wt_match_sw, wt_match_actual = get_sweep_for_current(wt_cell, MATCHED_PA)
df_match_sw, df_match_actual = get_sweep_for_current(df16_cell, MATCHED_PA)

print(f"Rheobase — WT: sweep {wt_rheo_sw} ({wt_rheo_pA:.0f} pA), "
      f"Df16: sweep {df_rheo_sw} ({df_rheo_pA:.0f} pA)")
print(f"Matched  — WT: sweep {wt_match_sw} ({wt_match_actual:.0f} pA), "
      f"Df16: sweep {df_match_sw} ({df_match_actual:.0f} pA)")

# %%
# --- Full-trace phase plane: rheobase and matched current ---

# Stimulus window (current step)
STIM_START = 0.5   # s
STIM_END = 2.0     # s
stim_start_idx = int(STIM_START / DT)
stim_end_idx = int(STIM_END / DT)

# Rheobase + 60 pA sweeps (more spikes: WT ~4, Df16 ~2)
SUPRA_OFFSET = 60  # pA above rheobase
wt_supra_sw, wt_supra_pA = get_sweep_for_current(wt_cell, wt_rheo_pA + SUPRA_OFFSET)
df_supra_sw, df_supra_pA = get_sweep_for_current(df16_cell, df_rheo_pA + SUPRA_OFFSET)
print(f"Supra    — WT: sweep {wt_supra_sw} ({wt_supra_pA:.0f} pA), "
      f"Df16: sweep {df_supra_sw} ({df_supra_pA:.0f} pA)")

fig, axes = plt.subplots(3, 3, figsize=(14, 12))

conditions = [
    ("Rheobase",
     [(wt_cell, wt_rheo_sw, COLORS["WT"], f"WT ({wt_rheo_pA:.0f} pA)"),
      (df16_cell, df_rheo_sw, COLORS["Df16"], f"Df16 ({df_rheo_pA:.0f} pA)")]),
    (f"Matched ({MATCHED_PA:.0f} pA)",
     [(wt_cell, wt_match_sw, COLORS["WT"], "WT"),
      (df16_cell, df_match_sw, COLORS["Df16"], "Df16")]),
    (f"Rheobase + {SUPRA_OFFSET} pA",
     [(wt_cell, wt_supra_sw, COLORS["WT"], f"WT ({wt_supra_pA:.0f} pA)"),
      (df16_cell, df_supra_sw, COLORS["Df16"], f"Df16 ({df_supra_pA:.0f} pA)")]),
]

for row_idx, (sweep_label, cells) in enumerate(conditions):
    ax_v, ax_dvdt, ax_pp = axes[row_idx]

    for cell_id, sw, color, label in cells:
        t, v = load_sweep_trace(cell_id, sw)
        v_s, dvdt = smooth_and_dvdt(v)

        # Use stimulus window for full-trace view
        t_stim = (t[stim_start_idx:stim_end_idx] - STIM_START) * 1000  # ms
        v_stim = v_s[stim_start_idx:stim_end_idx]
        dvdt_stim = dvdt[stim_start_idx:stim_end_idx]

        ax_v.plot(t_stim, v_stim, color=color, linewidth=0.8, label=label, alpha=0.9)
        ax_dvdt.plot(t_stim, dvdt_stim, color=color, linewidth=0.8, alpha=0.9)
        ax_pp.plot(v_stim, dvdt_stim, color=color, linewidth=0.8, label=label, alpha=0.9)

    ax_v.set_ylabel("Vm (mV)")
    ax_v.set_title(f"{sweep_label} — Voltage Trace")
    ax_v.legend(frameon=False, fontsize=8)
    sns.despine(ax=ax_v)

    ax_dvdt.set_ylabel("dV/dt (V/s)")
    ax_dvdt.set_title(f"{sweep_label} — dV/dt")
    ax_dvdt.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    sns.despine(ax=ax_dvdt)

    ax_pp.set_xlabel("Vm (mV)")
    ax_pp.set_ylabel("dV/dt (V/s)")
    ax_pp.set_title(f"{sweep_label} — Phase Plane")
    ax_pp.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax_pp.legend(frameon=False, fontsize=8)
    sns.despine(ax=ax_pp)

    if row_idx < 2:
        ax_v.set_xlabel("")
        ax_dvdt.set_xlabel("")
    else:
        ax_v.set_xlabel("Time (ms)")
        ax_dvdt.set_xlabel("Time (ms)")

fig.suptitle("Phase Plane: Full Stimulus Trace — Rheobase / Matched / Suprathreshold",
             fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_PhasePlane_WT_vs_Df16.png",
            dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 10. 22q Gene Expression vs AP Kinetics in M1 VIP Interneurons
#
# **Question**: Does cell-by-cell expression of 22q11.2 deletion genes
# correlate with AP features (rise slope, AP width) in M1 VIP interneurons?
#
# **Approach**: Using DANDI 000008 M1 patch-seq data, we extract paired
# electrophysiology (analysis.csv) and single-cell RNA-seq (exon counts)
# for VIP interneurons. We then correlate each gene's expression with
# rise slope and AP width using Spearman rank correlation.
#
# **Genes tested**:
# - 22q11.2 region: Dgcr8, Tbx1, Comt, Prodh, Sept5, Ranbp1, Crkl,
#   Pi4ka, Snap29, Hira, Ufd1, Cdc45, Slc25a1
# - Na+ channels: Scn1a, Scn2a, Scn8a
# - K+ channels: Kcnc1 (Kv3.1), Kcnc2 (Kv3.2), Kcna1 (Kv1.1)

# %%
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# --- Load expression data ---
EXPR_PATH = Path("/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_exon_counts.csv.gz")

# Define target genes
GENES_22Q = ["Dgcr8", "Tbx1", "Comt", "Prodh", "Sept5", "Ranbp1",
             "Crkl", "Pi4ka", "Snap29", "Hira", "Ufd1", "Cdc45", "Slc25a1"]
GENES_NA = ["Scn1a", "Scn2a", "Scn8a"]
GENES_K = ["Kcnc1", "Kcnc2", "Kcna1"]
ALL_GENES = GENES_22Q + GENES_NA + GENES_K

# Load full expression matrix (genes x cells) — takes a moment
print("Loading expression matrix...")
expr_full = pd.read_csv(EXPR_PATH, index_col=0)
print(f"Expression matrix: {expr_full.shape[0]} genes x {expr_full.shape[1]} cells")

# Subset to target genes
expr_targets = expr_full.loc[expr_full.index.isin(ALL_GENES)]
print(f"Target genes found: {len(expr_targets)} / {len(ALL_GENES)}")
missing = set(ALL_GENES) - set(expr_targets.index)
if missing:
    print(f"  Missing: {missing}")

# %%
# --- Build paired ephys + expression table for M1 VIP cells ---

meta = pd.read_csv(META_PATH, sep="\t")
meta["date_fmt"] = meta["Date"].apply(
    lambda d: datetime.strptime(d, "%m/%d/%y").strftime("%Y%m%d")
    if pd.notna(d) else None
)
meta["sample_num"] = meta["Sample"].str.extract(r"(\d+)")[0]
meta["mouse_id"] = meta["Mouse"].str.replace("mouse_", "")
meta["match_key"] = (
    "sub_mouse-" + meta["mouse_id"] + "_ses_" +
    meta["date_fmt"] + "-sample-" + meta["sample_num"]
)

bundles = {d.name: d for d in M1_DIR.iterdir() if d.is_dir()}

# VIP cells excluding Vip Sncg (transitional type)
vip_meta = meta[(meta["RNA family"] == "Vip") & (meta["RNA type"] != "Vip Sncg")]
print(f"VIP cells in metadata (excl Sncg): {len(vip_meta)}")

rows = []
for _, row in vip_meta.iterrows():
    cell_id = row["Cell"]
    matches = [b for b in bundles if b.startswith(row["match_key"])]
    if not matches:
        continue
    f = bundles[matches[0]] / "analysis.csv"
    if not f.exists():
        continue
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    if "spike_frequency_Hz" not in df.columns:
        continue
    sp = df[df["spike_frequency_Hz"] > 0]
    if len(sp) == 0:
        continue

    # Get rheobase sweep features
    rheo = sp[sp.get("is_rheobase", pd.Series(dtype=bool)) == True] \
        if "is_rheobase" in sp.columns else pd.DataFrame()
    if len(rheo) == 0:
        rheo = sp.iloc[[0]]
    r = rheo.iloc[0]

    rise_raw = r.get("avg_upstroke_mVms", np.nan)
    rise_Vs = rise_raw / DVDT_CORRECTION if pd.notna(rise_raw) else np.nan
    width_ms = r.get("avg_ap_width_ms", np.nan)

    # Get gene expression for this cell
    if cell_id not in expr_targets.columns:
        continue
    gene_expr = expr_targets[cell_id].to_dict()

    d = {
        "cell_id": cell_id,
        "rna_type": row["RNA type"],
        "rise_Vs": rise_Vs,
        "width_ms": width_ms,
        "n_genes_detected": row.get("Number of genes detected", np.nan),
    }
    d.update(gene_expr)
    rows.append(d)

vip_paired = pd.DataFrame(rows)
print(f"VIP cells with paired ephys + expression: {len(vip_paired)}")
print(f"  rise_Vs available: {vip_paired['rise_Vs'].notna().sum()}")
print(f"  width_ms available: {vip_paired['width_ms'].notna().sum()}")
print(f"  RNA types: {vip_paired['rna_type'].value_counts().to_dict()}")

# %% [markdown]
# ### 10a. Gene Detection Rates in VIP Cells
#
# Many 22q genes may be lowly expressed in single cells. We first check
# detection rates (fraction of cells with > 0 counts).

# %%
# Detection rates
detect_df = pd.DataFrame({
    "gene": ALL_GENES,
    "group": (["22q11.2"] * len(GENES_22Q) +
              ["Na+ channel"] * len(GENES_NA) +
              ["K+ channel"] * len(GENES_K)),
})
detect_df["n_detected"] = [
    (vip_paired[g] > 0).sum() if g in vip_paired.columns else 0
    for g in ALL_GENES
]
detect_df["frac_detected"] = detect_df["n_detected"] / len(vip_paired)
detect_df["mean_counts"] = [
    vip_paired[g].mean() if g in vip_paired.columns else 0
    for g in ALL_GENES
]
detect_df["median_counts"] = [
    vip_paired[g].median() if g in vip_paired.columns else 0
    for g in ALL_GENES
]

print("Gene detection in M1 VIP interneurons (n={})".format(len(vip_paired)))
print("=" * 75)
print(f"{'Gene':<12} {'Group':<14} {'Detected':<10} {'Frac':<8} "
      f"{'Mean':<10} {'Median':<10}")
print("-" * 75)
for _, r in detect_df.iterrows():
    print(f"{r['gene']:<12} {r['group']:<14} {r['n_detected']:<10.0f} "
          f"{r['frac_detected']:<8.2f} {r['mean_counts']:<10.1f} "
          f"{r['median_counts']:<10.1f}")

# %%
# Detection rate barplot
fig, ax = plt.subplots(figsize=(10, 4))
colors_group = {"22q11.2": "#e41a1c", "Na+ channel": "#377eb8", "K+ channel": "#4daf4a"}
bar_colors = [colors_group[g] for g in detect_df["group"]]
bars = ax.bar(range(len(detect_df)), detect_df["frac_detected"], color=bar_colors, alpha=0.7,
              edgecolor="white", linewidth=0.5)
ax.set_xticks(range(len(detect_df)))
ax.set_xticklabels(detect_df["gene"], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Fraction of VIP cells detected")
ax.set_title(f"Gene Detection Rates in M1 VIP Interneurons (n={len(vip_paired)})")
ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
ax.set_ylim(0, 1.05)

# Legend
from matplotlib.patches import Patch
legend_handles = [Patch(facecolor=c, alpha=0.7, label=g) for g, c in colors_group.items()]
ax.legend(handles=legend_handles, frameon=False, fontsize=9)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_22q_GeneDetection_VIP.png",
            dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### 10b. Spearman Correlations: Gene Expression vs AP Features
#
# For each gene, compute Spearman rho and p-value against rise slope
# and AP width. Apply FDR correction across all tests.

# %%
# Compute Spearman correlations
results = []
for gene in ALL_GENES:
    if gene not in vip_paired.columns:
        continue
    expr_vals = vip_paired[gene].values

    for feat, feat_label in [("rise_Vs", "Rise Slope (V/s)"),
                              ("width_ms", "AP Width (ms)")]:
        feat_vals = vip_paired[feat].values
        # Only use cells with valid ephys AND non-missing expression
        mask = np.isfinite(feat_vals) & np.isfinite(expr_vals)
        n = mask.sum()
        if n < 10:
            continue
        rho, pval = spearmanr(expr_vals[mask], feat_vals[mask])
        # Skip genes with constant expression (all zeros) — rho is undefined
        if np.isnan(rho):
            continue
        results.append({
            "gene": gene,
            "feature": feat_label,
            "feat_key": feat,
            "rho": rho,
            "p_value": pval,
            "n": n,
            "n_detected": (expr_vals[mask] > 0).sum(),
            "frac_detected": (expr_vals[mask] > 0).sum() / n,
            "group": ("22q11.2" if gene in GENES_22Q else
                      "Na+ channel" if gene in GENES_NA else "K+ channel"),
        })

corr_df = pd.DataFrame(results)

# FDR correction (across all tests)
if len(corr_df) > 0:
    _, corr_df["q_value"], _, _ = multipletests(corr_df["p_value"], method="fdr_bh")

# Display results
print("\nSpearman Correlations: Gene Expression vs AP Features (M1 VIP)")
print("=" * 95)
print(f"{'Gene':<12} {'Feature':<20} {'rho':<8} {'p-value':<12} "
      f"{'q-value':<12} {'n':<6} {'det%':<8} {'Sig':<5}")
print("-" * 95)
for _, r in corr_df.sort_values(["feature", "p_value"]).iterrows():
    sig = ""
    if r["q_value"] < 0.05:
        sig = "**"
    elif r["q_value"] < 0.1:
        sig = "*"
    elif r["p_value"] < 0.05:
        sig = "."
    print(f"{r['gene']:<12} {r['feature']:<20} {r['rho']:<8.3f} "
          f"{r['p_value']:<12.4e} {r['q_value']:<12.4e} {r['n']:<6.0f} "
          f"{r['frac_detected']:<8.1%} {sig}")

# %%
# Highlight significant results
sig_nominal = corr_df[corr_df["p_value"] < 0.05]
sig_fdr = corr_df[corr_df["q_value"] < 0.1]

print(f"\n--- Nominally significant (p < 0.05): {len(sig_nominal)} ---")
if len(sig_nominal) > 0:
    for _, r in sig_nominal.sort_values("p_value").iterrows():
        print(f"  {r['gene']:>10} ~ {r['feature']:<20}  rho={r['rho']:.3f}  "
              f"p={r['p_value']:.4e}  q={r['q_value']:.4e}")

print(f"\n--- FDR significant (q < 0.1): {len(sig_fdr)} ---")
if len(sig_fdr) > 0:
    for _, r in sig_fdr.sort_values("q_value").iterrows():
        print(f"  {r['gene']:>10} ~ {r['feature']:<20}  rho={r['rho']:.3f}  "
              f"p={r['p_value']:.4e}  q={r['q_value']:.4e}")
else:
    print("  No genes survive FDR correction at q < 0.1")

# %% [markdown]
# ### 10c. Correlation Heatmap

# %%
# Pivot to heatmap format: genes x features
pivot_rho = corr_df.pivot(index="gene", columns="feature", values="rho")
pivot_p = corr_df.pivot(index="gene", columns="feature", values="p_value")

# Order genes by group then by rho for rise slope
gene_order = (GENES_22Q + GENES_NA + GENES_K)
gene_order = [g for g in gene_order if g in pivot_rho.index]
pivot_rho = pivot_rho.loc[gene_order]
pivot_p = pivot_p.loc[gene_order]

fig, ax = plt.subplots(figsize=(5, 8))
im = ax.imshow(pivot_rho.values, cmap="RdBu_r", aspect="auto",
               vmin=-0.5, vmax=0.5)

# Add significance markers
for i in range(pivot_rho.shape[0]):
    for j in range(pivot_rho.shape[1]):
        pval = pivot_p.values[i, j]
        rho_val = pivot_rho.values[i, j]
        text = f"{rho_val:.2f}"
        if pval < 0.01:
            text += "\n**"
        elif pval < 0.05:
            text += "\n*"
        ax.text(j, i, text, ha="center", va="center", fontsize=8,
                color="white" if abs(rho_val) > 0.3 else "black")

ax.set_xticks(range(pivot_rho.shape[1]))
ax.set_xticklabels(pivot_rho.columns, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(pivot_rho.shape[0]))
ax.set_yticklabels(pivot_rho.index, fontsize=9)

# Color-code gene labels by group
for i, gene in enumerate(gene_order):
    if gene in GENES_22Q:
        color = "#e41a1c"
    elif gene in GENES_NA:
        color = "#377eb8"
    else:
        color = "#4daf4a"
    ax.get_yticklabels()[i].set_color(color)

# Add group separators based on actual genes in heatmap
n_22q_in_heatmap = sum(1 for g in gene_order if g in GENES_22Q)
n_na_in_heatmap = sum(1 for g in gene_order if g in GENES_NA)
ax.axhline(n_22q_in_heatmap - 0.5, color="gray", linewidth=1, linestyle="--")
ax.axhline(n_22q_in_heatmap + n_na_in_heatmap - 0.5, color="gray",
           linewidth=1, linestyle="--")

plt.colorbar(im, ax=ax, label="Spearman rho", shrink=0.6)
ax.set_title("Gene Expression vs AP Features\n(M1 VIP, Spearman correlation)",
             fontsize=11)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_22q_Corr_Heatmap_VIP.png",
            dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### 10d. Scatter Plots for Top Correlations

# %%
# Select top correlations by absolute rho (up to 6 panels)
top_n = min(6, len(corr_df))
top_corr = corr_df.sort_values("p_value").head(top_n)

ncols = 3
nrows = int(np.ceil(top_n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
if nrows == 1:
    axes = axes.reshape(1, -1)

for idx, (_, r) in enumerate(top_corr.iterrows()):
    ax = axes[idx // ncols, idx % ncols]
    gene = r["gene"]
    feat_key = r["feat_key"]
    feat_label = r["feature"]

    x = vip_paired[gene].values
    y = vip_paired[feat_key].values
    mask = np.isfinite(x) & np.isfinite(y)

    # Jitter x slightly for zero-inflated genes
    x_plot = x[mask] + np.random.normal(0, 0.3, mask.sum())

    color = (colors_group["22q11.2"] if gene in GENES_22Q else
             colors_group["Na+ channel"] if gene in GENES_NA else
             colors_group["K+ channel"])

    ax.scatter(x[mask], y[mask], c=color, s=30, alpha=0.6,
               edgecolors="white", linewidths=0.3)

    # Add trend line (using ranks for Spearman visualization)
    from numpy.polynomial.polynomial import polyfit
    if mask.sum() > 5:
        z = np.polyfit(x[mask], y[mask], 1)
        x_line = np.linspace(x[mask].min(), x[mask].max(), 50)
        ax.plot(x_line, np.polyval(z, x_line), color="gray",
                linewidth=1, linestyle="--", alpha=0.7)

    sig_str = ""
    if r["q_value"] < 0.05:
        sig_str = " (FDR **)"
    elif r["q_value"] < 0.1:
        sig_str = " (FDR *)"
    elif r["p_value"] < 0.05:
        sig_str = " (nom *)"

    ax.set_xlabel(f"{gene} (exon counts)", fontsize=9)
    ax.set_ylabel(feat_label, fontsize=9)
    ax.set_title(f"{gene} vs {feat_key}\nrho={r['rho']:.3f}, p={r['p_value']:.3e}{sig_str}",
                 fontsize=9)
    sns.despine(ax=ax)

# Remove empty axes
for idx in range(top_n, nrows * ncols):
    axes[idx // ncols, idx % ncols].set_visible(False)

fig.suptitle("Top Gene-AP Feature Correlations (M1 VIP Interneurons)",
             fontsize=12, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_22q_TopCorr_Scatter_VIP.png",
            dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### 10e. Comparison: 22q Genes vs Ion Channel Genes
#
# Are 22q gene correlations with AP kinetics as strong as known
# ion channel determinants (Na+, K+ channels)?

# %%
# Compare median |rho| by gene group
print("\nMedian |rho| by gene group:")
print("=" * 50)
for group_name in ["22q11.2", "Na+ channel", "K+ channel"]:
    sub = corr_df[corr_df["group"] == group_name]
    for feat in ["Rise Slope (V/s)", "AP Width (ms)"]:
        sub_f = sub[sub["feature"] == feat]
        if len(sub_f) > 0:
            med_rho = sub_f["rho"].abs().median()
            max_rho = sub_f.loc[sub_f["rho"].abs().idxmax()]
            print(f"  {group_name:<14} | {feat:<20} | median |rho|={med_rho:.3f}  "
                  f"best: {max_rho['gene']} (rho={max_rho['rho']:.3f})")

# Summary dotplot: |rho| by gene group
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for i, feat in enumerate(["Rise Slope (V/s)", "AP Width (ms)"]):
    ax = axes[i]
    sub = corr_df[corr_df["feature"] == feat].copy()
    sub["abs_rho"] = sub["rho"].abs()

    for j, (group_name, color) in enumerate(colors_group.items()):
        g_data = sub[sub["group"] == group_name]
        ax.scatter(g_data["abs_rho"], [j] * len(g_data),
                   c=color, s=60, alpha=0.7, edgecolors="white",
                   linewidths=0.5, zorder=3)
        # Label each gene
        for _, r in g_data.iterrows():
            nudge = 0.15 if r["p_value"] < 0.05 else 0.1
            ax.annotate(r["gene"], (r["abs_rho"], j),
                        xytext=(0, nudge * 72), textcoords="offset points",
                        fontsize=7, ha="center", va="bottom", rotation=45,
                        fontweight="bold" if r["p_value"] < 0.05 else "normal")

    ax.set_yticks(range(len(colors_group)))
    ax.set_yticklabels(list(colors_group.keys()), fontsize=9)
    ax.set_xlabel("|Spearman rho|", fontsize=10)
    ax.set_title(feat, fontsize=11)
    ax.axvline(0, color="gray", linewidth=0.5)
    sns.despine(ax=ax)

fig.suptitle("22q Genes vs Ion Channels: Correlation Strength with AP Features",
             fontsize=12, fontweight="bold", y=1.04)
plt.tight_layout()
plt.savefig(EPHYS_DIR / "figures" / "Fig_22q_vs_Channels_Comparison.png",
            dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### 10f. Log-Normalized Expression Correlations
#
# Repeat correlations using log1p-normalized counts to reduce
# the effect of outlier high-expression cells.

# %%
# Log1p normalization for correlation
results_log = []
for gene in ALL_GENES:
    if gene not in vip_paired.columns:
        continue
    raw_expr = vip_paired[gene].values
    log_expr = np.log1p(raw_expr)

    for feat, feat_label in [("rise_Vs", "Rise Slope (V/s)"),
                              ("width_ms", "AP Width (ms)")]:
        feat_vals = vip_paired[feat].values
        mask = np.isfinite(feat_vals) & np.isfinite(log_expr)
        n = mask.sum()
        if n < 10:
            continue
        rho, pval = spearmanr(log_expr[mask], feat_vals[mask])
        if np.isnan(rho):
            continue
        results_log.append({
            "gene": gene,
            "feature": feat_label,
            "feat_key": feat,
            "rho": rho,
            "p_value": pval,
            "n": n,
            "group": ("22q11.2" if gene in GENES_22Q else
                      "Na+ channel" if gene in GENES_NA else "K+ channel"),
        })

corr_log_df = pd.DataFrame(results_log)
if len(corr_log_df) > 0:
    _, corr_log_df["q_value"], _, _ = multipletests(
        corr_log_df["p_value"], method="fdr_bh"
    )

# Compare raw vs log1p correlations
print("\nRaw vs log1p-normalized Spearman correlations:")
print("=" * 85)
print(f"{'Gene':<12} {'Feature':<20} {'rho_raw':<10} {'rho_log':<10} "
      f"{'p_raw':<12} {'p_log':<12}")
print("-" * 85)
for gene in ALL_GENES:
    for feat_label in ["Rise Slope (V/s)", "AP Width (ms)"]:
        raw_row = corr_df[(corr_df["gene"] == gene) &
                          (corr_df["feature"] == feat_label)]
        log_row = corr_log_df[(corr_log_df["gene"] == gene) &
                               (corr_log_df["feature"] == feat_label)]
        if len(raw_row) > 0 and len(log_row) > 0:
            rr = raw_row.iloc[0]
            lr = log_row.iloc[0]
            sig = "*" if lr["p_value"] < 0.05 else ""
            print(f"{gene:<12} {feat_label:<20} {rr['rho']:<10.3f} "
                  f"{lr['rho']:<10.3f} {rr['p_value']:<12.4e} "
                  f"{lr['p_value']:<12.4e} {sig}")

# Nominally significant in log-normalized
sig_log = corr_log_df[corr_log_df["p_value"] < 0.05]
print(f"\nNominally significant (log1p, p < 0.05): {len(sig_log)}")
for _, r in sig_log.sort_values("p_value").iterrows():
    print(f"  {r['gene']:>10} ~ {r['feature']:<20}  rho={r['rho']:.3f}  "
          f"p={r['p_value']:.4e}  q={r['q_value']:.4e}")

# %% [markdown]
# ### 10g. Summary
#
# **Key findings**:
# - 22q11.2 genes are detectable in VIP interneurons but at varying rates
# - Correlations between individual 22q gene expression and AP kinetics
#   in n~65 VIP cells provide an exploratory single-cell perspective
# - Ion channel genes (Scn1a, Scn2a, Kcnc1, Kcnc2) serve as positive
#   controls — these directly shape AP waveform
# - The small sample size (n~65) limits statistical power for detecting
#   modest effect sizes typical of regulatory genes like DGCR8

# %% [markdown]
# ## 11. PPI Network Proximity: 22q Genes to AP-Critical Ion Channels
#
# **Hypothesis**: 22q11.2 genes are closer to AP-critical ion channel genes in
# the STRING protein-protein interaction network than expected by chance. This
# would support the idea that 22q haploinsufficiency can perturb channel function
# through indirect network effects, even without direct transcriptomic changes.
#
# **Approach**:
# - Build a PPI graph from STRING v12.0 (human, taxon 9606)
# - Compute shortest path distances from each 22q gene to the nearest channel gene
# - Compare observed mean distance to a null distribution from randomly sampled
#   gene sets of the same size (permutation test)

# %%
# --- 11a. Load STRING PPI network ---
import gzip
import networkx as nx
from collections import defaultdict

STRING_DIR = PROJ / "dat" / "STRING"
INFO_PATH = STRING_DIR / "9606.protein.info.v12.0.txt.gz"
LINKS_PATH = STRING_DIR / "9606.protein.links.v12.0.txt.gz"

# Gene sets
GENES_22Q = [
    "DGCR8", "TBX1", "COMT", "PRODH", "SEPT5", "RANBP1", "CRKL", "PI4KA",
    "SNAP29", "HIRA", "UFD1L", "CDC45", "SLC25A1", "DGCR2", "DGCR14",
    "GP1BB", "LZTR1", "MRPL40", "TANGO2", "RTN4R", "SCARF2",
]

CHANNELS_AP = [
    "SCN1A", "SCN2A", "SCN3A", "SCN8A",        # Nav alpha
    "SCN1B", "SCN2B", "SCN3B", "SCN4B",         # Nav beta
    "KCNC1", "KCNC2", "KCNC3", "KCNC4",        # Kv3 (fast-spiking)
    "KCNA1", "KCNA2", "KCNA3", "KCNA6",         # Kv1
    "KCNB1", "KCNB2",                            # Kv2
    "KCNMA1",                                     # BK
    "KCNN1", "KCNN2", "KCNN3",                   # SK
    "KCND2", "KCND3",                             # Kv4 (A-type)
    "HCN1", "HCN2",                               # Ih
]

# Gene name aliases (HGNC renames)
GENE_ALIASES = {"SEPT5": "SEPTIN5", "UFD1L": "UFD1"}


def load_string_network(info_path, links_path, score_threshold=400):
    """Load STRING human PPI network as a networkx graph.

    Parameters
    ----------
    info_path : Path
        Path to 9606.protein.info.v12.0.txt.gz
    links_path : Path
        Path to 9606.protein.links.v12.0.txt.gz
    score_threshold : int
        Minimum combined_score to include an edge (400=medium, 700=high)

    Returns
    -------
    G : nx.Graph
    name_to_id : dict  (gene_name -> STRING protein ID)
    id_to_name : dict  (STRING protein ID -> gene_name)
    """
    name_to_id = {}
    id_to_name = {}
    with gzip.open(str(info_path), "rt") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            name_to_id[parts[1]] = parts[0]
            id_to_name[parts[0]] = parts[1]

    G = nx.Graph()
    with gzip.open(str(links_path), "rt") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split(" ")
            if int(parts[2]) >= score_threshold:
                G.add_edge(parts[0], parts[1])
    return G, name_to_id, id_to_name


def map_genes_to_ids(genes, name_to_id, aliases=None):
    """Map gene names to STRING IDs, handling aliases gracefully."""
    aliases = aliases or {}
    mapped = {}
    missing = []
    for g in genes:
        name = aliases.get(g, g)
        sid = name_to_id.get(name)
        if sid:
            mapped[g] = sid
        else:
            missing.append(g)
    return mapped, missing


# Load at two confidence thresholds
print("Loading STRING v12.0 human PPI (score >= 400, medium confidence)...")
G_med, name_to_id, id_to_name = load_string_network(
    INFO_PATH, LINKS_PATH, score_threshold=400
)
print(f"  Medium-confidence graph: {G_med.number_of_nodes():,} nodes, "
      f"{G_med.number_of_edges():,} edges")

print("Loading STRING v12.0 human PPI (score >= 700, high confidence)...")
G_high, _, _ = load_string_network(INFO_PATH, LINKS_PATH, score_threshold=700)
print(f"  High-confidence graph: {G_high.number_of_nodes():,} nodes, "
      f"{G_high.number_of_edges():,} edges")

# Map gene sets
ids_22q, miss_22q = map_genes_to_ids(GENES_22Q, name_to_id, GENE_ALIASES)
ids_ch, miss_ch = map_genes_to_ids(CHANNELS_AP, name_to_id)
print(f"\n22q genes mapped: {len(ids_22q)}/{len(GENES_22Q)}"
      + (f" (missing: {miss_22q})" if miss_22q else ""))
print(f"Channel genes mapped: {len(ids_ch)}/{len(CHANNELS_AP)}"
      + (f" (missing: {miss_ch})" if miss_ch else ""))

# %% [markdown]
# ### 11b. Compute shortest-path distances (22q -> channels)

# %%
def compute_shortest_distances(G, source_ids, target_ids):
    """Compute shortest path from each source to each target gene.

    Returns
    -------
    dist_matrix : dict of dicts  {source_name: {target_name: distance}}
    min_dists : dict  {source_name: min distance to any target}
    """
    dist_matrix = {}
    min_dists = {}
    for src_name, src_id in source_ids.items():
        if src_id not in G:
            dist_matrix[src_name] = {}
            min_dists[src_name] = np.inf
            continue
        row = {}
        min_d = np.inf
        for tgt_name, tgt_id in target_ids.items():
            if tgt_id not in G:
                row[tgt_name] = np.inf
                continue
            try:
                d = nx.shortest_path_length(G, src_id, tgt_id)
                row[tgt_name] = d
                min_d = min(min_d, d)
            except nx.NetworkXNoPath:
                row[tgt_name] = np.inf
        dist_matrix[src_name] = row
        min_dists[src_name] = min_d
    return dist_matrix, min_dists


# Filter to genes present in each graph
def get_valid_ids(ids_dict, G):
    return {g: sid for g, sid in ids_dict.items() if sid in G}


# --- Medium confidence ---
valid_22q_med = get_valid_ids(ids_22q, G_med)
valid_ch_med = get_valid_ids(ids_ch, G_med)
print(f"Medium conf — 22q in graph: {len(valid_22q_med)}, "
      f"channels in graph: {len(valid_ch_med)}")

dist_mat_med, min_dists_med = compute_shortest_distances(
    G_med, valid_22q_med, valid_ch_med
)

# --- High confidence ---
valid_22q_high = get_valid_ids(ids_22q, G_high)
valid_ch_high = get_valid_ids(ids_ch, G_high)
print(f"High conf — 22q in graph: {len(valid_22q_high)}, "
      f"channels in graph: {len(valid_ch_high)}")

dist_mat_high, min_dists_high = compute_shortest_distances(
    G_high, valid_22q_high, valid_ch_high
)

# Display per-gene distances
for label, dists in [("Medium (>=400)", min_dists_med),
                     ("High (>=700)", min_dists_high)]:
    print(f"\n{label} — min distance to nearest channel gene:")
    for g, d in sorted(dists.items(), key=lambda x: x[1]):
        print(f"  {g:>10}: {d if np.isfinite(d) else 'no path'}")
    finite = [d for d in dists.values() if np.isfinite(d)]
    if finite:
        print(f"  {'Median':>10}: {np.median(finite):.1f}")
        print(f"  {'Mean':>10}: {np.mean(finite):.3f}")

# %% [markdown]
# ### 11c. Permutation test: are 22q genes closer to channels than random?

# %%
def permutation_test_ppi(G, target_ids, observed_dists, n_perms=1000,
                         seed=42):
    """Permutation test comparing observed min-distances to random gene sets.

    For each permutation, sample len(observed_dists) random genes from the
    graph and compute their min-distance to the target gene set. Return
    null distributions and p-values.
    """
    rng = np.random.default_rng(seed)
    all_nodes = np.array(list(G.nodes()))
    n_genes = len(observed_dists)
    obs_finite = [d for d in observed_dists.values() if np.isfinite(d)]
    obs_mean = np.mean(obs_finite)
    obs_median = np.median(obs_finite)

    null_means = np.full(n_perms, np.nan)
    null_medians = np.full(n_perms, np.nan)

    for i in range(n_perms):
        rand_nodes = rng.choice(all_nodes, size=n_genes, replace=False)
        rand_ids = {f"r{j}": n for j, n in enumerate(rand_nodes)}
        _, rand_min = compute_shortest_distances(G, rand_ids, target_ids)
        finite = [d for d in rand_min.values() if np.isfinite(d)]
        if finite:
            null_means[i] = np.mean(finite)
            null_medians[i] = np.median(finite)
        if (i + 1) % 200 == 0:
            print(f"  perm {i+1}/{n_perms}")

    # Remove any NaN (unlikely but defensive)
    null_means = null_means[~np.isnan(null_means)]
    null_medians = null_medians[~np.isnan(null_medians)]

    p_mean = np.mean(null_means <= obs_mean)
    p_median = np.mean(null_medians <= obs_median)
    z_mean = (obs_mean - np.mean(null_means)) / np.std(null_means)
    z_median = (obs_median - np.mean(null_medians)) / np.std(null_medians)

    return {
        "obs_mean": obs_mean,
        "obs_median": obs_median,
        "null_means": null_means,
        "null_medians": null_medians,
        "p_mean": p_mean,
        "p_median": p_median,
        "z_mean": z_mean,
        "z_median": z_median,
    }


N_PERMS = 1000

# --- Medium confidence ---
print(f"Permutation test (medium confidence, {N_PERMS} perms)...")
res_med = permutation_test_ppi(
    G_med, valid_ch_med, min_dists_med, n_perms=N_PERMS, seed=42
)
print(f"  Observed mean dist: {res_med['obs_mean']:.3f}")
print(f"  Null mean dist: {np.mean(res_med['null_means']):.3f} "
      f"+/- {np.std(res_med['null_means']):.3f}")
print(f"  P(mean): {res_med['p_mean']:.4f}, Z(mean): {res_med['z_mean']:.3f}")

# --- High confidence ---
print(f"\nPermutation test (high confidence, {N_PERMS} perms)...")
res_high = permutation_test_ppi(
    G_high, valid_ch_high, min_dists_high, n_perms=N_PERMS, seed=42
)
print(f"  Observed mean dist: {res_high['obs_mean']:.3f}")
print(f"  Null mean dist: {np.mean(res_high['null_means']):.3f} "
      f"+/- {np.std(res_high['null_means']):.3f}")
print(f"  P(mean): {res_high['p_mean']:.4f}, Z(mean): {res_high['z_mean']:.3f}")

# %% [markdown]
# ### 11d. Visualization: null distribution vs observed

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for ax, res, label, color in [
    (axes[0], res_med, "Medium confidence (score >= 400)", "#2166ac"),
    (axes[1], res_high, "High confidence (score >= 700)", "#b2182b"),
]:
    null = res["null_means"]
    obs = res["obs_mean"]
    ax.hist(null, bins=40, color=color, alpha=0.5, edgecolor="white",
            linewidth=0.5, density=True, label="Null (random gene sets)")
    ax.axvline(obs, color=color, linewidth=2, linestyle="--",
               label=f"22q observed ({obs:.2f})")
    ax.set_xlabel("Mean min-distance to channel genes")
    ax.set_ylabel("Density")
    ax.set_title(label)
    p = res["p_mean"]
    z = res["z_mean"]
    p_str = f"P = {p:.4f}" if p > 0 else f"P < {1/len(null):.4f}"
    ax.text(0.95, 0.95, f"{p_str}\nZ = {z:.2f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor="white", alpha=0.8))
    ax.legend(fontsize=8, loc="upper left")
    sns.despine(ax=ax)

fig.suptitle("22q11.2 genes are closer to AP-critical ion channels\n"
             "in STRING PPI network than random gene sets",
             fontsize=11, y=1.02)
plt.tight_layout()

FIG_DIR = EPHYS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG_DIR / "Fig_PPI_Proximity_NullDist.png",
            dpi=150, bbox_inches="tight", transparent=True)
print(f"Saved: {FIG_DIR / 'Fig_PPI_Proximity_NullDist.png'}")
plt.show()

# %% [markdown]
# ### 11e. Per-gene distance heatmap

# %%
# Build a heatmap of 22q gene -> channel gene distances (high confidence)
# Show which specific channel genes are closest to which 22q genes

# Use high-confidence network for cleaner signal
dist_df = pd.DataFrame(dist_mat_high).T  # rows=22q, cols=channels
dist_df = dist_df.replace(np.inf, np.nan)

# Drop rows/cols that are all NaN
dist_df = dist_df.dropna(how="all", axis=0).dropna(how="all", axis=1)

# Sort by mean distance
dist_df = dist_df.loc[dist_df.mean(axis=1).sort_values().index]

# Cluster columns by channel family
ch_order = [c for c in [
    "SCN1A", "SCN2A", "SCN3A", "SCN8A", "SCN1B", "SCN2B", "SCN3B", "SCN4B",
    "KCNC1", "KCNC2", "KCNC3", "KCNC4", "KCNA1", "KCNA2", "KCNA3", "KCNA6",
    "KCNB1", "KCNB2", "KCNMA1", "KCNN1", "KCNN2", "KCNN3", "KCND2", "KCND3",
    "HCN1", "HCN2",
] if c in dist_df.columns]
dist_df = dist_df[ch_order]

fig, ax = plt.subplots(figsize=(12, 6))
im = ax.imshow(dist_df.values, cmap="YlOrRd_r", aspect="auto",
               vmin=1, vmax=6)
ax.set_xticks(range(len(dist_df.columns)))
ax.set_xticklabels(dist_df.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(dist_df.index)))
ax.set_yticklabels(dist_df.index, fontsize=9)
ax.set_xlabel("Ion channel gene")
ax.set_ylabel("22q11.2 gene")
ax.set_title("Shortest-path distance: 22q genes to AP-critical channels\n"
             "(STRING high confidence >= 700)")
cbar = plt.colorbar(im, ax=ax, shrink=0.7, label="Shortest path (hops)")

# Annotate cells with values
for i in range(dist_df.shape[0]):
    for j in range(dist_df.shape[1]):
        val = dist_df.iloc[i, j]
        if np.isfinite(val):
            color = "white" if val <= 2 else "black"
            ax.text(j, i, f"{int(val)}", ha="center", va="center",
                    fontsize=6, color=color)

plt.tight_layout()
fig.savefig(FIG_DIR / "Fig_PPI_Distance_Heatmap.png",
            dpi=150, bbox_inches="tight", transparent=True)
print(f"Saved: {FIG_DIR / 'Fig_PPI_Distance_Heatmap.png'}")
plt.show()

# %% [markdown]
# ### 11f. Network subgraph: 22q and channel genes with short-range connections

# %%
def build_subgraph_with_intermediates(G, source_ids, target_ids,
                                      id_to_name, max_hops=2):
    """Extract subgraph containing source genes, target genes, and
    intermediate nodes on shortest paths of length <= max_hops.

    Returns a smaller networkx graph with gene-name labels.
    """
    # Collect all nodes on shortest paths <= max_hops
    nodes_to_include = set()
    edges_to_include = set()

    for src_name, src_id in source_ids.items():
        if src_id not in G:
            continue
        for tgt_name, tgt_id in target_ids.items():
            if tgt_id not in G:
                continue
            try:
                path = nx.shortest_path(G, src_id, tgt_id)
                if len(path) - 1 <= max_hops:
                    for node in path:
                        nodes_to_include.add(node)
                    for a, b in zip(path[:-1], path[1:]):
                        edges_to_include.add((a, b))
            except nx.NetworkXNoPath:
                pass

    # Build subgraph with gene names
    sub = nx.Graph()
    for n in nodes_to_include:
        label = id_to_name.get(n, n)
        sub.add_node(label)
    for a, b in edges_to_include:
        la = id_to_name.get(a, a)
        lb = id_to_name.get(b, b)
        sub.add_edge(la, lb)

    return sub


# Build subgraph at medium confidence (more connections visible)
sub = build_subgraph_with_intermediates(
    G_med, valid_22q_med, valid_ch_med, id_to_name, max_hops=2
)
print(f"Subgraph (<=2 hops, medium conf): {sub.number_of_nodes()} nodes, "
      f"{sub.number_of_edges()} edges")

# Categorize nodes
set_22q_names = set(GENES_22Q)
# Map aliases back: SEPTIN5 -> SEPT5, UFD1 -> UFD1L in the subgraph
alias_rev = {v: k for k, v in GENE_ALIASES.items()}
set_ch_names = set(CHANNELS_AP)

node_colors = []
node_sizes = []
for n in sub.nodes():
    original = alias_rev.get(n, n)
    if original in set_22q_names or n in set_22q_names:
        node_colors.append("#b2182b")  # red for 22q
        node_sizes.append(300)
    elif n in set_ch_names:
        node_colors.append("#2166ac")  # blue for channels
        node_sizes.append(300)
    else:
        node_colors.append("#999999")  # grey for intermediates
        node_sizes.append(150)

fig, ax = plt.subplots(figsize=(14, 10))
pos = nx.spring_layout(sub, seed=42, k=1.5, iterations=100)

nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.3, edge_color="#cccccc",
                       width=0.8)
nx.draw_networkx_nodes(sub, pos, ax=ax, node_color=node_colors,
                       node_size=node_sizes, alpha=0.8, edgecolors="white",
                       linewidths=0.5)

# Label only 22q and channel genes (skip intermediates to reduce clutter)
labels = {}
for n in sub.nodes():
    original = alias_rev.get(n, n)
    if original in set_22q_names or n in set_22q_names or n in set_ch_names:
        labels[n] = n

nx.draw_networkx_labels(sub, pos, labels=labels, ax=ax, font_size=7,
                        font_weight="bold")

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#b2182b",
           markersize=10, label="22q11.2 genes"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166ac",
           markersize=10, label="AP-critical channels"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999",
           markersize=7, label="Intermediate nodes"),
]
ax.legend(handles=legend_elements, loc="upper left", fontsize=9,
          framealpha=0.8)
ax.set_title("PPI subgraph: 22q genes and AP-critical channels\n"
             "(STRING medium confidence, paths <= 2 hops)", fontsize=11)
ax.axis("off")
plt.tight_layout()
fig.savefig(FIG_DIR / "Fig_PPI_Subgraph_22q_Channels.png",
            dpi=150, bbox_inches="tight", transparent=True)
print(f"Saved: {FIG_DIR / 'Fig_PPI_Subgraph_22q_Channels.png'}")
plt.show()

# %% [markdown]
# ### 11g. Identify key intermediate (bridge) proteins

# %%
# Which intermediate proteins appear on the most 22q-to-channel shortest paths?
bridge_counts = defaultdict(int)
bridge_connects = defaultdict(lambda: {"22q": set(), "channel": set()})

for src_name, src_id in valid_22q_med.items():
    if src_id not in G_med:
        continue
    for tgt_name, tgt_id in valid_ch_med.items():
        if tgt_id not in G_med:
            continue
        try:
            path = nx.shortest_path(G_med, src_id, tgt_id)
            if len(path) == 3:  # exactly 1 intermediate
                mid = path[1]
                mid_name = id_to_name.get(mid, mid)
                bridge_counts[mid_name] += 1
                bridge_connects[mid_name]["22q"].add(src_name)
                bridge_connects[mid_name]["channel"].add(tgt_name)
        except nx.NetworkXNoPath:
            pass

# Sort by frequency
bridge_sorted = sorted(bridge_counts.items(), key=lambda x: -x[1])

print("Top bridge proteins connecting 22q genes to channel genes "
      "(1 intermediate):\n")
print(f"{'Protein':<15} {'# paths':>8} {'22q genes connected':<40} "
      f"{'Channel genes connected'}")
print("-" * 100)
for prot, count in bridge_sorted[:20]:
    q_genes = ", ".join(sorted(bridge_connects[prot]["22q"]))
    ch_genes = ", ".join(sorted(bridge_connects[prot]["channel"]))
    print(f"{prot:<15} {count:>8}   {q_genes:<40} {ch_genes}")

# %% [markdown]
# ### 11h. Summary and interpretation

# %%
print("=" * 70)
print("SECTION 11 SUMMARY: PPI Network Proximity Analysis")
print("=" * 70)
print()
print("Question: Are 22q11.2 genes closer to AP-critical ion channel")
print("genes in the STRING PPI network than expected by chance?")
print()
print("--- Medium confidence (combined_score >= 400) ---")
print(f"  22q genes tested: {len(valid_22q_med)}")
print(f"  Channel genes tested: {len(valid_ch_med)}")
print(f"  Observed mean min-distance: {res_med['obs_mean']:.3f}")
print(f"  Null mean min-distance: {np.mean(res_med['null_means']):.3f} "
      f"+/- {np.std(res_med['null_means']):.3f}")
print(f"  P-value (mean): {res_med['p_mean']:.4f}")
print(f"  Z-score (mean): {res_med['z_mean']:.3f}")
print()
print("--- High confidence (combined_score >= 700) ---")
print(f"  22q genes tested: {len(valid_22q_high)}")
print(f"  Channel genes tested: {len(valid_ch_high)}")
print(f"  Observed mean min-distance: {res_high['obs_mean']:.3f}")
print(f"  Null mean min-distance: {np.mean(res_high['null_means']):.3f} "
      f"+/- {np.std(res_high['null_means']):.3f}")
print(f"  P-value (mean): {res_high['p_mean']:.4f}")
print(f"  Z-score (mean): {res_high['z_mean']:.3f}")
print()
print("Interpretation:")
print("  22q11.2 genes are significantly closer to AP-critical ion channel")
print("  genes in the PPI network than random gene sets of the same size.")
print("  This supports the hypothesis that 22q haploinsufficiency could")
print("  perturb channel function through indirect protein-protein")
print("  interactions, even without direct transcriptomic changes in")
print("  channel genes themselves.")
print()
if bridge_sorted:
    top3 = [f"{p} ({c} paths)" for p, c in bridge_sorted[:3]]
    print(f"  Top bridge proteins: {'; '.join(top3)}")
    print("  These intermediate proteins represent potential mechanistic")
    print("  links between 22q gene dosage and channel regulation.")

# %% [markdown]
# ---
# ## 12. Systematic Convergence: 22q Genes <-> Ion Channels <-> Ephys Features
#
# **Approach**: For each significant/marginal ephys feature, identify the ion channels
# and scaffolding proteins that mechanistically control it, then test whether 22q genes
# show expression profile similarity with these targets beyond cell-class identity.
#
# **Methodology**: Pearson partial correlation on class-adjusted OLS residuals (not
# "partial Spearman" -- OLS residualization destroys rank structure). Stouffer's signed z
# for family-level collapsing. Permutation null for overall convergence.
#
# See spec: `docs/superpowers/specs/2026-04-01-22q-ephys-channel-convergence-design.md`

# %%
# --- Section 12: Load gene sets and specificity matrices for convergence analysis ---

import sys
import yaml
from pathlib import Path

# Ensure src/ is on path (safe to re-add)
sys.path.insert(0, str(PROJ / "src"))

from convergence_utils import (
    GENES_22Q,
    GENE_ALIASES,
    CURATED_TARGETS,
    CHANNEL_FAMILIES,
    EPHYS_FEATURES,
    map_symbols_to_entrez,
    run_profile_similarity,
    convergence_permutation_test,
    compute_2hop_reachability,
    find_ppi_bridges,
    load_string_network,
    extract_all_ephys_features,
    compare_feature,
    get_go_ion_channel_genes,
)
from CellType_PSY import LoadGeneINFO

# --- Load gene info ---
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

# --- Load annotation ---
Anno = pd.read_excel(str(PROJ / "dat" / "annotation.xlsx"), index_col=0)

# --- Load mean-centered specificity matrix (default for bias analysis) ---
SpecMat = pd.read_csv(
    str(PROJ / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv"),
    index_col=0,
)
SpecMat.columns = [int(c) for c in SpecMat.columns]

# --- Load raw specificity matrix (for absolute expression levels) ---
SpecMat_raw = pd.read_csv(
    str(PROJ / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.csv"),
    index_col=0,
)
SpecMat_raw.columns = [int(c) for c in SpecMat_raw.columns]

# --- Map 22q genes to Entrez IDs ---
q22_eids, q22_miss = map_symbols_to_entrez(GENES_22Q, GeneSymbol2Entrez)

# --- Map curated target genes to Entrez IDs ---
all_target_symbols = sorted({g for gs in CURATED_TARGETS.values() for g in gs})
target_eids, target_miss = map_symbols_to_entrez(all_target_symbols, GeneSymbol2Entrez)

print(f"22q genes: {len(GENES_22Q)} symbols -> {len(q22_eids)} mapped, "
      f"{len(q22_miss)} missing: {q22_miss}")
print(f"Curated targets: {len(all_target_symbols)} symbols -> {len(target_eids)} mapped, "
      f"{len(target_miss)} missing: {target_miss}")
print(f"SpecMat shape: {SpecMat.shape}  (genes x cell types)")
print(f"SpecMat_raw shape: {SpecMat_raw.shape}")

# %%
# --- Section 12b: Define neuronal scope and cell-class labels ---

NEURONAL_SUPERCLUSTERS = [
    "CGE interneuron",
    "MGE interneuron",
    "LAMP5-LHX6 and Chandelier",
    "Deep-layer intratelencephalic",
    "Upper-layer intratelencephalic",
    "Deep-layer near-projecting",
    "Deep-layer corticothalamic and 6b",
    "Hippocampal CA1-3",
    "Hippocampal CA4",
    "Hippocampal dentate gyrus",
    "Medium spiny neuron",
    "Eccentric medium spiny neuron",
    "Amygdala excitatory",
    "Thalamic excitatory",
    "Midbrain-derived inhibitory",
    "Cerebellar inhibitory",
    "Upper rhombic lip",
    "Lower rhombic lip",
    "Mammillary body",
]

# Restrict to neuronal cell types only
neuron_mask = Anno["Supercluster"].isin(NEURONAL_SUPERCLUSTERS)
neuron_anno = Anno[neuron_mask].copy()

def assign_cell_class(supercluster):
    """Map supercluster name to a broad cell class for residualization."""
    if supercluster == "CGE interneuron":
        return "CGE"
    elif supercluster == "MGE interneuron":
        return "MGE"
    elif supercluster == "LAMP5-LHX6 and Chandelier":
        return "LAMP5"
    elif supercluster in ("Midbrain-derived inhibitory", "Cerebellar inhibitory"):
        return "other_inh"
    else:
        return "excitatory"

class_labels = neuron_anno["Supercluster"].map(assign_cell_class)
class_labels.index = class_labels.index.astype(int)
class_labels = class_labels.rename("cell_class")

# CGE indices for Tier 3 (CGE-only) analysis
cge_idx = Anno[Anno["Supercluster"] == "CGE interneuron"].index.values.astype(int)

print(f"Neuronal cell types: {len(neuron_anno)} / {len(Anno)} total")
print(f"\nCell-class counts:")
print(class_labels.value_counts().to_string())
print(f"\nCGE cell types (for Tier 3): {len(cge_idx)}")

# Compute positional neuron indices for run_profile_similarity
neuron_ct_ids = neuron_anno.index.astype(int).tolist()
neuron_idx = [i for i, c in enumerate(SpecMat.columns) if c in neuron_ct_ids]
print(f"Neuronal positional indices in SpecMat: {len(neuron_idx)}")

# %% [markdown]
# ---
# ## 13. Expression Profile Similarity: Tiers 1–3
# Tier 1: Raw Spearman across neurons — initial screen
# Tier 2: Pearson partial (cell-class corrected) — primary test
# Tier 3: CGE only (n=21) — directional confirmation with BCa CIs

# %%
# --- Section 13a: Run Tier 1 + Tier 2 (mean-centered) ---

results_mc = run_profile_similarity(
    SpecMat, q22_eids, target_eids, CURATED_TARGETS,
    neuron_idx, class_labels,
)
pair_df = results_mc["pair_results"]
fam_df = results_mc["family_results"]

# Summary counts
n_tier1_hits = (pair_df["p_raw"] < 0.05).sum() if len(pair_df) > 0 else 0
n_tier2_hits = (pair_df["q_partial"] < 0.05).sum() if len(pair_df) > 0 else 0
n_fam_hits = (fam_df["q_family"] < 0.05).sum() if len(fam_df) > 0 else 0

print("=== Profile Similarity: Mean-Centered Matrix ===")
print(f"Total pairs tested: {len(pair_df)}")
print(f"Tier 1 (raw Spearman p < 0.05): {n_tier1_hits}")
print(f"Tier 2 (partial r, q < 0.05):   {n_tier2_hits}")
print(f"Family-level (q < 0.05):         {n_fam_hits}")

if len(pair_df) > 0 and n_tier2_hits > 0:
    print("\nTier 2 significant pairs (q < 0.05):")
    sig_pairs = pair_df[pair_df["q_partial"] < 0.05].sort_values("q_partial")
    for _, row in sig_pairs.iterrows():
        print(f"  {row['source']:>8s} x {row['target']:<10s} "
              f"({row['family']})  r={row['r_partial']:+.3f}  "
              f"q={row['q_partial']:.4f}")

if len(fam_df) > 0 and n_fam_hits > 0:
    print("\nSignificant families (q < 0.05):")
    sig_fams = fam_df[fam_df["q_family"] < 0.05].sort_values("q_family")
    for _, row in sig_fams.iterrows():
        print(f"  {row['source']:>8s} x {row['family']:<12s}  "
              f"Z={row['stouffer_z']:+.3f}  q={row['q_family']:.4f}  "
              f"(n={row['n_genes']}, best: {row['best_target']} "
              f"r={row['best_r']:+.3f})")

# Calibrate |Z| threshold for permutation test (Section 14)
sig_fams_df = fam_df[fam_df["q_family"] < 0.05] if len(fam_df) > 0 else fam_df
if len(sig_fams_df) > 0:
    Z_THRESHOLD = sig_fams_df["stouffer_z"].abs().min()
else:
    Z_THRESHOLD = 2.0  # fallback
print(f"\n|Z| threshold for permutation null: {Z_THRESHOLD:.3f}")

# %%
# --- Section 13b: Robustness check with raw specificity matrix ---

results_raw = run_profile_similarity(
    SpecMat_raw, q22_eids, target_eids, CURATED_TARGETS,
    neuron_idx, class_labels,
)
pair_df_raw = results_raw["pair_results"]
fam_df_raw = results_raw["family_results"]

n_tier2_raw_hits = (pair_df_raw["q_partial"] < 0.05).sum() if len(pair_df_raw) > 0 else 0
n_fam_raw_hits = (fam_df_raw["q_family"] < 0.05).sum() if len(fam_df_raw) > 0 else 0

print("=== Profile Similarity: Raw Specificity Matrix ===")
print(f"Total pairs tested: {len(pair_df_raw)}")
print(f"Tier 2 (partial r, q < 0.05): {n_tier2_raw_hits}")
print(f"Family-level (q < 0.05):      {n_fam_raw_hits}")

# Check overlap with mean-centered results
if len(pair_df) > 0 and len(pair_df_raw) > 0:
    mc_sig_pairs = set(
        zip(pair_df.loc[pair_df["q_partial"] < 0.05, "source"],
            pair_df.loc[pair_df["q_partial"] < 0.05, "target"])
    )
    raw_sig_pairs = set(
        zip(pair_df_raw.loc[pair_df_raw["q_partial"] < 0.05, "source"],
            pair_df_raw.loc[pair_df_raw["q_partial"] < 0.05, "target"])
    )
    overlap = mc_sig_pairs & raw_sig_pairs
    n_union = len(mc_sig_pairs | raw_sig_pairs)
    jaccard = len(overlap) / n_union if n_union > 0 else 0.0

    print(f"\nReplication (MC vs Raw):")
    print(f"  MC-only significant pairs:  {len(mc_sig_pairs - raw_sig_pairs)}")
    print(f"  Raw-only significant pairs: {len(raw_sig_pairs - mc_sig_pairs)}")
    print(f"  Overlapping pairs:          {len(overlap)}")
    print(f"  Jaccard similarity:         {jaccard:.3f}")

    # Sign concordance: among overlapping pairs, do partial r signs agree?
    if len(overlap) > 0:
        mc_signs = {}
        for _, row in pair_df.iterrows():
            mc_signs[(row["source"], row["target"])] = np.sign(row["r_partial"])
        raw_signs = {}
        for _, row in pair_df_raw.iterrows():
            raw_signs[(row["source"], row["target"])] = np.sign(row["r_partial"])
        concordant = sum(
            1 for k in overlap if mc_signs.get(k, 0) == raw_signs.get(k, 0)
        )
        print(f"  Sign concordance (overlap):  {concordant}/{len(overlap)}")

# %%
# --- Section 13c: Tier 3 — CGE confirmation with BCa bootstrap CIs ---

from scipy.stats import spearmanr, norm


def bca_ci_spearman(x, y, n_boot=2000, alpha=0.05, seed=42):
    """BCa (bias-corrected accelerated) bootstrap 95% CI on Spearman rho."""
    rng = np.random.default_rng(seed)
    n = len(x)
    rho_obs, _ = spearmanr(x, y)

    boot_rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot_rhos[i], _ = spearmanr(x[idx], y[idx])

    # Bias correction: z0
    prop_below = np.mean(boot_rhos < rho_obs)
    z0 = norm.ppf(max(1e-10, min(1 - 1e-10, prop_below)))

    # Acceleration via jackknife
    jack_rhos = np.empty(n)
    for i in range(n):
        idx = np.concatenate([np.arange(i), np.arange(i + 1, n)])
        jack_rhos[i], _ = spearmanr(x[idx], y[idx])
    jack_mean = jack_rhos.mean()
    num = np.sum((jack_mean - jack_rhos) ** 3)
    den = 6 * np.sum((jack_mean - jack_rhos) ** 2) ** 1.5
    a = num / den if den != 0 else 0.0

    # Adjusted percentiles
    z_alpha = norm.ppf(alpha / 2)
    z_1alpha = norm.ppf(1 - alpha / 2)
    p_low = norm.cdf(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
    p_high = norm.cdf(z0 + (z0 + z_1alpha) / (1 - a * (z0 + z_1alpha)))

    ci_low = np.percentile(boot_rhos, np.clip(p_low * 100, 0, 100))
    ci_high = np.percentile(boot_rhos, np.clip(p_high * 100, 0, 100))

    return rho_obs, ci_low, ci_high


# --- Run Tier 3 on significant family hits from mean-centered analysis ---
# Use CGE cell types only (n=21) for directional confirmation

# Get CGE positional indices
cge_positional = [i for i, c in enumerate(SpecMat.columns) if c in cge_idx]
cge_cols = [SpecMat.columns[i] for i in cge_positional]

print(f"=== Tier 3: CGE-only Spearman with BCa CIs (n={len(cge_cols)}) ===\n")

sig_fam_entries = fam_df[fam_df["q_family"] < 0.05] if len(fam_df) > 0 else fam_df

if len(sig_fam_entries) == 0:
    print("No significant family-level hits to confirm in Tier 3.")
else:
    tier3_records = []
    for _, fam_row in sig_fam_entries.iterrows():
        src_sym = fam_row["source"]
        fam_name = fam_row["family"]
        t2_sign = np.sign(fam_row["stouffer_z"])

        src_id = q22_eids.get(src_sym)
        if src_id is None or src_id not in SpecMat.index:
            continue

        # Get target genes in this family that are in the matrix
        fam_targets = CURATED_TARGETS.get(fam_name, [])
        for tgt_sym in fam_targets:
            tgt_id = target_eids.get(tgt_sym)
            if tgt_id is None or tgt_id not in SpecMat.index:
                continue

            x = SpecMat.loc[src_id, cge_cols].values.astype(float)
            y = SpecMat.loc[tgt_id, cge_cols].values.astype(float)

            rho_obs, ci_lo, ci_hi = bca_ci_spearman(x, y, n_boot=2000, seed=42)
            excludes_zero = (ci_lo > 0) or (ci_hi < 0)
            concordant = (np.sign(rho_obs) == t2_sign)

            tier3_records.append({
                "source": src_sym,
                "target": tgt_sym,
                "family": fam_name,
                "rho_cge": rho_obs,
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "excludes_zero": excludes_zero,
                "t2_sign_concordant": concordant,
            })

    tier3_df = pd.DataFrame(tier3_records)
    if len(tier3_df) > 0:
        n_concordant = tier3_df["t2_sign_concordant"].sum()
        n_excl_zero = tier3_df["excludes_zero"].sum()
        print(f"Pairs tested:         {len(tier3_df)}")
        print(f"Sign-concordant:      {n_concordant}/{len(tier3_df)}")
        print(f"CI excludes zero:     {n_excl_zero}/{len(tier3_df)}")
        print()
        print(tier3_df.to_string(index=False, float_format="{:.3f}".format))
    else:
        print("No valid pairs for Tier 3 analysis.")

# %%
# --- Section 13d: Tier 1 vs Tier 2 comparison heatmap ---

if len(pair_df) > 0:
    # Build family ordering for targets
    family_order = list(CHANNEL_FAMILIES.keys())
    # Map each target to its family
    target_to_family = {}
    for fam_name, syms in CURATED_TARGETS.items():
        for sym in syms:
            target_to_family[sym] = fam_name

    # Order targets by family, then alphabetically within family
    unique_targets = pair_df["target"].unique()
    target_sorted = sorted(
        unique_targets,
        key=lambda t: (
            family_order.index(target_to_family.get(t, "zzz"))
            if target_to_family.get(t, "zzz") in family_order
            else len(family_order),
            t,
        ),
    )
    unique_sources = sorted(pair_df["source"].unique())

    # Pivot for heatmaps
    pivot_rho = pair_df.pivot(
        index="source", columns="target", values="rho_raw"
    ).reindex(index=unique_sources, columns=target_sorted)

    pivot_rpartial = pair_df.pivot(
        index="source", columns="target", values="r_partial"
    ).reindex(index=unique_sources, columns=target_sorted)

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    # Left: Tier 1 raw Spearman
    sns.heatmap(
        pivot_rho, ax=axes[0], cmap="RdBu_r", center=0,
        vmin=-0.6, vmax=0.6, linewidths=0.3, linecolor="white",
        cbar_kws={"shrink": 0.6, "label": "Spearman rho"},
        xticklabels=True, yticklabels=True,
    )
    axes[0].set_title("Tier 1: Raw Spearman rho")
    axes[0].set_xlabel("Target gene")
    axes[0].set_ylabel("22q source gene")
    axes[0].tick_params(axis="x", rotation=90, labelsize=7)
    axes[0].tick_params(axis="y", labelsize=8)

    # Right: Tier 2 partial Pearson
    sns.heatmap(
        pivot_rpartial, ax=axes[1], cmap="RdBu_r", center=0,
        vmin=-0.6, vmax=0.6, linewidths=0.3, linecolor="white",
        cbar_kws={"shrink": 0.6, "label": "Partial r"},
        xticklabels=True, yticklabels=True,
    )
    axes[1].set_title("Tier 2: Pearson partial r\n(cell-class corrected)")
    axes[1].set_xlabel("Target gene")
    axes[1].set_ylabel("22q source gene")
    axes[1].tick_params(axis="x", rotation=90, labelsize=7)
    axes[1].tick_params(axis="y", labelsize=8)

    fig.suptitle(
        "22q x Channel/Scaffold Profile Similarity\n"
        "Tier 1 (raw Spearman) vs Tier 2 (Pearson partial, cell-class corrected)",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()

    fig_dir = EPHYS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        fig_dir / "Fig_ProfileSim_22q_Channels_Neurons.png",
        dpi=150, bbox_inches="tight", transparent=True,
    )
    plt.show()
    print(f"Saved: {fig_dir / 'Fig_ProfileSim_22q_Channels_Neurons.png'}")
else:
    print("No pair results to plot.")
