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
