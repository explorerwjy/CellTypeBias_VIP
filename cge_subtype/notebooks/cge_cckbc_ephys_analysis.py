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
#     display_name: hh_sbi_py311
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # CGE Interneuron Subtype Electrophysiology Analysis
#
# The scANVI mapping (96.3% interneuron purity) identified several CGE interneuron
# subgroups with distinct electrophysiological profiles. This notebook documents:
#
# 1. **Summary statistics** — per-cluster ephys feature distributions
# 2. **Box/violin plots** — key features compared across cell type groups
# 3. **Statistical tests** — Mann-Whitney U between CCKBC candidates and other groups
# 4. **Voltage traces** — representative hero sweeps from DANDI 000636 recordings
#
# ### Cell Type Groupings
#
# | Group | Clusters | Hypothesis |
# |-------|----------|------------|
# | CCKBC candidates | 284, 289 | CCK+ basket cells — fast-spiking, low UD ratio |
# | High-rate VIP | 290 | Fast-adapting VIP interneurons |
# | Irregular VIP | 291, 292 | Irregularly spiking VIP cells |
# | Regular VIP | 279, 280 | Regular-spiking VIP interneurons |
# | LAMP5 | 278, 287, 288 | Late-spiking / neurogliaform-like |

# %% tags=["parameters"]
# ── Configuration ──────────────────────────────────────────────────────────────
SAVE_FIGURES = True
FIGURE_DPI = 200
FIGURE_FORMAT = "pdf"

# %%
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu

# Project paths
PROJECT_ROOT = Path("/home/jw3514/Work/NeurSim/TransEphys")
ATLAS_DIR = PROJECT_ROOT / "atlas_matching"
RESULTS_DIR = ATLAS_DIR / "results" / "human_interneuron"

# External data paths
EPHYS_FX_PATH = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")
CELL_META_PATH = Path("/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/cell_metadata.csv")
DANDI_BASE = Path("/mnt/data0/DANDI/Processed/000636")
DANDI_SUMMARY = DANDI_BASE / "summary" / "all_analysis.csv"

# Figure output
FIGURE_DIR = RESULTS_DIR / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Cluster groupings
CLUSTER_GROUPS = {
    "CCKBC": [284, 289],
    "High-rate VIP": [290],
    "Irregular VIP": [291, 292],
    "Regular VIP": [279, 280],
    "LAMP5": [278, 287, 288],
}
GROUP_ORDER = ["CCKBC", "High-rate VIP", "Irregular VIP", "Regular VIP", "LAMP5"]
GROUP_COLORS = {
    "CCKBC": "#d62728",
    "High-rate VIP": "#ff7f0e",
    "Irregular VIP": "#2ca02c",
    "Regular VIP": "#1f77b4",
    "LAMP5": "#9467bd",
}

# %% [markdown]
# ## Load scANVI Mapping Results + Ephys Features

# %%
# Load mapping results
mapping = pd.read_csv(
    RESULTS_DIR / "scvi_mapped" / "mapping_results.csv", index_col=0
)
mapping.index = mapping.index.astype(int)
print(f"Mapping results: {mapping.shape[0]} cells")

# Load ephys features
ephys = pd.read_csv(EPHYS_FX_PATH, index_col=0)
ephys.index = ephys.index.astype(int)
print(f"Ephys features: {ephys.shape[0]} cells, {ephys.shape[1]} features")

# Filter to CGE interneurons
cge_mask = mapping["assigned_supercluster"] == "CGE interneuron"
cge_mapping = mapping[cge_mask]
print(f"CGE interneurons: {len(cge_mapping)}")

# Join mapping with ephys
shared_ids = cge_mapping.index.intersection(ephys.index)
print(f"CGE cells with ephys features: {len(shared_ids)}")

df = cge_mapping.loc[shared_ids].copy()
df = df.join(ephys.loc[shared_ids])

# Add group labels
cluster_to_group = {}
for group, clusters in CLUSTER_GROUPS.items():
    for c in clusters:
        cluster_to_group[c] = group

df["group"] = df["assigned_cluster"].map(cluster_to_group)
df_grouped = df.dropna(subset=["group"])
print(f"CGE cells in defined groups: {len(df_grouped)}")
print(f"\nGroup distribution:")
print(df_grouped["group"].value_counts().reindex(GROUP_ORDER))

# %% [markdown]
# ## Build specimen_id → Session Directory Mapping
#
# Links cell metadata (`dandi_stem`) to DANDI 000636 processed session directories
# via the `all_analysis.csv` stems.

# %%
# Load cell metadata for dandi_stem
cell_meta = pd.read_csv(CELL_META_PATH)
cell_meta["specimen_id"] = cell_meta["specimen_id"].astype(int)
cell_meta = cell_meta.set_index("specimen_id")
print(f"Cell metadata: {len(cell_meta)} cells, {cell_meta['dandi_stem'].notna().sum()} with dandi_stem")

# Load DANDI session stems
aa = pd.read_csv(DANDI_SUMMARY, usecols=["subject", "stem"])
aa = aa.drop_duplicates(subset="stem")
print(f"DANDI sessions: {len(aa)}")

# Map dandi_stem (suffix) → full stem (sub-xxx_ses-xxx_icephys)
# dandi_stem format: "ses-636982248_icephys"
# DANDI stem format: "sub-1000181910_ses-1000781418_icephys"
suffix_to_full = {}
for _, row in aa.iterrows():
    # Extract the ses-..._icephys suffix from the full stem
    parts = row["stem"].split("_", 1)  # split on first underscore: ["sub-xxx", "ses-xxx_icephys"]
    if len(parts) == 2:
        suffix_to_full[parts[1]] = (row["subject"], row["stem"])

# Build specimen_id → session directory path
session_dir_map = {}
cge_ids = set(df_grouped.index)

for spec_id in cge_ids:
    if spec_id not in cell_meta.index:
        continue
    ds = cell_meta.loc[spec_id, "dandi_stem"]
    if pd.isna(ds):
        continue
    if ds in suffix_to_full:
        subject, stem = suffix_to_full[ds]
        session_dir = DANDI_BASE / subject / stem
        if session_dir.is_dir():
            session_dir_map[spec_id] = session_dir

print(f"CGE grouped cells with DANDI session dirs: {len(session_dir_map)} / {len(cge_ids)}")

# Per-group availability
for group in GROUP_ORDER:
    group_ids = set(df_grouped[df_grouped["group"] == group].index)
    n_mapped = len(group_ids & set(session_dir_map.keys()))
    print(f"  {group}: {n_mapped}/{len(group_ids)} cells with voltage traces")

# %% [markdown]
# ## Summary Statistics Table

# %%
# Key ephys features to summarize
# (column_name, display_label, scale_factor)
SUMMARY_FEATURES = [
    ("avg_rate_hero", "Firing rate (Hz)", 1),
    ("upstroke_downstroke_ratio_ramp", "UD ratio (ramp)", 1),
    ("width_ramp", "AP width (ms, ramp)", 1000),  # s → ms
    ("rheobase_i", "Rheobase (pA)", 1),
    ("isi_cv_hero", "ISI CV (hero)", 1),
    ("adapt_hero", "Adaptation (hero)", 1),
    ("sag", "Sag ratio", 1),
    ("tau", "Tau (s)", 1),
    ("latency_hero", "Latency (s, hero)", 1),
    ("input_resistance", "Input R (MOhm)", 1),
]

# Build summary table
summary_rows = []
for group in GROUP_ORDER:
    gdf = df_grouped[df_grouped["group"] == group]
    row = {"Group": group, "n": len(gdf)}
    for col, label, scale in SUMMARY_FEATURES:
        vals = gdf[col].dropna() * scale
        if len(vals) > 0:
            row[label] = f"{vals.median():.2f} ({vals.quantile(0.25):.2f}-{vals.quantile(0.75):.2f})"
        else:
            row[label] = "N/A"
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows).set_index("Group")
print("Median (IQR) per group:")
summary_df

# %% [markdown]
# ## Box/Violin Plot Comparison

# %%
plot_features = {
    "avg_rate_hero": "Firing Rate (Hz)",
    "upstroke_downstroke_ratio_ramp": "Upstroke/Downstroke Ratio",
    "rheobase_i": "Rheobase (pA)",
    "isi_cv_hero": "ISI CV",
}

fig, axes = plt.subplots(1, 4, figsize=(16, 4))
fig.patch.set_alpha(0)

for ax, (col, label) in zip(axes, plot_features.items()):
    ax.patch.set_alpha(0)
    data = []
    positions = []
    colors = []
    for i, group in enumerate(GROUP_ORDER):
        vals = df_grouped.loc[df_grouped["group"] == group, col].dropna().values
        if len(vals) > 0:
            data.append(vals)
            positions.append(i)
            colors.append(GROUP_COLORS[group])

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker="o", markersize=3, alpha=0.5),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Overlay individual points
    for i, (group, pos) in enumerate(zip(GROUP_ORDER, positions)):
        vals = df_grouped.loc[df_grouped["group"] == group, col].dropna().values
        if len(vals) > 0:
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(
                pos + jitter, vals, c=colors[i], s=10, alpha=0.5, zorder=3, edgecolors="none"
            )

    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([g.replace(" ", "\n") for g in GROUP_ORDER], fontsize=8)
    ax.set_ylabel(label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle("CGE Interneuron Subtypes — Key Electrophysiology Features", fontsize=12, y=1.02)
plt.tight_layout()

if SAVE_FIGURES:
    fig.savefig(FIGURE_DIR / f"cge_ephys_boxplots.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## Statistical Tests (Mann-Whitney U)
#
# Compare CCKBC candidates against each other group on key features.

# %%
test_features = ["avg_rate_hero", "upstroke_downstroke_ratio_ramp", "rheobase_i", "isi_cv_hero"]
cckbc_df = df_grouped[df_grouped["group"] == "CCKBC"]

results = []
for feat in test_features:
    cckbc_vals = cckbc_df[feat].dropna().values
    if len(cckbc_vals) < 2:
        print(f"Warning: CCKBC has <2 non-null values for {feat}, skipping tests")
        continue
    for other_group in GROUP_ORDER[1:]:
        other_vals = df_grouped.loc[df_grouped["group"] == other_group, feat].dropna().values
        if len(other_vals) < 2:
            continue
        stat, pval = mannwhitneyu(cckbc_vals, other_vals, alternative="two-sided")
        results.append({
            "Feature": feat,
            "Comparison": f"CCKBC vs {other_group}",
            "CCKBC median": np.median(cckbc_vals),
            "Other median": np.median(other_vals),
            "U statistic": stat,
            "p-value": pval,
            "Significant (p<0.05)": pval < 0.05,
        })

if results:
    stats_df = pd.DataFrame(results)
    stats_df["p-value"] = stats_df["p-value"].map(lambda x: f"{x:.4f}" if x >= 0.001 else f"{x:.2e}")
    display(stats_df)
else:
    print("Insufficient CCKBC data for statistical tests (n<2).")

# %% [markdown]
# ## Helper: Load Hero Sweep Voltage Trace

# %%
def load_hero_sweep(specimen_id, session_dir_map):
    """Load the hero sweep (highest spike frequency) voltage and current traces.

    Returns
    -------
    dict with keys: time_ms, voltage_mV, current_pA, sweep_num, stimulus_start_ms, stimulus_end_ms
    or None if data unavailable.
    """
    if specimen_id not in session_dir_map:
        return None

    session_dir = session_dir_map[specimen_id]

    # Load analysis to find hero sweep
    analysis_path = session_dir / "analysis.parquet"
    if not analysis_path.exists():
        return None
    analysis = pd.read_parquet(analysis_path)

    # Hero = sweep with highest spike frequency
    if "spike_frequency_Hz" in analysis.columns:
        spiking = analysis.dropna(subset=["spike_frequency_Hz"])
        spiking = spiking[spiking["spike_frequency_Hz"] > 0]
        if len(spiking) == 0:
            return None
        hero_row = spiking.loc[spiking["spike_frequency_Hz"].idxmax()]
    else:
        return None

    hero_sweep = int(hero_row["sweep"])

    # Load sweep config for timing windows
    config_path = session_dir / "sweep_config.json"
    if not config_path.exists():
        return None
    with open(config_path) as f:
        sweep_cfg = json.load(f)

    sweep_key = str(hero_sweep)
    if sweep_key not in sweep_cfg["sweeps"]:
        return None
    windows = sweep_cfg["sweeps"][sweep_key]["windows"]

    # Find mV and pA parquet files
    mv_files = list(session_dir.glob("mV_*.parquet"))
    pa_files = list(session_dir.glob("pA_*.parquet"))
    if not mv_files or not pa_files:
        return None

    # Load only the hero sweep from mV
    mv = pd.read_parquet(mv_files[0])
    mv_sweep = mv[mv["sweep"] == hero_sweep]
    if len(mv_sweep) == 0:
        return None

    pa = pd.read_parquet(pa_files[0])
    pa_sweep = pa[pa["sweep"] == hero_sweep]

    # Extract time window: 0.1s before stimulus start to 0.1s after stimulus end
    stim_start = windows["stimulus_start_s"]
    stim_end = windows["stimulus_end_s"]
    t_pad = 0.1  # seconds padding
    t_start = stim_start - t_pad
    t_end = stim_end + t_pad

    mask_v = (mv_sweep["t_s"] >= t_start) & (mv_sweep["t_s"] <= t_end)
    mask_c = (pa_sweep["t_s"] >= t_start) & (pa_sweep["t_s"] <= t_end)

    v_data = mv_sweep.loc[mask_v]
    c_data = pa_sweep.loc[mask_c]

    if len(v_data) == 0:
        return None

    # Normalize time to start at 0 and convert to ms
    t0 = v_data["t_s"].iloc[0]
    time_ms = (v_data["t_s"].values - t0) * 1000
    voltage_mV = v_data["value"].values  # already in mV despite "volts" label

    c_time_ms = (c_data["t_s"].values - t0) * 1000 if len(c_data) > 0 else None
    current_pA = c_data["value"].values if len(c_data) > 0 else None  # already in pA

    # Downsample for plotting: target ~2kHz
    n_samples = len(time_ms)
    duration_s = (time_ms[-1] - time_ms[0]) / 1000
    current_rate = n_samples / duration_s if duration_s > 0 else 50000
    target_rate = 2000
    step = max(1, int(current_rate / target_rate))

    time_ms = time_ms[::step]
    voltage_mV = voltage_mV[::step]
    if c_time_ms is not None:
        c_time_ms = c_time_ms[::step]
        current_pA = current_pA[::step]

    return {
        "time_ms": time_ms,
        "voltage_mV": voltage_mV,
        "current_time_ms": c_time_ms,
        "current_pA": current_pA,
        "sweep_num": hero_sweep,
        "stimulus_start_ms": (stim_start - t0) * 1000,
        "stimulus_end_ms": (stim_end - t0) * 1000,
        "firing_rate": hero_row.get("spike_frequency_Hz", np.nan),
    }


# Quick test
test_id = list(session_dir_map.keys())[0] if session_dir_map else None
if test_id is not None:
    test_result = load_hero_sweep(test_id, session_dir_map)
    if test_result:
        print(f"Test load for specimen {test_id}: sweep {test_result['sweep_num']}, "
              f"{len(test_result['time_ms'])} samples, "
              f"rate={test_result['firing_rate']:.1f} Hz")
    else:
        print(f"Test load for specimen {test_id}: failed")

# %% [markdown]
# ## Helper: Plot Voltage Trace Panel

# %%
def plot_trace_panel(ax_v, ax_c, trace_data, title="", color="k"):
    """Plot voltage (top) and current (bottom) traces on given axes."""
    if trace_data is None:
        ax_v.text(0.5, 0.5, "No data", transform=ax_v.transAxes, ha="center", va="center")
        ax_c.text(0.5, 0.5, "No data", transform=ax_c.transAxes, ha="center", va="center")
        return

    ax_v.patch.set_alpha(0)
    ax_c.patch.set_alpha(0)

    # Stimulus shading
    ax_v.axvspan(
        trace_data["stimulus_start_ms"], trace_data["stimulus_end_ms"],
        alpha=0.08, color="gray", zorder=0,
    )
    ax_c.axvspan(
        trace_data["stimulus_start_ms"], trace_data["stimulus_end_ms"],
        alpha=0.08, color="gray", zorder=0,
    )

    # Voltage trace
    ax_v.plot(trace_data["time_ms"], trace_data["voltage_mV"], color=color, linewidth=0.5)
    ax_v.set_ylabel("mV")
    ax_v.set_title(title, fontsize=9)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)

    # Current trace
    if trace_data["current_pA"] is not None:
        ax_c.plot(trace_data["current_time_ms"], trace_data["current_pA"], color=color, linewidth=0.5)
    ax_c.set_ylabel("pA")
    ax_c.set_xlabel("Time (ms)")
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)

# %% [markdown]
# ## CCKBC Candidate Traces (Clusters 284, 289)

# %%
cckbc_ids = df_grouped[df_grouped["group"] == "CCKBC"].index.tolist()
cckbc_with_traces = [sid for sid in cckbc_ids if sid in session_dir_map]
n_plot = min(3, len(cckbc_with_traces))

if n_plot > 0:
    fig, axes = plt.subplots(2, n_plot, figsize=(5 * n_plot, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        sid = cckbc_with_traces[i]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"Specimen {sid}\nCluster {cluster}, {rate_str}"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS["CCKBC"])

    fig.suptitle("CCKBC Candidates — Hero Sweep Traces", fontsize=12, y=1.02)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(FIGURE_DIR / f"traces_cckbc.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.show()
else:
    print("No CCKBC cells with available voltage traces.")

# %% [markdown]
# ## High-Rate VIP Traces (Cluster 290)

# %%
hr_ids = df_grouped[df_grouped["group"] == "High-rate VIP"].index.tolist()
hr_with_traces = [sid for sid in hr_ids if sid in session_dir_map]
n_plot = min(3, len(hr_with_traces))

if n_plot > 0:
    fig, axes = plt.subplots(2, n_plot, figsize=(5 * n_plot, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        sid = hr_with_traces[i]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"Specimen {sid}\nCluster {cluster}, {rate_str}"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS["High-rate VIP"])

    fig.suptitle("High-Rate VIP — Hero Sweep Traces", fontsize=12, y=1.02)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(FIGURE_DIR / f"traces_highrate_vip.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.show()
else:
    print("No High-rate VIP cells with available voltage traces.")

# %% [markdown]
# ## Irregular VIP Traces (Clusters 291, 292)

# %%
irr_ids = df_grouped[df_grouped["group"] == "Irregular VIP"].index.tolist()
irr_with_traces = [sid for sid in irr_ids if sid in session_dir_map]
n_plot = min(3, len(irr_with_traces))

if n_plot > 0:
    fig, axes = plt.subplots(2, n_plot, figsize=(5 * n_plot, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        sid = irr_with_traces[i]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"Specimen {sid}\nCluster {cluster}, {rate_str}"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS["Irregular VIP"])

    fig.suptitle("Irregular VIP — Hero Sweep Traces", fontsize=12, y=1.02)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(FIGURE_DIR / f"traces_irregular_vip.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.show()
else:
    print("No Irregular VIP cells with available voltage traces.")

# %% [markdown]
# ## Regular VIP Traces (Clusters 279, 280)

# %%
reg_ids = df_grouped[df_grouped["group"] == "Regular VIP"].index.tolist()
reg_with_traces = [sid for sid in reg_ids if sid in session_dir_map]
n_plot = min(3, len(reg_with_traces))

if n_plot > 0:
    fig, axes = plt.subplots(2, n_plot, figsize=(5 * n_plot, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        sid = reg_with_traces[i]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"Specimen {sid}\nCluster {cluster}, {rate_str}"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS["Regular VIP"])

    fig.suptitle("Regular VIP — Hero Sweep Traces", fontsize=12, y=1.02)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(FIGURE_DIR / f"traces_regular_vip.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.show()
else:
    print("No Regular VIP cells with available voltage traces.")

# %% [markdown]
# ## LAMP5 Traces (Clusters 278, 287, 288)

# %%
lamp5_ids = df_grouped[df_grouped["group"] == "LAMP5"].index.tolist()
lamp5_with_traces = [sid for sid in lamp5_ids if sid in session_dir_map]
n_plot = min(3, len(lamp5_with_traces))

if n_plot > 0:
    fig, axes = plt.subplots(2, n_plot, figsize=(5 * n_plot, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        sid = lamp5_with_traces[i]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"Specimen {sid}\nCluster {cluster}, {rate_str}"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS["LAMP5"])

    fig.suptitle("LAMP5 — Hero Sweep Traces", fontsize=12, y=1.02)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(FIGURE_DIR / f"traces_lamp5.{FIGURE_FORMAT}", dpi=FIGURE_DPI, transparent=True, bbox_inches="tight")
    plt.show()
else:
    print("No LAMP5 cells with available voltage traces.")

# %% [markdown]
# ## Multi-Panel Comparison: One Cell per Group

# %%
# Pick one representative cell per group (with available traces)
rep_cells = {}
for group in GROUP_ORDER:
    group_ids = df_grouped[df_grouped["group"] == group].index.tolist()
    candidates = [sid for sid in group_ids if sid in session_dir_map]
    if candidates:
        # Pick the cell closest to the group median firing rate
        rates = df_grouped.loc[candidates, "avg_rate_hero"].dropna()
        if len(rates) > 0:
            median_rate = rates.median()
            best = rates.sub(median_rate).abs().idxmin()
            rep_cells[group] = best
        else:
            rep_cells[group] = candidates[0]

n_groups = len(rep_cells)
if n_groups > 0:
    fig, axes = plt.subplots(2, n_groups, figsize=(4.5 * n_groups, 4), squeeze=False)
    fig.patch.set_alpha(0)

    for i, group in enumerate(GROUP_ORDER):
        if group not in rep_cells:
            axes[0, i].text(0.5, 0.5, f"No {group}\ndata", transform=axes[0, i].transAxes, ha="center")
            continue
        sid = rep_cells[group]
        trace = load_hero_sweep(sid, session_dir_map)
        cluster = df_grouped.loc[sid, "assigned_cluster"]
        rate = df_grouped.loc[sid, "avg_rate_hero"]
        rate_str = f"{rate:.0f} Hz" if pd.notna(rate) else "N/A"
        title = f"{group}\n(cl. {cluster}, {rate_str})"
        plot_trace_panel(axes[0, i], axes[1, i], trace, title=title, color=GROUP_COLORS[group])

    fig.suptitle("CGE Interneuron Subtypes — Representative Hero Sweeps", fontsize=12, y=1.04)
    plt.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(
            FIGURE_DIR / f"traces_comparison.{FIGURE_FORMAT}",
            dpi=FIGURE_DPI, transparent=True, bbox_inches="tight",
        )
    plt.show()
else:
    print("No cells with available voltage traces for comparison plot.")

# %% [markdown]
# ## Summary
#
# Key observations:
# - **CCKBC candidates** (clusters 284, 289): n=5, expected to show fast-spiking phenotype
#   with narrow AP width and high UD ratio
# - **High-rate VIP** (cluster 290): largest group (n=45), fast-adapting
# - **Irregular VIP** (clusters 291, 292): n=23, irregular firing patterns (high ISI CV)
# - **Regular VIP** (clusters 279, 280): n=16, regular spiking with lower adaptation
# - **LAMP5** (clusters 278, 287, 288): n=60, late-spiking / neurogliaform phenotype
