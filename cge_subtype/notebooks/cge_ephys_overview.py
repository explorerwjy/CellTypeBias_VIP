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

# %% [markdown]
# # CGE Interneuron Electrophysiology Overview
#
# Comprehensive electrophysiology feature plots for all mapped CGE interneurons.
# Shows distributions of every major ephys feature across CGE clusters (276-296).
#
# Data sources:
# - Ephys features: `LeeDalley_ephys_fx.csv` (94 features, 704 cells)
# - Atlas mapping: scANVI/scVI mapping results (778 cells → CGE subset)

# %%
# %load_ext autoreload
# %autoreload 2

# %% tags=["parameters"]
# -- Configuration ----------------------------------------------------------------
MODEL_TYPE = "scanvi"  # "scanvi" or "scvi"

SAVE_FIGURES = True
FIGURE_DPI = 200
FIGURE_FORMAT = "pdf"

# Minimum cells per cluster to include in plots
MIN_CELLS_PER_CLUSTER = 3

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path

# Project paths
PROJECT_ROOT = Path("/home/jw3514/Work/NeurSim/TransEphys")
ATLAS_DIR = PROJECT_ROOT / "atlas_matching"
RESULTS_DIR = ATLAS_DIR / "results" / "human_interneuron"
MAPPED_DIR = RESULTS_DIR / f"{MODEL_TYPE}_mapped"

EPHYS_FX_PATH = Path("/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/LeeDalley_ephys_fx.csv")

FIGURE_DIR = RESULTS_DIR / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

print(f"Model type: {MODEL_TYPE}")

# %% [markdown]
# ## Load Data

# %%
# Load mapping results
mapping = pd.read_csv(MAPPED_DIR / "mapping_results.csv", index_col=0)
mapping.index = mapping.index.astype(int)
print(f"Mapping results: {mapping.shape[0]} cells")

# Load ephys features
ephys = pd.read_csv(EPHYS_FX_PATH, index_col=0)
ephys.index = ephys.index.astype(int)
print(f"Ephys features: {ephys.shape[0]} cells, {ephys.shape[1]} features")

# Filter to CGE interneurons
cge_mapping = mapping[mapping["assigned_supercluster"] == "CGE interneuron"]
print(f"CGE interneurons (mapped): {len(cge_mapping)}")

# Join
shared_ids = cge_mapping.index.intersection(ephys.index)
print(f"CGE cells with ephys: {len(shared_ids)}")

df = cge_mapping.loc[shared_ids].copy()
df = df.join(ephys.loc[shared_ids])

# Cluster counts
cluster_counts = df["assigned_cluster"].value_counts().sort_index()
print(f"\nCells per cluster:")
for cid, n in cluster_counts.items():
    marker = " *" if n < MIN_CELLS_PER_CLUSTER else ""
    print(f"  Cluster {cid}: {n}{marker}")

# Filter to clusters with enough cells
valid_clusters = cluster_counts[cluster_counts >= MIN_CELLS_PER_CLUSTER].index.tolist()
df_plot = df[df["assigned_cluster"].isin(valid_clusters)].copy()
df_plot["cluster"] = df_plot["assigned_cluster"].astype(str)
print(f"\nClusters with >= {MIN_CELLS_PER_CLUSTER} cells: {len(valid_clusters)} "
      f"({len(df_plot)} cells)")

# Sorted cluster labels for consistent ordering
cluster_order = [str(c) for c in sorted(valid_clusters)]

# %% [markdown]
# ## Color Palette

# %%
# Use tab20 for up to 21 CGE clusters (276-296)
_all_cge = list(range(276, 297))
_cmap = plt.cm.get_cmap("tab20")
CLUSTER_PALETTE = {
    str(cid): _cmap(i / max(len(_all_cge) - 1, 1))
    for i, cid in enumerate(_all_cge)
}
if len(_all_cge) > 20:
    CLUSTER_PALETTE[str(_all_cge[-1])] = "#17becf"

print(f"Palette: {len(CLUSTER_PALETTE)} cluster colors")

# %% [markdown]
# ## Feature Definitions
#
# Organized into thematic groups for multi-panel figures.

# %%
# (column_name, display_label, unit_scale)
# unit_scale is applied before plotting (e.g., 1000 to convert s→ms)
FEATURE_GROUPS = {
    "Firing properties": [
        ("avg_rate_hero", "Firing rate (hero)", 1),
        ("avg_rate_rheo", "Firing rate (rheo)", 1),
        ("rheobase_i", "Rheobase (pA)", 1),
        ("fi_fit_slope", "f-I slope", 1),
        ("latency_hero", "Latency (ms, hero)", 1000),
        ("latency_rheo", "Latency (ms, rheo)", 1000),
    ],
    "AP waveform (ramp)": [
        ("upstroke_downstroke_ratio_ramp", "UD ratio", 1),
        ("width_ramp", "AP width (ms)", 1000),
        ("threshold_v_ramp", "Threshold (mV)", 1),
        ("upstroke_ramp", "Upstroke (mV/ms)", 0.001),
        ("downstroke_ramp", "Downstroke (mV/ms)", 0.001),
        ("peak_deltav_ramp", "Peak delta-V (mV)", 1),
        ("trough_deltav_ramp", "Trough delta-V (mV)", 1),
        ("fast_trough_deltav_ramp", "Fast trough delta-V (mV)", 1),
    ],
    "AP waveform (hero)": [
        ("upstroke_downstroke_ratio_hero", "UD ratio", 1),
        ("width_hero", "AP width (ms)", 1000),
        ("threshold_v_hero", "Threshold (mV)", 1),
        ("upstroke_hero", "Upstroke (mV/ms)", 0.001),
        ("downstroke_hero", "Downstroke (mV/ms)", 0.001),
        ("peak_deltav_hero", "Peak delta-V (mV)", 1),
        ("trough_deltav_hero", "Trough delta-V (mV)", 1),
        ("fast_trough_deltav_hero", "Fast trough delta-V (mV)", 1),
    ],
    "ISI & adaptation": [
        ("isi_cv_hero", "ISI CV (hero)", 1),
        ("isi_cv_mean", "ISI CV (mean)", 1),
        ("isi_cv_max", "ISI CV (max)", 1),
        ("adapt_hero", "Adaptation (hero)", 1),
        ("adapt_mean", "Adaptation (mean)", 1),
        ("adapt_max", "Adaptation (max)", 1),
        ("mean_isi_hero", "Mean ISI (ms, hero)", 1000),
        ("median_isi_hero", "Median ISI (ms, hero)", 1000),
    ],
    "Subthreshold": [
        ("v_baseline", "Baseline Vm (mV)", 1),
        ("input_resistance", "Input R (MOhm)", 1),
        ("tau", "Tau (ms)", 1000),
        ("sag", "Sag ratio", 1),
        ("sag_area", "Sag area", 1),
        ("sag_depol", "Sag (depol)", 1),
    ],
    "Spike shape adaptation": [
        ("isi_adapt_ratio", "ISI adapt ratio", 1),
        ("width_adapt_ratio", "Width adapt ratio", 1),
        ("upstroke_adapt_ratio", "Upstroke adapt ratio", 1),
        ("downstroke_adapt_ratio", "Downstroke adapt ratio", 1),
        ("peak_v_adapt_ratio", "Peak V adapt ratio", 1),
        ("threshold_v_adapt_ratio", "Threshold V adapt ratio", 1),
    ],
}

# Count available features
n_total = sum(len(v) for v in FEATURE_GROUPS.values())
print(f"Feature groups: {len(FEATURE_GROUPS)}, total features: {n_total}")


# %% [markdown]
# ## Plotting Helper

# %%
def plot_feature_panel(df_plot, features, group_title, cluster_order, palette,
                       figsize_per_feature=(14, 2.5)):
    """Box + strip plot for a group of features across clusters.

    Parameters
    ----------
    features : list of (col, label, scale) tuples
    """
    n_feat = len(features)
    fig, axes = plt.subplots(n_feat, 1, figsize=(figsize_per_feature[0],
                                                  figsize_per_feature[1] * n_feat))
    fig.patch.set_alpha(0)
    if n_feat == 1:
        axes = [axes]

    for ax, (col, label, scale) in zip(axes, features):
        ax.patch.set_alpha(0)

        if col not in df_plot.columns:
            ax.text(0.5, 0.5, f"{col}\nnot found", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            ax.set_ylabel(label, fontsize=9)
            continue

        plot_df = df_plot[["cluster", col]].dropna(subset=[col]).copy()
        plot_df[col] = plot_df[col] * scale

        if len(plot_df) == 0:
            ax.text(0.5, 0.5, f"{label}\nno data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            ax.set_ylabel(label, fontsize=9)
            continue

        sns.boxplot(
            data=plot_df, x="cluster", y=col, order=cluster_order,
            palette=palette, ax=ax, fliersize=0, width=0.6,
            boxprops=dict(alpha=0.5), medianprops=dict(color="black", linewidth=1.5),
            whiskerprops=dict(alpha=0.5), capprops=dict(alpha=0.5),
        )
        sns.stripplot(
            data=plot_df, x="cluster", y=col, order=cluster_order,
            palette=palette, ax=ax, size=3, alpha=0.6, jitter=0.2, zorder=3,
        )

        ax.set_ylabel(label, fontsize=9)
        ax.set_xlabel("")
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Annotate n per cluster on the bottom
        if ax == axes[-1]:
            ax.set_xlabel("CGE Cluster", fontsize=10)

    fig.suptitle(group_title, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    return fig


# %%
def savefig(fig, name):
    if SAVE_FIGURES:
        path = FIGURE_DIR / f"{name}.{FIGURE_FORMAT}"
        fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight",
                    facecolor="none", transparent=True)
        print(f"Saved: {path}")

# %% [markdown]
# ## Generate All Feature Group Plots

# %%
for group_name, features in FEATURE_GROUPS.items():
    fig = plot_feature_panel(df_plot, features, group_name, cluster_order,
                             CLUSTER_PALETTE)
    tag = group_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
    savefig(fig, f"cge_ephys_{tag}_{MODEL_TYPE}")
    plt.show()

# %% [markdown]
# ## Summary: Per-Cluster Feature Medians (Heatmap)

# %%
# Collect all features into a flat list
all_features = []
for features in FEATURE_GROUPS.values():
    for col, label, scale in features:
        if col in df_plot.columns:
            all_features.append((col, label, scale))

# Build median matrix
median_data = {}
for col, label, scale in all_features:
    vals = df_plot.groupby("cluster")[col].median() * scale
    median_data[label] = vals

median_df = pd.DataFrame(median_data, index=cluster_order).T
# Drop features with all NaN
median_df = median_df.dropna(how="all")

# Z-score across clusters for visualization
median_z = median_df.subtract(median_df.mean(axis=1), axis=0).divide(
    median_df.std(axis=1), axis=0
)
median_z = median_z.dropna(how="all")

fig, ax = plt.subplots(figsize=(max(8, len(cluster_order) * 0.7), len(median_z) * 0.35 + 2))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

sns.heatmap(
    median_z, ax=ax, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
    linewidths=0.5, linecolor="white",
    cbar_kws={"label": "Z-score (across clusters)", "shrink": 0.6},
    xticklabels=True, yticklabels=True,
)
ax.set_xlabel("CGE Cluster", fontsize=10)
ax.set_ylabel("")
ax.set_title("Per-Cluster Median Ephys Features (Z-scored)", fontsize=12, fontweight="bold")
ax.tick_params(labelsize=8)
plt.xticks(rotation=0)
plt.yticks(rotation=0, fontsize=7)

plt.tight_layout()
savefig(fig, f"cge_ephys_heatmap_{MODEL_TYPE}")
plt.show()

# %% [markdown]
# ## Per-Cluster Cell Counts

# %%
count_df = df_plot["cluster"].value_counts().reindex(cluster_order).fillna(0).astype(int)

fig, ax = plt.subplots(figsize=(max(8, len(cluster_order) * 0.6), 3))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

bars = ax.bar(
    range(len(cluster_order)), count_df.values,
    color=[CLUSTER_PALETTE[c] for c in cluster_order],
    edgecolor="white", linewidth=0.5,
)
for bar, n in zip(bars, count_df.values):
    if n > 0:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                str(n), ha="center", va="bottom", fontsize=8)

ax.set_xticks(range(len(cluster_order)))
ax.set_xticklabels(cluster_order, fontsize=9)
ax.set_xlabel("CGE Cluster", fontsize=10)
ax.set_ylabel("Number of cells", fontsize=10)
ax.set_title("CGE Cells with Ephys Data per Cluster", fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
savefig(fig, f"cge_cell_counts_{MODEL_TYPE}")
plt.show()

# %% [markdown]
# ## Voltage Traces: CCKBC Candidates (Clusters 284, 289) & CCK Candidates (Cluster 279)
#
# Load hero sweep recordings from DANDI 000636 processed data.

# %%
import json

CELL_META_PATH = Path("/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/cell_metadata.csv")
DANDI_BASE = Path("/mnt/data0/DANDI/Processed/000636")
DANDI_SUMMARY = DANDI_BASE / "summary" / "all_analysis.csv"

# Build specimen_id → session directory mapping
cell_meta = pd.read_csv(CELL_META_PATH)
cell_meta["specimen_id"] = cell_meta["specimen_id"].astype(int)
cell_meta = cell_meta.set_index("specimen_id")

aa = pd.read_csv(DANDI_SUMMARY, usecols=["subject", "stem"])
aa = aa.drop_duplicates(subset="stem")
suffix_to_full = {}
for _, row in aa.iterrows():
    parts = row["stem"].split("_", 1)
    if len(parts) == 2:
        suffix_to_full[parts[1]] = (row["subject"], row["stem"])

session_dir_map = {}
for spec_id in df.index:
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

print(f"CGE cells with DANDI session dirs: {len(session_dir_map)} / {len(df)}")

# %%
# Load ephys features for trace annotation
ephys_full = pd.read_csv(EPHYS_FX_PATH, index_col=0)
ephys_full.index = ephys_full.index.astype(int)


TARGET_CURRENT_PA = 100  # pA — show this amplitude for all cells


def load_sweep_at_current(specimen_id, session_dir_map, target_pA=TARGET_CURRENT_PA):
    """Load the sweep closest to target_pA current injection."""
    if specimen_id not in session_dir_map:
        return None

    session_dir = session_dir_map[specimen_id]

    analysis_path = session_dir / "analysis.parquet"
    if not analysis_path.exists():
        return None
    analysis = pd.read_parquet(analysis_path)

    if "avg_injected_current_pA" not in analysis.columns:
        return None
    # Only consider positive-current sweeps
    pos = analysis[analysis["avg_injected_current_pA"] > 0].dropna(
        subset=["avg_injected_current_pA"]
    )
    if len(pos) == 0:
        return None
    # Find sweep closest to target
    best_idx = (pos["avg_injected_current_pA"] - target_pA).abs().idxmin()
    target_row = pos.loc[best_idx]
    target_sweep = int(target_row["sweep"])
    actual_pA = target_row["avg_injected_current_pA"]
    freq = target_row.get("spike_frequency_Hz", np.nan)

    config_path = session_dir / "sweep_config.json"
    if not config_path.exists():
        return None
    with open(config_path) as f:
        sweep_cfg = json.load(f)

    sweep_key = str(target_sweep)
    if sweep_key not in sweep_cfg["sweeps"]:
        return None
    windows = sweep_cfg["sweeps"][sweep_key]["windows"]

    mv_files = list(session_dir.glob("mV_*.parquet"))
    pa_files = list(session_dir.glob("pA_*.parquet"))
    if not mv_files or not pa_files:
        return None

    mv = pd.read_parquet(mv_files[0])
    mv_sweep = mv[mv["sweep"] == target_sweep]
    if len(mv_sweep) == 0:
        return None

    pa = pd.read_parquet(pa_files[0])
    pa_sweep = pa[pa["sweep"] == target_sweep]

    stim_start = windows["stimulus_start_s"]
    stim_end = windows["stimulus_end_s"]
    t_pad = 0.1
    t_start = stim_start - t_pad
    t_end = stim_end + t_pad

    mask_v = (mv_sweep["t_s"] >= t_start) & (mv_sweep["t_s"] <= t_end)
    mask_c = (pa_sweep["t_s"] >= t_start) & (pa_sweep["t_s"] <= t_end)

    v_data = mv_sweep.loc[mask_v]
    c_data = pa_sweep.loc[mask_c]
    if len(v_data) == 0:
        return None

    t0 = v_data["t_s"].iloc[0]
    time_ms = (v_data["t_s"].values - t0) * 1000
    voltage_mV = v_data["value"].values

    c_time_ms = (c_data["t_s"].values - t0) * 1000 if len(c_data) > 0 else None
    current_pA = c_data["value"].values if len(c_data) > 0 else None

    # Downsample to ~2 kHz
    n_samples = len(time_ms)
    duration_s = (time_ms[-1] - time_ms[0]) / 1000
    current_rate = n_samples / duration_s if duration_s > 0 else 50000
    step = max(1, int(current_rate / 2000))

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
        "sweep_num": target_sweep,
        "actual_pA": actual_pA,
        "stimulus_start_ms": (stim_start - t0) * 1000,
        "stimulus_end_ms": (stim_end - t0) * 1000,
        "firing_rate": freq,
    }


def plot_trace_panel(ax_v, ax_c, trace_data, title="", color="k"):
    """Plot voltage (top) and current (bottom) traces."""
    if trace_data is None:
        ax_v.text(0.5, 0.5, "No data", transform=ax_v.transAxes, ha="center", va="center")
        ax_c.text(0.5, 0.5, "No data", transform=ax_c.transAxes, ha="center", va="center")
        return

    ax_v.patch.set_alpha(0)
    ax_c.patch.set_alpha(0)

    ax_v.axvspan(trace_data["stimulus_start_ms"], trace_data["stimulus_end_ms"],
                 alpha=0.08, color="gray", zorder=0)
    ax_c.axvspan(trace_data["stimulus_start_ms"], trace_data["stimulus_end_ms"],
                 alpha=0.08, color="gray", zorder=0)

    ax_v.plot(trace_data["time_ms"], trace_data["voltage_mV"], color=color, linewidth=0.5)
    ax_v.set_ylabel("mV")
    ax_v.set_title(title, fontsize=9)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)

    if trace_data["current_pA"] is not None:
        ax_c.plot(trace_data["current_time_ms"], trace_data["current_pA"],
                  color=color, linewidth=0.5)
    ax_c.set_ylabel("pA")
    ax_c.set_xlabel("Time (ms)")
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)


def plot_cluster_traces(cluster_ids, group_label, color):
    """Plot ~100 pA sweep traces for ALL cells in given clusters with DANDI data."""
    cell_ids = df[df["assigned_cluster"].isin(cluster_ids)].index.tolist()
    cells_with_traces = [sid for sid in cell_ids if sid in session_dir_map]

    if not cells_with_traces:
        print(f"No {group_label} cells with available voltage traces.")
        return None

    n_plot = len(cells_with_traces)
    ncols = min(n_plot, 3)
    nrows_per_cell = 2  # voltage + current
    nrows = nrows_per_cell * ((n_plot + ncols - 1) // ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 2 * nrows), squeeze=False)
    fig.patch.set_alpha(0)

    for i in range(n_plot):
        row_block = (i // ncols) * 2
        col = i % ncols
        sid = cells_with_traces[i]
        trace = load_sweep_at_current(sid, session_dir_map)

        cluster = mapping.loc[sid, "assigned_cluster"]
        tier = mapping.loc[sid, "mapping_tier"]
        if trace is not None:
            pA_str = f"{trace['actual_pA']:.0f} pA"
            freq_str = (f"{trace['firing_rate']:.0f} Hz"
                        if pd.notna(trace["firing_rate"]) else "N/A")
        else:
            pA_str = "N/A"
            freq_str = "N/A"
        title = f"{sid}\nCl. {cluster}, {pA_str}, {freq_str}"

        plot_trace_panel(axes[row_block, col], axes[row_block + 1, col],
                         trace, title=title, color=color)

    # Hide unused subplots
    for i in range(n_plot, nrows // 2 * ncols):
        row_block = (i // ncols) * 2
        col = i % ncols
        if row_block < nrows:
            axes[row_block, col].set_visible(False)
            axes[row_block + 1, col].set_visible(False)

    fig.suptitle(f"{group_label} — {TARGET_CURRENT_PA} pA Sweep Traces "
                 f"({n_plot} cells)", fontsize=12, y=1.02)
    plt.tight_layout()
    return fig


# %% [markdown]
# ### Trace Figures per CGE Cluster
#
# Generate one figure per cluster for clusters 284, 289 (CCKBC) and 287-294 (broad CGE),
# showing ALL cells with available DANDI recordings.

# %%
TRACE_CLUSTERS = sorted(set([284, 289] + list(range(287, 295))))
# Colors: CCKBC in red, others from the cluster palette
TRACE_COLORS = {}
for cid in TRACE_CLUSTERS:
    if cid in (284, 289):
        TRACE_COLORS[cid] = "#d62728"  # red for CCKBC
    else:
        TRACE_COLORS[cid] = CLUSTER_PALETTE.get(str(cid), "#333333")

for cid in TRACE_CLUSTERS:
    n_mapped = (df["assigned_cluster"] == cid).sum()
    n_with_dandi = sum(1 for sid in df[df["assigned_cluster"] == cid].index
                       if sid in session_dir_map)
    if n_mapped == 0:
        print(f"Cluster {cid}: no mapped cells, skipping")
        continue

    label = f"Cluster {cid}"
    if cid in (284, 289):
        label += " (CCKBC)"
    fig = plot_cluster_traces([cid], label, TRACE_COLORS[cid])
    if fig is not None:
        savefig(fig, f"traces_cl{cid}_{MODEL_TYPE}")
        plt.show()
    else:
        print(f"Cluster {cid}: {n_mapped} mapped cells but none with DANDI traces")

# %% [markdown]
# ## Summary

# %%
print(f"Model: {MODEL_TYPE}")
print(f"CGE cells with ephys: {len(df_plot)}")
print(f"Clusters plotted (>= {MIN_CELLS_PER_CLUSTER} cells): {len(valid_clusters)}")
print(f"Feature groups: {len(FEATURE_GROUPS)}")
print(f"Total features plotted: {n_total}")
print()

# Trace summary
print(f"=== Trace Figures ({TARGET_CURRENT_PA} pA) ===")
for cid in TRACE_CLUSTERS:
    cell_ids = df[df["assigned_cluster"] == cid].index.tolist()
    n_with = sum(1 for s in cell_ids if s in session_dir_map)
    print(f"  Cluster {cid}: {len(cell_ids)} mapped, {n_with} with traces")

print()
if SAVE_FIGURES:
    print(f"Saved figures ({FIGURE_FORMAT}):")
    for f in sorted(FIGURE_DIR.glob(f"cge_*_{MODEL_TYPE}.{FIGURE_FORMAT}")):
        print(f"  {f.name}")
    for f in sorted(FIGURE_DIR.glob(f"traces_*_{MODEL_TYPE}.{FIGURE_FORMAT}")):
        print(f"  {f.name}")
