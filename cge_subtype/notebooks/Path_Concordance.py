# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Module C: Path Concordance Validation
#
# This notebook validates that two distinct mapping strategies produce
# **concordant cell-type assignments** for mouse Patch-seq cells.
#
# ## Two Paths
#
# | Path | Steps |
# |------|-------|
# | **Direct** | Mouse Patch-seq → Human SC (Harmony label transfer, Module B) |
# | **Indirect** | Mouse Patch-seq → Mouse SC (known labels) → Human SC (RBH/AUROC, Module A) |
#
# Concordance is measured at three levels of granularity:
# - **Cluster** — finest grain, full Allen cluster ID
# - **Subclass** — e.g. "VIP", "SST", "Sncg"
# - **Supercluster** — broad groupings across species
#
# High concordance (> 0.7) and positive Cohen's kappa (> 0.4) support the
# robustness of the CGE interneuron mapping reported in the manuscript.

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import yaml

# Project paths
BIAS_PROJECT = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
CGE_DIR = BIAS_PROJECT / "cge_subtype"
ATLAS_MATCHING = Path("/home/jw3514/Work/NeurSim/TransEphys/atlas_matching")

sys.path.insert(0, str(BIAS_PROJECT))

from cge_subtype.src.concordance import (
    build_indirect_path,
    compute_concordance,
    compute_cohens_kappa,
    concordance_at_levels,
)

# Transparent backgrounds
matplotlib.rcParams["savefig.transparent"] = True

# %% [markdown]
# ## 1. Load Module A output — RBH pairs (mouse SC → human SC)

# %%
# TODO: replace with actual output path once Module A pipeline has been run
# Expected path: ATLAS_MATCHING / "results" / "rbh" / "mouse_to_human_rbh.csv"
# Columns required: mouse_cluster, human_cluster, correlation, is_rbh
RBH_PATH = ATLAS_MATCHING / "results" / "rbh" / "mouse_to_human_rbh.csv"

print(f"Module A RBH path: {RBH_PATH}")
if RBH_PATH.exists():
    rbh_df = pd.read_csv(RBH_PATH)
    print(f"  Loaded: {len(rbh_df)} rows, {rbh_df['is_rbh'].sum()} RBH pairs")
    # Filter to confirmed RBH pairs only
    mouse_to_human = rbh_df[rbh_df["is_rbh"]][["mouse_cluster", "human_cluster"]].copy()
    print(f"  RBH pairs retained: {len(mouse_to_human)}")
else:
    print("  [TODO] RBH file not found — using placeholder for demonstration")
    mouse_to_human = pd.DataFrame({
        "mouse_cluster": ["Sncg_1", "Sncg_2", "Vip_1", "Vip_2"],
        "human_cluster": ["CGE_Lamp5_1", "CGE_Lamp5_2", "CGE_VIP_1", "CGE_VIP_2"],
    })

print(mouse_to_human.head())

# %% [markdown]
# ## 2. Load Module B output — Harmony predictions (patch-seq → human SC, direct path)

# %%
# TODO: replace with actual output path once Module B pipeline has been run
# Expected path: ATLAS_MATCHING / "results" / "mapped" / "patchseq_mapping_results.csv"
# Columns required: cell_id, predicted_cluster, confidence, n_neighbors_agreeing
HARMONY_PATH = ATLAS_MATCHING / "results" / "mapped" / "patchseq_mapping_results.csv"

print(f"Module B Harmony predictions path: {HARMONY_PATH}")
if HARMONY_PATH.exists():
    harmony_df = pd.read_csv(HARMONY_PATH)
    print(f"  Loaded: {len(harmony_df)} cells")
    print(f"  Columns: {list(harmony_df.columns)}")
else:
    print("  [TODO] Harmony predictions not found — using placeholder for demonstration")
    harmony_df = pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(8)],
        "predicted_cluster": [
            "CGE_Lamp5_1", "CGE_Lamp5_1", "CGE_Lamp5_2",
            "CGE_VIP_1", "CGE_VIP_1", "CGE_VIP_2",
            "CGE_Lamp5_1", "CGE_VIP_1",
        ],
        "confidence": [0.9, 0.8, 0.85, 0.7, 0.75, 0.9, 0.6, 0.8],
    })

print(harmony_df.head())

# %% [markdown]
# ## 3. Load existing direct mapping — patch-seq → mouse cluster assignments

# %%
# TODO: replace with actual path to patch-seq mouse cluster labels
# This could come from the ephys-based or transcriptomic subtype assignment.
# Expected: DataFrame with columns cell_id, mouse_cluster (and optionally subtype info)
PATCHSEQ_MOUSE_PATH = CGE_DIR / "results" / "patchseq_mouse_cluster_labels.csv"

print(f"Patch-seq → mouse cluster labels: {PATCHSEQ_MOUSE_PATH}")
if PATCHSEQ_MOUSE_PATH.exists():
    patchseq_mouse = pd.read_csv(PATCHSEQ_MOUSE_PATH)
    print(f"  Loaded: {len(patchseq_mouse)} cells")
    print(f"  Mouse clusters: {patchseq_mouse['mouse_cluster'].nunique()} unique")
else:
    print("  [TODO] Patch-seq mouse labels not found — using placeholder for demonstration")
    patchseq_mouse = pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(8)],
        "mouse_cluster": [
            "Sncg_1", "Sncg_1", "Sncg_2",
            "Vip_1", "Vip_1", "Vip_2",
            "Sncg_1", "Vip_1",
        ],
    })

print(patchseq_mouse.head())

# %% [markdown]
# ## 4. Build indirect path: patch-seq → mouse SC → human SC

# %%
patchseq_indirect = build_indirect_path(patchseq_mouse, mouse_to_human)
print(f"Indirect path built: {len(patchseq_indirect)} cells")
print(f"  Cells with indirect human cluster: {patchseq_indirect['indirect_human_cluster'].notna().sum()}")
print(f"  Cells with no match (NaN): {patchseq_indirect['indirect_human_cluster'].isna().sum()}")
print(patchseq_indirect.head())

# %% [markdown]
# ## 5. Merge direct and indirect paths

# %%
# Merge Harmony (direct) predictions onto the indirect path DataFrame
merged = patchseq_indirect.merge(
    harmony_df[["cell_id", "predicted_cluster"]].rename(
        columns={"predicted_cluster": "direct_human_cluster"}
    ),
    on="cell_id",
    how="inner",
)
print(f"Merged: {len(merged)} cells with both direct and indirect assignments")
print(merged[["cell_id", "mouse_cluster", "direct_human_cluster", "indirect_human_cluster"]].head(10))

# %% [markdown]
# ## 6. Concordance at multiple levels

# %%
# TODO: load actual cluster → subclass / supercluster mapping tables once available
# Expected source: Allen Brain Atlas annotation tables or Siletti 2023 metadata
# Example format:
#   cluster_to_subclass = {"CGE_Lamp5_1": "Lamp5", "CGE_VIP_1": "VIP", ...}
#   cluster_to_supercluster = {"CGE_Lamp5_1": "CGE", "CGE_VIP_1": "CGE", ...}

# Placeholder mapping consistent with placeholder data above
cluster_to_subclass = {
    "CGE_Lamp5_1": "Lamp5",
    "CGE_Lamp5_2": "Lamp5",
    "CGE_VIP_1": "VIP",
    "CGE_VIP_2": "VIP",
}
cluster_to_supercluster = {
    "CGE_Lamp5_1": "CGE",
    "CGE_Lamp5_2": "CGE",
    "CGE_VIP_1": "CGE",
    "CGE_VIP_2": "CGE",
}

# %%
concordance_results = concordance_at_levels(
    direct_clusters=merged["direct_human_cluster"],
    indirect_clusters=merged["indirect_human_cluster"],
    cluster_to_subclass=cluster_to_subclass,
    cluster_to_supercluster=cluster_to_supercluster,
)

print("Path concordance results:")
print("-" * 45)
for key, val in concordance_results.items():
    label = key.replace("_", " ").title()
    if isinstance(val, float) and not (val != val):  # not NaN
        print(f"  {label:<32} {val:.3f}")
    else:
        print(f"  {label:<32} NaN")

# %% [markdown]
# ## 7. CCKBC-focused subsection
#
# CCKBCs (cholecystokinin basket cells) are a key CGE interneuron subtype
# implicated in the manuscript's findings.  Here we restrict the concordance
# analysis to cells with CCKBC-associated mouse cluster assignments (Sncg
# subclass in the Allen Mouse Brain Atlas).

# %%
# Identify CCKBC cells (Sncg subclass in mouse)
# TODO: adjust the CCKBC_CLUSTERS list to match the actual mouse cluster names
# in patchseq_mouse_cluster_labels.csv once real data is available
CCKBC_CLUSTERS = ["Sncg_1", "Sncg_2"]

cckbc_mask = merged["mouse_cluster"].isin(CCKBC_CLUSTERS)
cckbc_merged = merged[cckbc_mask].copy()
print(f"CCKBC cells: {len(cckbc_merged)} (out of {len(merged)} total)")

if len(cckbc_merged) > 0:
    cckbc_concordance = concordance_at_levels(
        direct_clusters=cckbc_merged["direct_human_cluster"],
        indirect_clusters=cckbc_merged["indirect_human_cluster"],
        cluster_to_subclass=cluster_to_subclass,
        cluster_to_supercluster=cluster_to_supercluster,
    )
    print("\nCCKBC-specific concordance:")
    print("-" * 45)
    for key, val in cckbc_concordance.items():
        label = key.replace("_", " ").title()
        if isinstance(val, float) and not (val != val):
            print(f"  {label:<32} {val:.3f}")
        else:
            print(f"  {label:<32} NaN")
else:
    print("  No CCKBC cells found with current placeholder data.")

# %% [markdown]
# ## 8. Summary bar chart: concordance by level

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

levels = ["cluster", "subclass", "supercluster"]
concordance_vals = [concordance_results[f"{lv}_concordance"] for lv in levels]
kappa_vals = [concordance_results[f"{lv}_kappa"] for lv in levels]

# Replace NaN with 0 for plotting
concordance_plot = [v if not (v != v) else 0.0 for v in concordance_vals]
kappa_plot = [v if not (v != v) else 0.0 for v in kappa_vals]

ax_conc = axes[0]
ax_conc.bar(levels, concordance_plot, color="#4C72B0", alpha=0.85)
ax_conc.axhline(0.7, color="red", linestyle="--", linewidth=1.0, label="0.7 threshold")
ax_conc.set_ylim(0, 1.05)
ax_conc.set_ylabel("Concordance (fraction)")
ax_conc.set_title("Path Concordance by Level")
ax_conc.legend(fontsize=9)
ax_conc.set_facecolor("none")

ax_kappa = axes[1]
ax_kappa.bar(levels, kappa_plot, color="#DD8452", alpha=0.85)
ax_kappa.axhline(0.4, color="red", linestyle="--", linewidth=1.0, label="κ=0.4 threshold")
ax_kappa.set_ylim(-0.1, 1.05)
ax_kappa.set_ylabel("Cohen's κ")
ax_kappa.set_title("Cohen's Kappa by Level")
ax_kappa.legend(fontsize=9)
ax_kappa.set_facecolor("none")

fig.patch.set_alpha(0)
plt.tight_layout()

# TODO: update output path when figure directory is established
FIG_OUT = CGE_DIR / "results" / "figures" / "path_concordance_summary.pdf"
print(f"[TODO] Save figure to: {FIG_OUT}")
plt.show()

# %% [markdown]
# ## 9. Sankey diagram placeholder
#
# A Sankey diagram showing direct vs. indirect label flow would visually
# summarize where paths agree and disagree.
#
# **Implementation note**: use `matplotlib-sankey` or `plotly.graph_objects.Sankey`
# once real data is available.  Key inputs:
# - Source nodes: unique `direct_human_cluster` values
# - Target nodes: unique `indirect_human_cluster` values
# - Link values: cell counts per (direct, indirect) pair

# %%
# Sankey placeholder: cross-tabulation of direct vs indirect assignments
if len(merged) > 0:
    cross_tab = pd.crosstab(
        merged["direct_human_cluster"].fillna("unassigned"),
        merged["indirect_human_cluster"].fillna("unassigned"),
    )
    print("Cross-tabulation (direct rows × indirect columns):")
    print(cross_tab.to_string())
else:
    print("[TODO] No merged data available for Sankey diagram.")

# TODO: implement full Sankey diagram once real data is loaded
# Example:
# import plotly.graph_objects as go
# fig = go.Figure(go.Sankey(
#     node=dict(label=all_labels),
#     link=dict(source=sources, target=targets, value=values),
# ))
# fig.write_html(CGE_DIR / "results" / "figures" / "path_concordance_sankey.html")
print("\n[TODO] Implement Sankey diagram using plotly once real data is available.")
print(f"[TODO] Output: {CGE_DIR / 'results' / 'figures' / 'path_concordance_sankey.html'}")

# %% [markdown]
# ## 10. Export concordance table

# %%
# Build a summary DataFrame for downstream reporting
summary_rows = []
for lv in ["cluster", "subclass", "supercluster"]:
    summary_rows.append({
        "level": lv,
        "concordance": concordance_results[f"{lv}_concordance"],
        "cohens_kappa": concordance_results[f"{lv}_kappa"],
        "n_cells": merged[["direct_human_cluster", "indirect_human_cluster"]].notna().all(axis=1).sum(),
    })

concordance_summary = pd.DataFrame(summary_rows)
print(concordance_summary.to_string(index=False))

# TODO: save to results once real data is available
CONCORDANCE_OUT = CGE_DIR / "results" / "path_concordance_summary.csv"
print(f"\n[TODO] Save concordance table to: {CONCORDANCE_OUT}")
# concordance_summary.to_csv(CONCORDANCE_OUT, index=False)
