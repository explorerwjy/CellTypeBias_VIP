# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Multi-Modal Cross-Species CCKBC/ISI VIP Analysis
#
# Integrates 4 datasets (mouse SC, mouse patch-seq, human SC, human patch-seq)
# to robustify CCKBC vs ISI VIP 22q bias analysis.
#
# **Key question:** Does the mouse finding (ISI VIP more affected than CCKBC by 22q deletion)
# replicate in human single-cell data?

# %%
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal, spearmanr

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path("..").resolve()
MAIN_ROOT = PROJECT_ROOT.parent
sys.path.insert(0, str(MAIN_ROOT))

# %% [markdown]
# ## 1. Load all data sources

# %%
# MetaNeighbor AUROC (Module A Pass 2)
auroc = pd.read_csv(PROJECT_ROOT / "results/cluster_bridge/metaneighbor_auroc_cge.csv", index_col=0)
print(f"MetaNeighbor AUROC: {auroc.shape} (mouse x human)")

# Previous convergent classification (Harmony, scVI, markers, ephys)
bias_summary = pd.read_csv(PROJECT_ROOT / "results/cckbc_convergent_bias_summary.csv")
conv_class = pd.read_csv(PROJECT_ROOT / "results/convergent_cckbc_classification.csv")
print(f"Convergent classification: {len(conv_class)} clusters")

# Updated 22q bias
updated_bias = pd.read_csv(PROJECT_ROOT / "results/updated_22q_bias/22q_bias_by_cluster.csv")
print(f"Updated bias: {len(updated_bias)} clusters")

# Ephys convergence
ephys_results = pd.read_csv(PROJECT_ROOT / "results/ephys_convergence/ephys_convergence_results.csv")
print(f"Ephys convergence: {len(ephys_results)} clusters")

# Allen mouse metadata (for subclass lookup)
mouse_meta = pd.read_csv(
    "/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/"
    "20230830/views/cell_metadata_with_cluster_annotation.csv",
    usecols=["cluster_alias", "subclass"],
)
cluster_to_subclass = mouse_meta.groupby("cluster_alias")["subclass"].first().to_dict()

# %% [markdown]
# ## 2. MetaNeighbor cross-species cluster correspondence
#
# For each human CGE cluster, find the best-matching mouse cluster and its subclass.

# %%
mn_results = []
for hc in sorted(auroc.columns, key=lambda x: int(float(x))):
    top5 = auroc[hc].nlargest(5)
    best_mc = int(float(top5.index[0]))
    best_auroc = top5.iloc[0]
    best_subclass = cluster_to_subclass.get(best_mc, "unknown")
    is_sncg = "Sncg" in best_subclass

    # Check if ALL top-5 are Sncg
    top5_subclasses = [cluster_to_subclass.get(int(float(m)), "") for m in top5.index]
    all_top5_sncg = all("Sncg" in s for s in top5_subclasses)

    mn_results.append({
        "human_cluster": int(float(hc)),
        "best_mouse_cluster": best_mc,
        "best_mouse_subclass": best_subclass,
        "best_auroc": best_auroc,
        "is_sncg": is_sncg,
        "all_top5_sncg": all_top5_sncg,
    })

mn_df = pd.DataFrame(mn_results)
print("MetaNeighbor top matches:")
print(mn_df[["human_cluster", "best_mouse_subclass", "best_auroc", "is_sncg", "all_top5_sncg"]].to_string(index=False))

# %% [markdown]
# ## 3. Integrated evidence table
#
# Combine MetaNeighbor (Module A), Harmony, scVI, marker expression,
# and 22q bias into one table per CGE cluster.

# %%
# Build master table
master = []
for _, row in bias_summary.iterrows():
    cid = int(row["Unnamed: 0"])
    mn_row = mn_df[mn_df["human_cluster"] == cid]
    bias_row = updated_bias[updated_bias["cluster_id"] == cid]

    entry = {
        "cluster_id": cid,
        "vip_status": "VIP+" if row["cluster_name"] == "INT-VIP" else "VIP-",
        # Harmony
        "harmony_cckbc_frac": row["cckbc_frac_harmony"],
        "harmony_imputed": row["is_imputed_cckbc"],
        # scVI
        "scvi_cckbc_frac": row["cckbc_frac_scvi"],
        "scvi_imputed": row["cckbc_frac_scvi"] > 0.5,
        # MetaNeighbor
        "mn_best_subclass": mn_row["best_mouse_subclass"].values[0] if len(mn_row) > 0 else "N/A",
        "mn_auroc": mn_row["best_auroc"].values[0] if len(mn_row) > 0 else np.nan,
        "mn_is_sncg": mn_row["is_sncg"].values[0] if len(mn_row) > 0 else False,
        "mn_all_top5_sncg": mn_row["all_top5_sncg"].values[0] if len(mn_row) > 0 else False,
        # Ephys
        "ephys_score": row["ephys_score"],
        "n_ephys": row["n_ephys"],
        # Bias
        "bias_22q_del": bias_row["22q_del"].values[0] if len(bias_row) > 0 else np.nan,
    }
    master.append(entry)

master_df = pd.DataFrame(master).sort_values("harmony_cckbc_frac", ascending=False)

# Consensus CCKBC: confirmed by >= 2 of (Harmony, scVI, MetaNeighbor)
master_df["n_methods_cckbc"] = (
    master_df["harmony_imputed"].astype(int)
    + master_df["scvi_imputed"].astype(int)
    + master_df["mn_is_sncg"].astype(int)
)
master_df["consensus"] = master_df["n_methods_cckbc"].map({
    0: "non-CCKBC", 1: "weak (1/3)", 2: "moderate (2/3)", 3: "strong (3/3)"
})

print("=== Integrated Evidence Table ===")
display_cols = [
    "cluster_id", "vip_status", "harmony_cckbc_frac", "scvi_cckbc_frac",
    "mn_is_sncg", "mn_auroc", "n_methods_cckbc", "consensus", "bias_22q_del",
]
print(master_df[display_cols].to_string(index=False))

# %% [markdown]
# ## 4. Three-way bias comparison
#
# Split CGE clusters into:
# 1. **VIP- putative-CCKBC** (277, 278): Harmony-imputed but NOT confirmed by MetaNeighbor
# 2. **VIP+ CCKBC** (279, 280, 281): Sncg-confirmed by MetaNeighbor
# 3. **VIP+ ISI** (remaining VIP+ clusters): non-CCKBC VIP interneurons

# %%
def assign_group(row):
    if row["harmony_imputed"] and row["vip_status"] == "VIP-":
        return "VIP- putative-CCKBC"
    elif row["mn_is_sncg"] and row["vip_status"] == "VIP+":
        return "VIP+ CCKBC (Sncg)"
    elif row["vip_status"] == "VIP+":
        return "VIP+ ISI"
    else:
        return "other"

master_df["group"] = master_df.apply(assign_group, axis=1)
print("Group assignments:")
print(master_df.groupby("group")["cluster_id"].apply(list).to_string())
print()
print("Group sizes:", master_df["group"].value_counts().to_dict())

# %%
# Statistical comparisons
groups = {}
for g in ["VIP- putative-CCKBC", "VIP+ CCKBC (Sncg)", "VIP+ ISI"]:
    vals = master_df.loc[master_df["group"] == g, "bias_22q_del"].dropna().values
    groups[g] = vals
    print(f"{g} (n={len(vals)}): mean={vals.mean():.4f}, median={np.median(vals):.4f}, values={[f'{v:.4f}' for v in sorted(vals)]}")

print()
# Pairwise Mann-Whitney
for g1, g2 in [
    ("VIP- putative-CCKBC", "VIP+ CCKBC (Sncg)"),
    ("VIP- putative-CCKBC", "VIP+ ISI"),
    ("VIP+ CCKBC (Sncg)", "VIP+ ISI"),
]:
    if len(groups[g1]) >= 2 and len(groups[g2]) >= 2:
        u, p = mannwhitneyu(groups[g1], groups[g2], alternative="two-sided")
        print(f"{g1} vs {g2}: U={u:.0f}, p={p:.4f}")

# Kruskal-Wallis
if all(len(v) >= 2 for v in groups.values()):
    h, p = kruskal(*groups.values())
    print(f"\nKruskal-Wallis (3-group): H={h:.3f}, p={p:.4f}")

# %% [markdown]
# ## 5. Concordance across mapping methods
#
# How well do Harmony, scVI, and MetaNeighbor agree on CCKBC classification?

# %%
concordance_df = master_df[["cluster_id", "vip_status", "harmony_imputed", "scvi_imputed",
                             "mn_is_sncg", "n_methods_cckbc", "consensus"]].copy()

# Method agreement matrix
print("=== Method Concordance ===")
print(f"Harmony-only CCKBC: {master_df['harmony_imputed'].sum():.0f} clusters")
print(f"scVI-only CCKBC:    {master_df['scvi_imputed'].sum():.0f} clusters")
print(f"MetaNeighbor Sncg:  {master_df['mn_is_sncg'].sum():.0f} clusters")
print(f"All 3 agree:        {(master_df['n_methods_cckbc'] == 3).sum():.0f} clusters -> {master_df.loc[master_df['n_methods_cckbc']==3, 'cluster_id'].tolist()}")
print(f"2/3 agree:          {(master_df['n_methods_cckbc'] == 2).sum():.0f} clusters -> {master_df.loc[master_df['n_methods_cckbc']==2, 'cluster_id'].tolist()}")
print(f"Harmony only (not confirmed): {master_df.loc[(master_df['harmony_imputed']) & (~master_df['scvi_imputed']) & (~master_df['mn_is_sncg']), 'cluster_id'].tolist()}")

# %% [markdown]
# ## 6. Visualization

# %%
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.patch.set_alpha(0)

# Panel A: 3-way bias boxplot
ax = axes[0]
group_order = ["VIP- putative-CCKBC", "VIP+ CCKBC (Sncg)", "VIP+ ISI"]
colors = ["#4ECDC4", "#FF6B6B", "#95A5A6"]
positions = [0, 1, 2]

for i, (g, c) in enumerate(zip(group_order, colors)):
    vals = groups[g]
    bp = ax.boxplot([vals], positions=[i], widths=0.5, patch_artist=True,
                    boxprops=dict(facecolor=c, alpha=0.7),
                    medianprops=dict(color="black", linewidth=2))
    ax.scatter(np.full_like(vals, i) + np.random.uniform(-0.1, 0.1, len(vals)),
               vals, c=c, edgecolors="black", s=50, zorder=3)

ax.set_xticks(positions)
ax.set_xticklabels(["VIP-\nputative-CCKBC", "VIP+\nCCKBC (Sncg)", "VIP+\nISI"], fontsize=9)
ax.set_ylabel("22q Deletion Bias")
ax.set_title("A. 22q Bias by VIP/CCKBC Group")
ax.patch.set_alpha(0)

# Significance annotations
# VIP- vs VIP+ ISI
if len(groups["VIP- putative-CCKBC"]) >= 2 and len(groups["VIP+ ISI"]) >= 2:
    _, p_val = mannwhitneyu(groups["VIP- putative-CCKBC"], groups["VIP+ ISI"])
    y_max = max(max(groups["VIP- putative-CCKBC"]), max(groups["VIP+ ISI"])) + 0.02
    ax.plot([0, 0, 2, 2], [y_max, y_max + 0.005, y_max + 0.005, y_max], "k-", lw=1)
    ax.text(1, y_max + 0.008, f"p={p_val:.3f}", ha="center", fontsize=9)

# Panel B: Concordance heatmap
ax = axes[1]
evidence_cols = ["harmony_imputed", "scvi_imputed", "mn_is_sncg"]
evidence_labels = ["Harmony", "scVI", "MetaNeighbor"]
sorted_df = master_df.sort_values(["n_methods_cckbc", "vip_status"], ascending=[False, True])
heatmap_data = sorted_df[evidence_cols].values.astype(float)

im = ax.imshow(heatmap_data.T, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(range(len(evidence_labels)))
ax.set_yticklabels(evidence_labels)
ax.set_xticks(range(len(sorted_df)))
ax.set_xticklabels(sorted_df["cluster_id"].values, rotation=90, fontsize=8)
ax.set_xlabel("Human CGE Cluster")
ax.set_title("B. CCKBC Evidence (Green=Yes)")

# Add VIP status bar
for i, (_, row) in enumerate(sorted_df.iterrows()):
    color = "#FF6B6B" if row["vip_status"] == "VIP+" else "#4ECDC4"
    ax.plot(i, -0.6, "s", color=color, markersize=6, clip_on=False)
ax.text(-1, -0.6, "VIP:", ha="right", va="center", fontsize=8)

ax.patch.set_alpha(0)

# Panel C: MetaNeighbor AUROC vs 22q bias
ax = axes[2]
for _, row in master_df.iterrows():
    if row["group"] == "VIP- putative-CCKBC":
        c, m = "#4ECDC4", "s"
    elif row["group"] == "VIP+ CCKBC (Sncg)":
        c, m = "#FF6B6B", "D"
    else:
        c, m = "#95A5A6", "o"
    ax.scatter(row["mn_auroc"], row["bias_22q_del"], c=c, marker=m,
               edgecolors="black", s=80, zorder=3)
    ax.annotate(str(row["cluster_id"]), (row["mn_auroc"], row["bias_22q_del"]),
                fontsize=7, ha="left", va="bottom")

ax.set_xlabel("MetaNeighbor Best AUROC")
ax.set_ylabel("22q Deletion Bias")
ax.set_title("C. AUROC vs 22q Bias")

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#4ECDC4", markersize=8, label="VIP- putative-CCKBC"),
    Line2D([0], [0], marker="D", color="w", markerfacecolor="#FF6B6B", markersize=8, label="VIP+ CCKBC (Sncg)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#95A5A6", markersize=8, label="VIP+ ISI"),
]
ax.legend(handles=legend_elements, fontsize=8, loc="lower right")
ax.patch.set_alpha(0)

plt.tight_layout()
plt.savefig(PROJECT_ROOT / "results" / "fig_cckbc_multimodal_analysis.pdf",
            transparent=True, bbox_inches="tight", dpi=150)
plt.show()
print("Saved: results/fig_cckbc_multimodal_analysis.pdf")

# %% [markdown]
# ## 7. Summary statistics for reviewer response

# %%
print("=" * 60)
print("SUMMARY FOR REVIEWER RESPONSE")
print("=" * 60)
print()
print("1. MetaNeighbor cross-species AUROC (CGE clusters):")
print(f"   - 589 mouse x 21 human clusters, 2000 HVGs")
print(f"   - Mean best-match AUROC: {mn_df['best_auroc'].mean():.3f}")
print(f"   - Clusters 279/280/281: ALL top-5 matches = Sncg Gaba (AUROC {mn_df[mn_df['human_cluster'].isin([279,280,281])]['best_auroc'].mean():.3f})")
print(f"   - Clusters 277/278: map to non-Sncg types (IA Mgp, OB Meis2)")
print()
print("2. Three-way 22q bias comparison:")
for g in group_order:
    vals = groups[g]
    print(f"   - {g} (n={len(vals)}): mean={vals.mean():.4f}")
_, p_kw = kruskal(*groups.values())
print(f"   - Kruskal-Wallis: p={p_kw:.4f}")
_, p_vn_vi = mannwhitneyu(groups["VIP- putative-CCKBC"], groups["VIP+ ISI"])
print(f"   - VIP- vs VIP+ ISI: p={p_vn_vi:.4f}")
print()
print("3. Method concordance:")
print(f"   - All 3 methods agree (CCKBC): clusters {master_df.loc[master_df['n_methods_cckbc']==3, 'cluster_id'].tolist()}")
print(f"   - Harmony-only (not confirmed): clusters {master_df.loc[(master_df['harmony_imputed']) & (master_df['n_methods_cckbc']==1), 'cluster_id'].tolist()}")
print()
print("4. Conclusion:")
print("   VIP expression status — not CCKBC identity — predicts 22q mutation bias.")
print("   The three Sncg-confirmed CCKBC homologs (279/280/281) are all VIP+ and show")
print("   bias comparable to VIP+ ISI. The apparent 'low-bias CCKBCs' (277/278) were")
print("   Harmony false positives corrected by MetaNeighbor.")
