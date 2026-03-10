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
# # Cross-Species CCKBC Mapping: Mouse → Human
#
# Map mouse Patch-seq cells with known CCKBC labels (Sncg subclass) onto the
# Siletti et al. 2023 human brain atlas to identify which human CGE clusters
# correspond to mouse CCKBCs. Connect assignments to 22q mutation bias.
#
# **Pipeline outputs used:**
# - `patchseq_mapping_results.csv` — mouse Patch-seq → human cluster assignments
# - `reference_with_latent.h5ad` — cross-species reference with UMAP
# - Human Patch-seq mapping results (for ephys-based CCKBC validation)
# - 22q mutation bias results from CellTypeBias_VIP pipeline

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Patch
from scipy import stats
import yaml

# Paths
ATLAS_MATCHING = Path("/home/jw3514/Work/NeurSim/TransEphys/atlas_matching")
BIAS_PROJECT = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")

sys.path.insert(0, str(ATLAS_MATCHING))

# Load config
with open(ATLAS_MATCHING / "configs" / "cross_species_config.yaml") as f:
    config = yaml.safe_load(f)

RESULTS_DIR = ATLAS_MATCHING / config["output"]["subdirs"]["mapped"]
PREPROC_DIR = ATLAS_MATCHING / config["output"]["subdirs"]["preprocessed"]
FIG_DIR = ATLAS_MATCHING / config["output"]["subdirs"]["figures"]
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Transparent backgrounds for all figures
matplotlib.rcParams["savefig.transparent"] = True

# %% [markdown]
# ## 1. Load mapping results

# %%
# Mouse Patch-seq mapping results
mapping_df = pd.read_csv(RESULTS_DIR / "patchseq_mapping_results.csv", index_col=0)
print(f"Mouse Patch-seq mapping results: {len(mapping_df)} cells")
print(f"Columns: {list(mapping_df.columns)}")

# Cross-species reference with latent+UMAP
adata_ref = sc.read_h5ad(PREPROC_DIR / "reference_with_latent.h5ad")
print(f"Reference: {adata_ref.shape}")

# Mouse Patch-seq query with latent
adata_query = sc.read_h5ad(RESULTS_DIR / "patchseq_query_mapped.h5ad")
print(f"Query: {adata_query.shape}")

# %% [markdown]
# ## 2. Mapping quality assessment

# %%
# Confidence distribution
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
fig.patch.set_alpha(0)

axes[0].hist(mapping_df["conf_supercluster"], bins=50, color="steelblue", edgecolor="white")
axes[0].set_xlabel("Supercluster confidence")
axes[0].set_ylabel("Count")
axes[0].set_title("Supercluster confidence")
axes[0].patch.set_alpha(0)

axes[1].hist(mapping_df["conf_cluster"], bins=50, color="darkorange", edgecolor="white")
axes[1].set_xlabel("Cluster confidence")
axes[1].set_title("Cluster confidence")
axes[1].patch.set_alpha(0)

tier_counts = mapping_df["mapping_tier"].value_counts()
axes[2].barh(tier_counts.index, tier_counts.values, color="mediumpurple", edgecolor="white")
axes[2].set_xlabel("Count")
axes[2].set_title("Mapping tier")
axes[2].patch.set_alpha(0)

plt.tight_layout()
plt.savefig(FIG_DIR / "mapping_confidence.pdf", transparent=True)
plt.show()

# %%
# Species mixing in UMAP
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Color by species
species_colors = {"human": "#1f77b4", "mouse": "#ff7f0e"}
ref_species = adata_ref.obs["species"].astype(str).map(species_colors).values
axes[0].scatter(
    adata_ref.obsm["X_umap"][:, 0], adata_ref.obsm["X_umap"][:, 1],
    c=ref_species, s=0.5, alpha=0.3, rasterized=True,
)
axes[0].set_title("Reference: species")
axes[0].set_xlabel("UMAP1")
axes[0].set_ylabel("UMAP2")
axes[0].legend(handles=[Patch(color=c, label=s) for s, c in species_colors.items()])
axes[0].patch.set_alpha(0)

# Color by cell class
if "cell_class" in adata_ref.obs.columns:
    class_colors = {"GABAergic": "#e41a1c", "Glutamatergic": "#4daf4a", "NonNeuronal": "#984ea3"}
    ref_class = adata_ref.obs["cell_class"].astype(str).map(class_colors).fillna("gray").values
    axes[1].scatter(
        adata_ref.obsm["X_umap"][:, 0], adata_ref.obsm["X_umap"][:, 1],
        c=ref_class, s=0.5, alpha=0.3, rasterized=True,
    )
    axes[1].set_title("Reference: cell class")
    axes[1].set_xlabel("UMAP1")
    axes[1].legend(handles=[Patch(color=c, label=s) for s, c in class_colors.items()])
    axes[1].patch.set_alpha(0)

fig.patch.set_alpha(0)
plt.tight_layout()
plt.savefig(FIG_DIR / "umap_species_class.pdf", transparent=True, dpi=150)
plt.show()

# %%
# Sanity checks: mouse subclass -> human supercluster concordance
print("=" * 60)
print("SANITY CHECK: Mouse subclass -> Human supercluster")
print("=" * 60)

# Expected: mouse Vip -> CGE, mouse Pvalb -> MGE, mouse Sst -> MGE
for meta_col in ["mouse_subclass", "mouse_corresponding_AIT2.3.1_alias", "mouse_RNA type"]:
    if meta_col in mapping_df.columns:
        cross = pd.crosstab(
            mapping_df[meta_col],
            mapping_df["assigned_human_supercluster"],
            margins=True,
        )
        print(f"\n{meta_col} x assigned_human_supercluster:")
        print(cross.to_string())
        break

# %% [markdown]
# ## 3. CCKBC cluster assignment

# %%
# Which human CGE clusters receive mouse CCKBC cells?
cckbc_mask = mapping_df["is_cckbc"].astype(bool) if "is_cckbc" in mapping_df.columns else pd.Series(False, index=mapping_df.index)
cckbc_results = mapping_df[cckbc_mask]
non_cckbc_results = mapping_df[~cckbc_mask]

print(f"CCKBC cells: {len(cckbc_results)}")
print(f"Non-CCKBC cells: {len(non_cckbc_results)}")

# Focus on CGE-assigned cells
cge_mask = mapping_df["assigned_human_supercluster"] == "CGE interneuron"
cge_clusters = mapping_df[cge_mask]

print(f"\nAll CGE-assigned cells: {len(cge_clusters)}")
print(f"CCKBC -> CGE: {(cckbc_mask & cge_mask).sum()}")

# %%
# Compute CCKBC fraction per human CGE cluster
cge_cluster_ids = sorted(cge_clusters["assigned_human_cluster"].unique())
print(f"\nCGE clusters receiving assignments: {len(cge_cluster_ids)}")

cluster_stats = []
for clust in cge_cluster_ids:
    mask_clust = cge_clusters["assigned_human_cluster"] == clust
    n_total = mask_clust.sum()
    n_cckbc = (cge_clusters.loc[mask_clust].index.isin(cckbc_results.index)).sum()
    cluster_stats.append({
        "cluster_id": clust,
        "n_total_assigned": n_total,
        "n_cckbc": n_cckbc,
        "cckbc_fraction": n_cckbc / n_total if n_total > 0 else 0,
    })

cluster_stats_df = pd.DataFrame(cluster_stats).sort_values("cckbc_fraction", ascending=False)
print("\nCCKBC fraction per CGE cluster:")
print(cluster_stats_df.to_string(index=False))

# %%
# Barplot: CCKBC fraction per CGE cluster
fig, ax = plt.subplots(figsize=(12, 5))
fig.patch.set_alpha(0)
ax.patch.set_alpha(0)

x = range(len(cluster_stats_df))
bars = ax.bar(x, cluster_stats_df["cckbc_fraction"],
              color="firebrick", edgecolor="white", alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(cluster_stats_df["cluster_id"].values, rotation=45, ha="right")
ax.set_xlabel("Human CGE cluster ID")
ax.set_ylabel("CCKBC fraction")
ax.set_title("Fraction of mouse CCKBC cells mapping to each human CGE cluster")

# Add count labels
for i, row in enumerate(cluster_stats_df.itertuples()):
    ax.text(i, row.cckbc_fraction + 0.01, f"n={row.n_cckbc}", ha="center", fontsize=7)

plt.tight_layout()
plt.savefig(FIG_DIR / "cckbc_fraction_per_cluster.pdf", transparent=True)
plt.show()

# %% [markdown]
# ## 4. Sankey diagram: mouse subclass -> human supercluster/cluster

# %%
# Alluvial flow: mouse subclass -> human supercluster
try:
    import plotly.graph_objects as go

    # Prepare data for Sankey
    meta_col = None
    for col in ["mouse_subclass", "mouse_corresponding_AIT2.3.1_alias", "mouse_RNA type"]:
        if col in mapping_df.columns:
            meta_col = col
            break

    if meta_col:
        flow_df = mapping_df[[meta_col, "assigned_human_supercluster"]].dropna()
        flow_counts = flow_df.groupby([meta_col, "assigned_human_supercluster"]).size().reset_index(name="count")
        flow_counts = flow_counts[flow_counts["count"] >= 5]  # Filter small flows

        source_labels = sorted(flow_counts[meta_col].unique())
        target_labels = sorted(flow_counts["assigned_human_supercluster"].unique())
        all_labels = source_labels + target_labels

        source_idx = [all_labels.index(s) for s in flow_counts[meta_col]]
        target_idx = [all_labels.index(t) for t in flow_counts["assigned_human_supercluster"]]

        fig_sankey = go.Figure(data=[go.Sankey(
            node=dict(label=all_labels, pad=15, thickness=20),
            link=dict(
                source=source_idx,
                target=target_idx,
                value=flow_counts["count"].values,
            ),
        )])
        fig_sankey.update_layout(title="Mouse subclass -> Human supercluster", font_size=10)
        fig_sankey.write_html(str(FIG_DIR / "sankey_subclass_supercluster.html"))
        fig_sankey.show()
        print("Sankey diagram saved")
    else:
        print("No suitable mouse subclass column found for Sankey diagram")
except ImportError:
    print("plotly not available; skipping Sankey diagram")

# %% [markdown]
# ## 5. Validation with human Patch-seq (ephys-based CCKBC)

# %%
# Load human Patch-seq mapping results
human_mapping_path = ATLAS_MATCHING / config["validation"]["human_patchseq_mapping"]
if human_mapping_path.exists():
    human_mapping_df = pd.read_csv(human_mapping_path, index_col=0)
    print(f"Human Patch-seq mapping: {len(human_mapping_df)} cells")

    # Load ephys features
    ephys_path = config["validation"]["human_patchseq_ephys"]
    if Path(ephys_path).exists():
        ephys_df = pd.read_csv(ephys_path, index_col=0)
        print(f"Ephys features: {ephys_df.shape}")
    else:
        ephys_df = None
        print(f"Ephys features not found at {ephys_path}")
else:
    human_mapping_df = None
    ephys_df = None
    print(f"Human Patch-seq mapping not found at {human_mapping_path}")

# %%
# Ephys-based CCKBC clusters from cge_cckbc_ephys_analysis
CLUSTER_GROUPS = config["validation"]["cluster_groups"]

print("Ephys-defined cluster groups:")
for group, clusters in CLUSTER_GROUPS.items():
    print(f"  {group}: {clusters}")

# Check: do ephys-CCKBC clusters (284, 289) receive mouse CCKBCs?
ephys_cckbc_clusters = set(str(c) for c in CLUSTER_GROUPS["CCKBC"])
mouse_cckbc_clusters = set(
    str(int(float(x))) for x in cluster_stats_df[cluster_stats_df["cckbc_fraction"] > 0.1]["cluster_id"]
)

print(f"\nEphys-defined CCKBC clusters: {ephys_cckbc_clusters}")
print(f"Mouse CCKBC-receiving clusters (>10%): {mouse_cckbc_clusters}")
print(f"Overlap: {ephys_cckbc_clusters & mouse_cckbc_clusters}")

# %% [markdown]
# ## 6. Convergent classification table

# %%
# Build convergent table: cluster_id, mouse_CCKBC_fraction, ephys_group, 22q_bias
convergent_rows = []

# Map cluster_id -> ephys group
cluster_to_ephys = {}
for group, clusters in CLUSTER_GROUPS.items():
    for c in clusters:
        cluster_to_ephys[str(c)] = group

# Load 22q bias results
bias_results_dir = Path(config["bias_data"]["bias_results_dir"])

# Search for 22q bias in preferred order: matched > random
bias_22q_path = None
for pattern in [
    "main_results/matched_*/Centering/22q_del_bias_addP.csv",
    "main_results/random/Centering/22q_del_bias_addP.csv",
    "set_level_matched*/Centering/22q_del_bias_addP.csv",
    "random/Centering/22q_del_bias_addP.csv",
]:
    hits = sorted(bias_results_dir.glob(pattern))
    if hits:
        bias_22q_path = hits[0]
        break

if bias_22q_path:
    bias_22q = pd.read_csv(bias_22q_path, index_col=0)
    print(f"22q bias: {bias_22q.shape} from {bias_22q_path.relative_to(bias_results_dir)}")
else:
    bias_22q = None
    print("22q bias results not found")

# Compute per-cluster human ephys stats (mean of cells mapped to each cluster)
EPHYS_FEATURES = ["avg_rate_hero", "width_ramp", "input_resistance",
                  "isi_cv_hero", "adapt_hero", "sag", "latency_rheo", "rheobase_i",
                  "upstroke_downstroke_ratio_ramp", "tau"]

cluster_ephys_stats = {}
if human_mapping_df is not None and ephys_df is not None:
    # Find cluster column in human mapping
    h_cluster_col = None
    for c in human_mapping_df.columns:
        if "cluster" in c.lower() and "super" not in c.lower():
            h_cluster_col = c
            break

    if h_cluster_col:
        common_ids = ephys_df.index.intersection(human_mapping_df.index)
        h_clusters = human_mapping_df.loc[common_ids, h_cluster_col].astype(float).astype(int)
        for clust_id in h_clusters.unique():
            mask = h_clusters == clust_id
            cell_ids = h_clusters[mask].index
            stats_dict = {}
            for feat in EPHYS_FEATURES:
                if feat in ephys_df.columns:
                    vals = ephys_df.loc[cell_ids, feat].dropna()
                    if len(vals) >= 2:
                        stats_dict[feat] = vals.mean()
            stats_dict["n_ephys_cells"] = int(mask.sum())
            cluster_ephys_stats[clust_id] = stats_dict
        print(f"Computed ephys stats for {len(cluster_ephys_stats)} clusters ({len(common_ids)} cells)")

# Build convergent table
for _, row in cluster_stats_df.iterrows():
    clust = str(int(float(row["cluster_id"])))
    entry = {
        "cluster_id": clust,
        "n_mouse_patchseq": row["n_total_assigned"],
        "n_cckbc": row["n_cckbc"],
        "mouse_cckbc_fraction": row["cckbc_fraction"],
        "ephys_group": cluster_to_ephys.get(clust, "Other"),
    }

    # Add 22q bias if available
    # The bias file index IS the cluster_id (0-461)
    if bias_22q is not None:
        try:
            clust_int = int(float(clust))
            if clust_int in bias_22q.index:
                entry["22q_bias"] = bias_22q.loc[clust_int, "EFFECT"]
                entry["22q_pvalue"] = bias_22q.loc[clust_int, "P-value"]
                entry["22q_zscore"] = bias_22q.loc[clust_int, "Z-score"]
                entry["22q_qvalue"] = bias_22q.loc[clust_int, "q-value"]
                entry["subtype"] = bias_22q.loc[clust_int, "Subtype"]
        except (ValueError, KeyError):
            pass

    # Add per-cluster ephys stats
    clust_int = int(float(clust))
    if clust_int in cluster_ephys_stats:
        for feat in EPHYS_FEATURES:
            if feat in cluster_ephys_stats[clust_int]:
                entry[f"ephys_{feat}"] = cluster_ephys_stats[clust_int][feat]
        entry["n_ephys_cells"] = cluster_ephys_stats[clust_int].get("n_ephys_cells", 0)

    convergent_rows.append(entry)

convergent_df = pd.DataFrame(convergent_rows)
print("\nConvergent classification table:")
# Show compact view with key columns
show_cols = ["cluster_id", "n_cckbc", "mouse_cckbc_fraction", "ephys_group",
             "22q_bias", "n_ephys_cells", "ephys_avg_rate_hero",
             "ephys_input_resistance", "ephys_isi_cv_hero", "ephys_adapt_hero"]
show_cols = [c for c in show_cols if c in convergent_df.columns]
print(convergent_df[show_cols].to_string(index=False))

# Save full table
convergent_df.to_csv(RESULTS_DIR / "convergent_cckbc_classification.csv", index=False)

# %% [markdown]
# ## 7. 22q bias comparison: CCKBC-receiving vs non-receiving CGE clusters

# %%
if "22q_bias" in convergent_df.columns and convergent_df["22q_bias"].notna().any():
    # Split into CCKBC-receiving (>10% fraction) vs non-receiving
    threshold = 0.10
    cckbc_receiving = convergent_df[convergent_df["mouse_cckbc_fraction"] >= threshold]
    non_receiving = convergent_df[convergent_df["mouse_cckbc_fraction"] < threshold]

    bias_receiving = cckbc_receiving["22q_bias"].dropna()
    bias_non_receiving = non_receiving["22q_bias"].dropna()

    print(f"CCKBC-receiving clusters (>={threshold}): n={len(bias_receiving)}")
    print(f"  Mean 22q bias: {bias_receiving.mean():.4f}")
    print(f"Non-receiving clusters: n={len(bias_non_receiving)}")
    print(f"  Mean 22q bias: {bias_non_receiving.mean():.4f}")

    if len(bias_receiving) > 1 and len(bias_non_receiving) > 1:
        stat, pval = stats.mannwhitneyu(bias_receiving, bias_non_receiving, alternative="two-sided")
        print(f"Mann-Whitney U: statistic={stat:.2f}, p={pval:.4f}")

    # Boxplot
    fig, ax = plt.subplots(figsize=(6, 5))
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    data = [bias_receiving.values, bias_non_receiving.values]
    labels = [f"CCKBC-receiving\n(n={len(bias_receiving)})",
              f"Non-receiving\n(n={len(bias_non_receiving)})"]
    bp = ax.boxplot(data, labels=labels, patch_artist=True,
                    boxprops=dict(facecolor="lightcoral", alpha=0.7),
                    medianprops=dict(color="black", linewidth=2))
    bp["boxes"][1].set_facecolor("lightblue")

    ax.set_ylabel("22q mutation bias (EFFECT)")
    ax.set_title("22q bias: CCKBC-receiving vs non-receiving CGE clusters")

    plt.tight_layout()
    plt.savefig(FIG_DIR / "22q_bias_cckbc_comparison.pdf", transparent=True)
    plt.show()

    # Spearman correlation: CCKBC fraction vs 22q bias across clusters
    valid = convergent_df[["mouse_cckbc_fraction", "22q_bias"]].dropna()
    if len(valid) >= 5:
        rho, pval = stats.spearmanr(valid["mouse_cckbc_fraction"], valid["22q_bias"])
        print(f"\nSpearman correlation (CCKBC fraction vs 22q bias): rho={rho:.3f}, p={pval:.4f}")

        fig, ax = plt.subplots(figsize=(6, 5))
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)
        ax.scatter(valid["mouse_cckbc_fraction"], valid["22q_bias"],
                   c="steelblue", s=50, edgecolors="white")

        # Label points with cluster_id
        for idx, row in convergent_df[["mouse_cckbc_fraction", "22q_bias", "cluster_id"]].dropna(subset=["22q_bias"]).iterrows():
            ax.annotate(str(row["cluster_id"]),
                        (row["mouse_cckbc_fraction"], row["22q_bias"]),
                        fontsize=7, alpha=0.7)

        ax.set_xlabel("Mouse CCKBC fraction")
        ax.set_ylabel("22q mutation bias")
        ax.set_title(f"CCKBC fraction vs 22q bias (Spearman rho={rho:.3f}, p={pval:.3f})")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "cckbc_fraction_vs_22q_bias.pdf", transparent=True)
        plt.show()
else:
    print("22q bias data not available for comparison")

# %% [markdown]
# ## 8. Cross-species ephys comparison
#
# Compare mouse CCKBC firing pattern (M1 Patch-seq) with human Patch-seq cells
# in CCKBC-receiving clusters vs ephys-defined CCKBC clusters.

# %%
# Load mouse M1 Patch-seq ephys features
m1_ephys_path = Path(config["mouse_patchseq"]["data_dir"]) / config["mouse_patchseq"]["M1"]["ephys_features"]
m1_meta_path = Path(config["mouse_patchseq"]["data_dir"]) / config["mouse_patchseq"]["M1"]["metadata"]

if m1_ephys_path.exists():
    m1_ephys = pd.read_csv(m1_ephys_path).set_index("cell id")
    m1_meta = pd.read_csv(m1_meta_path, sep="\t").set_index("Cell")
    m1_common = m1_ephys.index.intersection(m1_meta.index)
    m1_ephys = m1_ephys.loc[m1_common]
    m1_meta = m1_meta.loc[m1_common]

    # Identify mouse CCKBC using refined criteria
    sncg_family = m1_meta["RNA family"].astype(str) == "Sncg"
    rna_type = m1_meta["RNA type"].astype(str)
    m1_is_cckbc = sncg_family | rna_type.str.startswith("Vip Sncg") | rna_type.str.startswith("Vip Serpinf1")
    m1_is_vip = (m1_meta["RNA family"].astype(str) == "Vip") & ~m1_is_cckbc
    m1_is_pvalb = m1_meta["RNA family"].astype(str) == "Pvalb"
    m1_is_sst = m1_meta["RNA family"].astype(str) == "Sst"

    print(f"Mouse M1 ephys: {len(m1_ephys)} cells")
    print(f"  CCKBC: {m1_is_cckbc.sum()}, VIP: {m1_is_vip.sum()}, Pvalb: {m1_is_pvalb.sum()}, Sst: {m1_is_sst.sum()}")
else:
    m1_ephys = None
    print(f"M1 ephys not found at {m1_ephys_path}")

# %%
# Mouse CCKBC vs VIP ephys signature
if m1_ephys is not None:
    mouse_features = ["AP width", "AP amplitude", "AP threshold", "AP count",
                      "ISI CV", "ISI adapt", "R_input", "tau", "latency", "AHP"]
    mouse_groups = {"CCKBC": m1_is_cckbc, "VIP": m1_is_vip, "Pvalb": m1_is_pvalb, "Sst": m1_is_sst}

    print("Mouse M1 ephys by cell type:")
    print(f"{'Feature':<16}", end="")
    for g in mouse_groups:
        print(f" {g:>10}", end="")
    print(f" {'CCKBC vs VIP':>12}")
    print("-" * 72)
    for feat in mouse_features:
        if feat in m1_ephys.columns:
            print(f"{feat:<16}", end="")
            for g, mask in mouse_groups.items():
                v = m1_ephys.loc[mask, feat].dropna()
                print(f" {v.mean():>10.3f}" if len(v) > 0 else f" {'N/A':>10}", end="")
            # CCKBC vs VIP test
            c_vals = m1_ephys.loc[m1_is_cckbc, feat].dropna()
            v_vals = m1_ephys.loc[m1_is_vip, feat].dropna()
            if len(c_vals) >= 3 and len(v_vals) >= 3:
                _, pval = stats.mannwhitneyu(c_vals, v_vals, alternative="two-sided")
                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                print(f" p={pval:.4f}{sig}")
            else:
                print()

# %%
# Cross-species: mouse CCKBC profile vs human cluster groups
if m1_ephys is not None and ephys_df is not None and human_mapping_df is not None:
    h_clusters = human_mapping_df.loc[
        ephys_df.index.intersection(human_mapping_df.index), h_cluster_col
    ].astype(float).astype(int)

    human_groups = {
        "Mouse-CCKBC\ntargets(279-281)": h_clusters.isin({279, 280, 281}),
        "Ephys-CCKBC\n(284,289)": h_clusters.isin({284, 289}),
        "Irregular VIP\n(291,292)": h_clusters.isin({291, 292}),
        "LAMP5\n(278,287,288)": h_clusters.isin({278, 287, 288}),
    }

    # Compute mouse CCKBC reference values
    mouse_ref = {}
    for feat in ["AP width", "AP count", "R_input", "ISI CV", "tau", "ISI adapt"]:
        v = m1_ephys.loc[m1_is_cckbc, feat].dropna()
        if len(v) > 0:
            mouse_ref[feat] = v.mean()

    # Cross-species comparison table
    cross_features = [
        ("Firing rate", "AP count", "avg_rate_hero", "exp"),
        ("AP width", "AP width", "width_ramp", "ms_to_s"),
        ("R_input (MOhm)", "R_input", "input_resistance", None),
        ("ISI CV", "ISI CV", "isi_cv_hero", None),
        ("Adaptation", "ISI adapt", "adapt_hero", None),
        ("Tau", "tau", "tau", None),
        ("Sag", None, "sag", None),
        ("Rheobase", None, "rheobase_i", None),
        ("UD ratio", None, "upstroke_downstroke_ratio_ramp", None),
    ]

    print("\nCross-species ephys comparison:")
    print(f"{'Feature':<18} {'Mouse CCKBC':>12}", end="")
    for gname in human_groups:
        print(f" {gname.split(chr(10))[0]:>14}", end="")
    print()
    print("-" * (18 + 12 + 15 * len(human_groups)))

    for label, mouse_feat, human_feat, transform in cross_features:
        print(f"{label:<18}", end="")
        # Mouse value
        if mouse_feat and mouse_feat in mouse_ref:
            val = mouse_ref[mouse_feat]
            if transform == "exp":
                print(f" {np.exp(val):>10.1f}spk", end="")
            elif transform == "ms_to_s":
                print(f" {val:>10.3f}ms", end="")
            else:
                print(f" {val:>12.3f}", end="")
        else:
            print(f" {'N/A':>12}", end="")
        # Human values
        for gname, mask in human_groups.items():
            cell_ids = h_clusters[mask].index
            v = ephys_df.loc[cell_ids.intersection(ephys_df.index), human_feat].dropna()
            if len(v) > 0:
                if transform == "ms_to_s":
                    print(f" {v.mean()*1000:>12.3f}ms", end="")
                elif transform == "exp":
                    print(f" {v.mean():>12.1f}Hz", end="")
                else:
                    print(f" {v.mean():>14.3f}", end="")
            else:
                print(f" {'N/A':>14}", end="")
        print()

    print(f"\n{'N cells':<18} {m1_is_cckbc.sum():>12d}", end="")
    for gname, mask in human_groups.items():
        print(f" {mask.sum():>14d}", end="")
    print()

    # Barplot: key features across groups
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.patch.set_alpha(0)
    plot_features = [
        ("avg_rate_hero", "Firing rate (Hz)"),
        ("input_resistance", "R_input (MOhm)"),
        ("isi_cv_hero", "ISI CV"),
        ("adapt_hero", "Adaptation"),
        ("width_ramp", "AP width (s)"),
        ("sag", "Sag"),
    ]
    group_labels = [g.split("\n")[0] for g in human_groups]
    group_colors = ["firebrick", "steelblue", "goldenrod", "mediumpurple"]

    for ax, (feat, ylabel) in zip(axes.flat, plot_features):
        ax.patch.set_alpha(0)
        means, sems = [], []
        for gname, mask in human_groups.items():
            cell_ids = h_clusters[mask].index
            v = ephys_df.loc[cell_ids.intersection(ephys_df.index), feat].dropna()
            means.append(v.mean() if len(v) > 0 else 0)
            sems.append(v.sem() if len(v) > 1 else 0)
        x = range(len(means))
        ax.bar(x, means, yerr=sems, color=group_colors, edgecolor="white",
               capsize=3, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(group_labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(ylabel, fontsize=10)

    plt.suptitle("Human Patch-seq ephys by cluster group", fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "ephys_by_cluster_group.pdf", transparent=True, bbox_inches="tight")
    plt.show()

# %% [markdown]
# ## 9. Summary
#
# ### Which clusters are confidently CCKBC?
# Clusters receiving both mouse CCKBC mapping AND ephys-based CCKBC assignment
# provide convergent evidence for CCKBC identity.

# %%
# Final summary
print("=" * 60)
print("SUMMARY: CCKBC IDENTIFICATION")
print("=" * 60)

# Top CCKBC-receiving clusters
top_cckbc = cluster_stats_df[cluster_stats_df["cckbc_fraction"] > 0.05].copy()
top_cckbc["ephys_group"] = top_cckbc["cluster_id"].apply(lambda x: str(int(float(x)))).map(cluster_to_ephys).fillna("Other")

print(f"\nTop CCKBC-receiving CGE clusters (>5% mouse CCKBC):")
for _, row in top_cckbc.iterrows():
    print(f"  Cluster {row['cluster_id']}: "
          f"{row['cckbc_fraction']:.1%} CCKBC ({row['n_cckbc']}/{row['n_total_assigned']}), "
          f"ephys group: {row['ephys_group']}")

# Convergence check
convergent_clusters = []
for _, row in top_cckbc.iterrows():
    if row["ephys_group"] == "CCKBC":
        convergent_clusters.append(row["cluster_id"])

if convergent_clusters:
    print(f"\nConvergent CCKBC clusters (both mouse mapping + ephys): {convergent_clusters}")
else:
    print("\nNo convergent clusters found — mouse and ephys CCKBC definitions may point to different clusters")

# Does CCKBC 22q bias differ from VIP/other CGE?
if "22q_bias" in convergent_df.columns:
    print("\n22q bias by ephys group:")
    for group in ["CCKBC", "Regular VIP", "Irregular VIP", "High-rate VIP", "LAMP5"]:
        subset = convergent_df[convergent_df["ephys_group"] == group]
        bias_vals = subset["22q_bias"].dropna()
        if len(bias_vals) > 0:
            print(f"  {group}: mean 22q bias = {bias_vals.mean():.4f} (n={len(bias_vals)})")

# Ephys profile summary
if "ephys_avg_rate_hero" in convergent_df.columns:
    print("\nCross-species ephys summary:")
    print("  Mouse CCKBC signature: moderate firing, high R_input, high ISI irregularity")
    top_3 = convergent_df.nlargest(3, "mouse_cckbc_fraction")
    for _, row in top_3.iterrows():
        rate = row.get("ephys_avg_rate_hero", np.nan)
        rinput = row.get("ephys_input_resistance", np.nan)
        isicv = row.get("ephys_isi_cv_hero", np.nan)
        print(f"  Cluster {row['cluster_id']}: "
              f"CCKBC={row['mouse_cckbc_fraction']:.0%}, "
              f"rate={rate:.1f}Hz, R_in={rinput:.0f}MOhm, ISI_CV={isicv:.3f}")

print(f"\nResults saved to {RESULTS_DIR}")
print(f"Figures saved to {FIG_DIR}")
