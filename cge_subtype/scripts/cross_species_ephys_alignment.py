"""
Cross-species electrophysiology alignment: human vs mouse patch-seq.

Validates human cluster identities by aligning human and mouse cells
in shared electrophysiology feature space. Tests whether human cells
in "CCKBC clusters" actually have CCKBC-like electrophysiology.
"""

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
import umap
import warnings
warnings.filterwarnings("ignore")

OUTDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/results"

# --------------------------------------------------------------------------
# 1. Load data
# --------------------------------------------------------------------------
print("=" * 70)
print("1. LOADING DATA")
print("=" * 70)

# Mouse ephys
mouse_ephys = pd.read_csv(
    "/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_ephys_features.csv"
)

# Mouse metadata
mouse_meta = pd.read_csv(
    "/home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv", sep="\t"
)

# Human ephys
human_ephys = pd.read_csv(
    "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/dat/LeeDalley_ephys_fx.csv",
    index_col=0,
)

# Human cluster mapping
human_map = pd.read_csv(
    "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/cge_subtype/results/human_scvi_mapping_results.csv",
    index_col=0,
)

print(f"Mouse ephys: {mouse_ephys.shape}")
print(f"Mouse metadata: {mouse_meta.shape}")
print(f"Human ephys: {human_ephys.shape}")
print(f"Human cluster mapping: {human_map.shape}")

# --------------------------------------------------------------------------
# 2. Define cell type groups
# --------------------------------------------------------------------------
MOUSE_CCKBC_TYPES = [
    "Vip Sncg", "Vip Serpinf1_1", "Vip Serpinf1_2",
    "Sncg Npy2r", "Sncg Col14a1", "Sncg Calb1_1", "Sncg Calb1_2",
]
MOUSE_VIP_OTHER_TYPES = [
    "Vip Mybpc1_1", "Vip Mybpc1_2", "Vip Mybpc1_3",
    "Vip Gpc3", "Vip Chat_1", "Vip C1ql1", "Vip Htr1f", "Vip Col15a1",
]

HUMAN_VIP_NEG_CCKBC = [277, 278]        # low 22q bias
HUMAN_VIP_POS_CCKBC = [279, 280, 281]   # high 22q bias
HUMAN_VIP_POS_ISI = [276, 284, 285, 289, 290, 291, 292, 293, 294, 295, 296]

# Feature mapping: mouse name -> human name
FEATURE_MAP = {
    "AP width": "width_ramp",
    "AP threshold": "threshold_v_ramp",
    "ISI CV": "isi_cv_hero",
    "R_input": "input_resistance",
    "tau": "tau",
}
SHARED_FEATURES_MOUSE = list(FEATURE_MAP.keys())
SHARED_FEATURES_HUMAN = list(FEATURE_MAP.values())

# --------------------------------------------------------------------------
# 3. Prepare mouse data
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("2. PREPARING MOUSE DATA")
print("=" * 70)

# Merge ephys with metadata
mouse_merged = mouse_ephys.merge(
    mouse_meta[["Cell", "RNA type"]],
    left_on="cell id",
    right_on="Cell",
    how="inner",
)

# Filter to CCKBC and VIP-other types
mouse_cckbc = mouse_merged[mouse_merged["RNA type"].isin(MOUSE_CCKBC_TYPES)].copy()
mouse_vip_other = mouse_merged[mouse_merged["RNA type"].isin(MOUSE_VIP_OTHER_TYPES)].copy()

mouse_cckbc["group"] = "Mouse CCKBC"
mouse_vip_other["group"] = "Mouse VIP-other"

mouse_subset = pd.concat([mouse_cckbc, mouse_vip_other], axis=0)
mouse_feat = mouse_subset[SHARED_FEATURES_MOUSE].copy()
mouse_feat.columns = SHARED_FEATURES_HUMAN  # rename to common names
mouse_feat.index = mouse_subset["cell id"].values

print(f"Mouse CCKBC cells: {len(mouse_cckbc)}")
print(f"Mouse VIP-other cells: {len(mouse_vip_other)}")
print(f"  CCKBC subtypes: {dict(mouse_cckbc['RNA type'].value_counts())}")
print(f"  VIP-other subtypes: {dict(mouse_vip_other['RNA type'].value_counts())}")

# --------------------------------------------------------------------------
# 4. Prepare human data
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("3. PREPARING HUMAN DATA")
print("=" * 70)

# Merge ephys with mapping
human_merged = human_ephys.join(human_map[["assigned_cluster", "assigned_supercluster"]], how="inner")

# Classify into groups
def classify_human(row):
    c = row["assigned_cluster"]
    if c in HUMAN_VIP_NEG_CCKBC:
        return "Human VIP- CCKBC"
    elif c in HUMAN_VIP_POS_CCKBC:
        return "Human VIP+ CCKBC"
    elif c in HUMAN_VIP_POS_ISI:
        return "Human VIP+ ISI"
    else:
        return "Human Other"

human_merged["group"] = human_merged.apply(classify_human, axis=1)

# Filter to groups of interest
human_of_interest = human_merged[human_merged["group"] != "Human Other"].copy()
human_feat = human_of_interest[SHARED_FEATURES_HUMAN].copy()

print(f"Human cells of interest: {len(human_of_interest)}")
print(f"  VIP- CCKBC (clusters {HUMAN_VIP_NEG_CCKBC}): {(human_of_interest['group'] == 'Human VIP- CCKBC').sum()}")
print(f"  VIP+ CCKBC (clusters {HUMAN_VIP_POS_CCKBC}): {(human_of_interest['group'] == 'Human VIP+ CCKBC').sum()}")
print(f"  VIP+ ISI (clusters {HUMAN_VIP_POS_ISI}): {(human_of_interest['group'] == 'Human VIP+ ISI').sum()}")

# --------------------------------------------------------------------------
# 5. Combine, handle missing values, z-score within species
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("4. PREPROCESSING: Z-SCORING WITHIN SPECIES")
print("=" * 70)

# Drop rows with any NaN in shared features
mouse_feat_clean = mouse_feat.dropna()
human_feat_clean = human_feat.dropna()

# Map groups
mouse_groups = mouse_subset.set_index("cell id").loc[mouse_feat_clean.index, "group"]
human_groups = human_of_interest.loc[human_feat_clean.index, "group"]

print(f"After dropping NaN:")
print(f"  Mouse: {len(mouse_feat_clean)} (CCKBC={sum(mouse_groups=='Mouse CCKBC')}, VIP-other={sum(mouse_groups=='Mouse VIP-other')})")
print(f"  Human: {len(human_feat_clean)}")
for g in ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]:
    print(f"    {g}: {sum(human_groups == g)}")

# Z-score within species (so cross-species comparison is in relative space)
scaler_mouse = StandardScaler()
scaler_human = StandardScaler()

mouse_z = pd.DataFrame(
    scaler_mouse.fit_transform(mouse_feat_clean),
    index=mouse_feat_clean.index,
    columns=mouse_feat_clean.columns,
)
human_z = pd.DataFrame(
    scaler_human.fit_transform(human_feat_clean),
    index=human_feat_clean.index,
    columns=human_feat_clean.columns,
)

print("\nMouse feature stats (pre-z-score):")
print(mouse_feat_clean.describe().round(2))
print("\nHuman feature stats (pre-z-score):")
print(human_feat_clean.describe().round(2))

# Combine
combined_z = pd.concat([mouse_z, human_z], axis=0)
combined_groups = pd.concat([mouse_groups, human_groups], axis=0)
combined_species = pd.Series(
    ["Mouse"] * len(mouse_z) + ["Human"] * len(human_z),
    index=combined_z.index,
)

print(f"\nCombined matrix: {combined_z.shape}")

# --------------------------------------------------------------------------
# 6. PCA on combined z-scored data
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("5. PCA ON COMBINED Z-SCORED DATA")
print("=" * 70)

pca = PCA(n_components=min(5, combined_z.shape[1]))
pca_coords = pca.fit_transform(combined_z)

print(f"Explained variance ratio: {pca.explained_variance_ratio_.round(4)}")
print(f"Cumulative: {np.cumsum(pca.explained_variance_ratio_).round(4)}")
print(f"\nPC loadings:")
loadings = pd.DataFrame(
    pca.components_.T,
    index=SHARED_FEATURES_HUMAN,
    columns=[f"PC{i+1}" for i in range(pca.n_components_)],
)
print(loadings.round(3))

# --------------------------------------------------------------------------
# 7. UMAP on combined z-scored data
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("6. UMAP EMBEDDING")
print("=" * 70)

reducer = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42, n_components=2)
umap_coords = reducer.fit_transform(combined_z)
print(f"UMAP embedding computed: {umap_coords.shape}")

# --------------------------------------------------------------------------
# 8. Visualize PCA and UMAP
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("7. GENERATING FIGURES")
print("=" * 70)

# Color and marker scheme
GROUP_COLORS = {
    "Mouse CCKBC":       "#1f77b4",  # blue
    "Mouse VIP-other":   "#ff7f0e",  # orange
    "Human VIP- CCKBC":  "#2ca02c",  # green
    "Human VIP+ CCKBC":  "#d62728",  # red
    "Human VIP+ ISI":    "#9467bd",  # purple
}
GROUP_MARKERS = {
    "Mouse CCKBC":       "o",
    "Mouse VIP-other":   "s",
    "Human VIP- CCKBC":  "^",
    "Human VIP+ CCKBC":  "D",
    "Human VIP+ ISI":    "v",
}
# Plot order: mouse first (background), then human on top
PLOT_ORDER = ["Mouse VIP-other", "Mouse CCKBC", "Human VIP+ ISI", "Human VIP- CCKBC", "Human VIP+ CCKBC"]


def scatter_with_groups(ax, coords, groups, title, xlabel, ylabel):
    """Scatter plot colored/shaped by group."""
    for grp in PLOT_ORDER:
        mask = groups.values == grp
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=GROUP_COLORS[grp],
            marker=GROUP_MARKERS[grp],
            s=40 if "Mouse" in grp else 60,
            alpha=0.6 if "Mouse" in grp else 0.85,
            edgecolors="white",
            linewidths=0.3,
            label=f"{grp} (n={mask.sum()})",
            zorder=2 if "Mouse" in grp else 3,
        )
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.patch.set_alpha(0)


# --- Figure 1: PCA + UMAP side by side ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
fig.patch.set_alpha(0)

scatter_with_groups(
    axes[0], pca_coords[:, :2], combined_groups,
    "PCA: Cross-species Ephys Alignment",
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)",
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
)
scatter_with_groups(
    axes[1], umap_coords, combined_groups,
    "UMAP: Cross-species Ephys Alignment",
    "UMAP1", "UMAP2",
)

# Shared legend
handles = [
    Line2D([0], [0], marker=GROUP_MARKERS[g], color="w", markerfacecolor=GROUP_COLORS[g],
           markersize=9, label=g, markeredgecolor="grey", linewidth=0)
    for g in PLOT_ORDER
]
fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=10,
           frameon=True, fancybox=True, framealpha=0.8)
fig.subplots_adjust(bottom=0.18)
plt.tight_layout(rect=[0, 0.12, 1, 1])

figpath1 = f"{OUTDIR}/fig_cross_species_ephys_pca_umap.pdf"
fig.savefig(figpath1, dpi=150, transparent=True, bbox_inches="tight")
print(f"Saved: {figpath1}")

# --- Figure 2: PCA loadings biplot ---
fig2, ax2 = plt.subplots(figsize=(8, 7))
fig2.patch.set_alpha(0)
ax2.patch.set_alpha(0)

scatter_with_groups(
    ax2, pca_coords[:, :2], combined_groups,
    "PCA Biplot with Feature Loadings",
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)",
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
)

# Add loading arrows
scale = max(np.abs(pca_coords[:, :2]).max(axis=0)) * 0.8
for i, feat in enumerate(SHARED_FEATURES_HUMAN):
    ax2.annotate(
        feat,
        xy=(pca.components_[0, i] * scale, pca.components_[1, i] * scale),
        xytext=(pca.components_[0, i] * scale * 1.15, pca.components_[1, i] * scale * 1.15),
        fontsize=9, fontweight="bold", color="black",
        arrowprops=dict(arrowstyle="<-", color="black", lw=1.5),
        ha="center", va="center",
    )

ax2.legend(loc="upper right", fontsize=8, framealpha=0.8)
figpath2 = f"{OUTDIR}/fig_cross_species_pca_biplot.pdf"
fig2.savefig(figpath2, dpi=150, transparent=True, bbox_inches="tight")
print(f"Saved: {figpath2}")

# --------------------------------------------------------------------------
# 9. Nearest-Neighbor Analysis
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("8. NEAREST-NEIGHBOR ANALYSIS (k=5)")
print("=" * 70)

K = 5

# Fit NN on mouse data
nn = NearestNeighbors(n_neighbors=K, metric="euclidean")
nn.fit(mouse_z.values)

# For each human cell, find k nearest mouse neighbors
distances, indices = nn.kneighbors(human_z.values)

# Map mouse indices to cell types
mouse_labels = mouse_groups.values  # "Mouse CCKBC" or "Mouse VIP-other"
human_cell_ids = human_z.index
human_group_labels = human_groups.values

# For each human cell, compute fraction of CCKBC neighbors
cckbc_fracs = []
for i in range(len(human_z)):
    neighbor_labels = mouse_labels[indices[i]]
    frac_cckbc = np.mean(neighbor_labels == "Mouse CCKBC")
    cckbc_fracs.append(frac_cckbc)

nn_results = pd.DataFrame({
    "cell_id": human_cell_ids,
    "group": human_group_labels,
    "frac_cckbc_neighbors": cckbc_fracs,
    "mean_neighbor_dist": distances.mean(axis=1),
})

print("\nFraction of nearest mouse neighbors that are CCKBC, by human group:")
print("-" * 65)
for grp in ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]:
    sub = nn_results[nn_results["group"] == grp]
    if len(sub) == 0:
        continue
    print(f"  {grp} (n={len(sub)}):")
    print(f"    Mean CCKBC fraction:   {sub['frac_cckbc_neighbors'].mean():.3f}")
    print(f"    Median CCKBC fraction: {sub['frac_cckbc_neighbors'].median():.3f}")
    print(f"    Std:                   {sub['frac_cckbc_neighbors'].std():.3f}")
    print(f"    Cells with >50% CCKBC neighbors: {(sub['frac_cckbc_neighbors'] > 0.5).sum()}/{len(sub)} ({(sub['frac_cckbc_neighbors'] > 0.5).mean()*100:.1f}%)")
    print(f"    Mean neighbor distance: {sub['mean_neighbor_dist'].mean():.3f}")

# --- Figure 3: NN results boxplot ---
fig3, ax3 = plt.subplots(figsize=(8, 5))
fig3.patch.set_alpha(0)
ax3.patch.set_alpha(0)

human_groups_ordered = ["Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]
group_data = [nn_results.loc[nn_results["group"] == g, "frac_cckbc_neighbors"].values for g in human_groups_ordered]
group_colors = [GROUP_COLORS[g] for g in human_groups_ordered]

bp = ax3.boxplot(group_data, patch_artist=True, widths=0.5,
                 medianprops=dict(color="black", linewidth=2))

for patch, color in zip(bp["boxes"], group_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

# Overlay individual points (jittered)
rng = np.random.default_rng(42)
for i, (data, color) in enumerate(zip(group_data, group_colors)):
    jitter = rng.uniform(-0.12, 0.12, size=len(data))
    ax3.scatter(
        np.full(len(data), i + 1) + jitter,
        data, c=color, alpha=0.5, s=20, edgecolors="white", linewidths=0.3, zorder=3,
    )

ax3.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(human_groups_ordered, group_data)], fontsize=10)
ax3.set_ylabel("Fraction CCKBC among k=5 nearest mouse neighbors", fontsize=11)
ax3.set_title("Nearest-Neighbor Validation of Human Cluster Identity", fontsize=13, fontweight="bold")
ax3.axhline(0.5, ls="--", color="grey", alpha=0.5, label="Chance if equal CCKBC/VIP-other")
# Add mouse CCKBC baseline fraction
baseline = len(mouse_cckbc) / (len(mouse_cckbc) + len(mouse_vip_other))
ax3.axhline(baseline, ls=":", color="darkblue", alpha=0.5, label=f"Mouse CCKBC base rate ({baseline:.2f})")
ax3.legend(fontsize=9)
ax3.set_ylim(-0.05, 1.05)

figpath3 = f"{OUTDIR}/fig_cross_species_nn_validation.pdf"
fig3.savefig(figpath3, dpi=150, transparent=True, bbox_inches="tight")
print(f"\nSaved: {figpath3}")

# --------------------------------------------------------------------------
# 10. Classifier Analysis
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("9. CLASSIFIER ANALYSIS")
print("=" * 70)

# Train on mouse, predict on human
mouse_labels_binary = (mouse_groups == "Mouse CCKBC").astype(int).values
mouse_X = mouse_z.values
human_X = human_z.values

# Cross-validation on mouse data first
print("\nMouse-only cross-validation (5-fold):")
for name, clf in [
    ("Logistic Regression", LogisticRegression(max_iter=1000, random_state=42)),
    ("k-NN (k=5)", KNeighborsClassifier(n_neighbors=5)),
]:
    scores = cross_val_score(clf, mouse_X, mouse_labels_binary, cv=5, scoring="accuracy")
    print(f"  {name}: accuracy = {scores.mean():.3f} +/- {scores.std():.3f}")

# Train on ALL mouse data, predict on human
print("\nPredictions on human cells (trained on all mouse cells):")
print("-" * 65)

for name, clf in [
    ("Logistic Regression", LogisticRegression(max_iter=1000, random_state=42)),
    ("k-NN (k=5)", KNeighborsClassifier(n_neighbors=5)),
]:
    clf.fit(mouse_X, mouse_labels_binary)
    human_pred = clf.predict(human_X)
    human_prob = clf.predict_proba(human_X)[:, 1]  # P(CCKBC)

    print(f"\n  {name}:")
    for grp in human_groups_ordered:
        mask = human_groups.values == grp
        if mask.sum() == 0:
            continue
        pred_cckbc = human_pred[mask].sum()
        total = mask.sum()
        mean_prob = human_prob[mask].mean()
        print(f"    {grp} (n={total}):")
        print(f"      Predicted CCKBC: {pred_cckbc}/{total} ({pred_cckbc/total*100:.1f}%)")
        print(f"      Mean P(CCKBC):   {mean_prob:.3f}")

# --- Figure 4: Classifier probability distribution ---
fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))
fig4.patch.set_alpha(0)

for ax_idx, (name, clf) in enumerate([
    ("Logistic Regression", LogisticRegression(max_iter=1000, random_state=42)),
    ("k-NN (k=5)", KNeighborsClassifier(n_neighbors=5)),
]):
    clf.fit(mouse_X, mouse_labels_binary)
    human_prob = clf.predict_proba(human_X)[:, 1]

    ax = axes4[ax_idx]
    ax.patch.set_alpha(0)

    for grp in human_groups_ordered:
        mask = human_groups.values == grp
        if mask.sum() == 0:
            continue
        ax.hist(
            human_prob[mask], bins=np.linspace(0, 1, 21),
            alpha=0.5, color=GROUP_COLORS[grp], label=f"{grp} (n={mask.sum()})",
            edgecolor="white", linewidth=0.5,
        )

    ax.axvline(0.5, ls="--", color="grey", alpha=0.7)
    ax.set_xlabel("P(CCKBC)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(f"{name}", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

fig4.suptitle("Classifier Predicted CCKBC Probability for Human Cells",
              fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()

figpath4 = f"{OUTDIR}/fig_cross_species_classifier_probs.pdf"
fig4.savefig(figpath4, dpi=150, transparent=True, bbox_inches="tight")
print(f"\nSaved: {figpath4}")

# --------------------------------------------------------------------------
# 11. Feature-level comparison
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("10. FEATURE-LEVEL COMPARISON (Z-SCORED WITHIN SPECIES)")
print("=" * 70)

# Compare raw feature distributions across groups
all_groups = pd.concat([mouse_groups, human_groups])
all_z = pd.concat([mouse_z, human_z])

print(f"\n{'Feature':<22} {'Mouse CCKBC':>14} {'Mouse VIP-oth':>14} {'H VIP- CCKBC':>14} {'H VIP+ CCKBC':>14} {'H VIP+ ISI':>14}")
print("-" * 95)

for feat in SHARED_FEATURES_HUMAN:
    vals = []
    for grp in ["Mouse CCKBC", "Mouse VIP-other", "Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]:
        mask = all_groups == grp
        vals.append(f"{all_z.loc[mask, feat].mean():.2f}")
    print(f"{feat:<22} {'  '.join(f'{v:>12}' for v in vals)}")

# --- Figure 5: Feature comparison violin plots ---
fig5, axes5 = plt.subplots(1, len(SHARED_FEATURES_HUMAN), figsize=(4 * len(SHARED_FEATURES_HUMAN), 5))
fig5.patch.set_alpha(0)

all_plot_groups = ["Mouse CCKBC", "Mouse VIP-other", "Human VIP- CCKBC", "Human VIP+ CCKBC", "Human VIP+ ISI"]

for i, feat in enumerate(SHARED_FEATURES_HUMAN):
    ax = axes5[i]
    ax.patch.set_alpha(0)

    data_by_group = []
    for grp in all_plot_groups:
        mask = all_groups == grp
        data_by_group.append(all_z.loc[mask, feat].dropna().values)

    parts = ax.violinplot(data_by_group, positions=range(len(all_plot_groups)),
                          showmeans=True, showextrema=False)
    for j, (pc, grp) in enumerate(zip(parts["bodies"], all_plot_groups)):
        pc.set_facecolor(GROUP_COLORS[grp])
        pc.set_alpha(0.6)
    parts["cmeans"].set_color("black")

    ax.set_xticks(range(len(all_plot_groups)))
    ax.set_xticklabels(["M\nCCKBC", "M\nVIP-o", "H\nVIP-\nCCKBC", "H\nVIP+\nCCKBC", "H\nVIP+\nISI"],
                        fontsize=8)
    ax.set_title(feat, fontsize=10, fontweight="bold")
    if i == 0:
        ax.set_ylabel("Z-score (within species)", fontsize=10)

fig5.suptitle("Feature Distributions by Group (Z-scored within species)",
              fontsize=13, fontweight="bold")
plt.tight_layout()

figpath5 = f"{OUTDIR}/fig_cross_species_feature_violins.pdf"
fig5.savefig(figpath5, dpi=150, transparent=True, bbox_inches="tight")
print(f"\nSaved: {figpath5}")

# --------------------------------------------------------------------------
# 12. Statistical tests
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("11. STATISTICAL TESTS")
print("=" * 70)

from scipy import stats

print("\nMann-Whitney U tests on CCKBC neighbor fraction:")
print("-" * 65)

# VIP- CCKBC vs VIP+ ISI
vip_neg = nn_results.loc[nn_results["group"] == "Human VIP- CCKBC", "frac_cckbc_neighbors"]
vip_pos_cckbc = nn_results.loc[nn_results["group"] == "Human VIP+ CCKBC", "frac_cckbc_neighbors"]
vip_pos_isi = nn_results.loc[nn_results["group"] == "Human VIP+ ISI", "frac_cckbc_neighbors"]

comparisons = [
    ("VIP- CCKBC vs VIP+ ISI", vip_neg, vip_pos_isi),
    ("VIP+ CCKBC vs VIP+ ISI", vip_pos_cckbc, vip_pos_isi),
    ("VIP- CCKBC vs VIP+ CCKBC", vip_neg, vip_pos_cckbc),
    ("All CCKBC (VIP-/VIP+) vs VIP+ ISI",
     pd.concat([vip_neg, vip_pos_cckbc]), vip_pos_isi),
]

for label, g1, g2 in comparisons:
    if len(g1) == 0 or len(g2) == 0:
        print(f"  {label}: skipped (empty group)")
        continue
    stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
    effect_size = g1.mean() - g2.mean()
    print(f"  {label}:")
    print(f"    U={stat:.1f}, p={pval:.4e}, delta_mean={effect_size:+.3f}")

# Kruskal-Wallis across all 3 groups
if len(vip_neg) > 0 and len(vip_pos_cckbc) > 0 and len(vip_pos_isi) > 0:
    stat_kw, pval_kw = stats.kruskal(vip_neg, vip_pos_cckbc, vip_pos_isi)
    print(f"\n  Kruskal-Wallis (3-group): H={stat_kw:.2f}, p={pval_kw:.4e}")

# --------------------------------------------------------------------------
# 13. Summary
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

# Final NN summary
print("\nNearest-neighbor CCKBC fraction by human group:")
for grp in human_groups_ordered:
    sub = nn_results[nn_results["group"] == grp]
    if len(sub) > 0:
        print(f"  {grp}: {sub['frac_cckbc_neighbors'].mean():.3f} (n={len(sub)})")

# Classifier summary (LR)
lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(mouse_X, mouse_labels_binary)
human_pred_lr = lr.predict(human_X)
human_prob_lr = lr.predict_proba(human_X)[:, 1]

print("\nLogistic Regression CCKBC prediction rate by human group:")
for grp in human_groups_ordered:
    mask = human_groups.values == grp
    if mask.sum() > 0:
        rate = human_pred_lr[mask].mean()
        prob = human_prob_lr[mask].mean()
        print(f"  {grp}: {rate*100:.1f}% classified CCKBC, mean P(CCKBC)={prob:.3f} (n={mask.sum()})")

print(f"\nAll figures saved to: {OUTDIR}/")
print("Done.")
