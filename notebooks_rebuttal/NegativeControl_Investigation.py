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
#   kernelspec:
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Negative Control Investigation
#
# Systematically evaluate negative control traits to select the most appropriate
# set for the main paper figure. Tests all available genebass pLoF burden traits
# and the SCZ protective gene set across different weighting methods.
#
# **Gene sets (genebass pLoF burden, UKBB exome sequencing):**
#
# | Trait | Category | Expected pattern |
# |-------|----------|-----------------|
# | HDL Cholesterol | Non-brain metabolic | No neuronal enrichment |
# | Alanine AT | Liver enzyme | No neuronal enrichment |
# | HbA1c | Blood sugar | No neuronal enrichment |
# | T2D | Metabolic disease | No neuronal enrichment |
# | IBD | Gut inflammation | No neuronal enrichment |
# | Red blood cell count | Hematological | No neuronal enrichment |
# | Parkinson's | Neurological | Possible neuronal signal |
# | Alzheimer's | Neurological | Possible neuronal signal |
# | SCZ protective | Psychiatric (reversed) | Inverted CGE pattern |
#
# **Weighting methods:**
# 1. **Uniform (weight=1)** — all top-61 genes equal
# 2. **abs(BETA Burden)** — weighted by effect size
# 3. **Positive BETA only** — only genes with same-direction effect
#
# **Goal:** Identify which traits + weighting give a clean null (no neuronal enrichment)
# to use as negative controls in the main paper.

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import OrderedDict
import yaml

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'

FIG_DIR = PROJ_DIR / "results" / "figures" / "negative_controls"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Expression matrix
expression_matrix = _cfg['analysis_types']['Centering']
HumanCT_Z2 = pd.read_csv(str(PROJ_DIR / expression_matrix), index_col=0)
HumanCT_Z2.columns = HumanCT_Z2.columns.astype(int)

# %% [markdown]
# ## 1. Load all genebass trait data

# %%
CTRL_DIR = PROJ_DIR / "dat" / "CTRL"

TRAIT_FILES = OrderedDict([
    # Non-brain traits (expected: no neuronal enrichment)
    ("HDL Cholesterol",   CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30760-both_sexes--irnt_2025_11_25_15_49_42.csv"),
    ("Alanine AT",        CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30620-both_sexes--irnt_2025_12_16_20_18_07.csv"),
    ("HbA1c",             CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30750-both_sexes--irnt_2025_11_25_16_32_28.csv"),
    ("T2D",               CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-T2D_custom-both_sexes--custom_2025_11_25_16_14_20.csv"),
    ("IBD",               CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-IBD_custom2-both_sexes--custom_2025_11_25_16_13_49.csv"),
    ("Red Blood Cell",    CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30070-both_sexes--irnt_2025_12_16_20_41_11.csv"),
    # Neurological traits (may show some neuronal signal)
    ("Parkinson's",       CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131022-both_sexes--_2025_11_25_16_56_06.csv"),
    ("Alzheimer's",       CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131036-both_sexes--_2025_11_25_17_00_04.csv"),
])

# Annotate all traits with EntrezID
trait_data = {}
for trait_name, filepath in TRAIT_FILES.items():
    df = pd.read_csv(filepath, index_col=0)
    df["EntrezID"] = df.index.map(
        lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez.get(x)) else None
    )
    df = df[df["EntrezID"].notna()].copy()
    trait_data[trait_name] = df
    n_pos = (df.head(61)["BETA Burden"].dropna() > 0).sum()
    n_neg = (df.head(61)["BETA Burden"].dropna() < 0).sum()
    print(f"  {trait_name}: {len(df)} genes, top 61 BETA: {n_pos} pos / {n_neg} neg")

# %% [markdown]
# ## 2. Compute bias for all traits × weighting methods

# %%
def make_gene_weights(df, n_top=61, method="uniform"):
    """
    Create gene weight dict from annotated trait DataFrame.

    Methods:
    - 'uniform': weight=1 for all genes
    - 'abs_beta': weight=abs(BETA Burden)
    - 'pos_beta': only positive BETA genes, weight=BETA
    - 'neg_beta': only negative BETA genes, weight=abs(BETA)
    """
    top = df.head(n_top).copy()

    if method == "uniform":
        return dict(zip(top["EntrezID"].astype(int), [1] * len(top)))
    elif method == "abs_beta":
        valid = top[top["BETA Burden"].notna()]
        return dict(zip(valid["EntrezID"].astype(int), valid["BETA Burden"].abs()))
    elif method == "pos_beta":
        valid = top[(top["BETA Burden"].notna()) & (top["BETA Burden"] > 0)]
        return dict(zip(valid["EntrezID"].astype(int), valid["BETA Burden"]))
    elif method == "neg_beta":
        valid = top[(top["BETA Burden"].notna()) & (top["BETA Burden"] < 0)]
        return dict(zip(valid["EntrezID"].astype(int), valid["BETA Burden"].abs()))
    else:
        raise ValueError(f"Unknown method: {method}")


METHODS = ["uniform", "abs_beta"]
results = {}  # {(trait, method): bias_df}

for trait_name, df in trait_data.items():
    for method in METHODS:
        gw = make_gene_weights(df, n_top=61, method=method)
        if len(gw) == 0:
            print(f"  SKIP: {trait_name} / {method} — no genes")
            continue
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2, gw)
        bias = AnnotateCTDat(bias, Anno)
        results[(trait_name, method)] = bias

# Also add SCZ protective
SCZ_protect_gw = Fil2Dict(str(PROJ_DIR / "dat/GeneWeights/SCZ.top61.protect.gw"))
SCZ_protect_bias = AnnotateCTDat(HumanCT_AvgZ_Weighted(HumanCT_Z2, SCZ_protect_gw), Anno)
results[("SCZ protective", "case-control excess")] = SCZ_protect_bias

print(f"\nComputed {len(results)} bias profiles")

# %% [markdown]
# ## 3. Summary: Top supercluster and CGE interneuron rank

# %%
summary_rows = []
for (trait, method), bias_df in results.items():
    sc_mean = bias_df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
    top_sc = sc_mean.index[0]
    top_sc_val = sc_mean.iloc[0]

    # CGE interneuron rank and value
    cge_val = sc_mean.get("CGE interneuron", np.nan)
    cge_rank = list(sc_mean.index).index("CGE interneuron") + 1 if "CGE interneuron" in sc_mean.index else None

    # Is any neuronal supercluster in top 3?
    top3 = sc_mean.head(3).index.tolist()
    neur_superclusters = set(bias_df[bias_df["Class"].isin(["NEUR", "Neuron"])]["Supercluster"].unique())
    top3_neuronal = [s for s in top3 if s in neur_superclusters]

    summary_rows.append({
        "Trait": trait,
        "Method": method,
        "Top Supercluster": top_sc,
        "Top SC Value": top_sc_val,
        "CGE Rank": cge_rank,
        "CGE Value": cge_val,
        "Neuronal in Top 3": ", ".join(top3_neuronal) if top3_neuronal else "None",
    })

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))

# %% [markdown]
# ## 4. Supercluster bias boxplots — all traits (uniform weights)

# %%
# Plot all traits with uniform weights
uniform_traits = [k for k in results if k[1] == "uniform"]
n_traits = len(uniform_traits)
n_cols = 3
n_rows = (n_traits + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows), dpi=100)
axes = axes.flatten()

for i, (trait, method) in enumerate(uniform_traits):
    ax = axes[i]
    SuperClusterBias_BoxPlot(results[(trait, method)], trait, ax=ax)

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Negative Control Traits — Supercluster Bias (uniform weights, top 61 genes)",
             fontsize=16, fontweight="bold", y=1.01)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(str(FIG_DIR / "all_traits_uniform_boxplots.png"),
            dpi=200, bbox_inches='tight', transparent=True)
plt.show()

# %% [markdown]
# ## 5. Compare uniform vs abs(BETA) weighting

# %%
fig, axes = plt.subplots(len(trait_data), 2, figsize=(16, 5 * len(trait_data)), dpi=80)

for i, trait in enumerate(trait_data.keys()):
    for j, method in enumerate(["uniform", "abs_beta"]):
        ax = axes[i, j]
        key = (trait, method)
        if key in results:
            SuperClusterBias_BoxPlot(results[key], f"{trait}\n({method})", ax=ax)
        else:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center")
            ax.set_title(f"{trait} ({method})")

fig.suptitle("Uniform vs abs(BETA) weighting across all traits",
             fontsize=16, fontweight="bold", y=1.005)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(str(FIG_DIR / "uniform_vs_beta_comparison.png"),
            dpi=150, bbox_inches='tight', transparent=True)
plt.show()

# %% [markdown]
# ## 6. CGE interneuron bias: all traits compared to psychiatric disorders

# %%
# Load psychiatric disorder biases for comparison
Bias_Dir = str(PROJ_DIR / "results/main_results/random/Centering/") + "/"
psych_biases = {
    "ASD (all)": pd.read_csv(Bias_Dir + "ASD_All_bias_addP.csv", index_col=0),
    "ASD w/o ID": pd.read_csv(Bias_Dir + "ASD_HIQ_bias_addP.csv", index_col=0),
    "ASD with ID": pd.read_csv(Bias_Dir + "ASD_LIQ_bias_addP.csv", index_col=0),
    "SCZ": pd.read_csv(Bias_Dir + "SCZ_bias_addP.csv", index_col=0),
    "DDD": pd.read_csv(Bias_Dir + "DDD_61_bias_addP.csv", index_col=0),
    "22q11.2 del": pd.read_csv(Bias_Dir + "22q_del_bias_addP.csv", index_col=0),
}

# Compute CGE interneuron mean bias for all
cge_values = []
for label, bias_df in psych_biases.items():
    cge_mean = bias_df[bias_df["Supercluster"] == "CGE interneuron"]["EFFECT"].mean()
    cge_values.append({"Trait": label, "Category": "Psychiatric", "CGE Mean Bias": cge_mean})

for (trait, method), bias_df in results.items():
    if method != "uniform":
        continue
    cge_mean = bias_df[bias_df["Supercluster"] == "CGE interneuron"]["EFFECT"].mean()
    cat = "Neurological" if trait in ["Parkinson's", "Alzheimer's"] else "Non-brain"
    if trait == "SCZ protective":
        cat = "Protective"
    cge_values.append({"Trait": trait, "Category": cat, "CGE Mean Bias": cge_mean})

cge_df = pd.DataFrame(cge_values).sort_values("CGE Mean Bias", ascending=False)

# Plot
fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
colors = {"Psychiatric": "#E74C3C", "Non-brain": "#95A5A6", "Neurological": "#3498DB", "Protective": "#2ECC71"}
bar_colors = [colors[cat] for cat in cge_df["Category"]]

bars = ax.barh(range(len(cge_df)), cge_df["CGE Mean Bias"].values, color=bar_colors, edgecolor="white", lw=0.5)
ax.set_yticks(range(len(cge_df)))
ax.set_yticklabels(cge_df["Trait"].values, fontsize=11)
ax.set_xlabel("Mean CGE Interneuron Bias (EFFECT)", fontsize=12)
ax.set_title("CGE Interneuron Enrichment: Psychiatric vs Control Traits", fontsize=14, fontweight="bold")
ax.axvline(0, color="black", lw=0.8)
ax.invert_yaxis()

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=l) for l, c in colors.items()]
ax.legend(handles=legend_elements, loc="lower right", fontsize=10)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(str(FIG_DIR / "CGE_bias_psychiatric_vs_controls.png"),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

print(cge_df.to_string(index=False))

# %% [markdown]
# ## 7. Pipeline P-values: -log10P boxplots (matched null)
#
# Load bias+P results from the matched-null pipeline and plot -log10P by supercluster.
# This shows which cell types are significantly enriched after controlling for
# gene length, expression, and conservation.

# %%
# Load all negative control + psychiatric p-value results (matched null)
MATCHED_DIR = PROJ_DIR / "results" / "main_results" / "matched_WB_mean_phastCons_n_CDS_bases_Best1000" / "Centering"
RANDOM_DIR = PROJ_DIR / "results" / "main_results" / "random" / "Centering"

pval_traits = OrderedDict([
    # Psychiatric (for comparison)
    ("ASD w/o ID", "ASD_HIQ"),
    ("SCZ", "SCZ"),
    ("DDD", "DDD_61"),
    # Protective
    ("SCZ protective", "SCZ_protect"),
    # Non-brain
    ("HDL Cholesterol", "NegCtrl_HDL"),
    ("Alanine AT", "NegCtrl_Alanine"),
    ("HbA1c", "NegCtrl_HbA1c"),
    ("T2D", "NegCtrl_T2D"),
    ("IBD", "NegCtrl_IBD"),
    ("Red Blood Cell", "NegCtrl_RBC"),
    # Neurological
    ("Parkinson's", "NegCtrl_Parkinson"),
    ("Alzheimer's", "NegCtrl_Alzheimer"),
])

pval_results = {}
for label, key in pval_traits.items():
    fpath = MATCHED_DIR / f"{key}_bias_addP.csv"
    if fpath.exists():
        pval_results[label] = pd.read_csv(fpath, index_col=0)
    else:
        print(f"MISSING: {label} ({fpath.name})")

print(f"Loaded {len(pval_results)} bias+P results")

# %%
# Check for significant neuronal cell types (FDR < 0.1) across all traits
print("=" * 80)
print("SIGNIFICANT NEURONAL CELL TYPES (q < 0.1, matched null)")
print("=" * 80)
for label, df in pval_results.items():
    sig_neur = df[(df["q-value"] < 0.1) & (df["Class"].isin(["NEUR", "Neuron"]))]
    if len(sig_neur) > 0:
        print(f"\n  {label}: {len(sig_neur)} significant neuronal cell types")
        for _, row in sig_neur.iterrows():
            print(f"    {row['Supercluster']:<40s} Z={row['Z-score']:.2f}  q={row['q-value']:.4f}")
    else:
        print(f"  {label}: 0 significant neuronal cell types")

# %%
# Plot -log10P boxplots for all traits (matched null)
n_traits = len(pval_results)
n_cols = 3
n_rows = (n_traits + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows), dpi=100)
axes = axes.flatten()

for i, (label, df) in enumerate(pval_results.items()):
    ax = axes[i]
    SuperClusterBias_BoxPlot(df, label, EffectCol="-logP", ax=ax)
    # Add FDR threshold lines
    if "q-value" in df.columns:
        for fdr_cut, color, ls in [(0.05, "red", "--"), (0.1, "orange", ":")]:
            sig_mask = df["q-value"] < fdr_cut
            if sig_mask.any():
                pval_threshold = df.loc[sig_mask, "P-value"].max()
                ax.axhline(-np.log10(pval_threshold), color=color, ls=ls, lw=1.5, alpha=0.8,
                           label=f"FDR={fdr_cut}")
        ax.legend(fontsize=7, loc="upper right", framealpha=0.8)

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Cell-type mutation bias significance (-log10P, matched null)",
             fontsize=16, fontweight="bold", y=1.01)
fig.tight_layout()
fig.patch.set_alpha(0)
fig.savefig(str(FIG_DIR / "all_traits_matched_logP_boxplots.png"),
            dpi=200, bbox_inches='tight', transparent=True)
plt.show()

# %%
# Same for random null
fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows), dpi=100)
axes2 = axes2.flatten()

for i, (label, key) in enumerate(pval_traits.items()):
    ax = axes2[i]
    fpath = RANDOM_DIR / f"{key}_bias_addP.csv"
    if fpath.exists():
        df_rand = pd.read_csv(fpath, index_col=0)
        SuperClusterBias_BoxPlot(df_rand, f"{label} (random null)", EffectCol="-logP", ax=ax)
        if "q-value" in df_rand.columns:
            for fdr_cut, color, ls in [(0.05, "red", "--"), (0.1, "orange", ":")]:
                sig_mask = df_rand["q-value"] < fdr_cut
                if sig_mask.any():
                    pval_threshold = df_rand.loc[sig_mask, "P-value"].max()
                    ax.axhline(-np.log10(pval_threshold), color=color, ls=ls, lw=1.5, alpha=0.8,
                               label=f"FDR={fdr_cut}")
            ax.legend(fontsize=7, loc="upper right", framealpha=0.8)

for j in range(i + 1, len(axes2)):
    axes2[j].set_visible(False)

fig2.suptitle("Cell-type mutation bias significance (-log10P, random null)",
              fontsize=16, fontweight="bold", y=1.01)
fig2.tight_layout()
fig2.patch.set_alpha(0)
fig2.savefig(str(FIG_DIR / "all_traits_random_logP_boxplots.png"),
             dpi=200, bbox_inches='tight', transparent=True)
plt.show()

# %% [markdown]
# ## 8. Summary and Recommendation

# %%
print("=" * 70)
print("RECOMMENDATION FOR MAIN PAPER NEGATIVE CONTROLS")
print("=" * 70)
print()
print("Best negative controls (uniform weights, no significant neuronal enrichment):")
print()
for label, df in pval_results.items():
    if label in ["ASD w/o ID", "SCZ", "DDD"]:
        continue
    sig_neur = df[(df["q-value"] < 0.1) & (df["Class"].isin(["NEUR", "Neuron"]))]
    status = "CLEAN" if len(sig_neur) == 0 else f"HAS {len(sig_neur)} SIG NEURONAL"
    sc_mean = df.groupby("Supercluster")["EFFECT"].mean().sort_values(ascending=False)
    cge_rank = list(sc_mean.index).index("CGE interneuron") + 1 if "CGE interneuron" in sc_mean.index else None
    print(f"  [{status}] {label}: top={sc_mean.index[0]}, CGE rank={cge_rank}")

# %%
