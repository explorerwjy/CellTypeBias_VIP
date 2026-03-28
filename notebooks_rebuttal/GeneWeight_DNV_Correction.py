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

# %% [markdown]
# # Gene Weight Correction: Background Mutation Rate (BGMR) Subtraction
#
# **Motivation (Reviewer 3, R3.1):** Gene length is a confounder for de novo
# mutation counts — longer genes accumulate more mutations by chance.  We
# subtract expected de novo counts (Samocha et al. background mutation rates)
# from observed counts before computing gene weights.
#
# **Formula:**
# `corrected_weight = (n_obs_LGD − n_exp_LGD) × wLGD + (n_obs_Dmis − n_exp_Dmis) × wMis`
#
# where `n_expected = p_mutation_rate × 2 × N_probands` (factor of 2 for diploid).
#
# **Gene sets corrected here:**
# - ASD_All (EWS): pLI-stratified weights, N = 42,607
# - ASD_HIQ (IQ > 70): equal weights, N = 4,876
# - ASD_LIQ (IQ ≤ 70): equal weights, N = 2,619
# - DDD top-61: equal weights, N = 31,058 (was incorrectly 42,607)
# - DDD top-285: all 285 genome-wide significant genes (Kaplanis et al. 2020), N = 31,058
#
# SCZ is already corrected (case–control expected count subtraction).

# %% [markdown]
# ## A. Setup & BGMR

# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv

import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))

GW_DIR = PROJ_DIR / "dat" / "GeneWeights"

mpl.rcParams['figure.facecolor'] = 'none'
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['savefig.facecolor'] = 'none'
mpl.rcParams['font.size'] = 12

# %%
# Load BGMR (Samocha et al. per-gene mutation rates)
bgmr = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv", low_memory=False)
bgmr["entrez_id"] = bgmr["entrez_id"].astype(int)
bgmr = bgmr.set_index("entrez_id")
print(f"BGMR: {len(bgmr)} genes")
print(f"Key rate columns: p_LGD, p_misense, prevel_0.5")
print(bgmr[["GeneName", "p_LGD", "p_misense", "prevel_0.5"]].head())

# %%
# Load expression matrix gene list (for filtering to genes in our analysis)
spec = pd.read_csv(PROJ_DIR / "dat" / "ExpMats" / "HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
                    index_col=0, nrows=0)
expr_genes = set(spec.index.astype(int))
print(f"Expression matrix genes: {len(expr_genes)}")

# %% [markdown]
# ## B. ASD_All (EWS): Correct with BGMR
#
# Uses pLI-stratified weights (same as original `Aggregate_Gene_Weights2`):
# - Constrained (pLI ≥ 0.5): wLGD = 0.554, wDmis = 0.333
# - Unconstrained (pLI < 0.5): wLGD = 0.138, wDmis = 0.130

# %%
ews = pd.read_excel(PROJ_DIR / "dat" / "suppl.data" / "41588_2022_1148_MOESM4_ESM.xlsx",
                     sheet_name="Table S7", header=2)
print(f"EWS table: {len(ews)} genes")

# Filter to top EWS genes (same as original: pDenovoWEST_Meta ≤ 0.01)
ews["pDenovoWEST_Meta"] = pd.to_numeric(ews["pDenovoWEST_Meta"], errors="coerce")
ews_sig = ews[ews["pDenovoWEST_Meta"] <= 0.01].copy()
print(f"Significant EWS genes (p ≤ 0.01): {len(ews_sig)}")

# Load original weights to match the same gene set
old_asd_all = pd.read_csv(GW_DIR / "Spark_Meta_EWS.GeneWeight.csv", header=None,
                           names=["Entrez", "Weight"])
old_entrez = set(old_asd_all["Entrez"].values)
print(f"Original ASD_All gene weight file: {len(old_asd_all)} genes")

# %%
N_ASD_ALL = 42607

gene2weight_asd_all = {}
for _, row in ews_sig.iterrows():
    try:
        g = int(row["EntrezID"])
    except (ValueError, TypeError):
        continue

    # Only keep genes that were in original set
    if g not in old_entrez:
        continue

    obs_lof = float(row["AutismMerged_LoF"])
    obs_dmis = float(row["AutismMerged_Dmis_REVEL0.5"])

    # Get expected counts from BGMR
    if g in bgmr.index:
        exp_lof = bgmr.loc[g, "p_LGD"] * 2 * N_ASD_ALL
        exp_dmis = bgmr.loc[g, "prevel_0.5"] * 2 * N_ASD_ALL  # REVEL ≥ 0.5 rate
    else:
        print(f"  Gene {g} not in BGMR — skipping correction")
        exp_lof = 0
        exp_dmis = 0

    # pLI-stratified weights (same as Aggregate_Gene_Weights2)
    try:
        pLI = float(row["ExACpLI"])
    except (ValueError, TypeError):
        pLI = 0.0

    if pLI >= 0.5:
        wLGD, wDmis = 0.554, 0.333
    else:
        wLGD, wDmis = 0.138, 0.130

    weight = (obs_lof - exp_lof) * wLGD + (obs_dmis - exp_dmis) * wDmis
    gene2weight_asd_all[g] = weight

print(f"Corrected ASD_All: {len(gene2weight_asd_all)} genes")
print(f"Genes with negative weight: {sum(1 for v in gene2weight_asd_all.values() if v <= 0)}")

# Save (sorted by weight descending, matching original format)
out_asd_all = GW_DIR / "Spark_Meta_EWS.GeneWeight.bgmr.csv"
with open(out_asd_all, "w", newline="") as f:
    writer = csv.writer(f)
    for k, v in sorted(gene2weight_asd_all.items(), key=lambda x: x[1], reverse=True):
        writer.writerow([k, v])
print(f"Saved → {out_asd_all}")

# %% [markdown]
# ## C. ASD_HIQ & ASD_LIQ: Correct with BGMR
#
# Uses equal weights (wLGD = 1, wMis = 1), matching current HIQ/LIQ scheme.
# LGD = frameshift + splice_acceptor + splice_donor + stop_gained + stop_lost.
# Dmis = missense with REVEL > 0.5.
#
# Cohort sizes: N_HIQ = 4,876 (ASC 3,249 + SPARK 1,627),
# N_LIQ = 2,619 (ASC 1,572 + SPARK 1,047).

# %%
asd_iq = pd.read_csv(PROJ_DIR / "dat" / "ASD_IQ_Mut.csv")
print(f"ASD IQ mutations: {len(asd_iq)} variants")
print(f"IQ range: {asd_iq['IQ'].min():.0f} – {asd_iq['IQ'].max():.0f}")
print(f"GeneEff values: {asd_iq['GeneEff'].unique()}")

# %%
LGD_EFFECTS = {"frameshift", "splice_acceptor", "splice_donor", "stop_gained", "stop_lost"}
N_HIQ = 4876
N_LIQ = 2619

# Load original weights to match gene sets
old_hiq = pd.read_csv(GW_DIR / "HIQ.top61.nopLI.LGD_Dmis_SameWeight.gw",
                       header=None, names=["Entrez", "Weight"])
old_liq = pd.read_csv(GW_DIR / "LIQ.top61.nopLI.LGD_Dmis_SameWeight.gw",
                       header=None, names=["Entrez", "Weight"])

for label, N_proband, iq_cond, old_gw, out_name in [
    ("HIQ", N_HIQ, asd_iq["IQ"] > 70,  old_hiq, "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    ("LIQ", N_LIQ, asd_iq["IQ"] <= 70, old_liq, "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
]:
    sub = asd_iq[iq_cond].copy()
    print(f"\n--- {label} (N={N_proband}) ---")
    print(f"  Variants: {len(sub)}")

    # Aggregate per gene
    gene_counts = {}
    for _, row in sub.iterrows():
        try:
            g = int(row["Entrez"])
        except (ValueError, TypeError):
            continue
        if g not in gene_counts:
            gene_counts[g] = {"nLGD": 0, "nDmis": 0}
        eff = str(row["GeneEff"]).strip()
        if eff in LGD_EFFECTS:
            gene_counts[g]["nLGD"] += 1
        elif eff == "missense":
            try:
                revel = float(row["REVEL"])
                if revel > 0.5:
                    gene_counts[g]["nDmis"] += 1
            except (ValueError, TypeError):
                pass

    # Compute corrected weights (equal weights: wLGD=1, wMis=1)
    old_entrez_set = set(old_gw["Entrez"].values)
    gene2weight = {}
    for g, counts in gene_counts.items():
        if g not in old_entrez_set:
            continue
        if g in bgmr.index:
            exp_lgd = bgmr.loc[g, "p_LGD"] * 2 * N_proband
            exp_dmis = bgmr.loc[g, "prevel_0.5"] * 2 * N_proband
        else:
            exp_lgd = 0
            exp_dmis = 0
        weight = (counts["nLGD"] - exp_lgd) + (counts["nDmis"] - exp_dmis)
        gene2weight[g] = weight

    print(f"  Corrected genes: {len(gene2weight)}")
    print(f"  Genes with negative weight: {sum(1 for v in gene2weight.values() if v <= 0)}")

    # Save
    out_path = GW_DIR / out_name
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        for k, v in sorted(gene2weight.items(), key=lambda x: x[1], reverse=True):
            writer.writerow([k, v])
    print(f"  Saved → {out_path}")

# %% [markdown]
# ## D. DDD: Fix Nproband (31,058 instead of 42,607)
#
# The existing `DDD.top61.gw.csv` was generated using the default `Nproband=42,607`
# (ASD cohort size), but DDD has 31,058 trios (Kaplanis et al. 2020). This
# over-estimated expected counts by ~37%.

# %%
ddd = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
print(f"DDD table: {len(ddd)} genes")

N_DDD = 31058

# Gene name → Entrez mapping from BGMR
name2entrez = dict(zip(bgmr["GeneName"], bgmr.index))

# Filter significant genes (Bonferroni: 0.05/n_genes)
n_genes_ddd = len(ddd)
bf_thresh = 0.05 / n_genes_ddd
ddd_sig = ddd[ddd["denovoWEST_p_full"] <= bf_thresh].copy()
print(f"Significant DDD genes (p ≤ {bf_thresh:.2e}): {len(ddd_sig)}")

# Map to Entrez IDs
ddd_sig["EntrezID"] = ddd_sig["symbol"].map(name2entrez)
n_mapped = ddd_sig["EntrezID"].notna().sum()
print(f"Mapped to Entrez: {n_mapped}/{len(ddd_sig)}")

# Load old DDD weights
old_ddd = pd.read_csv(GW_DIR / "DDD.top61.gw.csv", header=None, names=["Entrez", "Weight"])
old_ddd_entrez = set(old_ddd["Entrez"].values)

# Compute corrected weights using Aggregate_Gene_Weights_NDD logic
# wLGD=1, wMis=1 (equal weights, matching existing scheme)
# IMPORTANT: keep the same 61 genes as original (selected by p_denovowest),
# only correct their weights for background mutation rate.
gene2weight_ddd = {}
for _, row in ddd_sig.iterrows():
    try:
        g = int(row["EntrezID"])
    except (ValueError, TypeError):
        continue

    # Only keep genes from original top-61 set (ranked by p_denovowest)
    if g not in old_ddd_entrez:
        continue

    nLGD = (row["frameshift_variant"] + row["splice_acceptor_variant"]
            + row["splice_donor_variant"] + row["stop_gained"] + row["stop_lost"])
    nMis = row["missense_variant"]

    if g in bgmr.index:
        exp_lgd = bgmr.loc[g, "p_LGD"] * 2 * N_DDD
        exp_mis = bgmr.loc[g, "p_misense"] * 2 * N_DDD  # all missense for DDD
    else:
        exp_lgd = 0
        exp_mis = 0

    weight = (nLGD - exp_lgd) * 1 + (nMis - exp_mis) * 1  # equal weights
    gene2weight_ddd[g] = weight

top_ddd = sorted(gene2weight_ddd.items(), key=lambda x: x[1], reverse=True)
print(f"\nCorrected DDD: {len(gene2weight_ddd)} genes (same membership as original)")
print(f"Genes with negative weight: {sum(1 for v in gene2weight_ddd.values() if v <= 0)}")

# Save
out_ddd = GW_DIR / "DDD.top61.gw.bgmr.csv"
with open(out_ddd, "w", newline="") as f:
    writer = csv.writer(f)
    for k, v in top_ddd:
        writer.writerow([k, v])
print(f"Saved → {out_ddd}")

# %% [markdown]
# ## D2. DDD top-285: All significant genes with BGMR correction
#
# Kaplanis et al. (2020) report 285 genome-wide significant genes.
# The previous `DDD.hc.gw` used `Nproband=42,607` (wrong — ASD default).
# Here we re-generate with correct `N=31,058` and take top 285 by significance.

# %%
# Use same ddd DataFrame loaded in section D
N_DDD_285 = 285

# Sort by significance, take top 285 (as reported in the paper)
ddd_top285 = ddd.sort_values("denovoWEST_p_full").head(N_DDD_285).copy()
ddd_top285["EntrezID"] = ddd_top285["symbol"].map(name2entrez)
n_mapped_285 = ddd_top285["EntrezID"].notna().sum()
print(f"DDD top 285 genes, mapped to Entrez: {n_mapped_285}/{N_DDD_285}")

# Compute BGMR-corrected weights (wLGD=1, wMis=1, N=31,058)
gene2weight_ddd285 = {}
for _, row in ddd_top285.iterrows():
    try:
        g = int(row["EntrezID"])
    except (ValueError, TypeError):
        continue

    nLGD = (row["frameshift_variant"] + row["splice_acceptor_variant"]
            + row["splice_donor_variant"] + row["stop_gained"] + row["stop_lost"])
    nMis = row["missense_variant"]

    if g in bgmr.index:
        exp_lgd = bgmr.loc[g, "p_LGD"] * 2 * N_DDD
        exp_mis = bgmr.loc[g, "p_misense"] * 2 * N_DDD
    else:
        exp_lgd = 0
        exp_mis = 0

    weight = (nLGD - exp_lgd) + (nMis - exp_mis)
    gene2weight_ddd285[g] = weight

print(f"Corrected DDD-285: {len(gene2weight_ddd285)} genes")
print(f"Genes with negative weight: {sum(1 for v in gene2weight_ddd285.values() if v <= 0)}")

# Save (sorted by weight descending)
out_ddd285 = GW_DIR / "DDD.top285.gw.bgmr.csv"
with open(out_ddd285, "w", newline="") as f:
    writer = csv.writer(f)
    for k, v in sorted(gene2weight_ddd285.items(), key=lambda x: x[1], reverse=True):
        writer.writerow([k, v])
print(f"Saved → {out_ddd285}")

# %% [markdown]
# ## E. Comparison: Old vs Corrected Weights

# %%
from scipy.stats import spearmanr

fig, axes = plt.subplots(1, 5, figsize=(24, 5))

comparisons = [
    ("ASD_All", "Spark_Meta_EWS.GeneWeight.csv", "Spark_Meta_EWS.GeneWeight.bgmr.csv"),
    ("ASD_HIQ", "HIQ.top61.nopLI.LGD_Dmis_SameWeight.gw", "HIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    ("ASD_LIQ", "LIQ.top61.nopLI.LGD_Dmis_SameWeight.gw", "LIQ.top61.nopLI.LGD_Dmis_SameWeight.bgmr.gw"),
    ("DDD_61", "DDD.top61.gw.csv", "DDD.top61.gw.bgmr.csv"),
    ("DDD_285", "DDD.hc.gw.csv", "DDD.top285.gw.bgmr.csv"),
]

for ax, (label, old_file, new_file) in zip(axes, comparisons):
    old = pd.read_csv(GW_DIR / old_file, header=None, names=["Entrez", "OldWeight"])
    new = pd.read_csv(GW_DIR / new_file, header=None, names=["Entrez", "NewWeight"])

    merged = old.merge(new, on="Entrez", how="inner")

    if len(merged) > 2:
        rho, pval = spearmanr(merged["OldWeight"], merged["NewWeight"])
    else:
        rho, pval = float("nan"), float("nan")

    ax.scatter(merged["OldWeight"], merged["NewWeight"], s=20, alpha=0.7)

    # Diagonal reference
    lo = min(merged["OldWeight"].min(), merged["NewWeight"].min())
    hi = max(merged["OldWeight"].max(), merged["NewWeight"].max())
    margin = (hi - lo) * 0.05
    ax.plot([lo - margin, hi + margin], [lo - margin, hi + margin],
            ls="--", lw=1, color="grey")

    ax.set_xlabel("Original weight")
    ax.set_ylabel("BGMR-corrected weight")
    ax.set_title(f"{label}\n(n={len(merged)}, $\\rho$={rho:.3f})", fontweight="bold")

fig.suptitle("Gene Weight Comparison: Original vs BGMR-Corrected", fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
plt.show()

# %%
# Summary table: rank changes
print("=" * 70)
print("SUMMARY: Gene Rank Changes After BGMR Correction")
print("=" * 70)

for label, old_file, new_file in comparisons:
    old = pd.read_csv(GW_DIR / old_file, header=None, names=["Entrez", "OldWeight"])
    new = pd.read_csv(GW_DIR / new_file, header=None, names=["Entrez", "NewWeight"])

    old["OldRank"] = range(1, len(old) + 1)
    new["NewRank"] = range(1, len(new) + 1)

    merged = old.merge(new, on="Entrez", how="inner")
    merged["RankChange"] = merged["OldRank"] - merged["NewRank"]
    merged["WeightPctChange"] = ((merged["NewWeight"] - merged["OldWeight"]) / merged["OldWeight"] * 100)

    rho, _ = spearmanr(merged["OldWeight"], merged["NewWeight"])

    print(f"\n{label}:")
    print(f"  Genes in common: {len(merged)}")
    print(f"  Spearman rho: {rho:.3f}")
    print(f"  Max rank change: {merged['RankChange'].abs().max()}")

    # Flag genes with > 50% weight reduction
    big_changes = merged[merged["WeightPctChange"].abs() > 50]
    if len(big_changes) > 0:
        print(f"  Genes with >50% weight change:")
        for _, r in big_changes.iterrows():
            print(f"    Entrez {int(r['Entrez'])}: {r['OldWeight']:.3f} → {r['NewWeight']:.3f} "
                  f"({r['WeightPctChange']:+.1f}%)")
    else:
        print(f"  No genes with >50% weight change")
