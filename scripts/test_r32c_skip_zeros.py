"""
R3.2c test: use inline-constructed weights (same as current notebook)
but completely remove zero-weight genes from the DataFrame before analysis.
Compare with the pipeline .gw file version.
"""
import sys
import os
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from joblib import Parallel, delayed
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

expression_matrix = _cfg['analysis_types']['Centering']
HumanCT_Z2_HCT = pd.read_csv(str(PROJ_DIR / expression_matrix), index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
print(f"Expression matrix: {HumanCT_Z2_HCT.shape}")

GeneWeightDIR = str(PROJ_DIR / "dat" / "GeneWeights") + "/"

# Load BGMR
BGMR = pd.read_csv("/home/jw3514/Work/Resources/BGMR.withEntrez.csv", low_memory=False)
BGMR["entrez_id"] = BGMR["entrez_id"].astype(int)
BGMR = BGMR.set_index("entrez_id")

# --- SCZ: same file as before (already all positive) ---
SCZ_GW = pd.read_csv(f"{GeneWeightDIR}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
                      index_col=0, header=None, names=["Weight"])

# --- DDD: construct inline (same as notebook), then drop zeros ---
DDD_Genes = pd.read_excel("/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx")
DDD_Genes = DDD_Genes.sort_values("denovoWEST_p_full")
entrez_ids = [int(GeneSymbol2Entrez.get(x, -1)) for x in DDD_Genes["symbol"].values]
DDD_Genes["EntrezID"] = entrez_ids

DDD_Genes_sub = DDD_Genes[["symbol", "EntrezID", "denovoWEST_p_full",
                            "missense_variant", "frameshift_variant",
                            "splice_acceptor_variant", "splice_donor_variant",
                            "stop_gained", "stop_lost"]].reset_index(drop=True)

N_DDD = 31058

def _ddd_gene_weights(df, Nproband, BGMR):
    valid = df[df["EntrezID"].isin(BGMR.index)]
    if len(valid) > 0:
        return Aggregate_Gene_Weights_NDD(valid, Nproband=Nproband, BGMR=BGMR, wLGD=1, wMis=1)
    return {}

# Compute weights for ALL DDD genes
_ddd_all_gw = _ddd_gene_weights(DDD_Genes_sub, Nproband=N_DDD, BGMR=BGMR)
DDD_GW_inline = pd.DataFrame([
    {"EntrezID": row["EntrezID"], "Weight": _ddd_all_gw.get(row["EntrezID"], 0)}
    for _, row in DDD_Genes_sub.iterrows()
    if row["EntrezID"] > 0
]).set_index("EntrezID")
DDD_GW_inline["Weight"] = DDD_GW_inline["Weight"].clip(lower=0)

print(f"\nDDD inline (before filtering): {len(DDD_GW_inline)} genes")
print(f"  zeros: {(DDD_GW_inline['Weight'] == 0).sum()}")
print(f"  negatives (clipped): {(DDD_GW_inline['Weight'] < 0).sum()}")
print(f"  positive: {(DDD_GW_inline['Weight'] > 0).sum()}")

# Drop zeros
DDD_GW_nozero = DDD_GW_inline[DDD_GW_inline["Weight"] > 0].copy()
print(f"DDD inline (after dropping zeros): {len(DDD_GW_nozero)} genes")
print(f"  min weight={DDD_GW_nozero['Weight'].min():.3f}, max={DDD_GW_nozero['Weight'].max():.3f}")

# For comparison, also load the pipeline file
DDD_GW_pipeline = pd.read_csv(f"{GeneWeightDIR}/DDD.top285.gw.bgmr.csv",
                               index_col=0, header=None, names=["Weight"])
print(f"\nDDD pipeline file: {len(DDD_GW_pipeline)} genes, all positive")

# --- R3.2c for both versions ---
n_add_values = np.arange(0, 150, 10)
n_reps = 100

configs = [
    ("SCZ (pipeline)",      SCZ_GW,          "#ff7f0e"),
    ("DDD (inline, skip 0)", DDD_GW_nozero,  "#2ca02c"),
    ("DDD (pipeline file)",  DDD_GW_pipeline, "#1f77b4"),
]

all_results = {}

for disorder_name, gw_df, color in configs:
    print(f"\n{'='*50}")
    print(f"Processing {disorder_name} ({len(gw_df)} genes available)")

    n_ref = min(61, len(gw_df))
    ref_gw = dict(zip(gw_df.index[:n_ref], gw_df["Weight"].values[:n_ref]))
    ref_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_gw)
    ref_bias = AnnotateCTDat(ref_bias, Anno)
    ref_effect = ref_bias["EFFECT"].values
    common_idx = ref_bias.index

    n_exclude = min(200, len(gw_df))
    exclude_set = set(gw_df.index[:n_exclude].astype(int))
    gene_pool = np.array([g for g in HumanCT_Z2_HCT.index.values if g not in exclude_set])

    all_weights = gw_df["Weight"].values

    # Real ranked genes
    real_corrs = []
    for n_add in n_add_values:
        n_total = min(61 + n_add, len(gw_df))
        gw = dict(zip(gw_df.index[:n_total], all_weights[:n_total]))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        real_corrs.append(r)

    # Random
    def _one_random_trial(n_add_trial, seed, top61_idx, top61_wt, gene_pool, rank_weights):
        rng = np.random.default_rng(seed)
        rand_genes = rng.choice(gene_pool, size=n_add_trial, replace=False)
        gw = dict(zip(top61_idx, top61_wt))
        gw.update(dict(zip(rand_genes, rank_weights[:n_add_trial])))
        bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, gw)
        bias = AnnotateCTDat(bias, Anno)
        r, _ = spearmanr(ref_effect, bias.loc[common_idx, "EFFECT"].values)
        return r

    top61_idx = gw_df.index[:n_ref]
    top61_wt = gw_df["Weight"].values[:n_ref]

    rand_mean = [1.0]
    rand_lo = [1.0]
    rand_hi = [1.0]

    for n_add in n_add_values[1:]:
        n_total = min(61 + n_add, len(gw_df))
        rank_weights = all_weights[61:n_total]
        if len(rank_weights) < n_add:
            rank_weights = np.concatenate([
                rank_weights,
                np.full(n_add - len(rank_weights), rank_weights[-1] if len(rank_weights) > 0 else 1.0)
            ])
        rs = Parallel(n_jobs=10)(
            delayed(_one_random_trial)(
                n_add, seed=hash((disorder_name, n_add, rep)) % (2**31),
                top61_idx=top61_idx, top61_wt=top61_wt,
                gene_pool=gene_pool, rank_weights=rank_weights,
            )
            for rep in range(n_reps)
        )
        rs = np.array(rs)
        rand_mean.append(rs.mean())
        rand_lo.append(np.percentile(rs, 2.5))
        rand_hi.append(np.percentile(rs, 97.5))

    all_results[disorder_name] = {
        "real": np.array(real_corrs),
        "rand_mean": np.array(rand_mean),
        "rand_lo": np.array(rand_lo),
        "rand_hi": np.array(rand_hi),
        "color": color,
    }
    print(f"  Done. Real ρ at N=201: {real_corrs[-1]:.3f}, Random mean: {rand_mean[-1]:.3f}")

# --- Plot ---
total_genes = 61 + n_add_values

fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=150, sharey=True)

for ax, (disorder_name, gw_df, color) in zip(axes, configs):
    res = all_results[disorder_name]

    ax.fill_between(total_genes, res["rand_lo"], res["rand_hi"],
                    color="#999999", alpha=0.25, label="Random genes (95% CI)")
    ax.plot(total_genes, res["rand_mean"], color="#999999", lw=2, ls="--",
            marker="s", markersize=3, label="Random genes (mean)")
    ax.plot(total_genes, res["real"], color=color, lw=2.5,
            marker="o", markersize=5, label=f"Ranked genes", zorder=5)

    ax.axvline(61, color="gray", ls=":", lw=1.5, alpha=0.5)
    ax.set_xlabel("Total number of genes", fontsize=12)
    ax.set_title(disorder_name, fontweight="bold", fontsize=14)
    ax.legend(fontsize=8, framealpha=0.8, loc="lower left")
    ax.set_ylim(0.3, 1.02)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

axes[0].set_ylabel("Spearman ρ with top-61 bias profile", fontsize=12)

fig.suptitle("R3.2c — SCZ (pipeline) vs DDD (inline skip-zeros vs pipeline file)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.patch.set_alpha(0)

outpath = str(PROJ_DIR / "results" / "figures" / "test_r32c_skip_zeros.png")
fig.savefig(outpath, dpi=150, bbox_inches='tight', transparent=True, facecolor='none')
print(f"\nSaved: {outpath}")

# Summary
for name, _, _ in configs:
    res = all_results[name]
    print(f"\n=== {name} ===")
    print(f"{'N total':>8s}  {'Real':>7s}  {'Rand mean':>10s}  {'Rand 95%CI':>14s}")
    for i, n_add in enumerate(n_add_values):
        nt = 61 + n_add
        if n_add == 0:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {'—':>10s}  {'—':>14s}")
        else:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {res['rand_mean'][i]:10.3f}  "
                  f"[{res['rand_lo'][i]:.3f}, {res['rand_hi'][i]:.3f}]")
