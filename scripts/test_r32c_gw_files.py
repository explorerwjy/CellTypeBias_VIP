"""
Quick test: reproduce R3.2c using pipeline gene weight files
SCZ: SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw (461 genes, all positive weights)
DDD: DDD.top285.gw.bgmr.csv (285 genes, all positive weights)

Compare with the current code's inline-constructed weights.
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

# Load expression matrix
expression_matrix = _cfg['analysis_types']['Centering']
HumanCT_Z2_HCT = pd.read_csv(str(PROJ_DIR / expression_matrix), index_col=0)
HumanCT_Z2_HCT.columns = HumanCT_Z2_HCT.columns.astype(int)
print(f"Expression matrix: {HumanCT_Z2_HCT.shape}")

GeneWeightDIR = str(PROJ_DIR / "dat" / "GeneWeights") + "/"

# Load gene weight files
SCZ_GW = pd.read_csv(f"{GeneWeightDIR}/SCZ.top500.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw",
                      index_col=0, header=None, names=["Weight"])
DDD_GW = pd.read_csv(f"{GeneWeightDIR}/DDD.top285.gw.bgmr.csv",
                      index_col=0, header=None, names=["Weight"])

print(f"SCZ: {len(SCZ_GW)} genes, min weight={SCZ_GW['Weight'].min():.3f}, max={SCZ_GW['Weight'].max():.3f}")
print(f"DDD: {len(DDD_GW)} genes, min weight={DDD_GW['Weight'].min():.3f}, max={DDD_GW['Weight'].max():.3f}")
print(f"SCZ zeros: {(SCZ_GW['Weight'] == 0).sum()}, DDD zeros: {(DDD_GW['Weight'] == 0).sum()}")
print(f"SCZ negatives: {(SCZ_GW['Weight'] < 0).sum()}, DDD negatives: {(DDD_GW['Weight'] < 0).sum()}")

# --- R3.2c analysis ---
n_add_values = np.arange(0, 150, 10)
n_reps = 100

disorder_configs = [
    ("SCZ", SCZ_GW, "#ff7f0e"),
    ("DDD", DDD_GW, "#2ca02c"),
]

all_results = {}

for disorder_name, gw_df, color in disorder_configs:
    print(f"\n{'='*50}")
    print(f"Processing {disorder_name} ({len(gw_df)} genes available)")

    n_ref = min(61, len(gw_df))
    ref_gw = dict(zip(gw_df.index[:n_ref], gw_df["Weight"].values[:n_ref]))
    ref_bias = HumanCT_AvgZ_Weighted(HumanCT_Z2_HCT, ref_gw)
    ref_bias = AnnotateCTDat(ref_bias, Anno)
    ref_effect = ref_bias["EFFECT"].values
    common_idx = ref_bias.index

    # Gene pool for random: exclude top-200 of this disorder
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

    # Random gene additions (weight-transferred)
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

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150, sharey=True)

for ax, (disorder_name, gw_df, color) in zip(axes, disorder_configs):
    res = all_results[disorder_name]

    ax.fill_between(total_genes, res["rand_lo"], res["rand_hi"],
                    color="#999999", alpha=0.25, label="Random genes (95% CI)")
    ax.plot(total_genes, res["rand_mean"], color="#999999", lw=2, ls="--",
            marker="s", markersize=3, label="Random genes (mean)")
    ax.plot(total_genes, res["real"], color=color, lw=2.5,
            marker="o", markersize=5, label=f"Ranked {disorder_name} genes", zorder=5)

    ax.axvline(61, color="gray", ls=":", lw=1.5, alpha=0.5)
    ax.set_xlabel("Total number of genes", fontsize=12)
    ax.set_title(disorder_name, fontweight="bold", fontsize=14)
    ax.legend(fontsize=8, framealpha=0.8, loc="lower left")
    ax.set_ylim(0.3, 1.02)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

axes[0].set_ylabel("Spearmans' R with top-61 bias profile", fontsize=12)

fig.suptitle("R3.2c test — Using pipeline gene weight files (all positive weights)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.patch.set_alpha(0)

outpath = str(PROJ_DIR / "results" / "figures" / "test_r32c_pipeline_gw.png")
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=150, bbox_inches='tight', transparent=True, facecolor='none')
print(f"\nSaved: {outpath}")

# Print summary table
for disorder_name, _, _ in disorder_configs:
    res = all_results[disorder_name]
    print(f"\n=== {disorder_name}: Correlation with top-61 reference ===")
    print(f"{'N total':>8s}  {'Real':>7s}  {'Rand mean':>10s}  {'Rand 95%CI':>14s}")
    for i, n_add in enumerate(n_add_values):
        nt = 61 + n_add
        if n_add == 0:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {'—':>10s}  {'—':>14s}")
        else:
            print(f"{nt:8d}  {res['real'][i]:7.3f}  {res['rand_mean'][i]:10.3f}  "
                  f"[{res['rand_lo'][i]:.3f}, {res['rand_hi'][i]:.3f}]")
