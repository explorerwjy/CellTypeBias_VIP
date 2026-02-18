#!/usr/bin/env python
"""
Downsampling Analysis for Mutation Bias Stability (TRA-18)

This script performs downsampling analysis to demonstrate the minimum sample size
needed for stable mutation bias estimates. At each downsampled fraction, it:
1. Subsamples raw mutation counts using binomial distribution
2. Re-runs gene discovery (Poisson test for de novo, Fisher exact for case-control)
3. Selects top-61 genes from subsampled data
4. Computes gene weights with BGMR correction
5. Computes cell-type bias and saves results

Usage:
    python scripts/script_downsampling.py \
        --disorder ASD \
        --fraction 0.5 \
        --n_iter 100 \
        --n_jobs 10 \
        --seed 42 \
        --outdir results/downsampling/
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from CellType_PSY import (
    HumanCT_AvgZ_Weighted,
    poisson_test_denovo,
    fisher_test_case_control,
)
from UNIMED import LoadGeneINFO


# ============================================================================
# Data Paths (absolute paths on lab workstation)
# ============================================================================
DATA_PATHS = {
    "ASD_mutations": "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/41588_2022_1148_MOESM4_ESM.xlsx",
    "SCZ_mutations": "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/SCZ.ALLGENE.MutCountModified.csv",
    "DDD_mutations": "/home/jw3514/Work/data/DDD/41586_2020_2832_MOESM4_ESM.xlsx",
    "BGMR": "/home/jw3514/Work/Resources/BGMR.withEntrez.csv",
    "expression_matrix": "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv",
    "annotation": "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/dat/annotation.xlsx",
}

# Sample sizes
SAMPLE_SIZES = {
    "ASD": 42607,
    "SCZ_case": 24248,
    "SCZ_ctrl": 97322,
    "DDD": 31058,
}

# Number of top genes to select for bias calculation
TOP_N_GENES = 61


def load_bgmr():
    """Load BGMR (background mutation rate) data."""
    bgmr = pd.read_csv(DATA_PATHS["BGMR"], index_col="entrez_id")
    return bgmr


def load_asd_data():
    """Load ASD mutation data from SPARK meta-analysis Table S7."""
    df = pd.read_excel(DATA_PATHS["ASD_mutations"], sheet_name="Table S7", header=2)
    # Required columns: EntrezID, AutismMerged_LoF, AutismMerged_Dmis_REVEL0.5
    df = df[["EntrezID", "AutismMerged_LoF", "AutismMerged_Dmis_REVEL0.5", "ExACpLI"]].copy()
    df = df.dropna(subset=["EntrezID"])
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")
    # Fill NaN mutation counts with 0
    df["AutismMerged_LoF"] = df["AutismMerged_LoF"].fillna(0).astype(int)
    df["AutismMerged_Dmis_REVEL0.5"] = df["AutismMerged_Dmis_REVEL0.5"].fillna(0).astype(int)
    return df


def load_scz_data():
    """Load SCZ case-control mutation data."""
    df = pd.read_csv(DATA_PATHS["SCZ_mutations"])
    # Required columns: Entrez, Case PTV, Ctrl PTV, Case mis3, Ctrl mis3
    df = df.rename(columns={"Entrez": "EntrezID"})
    df = df.dropna(subset=["EntrezID"])
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")
    # Fill NaN mutation counts with 0
    for col in ["Case PTV", "Ctrl PTV", "Case mis3", "Ctrl mis3"]:
        df[col] = df[col].fillna(0).astype(int)
    return df


def load_ddd_data(bgmr):
    """Load DDD mutation data from Kaplanis et al. 2020."""
    df = pd.read_excel(DATA_PATHS["DDD_mutations"], sheet_name="kaplanis_samocha_denovoWEST_res")

    # Map gene symbols to Entrez IDs using BGMR GeneName column
    symbol_to_entrez = dict(zip(bgmr["GeneName"].values, bgmr.index.values))

    df["EntrezID"] = df["symbol"].map(symbol_to_entrez)
    df = df.dropna(subset=["EntrezID"])
    df["EntrezID"] = df["EntrezID"].astype(int)
    df = df.set_index("EntrezID")

    # Compute LGD count
    lgd_cols = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "stop_lost",
    ]
    for col in lgd_cols:
        if col not in df.columns:
            df[col] = 0
    df["nLGD"] = df[lgd_cols].sum(axis=1).fillna(0).astype(int)
    df["nMis"] = df["missense_variant"].fillna(0).astype(int)

    return df


def load_expression_matrix():
    """Load the 461-cell-type expression matrix."""
    return pd.read_csv(DATA_PATHS["expression_matrix"], index_col=0)


def load_annotation():
    """Load cell type annotation."""
    return pd.read_excel(DATA_PATHS["annotation"], index_col=0)


# ============================================================================
# Downsampling Functions
# ============================================================================


def run_single_iteration_asd(
    asd_data, bgmr, expr_mat, fraction, seed, N_full=SAMPLE_SIZES["ASD"]
):
    """
    Run a single downsampling iteration for ASD.

    1. Subsample mutation counts using binomial distribution
    2. Run gene discovery using Poisson test
    3. Select top-61 genes
    4. Compute gene weights with BGMR correction
    5. Compute cell-type bias
    """
    rng = np.random.default_rng(seed)
    N_sub = int(fraction * N_full)

    results = []

    for entrez_id, row in asd_data.iterrows():
        if entrez_id not in bgmr.index:
            continue

        # 1. Subsample de novo counts (Binomial)
        lof_orig = int(row["AutismMerged_LoF"])
        dmis_orig = int(row["AutismMerged_Dmis_REVEL0.5"])

        lof_sub = rng.binomial(n=lof_orig, p=fraction) if lof_orig > 0 else 0
        dmis_sub = rng.binomial(n=dmis_orig, p=fraction) if dmis_orig > 0 else 0

        # 2. Gene discovery: Poisson test
        p_lgd_rate = bgmr.loc[entrez_id, "p_LGD"]
        p_dmis_rate = bgmr.loc[entrez_id, "prevel_0.5"]

        p_combined, p_lgd, p_dmis = poisson_test_denovo(
            lof_sub, dmis_sub, p_lgd_rate, p_dmis_rate, N_sub
        )

        # 3. Gene weight (BGMR-corrected, equal weights)
        exp_lof = p_lgd_rate * 2 * N_sub
        exp_dmis = p_dmis_rate * 2 * N_sub
        weight = (lof_sub - exp_lof) * 1.0 + (dmis_sub - exp_dmis) * 1.0

        results.append(
            {
                "EntrezID": entrez_id,
                "p_combined": p_combined,
                "p_lgd": p_lgd,
                "p_dmis": p_dmis,
                "weight": weight,
                "lof_sub": lof_sub,
                "dmis_sub": dmis_sub,
            }
        )

    # Convert to DataFrame and select top-61 genes
    df = pd.DataFrame(results)
    df = df.sort_values("p_combined", ascending=True).head(TOP_N_GENES)

    # Create gene weights dict (only positive weights)
    gene_weights = {
        int(row["EntrezID"]): row["weight"]
        for _, row in df.iterrows()
        if row["weight"] > 0
    }

    # If too few genes with positive weights, use all top genes with weight > 0
    if len(gene_weights) < 10:
        gene_weights = {
            int(row["EntrezID"]): max(row["weight"], 0.01)
            for _, row in df.iterrows()
        }

    # 5. Compute cell-type bias
    if len(gene_weights) == 0:
        return None

    bias = HumanCT_AvgZ_Weighted(expr_mat, gene_weights)
    return bias["EFFECT"]


def run_single_iteration_scz(
    scz_data,
    expr_mat,
    fraction,
    seed,
    N_case_full=SAMPLE_SIZES["SCZ_case"],
    N_ctrl=SAMPLE_SIZES["SCZ_ctrl"],
):
    """
    Run a single downsampling iteration for SCZ.

    Downsamples CASES ONLY, keeps controls fixed.
    """
    rng = np.random.default_rng(seed)
    N_case_sub = int(fraction * N_case_full)

    results = []

    for entrez_id, row in scz_data.iterrows():
        # 1. Subsample case counts only (Binomial)
        case_ptv_orig = int(row["Case PTV"])
        case_mis3_orig = int(row["Case mis3"])
        ctrl_ptv = int(row["Ctrl PTV"])
        ctrl_mis3 = int(row["Ctrl mis3"])

        case_ptv_sub = rng.binomial(n=case_ptv_orig, p=fraction) if case_ptv_orig > 0 else 0
        case_mis3_sub = rng.binomial(n=case_mis3_orig, p=fraction) if case_mis3_orig > 0 else 0

        # 2. Gene discovery: Fisher exact test
        p_combined, p_lgd, p_dmis = fisher_test_case_control(
            case_ptv_sub, case_mis3_sub, ctrl_ptv, ctrl_mis3, N_case_sub, N_ctrl
        )

        # 3. Gene weight (mutation count mode, equal weights)
        weight = case_ptv_sub * 1.0 + case_mis3_sub * 1.0

        results.append(
            {
                "EntrezID": entrez_id,
                "p_combined": p_combined,
                "p_lgd": p_lgd,
                "p_dmis": p_dmis,
                "weight": weight,
                "case_ptv_sub": case_ptv_sub,
                "case_mis3_sub": case_mis3_sub,
            }
        )

    # Convert to DataFrame and select top-61 genes
    df = pd.DataFrame(results)
    df = df.sort_values("p_combined", ascending=True).head(TOP_N_GENES)

    # Create gene weights dict
    gene_weights = {
        int(row["EntrezID"]): row["weight"]
        for _, row in df.iterrows()
        if row["weight"] > 0
    }

    if len(gene_weights) < 10:
        gene_weights = {
            int(row["EntrezID"]): max(row["weight"], 0.01)
            for _, row in df.iterrows()
        }

    if len(gene_weights) == 0:
        return None

    bias = HumanCT_AvgZ_Weighted(expr_mat, gene_weights)
    return bias["EFFECT"]


def run_single_iteration_ddd(
    ddd_data, bgmr, expr_mat, fraction, seed, N_full=SAMPLE_SIZES["DDD"]
):
    """
    Run a single downsampling iteration for DDD.

    Same approach as ASD but with different data source.
    """
    rng = np.random.default_rng(seed)
    N_sub = int(fraction * N_full)

    results = []

    for entrez_id, row in ddd_data.iterrows():
        if entrez_id not in bgmr.index:
            continue

        # 1. Subsample de novo counts (Binomial)
        lgd_orig = int(row["nLGD"])
        mis_orig = int(row["nMis"])

        lgd_sub = rng.binomial(n=lgd_orig, p=fraction) if lgd_orig > 0 else 0
        mis_sub = rng.binomial(n=mis_orig, p=fraction) if mis_orig > 0 else 0

        # 2. Gene discovery: Poisson test
        # Use p_misense for DDD (not prevel_0.5 since DDD uses all missense, not just REVEL-filtered)
        p_lgd_rate = bgmr.loc[entrez_id, "p_LGD"]
        p_mis_rate = bgmr.loc[entrez_id, "p_misense"]

        p_combined, p_lgd, p_dmis = poisson_test_denovo(
            lgd_sub, mis_sub, p_lgd_rate, p_mis_rate, N_sub
        )

        # 3. Gene weight (BGMR-corrected, equal weights)
        exp_lgd = p_lgd_rate * 2 * N_sub
        exp_mis = p_mis_rate * 2 * N_sub
        weight = (lgd_sub - exp_lgd) * 1.0 + (mis_sub - exp_mis) * 1.0

        results.append(
            {
                "EntrezID": entrez_id,
                "p_combined": p_combined,
                "p_lgd": p_lgd,
                "p_dmis": p_dmis,
                "weight": weight,
                "lgd_sub": lgd_sub,
                "mis_sub": mis_sub,
            }
        )

    # Convert to DataFrame and select top-61 genes
    df = pd.DataFrame(results)
    df = df.sort_values("p_combined", ascending=True).head(TOP_N_GENES)

    # Create gene weights dict
    gene_weights = {
        int(row["EntrezID"]): row["weight"]
        for _, row in df.iterrows()
        if row["weight"] > 0
    }

    if len(gene_weights) < 10:
        gene_weights = {
            int(row["EntrezID"]): max(row["weight"], 0.01)
            for _, row in df.iterrows()
        }

    if len(gene_weights) == 0:
        return None

    bias = HumanCT_AvgZ_Weighted(expr_mat, gene_weights)
    return bias["EFFECT"]


def run_downsampling(disorder, fraction, n_iter, n_jobs, seed, outdir):
    """
    Run downsampling analysis for a specific disorder and fraction.
    """
    print(f"Loading data for {disorder}...")

    # Load common data
    bgmr = load_bgmr()
    expr_mat = load_expression_matrix()

    # Load disorder-specific data
    if disorder == "ASD":
        data = load_asd_data()
        run_func = lambda s: run_single_iteration_asd(data, bgmr, expr_mat, fraction, s)
    elif disorder == "SCZ":
        data = load_scz_data()
        run_func = lambda s: run_single_iteration_scz(data, expr_mat, fraction, s)
    elif disorder == "DDD":
        data = load_ddd_data(bgmr)
        run_func = lambda s: run_single_iteration_ddd(data, bgmr, expr_mat, fraction, s)
    else:
        raise ValueError(f"Unknown disorder: {disorder}")

    print(f"Running {n_iter} iterations at fraction {fraction} with {n_jobs} workers...")

    # Generate seeds for each iteration
    seeds = [seed + i for i in range(n_iter)]

    # Run iterations in parallel
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(run_func)(s) for s in seeds
    )

    # Filter out None results
    valid_results = [r for r in results if r is not None]
    print(f"Got {len(valid_results)}/{n_iter} valid iterations")

    if len(valid_results) == 0:
        print("ERROR: No valid results. Check data and parameters.")
        return

    # Combine results into DataFrame
    bias_df = pd.DataFrame(valid_results).T
    bias_df.columns = [f"iter_{i}" for i in range(len(valid_results))]

    # Save results
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"{disorder}_f{fraction:.2f}_bias.csv"
    bias_df.to_csv(outfile)
    print(f"Results saved to: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Downsampling analysis for mutation bias stability"
    )
    parser.add_argument(
        "--disorder",
        type=str,
        required=True,
        choices=["ASD", "SCZ", "DDD"],
        help="Disorder to analyze",
    )
    parser.add_argument(
        "--fraction",
        type=float,
        required=True,
        help="Fraction of data to use (0.1 to 1.0)",
    )
    parser.add_argument(
        "--n_iter",
        type=int,
        default=100,
        help="Number of bootstrap iterations (default: 100)",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=10,
        help="Number of parallel workers (default: 10)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="results/downsampling/",
        help="Output directory (default: results/downsampling/)",
    )

    args = parser.parse_args()

    run_downsampling(
        disorder=args.disorder,
        fraction=args.fraction,
        n_iter=args.n_iter,
        n_jobs=args.n_jobs,
        seed=args.seed,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
