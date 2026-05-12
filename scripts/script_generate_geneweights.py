# Author: jywang	explorerwjy@gmail.com

# ========================================================================================================
# script_run_ctrl_sim.v3.py
# Run control simutations for expression biases (Cell Type; Structures)
# ========================================================================================================

import argparse
import sys
import os
from pathlib import Path
import multiprocessing
from multiprocessing import Pool

# Get project directory from script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ProjDIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(1, os.path.join(ProjDIR, 'src'))

from CellType_PSY import *

# Import matching functions from script_match_genes_by_variables.py
sys.path.insert(1, SCRIPT_DIR)
from script_match_genes_by_variables import (
    load_or_create_percentile_table,
    match_gene,
    convert_to_percentiles
)



def calculate_pairwise_distances(input_pct, candidate_pct, matching_vars):
    """
    Calculate pairwise distances between input genes and candidates.
    Returns: distance matrix (n_input × n_candidates)
    """
    import numpy as np

    pct_cols = [col + '_pct' for col in matching_vars]

    # Extract percentile values
    input_vals = input_pct[pct_cols].values  # (n_input, n_vars)
    candidate_vals = candidate_pct[pct_cols].values  # (n_candidates, n_vars)

    # Calculate pairwise Euclidean distances using broadcasting
    # input_vals[:, None, :] -> (n_input, 1, n_vars)
    # candidate_vals[None, :, :] -> (1, n_candidates, n_vars)
    # Result: (n_input, n_candidates, n_vars) -> sum and sqrt -> (n_input, n_candidates)
    diff_squared = (input_vals[:, None, :] - candidate_vals[None, :, :]) ** 2
    distances = np.sqrt(np.sum(diff_squared, axis=2))

    return distances


def apply_kernel_vectorized(distances, kernel, bandwidth):
    """
    Apply kernel to distance matrix.
    Returns: weight matrix (n_input × n_candidates)
    """
    import numpy as np

    if kernel == 'uniform':
        weights = (distances <= bandwidth).astype(float)
    elif kernel == 'tricubic':
        normalized_dist = distances / bandwidth
        mask = normalized_dist < 1.0
        weights = np.zeros_like(normalized_dist)
        weights[mask] = (1 - normalized_dist[mask]**3)**3
    else:
        raise ValueError(f"Unknown kernel: {kernel}")

    return weights


###########################################################################
## Worker Functions for Parallelization
###########################################################################

def worker_random_genes(sim_idx, valid_genes, n_genes, gene_probs=None, gene_pool=None):
    """Worker function for random gene sampling."""
    import numpy as np
    if gene_probs is not None and gene_pool is not None:
        genes = np.random.choice(gene_pool, size=n_genes, p=gene_probs, replace=False)
    else:
        genes = np.random.choice(valid_genes, size=n_genes, replace=False)
    return sim_idx, genes


def worker_matched_genes(sim_idx, entrez_ids, input_id_to_idx, weight_matrix,
                         candidate_ids_array, valid_genes, base_seed=None):
    """Worker function for gene-by-gene matched sampling."""
    import numpy as np

    # Set unique random seed for this simulation (fixes parallel RNG duplication)
    if base_seed is not None:
        np.random.seed(base_seed + sim_idx)

    # Track which candidates have been used in this simulation
    available_mask = np.ones(len(candidate_ids_array), dtype=bool)
    matched_genes = []

    for gene_idx, gene_id in enumerate(entrez_ids):
        if gene_id not in input_id_to_idx:
            # Gene not in master table, use random sampling
            matched_gene = np.random.choice(valid_genes)
        else:
            # Get pre-computed weights for this input gene
            weight_row_idx = input_id_to_idx[gene_id]
            weights = weight_matrix[weight_row_idx, :].copy()

            # Mask out already-used candidates
            weights[~available_mask] = 0

            # Check if we have valid candidates
            if weights.sum() > 0:
                # Normalize weights to probabilities
                probabilities = weights / weights.sum()

                # Sample one candidate
                matched_idx = np.random.choice(len(candidate_ids_array), p=probabilities)
                matched_gene = candidate_ids_array[matched_idx]

                # Mark as used
                available_mask[matched_idx] = False
            else:
                # No valid matches, fall back to random selection from remaining candidates
                remaining_candidates = candidate_ids_array[available_mask]
                if len(remaining_candidates) > 0:
                    matched_idx = np.random.choice(len(remaining_candidates))
                    matched_gene = remaining_candidates[matched_idx]
                    # Find original index and mark as used
                    orig_idx = np.where(candidate_ids_array == matched_gene)[0][0]
                    available_mask[orig_idx] = False
                else:
                    # Fallback to any valid gene
                    matched_gene = np.random.choice(valid_genes)

        matched_genes.append(matched_gene)

    return sim_idx, np.array(matched_genes)


def worker_set_level_matched(sim_idx, n_genes, candidate_pool_ids, master_table_complete,
                              matching_vars, target_stats, use_ks_test, threshold):
    """Worker function for set-level matched sampling with rejection sampling."""
    import numpy as np
    from scipy import stats

    def calculate_distance(sampled_genes, use_ks=use_ks_test):
        """Calculate how well sampled genes match target distribution."""
        sampled_data = master_table_complete.loc[sampled_genes]

        if use_ks:
            # Use KS test for each variable
            ks_stats = []
            for var in matching_vars:
                ks_stat, _ = stats.ks_2samp(target_stats[var]['values'],
                                            sampled_data[var].values)
                ks_stats.append(ks_stat)
            return np.mean(ks_stats)  # Average KS statistic
        else:
            # Use normalized mean and std differences
            distances = []
            for var in matching_vars:
                sampled_mean = sampled_data[var].mean()
                sampled_std = sampled_data[var].std()

                target_mean = target_stats[var]['mean']
                target_std = target_stats[var]['std']

                # Normalize by target std to make comparable across variables
                if target_std > 0:
                    mean_diff = abs(sampled_mean - target_mean) / target_std
                    std_diff = abs(sampled_std - target_std) / target_std
                    distances.append(mean_diff + std_diff)
                else:
                    distances.append(0)

            return np.mean(distances)

    best_genes = None
    best_distance = float('inf')
    max_attempts_per_sim = 100

    # Try to find acceptable gene set
    for attempt in range(max_attempts_per_sim):
        # Sample random genes from candidate pool
        sampled_genes = np.random.choice(candidate_pool_ids, size=n_genes, replace=False)

        # Calculate distance from target distribution
        distance = calculate_distance(sampled_genes)

        # Track best attempt
        if distance < best_distance:
            best_distance = distance
            best_genes = sampled_genes.copy()

        # Accept if good enough
        if distance <= threshold:
            break

    return sim_idx, best_genes, best_distance


def worker_best_of_n_percentile(sim_idx, n_genes, candidate_pool_ids, candidate_values_pct,
                                  matching_vars, target_stats_pct, n_candidates=100, base_seed=None):
    """Worker function for best-of-N sampling in percentile space.

    For each simulation:
    1. Generate N candidate samples (default 100)
    2. Calculate distance for each candidate
    3. Select the best matching candidate

    This works in percentile space where all variables are on [0,100] scale.
    """
    import numpy as np

    # Set unique random seed for this simulation
    if base_seed is not None:
        np.random.seed(base_seed + sim_idx)

    def calculate_distance_percentile(sampled_indices):
        """Calculate Euclidean distance in percentile space.

        All variables are on [0, 100] percentile scale, so equal weighting
        is appropriate. Distance is computed on mean percentile values only
        (ignoring std) for flexibility.
        """
        sum_sq = 0.0
        for var in matching_vars:
            sampled_vals = candidate_values_pct[var][sampled_indices]
            sampled_mean = sampled_vals.mean()
            target_mean = target_stats_pct[var]['mean']
            sum_sq += (sampled_mean - target_mean) ** 2

        return np.sqrt(sum_sq)

    best_genes = None
    best_distance = float('inf')
    best_indices = None

    # Total number of available candidates
    n_available = len(candidate_pool_ids)

    # Generate N candidate samples and find best
    for attempt in range(n_candidates):
        # Sample random gene indices from candidate pool
        sampled_indices = np.random.choice(n_available, size=n_genes, replace=False)

        # Calculate distance from target distribution
        distance = calculate_distance_percentile(sampled_indices)

        # Track best candidate
        if distance < best_distance:
            best_distance = distance
            best_indices = sampled_indices.copy()

    # Convert best indices to gene IDs
    best_genes = candidate_pool_ids[best_indices]

    return sim_idx, best_genes, best_distance


def worker_stratified_sampling(sim_idx, n_genes, candidate_pool_ids, strata_info,
                                input_gene_strata, base_seed=None):
    """Worker function for stratified sampling.

    Samples genes to match the distribution of input genes across strata (bins).
    Each input gene is assigned to a stratum, and we sample matched genes from
    the same strata in the same proportions.

    This preserves the FULL distribution including tails, unlike propensity
    weighting which only matches the mean.
    """
    import numpy as np

    # Set unique random seed for this simulation
    if base_seed is not None:
        np.random.seed(base_seed + sim_idx)

    sampled_genes = []

    # For each input gene, sample a matched gene from the same stratum
    for input_stratum_key in input_gene_strata:
        # Get candidates in this stratum
        candidates_in_stratum = strata_info[input_stratum_key]['candidates']

        if len(candidates_in_stratum) > 0:
            # Sample one gene from this stratum
            # Use sampling without replacement within this simulation
            available_candidates = [c for c in candidates_in_stratum if c not in sampled_genes]

            if len(available_candidates) > 0:
                sampled_gene = np.random.choice(available_candidates)
            else:
                # Fallback: allow replacement if we've exhausted the stratum
                sampled_gene = np.random.choice(candidates_in_stratum)
        else:
            # Fallback: if stratum is empty, sample from neighboring stratum
            # This is rare but can happen with small datasets
            sampled_gene = np.random.choice(candidate_pool_ids)

        sampled_genes.append(sampled_gene)

    return sim_idx, np.array(sampled_genes)


def worker_propensity_weighted(sim_idx, n_genes, candidate_pool_ids, propensity_probs, base_seed=None):
    """Worker function for propensity score weighted sampling.

    Each candidate gene has a pre-computed probability based on how well it matches
    the target distribution. We simply sample genes according to these probabilities.

    This is much faster than best-of-N since probabilities are computed once upfront.
    """
    import numpy as np

    # Set unique random seed for this simulation
    if base_seed is not None:
        np.random.seed(base_seed + sim_idx)

    # Sample genes according to propensity probabilities
    sampled_genes = np.random.choice(
        candidate_pool_ids,
        size=n_genes,
        p=propensity_probs,
        replace=False  # Each gene can only be selected once per simulation
    )

    return sim_idx, sampled_genes


def worker_sis_matched_percentile(sim_idx, n_genes, candidate_pool_ids, candidate_values_pct,
                                   matching_vars, target_stats_pct, temperature, adaptive_temp):
    """Worker function for Sequential Importance Sampling using PERCENTILES.

    This version uses percentile values which are naturally on the same scale (0-100),
    avoiding distribution and scaling issues.
    """
    import numpy as np

    # Initialize for this simulation
    selected_genes = []
    selected_values_pct = {var: [] for var in matching_vars}
    available_mask = np.ones(len(candidate_pool_ids), dtype=bool)

    # Sequentially select genes
    for i in range(n_genes):
        n_selected = len(selected_genes)
        remaining_slots = n_genes - n_selected

        if remaining_slots == 0:
            break

        # Calculate what we need from remaining genes to hit target
        needed_values_pct = {}
        for var in matching_vars:
            target_mean_pct = target_stats_pct[var]['mean']

            if n_selected > 0:
                current_sum = sum(selected_values_pct[var])
                needed_sum = target_mean_pct * n_genes - current_sum
                needed_mean_pct = needed_sum / remaining_slots
            else:
                needed_mean_pct = target_mean_pct

            needed_values_pct[var] = needed_mean_pct

        # Calculate weights for each available candidate
        # In percentile space, all variables are on [0,100] scale
        all_similarities = []

        for var in matching_vars:
            # Distance in percentile space (0-100 scale)
            distances_pct = np.abs(candidate_values_pct[var] - needed_values_pct[var])

            # Gaussian kernel in percentile space
            # bandwidth=15 percentile points is reasonable (within ~15% of distribution)
            bandwidth_pct = 15.0
            similarities = np.exp(-(distances_pct**2) / (2 * bandwidth_pct**2))

            all_similarities.append(similarities)

        # Combine similarities (sum works well since all are [0,1])
        weights = np.sum(all_similarities, axis=0)

        # Apply temperature for softmax
        if adaptive_temp:
            current_temp = temperature * (1 + n_selected / n_genes)
        else:
            current_temp = temperature

        # Mask out already-selected genes
        weights = weights * available_mask

        # Convert to probabilities with temperature
        if weights.sum() > 0:
            weights = weights ** (1.0 / current_temp)
            probabilities = weights / weights.sum()

            # Sample one gene
            selected_idx = np.random.choice(len(candidate_pool_ids), p=probabilities)
            selected_gene = candidate_pool_ids[selected_idx]

            # Update state
            selected_genes.append(selected_gene)
            for var in matching_vars:
                selected_values_pct[var].append(candidate_values_pct[var][selected_idx])
            available_mask[selected_idx] = False
        else:
            # Fallback: sample from remaining available candidates
            available_indices = np.where(available_mask)[0]
            if len(available_indices) > 0:
                selected_idx = np.random.choice(available_indices)
            else:
                selected_idx = np.random.randint(len(candidate_pool_ids))
            selected_gene = candidate_pool_ids[selected_idx]
            selected_genes.append(selected_gene)
            for var in matching_vars:
                selected_values_pct[var].append(candidate_values_pct[var][selected_idx])
            available_mask[selected_idx] = False

    # Calculate final distance for this simulation
    distance = 0
    for var in matching_vars:
        if len(selected_values_pct[var]) > 0:
            sampled_mean = np.mean(selected_values_pct[var])
            target_mean = target_stats_pct[var]['mean']
            # Distance in percentile units
            distance += abs(sampled_mean - target_mean)

    final_distance = distance / len(matching_vars)

    return sim_idx, np.array(selected_genes), final_distance


def worker_sis_matched(sim_idx, n_genes, candidate_pool_ids, candidate_values,
                       matching_vars, target_stats, temperature, adaptive_temp):
    """Worker function for Sequential Importance Sampling."""
    import numpy as np

    # Initialize for this simulation
    selected_genes = []
    selected_values = {var: [] for var in matching_vars}
    available_mask = np.ones(len(candidate_pool_ids), dtype=bool)

    # Sequentially select genes
    for i in range(n_genes):
        # How many genes have we selected?
        n_selected = len(selected_genes)

        # How many genes left to select?
        remaining_slots = n_genes - n_selected

        if remaining_slots == 0:
            break

        # Calculate what we need from remaining genes to hit target
        needed_values = {}
        for var in matching_vars:
            target_mean = target_stats[var]['mean']

            if n_selected > 0:
                current_sum = sum(selected_values[var])
                needed_sum = target_mean * n_genes - current_sum
                needed_mean = needed_sum / remaining_slots
            else:
                needed_mean = target_mean

            needed_values[var] = needed_mean

        # Calculate weights for each available candidate
        # Use Gaussian kernel which is naturally bounded to [0,1]
        all_similarities = []

        for var in matching_vars:
            # Distance from needed value
            distances = np.abs(candidate_values[var] - needed_values[var])

            # Normalize by std to make comparable across variables
            std = target_stats[var]['std']
            if std > 0:
                normalized_distances = distances / std
            else:
                normalized_distances = distances

            # Gaussian kernel: exp(-d^2 / (2*bandwidth^2))
            # This naturally bounds similarities to [0, 1]
            # Using bandwidth=1.5 std for good balance
            bandwidth = 1.5
            similarities = np.exp(-(normalized_distances**2) / (2 * bandwidth**2))

            all_similarities.append(similarities)

        # Combine similarities (each variable contributes [0,1])
        # Sum is more forgiving than geometric mean
        weights = np.sum(all_similarities, axis=0)

        # Apply temperature for softmax
        if adaptive_temp:
            # Increase temperature as we fill up
            current_temp = temperature * (1 + n_selected / n_genes)
        else:
            current_temp = temperature

        # Mask out already-selected genes
        weights = weights * available_mask

        # Convert to probabilities with temperature
        if weights.sum() > 0:
            # Apply temperature scaling
            weights = weights ** (1.0 / current_temp)
            probabilities = weights / weights.sum()

            # Sample one gene
            selected_idx = np.random.choice(len(candidate_pool_ids), p=probabilities)
            selected_gene = candidate_pool_ids[selected_idx]

            # Update state
            selected_genes.append(selected_gene)
            for var in matching_vars:
                selected_values[var].append(candidate_values[var][selected_idx])
            available_mask[selected_idx] = False
        else:
            # Fallback: sample from remaining available candidates
            available_indices = np.where(available_mask)[0]
            if len(available_indices) > 0:
                selected_idx = np.random.choice(available_indices)
            else:
                selected_idx = np.random.randint(len(candidate_pool_ids))
            selected_gene = candidate_pool_ids[selected_idx]
            selected_genes.append(selected_gene)
            for var in matching_vars:
                selected_values[var].append(candidate_values[var][selected_idx])
            available_mask[selected_idx] = False

    # Calculate final distance for this simulation
    distance = 0
    for var in matching_vars:
        sampled_mean = np.mean(selected_values[var])
        target_mean = target_stats[var]['mean']

        std = target_stats[var]['std']
        if std > 0:
            distance += abs(sampled_mean - target_mean) / std

    final_distance = distance / len(matching_vars)

    return sim_idx, np.array(selected_genes), final_distance


def MatchedGenes(ExpMat, WeightDF, outfile, n_sims=10000,
                 master_table_path=None, percentile_table_path=None,
                 kernel='tricubic', bandwidth=10.0, random_seed=None,
                 matching_vars=None, n_processes=10):
    """
    Generate matched gene sets based on selected variables.
    For each simulation, match each input gene to a similar gene.
    Optimized version: pre-computes distances once.

    Parameters:
    -----------
    matching_vars : list of str, optional
        Variables to match on. Options: 'CDS_length', 'WB', 'LOEUF'
        Default: ['CDS_length', 'WB', 'LOEUF']
    """
    import pandas as pd
    import numpy as np

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    print(f"Input genes: {len(Gene_Weights)}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("MATCHED GENE SAMPLING PARAMETERS")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Kernel type: {kernel}")
    print(f"Bandwidth: {bandwidth}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print("="*60 + "\n")

    # Load or create percentile table
    print(f"Loading percentile table from: {percentile_table_path}")
    master_table_pct = load_or_create_percentile_table(
        master_table_path,
        percentile_table_path,
        matching_vars
    )

    # Filter to genes with complete data
    master_table_complete = master_table_pct[matching_vars + [col + '_pct' for col in matching_vars]].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # Create candidate pool (all valid genes except input genes)
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids)
    print(f"Candidate pool size: {len(candidate_pool_ids)}")

    # PRE-COMPUTE distances and weights for all input genes vs all candidates
    print("Pre-computing distances and weights...")
    input_pct = master_table_complete.loc[input_genes_in_master]
    candidate_pct = master_table_complete.loc[candidate_pool_ids]

    # Calculate distance matrix: (n_input × n_candidates)
    distance_matrix = calculate_pairwise_distances(input_pct, candidate_pct, matching_vars)

    # Apply kernel to get weight matrix: (n_input × n_candidates)
    weight_matrix = apply_kernel_vectorized(distance_matrix, kernel, bandwidth)

    print(f"Distance matrix shape: {distance_matrix.shape}")
    print(f"Weight matrix shape: {weight_matrix.shape}")

    # Create mapping from entrez_id to index in weight matrix
    input_id_to_idx = {gene_id: idx for idx, gene_id in enumerate(input_genes_in_master)}
    candidate_ids_array = candidate_pool_ids.values

    # Run simulations in parallel
    print(f"Running {n_sims} matched simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_matched_genes,
                          [(i, entrez_ids, input_id_to_idx, weight_matrix,
                            candidate_ids_array, valid_genes, random_seed)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    for sim_idx, genes in results:
        sim_matrix[:, sim_idx] = genes

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # Print matching quality summary
    print("\n" + "="*60)
    print("MATCHING QUALITY SUMMARY")
    print("="*60)

    # Calculate some statistics about the matches
    for i, var in enumerate(matching_vars):
        input_vals = input_pct[var].values
        # Get matched genes from first simulation as example
        example_sim = sim_matrix[:len(input_genes_in_master), 0]
        matched_vals = candidate_pct.loc[example_sim][var].values if len(example_sim) > 0 else []

        if len(matched_vals) > 0:
            diff = input_vals - matched_vals
            print(f"{var}:")
            print(f"  Mean absolute difference: {np.abs(diff).mean():.3f}")
            print(f"  Std of difference: {diff.std():.3f}")

    # Show weight statistics
    non_zero_weights = weight_matrix[weight_matrix > 0]
    print(f"\nWeight matrix statistics:")
    print(f"  Total possible matches: {weight_matrix.size}")
    print(f"  Non-zero weights: {len(non_zero_weights)} ({100*len(non_zero_weights)/weight_matrix.size:.2f}%)")
    print(f"  Mean candidates per gene: {(weight_matrix > 0).sum(axis=1).mean():.1f}")
    print("="*60 + "\n")

    # make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} matched simulations to {outfile}")


def SetLevelMatchedGenes_SIS(ExpMat, WeightDF, outfile, n_sims=10000,
                             master_table_path=None, percentile_table_path=None,
                             random_seed=None, matching_vars=None,
                             temperature=1.0, adaptive_temp=True, n_processes=10):
    """
    Generate matched gene sets using SET-LEVEL matching with Sequential Importance Sampling.

    This algorithm builds up the gene set iteratively, at each step selecting genes
    that help guide the overall distribution toward the target.

    Much more efficient than rejection sampling!

    Parameters:
    -----------
    temperature : float
        Temperature parameter for softmax weighting (default 1.0)
        Lower = more greedy (picks best matches), Higher = more random
    adaptive_temp : bool
        If True, increase temperature as we fill the set to maintain diversity
    """
    import pandas as pd
    import numpy as np
    from scipy import stats

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("SET-LEVEL MATCHED GENE SAMPLING (Sequential Importance Sampling)")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Temperature: {temperature}")
    print(f"Adaptive temperature: {adaptive_temp}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print("="*60 + "\n")

    # Load master table
    print(f"Loading master table from: {master_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)

    # Filter to genes with complete data for matching variables
    master_table_complete = master_table[matching_vars].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # Create candidate pool (all valid genes except input genes)
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size: {len(candidate_pool_ids)}")

    # Calculate target distribution statistics from input genes
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats = {}
    for var in matching_vars:
        target_stats[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
        }

    print("\nTarget distribution statistics:")
    for var in matching_vars:
        print(f"  {var}: mean={target_stats[var]['mean']:.3f}, std={target_stats[var]['std']:.3f}")

    # Pre-extract candidate values for efficiency
    candidate_data = master_table_complete.loc[candidate_pool_ids]
    candidate_values = {var: candidate_data[var].values for var in matching_vars}

    # Run simulations in parallel
    print(f"\nRunning {n_sims} SIS simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_sis_matched,
                          [(i, n_genes, candidate_pool_ids, candidate_values,
                            matching_vars, target_stats, temperature, adaptive_temp)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix and collect distances from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    final_distances = []
    for sim_idx, genes, distance in results:
        sim_matrix[:, sim_idx] = genes
        final_distances.append(distance)

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # Print matching quality summary
    print("\n" + "="*60)
    print("SET-LEVEL MATCHING QUALITY SUMMARY (SIS)")
    print("="*60)
    print(f"Mean normalized distance: {np.mean(final_distances):.4f}")
    print(f"Median normalized distance: {np.median(final_distances):.4f}")
    print(f"Max normalized distance: {np.max(final_distances):.4f}")
    print(f"Min normalized distance: {np.min(final_distances):.4f}")

    # Validate on first few simulations
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx}:")
        for var in matching_vars:
            sampled_mean = sampled_data[var].mean()
            sampled_std = sampled_data[var].std()
            target_mean = target_stats[var]['mean']
            target_std = target_stats[var]['std']

            print(f"    {var}: target=({target_mean:.3f}, {target_std:.3f}), "
                  f"sampled=({sampled_mean:.3f}, {sampled_std:.3f})")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} SIS-matched simulations to {outfile}")


def SetLevelMatchedGenes_SIS_Percentile(ExpMat, WeightDF, outfile, n_sims=10000,
                                        master_table_path=None, percentile_table_path=None,
                                        random_seed=None, matching_vars=None,
                                        temperature=1.0, adaptive_temp=True, n_processes=10):
    """
    Generate matched gene sets using SET-LEVEL matching with SIS on PERCENTILES.

    This version uses percentile-transformed values which are naturally on the same
    scale (0-100), avoiding distribution issues from skewed variables.

    Advantages over raw value matching:
    - All variables on same scale (percentiles 0-100)
    - Handles skewed distributions naturally (e.g., CDS_length has skewness=13.6)
    - Distance of 10 percentile points means the same thing for all variables
    - No need to normalize by std

    Parameters:
    -----------
    temperature : float
        Temperature parameter for softmax weighting (default 1.0)
    adaptive_temp : bool
        If True, increase temperature as we fill the set to maintain diversity
    """
    import pandas as pd
    import numpy as np

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("SET-LEVEL MATCHED GENE SAMPLING (SIS with PERCENTILES)")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Temperature: {temperature}")
    print(f"Adaptive temperature: {adaptive_temp}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print(f"Using PERCENTILE matching (bandwidth=15 percentile points)")
    print("="*60 + "\n")

    # Load master table AND percentile table
    print(f"Loading percentile table from: {percentile_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)
    master_table_pct = pd.read_csv(percentile_table_path, index_col=0)

    # Get percentile column names
    pct_cols = [var + '_pct' for var in matching_vars]

    # Filter to genes with complete data for both raw and percentile values
    required_cols = matching_vars + pct_cols
    master_table_complete = master_table_pct[required_cols].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # Create candidate pool (all valid genes except input genes)
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size: {len(candidate_pool_ids)}")

    # Calculate target distribution statistics from input genes (PERCENTILES)
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats_pct = {}
    target_stats_raw = {}

    for var in matching_vars:
        pct_col = var + '_pct'
        target_stats_pct[var] = {
            'mean': input_data[pct_col].mean(),
            'std': input_data[pct_col].std(),
        }
        target_stats_raw[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
        }

    print("\nTarget distribution (PERCENTILES):")
    for var in matching_vars:
        print(f"  {var}_pct: mean={target_stats_pct[var]['mean']:.1f}, std={target_stats_pct[var]['std']:.1f}")

    print("\nTarget distribution (RAW VALUES for reference):")
    for var in matching_vars:
        print(f"  {var}: mean={target_stats_raw[var]['mean']:.3f}, std={target_stats_raw[var]['std']:.3f}")

    # Pre-extract candidate percentile values for efficiency
    candidate_data = master_table_complete.loc[candidate_pool_ids]
    candidate_values_pct = {}
    for var in matching_vars:
        pct_col = var + '_pct'
        candidate_values_pct[var] = candidate_data[pct_col].values

    # Run simulations in parallel
    print(f"\nRunning {n_sims} SIS-PERCENTILE simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_sis_matched_percentile,
                          [(i, n_genes, candidate_pool_ids, candidate_values_pct,
                            matching_vars, target_stats_pct, temperature, adaptive_temp)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix and collect distances from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    final_distances = []
    for sim_idx, genes, distance in results:
        sim_matrix[:, sim_idx] = genes
        final_distances.append(distance)

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # Print matching quality summary
    print("\n" + "="*60)
    print("SET-LEVEL MATCHING QUALITY SUMMARY (SIS-PERCENTILE)")
    print("="*60)
    print(f"Mean distance (percentile units): {np.mean(final_distances):.2f}")
    print(f"Median distance: {np.median(final_distances):.2f}")
    print(f"Max distance: {np.max(final_distances):.2f}")
    print(f"Min distance: {np.min(final_distances):.2f}")

    # Validate on first few simulations - compare BOTH percentiles AND raw values
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx}:")
        print(f"    PERCENTILES:")
        for var in matching_vars:
            pct_col = var + '_pct'
            sampled_mean = sampled_data[pct_col].mean()
            target_mean = target_stats_pct[var]['mean']
            diff = sampled_mean - target_mean
            print(f"      {pct_col}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, diff={diff:+.1f}")

        print(f"    RAW VALUES:")
        for var in matching_vars:
            sampled_mean = sampled_data[var].mean()
            target_mean = target_stats_raw[var]['mean']
            diff = sampled_mean - target_mean
            pct_error = abs(diff) / target_mean * 100
            print(f"      {var}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, error={pct_error:.1f}%")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} SIS-PERCENTILE matched simulations to {outfile}")


def SetLevelMatchedGenes(ExpMat, WeightDF, outfile, n_sims=10000,
                         master_table_path=None, percentile_table_path=None,
                         max_distance=0.15, max_attempts_per_sim=100, random_seed=None,
                         matching_vars=None, use_ks_test=False, ks_threshold=0.3, n_processes=10):
    """
    Generate matched gene sets using SET-LEVEL matching with REJECTION SAMPLING.
    Instead of matching each gene individually, this samples gene sets that have
    similar OVERALL distributions of matching variables.

    This allows much more flexibility in sampling while still controlling for confounders.

    Parameters:
    -----------
    max_distance : float
        Maximum acceptable normalized distance between distributions (default 0.15)
        Distance = mean(|mean_diff|/std_input + |std_diff|/std_input) across variables
    max_attempts_per_sim : int
        Maximum attempts per simulation before accepting best attempt (default 100)
    use_ks_test : bool
        If True, use Kolmogorov-Smirnov test instead of mean/std matching
    ks_threshold : float
        Maximum acceptable KS statistic if use_ks_test=True (default 0.3)
    """
    import pandas as pd
    import numpy as np
    from scipy import stats

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("SET-LEVEL MATCHED GENE SAMPLING PARAMETERS")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Matching criterion: {'KS test' if use_ks_test else 'Mean/Std distance'}")
    if use_ks_test:
        print(f"KS threshold: {ks_threshold}")
    else:
        print(f"Max distance: {max_distance}")
    print(f"Max attempts per simulation: {max_attempts_per_sim}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print("="*60 + "\n")

    # Load master table
    print(f"Loading master table from: {master_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)

    # Filter to genes with complete data for matching variables
    master_table_complete = master_table[matching_vars].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # Create candidate pool (all valid genes except input genes)
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size: {len(candidate_pool_ids)}")

    # Calculate target distribution statistics from input genes
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats = {}
    for var in matching_vars:
        target_stats[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
            'values': input_data[var].values  # For KS test
        }

    print("\nTarget distribution statistics:")
    for var in matching_vars:
        print(f"  {var}: mean={target_stats[var]['mean']:.3f}, std={target_stats[var]['std']:.3f}")

    # Run simulations in parallel
    threshold = ks_threshold if use_ks_test else max_distance
    print(f"\nRunning {n_sims} set-level matched simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_set_level_matched,
                          [(i, n_genes, candidate_pool_ids, master_table_complete,
                            matching_vars, target_stats, use_ks_test, threshold)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix and collect distances from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    accepted_distances = []
    relaxed_count = 0
    for sim_idx, genes, distance in results:
        sim_matrix[:, sim_idx] = genes
        accepted_distances.append(distance)
        if distance > threshold:
            relaxed_count += 1

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # Print matching quality summary
    print("\n" + "="*60)
    print("SET-LEVEL MATCHING QUALITY SUMMARY")
    print("="*60)
    print(f"Simulations where criteria relaxed: {relaxed_count}/{n_sims} ({100*relaxed_count/n_sims:.1f}%)")
    print(f"Distance statistics:")
    print(f"  Mean: {np.mean(accepted_distances):.4f}")
    print(f"  Median: {np.median(accepted_distances):.4f}")
    print(f"  Max: {np.max(accepted_distances):.4f}")

    # Validate on first few simulations
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx}:")
        for var in matching_vars:
            sampled_mean = sampled_data[var].mean()
            sampled_std = sampled_data[var].std()
            target_mean = target_stats[var]['mean']
            target_std = target_stats[var]['std']

            print(f"    {var}: target=({target_mean:.3f}, {target_std:.3f}), "
                  f"sampled=({sampled_mean:.3f}, {sampled_std:.3f})")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} set-level matched simulations to {outfile}")


def SetLevelMatchedGenes_BestOfN_Percentile(ExpMat, WeightDF, outfile, n_sims=10000,
                                             master_table_path=None, percentile_table_path=None,
                                             random_seed=None, matching_vars=None,
                                             n_candidates=100, n_processes=10):
    """
    Generate matched gene sets using best-of-N sampling in PERCENTILE SPACE.

    For each simulation:
    1. Generate N candidate samples (default 100)
    2. Calculate similarity metric for each candidate in percentile space
    3. Select the best matching candidate among the N

    This is simpler than SIS and should give better matching than pure rejection sampling.

    Parameters:
    -----------
    n_candidates : int
        Number of candidate samples to generate per simulation (default: 100)
        Higher values = better matching but slower
    """
    import pandas as pd
    import numpy as np

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("SET-LEVEL MATCHED GENE SAMPLING (Best-of-N with PERCENTILES)")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Candidates per simulation: {n_candidates}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print(f"Parallel processes: {n_processes}")
    print("="*60 + "\n")

    # Load master table AND percentile table
    print(f"Loading percentile table from: {percentile_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)
    master_table_pct = pd.read_csv(percentile_table_path, index_col=0)

    # Get percentile column names
    pct_cols = [var + '_pct' for var in matching_vars]

    # Filter to genes with complete data for both raw and percentile values
    required_cols = matching_vars + pct_cols
    master_table_complete = master_table_pct[required_cols].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # Create candidate pool (all valid genes except input genes)
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size: {len(candidate_pool_ids)}")

    # Calculate target distribution statistics from input genes (PERCENTILES)
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats_pct = {}
    target_stats_raw = {}

    for var in matching_vars:
        pct_col = var + '_pct'
        target_stats_pct[var] = {
            'mean': input_data[pct_col].mean(),
            'std': input_data[pct_col].std(),
        }
        target_stats_raw[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
        }

    print("\nTarget distribution (PERCENTILES):")
    for var in matching_vars:
        print(f"  {var}_pct: mean={target_stats_pct[var]['mean']:.1f}, std={target_stats_pct[var]['std']:.1f}")

    print("\nTarget distribution (RAW VALUES for reference):")
    for var in matching_vars:
        print(f"  {var}: mean={target_stats_raw[var]['mean']:.3f}, std={target_stats_raw[var]['std']:.3f}")

    # Pre-extract candidate percentile values for efficiency
    candidate_data = master_table_complete.loc[candidate_pool_ids]
    candidate_values_pct = {}
    for var in matching_vars:
        pct_col = var + '_pct'
        candidate_values_pct[var] = candidate_data[pct_col].values

    # Run simulations in parallel
    print(f"\nRunning {n_sims} best-of-{n_candidates} simulations using {n_processes} processes...")
    print(f"Total candidate evaluations: {n_sims * n_candidates:,}")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_best_of_n_percentile,
                          [(i, n_genes, candidate_pool_ids, candidate_values_pct,
                            matching_vars, target_stats_pct, n_candidates, random_seed)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix and collect distances from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    final_distances = []
    for sim_idx, genes, distance in results:
        sim_matrix[:, sim_idx] = genes
        final_distances.append(distance)

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # Print matching quality summary
    print("\n" + "="*60)
    print("SET-LEVEL MATCHING QUALITY SUMMARY (Best-of-N Percentile)")
    print("="*60)
    print(f"Distance statistics (lower is better):")
    print(f"  Mean distance: {np.mean(final_distances):.2f}")
    print(f"  Median distance: {np.median(final_distances):.2f}")
    print(f"  Std distance: {np.std(final_distances):.2f}")
    print(f"  Max distance: {np.max(final_distances):.2f}")
    print(f"  Min distance: {np.min(final_distances):.2f}")

    # Validate on first few simulations - compare BOTH percentiles AND raw values
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx} (distance={final_distances[sim_idx]:.2f}):")
        print(f"    PERCENTILES:")
        for var in matching_vars:
            pct_col = var + '_pct'
            sampled_mean = sampled_data[pct_col].mean()
            target_mean = target_stats_pct[var]['mean']
            diff = sampled_mean - target_mean
            print(f"      {pct_col}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, diff={diff:+.1f}")

        print(f"    RAW VALUES:")
        for var in matching_vars:
            sampled_mean = sampled_data[var].mean()
            target_mean = target_stats_raw[var]['mean']
            diff = sampled_mean - target_mean
            pct_error = abs(diff) / target_mean * 100 if target_mean != 0 else 0
            print(f"      {var}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, error={pct_error:.1f}%")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} best-of-{n_candidates} matched simulations to {outfile}")


def PropensityWeightedGenes(ExpMat, WeightDF, outfile, n_sims=10000,
                            master_table_path=None, percentile_table_path=None,
                            random_seed=None, matching_vars=None,
                            bandwidth=15.0, kernel='tricubic', n_processes=10,
                            add_noise=False, noise_scale=5.0, relaxed_matching=False,
                            variable_weights=None):
    """
    Generate matched gene sets using PROPENSITY SCORE WEIGHTING.

    Key idea:
    1. Pre-compute matching probability for ALL candidate genes relative to target set
    2. Sample genes with probability proportional to their matching likelihood
    3. Much faster than best-of-N since probabilities computed once upfront

    This is similar to propensity score matching in causal inference.

    Parameters:
    -----------
    bandwidth : float
        Bandwidth for kernel in percentile space (default: 15.0)
        Smaller = more selective (prefer closer matches)
        Larger = more diverse (broader pool of genes)
    kernel : str
        Kernel function for converting distances to probabilities
        Options:
        - 'gaussian': exp(-(d/h)^2) - Infinite support, smooth decay
        - 'tricubic': (1-(d/h)^3)^3 for d<h, 0 otherwise - Compact support, most diversity
        - 'epanechnikov': (1-(d/h)^2) for d<h, 0 otherwise - Compact support, smooth
        - 'uniform': 1 for d<h, 0 otherwise - Compact support, no weighting within
        - 'linear': (1-d/h) for d<h, 0 otherwise - Compact support, linear decay
    add_noise : bool
        If True, add random noise to target distribution to prevent over-matching
        This preserves variance and statistical power (default: False)
    noise_scale : float
        Scale of noise to add (in percentile units) if add_noise=True (default: 5.0)
        Larger values = more variance preserved = better statistical power
    relaxed_matching : bool
        If True, use relaxed matching that intentionally allows some mismatch
        to preserve statistical power. Adds noise and increases effective bandwidth (default: False)
    variable_weights : dict, optional
        Dictionary mapping variable names to weights for distance calculation.
        Default: {'CDS_length': 1.0, 'WB': 1.0, 'LOEUF': 0.5}
        Lower LOEUF weight helps prevent it from dominating distance calculation.
        If None, uses default weights.
    """
    import pandas as pd
    import numpy as np

    # Set random seed if provided (for reproducibility)
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")
    
    # Set variable weights (default: LOEUF weighted less to prevent domination)
    if variable_weights is None:
        variable_weights = {
            'CDS_length': 1.0,
            'WB': 1.0,
            'LOEUF': 0.5,  # Weight LOEUF less - it often has larger differences and can dominate
            'mean_phastCons': 1.0,
            'n_CDS_bases': 1.0
        }
    
    # Ensure all matching variables have weights
    for var in matching_vars:
        if var not in variable_weights:
            variable_weights[var] = 1.0  # Default weight

    # Print matching parameters
    print("\n" + "="*60)
    print("PROPENSITY SCORE WEIGHTED GENE SAMPLING")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Variable weights: {', '.join([f'{v}={variable_weights.get(v, 1.0):.2f}' for v in matching_vars])}")
    print(f"Kernel function: {kernel}")
    print(f"Kernel bandwidth: {bandwidth} percentile points")
    print(f"Add noise: {add_noise} (scale={noise_scale if add_noise else 'N/A'})")
    print(f"Relaxed matching: {relaxed_matching}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print(f"Parallel processes: {n_processes}")
    print("="*60 + "\n")

    # Load master table AND percentile table
    print(f"Loading percentile table from: {percentile_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)
    master_table_pct = pd.read_csv(percentile_table_path, index_col=0)

    # Get percentile column names
    pct_cols = [var + '_pct' for var in matching_vars]

    # Filter to genes with complete data for both raw and percentile values
    required_cols = matching_vars + pct_cols
    master_table_complete = master_table_pct[required_cols].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # IMPORTANT: Create candidate pool EXCLUDING input genes
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size (before filtering): {len(candidate_pool_ids)}")

    # If relaxed matching, increase effective bandwidth
    effective_bandwidth = bandwidth
    if relaxed_matching:
        effective_bandwidth = bandwidth * 1.5  # Increase bandwidth by 50%
        print(f"\nRelaxed matching enabled: bandwidth increased from {bandwidth} to {effective_bandwidth:.1f}")

    # Calculate target distribution statistics from input genes (PERCENTILES)
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats_pct = {}
    target_stats_raw = {}

    for var in matching_vars:
        pct_col = var + '_pct'
        target_mean = input_data[pct_col].mean()
        target_std = input_data[pct_col].std()
        
        # Apply noise/jitter if requested to prevent over-matching
        if add_noise or relaxed_matching:
            # Add noise to target mean to preserve variance
            # This prevents the null distribution from being too similar to real
            noise = np.random.normal(0, noise_scale, size=1)[0]
            target_mean_adjusted = target_mean + noise
            # Clamp to valid percentile range [0, 100]
            target_mean_adjusted = np.clip(target_mean_adjusted, 0, 100)
        else:
            target_mean_adjusted = target_mean
        
        target_stats_pct[var] = {
            'mean': target_mean_adjusted,
            'mean_original': target_mean,  # Keep original for reporting
            'std': target_std,
        }
        target_stats_raw[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
        }

    print("\nTarget distribution (PERCENTILES):")
    for var in matching_vars:
        mean_val = target_stats_pct[var]['mean']
        mean_orig = target_stats_pct[var].get('mean_original', mean_val)
        if add_noise or relaxed_matching:
            noise_applied = mean_val - mean_orig
            print(f"  {var}_pct: mean={mean_val:.1f} (original={mean_orig:.1f}, noise={noise_applied:+.1f}), std={target_stats_pct[var]['std']:.1f}")
        else:
            print(f"  {var}_pct: mean={mean_val:.1f}, std={target_stats_pct[var]['std']:.1f}")

    print("\nTarget distribution (RAW VALUES for reference):")
    for var in matching_vars:
        print(f"  {var}: mean={target_stats_raw[var]['mean']:.3f}, std={target_stats_raw[var]['std']:.3f}")

    # ====================================================================
    # STEP 1: PRE-COMPUTE PROPENSITY SCORES FOR ALL CANDIDATE GENES
    # ====================================================================
    print("\nComputing propensity scores for all candidate genes...")

    candidate_data = master_table_complete.loc[candidate_pool_ids]
    # Use actual length of candidate_data in case some IDs were filtered out
    n_candidates = len(candidate_data)
    print(f"Candidate pool size (after filtering): {n_candidates}")

    # Compute distance for each candidate gene to target distribution
    distances = np.zeros(n_candidates)

    for var in matching_vars:
        pct_col = var + '_pct'
        candidate_pcts = candidate_data[pct_col].values
        target_mean_pct = target_stats_pct[var]['mean']

        # Squared distance in percentile space
        var_distances = (candidate_pcts - target_mean_pct) ** 2

        # Variable-specific weighting
        # Use provided weights (default: LOEUF weighted less to prevent domination)
        weight = variable_weights.get(var, 1.0)

        distances += weight * var_distances

    # Take square root to get Euclidean distance
    distances = np.sqrt(distances)

    # Convert distances to probabilities using selected kernel
    # Smaller distance = higher probability
    # Use effective_bandwidth (may be adjusted for relaxed matching)
    if kernel == 'gaussian':
        # Gaussian: infinite support, smooth decay
        # exp(-(d/h)^2/2)
        propensity_scores = np.exp(-(distances**2) / (2 * effective_bandwidth**2))

    elif kernel == 'tricubic':
        # Tricubic: compact support (d < bandwidth), smooth at boundaries
        # (1 - (d/h)^3)^3 for d < h, 0 otherwise
        normalized_dist = distances / effective_bandwidth
        propensity_scores = np.zeros_like(distances)
        mask = normalized_dist < 1.0
        propensity_scores[mask] = (1 - normalized_dist[mask]**3)**3

    elif kernel == 'epanechnikov':
        # Epanechnikov: compact support, parabolic shape
        # (1 - (d/h)^2) for d < h, 0 otherwise
        normalized_dist = distances / effective_bandwidth
        propensity_scores = np.zeros_like(distances)
        mask = normalized_dist < 1.0
        propensity_scores[mask] = 1 - normalized_dist[mask]**2

    elif kernel == 'uniform':
        # Uniform: all genes within bandwidth have equal probability
        # 1 for d < h, 0 otherwise
        propensity_scores = (distances < effective_bandwidth).astype(float)

    elif kernel == 'linear':
        # Linear: compact support, linear decay
        # (1 - d/h) for d < h, 0 otherwise
        normalized_dist = distances / effective_bandwidth
        propensity_scores = np.zeros_like(distances)
        mask = normalized_dist < 1.0
        propensity_scores[mask] = 1 - normalized_dist[mask]

    else:
        raise ValueError(f"Unknown kernel: {kernel}. Options: gaussian, tricubic, epanechnikov, uniform, linear")

    # Normalize to sum to 1 (to use as probabilities)
    propensity_probs = propensity_scores / propensity_scores.sum()

    # Print statistics about propensity scores
    print(f"\nPropensity score statistics:")
    print(f"  Min distance: {distances.min():.2f} percentile units")
    print(f"  Mean distance: {distances.mean():.2f}")
    print(f"  Max distance: {distances.max():.2f}")
    print(f"  Min probability: {propensity_probs.min():.2e}")
    print(f"  Max probability: {propensity_probs.max():.2e}")
    print(f"  Effective pool size (genes with >1% of max prob): "
          f"{(propensity_probs > 0.01 * propensity_probs.max()).sum()}")

    # Show top 10 most likely genes
    top_indices = np.argsort(propensity_probs)[-10:][::-1]
    print(f"\nTop 10 most likely genes to be sampled:")

    # Build header dynamically based on matching variables
    header = f"  {'Gene ID':<10} {'Probability':<12} {'Distance':<10}"
    for var in matching_vars:
        header += f" {var+'_pct':<10}"
    print(header)

    # Print top genes with only the matching variables
    for idx in top_indices:
        gene_id = candidate_data.index[idx]  # Use filtered candidate_data.index instead
        prob = propensity_probs[idx]
        dist = distances[idx]

        # Build row dynamically
        row = f"  {gene_id:<10} {prob:<12.6f} {dist:<10.2f}"
        for var in matching_vars:
            pct_col = var + '_pct'
            var_pct = candidate_data.loc[gene_id, pct_col]
            row += f" {var_pct:<10.1f}"
        print(row)

    # ====================================================================
    # STEP 2: RUN SIMULATIONS IN PARALLEL
    # ====================================================================
    print(f"\nRunning {n_sims} propensity-weighted simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    # Use candidate_data.index.values to ensure we use the filtered candidate list
    filtered_candidate_ids = candidate_data.index.values
    results = pool.starmap(worker_propensity_weighted,
                          [(i, n_genes, filtered_candidate_ids, propensity_probs, random_seed)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    for sim_idx, genes in results:
        sim_matrix[:, sim_idx] = genes

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # ====================================================================
    # STEP 3: VALIDATE MATCHING QUALITY
    # ====================================================================
    print("\n" + "="*60)
    print("PROPENSITY WEIGHTED MATCHING QUALITY SUMMARY")
    print("="*60)

    # Validate on first few simulations
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx}:")
        print(f"    PERCENTILES:")
        for var in matching_vars:
            pct_col = var + '_pct'
            sampled_mean = sampled_data[pct_col].mean()
            target_mean = target_stats_pct[var]['mean']
            diff = sampled_mean - target_mean
            print(f"      {pct_col}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, diff={diff:+.1f}")

        print(f"    RAW VALUES:")
        for var in matching_vars:
            sampled_mean = sampled_data[var].mean()
            target_mean = target_stats_raw[var]['mean']
            diff = sampled_mean - target_mean
            pct_error = abs(diff) / target_mean * 100 if target_mean != 0 else 0
            print(f"      {var}: target={target_mean:.1f}, sampled={sampled_mean:.1f}, error={pct_error:.1f}%")

    # Check gene diversity across simulations
    print(f"\n  Gene diversity across all {n_sims} simulations:")
    unique_genes_used = set()
    for sim_idx in range(n_sims):
        unique_genes_used.update(sim_matrix[:, sim_idx])
    print(f"    Unique genes used: {len(unique_genes_used)}/{len(candidate_pool_ids)} "
          f"({100*len(unique_genes_used)/len(candidate_pool_ids):.1f}% of pool)")

    # Most frequently sampled genes
    gene_counts = {}
    for sim_idx in range(n_sims):
        for gene in sim_matrix[:, sim_idx]:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1

    top_frequent = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\n  Top 10 most frequently sampled genes:")
    print(f"    {'Gene ID':<10} {'Count':<10} {'Frequency':<10}")
    for gene_id, count in top_frequent:
        freq = count / n_sims * 100
        print(f"    {gene_id:<10} {count:<10} {freq:<10.1f}%")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} propensity-weighted simulations to {outfile}")


def StratifiedSamplingGenes(ExpMat, WeightDF, outfile, n_sims=10000,
                             master_table_path=None, percentile_table_path=None,
                             random_seed=None, matching_vars=None,
                             n_bins=10, n_processes=10):
    """
    Generate matched gene sets using STRATIFIED SAMPLING.

    Key idea:
    1. Divide genes into strata (bins) based on percentiles for each matching variable
    2. For each input gene, determine its stratum
    3. Sample matched genes from the same strata in the same proportions
    4. This PRESERVES THE FULL DISTRIBUTION including tails!

    Unlike propensity weighting (which matches the mean), stratified sampling
    matches the entire distribution shape.

    Parameters:
    -----------
    n_bins : int
        Number of bins per variable (default: 10 = deciles)
        More bins = finer matching, fewer candidates per bin
        Fewer bins = coarser matching, more candidates per bin
    """
    import pandas as pd
    import numpy as np

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load expression matrix and get valid genes
    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Set up paths for master table
    if master_table_path is None:
        proj_dir = Path(__file__).parent.parent
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    if percentile_table_path is None:
        proj_dir = Path(__file__).parent.parent
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'

    # Set matching variables (use default if not provided)
    if matching_vars is None:
        matching_vars = ['CDS_length', 'WB', 'LOEUF']

    # Validate matching variables
    valid_vars = ['CDS_length', 'WB', 'LOEUF', 'mean_phastCons', 'n_CDS_bases']
    for var in matching_vars:
        if var not in valid_vars:
            raise ValueError(f"Invalid matching variable: {var}. Must be one of {valid_vars}")

    if len(matching_vars) == 0:
        raise ValueError("At least one matching variable must be specified")

    # Print matching parameters
    print("\n" + "="*60)
    print("STRATIFIED SAMPLING (Distribution-Preserving)")
    print("="*60)
    print(f"Matching variables: {', '.join(matching_vars)}")
    print(f"Number of bins per variable: {n_bins}")
    print(f"Total strata: {n_bins ** len(matching_vars)}")
    print(f"Random seed: {random_seed if random_seed is not None else 'None (random)'}")
    print(f"Number of simulations: {n_sims}")
    print(f"Parallel processes: {n_processes}")
    print("="*60 + "\n")

    # Load master table AND percentile table
    print(f"Loading percentile table from: {percentile_table_path}")
    master_table = pd.read_csv(master_table_path, index_col=0)
    master_table_pct = pd.read_csv(percentile_table_path, index_col=0)

    # Get percentile column names
    pct_cols = [var + '_pct' for var in matching_vars]

    # Filter to genes with complete data
    required_cols = matching_vars + pct_cols
    master_table_complete = master_table_pct[required_cols].dropna()

    # Filter to valid genes in expression matrix
    master_table_complete = master_table_complete[master_table_complete.index.isin(valid_genes)]

    # Remove duplicate index entries to prevent dimension mismatch
    master_table_complete = master_table_complete[~master_table_complete.index.duplicated(keep='first')]

    # Find input genes in master table
    input_genes_in_master = [g for g in entrez_ids if g in master_table_complete.index]
    print(f"Input genes in master table: {len(input_genes_in_master)}/{len(entrez_ids)}")

    # IMPORTANT: Create candidate pool EXCLUDING input genes
    candidate_pool_ids = master_table_complete.index.difference(entrez_ids).values
    print(f"Candidate pool size (before filtering): {len(candidate_pool_ids)}")

    # Calculate target distribution statistics (for reporting only)
    input_data = master_table_complete.loc[input_genes_in_master]
    target_stats_pct = {}
    target_stats_raw = {}

    for var in matching_vars:
        pct_col = var + '_pct'
        target_stats_pct[var] = {
            'mean': input_data[pct_col].mean(),
            'std': input_data[pct_col].std(),
            'min': input_data[pct_col].min(),
            'max': input_data[pct_col].max(),
        }
        target_stats_raw[var] = {
            'mean': input_data[var].mean(),
            'std': input_data[var].std(),
            'min': input_data[var].min(),
            'max': input_data[var].max(),
        }

    print("\nTarget distribution (PERCENTILES):")
    for var in matching_vars:
        stats = target_stats_pct[var]
        print(f"  {var}_pct: mean={stats['mean']:.1f}, std={stats['std']:.1f}, "
              f"range=[{stats['min']:.1f}, {stats['max']:.1f}]")

    print("\nTarget distribution (RAW VALUES for reference):")
    for var in matching_vars:
        stats = target_stats_raw[var]
        print(f"  {var}: mean={stats['mean']:.1f}, std={stats['std']:.1f}, "
              f"range=[{stats['min']:.1f}, {stats['max']:.1f}]")

    # ====================================================================
    # STEP 1: CREATE STRATA (BINS) FOR EACH VARIABLE
    # ====================================================================
    print("\nCreating strata...")

    def assign_to_bin(percentile_value, n_bins):
        """Assign a percentile value to a bin (0 to n_bins-1)."""
        bin_idx = int(percentile_value / (100.0 / n_bins))
        # Handle edge case for 100th percentile
        if bin_idx >= n_bins:
            bin_idx = n_bins - 1
        return bin_idx

    # Assign each gene (input and candidate) to strata
    def get_stratum_key(gene_id, data):
        """Get multi-dimensional stratum key for a gene."""
        keys = []
        for var in matching_vars:
            pct_col = var + '_pct'
            pct_val = data.loc[gene_id, pct_col]
            bin_idx = assign_to_bin(pct_val, n_bins)
            keys.append(f"{var}:{bin_idx}")
        return tuple(keys)

    # Assign input genes to strata
    input_gene_strata = []
    for gene_id in input_genes_in_master:
        stratum_key = get_stratum_key(gene_id, master_table_complete)
        input_gene_strata.append(stratum_key)

    # Build strata_info: for each stratum, list of candidate genes
    strata_info = {}
    for gene_id in candidate_pool_ids:
        stratum_key = get_stratum_key(gene_id, master_table_complete)
        if stratum_key not in strata_info:
            strata_info[stratum_key] = {'candidates': []}
        strata_info[stratum_key]['candidates'].append(gene_id)

    # Print statistics about strata
    stratum_counts_input = {}
    for stratum_key in input_gene_strata:
        stratum_counts_input[stratum_key] = stratum_counts_input.get(stratum_key, 0) + 1

    stratum_counts_candidates = {k: len(v['candidates']) for k, v in strata_info.items()}

    print(f"\nStrata statistics:")
    print(f"  Total possible strata: {n_bins ** len(matching_vars)}")
    print(f"  Strata occupied by input genes: {len(stratum_counts_input)}")
    print(f"  Strata occupied by candidates: {len(stratum_counts_candidates)}")

    # Check for empty strata (strata with input genes but no candidates)
    empty_strata = []
    for stratum_key in stratum_counts_input:
        if stratum_key not in strata_info or len(strata_info[stratum_key]['candidates']) == 0:
            empty_strata.append(stratum_key)

    if len(empty_strata) > 0:
        print(f"  WARNING: {len(empty_strata)} strata have input genes but no candidates!")
        print(f"    These will fall back to nearest-neighbor sampling")
    else:
        print(f"  All input strata have candidate genes ✓")

    # Show a few example strata
    print(f"\nExample strata (showing first 5):")
    for i, (stratum_key, count) in enumerate(list(stratum_counts_input.items())[:5]):
        n_candidates = len(strata_info.get(stratum_key, {}).get('candidates', []))
        print(f"  {stratum_key}: {count} input genes, {n_candidates} candidates")

    # ====================================================================
    # STEP 2: RUN SIMULATIONS IN PARALLEL
    # ====================================================================
    print(f"\nRunning {n_sims} stratified simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_stratified_sampling,
                          [(i, n_genes, candidate_pool_ids, strata_info,
                            input_gene_strata, random_seed)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix from results
    sim_matrix = np.empty((len(entrez_ids), n_sims), dtype=object)
    for sim_idx, genes in results:
        sim_matrix[:, sim_idx] = genes

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # ====================================================================
    # STEP 3: VALIDATE MATCHING QUALITY
    # ====================================================================
    print("\n" + "="*60)
    print("STRATIFIED SAMPLING QUALITY SUMMARY")
    print("="*60)

    # Validate on first few simulations
    print(f"\nValidation on first 5 simulations:")
    for sim_idx in range(min(5, n_sims)):
        sampled_genes = sim_matrix[:, sim_idx]
        sampled_data = master_table_complete.loc[sampled_genes]

        print(f"\n  Simulation {sim_idx}:")
        print(f"    PERCENTILES:")
        for var in matching_vars:
            pct_col = var + '_pct'
            sampled_mean = sampled_data[pct_col].mean()
            sampled_min = sampled_data[pct_col].min()
            sampled_max = sampled_data[pct_col].max()
            target_mean = target_stats_pct[var]['mean']
            target_min = target_stats_pct[var]['min']
            target_max = target_stats_pct[var]['max']
            diff_mean = sampled_mean - target_mean
            diff_min = sampled_min - target_min
            diff_max = sampled_max - target_max
            print(f"      {pct_col}:")
            print(f"        mean: target={target_mean:.1f}, sampled={sampled_mean:.1f}, diff={diff_mean:+.1f}")
            print(f"        range: target=[{target_min:.1f},{target_max:.1f}], sampled=[{sampled_min:.1f},{sampled_max:.1f}]")

    # Check gene diversity
    print(f"\n  Gene diversity across all {n_sims} simulations:")
    unique_genes_used = set()
    for sim_idx in range(n_sims):
        unique_genes_used.update(sim_matrix[:, sim_idx])
    print(f"    Unique genes used: {len(unique_genes_used)}/{len(candidate_pool_ids)} "
          f"({100*len(unique_genes_used)/len(candidate_pool_ids):.1f}% of pool)")

    # Most frequently sampled genes
    gene_counts = {}
    for sim_idx in range(n_sims):
        for gene in sim_matrix[:, sim_idx]:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1

    top_frequent = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\n  Top 10 most frequently sampled genes:")
    print(f"    {'Gene ID':<10} {'Count':<10} {'Frequency':<10}")
    for gene_id, count in top_frequent:
        freq = count / n_sims * 100
        print(f"    {gene_id:<10} {count:<10} {freq:<10.1f}%")

    print("="*60 + "\n")

    # Make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} stratified simulations to {outfile}")


def RandomGenes(ExpMat, WeightDF, outfile, GeneProb, n_sims=10000, n_processes=10):
    import pandas as pd
    import numpy as np
    import os

    ExpMat = pd.read_csv(ExpMat, index_col=0)
    valid_genes = ExpMat.index.values

    # Load gene weights and filter to valid genes
    WeightDF = pd.read_csv(WeightDF, header=None)
    ValidWeightDF = WeightDF[WeightDF[0].isin(valid_genes)]
    entrez_ids = ValidWeightDF[0].values
    Gene_Weights = ValidWeightDF[1].values
    n_genes = len(Gene_Weights)
    print(f"Input genes: {n_genes}")

    # Exclude input genes from sampling pool
    valid_genes = np.setdiff1d(valid_genes, entrez_ids)
    print(f"Candidate pool (excluding input genes): {len(valid_genes)}")

    # Prepare probability distribution if provided
    gene_pool = None
    gene_probs = None
    if GeneProb is not None:
        Gene2Prob = pd.read_csv(GeneProb, index_col=0)
        # Only keep genes in Gene2Prob that are also in valid_genes
        filtered_genes = [g for g in Gene2Prob.index.values if g in valid_genes]
        Gene2Prob = Gene2Prob.loc[filtered_genes, :]
        probs = Gene2Prob["Prob"].values
        total = np.sum(probs)
        probs = probs / total
        # Ensure sum to 1
        if len(probs) > 1:
            probs[-1] = 1 - np.sum(probs[:-1])
        Gene2Prob["Prob"] = probs

        gene_pool = Gene2Prob.index.values
        gene_probs = Gene2Prob["Prob"].values

    # Run simulations in parallel
    print(f"Running {n_sims} random simulations using {n_processes} processes...")
    pool = multiprocessing.Pool(processes=n_processes)
    results = pool.starmap(worker_random_genes,
                          [(i, valid_genes, n_genes, gene_probs, gene_pool)
                           for i in range(n_sims)])
    pool.close()
    pool.join()

    # Build simulation matrix from results
    sim_matrix = np.empty((n_genes, n_sims), dtype=object)
    for sim_idx, genes in results:
        sim_matrix[:, sim_idx] = genes

    # Build output DataFrame
    out_df = pd.DataFrame(sim_matrix, index=entrez_ids, columns=[str(i) for i in range(n_sims)])
    out_df.insert(0, "GeneWeight", Gene_Weights)

    # make dir if not exist
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(outfile)
    print(f"Saved all {n_sims} simulations to {outfile}")


###########################################################################
## Args and Main Functions
###########################################################################
def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', type=str, help='Output file')
    #parser.add_argument('-o', '--outdir', type=str, help='Output directory')
    parser.add_argument('-w', '--WeightDF', type=str, help='Weight DF for control geneset')
    parser.add_argument('-p', '--GeneProb', default=None, help='GeneProb Filname or None if dont use')
    parser.add_argument('--n_sims', type=int, default=10000, help='Number of ctrl simulations')
    parser.add_argument('--GW_Dir', type=str, help="dirctory of ctrl gene weights")
    parser.add_argument('--SpecMat', type=str, help="Filename of bias matrix")
    parser.add_argument('--sampling_mode', type=str, default='random',
                       choices=['random', 'matched', 'set_level_matched'],
                       help='Sampling mode: random (default), matched (gene-by-gene matching), or set_level_matched (distribution matching)')
    parser.add_argument('--master_table', type=str, default=None,
                       help='Path to master table for matched mode (default: dat/Variable_2_Match_master_table.csv)')
    parser.add_argument('--percentile_table', type=str, default=None,
                       help='Path to percentile table for matched mode (default: dat/Variable_2_Match_master_table_pct.csv)')
    parser.add_argument('--kernel', type=str, default='tricubic',
                       choices=['uniform', 'tricubic'],
                       help='Kernel type for gene-by-gene matched mode (default: tricubic)')
    parser.add_argument('--bandwidth', type=float, default=10.0,
                       help='Bandwidth for kernel in gene-by-gene matched mode (default: 10.0)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for reproducibility')
    parser.add_argument('--matching_variables', type=str, default=None,
                       help='Comma-separated list of variables to match on. Available: CDS_length, WB, LOEUF, mean_phastCons, n_CDS_bases (e.g., "CDS_length,WB,LOEUF")')

    # Parameters for set-level matching (rejection sampling)
    parser.add_argument('--max_distance', type=float, default=0.15,
                       help='Max normalized distance for set-level matching (default: 0.15)')
    parser.add_argument('--max_attempts', type=int, default=100,
                       help='Max attempts per simulation for set-level matching (default: 100)')
    parser.add_argument('--use_ks_test', action='store_true',
                       help='Use KS test instead of mean/std for set-level matching')
    parser.add_argument('--ks_threshold', type=float, default=0.3,
                       help='KS statistic threshold if using KS test (default: 0.3)')

    # Parameters for Sequential Importance Sampling (SIS)
    parser.add_argument('--use_sis', action='store_true',
                       help='Use Sequential Importance Sampling for set-level matching (more efficient)')
    parser.add_argument('--use_percentile', action='store_true',
                       help='Use percentile-based matching (recommended - handles skewed distributions better)')
    parser.add_argument('--temperature', type=float, default=1.0,
                       help='Temperature for SIS sampling (default: 1.0, lower=greedy, higher=random)')
    parser.add_argument('--adaptive_temp', action='store_true', default=True,
                       help='Use adaptive temperature in SIS (default: True)')

    # Parameters for Best-of-N sampling
    parser.add_argument('--use_best_of_n', action='store_true',
                       help='Use Best-of-N sampling for set-level matching (simple and effective)')
    parser.add_argument('--n_candidates', type=int, default=100,
                       help='Number of candidates to evaluate per simulation in Best-of-N mode (default: 100)')

    # Parameters for Propensity Score Weighted sampling
    parser.add_argument('--use_propensity', action='store_true',
                       help='Use propensity score weighted sampling (fast and effective)')
    parser.add_argument('--propensity_bandwidth', type=float, default=15.0,
                       help='Bandwidth for kernel in propensity weighting (default: 15.0 percentile units)')
    parser.add_argument('--propensity_kernel', type=str, default='tricubic',
                       choices=['gaussian', 'tricubic', 'epanechnikov', 'uniform', 'linear'],
                       help='Kernel function for propensity weighting (default: tricubic)')
    parser.add_argument('--add_noise', action='store_true',
                       help='Add noise to target distribution to prevent over-matching (preserves statistical power)')
    parser.add_argument('--noise_scale', type=float, default=5.0,
                       help='Scale of noise to add in percentile units if --add_noise is used (default: 5.0)')
    parser.add_argument('--relaxed_matching', action='store_true',
                       help='Use relaxed matching mode: adds noise and increases bandwidth to preserve variance')
    parser.add_argument('--loeuf_weight', type=float, default=0.5,
                       help='Weight for LOEUF in distance calculation (default: 0.5). Lower values reduce LOEUF influence.')

    # Parallelization parameters
    parser.add_argument('--n_processes', type=int, default=10,
                       help='Number of parallel processes to use (default: 10)')

    args = parser.parse_args()
    return args

def main():
    args = GetOptions()
    SpecMat = args.SpecMat
    WeightDF = args.WeightDF
    outfile = args.outfile
    n_sims = args.n_sims
    sampling_mode = args.sampling_mode

    print("\n" + "="*70)
    print(f"GENE WEIGHT GENERATION - MODE: {sampling_mode.upper()}")
    print("="*70)

    if sampling_mode == 'matched':
        # Parse matching variables from comma-separated string
        matching_vars = None
        if args.matching_variables:
            matching_vars = [v.strip() for v in args.matching_variables.split(',')]

        # Use gene-by-gene matched sampling
        MatchedGenes(
            SpecMat,
            WeightDF,
            outfile,
            n_sims=n_sims,
            master_table_path=args.master_table,
            percentile_table_path=args.percentile_table,
            kernel=args.kernel,
            bandwidth=args.bandwidth,
            random_seed=args.seed,
            matching_vars=matching_vars,
            n_processes=args.n_processes
        )
    elif sampling_mode == 'set_level_matched':
        # Parse matching variables from comma-separated string
        matching_vars = None
        if args.matching_variables:
            matching_vars = [v.strip() for v in args.matching_variables.split(',')]

        # Check if using Propensity Score Weighted sampling (NEW - FASTEST)
        if args.use_propensity:
            print("Using Propensity Score Weighted sampling - Fast and Effective!")
            # Build variable weights dictionary
            var_weights = {
                'CDS_length': 1.0,
                'WB': 1.0,
                'LOEUF': args.loeuf_weight  # Default 0.5 to reduce LOEUF influence
            }
            PropensityWeightedGenes(
                SpecMat,
                WeightDF,
                outfile,
                n_sims=n_sims,
                master_table_path=args.master_table,
                percentile_table_path=args.percentile_table,
                random_seed=args.seed,
                matching_vars=matching_vars,
                bandwidth=args.propensity_bandwidth,
                kernel=args.propensity_kernel,
                n_processes=args.n_processes,
                add_noise=args.add_noise,
                noise_scale=args.noise_scale,
                relaxed_matching=args.relaxed_matching,
                variable_weights=var_weights
            )
        # Check if using Best-of-N sampling (RECOMMENDED)
        elif args.use_best_of_n:
            print("Using Best-of-N sampling with PERCENTILES - Simple and Effective!")
            SetLevelMatchedGenes_BestOfN_Percentile(
                SpecMat,
                WeightDF,
                outfile,
                n_sims=n_sims,
                master_table_path=args.master_table,
                percentile_table_path=args.percentile_table,
                random_seed=args.seed,
                matching_vars=matching_vars,
                n_candidates=args.n_candidates,
                n_processes=args.n_processes
            )
        # Check if using Sequential Importance Sampling
        elif args.use_sis:
            # Check if using percentile-based matching
            if args.use_percentile:
                print("Using Sequential Importance Sampling (SIS) with PERCENTILES - RECOMMENDED!")
                SetLevelMatchedGenes_SIS_Percentile(
                    SpecMat,
                    WeightDF,
                    outfile,
                    n_sims=n_sims,
                    master_table_path=args.master_table,
                    percentile_table_path=args.percentile_table,
                    random_seed=args.seed,
                    matching_vars=matching_vars,
                    temperature=args.temperature,
                    adaptive_temp=args.adaptive_temp,
                    n_processes=args.n_processes
                )
            else:
                print("Using Sequential Importance Sampling (SIS) with raw values")
                SetLevelMatchedGenes_SIS(
                    SpecMat,
                    WeightDF,
                    outfile,
                    n_sims=n_sims,
                    master_table_path=args.master_table,
                    percentile_table_path=args.percentile_table,
                    random_seed=args.seed,
                    matching_vars=matching_vars,
                    temperature=args.temperature,
                    adaptive_temp=args.adaptive_temp,
                    n_processes=args.n_processes
                )
        else:
            print("Using Rejection Sampling")
            # Use set-level matched sampling (distribution matching with rejection sampling)
            SetLevelMatchedGenes(
                SpecMat,
                WeightDF,
                outfile,
                n_sims=n_sims,
                master_table_path=args.master_table,
                percentile_table_path=args.percentile_table,
                max_distance=args.max_distance,
                max_attempts_per_sim=args.max_attempts,
                random_seed=args.seed,
                matching_vars=matching_vars,
                use_ks_test=args.use_ks_test,
                ks_threshold=args.ks_threshold,
                n_processes=args.n_processes
            )
    else:
        # Use random gene sampling
        GeneProb = args.GeneProb
        RandomGenes(SpecMat, WeightDF, outfile, GeneProb, n_sims, n_processes=args.n_processes)

    return

if __name__ == '__main__':
    main()
