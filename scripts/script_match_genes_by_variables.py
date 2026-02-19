#!/usr/bin/env python3
"""
Match genes based on CDS_length, WB, and LOEUF using percentile-based kernel weighting.

This script takes an input gene list and produces a matched set of genes based on
percentile distances in CDS_length, WB (BrainSpan expression), and LOEUF space.
Uses uniform or tricubic kernel weighting for matching.
"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

import sys
import os

# Derive ProjDIR from script location
script_dir = Path(__file__).parent
proj_dir = script_dir.parent
ProjDIR = str(proj_dir) + "/"
sys.path.insert(1, str(proj_dir / 'src'))
from CellType_PSY import *


def convert_to_percentiles(df, columns):
    """Convert specified columns to percentiles (0-100 scale)."""
    df_percentiles = df.copy()
    for col in columns:
        if col in df.columns:
            # Calculate percentiles (0-100 scale)
            df_percentiles[col + '_pct'] = df[col].rank(pct=True) * 100
    return df_percentiles


def uniform_kernel(distance, bandwidth):
    """Uniform kernel: weight = 1 if distance <= bandwidth, else 0."""
    return (distance <= bandwidth).astype(float)


def tricubic_kernel(distance, bandwidth):
    """Tricubic kernel: weight = (1 - (d/h)^3)^3 if d < h, else 0."""
    normalized_dist = distance / bandwidth
    mask = normalized_dist < 1.0
    weights = np.zeros_like(normalized_dist)
    weights[mask] = (1 - normalized_dist[mask]**3)**3
    return weights


def calculate_distance_percentile(gene_percentiles, candidate_percentiles):
    """Calculate Euclidean distance in percentile space."""
    # Get percentile columns (gene_percentiles is a Series, candidate_percentiles is a DataFrame)
    pct_cols = [col for col in candidate_percentiles.columns if col.endswith('_pct')]
    
    # Extract values - gene_percentiles is a Series, candidate_percentiles is a DataFrame
    gene_pct = gene_percentiles[pct_cols].values  # 1D array
    candidate_pct = candidate_percentiles[pct_cols].values  # 2D array (n_genes x n_vars)
    
    # Calculate Euclidean distance: sqrt(sum((candidate - gene)^2))
    # Broadcasting: candidate_pct (n_genes, n_vars) - gene_pct (n_vars,) -> (n_genes, n_vars)
    distances = np.sqrt(np.sum((candidate_pct - gene_pct)**2, axis=1))
    return distances


def match_gene(gene_id, candidate_table_pct, master_table_pct, 
               kernel='tricubic', bandwidth=10.0, min_candidates=10, random_seed=None):
    """
    Match a gene to candidates based on percentile distance and kernel weighting.
    
    Parameters:
    -----------
    gene_id : int
        Entrez ID of gene to match
    candidate_table_pct : DataFrame
        Candidate genes with percentile values
    master_table_pct : DataFrame
        All genes with percentile values (for getting target gene values)
    kernel : str
        'uniform' or 'tricubic'
    bandwidth : float
        Bandwidth for kernel (in percentile units)
    min_candidates : int
        Minimum number of candidates required
    random_seed : int, optional
        Random seed for reproducibility
    
    Returns:
    --------
    matched_gene_id : int or None
        Matched gene ID, or None if no suitable match found
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    if gene_id not in master_table_pct.index:
        return None
    
    # Get target gene percentiles
    target_percentiles = master_table_pct.loc[gene_id]
    
    # Calculate distances to all candidates
    distances = calculate_distance_percentile(
        target_percentiles, 
        candidate_table_pct
    )
    
    # Apply kernel
    if kernel == 'uniform':
        weights = uniform_kernel(distances, bandwidth)
    elif kernel == 'tricubic':
        weights = tricubic_kernel(distances, bandwidth)
    else:
        raise ValueError(f"Unknown kernel: {kernel}. Use 'uniform' or 'tricubic'.")
    
    # Get candidate gene IDs
    candidate_ids = candidate_table_pct.index.values
    
    # Filter to candidates with non-zero weight
    valid_mask = weights > 0
    valid_candidates = candidate_ids[valid_mask]
    valid_weights = weights[valid_mask]
    
    if len(valid_candidates) < min_candidates:
        # If not enough candidates, increase bandwidth and try again
        if kernel == 'tricubic':
            weights = tricubic_kernel(distances, bandwidth * 2)
            valid_mask = weights > 0
            valid_candidates = candidate_ids[valid_mask]
            valid_weights = weights[valid_mask]
        
        if len(valid_candidates) < min_candidates:
            return None
    
    # Normalize weights to probabilities
    if valid_weights.sum() > 0:
        probabilities = valid_weights / valid_weights.sum()
        # Sample one gene based on weights
        matched_idx = np.random.choice(len(valid_candidates), p=probabilities)
        return valid_candidates[matched_idx]
    else:
        return None


def load_master_table(master_table_path):
    """Load the master table with matching variables."""
    if not os.path.exists(master_table_path):
        raise FileNotFoundError(f"Master table not found: {master_table_path}")
    
    master_table = pd.read_csv(master_table_path, index_col=0)
    master_table.index = master_table.index.astype(int)
    return master_table


def load_or_create_percentile_table(master_table_path, percentile_table_path, matching_vars):
    """Load or create percentile version of master table."""
    if os.path.exists(percentile_table_path):
        print(f"Loading percentile table from: {percentile_table_path}")
        master_table_pct = pd.read_csv(percentile_table_path, index_col=0)
        master_table_pct.index = master_table_pct.index.astype(int)
        return master_table_pct
    
    print(f"Creating percentile table and saving to: {percentile_table_path}")
    master_table = load_master_table(master_table_path)
    
    # Filter to genes with complete data
    master_table_complete = master_table[matching_vars].dropna()
    
    # Convert to percentiles
    master_table_pct = convert_to_percentiles(master_table_complete, matching_vars)
    
    # Save percentile table
    master_table_pct.to_csv(percentile_table_path)
    print(f"Percentile table saved: {len(master_table_pct)} genes")
    
    return master_table_pct


def load_gene_list(gene_list_path):
    """Load gene list from file (one gene per line, entrez IDs)."""
    GW = Fil2Dict(gene_list_path)
    genes = list(GW.keys())
    return genes


def match_gene_set(input_genes, master_table_pct, matching_vars, 
                  kernel='tricubic', bandwidth=10.0, exclude_genes=None,
                  random_seed=None):
    """
    Match a set of genes to candidates.
    
    Parameters:
    -----------
    input_genes : list of int
        List of entrez IDs to match
    master_table_pct : DataFrame
        Master table with percentile values
    matching_vars : list of str
        Variable names to match on
    kernel : str
        'uniform' or 'tricubic'
    bandwidth : float
        Bandwidth for kernel
    exclude_genes : set, optional
        Genes to exclude from candidate pool (e.g., input genes)
    random_seed : int, optional
        Random seed for reproducibility
    
    Returns:
    --------
    dict
        Dictionary mapping input gene IDs to matched gene IDs
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Filter master table to complete data
    master_table_complete = master_table_pct[matching_vars + [col + '_pct' for col in matching_vars]].dropna()
    
    # Find input genes in master table
    input_genes_in_master = [g for g in input_genes if g in master_table_complete.index]
    
    if exclude_genes is None:
        exclude_genes = set()
    
    # Create candidate pool (exclude input genes and specified exclusions)
    exclude_set = set(input_genes_in_master) | set(exclude_genes)
    candidate_genes = master_table_complete.index.difference(exclude_set).tolist()
    candidate_table_pct = master_table_complete.loc[candidate_genes]
    
    print(f"\nMatching {len(input_genes_in_master)}/{len(input_genes)} genes")
    print(f"Candidate pool size: {len(candidate_genes)}")
    
    # Match each gene
    matched_genes = {}
    candidate_pool = candidate_genes.copy()
    
    for gene_id in input_genes_in_master:
        if gene_id not in master_table_pct.index:
            continue
        
        # Update candidate table from current pool
        candidate_table_pct = master_table_complete.loc[candidate_pool]
        
        matched = match_gene(gene_id, candidate_table_pct, master_table_pct,
                            kernel=kernel, bandwidth=bandwidth, random_seed=None)
        
        if matched:
            matched_genes[gene_id] = matched
            # Remove matched gene from candidate pool to avoid duplicates
            if matched in candidate_pool:
                candidate_pool.remove(matched)
        else:
            print(f"Warning: No match found for gene {gene_id}")
    
    return matched_genes


def GetOptions():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Match genes based on CDS_length, WB, and LOEUF using percentile-based kernel weighting.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, type=str,
                       help='Input gene list file (one entrez ID per line)')
    parser.add_argument('-o', '--output', required=True, type=str,
                       help='Output directory for matched genes')
    parser.add_argument('-m', '--master_table', type=str,
                       default=None,
                       help='Path to master table CSV (default: ../dat/Variable_2_Match_master_table.csv)')
    parser.add_argument('-p', '--percentile_table', type=str,
                       default=None,
                       help='Path to percentile table CSV (default: ../dat/Variable_2_Match_master_table_pct.csv)')
    parser.add_argument('--kernel', type=str, default='tricubic',
                       choices=['uniform', 'tricubic'],
                       help='Kernel type for weighting (default: tricubic)')
    parser.add_argument('--bandwidth', type=float, default=10.0,
                       help='Bandwidth for kernel in percentile units (default: 10.0)')
    parser.add_argument('--min_candidates', type=int, default=10,
                       help='Minimum number of candidates required (default: 10)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for reproducibility')
    parser.add_argument('--exclude', type=str, default=None,
                       help='File with genes to exclude from candidate pool (one per line)')
    
    args = parser.parse_args()
    return args


def main():
    """Main function."""
    args = GetOptions()
    
    # Set up paths
    script_dir = Path(__file__).parent
    proj_dir = script_dir.parent
    
    if args.master_table is None:
        master_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table.csv'
    else:
        master_table_path = Path(args.master_table)
    
    if args.percentile_table is None:
        percentile_table_path = proj_dir / 'dat' / 'Variable_2_Match_master_table_pct.csv'
    else:
        percentile_table_path = Path(args.percentile_table)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Matching variables
    matching_vars = ['CDS_length', 'WB', 'LOEUF']
    
    # Load or create percentile table
    print("Loading percentile table...")
    master_table_pct = load_or_create_percentile_table(
        master_table_path, 
        percentile_table_path, 
        matching_vars
    )
    
    # Load input gene list
    print(f"Loading input gene list from: {args.input}")
    input_genes = load_gene_list(args.input)
    print(f"Loaded {len(input_genes)} genes")
    
    # Load exclude genes if provided
    exclude_genes = None
    if args.exclude:
        exclude_genes = load_gene_list(args.exclude)
        print(f"Excluding {len(exclude_genes)} genes from candidate pool")
    
    # Match genes
    print(f"\nMatching parameters:")
    print(f"  Kernel: {args.kernel}")
    print(f"  Bandwidth: {args.bandwidth}")
    print(f"  Min candidates: {args.min_candidates}")
    if args.seed is not None:
        print(f"  Random seed: {args.seed}")
    
    matched_genes = match_gene_set(
        input_genes,
        master_table_pct,
        matching_vars,
        kernel=args.kernel,
        bandwidth=args.bandwidth,
        exclude_genes=exclude_genes,
        random_seed=args.seed
    )
    
    print(f"\nMatching complete: {len(matched_genes)}/{len(input_genes)} genes matched")
    
    # Create output DataFrame
    if matched_genes:
        # Load master table for comparison
        master_table = load_master_table(master_table_path)
        
        matching_df = pd.DataFrame({
            'Input_Gene': list(matched_genes.keys()),
            'Matched_Gene': list(matched_genes.values())
        })
        
        # Add variable values for comparison
        for var in matching_vars:
            matching_df[f'Input_{var}'] = matching_df['Input_Gene'].map(master_table[var])
            matching_df[f'Matched_{var}'] = matching_df['Matched_Gene'].map(master_table[var])
            matching_df[f'{var}_diff'] = matching_df[f'Input_{var}'] - matching_df[f'Matched_{var}']
        
        # Save matching table
        matching_output = output_dir / 'matched_genes.csv'
        matching_df.to_csv(matching_output, index=False)
        print(f"Saved matching table to: {matching_output}")
        
        # Save matched gene list
        matched_list_output = output_dir / 'matched_genes_list.txt'
        with open(matched_list_output, 'w') as f:
            for gene in matched_genes.values():
                f.write(f"{gene}\n")
        print(f"Saved matched gene list to: {matched_list_output}")
        
        # Print summary
        print("\nMatching summary:")
        print(matching_df[[col for col in matching_df.columns if 'diff' in col]].describe())
    else:
        print("No genes were matched!")


if __name__ == '__main__':
    main()

