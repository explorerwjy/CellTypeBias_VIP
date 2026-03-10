#!/usr/bin/env python3
"""
Step 1: Preprocess human reference from Siletti et al. 2023.

Memory-safe approach: loads one h5ad at a time, subsamples immediately,
subsets to shared gene set, then concatenates small chunks.

Output: results/cross_species/preprocessed/human_reference.h5ad
"""

import sys
import gc
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import psutil
import yaml
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_MATCHING_DIR = SCRIPT_DIR.parent.parent


def load_config():
    config_path = ATLAS_MATCHING_DIR / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def print_memory(label=""):
    rss_gb = psutil.Process().memory_info().rss / 1e9
    print(f"  [MEM] {label}: {rss_gb:.1f} GB RSS")


def load_ortholog_gene_set(ortho_dir):
    """Load the set of human symbols that have mouse orthologs."""
    ortho_path = Path(ortho_dir) / "ortholog_mapping.csv"
    df = pd.read_csv(ortho_path)
    human_symbols = set(df["human_symbol"].dropna().unique())
    print(f"  Ortholog gene set: {len(human_symbols)} human symbols with mouse orthologs")
    return human_symbols


def process_one_supercluster(h5ad_dir, file_info, subsample_n, gene_set, species_label,
                              cell_class_label, seed=42):
    """Load one Siletti h5ad, subsample, convert genes, subset to gene_set."""
    h5ad_path = Path(h5ad_dir) / file_info["file"]
    total_cells = file_info.get("total_cells", None)

    print(f"\n  Loading {h5ad_path.name}...")
    print_memory("before load")

    # Load in backed mode to avoid full memory allocation
    adata_backed = sc.read_h5ad(h5ad_path, backed="r")
    n_cells = adata_backed.n_obs
    print(f"    Cells: {n_cells:,}, Genes: {adata_backed.n_vars}")

    # Subsample indices
    np.random.seed(seed)
    n_take = min(subsample_n, n_cells)
    if n_take < n_cells:
        indices = np.sort(np.random.choice(n_cells, n_take, replace=False))
    else:
        indices = np.arange(n_cells)

    # Build ensembl->symbol mapping from .var
    var_df = adata_backed.var.copy()
    gene_symbol_col = "Gene"
    if gene_symbol_col not in var_df.columns:
        # Try alternative column names
        for col in ["gene_symbol", "gene_name", "Symbol"]:
            if col in var_df.columns:
                gene_symbol_col = col
                break

    # Map ensembl -> symbol, find which columns to keep
    ensembl_to_symbol = {}
    for ensembl_id, row in var_df.iterrows():
        sym = row.get(gene_symbol_col, None)
        if sym and pd.notna(sym):
            ensembl_to_symbol[ensembl_id] = str(sym).strip()

    # Find gene indices that map to our gene_set
    keep_gene_indices = []
    keep_gene_symbols = []
    keep_gene_ensembl = []
    seen_symbols = set()
    for i, ensembl_id in enumerate(var_df.index):
        sym = ensembl_to_symbol.get(ensembl_id, None)
        if sym and sym in gene_set and sym not in seen_symbols:
            keep_gene_indices.append(i)
            keep_gene_symbols.append(sym)
            keep_gene_ensembl.append(ensembl_id)
            seen_symbols.add(sym)

    keep_gene_indices = np.array(keep_gene_indices)
    print(f"    Keeping {len(keep_gene_indices)} genes (of {len(gene_set)} in gene set)")

    # Extract subsampled + gene-subsetted matrix
    print(f"    Extracting {n_take:,} cells x {len(keep_gene_indices)} genes...")
    X_sub = adata_backed.X[indices][:, keep_gene_indices]
    if not sparse.issparse(X_sub):
        X_sub = sparse.csr_matrix(X_sub)
    else:
        X_sub = X_sub.copy()

    # Extract obs
    obs_sub = adata_backed.obs.iloc[indices].copy()
    adata_backed.file.close()
    del adata_backed
    gc.collect()

    # Build new AnnData with human gene symbols as var_names
    var_new = pd.DataFrame(
        {"ensembl_id": keep_gene_ensembl},
        index=keep_gene_symbols,
    )

    adata_sub = ad.AnnData(X=X_sub, obs=obs_sub, var=var_new)
    adata_sub.layers["counts"] = adata_sub.X.copy()

    # Add metadata columns
    adata_sub.obs["species"] = "human"
    adata_sub.obs["cell_class"] = cell_class_label
    # sample_species = sample_id + "_human"
    if "sample_id" in adata_sub.obs.columns:
        adata_sub.obs["sample_species"] = (
            adata_sub.obs["sample_id"].astype(str) + "_human"
        )
    else:
        adata_sub.obs["sample_species"] = "unknown_human"

    # Keep only essential obs columns
    essential_cols = [
        "sample_id", "supercluster_term", "cluster_id", "subcluster_id",
        "species", "cell_class", "sample_species",
    ]
    cols_to_keep = [c for c in essential_cols if c in adata_sub.obs.columns]
    adata_sub.obs = adata_sub.obs[cols_to_keep].copy()

    print(f"    Result: {adata_sub.shape}")
    print_memory("after processing")

    return adata_sub


def main():
    parser = argparse.ArgumentParser(description="Preprocess human reference")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    config = load_config()
    output_dir = args.output_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["preprocessed"])
    output_dir.mkdir(parents=True, exist_ok=True)
    ortho_dir = ATLAS_MATCHING_DIR / config["output"]["subdirs"]["orthologs"]

    print("=" * 70)
    print("STEP 1: PREPROCESS HUMAN REFERENCE (SILETTI ET AL. 2023)")
    print("=" * 70)
    print_memory("start")

    # Load ortholog gene set to know which genes to keep
    print("\n[GENES] Loading ortholog gene set...")
    ortho_genes = load_ortholog_gene_set(ortho_dir)

    # Process each supercluster file sequentially
    h5ad_dir = config["human_reference"]["h5ad_dir"]
    chunks = []

    # Interneuron files
    print("\n[INTERNEURONS] Processing interneuron superclusters...")
    for name, info in config["human_reference"]["interneuron_files"].items():
        adata_chunk = process_one_supercluster(
            h5ad_dir, info, info["subsample_n"],
            gene_set=ortho_genes,
            species_label="human",
            cell_class_label="GABAergic",
            seed=args.seed,
        )
        chunks.append(adata_chunk)
        del adata_chunk
        gc.collect()

    # Control files
    print("\n[CONTROLS] Processing control superclusters...")
    for name, info in config["human_reference"]["control_files"].items():
        if "Glut" in name:
            cell_class = "Glutamatergic"
        else:
            cell_class = "NonNeuronal"

        adata_chunk = process_one_supercluster(
            h5ad_dir, info, info["subsample_n"],
            gene_set=ortho_genes,
            species_label="human",
            cell_class_label=cell_class,
            seed=args.seed,
        )
        chunks.append(adata_chunk)
        del adata_chunk
        gc.collect()

    # Concatenate all chunks
    print(f"\n[CONCAT] Concatenating {len(chunks)} chunks...")
    print_memory("before concat")

    # Ensure all chunks have the same gene set (intersect)
    common_genes = set(chunks[0].var_names)
    for chunk in chunks[1:]:
        common_genes &= set(chunk.var_names)
    common_genes = sorted(common_genes)
    print(f"  Common genes across all chunks: {len(common_genes)}")

    # Subset each chunk to common genes
    for i in range(len(chunks)):
        chunks[i] = chunks[i][:, common_genes].copy()

    adata_human = ad.concat(chunks, join="inner")
    del chunks
    gc.collect()

    # Ensure counts layer
    if "counts" not in adata_human.layers:
        adata_human.layers["counts"] = adata_human.X.copy()

    print(f"  Final human reference: {adata_human.shape}")
    print_memory("after concat")

    # Summary stats
    print(f"\n{'=' * 70}")
    print("HUMAN REFERENCE SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total cells: {adata_human.n_obs:,}")
    print(f"  Total genes: {adata_human.n_vars}")
    if "supercluster_term" in adata_human.obs.columns:
        print(f"\n  Cells per supercluster:")
        for sc_name, count in adata_human.obs["supercluster_term"].value_counts().items():
            print(f"    {sc_name}: {count:,}")
    if "cell_class" in adata_human.obs.columns:
        print(f"\n  Cells per class:")
        for cls, count in adata_human.obs["cell_class"].value_counts().items():
            print(f"    {cls}: {count:,}")

    # Save
    out_path = output_dir / "human_reference.h5ad"
    print(f"\n  Saving to {out_path}...")
    adata_human.write(out_path)
    print(f"  Done. File size: {out_path.stat().st_size / 1e9:.2f} GB")
    print_memory("end")


if __name__ == "__main__":
    main()
