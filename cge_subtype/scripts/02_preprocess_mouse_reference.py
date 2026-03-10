#!/usr/bin/env python3
"""
Step 2: Preprocess mouse reference from Allen Brain Cell Atlas (WMB-10Xv3).

Memory-safe approach: load metadata first to identify target cells, subsample
barcodes BEFORE loading expression, then load only needed cells from partitioned
h5ad files.

Output: results/cross_species/preprocessed/mouse_reference.h5ad
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


def load_ortholog_mapping(ortho_dir):
    """Load mouse->human symbol mapping from ortholog CSV."""
    ortho_path = Path(ortho_dir) / "ortholog_mapping.csv"
    df = pd.read_csv(ortho_path)

    # Build mouse_symbol -> human_symbol dict (prefer 1:1, take first for many)
    # Group by mouse_symbol, pick best human mapping
    mapping = {}
    for _, row in df.iterrows():
        ms = row["mouse_symbol"]
        hs = row["human_symbol"]
        if pd.isna(ms) or pd.isna(hs):
            continue
        if ms not in mapping:
            mapping[ms] = hs
        # If already exists, prefer 1:1 mapping
        elif row.get("is_one_to_one", False):
            mapping[ms] = hs

    print(f"  Loaded {len(mapping)} mouse->human symbol mappings")
    return mapping


def filter_metadata(config):
    """Load and filter Allen metadata to target regions and cell classes."""
    meta_path = config["mouse_reference"]["metadata_csv"]
    print(f"  Loading metadata from {meta_path}...")
    meta = pd.read_csv(meta_path, low_memory=False)
    print(f"  Total cells in metadata: {len(meta):,}")

    # Filter by dataset (WMB-10Xv3 only — expression files are from this dataset)
    dataset_filter = config["mouse_reference"].get("dataset_filter", None)
    if dataset_filter and "dataset_label" in meta.columns:
        meta = meta[meta["dataset_label"] == dataset_filter].copy()
        print(f"  After dataset filter ({dataset_filter}): {len(meta):,} cells")

    # Filter by region
    region_prefixes = config["mouse_reference"]["region_filter"]["include_prefixes"]
    if "region_of_interest_acronym" in meta.columns:
        roi_col = "region_of_interest_acronym"
    else:
        # Try alternative column names
        for col in meta.columns:
            if "region" in col.lower():
                roi_col = col
                break
        else:
            print("  WARNING: No region column found, skipping region filter")
            roi_col = None

    if roi_col:
        region_mask = pd.Series(False, index=meta.index)
        for prefix in region_prefixes:
            region_mask |= meta[roi_col].astype(str).str.startswith(prefix)
        meta = meta[region_mask].copy()
        print(f"  After region filter: {len(meta):,} cells")

    # Filter by cell class
    include_classes = config["mouse_reference"]["cell_class_filter"]["include_classes"]
    if "class" in meta.columns:
        class_mask = meta["class"].isin(include_classes)
        meta = meta[class_mask].copy()
        print(f"  After class filter: {len(meta):,} cells")

    return meta


def subsample_metadata(meta, config, seed=42):
    """Subsample metadata to target N per cell class."""
    subsampling = config["mouse_reference"]["subsampling"]
    np.random.seed(seed)

    # Classify cells into broad categories
    gaba_classes = [c for c in meta["class"].unique()
                    if any(k in c.lower() for k in ["gaba", "lamp5", "sncg", "vip", "pvalb", "sst"])]
    glut_classes = [c for c in meta["class"].unique()
                    if any(k in c.lower() for k in ["glut", "it-et", "np-ct"])]
    nonneuronal_classes = [c for c in meta["class"].unique()
                           if any(k in c.lower() for k in ["astro", "oligo", "opc"])]

    selected = []

    for classes, label, max_n in [
        (gaba_classes, "GABAergic", subsampling["GABAergic"]),
        (glut_classes, "Glutamatergic", subsampling["Glutamatergic"]),
        (nonneuronal_classes, "NonNeuronal", subsampling["NonNeuronal"]),
    ]:
        mask = meta["class"].isin(classes)
        sub = meta[mask]
        n_available = len(sub)
        n_take = min(max_n, n_available)

        if n_take < n_available:
            # Stratified subsample by subclass
            if "subclass" in sub.columns:
                sampled = sub.groupby("subclass", group_keys=False).apply(
                    lambda x: x.sample(
                        n=max(1, int(n_take * len(x) / n_available)),
                        random_state=seed,
                    )
                )
                # Trim to exact target
                if len(sampled) > n_take:
                    sampled = sampled.sample(n=n_take, random_state=seed)
            else:
                sampled = sub.sample(n=n_take, random_state=seed)
        else:
            sampled = sub

        sampled = sampled.copy()
        sampled["cell_class"] = label
        selected.append(sampled)
        print(f"  {label}: {n_available:,} available -> {len(sampled):,} selected")

    meta_sub = pd.concat(selected, ignore_index=False)
    print(f"  Total subsampled: {len(meta_sub):,} cells")
    return meta_sub


def load_expression_for_cells(config, cell_labels, mouse_to_human, human_gene_set):
    """Load expression data only for selected cells from partitioned h5ad files."""
    expr_dir = Path(config["mouse_reference"]["expression_dir"])
    expr_files = config["mouse_reference"]["expression_files"]

    cell_labels_set = set(cell_labels)
    chunks = []

    for fname in expr_files:
        fpath = expr_dir / fname
        if not fpath.exists():
            print(f"  WARNING: {fpath} not found, skipping")
            continue

        print(f"\n  Loading {fname}...")
        print_memory(f"before {fname}")

        # Load in backed mode
        adata_backed = sc.read_h5ad(fpath, backed="r")
        n_total = adata_backed.n_obs

        # Find which cells in this file match our target set
        obs_names = adata_backed.obs_names.tolist()
        # Also check cell_label column if it exists
        if "cell_label" in adata_backed.obs.columns:
            file_cell_labels = adata_backed.obs["cell_label"].values
        else:
            file_cell_labels = np.array(obs_names)

        keep_indices = []
        keep_cell_labels = []
        for i, cl in enumerate(file_cell_labels):
            if cl in cell_labels_set:
                keep_indices.append(i)
                keep_cell_labels.append(cl)

        if len(keep_indices) == 0:
            print(f"    No matching cells in {fname}")
            adata_backed.file.close()
            del adata_backed
            continue

        keep_indices = np.array(keep_indices)
        print(f"    Found {len(keep_indices):,} matching cells (of {n_total:,})")

        # Get mouse gene symbols from .var
        # Allen WMB-10Xv3 h5ad files use Ensembl IDs as .var index,
        # with gene symbols in .var['gene_symbol'] column
        var_df = adata_backed.var.copy()
        if "gene_symbol" in var_df.columns:
            mouse_genes = var_df["gene_symbol"].tolist()
        elif "Gene" in var_df.columns:
            mouse_genes = var_df["Gene"].tolist()
        else:
            mouse_genes = var_df.index.tolist()
        print(f"    Gene ID sample: {mouse_genes[:3]} (from {'gene_symbol col' if 'gene_symbol' in var_df.columns else 'index'})")

        # Map mouse genes to human symbols
        gene_mapping = {}  # col_index -> human_symbol
        for j, mg in enumerate(mouse_genes):
            hs = mouse_to_human.get(mg, None)
            if hs and hs in human_gene_set:
                gene_mapping[j] = hs

        keep_gene_indices = sorted(gene_mapping.keys())
        mapped_symbols = [gene_mapping[j] for j in keep_gene_indices]
        print(f"    Mapped {len(keep_gene_indices)} mouse genes -> human symbols")

        # Handle many-to-one: if multiple mouse genes map to same human gene,
        # we'll sum their counts later
        keep_gene_indices = np.array(keep_gene_indices)

        # Extract expression matrix (cells x genes)
        X_sub = adata_backed.X[keep_indices][:, keep_gene_indices]
        if not sparse.issparse(X_sub):
            X_sub = sparse.csr_matrix(X_sub)
        else:
            X_sub = X_sub.copy()

        adata_backed.file.close()
        del adata_backed
        gc.collect()

        # Build AnnData with human gene symbols
        var_new = pd.DataFrame(index=mapped_symbols)
        obs_new = pd.DataFrame(index=keep_cell_labels)

        adata_chunk = ad.AnnData(X=X_sub, obs=obs_new, var=var_new)

        # Handle duplicate human gene symbols (many-to-one): sum columns
        dup_genes = adata_chunk.var_names[adata_chunk.var_names.duplicated(keep=False)]
        if len(dup_genes) > 0:
            unique_dup = dup_genes.unique()
            print(f"    Summing {len(unique_dup)} duplicate human symbols (many-to-one)")

            # For each duplicate, sum the columns
            X_dense_parts = []
            non_dup_mask = ~adata_chunk.var_names.duplicated(keep=False)
            X_non_dup = adata_chunk.X[:, non_dup_mask]
            non_dup_names = adata_chunk.var_names[non_dup_mask].tolist()

            summed_names = []
            summed_cols = []
            for gene in unique_dup:
                gene_mask = adata_chunk.var_names == gene
                col_sum = np.asarray(adata_chunk.X[:, gene_mask].sum(axis=1)).ravel()
                summed_cols.append(col_sum)
                summed_names.append(gene)

            # Reconstruct
            if len(summed_cols) > 0:
                X_summed = np.column_stack(summed_cols)
                if sparse.issparse(X_non_dup):
                    X_combined = sparse.hstack([X_non_dup, sparse.csr_matrix(X_summed)])
                else:
                    X_combined = np.hstack([X_non_dup, X_summed])
                all_names = non_dup_names + summed_names
                var_combined = pd.DataFrame(index=all_names)
                adata_chunk = ad.AnnData(
                    X=X_combined if sparse.issparse(X_combined) else sparse.csr_matrix(X_combined),
                    obs=obs_new,
                    var=var_combined,
                )

        chunks.append(adata_chunk)
        print(f"    Chunk: {adata_chunk.shape}")
        print_memory(f"after {fname}")

    if not chunks:
        raise RuntimeError("No expression data loaded!")

    # Concatenate chunks from different files
    print(f"\n  Concatenating {len(chunks)} expression chunks...")
    # Find common genes across chunks
    common_genes = set(chunks[0].var_names)
    for chunk in chunks[1:]:
        common_genes &= set(chunk.var_names)
    common_genes = sorted(common_genes)
    print(f"  Common genes across files: {len(common_genes)}")

    for i in range(len(chunks)):
        chunks[i] = chunks[i][:, common_genes].copy()

    adata_expr = ad.concat(chunks, join="inner")
    adata_expr.var_names_make_unique()
    del chunks
    gc.collect()

    return adata_expr, common_genes


def main():
    parser = argparse.ArgumentParser(description="Preprocess mouse reference")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    config = load_config()
    output_dir = args.output_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["preprocessed"])
    output_dir.mkdir(parents=True, exist_ok=True)
    ortho_dir = ATLAS_MATCHING_DIR / config["output"]["subdirs"]["orthologs"]

    print("=" * 70)
    print("STEP 2: PREPROCESS MOUSE REFERENCE (ALLEN WMB-10Xv3)")
    print("=" * 70)
    print_memory("start")

    # Load ortholog mapping
    print("\n[1/5] Loading ortholog mapping...")
    mouse_to_human = load_ortholog_mapping(ortho_dir)

    # Load human reference gene set (to ensure overlap)
    human_ref_path = output_dir / "human_reference.h5ad"
    if human_ref_path.exists():
        print(f"\n[2/5] Loading human reference gene set from {human_ref_path.name}...")
        adata_human_backed = sc.read_h5ad(human_ref_path, backed="r")
        human_gene_set = set(adata_human_backed.var_names.tolist())
        adata_human_backed.file.close()
        del adata_human_backed
        print(f"  Human reference genes: {len(human_gene_set)}")
    else:
        # Fall back to all orthologs
        print(f"\n[2/5] Human reference not found, using all ortholog genes...")
        human_gene_set = set(mouse_to_human.values())

    # Filter and subsample metadata
    print(f"\n[3/5] Filtering and subsampling metadata...")
    meta = filter_metadata(config)
    meta_sub = subsample_metadata(meta, config, seed=args.seed)
    del meta
    gc.collect()

    # Get cell labels for expression loading
    if "cell_label" in meta_sub.columns:
        cell_labels = meta_sub["cell_label"].tolist()
    else:
        cell_labels = meta_sub.index.tolist()

    # Load expression for selected cells
    print(f"\n[4/5] Loading expression data for {len(cell_labels):,} cells...")
    print_memory("before expr load")
    adata_expr, common_genes = load_expression_for_cells(
        config, cell_labels, mouse_to_human, human_gene_set
    )
    print_memory("after expr load")

    # Merge metadata into expression AnnData
    print(f"\n[5/5] Merging metadata...")

    # Match metadata to expression cells
    if "cell_label" in meta_sub.columns:
        meta_indexed = meta_sub.set_index("cell_label", drop=False)
    else:
        meta_indexed = meta_sub

    common_cells = set(adata_expr.obs_names) & set(meta_indexed.index)
    print(f"  Matched {len(common_cells):,} cells between expression and metadata")

    if len(common_cells) < len(adata_expr):
        adata_expr = adata_expr[list(common_cells)].copy()

    # Add metadata columns
    essential_meta_cols = [
        "class", "subclass", "supertype", "cluster",
        "neurotransmitter", "library_label", "region_of_interest_acronym",
    ]
    for col in essential_meta_cols:
        if col in meta_indexed.columns:
            adata_expr.obs[col] = meta_indexed.loc[adata_expr.obs_names, col].values

    # Add cell_class from subsampling
    if "cell_class" in meta_indexed.columns:
        adata_expr.obs["cell_class"] = meta_indexed.loc[adata_expr.obs_names, "cell_class"].values

    # Add species and sample_species
    adata_expr.obs["species"] = "mouse"
    if "library_label" in adata_expr.obs.columns:
        adata_expr.obs["sample_species"] = (
            adata_expr.obs["library_label"].astype(str) + "_mouse"
        )
    else:
        adata_expr.obs["sample_species"] = "allen_mouse"

    # Ensure counts layer
    adata_expr.layers["counts"] = adata_expr.X.copy()

    # Subset to genes present in human reference
    genes_in_human = [g for g in adata_expr.var_names if g in human_gene_set]
    if len(genes_in_human) < adata_expr.n_vars:
        print(f"  Subsetting to {len(genes_in_human)} genes present in human reference")
        adata_expr = adata_expr[:, genes_in_human].copy()
        adata_expr.layers["counts"] = adata_expr.X.copy()

    # Summary
    print(f"\n{'=' * 70}")
    print("MOUSE REFERENCE SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total cells: {adata_expr.n_obs:,}")
    print(f"  Total genes (human symbols): {adata_expr.n_vars}")
    if "cell_class" in adata_expr.obs.columns:
        print(f"\n  Cells per class:")
        for cls, count in adata_expr.obs["cell_class"].value_counts().items():
            print(f"    {cls}: {count:,}")
    if "subclass" in adata_expr.obs.columns:
        print(f"\n  Top subclasses:")
        for sc_name, count in adata_expr.obs["subclass"].value_counts().head(15).items():
            print(f"    {sc_name}: {count:,}")

    # Save
    out_path = output_dir / "mouse_reference.h5ad"
    print(f"\n  Saving to {out_path}...")
    adata_expr.write(out_path)
    print(f"  Done. File size: {out_path.stat().st_size / 1e9:.2f} GB")
    print_memory("end")


if __name__ == "__main__":
    main()
