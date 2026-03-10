#!/usr/bin/env python3
"""
Step 1H: Preprocess Human Interneuron Data for scVI Atlas Matching

This script:
1. Loads Siletti et al. (2023) human brain reference atlas (interneuron + decoy superclusters)
2. Loads Lee & Dalley human Patch-seq GABAergic interneuron data (DANDI 000636)
3. Converts Patch-seq log2(CPM+1) expression to integer pseudo-counts
4. Maps gene IDs (Ensembl → symbols), intersects gene spaces
5. Selects HVGs from reference only, force-includes ion channel genes
6. Saves preprocessed h5ad files ready for scVI training

Config: atlas_matching/configs/human_interneuron_config.yaml
Linear: TRA-30
"""

import sys
import argparse
import warnings
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "atlas_matching"))


# =============================================================================
# Reference Atlas Loading
# =============================================================================

def load_reference_interneurons(
    h5ad_dir: Path,
    annotation_xlsx: Path,
    interneuron_files: dict,
    subsample: int | None = None,
    random_seed: int = 42,
) -> ad.AnnData:
    """Load and concatenate interneuron superclusters from Siletti et al. (2023).

    Each h5ad file has raw UMI counts in .X (CSR sparse, float32 int-valued).
    Gene IDs are Ensembl in .var index; symbols in .var["Gene"].
    """
    print("[REF] Loading interneuron reference superclusters...")
    rng = np.random.default_rng(random_seed)
    adatas = []

    var_df = None  # Save var info from first file (shared across all)

    for name, info in interneuron_files.items():
        path = h5ad_dir / info["file"]
        print(f"  Loading {name}: {path.name} ...", end=" ", flush=True)
        adata = ad.read_h5ad(path)
        print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")

        # Save var info from first file (all files share the same 59,236 genes)
        if var_df is None:
            var_df = adata.var.copy()

        # Optional subsampling per supercluster
        if subsample is not None and adata.n_obs > subsample:
            idx = rng.choice(adata.n_obs, size=subsample, replace=False)
            adata = adata[idx].copy()
            print(f"    Subsampled to {adata.n_obs:,} cells")

        adatas.append(adata)

    # Concatenate interneuron superclusters
    print(f"\n  Concatenating {len(adatas)} interneuron files...")
    adata_ref = ad.concat(adatas, join="outer")
    del adatas
    gc.collect()

    # Restore var info (ad.concat drops var columns)
    adata_ref.var = var_df.loc[adata_ref.var_names].copy()

    print(f"  Combined interneurons: {adata_ref.n_obs:,} cells x {adata_ref.n_vars:,} genes")
    return adata_ref


def load_decoy_cells(
    h5ad_dir: Path,
    decoy_config: dict,
    random_seed: int = 42,
) -> ad.AnnData:
    """Load subsampled non-interneuron cells as negative controls."""
    print("\n[REF] Loading decoy (non-interneuron) cells...")
    rng = np.random.default_rng(random_seed)
    adatas = []
    var_df = None

    for category, file_list in decoy_config.items():
        for entry in file_list:
            path = h5ad_dir / entry["file"]
            sample_n = entry["sample_n"]
            print(f"  Loading {category}: {path.name} (sampling {sample_n:,}) ...", end=" ", flush=True)

            adata = ad.read_h5ad(path, backed="r")
            n_total = adata.n_obs

            if var_df is None:
                var_df = adata.var.to_memory() if hasattr(adata.var, "to_memory") else adata.var.copy()

            if sample_n < n_total:
                idx = sorted(rng.choice(n_total, size=sample_n, replace=False))
                adata_sub = adata[idx].to_memory()
            else:
                adata_sub = adata.to_memory()

            adata.file.close()
            print(f"{adata_sub.n_obs:,} / {n_total:,} cells")
            adatas.append(adata_sub)

    adata_decoy = ad.concat(adatas, join="outer")
    del adatas
    gc.collect()

    # Restore var info
    if var_df is not None:
        adata_decoy.var = var_df.loc[adata_decoy.var_names].copy()

    print(f"  Total decoy cells: {adata_decoy.n_obs:,}")
    return adata_decoy


def convert_ensembl_to_symbols(adata: ad.AnnData) -> ad.AnnData:
    """Convert .var index from Ensembl IDs to gene symbols using .var['Gene'] column."""
    print("\n[GENES] Converting Ensembl IDs to gene symbols...")

    if "Gene" not in adata.var.columns:
        raise ValueError("Expected .var['Gene'] column with gene symbols")

    # Save original Ensembl IDs
    adata.var["ensembl_id"] = adata.var.index.copy()

    # Map to symbols, handle duplicates by keeping first occurrence
    symbols = adata.var["Gene"].astype(str)
    n_before = len(symbols)

    # Drop genes with missing/empty symbols
    valid_mask = (symbols != "") & (symbols != "nan") & symbols.notna()
    adata = adata[:, valid_mask].copy()
    symbols = adata.var["Gene"].astype(str)
    print(f"  Valid gene symbols: {len(symbols):,} / {n_before:,}")

    # Handle duplicates: keep first occurrence
    dup_mask = symbols.duplicated(keep="first")
    if dup_mask.any():
        n_dup = dup_mask.sum()
        print(f"  Removing {n_dup} duplicate gene symbols")
        adata = adata[:, ~dup_mask].copy()
        symbols = adata.var["Gene"].astype(str)

    # Set gene symbols as index
    adata.var.index = symbols.values
    adata.var_names_make_unique()
    print(f"  Final: {adata.n_vars:,} genes with unique symbols")

    return adata


# =============================================================================
# Patch-seq Query Loading
# =============================================================================

def load_human_patchseq(
    expression_csv: Path,
    metadata_csv: Path,
) -> ad.AnnData:
    """Load human Patch-seq GABAergic interneurons and convert to pseudo-counts.

    Expression is log2(CPM+1). We inverse-transform: round(2^x - 1) to get
    integer pseudo-counts (~1M total per cell). scVI's per-cell size factor
    handles the scaling difference vs 10x reference.
    """
    print("\n[QUERY] Loading human Patch-seq data...")

    # Load expression
    print(f"  Reading expression: {expression_csv.name} ...", end=" ", flush=True)
    expr_df = pd.read_csv(expression_csv, index_col="specimen_id")
    print(f"{expr_df.shape[0]} cells x {expr_df.shape[1]} genes")

    # Convert specimen_id to string for consistent indexing
    expr_df.index = expr_df.index.astype(str)

    # Inverse-transform log2(CPM+1) → integer pseudo-counts
    print("  Converting log2(CPM+1) → pseudo-counts: round(2^x - 1) ...", end=" ", flush=True)
    cpm_values = np.power(2.0, expr_df.values) - 1.0
    pseudo_counts = np.round(cpm_values).astype(np.int32)
    # Clamp any negative values (numerical noise) to 0
    pseudo_counts = np.maximum(pseudo_counts, 0)
    print(f"done. Row sums: {pseudo_counts.sum(axis=1).mean():.0f} mean")

    # Create AnnData
    adata = ad.AnnData(
        X=sparse.csr_matrix(pseudo_counts),
        obs=pd.DataFrame(index=expr_df.index),
        var=pd.DataFrame(index=expr_df.columns),
    )
    adata.obs_names = expr_df.index.tolist()
    adata.var_names = expr_df.columns.tolist()
    adata.layers["counts"] = adata.X.copy()

    del expr_df, cpm_values, pseudo_counts
    gc.collect()

    # Load and merge metadata
    print(f"  Reading metadata: {metadata_csv.name} ...", end=" ", flush=True)
    meta = pd.read_csv(metadata_csv)
    meta["specimen_id"] = meta["specimen_id"].astype(str)
    meta = meta.set_index("specimen_id")
    print(f"{len(meta)} cells")

    # Match metadata to expression
    common_ids = adata.obs_names.intersection(meta.index)
    print(f"  Matched: {len(common_ids)} / {adata.n_obs} cells")

    if len(common_ids) < adata.n_obs:
        adata = adata[common_ids].copy()

    # Add metadata columns to .obs
    meta_cols = [
        "subclass_label",
        "transcriptomic_type",
        "tree_call",
        "donor",
        "condition",
        "brain_region",
        "quality_score",
        "nwb_id",
    ]
    for col in meta_cols:
        if col in meta.columns:
            adata.obs[col] = meta.loc[adata.obs_names, col].values

    # Add tech_batch
    adata.obs["tech_batch"] = "PatchSeq_Human_GABA"
    adata.obs["tech_batch"] = adata.obs["tech_batch"].astype("category")

    # Add supercluster_term for compatibility with reference
    adata.obs["supercluster_term"] = "PatchSeq_query"

    print(f"  Query AnnData: {adata.shape}")
    print(f"  Subclass distribution:")
    if "subclass_label" in adata.obs.columns:
        for lbl, cnt in adata.obs["subclass_label"].value_counts().items():
            print(f"    {lbl}: {cnt}")

    return adata


# =============================================================================
# Ion Channel Genes (Human)
# =============================================================================

HUMAN_ION_CHANNEL_GENES = [
    # Sodium channels
    "SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A",
    # Potassium Kv1
    "KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6",
    # Potassium Kv2
    "KCNB1", "KCNB2",
    # Potassium Kv3
    "KCNC1", "KCNC2", "KCNC3", "KCNC4",
    # Potassium Kv4
    "KCND1", "KCND2", "KCND3",
    # Calcium
    "CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E",
    # HCN
    "HCN1", "HCN2", "HCN3", "HCN4",
]

HUMAN_FORCED_GENES = [
    # AIS anchoring / clustering
    "ANK3", "SPTAN1", "SPTBN4", "NFASC", "CNTNAP2",
    # Nav modulators
    "FGF14", "FGF13",
    # Neuronal markers (neuron vs glia discriminators)
    "SNAP25", "RBFOX3", "SYP", "SYT1",
    # Glial markers (help model learn neuron-glia boundary)
    "GFAP", "ALDOC", "MOG", "PDGFRA",
    # Interneuron subclass markers
    "GAD1", "GAD2", "SLC32A1",  # GABAergic
    "PVALB", "SST", "VIP", "LAMP5", "LHX6",  # Subclass markers
    "CCK", "RELN", "ADARB2",  # CGE markers
    # Excitatory markers (for decoy separation)
    "SLC17A7", "SLC17A6",  # VGLUT1, VGLUT2
]


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Step 1H: Preprocess human interneuron data for scVI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "atlas_matching" / "results" / "human_interneuron" / "preprocessed",
        help="Output directory for preprocessed h5ad files",
    )
    parser.add_argument(
        "--n-hvgs", type=int, default=3000, help="Number of HVGs from reference"
    )
    parser.add_argument(
        "--hvg-ncells",
        type=int,
        default=200_000,
        help="Number of reference cells to subsample for HVG selection",
    )
    parser.add_argument(
        "--subsample-interneurons",
        type=int,
        default=None,
        help="Subsample each interneuron supercluster to N cells (for testing)",
    )
    parser.add_argument(
        "--no-decoys",
        action="store_true",
        help="Skip loading decoy (non-interneuron) cells",
    )
    parser.add_argument(
        "--random-seed", type=int, default=42, help="Random seed"
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Paths from config
    H5AD_DIR = Path("/mnt/data0/HumanBrainCellType/SuperTypeRawDat")
    ANNOTATION_XLSX = Path("/mnt/data0/HumanBrainCellType/annotation.xlsx")
    EXPRESSION_CSV = Path(
        "/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/expression_by_cell.csv.gz"
    )
    METADATA_CSV = Path(
        "/home/jw3514/Work/NeurSim/EphysSumStats/workspace/HumanCortexGaba/cell_metadata.csv"
    )

    interneuron_files = {
        "CGE": {"file": "Supercluster_CGE_interneuron.h5ad"},
        "MGE": {"file": "Supercluster_MGE_interneuron.h5ad"},
        "LAMP5_LHX6_Chandelier": {"file": "Supercluster_LAMP5_LHX6_and_Chandelier.h5ad"},
    }

    decoy_config = {
        "glutamatergic": [
            {"file": "Supercluster_Upper_layer_intratelencephalic.h5ad", "sample_n": 5000},
            {"file": "Supercluster_Deep_layer_intratelencephalic.h5ad", "sample_n": 5000},
        ],
        "non_neuronal": [
            {"file": "Supercluster_Astrocyte.h5ad", "sample_n": 3000},
            {"file": "Supercluster_Oligodendrocyte.h5ad", "sample_n": 3000},
            {"file": "Supercluster_Microglia.h5ad", "sample_n": 2000},
        ],
    }

    print("=" * 80)
    print("STEP 1H: PREPROCESS HUMAN INTERNEURON DATA FOR SCVI")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  HVGs: {args.n_hvgs}")
    print(f"  HVG subsample: {args.hvg_ncells:,}")
    print(f"  Subsample interneurons: {args.subsample_interneurons}")
    print(f"  Include decoys: {not args.no_decoys}")
    print(f"  Output: {args.output_dir}")
    print("=" * 80)

    # ---- Phase 1.1: Load reference ----
    adata_interneurons = load_reference_interneurons(
        H5AD_DIR,
        ANNOTATION_XLSX,
        interneuron_files,
        subsample=args.subsample_interneurons,
        random_seed=args.random_seed,
    )

    if not args.no_decoys:
        adata_decoys = load_decoy_cells(
            H5AD_DIR,
            decoy_config,
            random_seed=args.random_seed,
        )
        # Concatenate interneurons + decoys
        print("\n[REF] Concatenating interneurons + decoys...")
        var_backup = adata_interneurons.var.copy()
        adata_ref = ad.concat([adata_interneurons, adata_decoys], join="outer")
        adata_ref.var = var_backup.loc[adata_ref.var_names].copy()
        del adata_interneurons, adata_decoys, var_backup
    else:
        adata_ref = adata_interneurons
        del adata_interneurons

    gc.collect()
    print(f"  Reference total: {adata_ref.n_obs:,} cells x {adata_ref.n_vars:,} genes")

    # ---- Convert Ensembl → gene symbols ----
    adata_ref = convert_ensembl_to_symbols(adata_ref)

    # Store raw counts in layers
    if "counts" not in adata_ref.layers:
        adata_ref.layers["counts"] = adata_ref.X.copy()

    # Add tech_batch for reference
    adata_ref.obs["tech_batch"] = "Siletti_Atlas"
    adata_ref.obs["tech_batch"] = adata_ref.obs["tech_batch"].astype("category")

    # Keep essential obs columns
    essential_ref_cols = [
        "tech_batch", "supercluster_term", "cluster_id", "subcluster_id",
        "donor_id", "sample_id", "roi", "tissue",
    ]
    cols_to_keep = [c for c in essential_ref_cols if c in adata_ref.obs.columns]
    adata_ref.obs = adata_ref.obs[cols_to_keep].copy()

    # Clear unnecessary attributes
    adata_ref.uns.clear()
    adata_ref.obsm.clear()
    adata_ref.varm.clear()
    adata_ref.obsp.clear()
    adata_ref.varp.clear()
    gc.collect()

    print(f"\n  Reference after cleanup: {adata_ref.shape}")
    print(f"  Obs columns: {list(adata_ref.obs.columns)}")
    if "supercluster_term" in adata_ref.obs.columns:
        print(f"  Supercluster distribution:")
        for lbl, cnt in adata_ref.obs["supercluster_term"].value_counts().items():
            print(f"    {lbl}: {cnt:,}")

    # ---- Phase 1.2: Load query ----
    adata_query = load_human_patchseq(EXPRESSION_CSV, METADATA_CSV)

    # ---- Phase 1.3: Gene intersection + HVG selection ----
    print("\n[GENES] Intersecting gene spaces...")
    ref_genes = set(adata_ref.var_names)
    query_genes = set(adata_query.var_names)
    common_genes = sorted(ref_genes & query_genes)
    print(f"  Reference genes: {len(ref_genes):,}")
    print(f"  Query genes: {len(query_genes):,}")
    print(f"  Common genes: {len(common_genes):,}")

    # Check ion channel gene coverage
    ion_in_common = [g for g in HUMAN_ION_CHANNEL_GENES if g in common_genes]
    print(f"  Ion channel genes in common: {len(ion_in_common)} / {len(HUMAN_ION_CHANNEL_GENES)}")
    missing_ion = [g for g in HUMAN_ION_CHANNEL_GENES if g not in common_genes]
    if missing_ion:
        print(f"  WARNING: Missing ion channel genes: {missing_ion}")

    # Subset both to common genes
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()
    gc.collect()

    # Ensure counts layers are present after subsetting
    if "counts" not in adata_ref.layers:
        adata_ref.layers["counts"] = adata_ref.X.copy()
    if "counts" not in adata_query.layers:
        adata_query.layers["counts"] = adata_query.X.copy()

    print(f"  Reference after intersection: {adata_ref.shape}")
    print(f"  Query after intersection: {adata_query.shape}")

    # ---- HVG selection from reference only ----
    print(f"\n[HVG] Selecting {args.n_hvgs} HVGs from reference (subsample={args.hvg_ncells:,})...")
    rng = np.random.default_rng(args.random_seed)
    n_ref = adata_ref.n_obs
    n_sub = min(args.hvg_ncells, n_ref)

    if n_sub < n_ref:
        sub_idx = rng.choice(n_ref, size=n_sub, replace=False)
        ref_hvg = adata_ref[sub_idx, :].copy()
    else:
        ref_hvg = adata_ref.copy()

    # seurat_v3 expects counts in .X
    ref_hvg.X = ref_hvg.layers["counts"]
    sc.pp.highly_variable_genes(ref_hvg, n_top_genes=args.n_hvgs, flavor="seurat_v3", subset=False)
    hvg_genes = ref_hvg.var_names[ref_hvg.var["highly_variable"]].tolist()
    del ref_hvg
    gc.collect()
    print(f"  Selected {len(hvg_genes)} HVGs")

    # Force-include ion channel + marker genes
    forced_genes = set(HUMAN_ION_CHANNEL_GENES) | set(HUMAN_FORCED_GENES)
    forced_present = [g for g in forced_genes if g in adata_ref.var_names]
    final_genes = sorted(set(hvg_genes) | set(forced_present))
    print(f"  Forced genes: {len(forced_present)} / {len(forced_genes)} present")
    print(f"  Final gene set: {len(final_genes)}")

    # Subset both datasets to final genes
    print("[SUBSET] Subsetting to final gene set...")
    gene_mask = adata_ref.var_names.isin(final_genes)

    adata_ref_final = ad.AnnData(
        X=adata_ref.layers["counts"][:, gene_mask].copy(),
        obs=adata_ref.obs.copy(),
        var=adata_ref.var.loc[gene_mask].copy(),
    )
    adata_ref_final.layers["counts"] = adata_ref_final.X.copy()

    adata_query_final = ad.AnnData(
        X=adata_query.layers["counts"][:, gene_mask].copy(),
        obs=adata_query.obs.copy(),
        var=adata_query.var.loc[gene_mask].copy(),
    )
    adata_query_final.layers["counts"] = adata_query_final.X.copy()

    del adata_ref, adata_query
    gc.collect()

    print(f"  Reference final: {adata_ref_final.shape}")
    print(f"  Query final: {adata_query_final.shape}")

    # ---- Save outputs ----
    print(f"\n[SAVE] Saving preprocessed data to {args.output_dir}")
    ref_path = args.output_dir / "human_reference_preprocessed.h5ad"
    query_path = args.output_dir / "human_patchseq_query_preprocessed.h5ad"
    genes_path = args.output_dir / "human_genes.txt"

    adata_ref_final.write(ref_path)
    print(f"  Reference: {ref_path} ({adata_ref_final.shape})")

    adata_query_final.write(query_path)
    print(f"  Query: {query_path} ({adata_query_final.shape})")

    with open(genes_path, "w") as f:
        f.write("\n".join(final_genes))
    print(f"  Gene list: {genes_path} ({len(final_genes)} genes)")

    # Summary
    print(f"\n{'=' * 80}")
    print("PREPROCESSING COMPLETE")
    print(f"{'=' * 80}")
    print(f"  Reference: {adata_ref_final.n_obs:,} cells x {adata_ref_final.n_vars:,} genes")
    print(f"  Query: {adata_query_final.n_obs:,} cells x {adata_query_final.n_vars:,} genes")
    if "supercluster_term" in adata_ref_final.obs.columns:
        print(f"\n  Reference composition:")
        for lbl, cnt in adata_ref_final.obs["supercluster_term"].value_counts().items():
            print(f"    {lbl}: {cnt:,}")
    if "subclass_label" in adata_query_final.obs.columns:
        print(f"\n  Query composition:")
        for lbl, cnt in adata_query_final.obs["subclass_label"].value_counts().items():
            print(f"    {lbl}: {cnt}")
    print(f"\nNext step: python scripts/02h_train_scvi_human.py")


if __name__ == "__main__":
    main()
