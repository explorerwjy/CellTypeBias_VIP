#!/usr/bin/env python3
"""
Step 3: Preprocess mouse Patch-seq query data (V1 + M1).

Loads V1 and M1 Patch-seq datasets, converts mouse gene symbols to human
symbols using ortholog mapping, flags CCKBC cells, and saves for scArches mapping.

Output: results/cross_species/preprocessed/mouse_patchseq_query.h5ad
"""

import sys
import gc
import gzip
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import yaml
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_MATCHING_DIR = SCRIPT_DIR.parent.parent
PROJECT_ROOT = ATLAS_MATCHING_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(ATLAS_MATCHING_DIR))


def load_config():
    config_path = ATLAS_MATCHING_DIR / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_ortholog_mapping(ortho_dir):
    """Load mouse->human symbol mapping from ortholog CSV."""
    ortho_path = Path(ortho_dir) / "ortholog_mapping.csv"
    df = pd.read_csv(ortho_path)
    mapping = {}
    for _, row in df.iterrows():
        ms = row["mouse_symbol"]
        hs = row["human_symbol"]
        if pd.isna(ms) or pd.isna(hs):
            continue
        if ms not in mapping:
            mapping[ms] = hs
        elif row.get("is_one_to_one", False):
            mapping[ms] = hs
    return mapping


def load_v1_patchseq(data_dir, mouse_to_human, target_genes):
    """Load V1 Patch-seq, convert gene symbols, subset to target genes."""
    print("\n  [V1] Loading V1 Patch-seq...")
    counts_path = Path(data_dir) / "V1_patchseq_counts.csv"
    meta_path = Path(data_dir) / "V1_patchseq_metadata.csv"

    # Load counts (genes x cells)
    counts_df = pd.read_csv(counts_path, index_col=0)
    counts_df = counts_df.T  # -> cells x genes
    print(f"    Expression: {counts_df.shape[0]} cells x {counts_df.shape[1]} genes")

    # Load metadata
    metadata = pd.read_csv(meta_path, sep=",")

    # Match metadata to expression
    if "transcriptomics_sample_id" in metadata.columns:
        metadata = metadata.set_index("transcriptomics_sample_id", drop=False)

    common_cells = sorted(set(counts_df.index) & set(metadata.index))
    counts_df = counts_df.loc[common_cells]
    metadata = metadata.loc[common_cells]
    print(f"    Matched cells: {len(common_cells)}")

    # Convert mouse gene symbols to human
    mouse_genes = counts_df.columns.tolist()
    mapped_genes = {}
    for mg in mouse_genes:
        hs = mouse_to_human.get(mg)
        if hs and hs in target_genes:
            if hs not in mapped_genes:
                mapped_genes[hs] = []
            mapped_genes[hs].append(mg)

    # Build human-symbol expression matrix (sum many-to-one)
    human_expr = pd.DataFrame(index=counts_df.index)
    for hs, mouse_cols in mapped_genes.items():
        if len(mouse_cols) == 1:
            human_expr[hs] = counts_df[mouse_cols[0]].values
        else:
            human_expr[hs] = counts_df[mouse_cols].sum(axis=1).values

    print(f"    Mapped {len(mapped_genes)} human genes (from {sum(len(v) for v in mapped_genes.values())} mouse genes)")

    # Create AnnData
    X = sparse.csr_matrix(human_expr.values.astype(np.float32))
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=human_expr.index),
        var=pd.DataFrame(index=human_expr.columns),
    )
    adata.layers["counts"] = adata.X.copy()

    # Add metadata
    for col in metadata.columns:
        if col in ["transcriptomics_sample_id"]:
            continue
        adata.obs[col] = metadata.loc[adata.obs_names, col].values

    adata.obs["dataset"] = "V1"
    adata.obs["species"] = "mouse_patchseq"
    adata.obs["sample_species"] = "patchseq_V1_mouse"

    return adata


def load_m1_patchseq(data_dir, mouse_to_human, target_genes):
    """Load M1 Patch-seq, convert gene symbols, subset to target genes."""
    print("\n  [M1] Loading M1 Patch-seq...")
    counts_path = Path(data_dir) / "M1_patchseq_counts.csv.gz"
    meta_path = Path(data_dir) / "M1_patchseq_metadata.csv"

    # Load counts (genes x cells, gzipped)
    counts_df = pd.read_csv(counts_path, index_col=0, compression="gzip")
    counts_df = counts_df.T  # -> cells x genes
    print(f"    Expression: {counts_df.shape[0]} cells x {counts_df.shape[1]} genes")

    # Load metadata (tab-separated)
    try:
        metadata = pd.read_csv(meta_path, sep="\t")
    except Exception:
        metadata = pd.read_csv(meta_path, sep=",")

    # Match metadata to expression
    if "Cell" in metadata.columns:
        metadata = metadata.set_index("Cell", drop=False)

    common_cells = sorted(set(counts_df.index) & set(metadata.index))
    counts_df = counts_df.loc[common_cells]
    metadata = metadata.loc[common_cells]
    print(f"    Matched cells: {len(common_cells)}")

    # Convert mouse gene symbols to human
    mouse_genes = counts_df.columns.tolist()
    mapped_genes = {}
    for mg in mouse_genes:
        hs = mouse_to_human.get(mg)
        if hs and hs in target_genes:
            if hs not in mapped_genes:
                mapped_genes[hs] = []
            mapped_genes[hs].append(mg)

    # Build human-symbol expression matrix
    human_expr = pd.DataFrame(index=counts_df.index)
    for hs, mouse_cols in mapped_genes.items():
        if len(mouse_cols) == 1:
            human_expr[hs] = counts_df[mouse_cols[0]].values
        else:
            human_expr[hs] = counts_df[mouse_cols].sum(axis=1).values

    print(f"    Mapped {len(mapped_genes)} human genes")

    # Create AnnData
    X = sparse.csr_matrix(human_expr.values.astype(np.float32))
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=human_expr.index),
        var=pd.DataFrame(index=human_expr.columns),
    )
    adata.layers["counts"] = adata.X.copy()

    # Add metadata
    for col in metadata.columns:
        if col in ["Cell"]:
            continue
        adata.obs[col] = metadata.loc[adata.obs_names, col].values

    adata.obs["dataset"] = "M1"
    adata.obs["species"] = "mouse_patchseq"
    adata.obs["sample_species"] = "patchseq_M1_mouse"

    return adata


def flag_cckbc_cells(adata, config):
    """Flag cells as CCKBCs based on Gouwens et al. 2020 classification.

    CCKBC types (see docs/cckbc_definition.md):
      1. Sncg subclass (all types) — large CCK basket cells
      2. Vip Sncg supertype — small CCK basket cells
      3. Vip Serpinf1 supertype — deep CCK+ cells

    Excludes other Vip types (Mybpc1, Col15a1, Chat, etc.) and all
    Lamp5/Sst/Pvalb types.
    """
    cckbc_def = config["mouse_patchseq"]["cckbc_definition"]
    sncg_prefixes = cckbc_def["sncg_subclass_prefixes"]       # ["Sncg"]
    vip_supertypes = cckbc_def["vip_cckbc_supertypes"]        # ["Vip Sncg", "Vip Serpinf1"]
    m1_families = cckbc_def["m1_family_include"]               # ["Sncg"]
    exclude_types = cckbc_def["exclude_types"]                 # ["Mybpc1", ...]

    adata.obs["is_cckbc"] = False

    def _is_excluded(val):
        """Check if a type name contains any excluded marker."""
        for exc in exclude_types:
            if exc.lower() in val.lower():
                return True
        return False

    # --- V1: uses "corresponding_AIT2.3.1_alias" column ---
    # Format: "Sncg Vip Itih5", "Sncg Gpr50", "Vip Col15a1 Pde1a", etc.
    alias_col = None
    # Prefer columns with "alias" in name; fall back to "ait" prefix columns
    for col in adata.obs.columns:
        if "alias" in col.lower():
            alias_col = col
            break
    if alias_col is None:
        for col in adata.obs.columns:
            if "ait" in col.lower() and "id" not in col.lower():
                alias_col = col
                break

    if alias_col is not None:
        values = adata.obs[alias_col].astype(str)
        for prefix in sncg_prefixes:
            # Types starting with "Sncg" = Sncg subclass
            mask = values.str.startswith(prefix, na=False) & ~values.apply(_is_excluded)
            if mask.any():
                adata.obs.loc[mask, "is_cckbc"] = True
                print(f"    Flagged {mask.sum()} cells as CCKBC via {alias_col} starting with '{prefix}'")
        # Also check for Serpinf1 in V1 alias (e.g. "Serpinf1 Clrn1", "Serpinf1 Aqp5 Vip")
        serpinf1_mask = values.str.contains("Serpinf1", case=False, na=False) & ~values.apply(_is_excluded)
        if serpinf1_mask.any():
            adata.obs.loc[serpinf1_mask, "is_cckbc"] = True
            print(f"    Flagged {serpinf1_mask.sum()} cells as CCKBC via {alias_col} containing 'Serpinf1'")

    # --- M1: uses "RNA type" and "RNA family" columns ---
    # RNA family: "Vip", "Sncg", "Lamp5", etc.
    # RNA type: "Vip Sncg", "Sncg Npy2r", "Vip Mybpc1_2", etc.
    family_col = None
    type_col = None
    for col in adata.obs.columns:
        col_lower = col.lower().replace(" ", "_")
        if col_lower in ("rna_family", "rnafamily"):
            family_col = col
        elif col_lower in ("rna_type", "rnatype"):
            type_col = col

    if family_col is not None:
        families = adata.obs[family_col].astype(str)
        for fam in m1_families:
            mask = (families == fam)
            if mask.any():
                adata.obs.loc[mask, "is_cckbc"] = True
                print(f"    Flagged {mask.sum()} cells as CCKBC via {family_col} == '{fam}'")

    if type_col is not None:
        types = adata.obs[type_col].astype(str)
        for vip_st in vip_supertypes:
            # Match types starting with the supertype name (e.g. "Vip Sncg", "Vip Serpinf1_1")
            mask = types.str.startswith(vip_st, na=False) & ~types.apply(_is_excluded)
            if mask.any():
                # Don't double-count cells already flagged by family
                new_flags = mask & ~adata.obs["is_cckbc"]
                if new_flags.any():
                    adata.obs.loc[new_flags, "is_cckbc"] = True
                    print(f"    Flagged {new_flags.sum()} additional cells as CCKBC via {type_col} starting with '{vip_st}'")

    n_cckbc = adata.obs["is_cckbc"].sum()
    print(f"  Total CCKBC cells: {n_cckbc} / {len(adata)}")

    # Print breakdown by type for verification
    for col in [alias_col, type_col]:
        if col is not None and col in adata.obs.columns:
            cckbc_types = adata.obs.loc[adata.obs["is_cckbc"], col].value_counts()
            if len(cckbc_types) > 0:
                print(f"\n  CCKBC breakdown by {col}:")
                for typ, cnt in cckbc_types.items():
                    print(f"    {typ}: {cnt}")

    return adata


def main():
    parser = argparse.ArgumentParser(description="Preprocess mouse Patch-seq query")
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    config = load_config()
    output_dir = args.output_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["preprocessed"])
    output_dir.mkdir(parents=True, exist_ok=True)
    ortho_dir = ATLAS_MATCHING_DIR / config["output"]["subdirs"]["orthologs"]

    print("=" * 70)
    print("STEP 3: PREPROCESS MOUSE PATCH-SEQ QUERY (V1 + M1)")
    print("=" * 70)

    # Load ortholog mapping
    print("\n[1/4] Loading ortholog mapping...")
    mouse_to_human = load_ortholog_mapping(ortho_dir)
    print(f"  {len(mouse_to_human)} mouse->human mappings")

    # Load target gene set from human reference
    human_ref_path = output_dir / "human_reference.h5ad"
    if human_ref_path.exists():
        print(f"\n[2/4] Loading target gene set from human reference...")
        adata_human = sc.read_h5ad(human_ref_path, backed="r")
        target_genes = set(adata_human.var_names.tolist())
        adata_human.file.close()
        del adata_human
    else:
        print(f"\n[2/4] Human reference not found, using all ortholog genes...")
        target_genes = set(mouse_to_human.values())
    print(f"  Target gene set: {len(target_genes)} genes")

    # Load V1 and M1 Patch-seq
    data_dir = config["mouse_patchseq"]["data_dir"]

    print(f"\n[3/4] Loading Patch-seq datasets...")
    adata_v1 = load_v1_patchseq(data_dir, mouse_to_human, target_genes)
    adata_m1 = load_m1_patchseq(data_dir, mouse_to_human, target_genes)

    # Concatenate V1 + M1
    print(f"\n  Concatenating V1 ({adata_v1.n_obs} cells) + M1 ({adata_m1.n_obs} cells)...")

    # Intersect genes
    common_genes = sorted(set(adata_v1.var_names) & set(adata_m1.var_names))
    print(f"  Common genes: {len(common_genes)}")
    adata_v1 = adata_v1[:, common_genes].copy()
    adata_m1 = adata_m1[:, common_genes].copy()

    # Save obs metadata before concat (concat may drop non-shared obs columns)
    obs_v1 = adata_v1.obs.copy()
    obs_m1 = adata_m1.obs.copy()

    adata_query = ad.concat([adata_v1, adata_m1], join="inner")
    adata_query.var_names_make_unique()

    # Restore full obs metadata using pd.concat (keeps all columns, fills NaN)
    obs_combined = pd.concat([obs_v1, obs_m1], axis=0)
    for col in obs_combined.columns:
        if col not in adata_query.obs.columns:
            adata_query.obs[col] = obs_combined.loc[adata_query.obs_names, col].values
    print(f"  Obs columns after merge: {len(adata_query.obs.columns)}")
    del adata_v1, adata_m1
    gc.collect()

    # Ensure counts layer
    adata_query.layers["counts"] = adata_query.X.copy()

    # Flag CCKBC cells
    print(f"\n[4/4] Flagging CCKBC cells...")
    adata_query = flag_cckbc_cells(adata_query, config)

    # Summary
    print(f"\n{'=' * 70}")
    print("MOUSE PATCH-SEQ QUERY SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total cells: {adata_query.n_obs}")
    print(f"  Total genes (human symbols): {adata_query.n_vars}")
    if "dataset" in adata_query.obs.columns:
        print(f"\n  Cells per dataset:")
        for ds, count in adata_query.obs["dataset"].value_counts().items():
            print(f"    {ds}: {count}")
    if "is_cckbc" in adata_query.obs.columns:
        n_cckbc = adata_query.obs["is_cckbc"].sum()
        print(f"\n  CCKBC cells: {n_cckbc} ({100*n_cckbc/len(adata_query):.1f}%)")

    # Save
    out_path = output_dir / "mouse_patchseq_query.h5ad"
    print(f"\n  Saving to {out_path}...")
    adata_query.write(out_path)
    print(f"  Done.")


if __name__ == "__main__":
    main()
