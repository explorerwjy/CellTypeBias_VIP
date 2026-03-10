#!/usr/bin/env python3
"""
Step 0: Prepare ortholog mapping between mouse and human gene symbols.

Builds a comprehensive mouse_symbol -> human_symbol -> human_ensembl mapping
using multiple sources:
1. Pre-built Mouse2Human_Symbol.pk pickle
2. MGI HOM_MouseHumanSequence.rpt (supplementary 1:1 orthologs)
3. Siletti .var for human_symbol -> human_ensembl

Output: results/cross_species/orthologs/ortholog_mapping.csv
"""

import sys
import pickle
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

# Project paths
SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_MATCHING_DIR = SCRIPT_DIR.parent.parent
PROJECT_ROOT = ATLAS_MATCHING_DIR.parent


def load_config():
    config_path = ATLAS_MATCHING_DIR / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_pickle_mapping(pickle_path):
    """Load pre-built mouse->human symbol mapping from pickle.

    Pickle format: {mouse_symbol: {'Entrez': int, 'humanHomo': [(human_symbol, human_entrez), ...]}}
    Returns: {mouse_symbol: [human_symbol, ...]}
    """
    with open(pickle_path, "rb") as f:
        raw = pickle.load(f)

    mapping = {}
    for mouse_sym, val in raw.items():
        if isinstance(val, dict) and "humanHomo" in val:
            human_syms = [hs for hs, _ in val["humanHomo"]]
            mapping[mouse_sym] = human_syms
        elif isinstance(val, str):
            mapping[mouse_sym] = [val]
        elif isinstance(val, list):
            mapping[mouse_sym] = val
        else:
            # Skip unexpected formats
            continue

    print(f"  Loaded {len(mapping)} mouse->human symbol mappings from pickle")
    # Sample
    sample_keys = list(mapping.keys())[:3]
    for k in sample_keys:
        print(f"    {k} -> {mapping[k]}")
    return mapping


def load_hom_table(hom_path):
    """Load MGI HOM_MouseHumanSequence.rpt for supplementary orthologs."""
    df = pd.read_csv(hom_path, sep="\t", low_memory=False)
    print(f"  Loaded HOM table: {len(df)} rows")

    # Group by DB Class Key to pair mouse/human orthologs
    mouse_rows = df[df["Common Organism Name"] == "mouse, laboratory"].copy()
    human_rows = df[df["Common Organism Name"] == "human"].copy()

    # Build mapping: DB Class Key -> (mouse_symbol, human_symbol)
    mouse_by_key = mouse_rows.groupby("DB Class Key")["Symbol"].apply(list).to_dict()
    human_by_key = human_rows.groupby("DB Class Key")["Symbol"].apply(list).to_dict()

    hom_mapping = {}
    for key in mouse_by_key:
        if key in human_by_key:
            for ms in mouse_by_key[key]:
                for hs in human_by_key[key]:
                    if ms not in hom_mapping:
                        hom_mapping[ms] = []
                    hom_mapping[ms].append(hs)

    # Deduplicate
    for ms in hom_mapping:
        hom_mapping[ms] = list(set(hom_mapping[ms]))

    n_one_to_one = sum(1 for v in hom_mapping.values() if len(v) == 1)
    print(f"  HOM table mappings: {len(hom_mapping)} mouse genes -> human")
    print(f"  1:1 orthologs: {n_one_to_one}")
    return hom_mapping


def build_human_symbol_to_ensembl(h5ad_dir, h5ad_file):
    """Build human_symbol -> human_ensembl from a Siletti h5ad .var."""
    h5ad_path = Path(h5ad_dir) / h5ad_file
    print(f"  Loading .var from {h5ad_path.name}...")
    adata = sc.read_h5ad(h5ad_path, backed="r")
    var_df = adata.var.copy()
    adata.file.close()

    # .var index = Ensembl IDs, .var['Gene'] = symbols
    symbol_to_ensembl = {}
    for ensembl_id, row in var_df.iterrows():
        symbol = row.get("Gene", None)
        if symbol and pd.notna(symbol):
            symbol = str(symbol).strip()
            if symbol and symbol not in symbol_to_ensembl:
                symbol_to_ensembl[symbol] = ensembl_id

    print(f"  Built {len(symbol_to_ensembl)} human symbol->ensembl mappings")
    return symbol_to_ensembl


def main():
    parser = argparse.ArgumentParser(description="Prepare ortholog mapping")
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    config = load_config()
    output_dir = args.output_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["orthologs"])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("STEP 0: PREPARE ORTHOLOG MAPPING")
    print("=" * 70)

    # 1. Load pre-built pickle mapping (mouse symbol -> human symbol)
    print("\n[1/4] Loading pre-built mouse->human symbol pickle...")
    pickle_path = config["orthologs"]["mouse2human_symbol_pickle"]
    pkl_mapping = load_pickle_mapping(pickle_path)

    # 2. Load HOM table for supplementary mappings
    print("\n[2/4] Loading HOM_MouseHumanSequence.rpt...")
    hom_path = config["orthologs"]["hom_table"]
    hom_mapping = load_hom_table(hom_path)

    # 3. Build human symbol -> ensembl from Siletti reference
    print("\n[3/4] Building human symbol->ensembl from Siletti .var...")
    h5ad_dir = config["human_reference"]["h5ad_dir"]
    # Use CGE interneuron file (any will work — they share the same .var)
    h5ad_file = config["human_reference"]["interneuron_files"]["CGE"]["file"]
    symbol_to_ensembl = build_human_symbol_to_ensembl(h5ad_dir, h5ad_file)

    # 4. Compose: mouse_symbol -> human_symbol -> human_ensembl
    print("\n[4/4] Composing full ortholog mapping...")

    rows = []

    # Start with pickle mapping (primary source)
    for mouse_sym, human_sym in pkl_mapping.items():
        if isinstance(human_sym, list):
            for hs in human_sym:
                ensembl = symbol_to_ensembl.get(hs, None)
                rows.append({
                    "mouse_symbol": mouse_sym,
                    "human_symbol": hs,
                    "human_ensembl": ensembl,
                    "mapping_source": "pickle",
                })
        else:
            ensembl = symbol_to_ensembl.get(human_sym, None)
            rows.append({
                "mouse_symbol": mouse_sym,
                "human_symbol": human_sym,
                "human_ensembl": ensembl,
                "mapping_source": "pickle",
            })

    # Add supplementary HOM table mappings (only for genes not already in pickle)
    pkl_mouse_genes = set(pkl_mapping.keys())
    n_added_from_hom = 0
    for mouse_sym, human_syms in hom_mapping.items():
        if mouse_sym not in pkl_mouse_genes:
            for hs in human_syms:
                ensembl = symbol_to_ensembl.get(hs, None)
                rows.append({
                    "mouse_symbol": mouse_sym,
                    "human_symbol": hs,
                    "human_ensembl": ensembl,
                    "mapping_source": "HOM_table",
                })
                n_added_from_hom += 1

    print(f"  Added {n_added_from_hom} supplementary mappings from HOM table")

    # Build DataFrame
    df = pd.DataFrame(rows)

    # Flag 1:1 vs many-to-one
    mouse_counts = df.groupby("mouse_symbol")["human_symbol"].nunique()
    human_counts = df.groupby("human_symbol")["mouse_symbol"].nunique()

    df["is_one_to_one"] = df["mouse_symbol"].map(
        lambda x: mouse_counts.get(x, 0) == 1
    ) & df["human_symbol"].map(
        lambda x: human_counts.get(x, 0) == 1
    )

    # Flag many-to-one (multiple mouse genes -> same human gene)
    df["many_mouse_to_one_human"] = df["human_symbol"].map(
        lambda x: human_counts.get(x, 0) > 1
    )

    # Save
    out_path = output_dir / "ortholog_mapping.csv"
    df.to_csv(out_path, index=False)

    # Report coverage stats
    n_mouse = df["mouse_symbol"].nunique()
    n_human = df["human_symbol"].nunique()
    n_with_ensembl = df["human_ensembl"].notna().sum()
    n_one_to_one = df["is_one_to_one"].sum()
    n_many_to_one = df["many_mouse_to_one_human"].sum()

    print(f"\n{'=' * 70}")
    print("ORTHOLOG MAPPING SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total mapping rows: {len(df)}")
    print(f"  Unique mouse genes: {n_mouse}")
    print(f"  Unique human genes: {n_human}")
    print(f"  With Ensembl ID:    {n_with_ensembl} ({100*n_with_ensembl/len(df):.1f}%)")
    print(f"  1:1 orthologs:      {n_one_to_one}")
    print(f"  Many-mouse-to-one:  {n_many_to_one}")
    print(f"\n  Saved to: {out_path}")
    print(f"  Total human genes in Siletti .var: {len(symbol_to_ensembl)}")


if __name__ == "__main__":
    main()
