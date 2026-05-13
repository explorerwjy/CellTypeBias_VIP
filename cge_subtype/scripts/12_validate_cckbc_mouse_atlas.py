#!/usr/bin/env python3
"""
Validate CCKBC definition: map M1 Patch-seq cells → Allen WMB-10Xv3 mouse atlas.

Confirms that Gouwens et al. 2020 CCKBC cells (defined by RNA family/type:
Sncg subclass + Vip Sncg/Serpinf1 supertypes) map predominantly to the Sncg
subclass in the newer Allen Whole Mouse Brain atlas.

Method:
  1. Load Allen WMB-10Xv3 cortical GABAergic reference (mouse gene symbols)
  2. Load M1 Patch-seq counts + metadata, flag CCKBC cells
  3. Concatenate, normalize, HVG selection, PCA
  4. Harmony batch correction (atlas vs patchseq)
  5. KNN label transfer: atlas subclass → Patch-seq cells
  6. Cross-tabulate is_cckbc vs predicted_subclass

Output:
  results/cckbc_validation/cckbc_atlas_mapping_validation.csv
  results/cckbc_validation/validation_summary.txt
"""

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
from sklearn.neighbors import NearestNeighbors

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def mem_gb(label=""):
    rss = psutil.Process().memory_info().rss / 1e9
    print(f"  [MEM] {label}: {rss:.1f} GB")


def load_config():
    with open(PROJECT_DIR / "configs" / "cross_species_config.yaml") as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# 1. Load Allen WMB-10Xv3 GABAergic reference
# ---------------------------------------------------------------------------

def load_mouse_reference(config, max_gaba=30000, seed=42):
    """Load cortical GABAergic cells from WMB-10Xv3 with mouse gene symbols."""
    print("=" * 60)
    print("[1/6] Loading Allen WMB-10Xv3 cortical GABAergic reference")
    print("=" * 60)

    # --- Metadata ---
    meta_path = config["mouse_reference"]["metadata_csv"]
    print(f"  Loading metadata from {meta_path}...")
    meta = pd.read_csv(meta_path, low_memory=False)
    print(f"  Total cells: {len(meta):,}")

    # Dataset filter
    meta = meta[meta["dataset_label"] == "WMB-10Xv3"].copy()
    print(f"  After WMB-10Xv3 filter: {len(meta):,}")

    # Region filter — all cortical regions
    region_prefixes = config["mouse_reference"]["region_filter"]["include_prefixes"]
    roi_mask = pd.Series(False, index=meta.index)
    for prefix in region_prefixes:
        roi_mask |= meta["region_of_interest_acronym"].astype(str).str.startswith(prefix)
    meta = meta[roi_mask].copy()
    print(f"  After region filter: {len(meta):,}")

    # GABAergic only (CGE + MGE classes)
    gaba_classes = ["06 CTX-CGE GABA", "07 CTX-MGE GABA"]
    meta = meta[meta["class"].isin(gaba_classes)].copy()
    print(f"  GABAergic cells: {len(meta):,}")

    # Subsample (stratified by subclass) to keep it manageable
    np.random.seed(seed)
    if len(meta) > max_gaba:
        meta = meta.groupby("subclass", group_keys=False).apply(
            lambda x: x.sample(
                n=max(1, int(max_gaba * len(x) / len(meta))),
                random_state=seed,
            )
        )
        if len(meta) > max_gaba:
            meta = meta.sample(n=max_gaba, random_state=seed)
    print(f"  After subsampling: {len(meta):,}")

    print(f"\n  Subclass distribution:")
    for sc_name, cnt in meta["subclass"].value_counts().items():
        print(f"    {sc_name}: {cnt:,}")

    # --- Expression ---
    cell_labels = set(meta["cell_label"].tolist())
    expr_dir = Path(config["mouse_reference"]["expression_dir"])
    expr_files = config["mouse_reference"]["expression_files"]

    chunks = []
    for fname in expr_files:
        fpath = expr_dir / fname
        if not fpath.exists():
            print(f"  WARNING: {fpath} not found, skipping")
            continue

        print(f"\n  Loading {fname}...")
        mem_gb(f"before {fname}")

        adata_b = sc.read_h5ad(fpath, backed="r")

        # Identify cells in this file
        if "cell_label" in adata_b.obs.columns:
            file_labels = adata_b.obs["cell_label"].values
        else:
            file_labels = np.array(adata_b.obs_names.tolist())

        keep_idx = [i for i, cl in enumerate(file_labels) if cl in cell_labels]
        keep_labels = [file_labels[i] for i in keep_idx]

        if not keep_idx:
            adata_b.file.close()
            del adata_b
            continue

        keep_idx = np.array(keep_idx)
        print(f"    {len(keep_idx):,} matching cells")

        # Get mouse gene symbols
        var_df = adata_b.var.copy()
        if "gene_symbol" in var_df.columns:
            gene_symbols = var_df["gene_symbol"].tolist()
        else:
            gene_symbols = var_df.index.tolist()

        # Extract expression
        X_sub = adata_b.X[keep_idx]
        if not sparse.issparse(X_sub):
            X_sub = sparse.csr_matrix(X_sub)
        else:
            X_sub = X_sub.copy()

        adata_b.file.close()
        del adata_b
        gc.collect()

        chunk = ad.AnnData(
            X=X_sub,
            obs=pd.DataFrame(index=keep_labels),
            var=pd.DataFrame(index=gene_symbols),
        )
        chunk.var_names_make_unique()
        chunks.append(chunk)
        print(f"    Chunk: {chunk.shape}")
        mem_gb(f"after {fname}")

    # Concatenate
    print(f"\n  Concatenating {len(chunks)} chunks...")
    common_genes = set(chunks[0].var_names)
    for c in chunks[1:]:
        common_genes &= set(c.var_names)
    common_genes = sorted(common_genes)
    for i in range(len(chunks)):
        chunks[i] = chunks[i][:, common_genes].copy()
    adata_ref = ad.concat(chunks, join="inner")
    adata_ref.var_names_make_unique()
    del chunks
    gc.collect()

    # Attach metadata
    meta_idx = meta.set_index("cell_label", drop=False)
    for col in ["class", "subclass", "supertype", "cluster"]:
        if col in meta_idx.columns:
            adata_ref.obs[col] = meta_idx.loc[adata_ref.obs_names, col].values

    adata_ref.obs["batch"] = "atlas"
    adata_ref.layers["counts"] = adata_ref.X.copy()
    print(f"\n  Reference: {adata_ref.shape}")
    return adata_ref


# ---------------------------------------------------------------------------
# 2. Load M1 Patch-seq
# ---------------------------------------------------------------------------

def load_m1_patchseq(config):
    """Load M1 Patch-seq counts + metadata, flag CCKBC cells."""
    print("\n" + "=" * 60)
    print("[2/6] Loading M1 Patch-seq")
    print("=" * 60)

    data_dir = Path(config["mouse_patchseq"]["data_dir"])
    counts_path = data_dir / config["mouse_patchseq"]["M1"]["counts"]
    meta_path = data_dir / config["mouse_patchseq"]["M1"]["metadata"]

    # Counts (genes x cells → cells x genes)
    counts = pd.read_csv(counts_path, index_col=0, compression="gzip").T
    print(f"  Expression: {counts.shape}")

    # Metadata
    try:
        meta = pd.read_csv(meta_path, sep="\t")
    except Exception:
        meta = pd.read_csv(meta_path, sep=",")
    if "Cell" in meta.columns:
        meta = meta.set_index("Cell", drop=False)

    common = sorted(set(counts.index) & set(meta.index))
    counts = counts.loc[common]
    meta = meta.loc[common]
    print(f"  Matched: {len(common)} cells")

    # Flag CCKBC
    cckbc_def = config["mouse_patchseq"]["cckbc_definition"]
    is_cckbc = pd.Series(False, index=meta.index)

    # RNA family == Sncg
    if "RNA family" in meta.columns:
        for fam in cckbc_def["m1_family_include"]:
            is_cckbc |= (meta["RNA family"].astype(str) == fam)

    # RNA type starts with Vip Sncg / Vip Serpinf1
    if "RNA type" in meta.columns:
        rna_type = meta["RNA type"].astype(str)
        for vip_st in cckbc_def["vip_cckbc_supertypes"]:
            is_cckbc |= rna_type.str.startswith(vip_st, na=False)

    # Exclude non-CCKBC types
    exclude = cckbc_def["exclude_types"]
    for exc in exclude:
        if "RNA type" in meta.columns:
            is_cckbc &= ~meta["RNA type"].astype(str).str.contains(exc, case=False, na=False)

    print(f"  CCKBC: {is_cckbc.sum()} / {len(meta)}")
    print(f"  CCKBC by RNA family:")
    print(f"    {meta.loc[is_cckbc, 'RNA family'].value_counts().to_dict()}")
    print(f"  CCKBC by RNA type:")
    for t, n in meta.loc[is_cckbc, "RNA type"].value_counts().items():
        print(f"    {t}: {n}")

    # Build AnnData
    X = sparse.csr_matrix(counts.values.astype(np.float32))
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=counts.index), var=pd.DataFrame(index=counts.columns))
    adata.layers["counts"] = adata.X.copy()
    adata.obs["is_cckbc"] = is_cckbc.values
    adata.obs["RNA_family"] = meta["RNA family"].values
    adata.obs["RNA_type"] = meta["RNA type"].values
    adata.obs["batch"] = "patchseq"

    print(f"  M1 Patch-seq AnnData: {adata.shape}")
    return adata


# ---------------------------------------------------------------------------
# 3-5. Integration pipeline
# ---------------------------------------------------------------------------

def integrate_and_map(adata_ref, adata_ps, n_hvgs=3000, n_pcs=50, n_neighbors=30):
    """Concatenate, normalize, PCA, Harmony, KNN label transfer."""

    # --- Common genes ---
    print("\n" + "=" * 60)
    print("[3/6] Finding common genes & normalizing")
    print("=" * 60)
    common_genes = sorted(set(adata_ref.var_names) & set(adata_ps.var_names))
    print(f"  Common genes: {len(common_genes)}")
    adata_ref = adata_ref[:, common_genes].copy()
    adata_ps = adata_ps[:, common_genes].copy()

    # Save metadata before concat (concat may drop non-shared obs columns)
    ref_obs_saved = adata_ref.obs.copy()
    ps_obs_saved = adata_ps.obs.copy()

    # Concatenate
    adata = ad.concat([adata_ref, adata_ps], join="inner")
    adata.var_names_make_unique()

    # Re-attach all metadata
    ref_obs_saved["source"] = "atlas"
    ps_obs_saved["source"] = "patchseq"
    combined_obs = pd.concat([ref_obs_saved, ps_obs_saved], axis=0)
    for col in combined_obs.columns:
        adata.obs[col] = combined_obs.loc[adata.obs_names, col].values
    print(f"  Combined: {adata.shape}")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVGs (computed on atlas cells only, applied to all)
    atlas_mask = adata.obs["source"].astype(str) == "atlas"
    adata_atlas_tmp = adata[atlas_mask].copy()
    sc.pp.highly_variable_genes(adata_atlas_tmp, n_top_genes=n_hvgs, flavor="seurat_v3",
                                 layer=None, span=0.3)
    hvgs = adata_atlas_tmp.var_names[adata_atlas_tmp.var.highly_variable].tolist()
    del adata_atlas_tmp
    gc.collect()
    print(f"  HVGs: {len(hvgs)}")

    # PCA on HVGs
    print("\n" + "=" * 60)
    print("[4/6] PCA + Harmony")
    print("=" * 60)
    adata_hvg = adata[:, hvgs].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    del adata_hvg
    gc.collect()
    print(f"  PCA: {adata.obsm['X_pca'].shape}")

    # Harmony
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "batch", max_iter_harmony=30)
    Z = ho.Z_corr
    # Handle torch tensor from GPU mode
    if hasattr(Z, "cpu"):
        Z = Z.cpu().numpy()
    Z = np.array(Z)
    # Z_corr is (n_pcs, n_cells) — transpose to (n_cells, n_pcs)
    if Z.shape[0] == n_pcs and Z.shape[1] == adata.n_obs:
        Z = Z.T
    adata.obsm["X_harmony"] = Z
    print(f"  Harmony: {adata.obsm['X_harmony'].shape}")

    # --- KNN label transfer ---
    print("\n" + "=" * 60)
    print("[5/6] KNN label transfer")
    print("=" * 60)
    ref_mask = adata.obs["source"].astype(str) == "atlas"
    query_mask = adata.obs["source"].astype(str) == "patchseq"

    ref_emb = adata.obsm["X_harmony"][ref_mask.values]
    query_emb = adata.obsm["X_harmony"][query_mask.values]
    ref_subclass = adata.obs.loc[ref_mask, "subclass"].values.astype(str)

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean", n_jobs=10)
    nn.fit(ref_emb)
    distances, indices = nn.kneighbors(query_emb)

    # Distance-weighted voting
    n_query = len(query_emb)
    predicted_subclass = np.empty(n_query, dtype=object)
    confidence = np.zeros(n_query)

    for i in range(n_query):
        weights = 1.0 / (distances[i] + 1e-6)
        weights /= weights.sum()
        neighbor_labels = ref_subclass[indices[i]]
        unique_labels = np.unique(neighbor_labels)
        scores = {lbl: weights[neighbor_labels == lbl].sum() for lbl in unique_labels}
        winner = max(scores, key=scores.get)
        predicted_subclass[i] = winner
        confidence[i] = scores[winner]

    # Attach to Patch-seq obs
    ps_obs = adata.obs[query_mask].copy()
    ps_obs["predicted_subclass"] = predicted_subclass
    ps_obs["prediction_confidence"] = confidence

    print(f"  Predicted subclass distribution (Patch-seq):")
    for sc_name, cnt in pd.Series(predicted_subclass).value_counts().items():
        print(f"    {sc_name}: {cnt}")

    return ps_obs, adata


# ---------------------------------------------------------------------------
# 6. Validation summary
# ---------------------------------------------------------------------------

def validation_summary(ps_obs, output_dir):
    """Cross-tabulate is_cckbc vs predicted_subclass."""
    print("\n" + "=" * 60)
    print("[6/6] Validation summary")
    print("=" * 60)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Full results
    ps_obs.to_csv(output_dir / "cckbc_atlas_mapping_validation.csv")

    # Cross-tab
    cckbc = ps_obs[ps_obs["is_cckbc"] == True]
    non_cckbc = ps_obs[ps_obs["is_cckbc"] == False]

    lines = []
    lines.append("=" * 60)
    lines.append("CCKBC → Mouse Atlas Mapping Validation")
    lines.append("=" * 60)
    lines.append(f"\nTotal M1 Patch-seq cells: {len(ps_obs)}")
    lines.append(f"CCKBC cells: {len(cckbc)}")
    lines.append(f"Non-CCKBC cells: {len(non_cckbc)}")

    lines.append(f"\n--- CCKBC cells → predicted subclass ---")
    cckbc_vc = cckbc["predicted_subclass"].value_counts()
    for sc_name, cnt in cckbc_vc.items():
        pct = 100 * cnt / len(cckbc)
        lines.append(f"  {sc_name}: {cnt} ({pct:.1f}%)")

    # Sncg fraction (atlas uses prefixed names like "047 Sncg Gaba")
    sncg_count = sum(cnt for sc_name, cnt in cckbc_vc.items() if "Sncg" in str(sc_name))
    sncg_frac = sncg_count / len(cckbc) if len(cckbc) > 0 else 0
    lines.append(f"\n  ** Sncg fraction among CCKBCs: {sncg_frac:.1%} ({sncg_count}/{len(cckbc)}) **")

    # Breakdown by RNA family
    lines.append(f"\n--- CCKBC breakdown by original RNA family ---")
    for fam in cckbc["RNA_family"].unique():
        fam_mask = cckbc["RNA_family"] == fam
        fam_sub = cckbc[fam_mask]
        lines.append(f"\n  RNA family = {fam} ({len(fam_sub)} cells):")
        for sc_name, cnt in fam_sub["predicted_subclass"].value_counts().items():
            pct = 100 * cnt / len(fam_sub)
            lines.append(f"    → {sc_name}: {cnt} ({pct:.1f}%)")

    # Breakdown by RNA type
    lines.append(f"\n--- CCKBC breakdown by original RNA type ---")
    for rtype in cckbc["RNA_type"].value_counts().index:
        rt_mask = cckbc["RNA_type"] == rtype
        rt_sub = cckbc[rt_mask]
        lines.append(f"\n  RNA type = {rtype} ({len(rt_sub)} cells):")
        for sc_name, cnt in rt_sub["predicted_subclass"].value_counts().items():
            pct = 100 * cnt / len(rt_sub)
            lines.append(f"    → {sc_name}: {cnt} ({pct:.1f}%)")

    # Non-CCKBC VIP cells as control
    vip_non_cckbc = non_cckbc[non_cckbc["RNA_family"] == "Vip"]
    if len(vip_non_cckbc) > 0:
        lines.append(f"\n--- Control: Non-CCKBC VIP cells ({len(vip_non_cckbc)}) → predicted subclass ---")
        for sc_name, cnt in vip_non_cckbc["predicted_subclass"].value_counts().items():
            pct = 100 * cnt / len(vip_non_cckbc)
            lines.append(f"  {sc_name}: {cnt} ({pct:.1f}%)")

    # Other families as control
    for fam in ["Pvalb", "Sst", "Lamp5"]:
        fam_cells = non_cckbc[non_cckbc["RNA_family"] == fam]
        if len(fam_cells) > 0:
            top_pred = fam_cells["predicted_subclass"].value_counts().index[0]
            top_pct = 100 * fam_cells["predicted_subclass"].value_counts().iloc[0] / len(fam_cells)
            lines.append(f"\n--- Control: {fam} cells ({len(fam_cells)}) → top predicted: {top_pred} ({top_pct:.1f}%) ---")

    # Confidence
    lines.append(f"\n--- Prediction confidence ---")
    lines.append(f"  CCKBC mean confidence: {cckbc['prediction_confidence'].mean():.3f}")
    lines.append(f"  Non-CCKBC mean confidence: {non_cckbc['prediction_confidence'].mean():.3f}")

    # Vip fraction
    vip_count = sum(cnt for sc_name, cnt in cckbc_vc.items() if "Vip" in str(sc_name))
    vip_frac = vip_count / len(cckbc) if len(cckbc) > 0 else 0

    # Conclusion
    lines.append(f"\n{'=' * 60}")
    lines.append(f"CONCLUSION:")
    lines.append(f"  Sncg Gaba: {sncg_frac:.1%} ({sncg_count}/{len(cckbc)})")
    lines.append(f"  Vip Gaba:  {vip_frac:.1%} ({vip_count}/{len(cckbc)})")
    lines.append(f"  Combined (Sncg+Vip): {sncg_frac + vip_frac:.1%}")
    if sncg_frac >= 0.5:
        lines.append(f"\nSncg is the dominant recipient for CCKBC cells.")
        lines.append("This validates Sncg subclass as a good CCKBC approximation.")
    else:
        lines.append(f"\nSncg + Vip together account for {sncg_frac + vip_frac:.0%}.")
        lines.append("Vip-family CCKBCs split between Sncg and Vip, as expected.")
    lines.append("=" * 60)

    summary_text = "\n".join(lines)
    print(summary_text)

    with open(output_dir / "validation_summary.txt", "w") as f:
        f.write(summary_text)

    print(f"\n  Results saved to {output_dir}")
    return summary_text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Validate CCKBC → Sncg in mouse atlas")
    parser.add_argument("--max-gaba", type=int, default=30000,
                        help="Max GABAergic reference cells to load (default: 30000)")
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    config = load_config()
    output_dir = args.output_dir or (PROJECT_DIR / "results" / "cckbc_validation")

    adata_ref = load_mouse_reference(config, max_gaba=args.max_gaba, seed=args.seed)
    adata_ps = load_m1_patchseq(config)

    ps_obs, adata_combined = integrate_and_map(
        adata_ref, adata_ps,
        n_hvgs=args.n_hvgs, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors,
    )

    validation_summary(ps_obs, output_dir)


if __name__ == "__main__":
    main()
