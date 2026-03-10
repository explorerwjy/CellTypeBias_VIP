#!/usr/bin/env python3
"""
Step 4: Train cross-species scVI model on joint human + mouse reference.

Trains scVI jointly on human (Siletti atlas) + mouse (Allen WMB-10Xv3) cells
with batch correction via sample_species key. Both species contribute equally
to the latent space.

Output:
  - results/cross_species/models/scvi_cross_species/ (trained model)
  - results/cross_species/preprocessed/reference_with_latent.h5ad
"""

import sys
import gc
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
import yaml
from scipy import sparse

warnings.filterwarnings("ignore")

# Fix: PyTorch 2.6+ defaults weights_only=True; scvi-tools 1.2.0 predates this
_orig_torch_load = torch.load
torch.load = lambda *a, **kw: _orig_torch_load(*a, **{**kw, "weights_only": False})

SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_MATCHING_DIR = SCRIPT_DIR.parent.parent
PROJECT_ROOT = ATLAS_MATCHING_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(ATLAS_MATCHING_DIR))


def load_config():
    config_path = ATLAS_MATCHING_DIR / "configs" / "cross_species_config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_forced_genes(config):
    """Collect all forced genes from config."""
    forced = set()
    gene_sel = config["gene_selection"]["forced_genes"]

    # Interneuron markers
    for g in gene_sel.get("interneuron_markers", []):
        forced.add(g)

    # Ion channels
    for family_genes in gene_sel.get("ion_channels", {}).values():
        if isinstance(family_genes, list):
            forced.update(family_genes)

    return forced


def select_hvgs_per_species(adata, n_hvgs, species_col="species", seed=42):
    """Select top N HVGs per species using seurat_v3, return union."""
    all_hvgs = set()

    for sp in adata.obs[species_col].unique():
        mask = adata.obs[species_col] == sp
        adata_sp = adata[mask].copy()

        # Subsample for HVG if too large
        if adata_sp.n_obs > 200_000:
            np.random.seed(seed)
            idx = np.random.choice(adata_sp.n_obs, 200_000, replace=False)
            adata_sp = adata_sp[idx].copy()

        adata_sp.X = adata_sp.layers["counts"]
        sc.pp.highly_variable_genes(
            adata_sp, n_top_genes=n_hvgs, flavor="seurat_v3", subset=False
        )
        hvgs = adata_sp.var_names[adata_sp.var["highly_variable"]].tolist()
        all_hvgs.update(hvgs)
        print(f"  {sp}: {len(hvgs)} HVGs (from {adata_sp.n_obs:,} cells)")
        del adata_sp
        gc.collect()

    return all_hvgs


def main():
    parser = argparse.ArgumentParser(description="Train cross-species scVI")
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--model-dir", type=Path, default=None)
    parser.add_argument("--no-gpu", action="store_true")
    parser.add_argument("--force-retrain", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    config = load_config()
    preproc_dir = args.preprocessed_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["preprocessed"])
    model_dir = args.model_dir or (ATLAS_MATCHING_DIR / config["output"]["subdirs"]["models"])
    model_dir.mkdir(parents=True, exist_ok=True)
    use_gpu = not args.no_gpu

    scvi_config = config["scvi"]
    model_name = f"scvi_cross_species_lat{scvi_config['n_latent']}_hid{scvi_config['n_hidden']}"
    model_path = model_dir / model_name

    print("=" * 70)
    print("STEP 4: TRAIN CROSS-SPECIES SCVI MODEL")
    print("=" * 70)
    print(f"  n_latent={scvi_config['n_latent']}, n_hidden={scvi_config['n_hidden']}")
    print(f"  n_layers={scvi_config['n_layers']}, gene_likelihood={scvi_config['gene_likelihood']}")
    print(f"  batch_key={scvi_config['batch_key']}")
    print(f"  GPU: {use_gpu}")

    # Check if model already exists
    if model_path.exists() and (model_path / "model.pt").exists() and not args.force_retrain:
        print(f"\n  Model already exists at {model_path}")
        print(f"  Use --force-retrain to retrain")
        print(f"  Skipping to latent computation...")

        # Load model and compute latent
        from scvi.model import SCVI
        adata_ref = sc.read_h5ad(preproc_dir / "cross_species_reference_hvg.h5ad")
        if "counts" not in adata_ref.layers:
            adata_ref.layers["counts"] = adata_ref.X.copy()

        # Recompute latent and UMAP
        SCVI.setup_anndata(adata_ref, batch_key="_scvi_batch", layer="counts")
        model = SCVI.load(str(model_path), adata=adata_ref)
        adata_ref.obsm["X_latent"] = model.get_latent_representation(adata_ref)

        print("  Computing UMAP on latent space...")
        sc.pp.neighbors(adata_ref, use_rep="X_latent", n_neighbors=30)
        sc.tl.umap(adata_ref)

        out_path = preproc_dir / "reference_with_latent.h5ad"
        adata_ref.write(out_path)
        print(f"  Saved: {out_path}")
        return

    # ================================================================
    # Load preprocessed data
    # ================================================================
    print(f"\n[1/6] Loading preprocessed reference data...")
    human_path = preproc_dir / "human_reference.h5ad"
    mouse_path = preproc_dir / "mouse_reference.h5ad"

    adata_human = sc.read_h5ad(human_path)
    adata_mouse = sc.read_h5ad(mouse_path)
    print(f"  Human: {adata_human.shape}")
    print(f"  Mouse: {adata_mouse.shape}")

    # ================================================================
    # Intersect gene sets
    # ================================================================
    print(f"\n[2/6] Intersecting gene sets...")
    common_genes = sorted(set(adata_human.var_names) & set(adata_mouse.var_names))
    print(f"  Human genes: {adata_human.n_vars}")
    print(f"  Mouse genes: {adata_mouse.n_vars}")
    print(f"  Common genes: {len(common_genes)}")

    adata_human = adata_human[:, common_genes].copy()
    adata_mouse = adata_mouse[:, common_genes].copy()

    # Ensure counts layers
    if "counts" not in adata_human.layers:
        adata_human.layers["counts"] = adata_human.X.copy()
    if "counts" not in adata_mouse.layers:
        adata_mouse.layers["counts"] = adata_mouse.X.copy()

    # ================================================================
    # Concatenate
    # ================================================================
    print(f"\n[3/6] Concatenating human + mouse reference...")
    # Save obs before concat (join="inner" on .var drops non-shared .obs columns)
    obs_human = adata_human.obs.copy()
    obs_mouse = adata_mouse.obs.copy()

    adata_ref = ad.concat([adata_human, adata_mouse], join="inner",
                          index_unique="-")

    # Restore full obs (with species-specific columns; NaN for missing)
    obs_combined = pd.concat([obs_human, obs_mouse], axis=0)
    # Update index to match concat's index_unique suffix
    obs_combined.index = adata_ref.obs_names
    adata_ref.obs = obs_combined
    del obs_human, obs_mouse, obs_combined
    del adata_human, adata_mouse
    gc.collect()

    # Ensure counts layer survived concatenation
    if "counts" not in adata_ref.layers:
        adata_ref.layers["counts"] = adata_ref.X.copy()

    print(f"  Combined reference: {adata_ref.shape}")

    # ================================================================
    # HVG selection
    # ================================================================
    print(f"\n[4/6] Selecting HVGs (per-species union + forced genes)...")
    n_hvgs_per = config["gene_selection"]["n_hvgs_per_species"]

    hvg_union = select_hvgs_per_species(adata_ref, n_hvgs_per, species_col="species", seed=args.seed)

    # Add forced genes
    forced_genes = get_forced_genes(config)
    forced_present = {g for g in forced_genes if g in adata_ref.var_names}
    final_genes = sorted(hvg_union | forced_present)

    print(f"  HVG union: {len(hvg_union)}")
    print(f"  Forced genes (present): {len(forced_present)}")
    print(f"  Final gene set: {len(final_genes)}")

    # Subset to final genes
    adata_ref = adata_ref[:, final_genes].copy()
    adata_ref.layers["counts"] = adata_ref.X.copy()

    # Save HVG-subset reference (for model loading later)
    hvg_ref_path = preproc_dir / "cross_species_reference_hvg.h5ad"
    print(f"  Saving HVG reference: {hvg_ref_path}")
    adata_ref.write(hvg_ref_path)

    # ================================================================
    # Setup and train scVI
    # ================================================================
    print(f"\n[5/6] Training scVI...")
    from scvi.model import SCVI

    # Setup batch key (string categorical)
    batch_key = scvi_config["batch_key"]
    if batch_key not in adata_ref.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in obs columns: {list(adata_ref.obs.columns)}")

    batch_values = [str(x) for x in adata_ref.obs[batch_key].values]
    adata_ref.obs["_scvi_batch"] = pd.Categorical(batch_values)
    n_batches = adata_ref.obs["_scvi_batch"].nunique()
    print(f"  Batch key: {batch_key} ({n_batches} batches)")

    # Verify all batch categories are strings
    cats = adata_ref.obs["_scvi_batch"].cat.categories
    cat_types = set(type(x).__name__ for x in cats)
    if cat_types != {"str"}:
        raise TypeError(f"Batch categories must be strings, found: {cat_types}")

    SCVI.setup_anndata(adata_ref, batch_key="_scvi_batch", layer="counts")

    model = SCVI(
        adata_ref,
        n_latent=scvi_config["n_latent"],
        n_layers=scvi_config["n_layers"],
        n_hidden=scvi_config["n_hidden"],
        gene_likelihood=scvi_config["gene_likelihood"],
    )

    train_kwargs = {
        "max_epochs": scvi_config["max_epochs"],
        "train_size": scvi_config["train_size"],
        "batch_size": scvi_config["batch_size"],
        "early_stopping": scvi_config["early_stopping"],
        "early_stopping_patience": scvi_config["early_stopping_patience"],
        "check_val_every_n_epoch": 5,
        "enable_progress_bar": True,
    }
    if use_gpu:
        train_kwargs["accelerator"] = "gpu"
        train_kwargs["devices"] = "auto"
    else:
        train_kwargs["accelerator"] = "cpu"

    model.train(**train_kwargs)

    # Print training summary
    history = model.history
    if history is not None and hasattr(history, "columns"):
        if "elbo_validation" in history.columns:
            final_train = float(history["elbo_train"].iloc[-1])
            final_val = float(history["elbo_validation"].iloc[-1])
            n_epochs = len(history["elbo_train"])
            print(f"  Final ELBO - Train: {final_train:.2f}, Val: {final_val:.2f}")
            print(f"  Epochs: {n_epochs}")

    # Save model
    print(f"\n  Saving model to {model_path}...")
    model_path.mkdir(parents=True, exist_ok=True)
    model.save(str(model_path), overwrite=True, save_anndata=True)

    # ================================================================
    # Latent representations + UMAP
    # ================================================================
    print(f"\n[6/6] Computing latent representations and UMAP...")
    adata_ref.obsm["X_latent"] = model.get_latent_representation(adata_ref)
    print(f"  Latent shape: {adata_ref.obsm['X_latent'].shape}")

    sc.pp.neighbors(adata_ref, use_rep="X_latent", n_neighbors=30)
    sc.tl.umap(adata_ref)
    print(f"  UMAP computed")

    # Save
    out_path = preproc_dir / "reference_with_latent.h5ad"
    adata_ref.write(out_path)
    print(f"  Saved: {out_path}")

    # Summary
    print(f"\n{'=' * 70}")
    print("CROSS-SPECIES SCVI TRAINING COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Model: {model_path}")
    print(f"  Reference + latent: {out_path}")
    print(f"  Reference shape: {adata_ref.shape}")


if __name__ == "__main__":
    main()
