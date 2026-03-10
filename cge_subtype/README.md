# CGE Subtype Analysis: CCKBC vs ISI in Human Data

**Question**: Do human CGE interneuron subtypes show differential mutation bias, paralleling the mouse finding that CCK basket cells (CCKBCs) are less affected than irregular-spiking (ISI) VIP+ interneurons?

**Reviewer context**: Both Reviewer #2 and #3 ask whether the mouse VIP/CCKBC distinction holds in human transcriptomic data.

## Status & Key Finding

**Previous result**: No convergent CCKBC definition across three approaches:
- **Marker-based** (CR-, SNCG+, CCK+, M2R-): clusters 276, 279
- **Ephys-based** (fast-spiking, low UD ratio via scANVI-mapped human patch-seq): clusters 284, 289
- **Mouse mapping** (cross-species scVI, Sncg subclass → human): clusters 279, 280, 281 (high fraction)

Mouse SNCG markers do NOT predict human CCKBC identity — approaches disagree on which human clusters are CCKBC homologs.

## Directory Structure

```
cge_subtype/
├── notebooks/
│   ├── CGE_Subtype_22q_Bias.py       # Main analysis: marker classification + multi-disorder bias
│   ├── CCKBC_Human_Mapping.py        # Cross-species mouse→human CCKBC mapping
│   ├── cge_cckbc_ephys_analysis.py   # Human patch-seq ephys validation
│   ├── cge_ephys_overview.py         # Ephys overview across CGE clusters
│   ├── AllenHumanCGE.py              # Siletti atlas CGE marker extraction
│   └── AllenMouseCGE.py              # Allen WMB-10Xv3 CGE marker extraction
├── scripts/
│   ├── 00_prepare_orthologs.py       # Mouse↔human gene ortholog mapping
│   ├── 01_preprocess_human_reference.py
│   ├── 01h_preprocess_human_interneuron.py
│   ├── 02_preprocess_mouse_reference.py
│   ├── 02h_train_and_map_human_scanvi.py  # scANVI for human patch-seq → atlas
│   ├── 03_preprocess_mouse_patchseq.py
│   ├── 04_train_cross_species_scvi.py     # Joint human+mouse scVI
│   └── 05_map_mouse_patchseq.py           # Map mouse CCKBC → human clusters
├── configs/
│   ├── cross_species_config.yaml     # Main config (atlas paths, scVI params, CCKBC def)
│   ├── human_interneuron_config.yaml
│   ├── scvi_model_config.yaml
│   └── species_config.yaml
├── src/
│   └── scvi_integration.py           # scVI train/map helper functions
├── dat/
│   └── LeeDalley_ephys_fx.csv        # Human patch-seq ephys features
└── results/
    ├── convergent_cckbc_classification.csv  # Key output: per-cluster convergent table
    ├── human_scvi_mapping_results.csv       # Human patch-seq → atlas mapping
    ├── patchseq_mapping_results.csv         # Mouse patch-seq → human mapping
    ├── human_scvi_reference_obs.csv → ...   # (symlink) large reference metadata
    ├── cross_species_models → ...           # (symlink) trained scVI models (2.3G)
    ├── human_interneuron_models → ...       # (symlink) scANVI models (31G)
    ├── cross_species_preprocessed → ...     # (symlink) preprocessed h5ad
    ├── human_interneuron_preprocessed → ... # (symlink) preprocessed h5ad
    └── orthologs → ...                      # (symlink) ortholog tables
```

## Provenance

Code migrated from `/home/jw3514/Work/NeurSim/TransEphys/atlas_matching/` on 2026-03-10.
Large data (models, preprocessed h5ad) are symlinked to original locations.
Small result CSVs and ephys features are copied.
