# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

- **Conda env**: `gencic` — activate before running any Python, snakemake, or pip commands

## Overview

CellTypeBias_VIP is a computational pipeline for analyzing cell type-specific mutation bias in psychiatric and neurodevelopmental disorders. The codebase processes single-cell RNA-seq expression data alongside genetic mutation data to identify which cell types are preferentially affected by disease-associated mutations.

**Key Publication**: "Rare mutations implicate CGE interneurons as a vulnerable axis of cognitive deficits across psychiatric disorders"

## Core Workflow

The pipeline follows a three-step workflow managed by Snakemake:

1. **Generate Null Gene Weights** (`script_generate_geneweights.py`)
   - Creates null distributions for statistical testing
   - Supports three sampling modes:
     - `random`: Pure random gene sampling
     - `matched`: Gene-by-gene kernel-weighted matching
     - `set_level_matched`: Distribution-level matching (recommended)
   - Produces CSV files with N simulations (default: 10,000)

2. **Compute Null Bias** (`script_run_ctrl_sim.py`)
   - Calculates cell type bias for each null gene set using multiprocessing
   - Applies the `HumanCT_AvgZ_Weighted` function from `CellType_PSY.py`
   - Generates null distribution of bias scores across 461 cell types

3. **Compute Bias + P-values** (`script_bias_cal.py`)
   - Calculates observed bias for the real gene set
   - Compares observed bias against null distribution (vectorized permutation test)
   - Adds FDR-corrected q-values and z-scores

## Running the Pipeline

### Full Pipeline Execution

```bash
# Activate environment first
conda activate gencic

# Run entire workflow with N cores
snakemake --cores 20

# Dry run to see what will be executed
snakemake -n

# Run specific target
snakemake results/Centering/ASD_All_bias_addP.csv --cores 10
```

### Configuration

Edit `config/config.yaml` before running:
- **Gene sets**: Define paths to gene weight files (ASD, SCZ, DDD, 22q, etc.)
- **Expression matrices**: Specify expression data path (Raw, Centering, QuantileNorm)
- **Matching configs**: Define one or more matching strategies to run

### Individual Script Usage

**Generate null gene weights:**

```bash
# Random sampling
python scripts/script_generate_geneweights.py \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --n_sims 10000 \
    --sampling_mode random \
    --outfile results/random/ASD_geneweights.csv

# Set-level matching with best-of-N (recommended)
python scripts/script_generate_geneweights.py \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --n_sims 10000 \
    --sampling_mode set_level_matched \
    --matching_variables n_CDS_bases,WB,LOEUF \
    --use_best_of_n \
    --n_candidates 1000 \
    --use_percentile \
    --seed 42 \
    --outfile results/set_level_matched/ASD_geneweights.csv
```

**Compute null bias:**
```bash
python scripts/script_run_ctrl_sim.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --Ctrl_Genes_Fil results/random/ASD_geneweights.csv \
    --outfile results/random/ASD_null_bias.csv \
    --n_processes 20
```

**Add p-values to bias results:**
```bash
python scripts/script_bias_cal.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --gw dat/GeneWeights/ASD.gw.csv \
    --Bias_Null results/random/ASD_null_bias.csv \
    --Bias_Out results/ASD_bias_addP.csv
```

## Project Structure

```
CellTypeBias_VIP/
├── config/
│   └── config.yaml              # Snakemake configuration
├── dat/
│   ├── GeneWeights/             # Gene weight files by disorder
│   ├── ExpMats/                 # Expression specificity matrices
│   ├── Variable_2_Match_master_table.csv      # Gene annotations for matching
│   └── annotation.xlsx          # Cell type annotations
├── scripts/
│   ├── script_generate_geneweights.py  # Step 1: Null gene sets
│   ├── script_run_ctrl_sim.py          # Step 2: Null bias calculation
│   └── script_bias_cal.py              # Step 3: Observed bias + p-values
├── src/
│   ├── CellType_PSY.py          # Core bias calculation functions
│   ├── ASD_Circuits.py          # ASD circuit analysis
│   ├── Neurotransmitter.py      # Neurotransmitter annotations
│   └── SA.py                    # Statistical analysis utilities
├── notebooks/                   # Analysis notebooks (see below)
├── docs/                        # Documentation for matching algorithms
├── results/                     # Output directory (auto-generated)
└── Snakefile                    # Workflow definition
```

## Key Concepts

### Expression Specificity Matrix
- **Structure**: Genes (rows) × Cell Types (columns), ~18,000 genes × 461 cell types
- **Values**: Z-scored expression specificity (optionally mean-centered or quantile-normalized)
- **Files**: Located in `dat/ExpMats/`
  - `*.csv` — raw specificity
  - `*.mean_centered.csv` — mean-centered (default)
  - `*.qn.csv` — quantile-normalized

### Gene Weights
- **Format**: Two-column CSV: `Entrez_ID,Weight` (no header)
- **Sources**: De novo mutations (ASD, DDD), case-control excess (SCZ), GWAS scores (UKBB)
- **Location**: `dat/GeneWeights/`

### Bias Calculation
Weighted average specificity per cell type:
```
Bias(cell_type) = Σ(weight × specificity) / Σ(weight)
```

### Matching Modes for Null Distributions

See `docs/MATCHING_MODES.md` for details.

| Mode | Description | Use Case |
|------|-------------|----------|
| `random` | Pure random sampling | Baseline, exploratory |
| `matched` | Gene-by-gene kernel matching | Very restrictive |
| `set_level_matched` | Distribution matching | **Recommended** |

**Set-level matching options:**
- `use_best_of_n`: Sample N candidates, pick best match (recommended)
- `use_propensity`: Propensity-weighted sampling
- Matching variables: `n_CDS_bases`, `WB`, `LOEUF`, `mean_phastCons`

## Analysis Notebooks

| Notebook | Purpose |
|----------|---------|
| `Figures.spec.ipynb` | Reproduce manuscript figures |
| `notebook_preprocessing.ipynb` | Expression data preprocessing |
| `Bias_Mutation_Weights.ipynb` | Generate gene weights from mutation data |
| `notebook_SpecBias.ipynb` | Interactive bias calculation |
| `Similarity_ASD_SCZ.spec.ipynb` | ASD-SCZ bias correlation analysis |
| `Number_Gene_Effect.ipynb` | Gene set size robustness |
| `Phenotype.UKBB.ipynb` | UKBB phenotype analysis |
| `notebook_PBS_FSIQ.ipynb` | IQ-bias correlation |
| `VIP_OtherMarkers.ipynb` | VIP interneuron marker analysis |
| `BiasSignificance.spec.ipynb` | Statistical significance testing |

## Key Functions in `src/CellType_PSY.py`

| Function | Purpose |
|----------|---------|
| `HumanCT_AvgZ_Weighted()` | Core bias calculation |
| `AddPvalue_optimized()` | Vectorized permutation p-values |
| `AnnotateCTDat()` | Add cell type names to results |
| `SuperClusterBias_BoxPlot()` | Visualize bias by supercluster |
| `Z1Conversion_V2()` | Convert expression to Z-scores |
| `quantileNormalize()` | Quantile normalization |

## Output Interpretation

**Key columns in bias results:**
- `EFFECT`: Weighted mean specificity (bias score)
- `P-value`: Empirical permutation p-value
- `Z-score`: Standardized effect size
- `q-value`: FDR-corrected p-value (BH method)
- `-logP`: -log10(p-value) for visualization

**Interpretation:**
- Positive bias = cell type enriched for high-specificity genes with high mutation burden
- Significant: typically q-value < 0.1, Z-score > 2.0

## Common Issues

**"Gene not in Expression dataset"**
- Check gene IDs match (Entrez IDs, string vs int)

**Memory errors during null bias**
- Reduce `--n_sims` or `--n_processes`

**Snakemake can't find files**
- Verify paths in `config.yaml` are absolute
- Run `snakemake -n` to debug

**Matched sampling "No valid matches"**
- Increase bandwidth or use `set_level_matched` mode
- Verify genes exist in master table
