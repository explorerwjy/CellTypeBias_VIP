# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

- **Conda env**: `gencic` — activate before running any Python, snakemake, or pip commands

## Overview

CellTypeBias_VIP is a computational pipeline for analyzing cell type-specific mutation bias in psychiatric and neurodevelopmental disorders. The codebase processes single-cell RNA-seq expression data alongside genetic mutation data to identify which cell types are preferentially affected by disease-associated mutations.

**Key Publication**: "Rare mutations implicate CGE interneurons as a vulnerable axis of cognitive deficits across psychiatric disorders"

## Core Workflow

The pipeline follows a three-step workflow managed by Snakemake:

1. **Generate Random Gene Weights** (`script_generate_geneweights.py`)
   - Creates null distributions for statistical testing
   - Supports two sampling modes:
     - `random`: Pure random gene sampling
     - `matched`: Kernel-weighted matching on CDS_length, WB (whole brain expression), and LOEUF (constraint score)
   - Produces CSV files with N simulations (default: 10,000) of matched/random gene sets

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
# Run entire workflow with N cores
snakemake --cores 20

# Dry run to see what will be executed
snakemake -n

# Run specific target
snakemake results/Centering/ASD_All_bias_addP.csv --cores 10
```

### Configuration

Edit `config/config.yaml` before running:
- **Gene sets**: Define paths to gene weight files for each disorder/phenotype (ASD, SCZ, DDD, UKBB, etc.)
- **Expression matrices**: Specify paths to expression data (supports Raw, Centering, QuantileNorm)
- **Sampling mode**: Choose between `random` or `matched` null gene generation
- **Matched sampling parameters**: Configure kernel type, bandwidth, and matching variables if using matched mode

### Individual Script Usage

**Generate gene weights (null distributions):**

Random mode:
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --n_sims 10000 \
    --sampling_mode random \
    --outfile results/null_weights/ASD_random_geneweights.csv
```

Gene-by-gene matched mode:
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --n_sims 10000 \
    --sampling_mode matched \
    --kernel tricubic \
    --bandwidth 100.0 \
    --matching_variables CDS_length,WB,LOEUF \
    --seed 42 \
    --outfile results/null_weights/ASD_matched_geneweights.csv
```

Set-level matched mode (recommended):
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --n_sims 10000 \
    --sampling_mode set_level_matched \
    --matching_variables CDS_length,WB,LOEUF \
    --max_distance 0.15 \
    --max_attempts 100 \
    --seed 42 \
    --outfile results/null_weights/ASD_setlevel_geneweights.csv
```

**Compute null bias:**
```bash
python scripts/script_run_ctrl_sim.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --Ctrl_Genes_Fil results/null_weights/ASD_random_geneweights.csv \
    --outfile results/null_bias/ASD_null_bias.csv \
    --mode human_ct_bias \
    --n_processes 20
```

**Add p-values to bias results:**
```bash
python scripts/script_bias_cal.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --gw dat/GeneWeights/ASD.gw.csv \
    --Bias_Null results/null_bias/ASD_null_bias.csv \
    --Bias_Out results/ASD_bias_addP.csv
```

## Architecture & Key Concepts

### Expression Specificity Matrix
- **Structure**: Genes (rows) × Cell Types (columns), typically ~18,000 genes × 461 cell types
- **Values**: Z-scored expression specificity, optionally mean-centered or quantile-normalized
- **Format**: CSV with gene Entrez IDs as row index, cell type cluster IDs as column names

### Gene Weights
- **Format**: Two-column CSV: `[Entrez_ID, Weight]` (no header)
- **Weight interpretation**: Mutation burden, effect size, or other gene-level statistic
- **Sources**:
  - De novo mutation counts (ASD, DDD)
  - Case-control excess mutations (SCZ)
  - GWAS-derived polygenic scores (UKBB phenotypes)

### Bias Calculation (Core Algorithm)
The weighted average specificity is computed as:

```
Bias(cell_type) = Σ(weight_gene × specificity_gene,cell_type) / Σ(weight_gene)
```

Implemented in `CellType_PSY.py` as `HumanCT_AvgZ_Weighted()` which:
1. Filters genes to those present in expression matrix
2. Computes weighted mean specificity per cell type
3. Returns DataFrame with EFFECT column (bias score)

### Gene Matching Modes

The pipeline supports three sampling modes for null distribution generation (see `docs/MATCHING_MODES.md` for details):

**1. Random Sampling (`random`)**
- Pure random gene sampling without matching
- Fastest, but doesn't control for confounders
- Use for exploratory analysis or baselines

**2. Gene-by-Gene Matching (`matched`)**
- Matches each input gene individually on all selected variables
- Uses kernel-weighted matching in percentile space
- **Very restrictive** - can severely limit sampling space
- Process:
  1. Convert all variables to percentiles
  2. Compute pairwise Euclidean distances in percentile space
  3. Apply kernel function to get sampling weights
  4. Sample matched genes without replacement (each gene gets individual match)

**3. Set-Level Distribution Matching (`set_level_matched`)** - **RECOMMENDED**
- Matches the overall distribution of the gene set, not individual genes
- Uses rejection sampling with mean/std or Kolmogorov-Smirnov distance
- **Much larger sampling space** while still controlling confounders
- Process:
  1. Calculate target distribution statistics (mean, std) from input genes
  2. For each simulation, sample N random genes
  3. Calculate distribution statistics of sampled set
  4. Accept if distance from target is below threshold, otherwise retry
  5. Track matching quality across simulations

**Key difference:** Gene-by-gene matching requires finding similar genes for each input gene on ALL variables simultaneously, which is very restrictive. Set-level matching only requires the sampled set to have similar aggregate properties, allowing much more flexibility.

### Statistical Testing
- **Permutation-based p-values**: Vectorized comparison of observed bias vs. null distribution
- **FDR correction**: Benjamini-Hochberg procedure (α = 0.1)
- **Z-scores**: Standardized effect sizes relative to null mean/SD

## Project Structure

```
CellTypeBias_VIP/
├── config/
│   └── config.yaml              # Main configuration for Snakemake
├── dat/
│   ├── GeneWeights/             # Input gene weight files for each disorder
│   ├── ExpMats/                 # Expression specificity matrices
│   └── Variable_2_Match_master_table.csv  # Gene annotations for matching
├── scripts/
│   ├── script_generate_geneweights.py     # Step 1: Null gene sets
│   ├── script_run_ctrl_sim.py             # Step 2: Null bias calculation
│   └── script_bias_cal.py                 # Step 3: Observed bias + p-values
├── src/
│   └── CellType_PSY.py          # Core functions for bias calculation
├── notebooks/
│   ├── notebook_preprocessing.ipynb       # Expression data preprocessing
│   ├── Bias_Mutation_Weights.ipynb        # Generate gene weights
│   ├── Figures.spec.ipynb                 # Reproduce paper figures
│   └── Similarity_ASD_SCZ.spec.ipynb      # Constraint score analysis
├── results/                     # Output directory (generated by pipeline)
└── Snakefile                    # Workflow definition
```

## Important Implementation Details

### Path Configuration
Many scripts contain hardcoded project directory paths (e.g., `/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/`). When modifying scripts:
- Update `ProjDIR` variable at the top of scripts
- Ensure paths in `config.yaml` are absolute
- Check annotation file path in `src/CellType_PSY.py` (line 17: `/mnt/data0/HumanBrainCellAtlas/annotation.xlsx`)

### Cell Type Annotations
The global `Anno` DataFrame (loaded in `CellType_PSY.py`) provides:
- Cluster ID → Cell type name mapping
- Neurotransmitter annotations
- Supercluster groupings (used for aggregating results)

Access via `AnnotateCTDat(bias_df, Anno)` to add readable labels to results.

### Multiprocessing
- `script_run_ctrl_sim.py` uses `multiprocessing.Pool` (default: 20 processes)
- Ensure sufficient memory when running large simulations (10,000 iterations × 18,000 genes)
- Control parallelism via `--n_processes` flag or Snakemake `--cores`

### File Format Expectations

**Gene weights input:**
```
1234,0.85
5678,1.23
9012,0.42
```

**Expression matrix:**
```
gene_id,0,1,2,...,460
1234,0.45,-0.23,0.12,...,-0.34
5678,-0.12,1.23,-0.45,...,0.67
```

**Bias output (with p-values):**
```
cluster_id,EFFECT,P-value,Z-score,q-value,-logP,Cell_Type_Name,...
0,1.234,0.001,3.45,0.015,3.0,Excitatory_Neuron_L2/3
1,0.567,0.023,2.12,0.089,1.64,Inhibitory_Neuron_VIP
```

## Data Preprocessing

Before running the pipeline, prepare input files using the notebooks:

1. **Expression matrices** (`notebooks/notebook_preprocessing.ipynb`):
   - Load single-cell data (CSV, TSV, or HDF5)
   - Filter low-quality cells and genes
   - Calculate specificity scores (z-scores within each cell type)
   - Apply centering or quantile normalization if desired

2. **Gene weights** (`notebooks/Bias_Mutation_Weights.ipynb`):
   - Aggregate mutation counts (LGD, missense) with weights
   - Apply background mutation rate corrections
   - Filter to genes present in expression matrix
   - Export as two-column CSV (Entrez ID, weight)

## Simulation Framework

The repository includes specifications for forward genetic simulations (`Genetic_simulation_pipeline_spec.md`) to validate the bias metric:

**Key concepts:**
- Define ground truth causal cell type (e.g., CGE interneurons)
- Simulate de novo or case-control mutation data with known effect sizes
- Test whether the pipeline correctly recovers the causal cell type
- Evaluate stability across sample sizes and genetic architectures (ASD-like, SCZ-like, NDD-like)

**Simulation parameters:**
- Sample sizes: 1K - 100K
- Relative risk distributions: Gamma-distributed, fitted from real data
- Causal gene selection: Top 500 genes by cell type specificity
- Null models: Random or matched gene sampling

## Python Environment

The repository uses conda/pip with extensive dependencies (see `requirements.txt`). Key packages:
- **Data**: pandas, numpy, scipy
- **Statistics**: statsmodels, scikit-learn, pingouin
- **Visualization**: matplotlib, seaborn, plotly
- **Single-cell**: scanpy, anndata
- **Workflow**: snakemake (install separately if not in requirements)

Install dependencies:
```bash
pip install -r requirements.txt
```

Note: Some packages have pinned versions from conda builds. If encountering conflicts, use conda:
```bash
conda env create -f environment.yml  # If available
# OR manually install core packages:
conda install pandas numpy scipy statsmodels matplotlib seaborn scanpy snakemake -c conda-forge
```

## Common Issues & Solutions

**Issue: "Gene not in Expression dataset"**
- Ensure gene weight file uses Entrez IDs matching expression matrix row names
- Check for type mismatches (string vs integer IDs)

**Issue: Memory errors during null bias calculation**
- Reduce `--n_sims` parameter (e.g., 1000 instead of 10,000)
- Decrease `--n_processes` to reduce memory footprint
- Use smaller expression matrix (filter to fewer genes/cell types)

**Issue: Snakemake fails to find input files**
- Verify all paths in `config.yaml` are absolute
- Check that `ProjDIR` in config matches actual project location
- Run `snakemake -n` to debug missing files before execution

**Issue: Matched sampling produces warnings about "No valid matches"**
- Increase `--bandwidth` parameter (e.g., from 10 to 50 or 100)
- Check that `Variable_2_Match_master_table.csv` contains matching variables
- Verify genes in weight file exist in master table

## Analysis Notebooks

The `notebooks/` directory contains Jupyter notebooks for analysis and visualization:

- **Figures.spec.ipynb**: Reproduces all main figures from the manuscript
- **Similarity_ASD_SCZ.spec.ipynb**: Analyzes constraint score effects on bias correlation
- **Number_Gene_Effect.ipynb**: Tests robustness to gene set size
- **notebook_SpecBias.ipynb**: Interactive bias calculation and visualization

These notebooks import functions from `src/CellType_PSY.py` and expect data in `dat/` and `results/` directories.

## Code Style Notes

- Functions in `CellType_PSY.py` use uppercase names (e.g., `HumanCT_AvgZ_Weighted`)
- Gene dictionaries use format: `{entrez_id: weight}`
- DataFrame operations are vectorized where possible (see `AddPvalue_vec` function)
- Statistical functions return DataFrames with standardized column names: `EFFECT`, `P-value`, `Z-score`, `q-value`

## Output Interpretation

**Key output columns:**
- `EFFECT`: Weighted mean specificity (bias score)
- `P-value`: Empirical p-value from permutation test
- `Z-score`: Standardized effect size
- `q-value`: FDR-corrected p-value
- `-logP`: -log10(p-value) for visualization

**Positive bias** indicates cell type is enriched for high-specificity genes with high mutation burden.

**Significant hits**: Typically q-value < 0.1, Z-score > 2.0
