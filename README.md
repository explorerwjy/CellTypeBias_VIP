# CellTypeBias_VIP

## Introduction

CellTypeBias_VIP is the official codebase accompanying the paper ["Rare mutations implicate CGE interneurons as a vulnerable axis of cognitive deficits across psychiatric disorders"](https://www.biorxiv.org/content/10.1101/2025.03.28.645799v1.abstract). This toolkit provides scripts and notebooks to process genetics and cell type transcriptomics data, perform cell type mutation bias analysis, and generate figures reported in the paper.

## Dataset Requirements

To run the CellTypeBias_VIP pipeline, you will need the following datasets:

- **Expression Data:**  
  Single-cell RNA-seq expression matrices, preferably in CSV, TSV, or HDF5 format.  
  - Example resource: [Human Brain Cell Atlas](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443)
  - For detailed preprocessing instructions, refer to the Jupyter notebook:  
    [`notebooks/notebook_preprocessing.ipynb`](notebooks/notebook_preprocessing.ipynb)  
    This notebook guides you through cleaning, filtering, and formatting your raw data into suitable expression matrices.

- **Metadata:**  
  Cell annotation files containing cell type labels and sample information.  
  - These files should correspond to your expression matrices and provide necessary cell-level metadata.  
  - See the preprocessing notebook above for guidance on generating or formatting metadata.

- **Genetic Data:**  
  Gene weight files representing mutation burden or other genetic metrics for your gene sets of interest.  
  - For instructions on generating these files, see:  
    [`notebooks/Bias_Mutation_Weights.ipynb`](notebooks/Bias_Mutation_Weights.ipynb)  
    This notebook details how to create gene weight files required for downstream bias analysis.

Please ensure all input files are properly formatted and referenced in your configuration file before running the pipeline.

## Installation

Clone this repository and set up a Python environment with the required dependencies:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/CellTypeBias_VIP.git
   cd CellTypeBias_VIP
   ```

2. **(Recommended) Create and activate a virtual environment:**
   - Using `venv`:
     ```bash
     python3 -m venv env
     source env/bin/activate
     ```
   - Or using `conda`:
     ```bash
     conda create -n celltypebias_vip python=3.9
     conda activate celltypebias_vip
     ```

3. **Install the required dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

This will ensure all necessary packages are installed in an isolated environment.

## Usage Example

### 1. Configure the Pipeline

Before running the analysis, make sure your configuration file (`config/config.yaml`) is properly set up. This file should include:

- Paths to your input expression matrices
- Metadata files with cell type annotations
- Any other required parameters for your analysis

You can use the example `config.yaml` provided in the repository as a template.

**Preparing Input Files**

To generate all necessary input files for the mutation bias pipeline, please follow these Jupyter notebooks:

- [`notebooks/notebook_preprocessing.ipynb`](notebooks/notebook_preprocessing.ipynb): Preprocesses your raw data and prepares expression matrices.
- [`notebooks/Bias_Mutation_Weights.ipynb`](notebooks/Bias_Mutation_Weights.ipynb): Generates gene weight files required for the bias analysis.

These notebooks will guide you through the steps needed to produce all required files for running the pipeline.

### 2. Run Bias Calculations with Snakemake

Once your configuration is ready, you can execute the full bias calculation workflow using [Snakemake](https://snakemake.readthedocs.io/en/stable/):
snakemake --cores N_CORES 

### 3. Reproducing Figures and Analyses from the Manuscript

You can reproduce all main figures from our manuscript using the following Jupyter notebook:

- [`notebooks/Figures.spec.ipynb`](notebooks/Figures.spec.ipynb): Recreates all primary figures presented in the paper.

To explore the effect of gene constraint scores on bias correlations, please refer to:

- [`notebooks/Similarity_ASD_SCZ.spec.ipynb`](notebooks/Similarity_ASD_SCZ.spec.ipynb): Analyzes how gene constraint scores influence bias correlations between ASD and SCZ.

Simply open these notebooks and follow the instructions within to reproduce the analyses and figures.
