# CellTypeBias_VIP

## Introduction

CellTypeBias_VIP is the official codebase accompanying the paper ["Rare mutations implicate CGE interneurons as a vulnerable axis of cognitive deficits across psychiatric disorders"](https://www.biorxiv.org/content/10.1101/2025.03.28.645799v1.abstract). This toolkit provides scripts and notebooks to process genetics and cell type transcriptomics data, perform cell type mutation bias analysis, and generate figures reported in the paper.

## Dataset Requirements

- **Input Data:** Single-cell RNA-seq expression matrices (e.g., in CSV, TSV, or HDF5 format).
- **Metadata:** Cell annotation files containing cell type labels and sample information.
- **Recommended:** Preprocessed and quality-controlled datasets.

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

Before running the analysis, ensure that your configuration file (`config.yaml`) is correctly set up. This file should specify the paths to your input expression matrices, metadata files, and any other required parameters. Refer to the example `config.yaml` provided in the repository for guidance.

### 2. Run Bias Calculations with Snakemake

Once your configuration is ready, you can execute the full bias calculation workflow using [Snakemake](https://snakemake.readthedocs.io/en/stable/):
