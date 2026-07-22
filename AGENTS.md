# Codex Project Instructions

These instructions are for Codex agents working in this repository. They are derived from `CLAUDE.md` plus a direct scan of the current codebase. When behavior or parameters matter scientifically, read the implementation and config before answering or editing.

## Repository Purpose

`CellTypeBias_VIP` is the analysis code for the study "Rare mutations implicate CGE interneurons as a vulnerable axis of cognitive deficits across psychiatric disorders." The repository combines single-cell expression specificity matrices with mutation or burden-derived gene weights to estimate disease-associated cell type bias.

Core bias definition:

```text
Bias(cell_type) = sum(weight * specificity) / sum(weight)
```

Positive bias means high-weight genes are enriched for high specificity in that cell type. Main result tables usually include `EFFECT`, `P-value`, `Z-score`, `q-value`, and `-logP`.

## Environment

- Use the `gencic` conda environment before running Python, Snakemake, Jupytext, or pip commands:

```bash
conda activate gencic
```

- Prefer running commands from the repository root:

```bash
cd /home/jw3514/Work/CellType_Psy/CellTypeBias_VIP
```

- Many scripts use absolute project/data paths and some imports load external resources at import time. For example, `src/CellType_PSY.py` reads `/mnt/data0/HumanBrainCellType/annotation.xlsx` when imported. If a command fails because a mounted resource is missing, report that explicitly instead of rewriting paths blindly.

## Worktree Rules

- The working tree may already contain user edits, generated notebooks, figures, logs, and result files. Do not revert or overwrite unrelated changes.
- Treat `dat/`, `results/`, `logs/`, `.snakemake/`, and large symlinked subproject outputs as data/artifacts. Modify them only when the task explicitly requires regeneration or data updates.
- Keep code changes scoped. This is a scientific analysis repo, so small unintended changes to config, random seeds, gene sets, or matching parameters can change manuscript conclusions.
- Do not describe computational methods from memory. Verify against the relevant function, script, and `config/config.yaml`.

## Project Map

- `Snakefile` orchestrates the main mutation-bias pipeline.
- `config/config.yaml` holds production paths, gene sets, expression matrices, random-weight count, process count, and matching configurations.
- `scripts/` contains pipeline and analysis entry points:
  - `script_generate_geneweights.py` samples null gene-weight sets.
  - `script_run_ctrl_sim.py` computes null bias scores.
  - `script_bias_cal.py` computes observed bias, permutation p-values, and supercluster summaries.
  - `build_negctrl_geneweights.py` builds beta-weighted negative-control gene weights.
- `src/CellType_PSY.py` contains the core bias/statistical plotting utilities, including `HumanCT_AvgZ_Weighted()`, `AddPvalue_optimized()`, `AnnotateCTDat()`, `SuperClusterBias_BoxPlot()`, `Z1Conversion_V2()`, and `quantileNormalize()`.
- `src/convergence_utils.py` supports PPI/GO/ephys convergence analyses and has focused unit tests in `tests/test_convergence_utils.py`.
- `notebooks/` contains main and supplementary analyses.
- `notebooks_rebuttal/` contains rebuttal-cycle analyses.
- `cge_subtype/` is a subproject for human CGE subtype / CCKBC analyses. See `cge_subtype/README.md` before editing that area.
- `docs/MS/` contains manuscript files and markdown conversions.
- `docs/plans/` contains algorithm notes and implementation plans.

## Main Pipeline

The Snakemake workflow has three main steps:

1. `scripts/script_generate_geneweights.py`
   Samples null gene sets. Supported modes include `random`, `matched`, and `set_level_matched`.
2. `scripts/script_run_ctrl_sim.py`
   Computes bias for null sets, calling `HumanCT_AvgZ_Weighted()` from `src/CellType_PSY.py`.
3. `scripts/script_bias_cal.py`
   Computes observed bias, permutation p-values, and supercluster-aggregated output.

Common commands:

```bash
conda activate gencic
snakemake -n
snakemake --cores 20
snakemake results/main_results/matched_WB_mean_phastCons_n_CDS_bases_Best1000/Centering/ASD_All_bias_addP.csv --cores 10
```

For script-level use, inspect the matching Snakemake rule and run:

```bash
python scripts/<script_name>.py --help
```

## Current Configuration Notes

Always verify `config/config.yaml` before assuming these are still current.

- `results_base` is `results/main_results`.
- The default enabled analysis type is `Centering`, using `dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip.lowexp.cut1e4.mean_centered.csv`.
- `n_random_weights` is `10000`.
- `n_processes` is `10`.
- Current matching configs include:
  - `random`, writing to `results/main_results/random`
  - `matched`: `set_level_matched`, variables `n_CDS_bases`, `WB`, `mean_phastCons`, `use_best_of_n=true`, `n_candidates=1000`, `use_percentile=true`, `seed=42`, writing to `results/main_results/matched_WB_mean_phastCons_n_CDS_bases_Best1000`
  - `matched_500`: same variables with `n_candidates=500`, writing to `results/main_results/matched_500_WB_mean_phastCons_n_CDS_bases_Best500`
  - `matched_2000`: same variables with `n_candidates=2000`, writing to `results/main_results/matched_2000_WB_mean_phastCons_n_CDS_bases_Best2000`
- Production manuscript descriptions should usually refer to the `matched` best-of-1000 set-level configuration unless the code/config says otherwise.

## Data Formats

- Expression specificity matrices are genes by cell types, usually in `dat/ExpMats/`.
- Gene-weight files live in `dat/GeneWeights/`.
- Gene-weight format is two columns with no header:

```text
Entrez_ID,Weight
```

- Bias results are usually written under auto-generated directories below `results/main_results/<matching_config_and_parameters>/<analysis>/`.

## Negative Controls

The manuscript negative-control traits `NegCtrl_HDL`, `NegCtrl_Alanine`, `NegCtrl_RBC`, and `NegCtrl_IBD` use weights equal to `abs(BETA Burden)` from genebass burden statistics. They preserve the same Entrez IDs as the archived uniform-weighted versions in `dat/GeneWeights/_archive_pre_beta/`. The generating script is `scripts/build_negctrl_geneweights.py`.

Other genebass negative controls such as `NegCtrl_HbA1c`, `NegCtrl_T2D`, `NegCtrl_Parkinson`, and `NegCtrl_Alzheimer` may remain uniform-weighted for reproducibility but are not main-figure negative controls. Check `docs/plans/2026-05-13-negctrl-beta-vs-uniform.md` before changing or describing this logic.

## Notebooks

- Most notebooks are paired with `.py` files using Jupytext.
- Edit the `.py` file first, then sync:

```bash
conda activate gencic
jupytext --sync path/to/notebook_pair.py
```

- Preserve the autoreload/setup cell at the top of paired notebooks.
- When editing figure notebooks, expect generated `.ipynb`, image, and result artifacts to change. Only commit or summarize artifacts that the task actually required.

## Testing

Use focused tests when possible:

```bash
conda activate gencic
pytest tests/test_convergence_utils.py
pytest cge_subtype/tests
```

Some tests are data-dependent or marked slow. If data mounts are unavailable or a test skips because optional data are missing, report that in the final response. Avoid running the full pipeline or all notebooks unless the task explicitly calls for it.

## Manuscript Rules

All manuscript text, Methods descriptions, and response-letter claims must be checked against the actual code and config before finalizing. Confirm exact algorithms, parameters, randomization strategy, weighting, and matching variables.

For manuscript `.docx` to markdown conversion, work in `docs/MS/` and accept tracked changes:

```bash
conda activate gencic
pandoc --track-changes=accept \
       --extract-media=docs/MS/media_<tag> \
       docs/MS/<basename>.docx \
       -o docs/MS/<basename>.md
```

Do not write manuscript conversions to `/tmp`, `/var/tmp`, or another temporary directory. If newer-dated manuscript files replace older versions, archive superseded `.docx`, `.md`, and `media_<tag>/` files under `docs/MS/archive/YYYY-MM-DD/` using the old file's modification date. Do not delete manuscript history except for macOS `._*` metadata or zero-byte files.

## CGE Subtype Subproject

`cge_subtype/` analyzes whether human CGE interneuron subtypes show differential mutation bias comparable to mouse CCKBC/ISI findings. It has its own configs, scripts, tests, notebooks, and symlinked large outputs.

Before editing this area, read:

```bash
cge_subtype/README.md
```

Key tested modules include `cge_subtype/src/concordance.py`, `harmony_mapping.py`, `cluster_correspondence.py`, `ephys_harmonization.py`, and `multimodal_classification.py`.

## Coding Style

- Follow the existing script/notebook style. This repo is not packaged as a conventional installable Python package.
- Existing code often inserts `src/` onto `sys.path`; use the local pattern instead of doing a broad packaging refactor.
- Prefer pandas, NumPy, SciPy, statsmodels, seaborn, matplotlib, Snakemake, and Jupytext already present in the project.
- Avoid changing random seeds, matching variables, output directory names, or gene-set definitions unless the user explicitly asks for an analysis change.
- Add concise comments only where they clarify non-obvious scientific or statistical logic.

## When Answering Questions

- Cite file paths and function names so the user can audit the answer.
- For current methods, read `Snakefile`, `config/config.yaml`, and the relevant `scripts/` or `src/` implementation.
- For figure or manuscript claims, cross-check notebooks, generated result tables, and `docs/MS/Code_Method_Crossref.md` when relevant.
- If an answer depends on generated outputs that are not present or not reproducible in the current environment, say so directly.
