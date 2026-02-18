# Multi-Configuration Matching Guide

## Overview

The pipeline now supports running multiple matching strategies in parallel. Each configuration generates results in its own auto-generated directory, making it easy to compare different matching approaches.

## Quick Start

### Single Configuration

```yaml
matching_configs:
  full_model:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]
      use_propensity: true
      # ... other parameters
```

Run: `snakemake --cores 20`

Results: `results/full_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/`

### Multiple Configurations

```yaml
matching_configs:
  phastCons_only:
    sampling_mode: "matched"
    matched_sampling:
      matching_variables: ["mean_phastCons"]
      kernel: "tricubic"
      bandwidth: 100.0

  full_model:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]
      use_propensity: true
```

Run: `snakemake --cores 20`

Results:
- `results/phastCons_only_mean_phastCons_Tricubic/`
- `results/full_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/`

## Common Use Cases

### 1. Compare Different Variable Combinations

Test which variables matter most:

```yaml
matching_configs:
  conservation_only:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["mean_phastCons"]
      use_propensity: true

  expression_only:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["WB"]
      use_propensity: true

  conservation_expression:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["mean_phastCons", "WB"]
      use_propensity: true

  full_model:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]
      use_propensity: true
```

**Output directories:**
- `results/conservation_only_mean_phastCons_PropWeight_Tricubic/`
- `results/expression_only_WB_PropWeight_Tricubic/`
- `results/conservation_expression_WB_mean_phastCons_PropWeight_Tricubic/`
- `results/full_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/`

### 2. Compare Matching Algorithms

Test propensity vs best-of-N vs random:

```yaml
matching_configs:
  random_baseline:
    sampling_mode: "random"

  propensity_matching:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]
      use_propensity: true
      propensity_kernel: "tricubic"
      propensity_bandwidth: 60.0

  best_of_1000:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]
      use_best_of_n: true
      n_candidates: 1000
```

**Output directories:**
- `results/random_baseline/`
- `results/propensity_matching_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/`
- `results/best_of_1000_WB_mean_phastCons_n_CDS_bases_Best1000/`

### 3. Test Different Kernel Bandwidths

Find optimal kernel bandwidth for propensity weighting:

```yaml
matching_configs:
  bw_40:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["WB", "mean_phastCons"]
      use_propensity: true
      propensity_bandwidth: 40.0

  bw_60:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["WB", "mean_phastCons"]
      use_propensity: true
      propensity_bandwidth: 60.0

  bw_80:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["WB", "mean_phastCons"]
      use_propensity: true
      propensity_bandwidth: 80.0
```

**Note:** All configs use same variable set, so directory names will differ only in method details. Add descriptive config names to distinguish them.

### 4. Gene-by-Gene vs Set-Level Matching

Compare restrictive vs flexible matching:

```yaml
matching_configs:
  gene_by_gene:
    sampling_mode: "matched"
    matched_sampling:
      matching_variables: ["WB", "LOEUF"]
      kernel: "tricubic"
      bandwidth: 100.0

  set_level:
    sampling_mode: "set_level_matched"
    set_level_matched_sampling:
      matching_variables: ["WB", "LOEUF"]
      use_propensity: true
      propensity_bandwidth: 60.0
```

**Output directories:**
- `results/gene_by_gene_LOEUF_WB_Tricubic/`
- `results/set_level_LOEUF_WB_PropWeight_Tricubic/`

## Directory Naming Convention

Directories are auto-generated using this pattern:

`results/{config_name}_{sorted_variables}_{method}`

**Examples:**
- Random: `results/{config_name}/`
- Matched: `results/{config_name}_{vars}_{Kernel}/`
  - `results/test_CDS_length_WB_Tricubic/`
- Set-level (propensity): `results/{config_name}_{vars}_PropWeight_{Kernel}/`
  - `results/full_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/`
- Set-level (best-of-N): `results/{config_name}_{vars}_Best{N}/`
  - `results/test_WB_mean_phastCons_Best1000/`
- Set-level (SIS): `results/{config_name}_{vars}_SIS/`
- Set-level (rejection): `results/{config_name}_{vars}_Rejection/`

**Tips:**
- Use descriptive config names to identify results easily
- Variable names are sorted alphabetically
- Config name is always included in the path

## Running the Pipeline

### Dry run (see what will execute)
```bash
snakemake -n
```

### Run all configs
```bash
snakemake --cores 20
```

### Run specific target
```bash
# Run just one gene set with one config
snakemake results/full_model_WB_mean_phastCons_n_CDS_bases_PropWeight_Tricubic/Centering/ASD_All_bias_addP.csv --cores 10
```

### Check job counts
With 1 config, 7 gene sets, 1 analysis type:
- Jobs = 1 + (7 × 3) = **22 jobs** (1 all + 7 geneweights + 7 null_bias + 7 bias_pvalue)

With 2 configs, 7 gene sets, 1 analysis type:
- Jobs = 1 + (2 × 7 × 3) = **43 jobs**

With 4 configs, 7 gene sets, 1 analysis type:
- Jobs = 1 + (4 × 7 × 3) = **85 jobs**

## Complete Configuration Template

Here's a template showing all available options:

```yaml
matching_configs:

  # Your config name (used in output directory)
  my_config_name:

    # Choose sampling mode
    sampling_mode: "set_level_matched"  # Options: random, matched, set_level_matched

    # For sampling_mode = "matched" (gene-by-gene matching)
    matched_sampling:
      matching_variables: ["WB", "LOEUF"]  # Available: CDS_length, WB, LOEUF, mean_phastCons, n_CDS_bases
      kernel: "tricubic"                    # Options: uniform, tricubic
      bandwidth: 100.0                      # Kernel bandwidth in percentile units
      seed: 42                              # Random seed (optional)

    # For sampling_mode = "set_level_matched" (distribution matching)
    set_level_matched_sampling:
      matching_variables: ["n_CDS_bases", "WB", "mean_phastCons"]

      # Choose ONE algorithm:
      use_propensity: true       # Propensity weighted (FASTEST, recommended)
      use_best_of_n: false       # Best-of-N sampling
      use_sis: false             # Sequential importance sampling
      # If all false, uses rejection sampling

      # Propensity weighted parameters
      propensity_kernel: "tricubic"     # Options: gaussian, tricubic, epanechnikov, uniform, linear
      propensity_bandwidth: 60.0        # Kernel bandwidth (40-80 typical)
      add_noise: true                   # Add noise to prevent over-matching
      noise_scale: 10.0                 # Noise scale in percentile units
      relaxed_matching: true            # Use relaxed matching (increases bandwidth)
      loeuf_weight: 2                   # Weight for LOEUF variable (default: 0.5)

      # Best-of-N parameters
      n_candidates: 1000                # Number of candidates per simulation

      # SIS parameters
      temperature: 1.0                  # SIS temperature (0.5-2.0)
      adaptive_temp: true               # Increase temperature during sampling
      use_percentile: true              # Use percentile-based matching

      # Rejection sampling parameters
      max_distance: 0.15                # Max normalized distance
      max_attempts: 100                 # Max attempts per simulation
      use_ks_test: false                # Use KS test instead of mean/std
      ks_threshold: 0.3                 # KS statistic threshold

      seed: 42                          # Random seed (optional)
```

## Tips and Best Practices

1. **Start small**: Test with 1-2 configs before scaling up
2. **Use descriptive names**: Config names like `WB_only` are clearer than `test1`
3. **Check dry run**: Always run `snakemake -n` first to see what will execute
4. **Monitor resources**: More configs = more parallel jobs = more memory/CPU needed
5. **Compare systematically**: Change one variable at a time to isolate effects
6. **Save old results**: Output directories are overwritten if config parameters change

## Troubleshooting

**Problem**: "Could not determine config from path"
- **Solution**: Make sure all config names are unique and don't start with '#'

**Problem**: Too many jobs running at once
- **Solution**: Reduce `--cores` or run fewer configs simultaneously

**Problem**: Can't find results
- **Solution**: Check output directory names with `snakemake -n` to see auto-generated paths

**Problem**: Want to change output directory
- **Solution**: Modify the config name, which changes the directory prefix
