# Gene Matching Modes for Null Distribution Generation

This document explains the three different sampling modes available for generating null gene sets in the CellTypeBias_VIP pipeline.

## Overview

When calculating cell type mutation bias, we need to generate null distributions by creating random gene sets that can serve as controls. However, purely random sampling may not account for confounding variables like gene length, expression level, or evolutionary constraint. Different matching strategies offer different trade-offs between controlling confounders and maintaining sampling flexibility.

## Available Modes

### 1. Random Sampling (`random`)

**How it works:**
- Randomly samples genes from the gene pool without any matching criteria
- Each simulation draws N genes uniformly at random

**Pros:**
- Fastest method
- Maximum sampling space (all genes available)
- Simple and transparent

**Cons:**
- Does not control for confounders (gene length, expression, constraint)
- May lead to biased results if input genes have unusual properties

**When to use:**
- Quick exploratory analysis
- When confounders are not a concern
- As a baseline comparison

**Configuration:**
```yaml
sampling_mode: "random"
```

---

### 2. Gene-by-Gene Matching (`matched`)

**How it works:**
- For **each input gene**, finds a matched gene with similar properties
- Matching is done on ALL selected variables simultaneously (CDS_length, WB, LOEUF)
- Uses kernel-weighted distance in percentile space
- Gene X with (CDS_length=1000, WB=5, LOEUF=0.3) is matched to gene Y with similar values on ALL three dimensions

**Pros:**
- Strongest confounder control
- Each gene has an individually matched counterpart
- Precise matching on selected variables

**Cons:**
- **Very restrictive** - sampling space can be severely limited
- May fail to find matches for genes with unusual combinations of features
- Requires larger bandwidth parameters to find enough candidates
- Slowest of the three methods

**When to use:**
- When you need the strictest confounder control
- When you have abundant genes with diverse property combinations
- When you can afford to increase bandwidth to find matches

**Configuration:**
```yaml
sampling_mode: "matched"
matched_sampling:
  matching_variables: ["CDS_length", "WB", "LOEUF"]
  kernel: "tricubic"
  bandwidth: 100.0
  seed: 42
```

**Limitations:**
If your input genes have:
- Long genes (high CDS_length)
- High brain expression (high WB)
- Low constraint (high LOEUF)

Then gene-by-gene matching requires finding genes that match on ALL three criteria simultaneously, which severely limits the pool.

---

### 3. Set-Level Distribution Matching (`set_level_matched`) - **RECOMMENDED**

**How it works:**
- Instead of matching each gene individually, matches the **overall distribution** of the gene set
- Samples N random genes, then checks if their aggregate properties match the input set
- Uses rejection sampling with mean/std or Kolmogorov-Smirnov distance metrics
- Accepts gene sets where the distributions are "close enough" to the target

**Pros:**
- **Much larger sampling space** - genes don't need individual matches
- Still controls for confounders at the set level
- More robust - easier to find acceptable gene sets
- Faster than gene-by-gene matching
- More realistic null model (preserves overall properties)

**Cons:**
- Weaker individual gene matching than gene-by-gene mode
- May require tuning distance thresholds

**When to use:**
- **Recommended for most applications**
- When you want to control for confounders without over-restricting the sampling space
- When gene-by-gene matching is too restrictive
- When you care about aggregate properties more than individual gene matches

**Configuration (Mean/Std matching):**
```yaml
sampling_mode: "set_level_matched"
set_level_matched_sampling:
  matching_variables: ["CDS_length", "WB", "LOEUF"]
  max_distance: 0.15  # Normalized distance threshold
  max_attempts: 100   # Attempts per simulation
  use_ks_test: false
  seed: 42
```

**Configuration (Kolmogorov-Smirnov test):**
```yaml
sampling_mode: "set_level_matched"
set_level_matched_sampling:
  matching_variables: ["CDS_length", "WB", "LOEUF"]
  use_ks_test: true
  ks_threshold: 0.3  # KS statistic threshold
  max_attempts: 100
  seed: 42
```

---

## Comparison Example

Suppose your input gene set has:
- Mean CDS_length: 1200 bp
- Mean WB expression: 4.5
- Mean LOEUF: 0.35

### Random Mode
- Sampled sets might have mean CDS_length: 800, WB: 3.2, LOEUF: 0.6
- **No guarantee** of matching these properties

### Gene-by-Gene Matched Mode
- Each input gene is individually matched
- Input gene: (CDS=1100, WB=4.2, LOEUF=0.30) → Matched gene: (CDS=1095, WB=4.25, LOEUF=0.32)
- **Problem**: Very few genes might satisfy these constraints for all input genes
- May need bandwidth=100+ to find enough candidates

### Set-Level Matched Mode
- Sampled sets will have similar **aggregate** properties:
  - Mean CDS_length: ~1190-1210 bp
  - Mean WB: ~4.4-4.6
  - Mean LOEUF: ~0.33-0.37
- **Individual genes** don't need to match, but the **overall distribution** does
- Much larger sampling space while still controlling confounders

---

## Tuning Parameters

### For `matched` mode:

**`bandwidth`**: Controls how strict the matching is
- Lower bandwidth (e.g., 10) = stricter matching, fewer candidates
- Higher bandwidth (e.g., 100) = looser matching, more candidates
- **Recommended**: Start with 100, increase if you see warnings about "no valid matches"

**`kernel`**: Weighting function
- `tricubic`: Smooth weighting (recommended)
- `uniform`: Hard threshold

### For `set_level_matched` mode:

**`max_distance`**: Controls how close distributions must be
- Lower (e.g., 0.10) = stricter matching
- Higher (e.g., 0.20) = more lenient
- **Recommended**: 0.15 is a good starting point
- Distance = mean(|mean_diff|/std + |std_diff|/std) across variables

**`max_attempts`**: How many tries per simulation
- More attempts = better chance of finding good match
- **Recommended**: 100 attempts is usually sufficient

**`use_ks_test`**: Use KS test instead of mean/std
- `false`: Uses mean and std differences (faster, recommended)
- `true`: Uses Kolmogorov-Smirnov two-sample test (more rigorous)

**`ks_threshold`**: KS statistic threshold (if use_ks_test=true)
- Lower (e.g., 0.2) = stricter
- Higher (e.g., 0.4) = more lenient
- **Recommended**: 0.3

---

## Choosing Matching Variables

You can match on any subset of:
- **`CDS_length`**: Coding sequence length (controls for mutation rate)
- **`WB`**: Whole brain expression level (controls for expression)
- **`LOEUF`**: Loss-of-function observed/expected upper bound (controls for constraint)

**Recommendations:**
- Use all three `["CDS_length", "WB", "LOEUF"]` for comprehensive matching
- If mainly concerned about constraint: `["LOEUF"]`
- If mainly concerned about expression: `["WB", "LOEUF"]`

---

## Example Use Cases

### Case 1: Exploratory Analysis
```yaml
sampling_mode: "random"
```

### Case 2: Standard Analysis with Confounder Control (Recommended)
```yaml
sampling_mode: "set_level_matched"
set_level_matched_sampling:
  matching_variables: ["CDS_length", "WB", "LOEUF"]
  max_distance: 0.15
  max_attempts: 100
  seed: 42
```

### Case 3: Strict Individual Gene Matching
```yaml
sampling_mode: "matched"
matched_sampling:
  matching_variables: ["CDS_length", "WB", "LOEUF"]
  kernel: "tricubic"
  bandwidth: 100.0
  seed: 42
```

### Case 4: Matching on Constraint Only
```yaml
sampling_mode: "set_level_matched"
set_level_matched_sampling:
  matching_variables: ["LOEUF"]
  max_distance: 0.15
  max_attempts: 100
  seed: 42
```

---

## Validation

All modes print diagnostic information:

### Random mode:
- Reports number of genes sampled

### Matched mode:
- Shows mean candidates per gene
- Reports percentage of non-zero weights
- Displays mean absolute differences for each variable

### Set-level matched mode:
- Shows target distribution statistics
- Reports how many simulations needed relaxed criteria
- Displays distance statistics (mean, median, max)
- Validates first 5 simulations with actual matched values

**Check these outputs to ensure matching quality is acceptable!**

---

## Command Line Usage

You can also run the script directly:

**Random:**
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF data/gene_weights.csv \
    --SpecMat data/expression_matrix.csv \
    --n_sims 10000 \
    --sampling_mode random \
    --outfile results/null_weights.csv
```

**Gene-by-gene matched:**
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF data/gene_weights.csv \
    --SpecMat data/expression_matrix.csv \
    --n_sims 10000 \
    --sampling_mode matched \
    --matching_variables CDS_length,WB,LOEUF \
    --kernel tricubic \
    --bandwidth 100.0 \
    --seed 42 \
    --outfile results/null_weights.csv
```

**Set-level matched:**
```bash
python scripts/script_generate_geneweights.py \
    --WeightDF data/gene_weights.csv \
    --SpecMat data/expression_matrix.csv \
    --n_sims 10000 \
    --sampling_mode set_level_matched \
    --matching_variables CDS_length,WB,LOEUF \
    --max_distance 0.15 \
    --max_attempts 100 \
    --seed 42 \
    --outfile results/null_weights.csv
```

---

## Summary

| Feature | Random | Gene-by-Gene | Set-Level (Recommended) |
|---------|--------|-------------|------------------------|
| **Confounder Control** | None | Strongest | Good |
| **Sampling Space** | Maximum | Very Limited | Large |
| **Speed** | Fastest | Slowest | Fast |
| **Complexity** | Simplest | Most Complex | Moderate |
| **Recommended Use** | Baseline | Strict matching needed | Most applications |

**For most users: Use `set_level_matched` mode.** It provides a good balance between controlling for confounders and maintaining sufficient sampling diversity.
