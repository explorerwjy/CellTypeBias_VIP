# Best-of-N Sampling for Gene Matching

## Overview

The Best-of-N sampling approach provides a simple and effective method for generating matched null gene sets. For each null simulation, it:

1. Generates N candidate samples (default: 100)
2. Calculates a similarity metric for each candidate in percentile space
3. Selects the best matching candidate among the N

This approach is **simpler than Sequential Importance Sampling (SIS)** and provides **better control over matching quality** compared to pure rejection sampling.

## Key Advantages

1. **Simplicity**: Easy to understand and implement - no complex sequential selection logic
2. **Predictable matching quality**: Always selects the best match among N candidates
3. **Works in percentile space**: Handles skewed distributions naturally (e.g., CDS_length, LOEUF)
4. **Parallelizable**: Each simulation is independent and can run in parallel
5. **Tunable**: Increase `n_candidates` for better matching (at cost of computation time)

## Usage

### Basic Command

```bash
python scripts/script_generate_geneweights.py \
    --SpecMat dat/ExpMats/HumanCT.TPM.0.1.Filt.csv \
    --WeightDF dat/GeneWeights/DDD.top61.gw.csv \
    --outfile results/null_weights/DDD_bestofn_geneweights.csv \
    --sampling_mode set_level_matched \
    --use_best_of_n \
    --n_candidates 100 \
    --n_sims 10000 \
    --seed 42 \
    --matching_variables CDS_length,WB,LOEUF \
    --n_processes 20
```

### Parameters

- `--use_best_of_n`: Enable best-of-N sampling mode
- `--n_candidates`: Number of candidate samples to evaluate per simulation (default: 100)
  - Higher values = better matching but slower
  - Recommended: 100-500 for good balance
- `--n_sims`: Number of null simulations to generate (default: 10,000)
- `--n_processes`: Number of parallel processes (default: 10)
- `--matching_variables`: Variables to match on (default: `CDS_length,WB,LOEUF`)
- `--seed`: Random seed for reproducibility

### Performance Considerations

**Computational cost**: With best-of-N sampling, you evaluate `n_sims × n_candidates` total samples.
- Example: 10,000 simulations × 100 candidates = 1,000,000 evaluations
- With 20 parallel processes: typically completes in 1-5 minutes for ~60 input genes

**Scaling with n_candidates**:
- 50 candidates: Faster, good matching
- 100 candidates: Balanced (recommended)
- 200 candidates: Better matching, 2× slower
- 500 candidates: Best matching, 5× slower

## Output and Quality Assessment

The script outputs detailed matching quality statistics:

```
============================================================
SET-LEVEL MATCHING QUALITY SUMMARY (Best-of-N Percentile)
============================================================
Distance statistics (lower is better):
  Mean distance: 39.67
  Median distance: 40.00
  Std distance: 1.81
  Max distance: 42.60
  Min distance: 36.11
```

### Interpreting Distance Metrics

The distance metric combines:
1. **Mean difference** in percentile space (0-100 scale)
2. **Std difference** (relative, scaled by 50 for comparability)

**Good matching**:
- Mean distance: < 50 (excellent), < 80 (good), < 120 (acceptable)
- Low std of distances indicates consistent matching quality across simulations

### Validation Output

For each simulation, you'll see both percentile and raw value comparisons:

```
Simulation 0 (distance=39.54):
  PERCENTILES:
    CDS_length_pct: target=72.5, sampled=56.7, diff=-15.9
    WB_pct: target=71.6, sampled=49.3, diff=-22.3
    LOEUF_pct: target=11.0, sampled=38.1, diff=+27.1
  RAW VALUES:
    CDS_length: target=2859.3, sampled=2128.3, error=25.6%
    WB: target=3.7, sampled=2.3, error=38.9%
    LOEUF: target=0.3, sampled=0.7, error=131.6%
```

**Note**: High LOEUF error is expected when input genes are highly constrained (low percentile). The percentile matching is more important than raw value matching for skewed distributions.

## Comparison with Other Methods

| Method | Matching Quality | Speed | Complexity | Sampling Space |
|--------|-----------------|-------|------------|----------------|
| Random | Poor | Fastest | Simplest | Largest |
| Rejection | Variable | Fast | Simple | Large |
| **Best-of-N** | **Good** | **Moderate** | **Simple** | **Large** |
| SIS | Good | Slow | Complex | Moderate |
| Gene-by-Gene | Excellent | Slowest | Moderate | Smallest |

**Recommended**: Use Best-of-N for most analyses. It provides good matching quality with reasonable computation time and simple, interpretable behavior.

## Example Workflow

1. **Generate null gene weights** using best-of-N:
```bash
python scripts/script_generate_geneweights.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --outfile results/null_weights/ASD_bestofn100_geneweights.csv \
    --sampling_mode set_level_matched \
    --use_best_of_n \
    --n_candidates 100 \
    --n_sims 10000 \
    --seed 42 \
    --n_processes 20
```

2. **Compute null bias**:
```bash
python scripts/script_run_ctrl_sim.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --Ctrl_Genes_Fil results/null_weights/ASD_bestofn100_geneweights.csv \
    --outfile results/null_bias/ASD_null_bias.csv \
    --mode human_ct_bias \
    --n_processes 20
```

3. **Add p-values**:
```bash
python scripts/script_bias_cal.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --gw dat/GeneWeights/ASD.gw.csv \
    --Bias_Null results/null_bias/ASD_null_bias.csv \
    --Bias_Out results/ASD_bias_addP.csv
```

## Technical Details

### Distance Calculation in Percentile Space

For each candidate sample, the distance is calculated as:

```python
for each matching variable:
    mean_diff = |sampled_mean_pct - target_mean_pct|
    std_diff = |sampled_std_pct - target_std_pct| / target_std_pct
    variable_distance = mean_diff + std_diff × 50

total_distance = mean(variable_distance across all variables)
```

The scaling factor of 50 for std_diff makes it comparable to mean_diff in percentile space.

### Why Percentile Space?

1. **Uniform scale**: All variables on [0, 100] scale
2. **Handles skewness**: CDS_length has skewness=13.6, but percentiles are uniform
3. **Interpretable distances**: 10 percentile points has the same meaning for all variables
4. **No normalization issues**: Don't need to worry about different units or scales

### Reproducibility

Each simulation gets a unique random seed: `base_seed + sim_idx`
- With `--seed 42`: simulation 0 uses seed 42, simulation 1 uses seed 43, etc.
- This ensures reproducibility while generating diverse null samples
