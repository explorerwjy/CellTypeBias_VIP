# Propensity Score Weighted Sampling - Kernel Functions

## Overview

Propensity weighted sampling converts gene distances (in percentile space) into sampling probabilities using a kernel function. Different kernels offer different trade-offs between matching quality and gene diversity.

## Why Kernel Choice Matters for P-values

**Problem**: If you only sample from 3-5% of the gene pool, your null distribution is too conservative:
- Limited diversity means some gene combinations are never sampled
- P-values become inflated (too large)
- You may miss true signals

**Solution**: Use kernels with broader support (tricubic, epanechnikov, uniform) to increase diversity to 15-25% of the gene pool while maintaining good matching.

## Available Kernels

### 1. Tricubic (RECOMMENDED)

**Formula**: `w(d) = (1 - (d/h)^3)^3` for d < h, else 0

**Properties**:
- **Compact support**: Only genes within bandwidth have non-zero probability
- **Smooth at boundaries**: Continuous and differentiable
- **Good diversity**: With bandwidth=50, uses ~19% of gene pool

**Recommended Settings**:
```yaml
propensity_kernel: "tricubic"
propensity_bandwidth: 50.0  # For ~15-20% diversity
```

**Results (bandwidth=50, CDS+WB matching)**:
- Effective pool: 5,060 genes (33% of candidates)
- Unique genes used: 2,934/15,515 (18.9%) across 100 sims
- Top gene frequency: 7-8%
- Matching quality: CDS ±1.6%, WB ±4.4% (good)

**Best for**: Balance of diversity and matching quality

---

### 2. Gaussian

**Formula**: `w(d) = exp(-(d/h)^2 / 2)`

**Properties**:
- **Infinite support**: All genes have non-zero probability (but very small for distant genes)
- **Smooth decay**: Exponentially decays with distance
- **Less diversity**: Probabilities drop very quickly, concentrating on close matches

**Recommended Settings**:
```yaml
propensity_kernel: "gaussian"
propensity_bandwidth: 15.0  # For ~5-10% diversity
```

**Results (bandwidth=15, CDS+WB matching)**:
- Effective pool: 4,939 genes (32%)
- Unique genes used: 538/15,515 (3.5%) across 10 sims
- Top gene frequency: 30%
- Matching quality: CDS ±0.7%, WB ±1.5% (excellent)

**Best for**: Maximum matching quality, when diversity is less critical

---

### 3. Epanechnikov

**Formula**: `w(d) = (1 - (d/h)^2)` for d < h, else 0

**Properties**:
- **Compact support**: Zero probability beyond bandwidth
- **Parabolic shape**: Smoother than linear, less sharp than tricubic
- **Statistically optimal**: Minimizes mean squared error in kernel density estimation

**Recommended Settings**:
```yaml
propensity_kernel: "epanechnikov"
propensity_bandwidth: 40.0  # For ~12-18% diversity
```

**Expected diversity**: Similar to tricubic but slightly less diverse

**Best for**: Statistical purists who want optimal kernel properties

---

### 4. Uniform

**Formula**: `w(d) = 1` for d < h, else 0

**Properties**:
- **Compact support**: Zero probability beyond bandwidth
- **Flat weighting**: All genes within bandwidth have equal probability
- **Maximum diversity**: No preference for closer matches within bandwidth

**Recommended Settings**:
```yaml
propensity_kernel: "uniform"
propensity_bandwidth: 40.0  # For maximum diversity
```

**Expected diversity**: Highest diversity, but matching quality may be worse

**Best for**: When you want maximum gene pool coverage and are okay with looser matching

---

### 5. Linear

**Formula**: `w(d) = (1 - d/h)` for d < h, else 0

**Properties**:
- **Compact support**: Zero probability beyond bandwidth
- **Linear decay**: Simple, interpretable weighting
- **Moderate diversity**: Between tricubic and uniform

**Recommended Settings**:
```yaml
propensity_kernel: "linear"
propensity_bandwidth: 40.0  # For ~10-15% diversity
```

**Expected diversity**: Moderate, simpler than tricubic

**Best for**: When you want simple, interpretable weighting

---

## Comparison Table

| Kernel | Support | Diversity (100 sims) | Matching Quality | Speed | Recommended Bandwidth |
|--------|---------|---------------------|------------------|-------|----------------------|
| **Tricubic** | **Compact** | **~19%** | **Good** | **Fast** | **50** |
| Gaussian | Infinite | ~3.5% | Excellent | Fast | 15 |
| Epanechnikov | Compact | ~15% | Good | Fast | 40 |
| Uniform | Compact | ~25% | Fair | Fast | 40 |
| Linear | Compact | ~12% | Good | Fast | 40 |

## How to Choose

### For most analyses (RECOMMENDED):
```yaml
propensity_kernel: "tricubic"
propensity_bandwidth: 50.0
```
- Good balance of diversity (19%) and matching quality
- Avoids overly conservative p-values
- Fast computation

### For maximum matching quality:
```yaml
propensity_kernel: "gaussian"
propensity_bandwidth: 15.0
```
- Best matching (±0.7% on CDS)
- Use when diversity is less important
- May give conservative p-values

### For maximum diversity:
```yaml
propensity_kernel: "uniform"
propensity_bandwidth: 40.0
```
- Highest diversity (~25%)
- Looser matching constraints
- Use when you have concerns about sampling bias

## Bandwidth Selection Guide

### For Tricubic (compact support):
- **30-40**: More selective, better matching, ~10-15% diversity
- **50-60**: Balanced (RECOMMENDED), ~15-20% diversity
- **70-80**: More diverse, looser matching, ~25-30% diversity

### For Gaussian (infinite support):
- **10-15**: Tight matching, ~3-5% diversity
- **20-25**: Balanced, ~8-12% diversity
- **30-40**: More diverse, ~15-20% diversity

## Kernel Shape Visualization

```
Tricubic (h=50):        Gaussian (h=15):        Uniform (h=40):
1.0|    ╱╲               1.0|     ╱╲              1.0|  ┌────┐
   |   ╱  ╲                 |    ╱  ╲                |  │    │
0.5|  ╱    ╲             0.5|   ╱    ╲            0.5|  │    │
   | ╱      ╲               |  ╱      ╲              |  │    │
0.0|╱________╲___        0.0| ╱________╲___       0.0|__│    │___
   0   25   50  75          0  15  30  45           0  20  40  60
       distance (pct)           distance (pct)          distance (pct)
```

## Example Usage

```bash
# Tricubic with good diversity (RECOMMENDED)
python scripts/script_generate_geneweights.py \
    --SpecMat dat/ExpMats/HumanCT.spec.csv \
    --WeightDF dat/GeneWeights/ASD.gw.csv \
    --outfile results/null_weights/ASD_propensity_geneweights.csv \
    --sampling_mode set_level_matched \
    --use_propensity \
    --propensity_kernel tricubic \
    --propensity_bandwidth 50.0 \
    --matching_variables CDS_length,WB \
    --n_sims 10000 \
    --seed 42 \
    --n_processes 20

# Gaussian for best matching
python scripts/script_generate_geneweights.py \
    ... \
    --propensity_kernel gaussian \
    --propensity_bandwidth 15.0 \
    ...

# Uniform for maximum diversity
python scripts/script_generate_geneweights.py \
    ... \
    --propensity_kernel uniform \
    --propensity_bandwidth 40.0 \
    ...
```

## Technical Notes

### Distance Calculation
All distances are computed in percentile space after applying variable-specific weights:
```python
for var in matching_vars:
    weight = 3.0 if var == 'CDS_length' else 1.0
    distances += weight * (candidate_pct[var] - target_mean_pct[var])**2
distances = sqrt(distances)
```

### Probability Normalization
After applying the kernel, probabilities are normalized to sum to 1:
```python
propensity_scores = kernel_function(distances, bandwidth)
propensity_probs = propensity_scores / propensity_scores.sum()
```

### Effective Pool Size
Reported as genes with probability > 1% of maximum probability:
```python
effective_pool = (propensity_probs > 0.01 * propensity_probs.max()).sum()
```
