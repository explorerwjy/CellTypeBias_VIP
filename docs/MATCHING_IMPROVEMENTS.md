# Matching Framework Improvements

## Problem Statement

The original matching framework was **over-matching**, which caused:
- Null distributions that were too similar to real distributions
- Loss of statistical power (p-values becoming non-significant)
- Inability to detect real biological differences

This is a classic problem in statistical matching: you want to control for confounders but NOT remove the biological signal you're trying to detect.

## Solution: Relaxed Matching Options

Three new parameters have been added to prevent over-matching:

### 1. `add_noise` (boolean)
- **Purpose**: Add random noise to the target distribution to prevent perfect matching
- **How it works**: Instead of matching exactly to the mean, adds Gaussian noise (scale controlled by `noise_scale`)
- **When to use**: When you see matching is too perfect (CDS_length differences < 1%)
- **Default**: `false`

### 2. `noise_scale` (float)
- **Purpose**: Controls the magnitude of noise added (in percentile units)
- **Range**: Typically 3-10 percentile units
- **Larger values**: More variance preserved = better statistical power
- **Default**: `5.0`

### 3. `relaxed_matching` (boolean) ⭐ **RECOMMENDED**
- **Purpose**: Combines noise addition with bandwidth increase for maximum variance preservation
- **How it works**: 
  - Adds noise to target distribution (`noise_scale` percentile units)
  - Increases effective bandwidth by 50% (more diverse gene selection)
- **When to use**: When p-values are being killed by over-matching
- **Default**: `false`

## Usage Examples

### Example 1: Enable Relaxed Matching (Recommended)
```yaml
set_level_matched_sampling:
  use_propensity: true
  propensity_kernel: "tricubic"
  propensity_bandwidth: 50.0
  relaxed_matching: true  # Enable relaxed matching
  noise_scale: 5.0        # 5 percentile units of noise
```

### Example 2: Add Noise Only
```yaml
set_level_matched_sampling:
  use_propensity: true
  propensity_kernel: "tricubic"
  propensity_bandwidth: 50.0
  add_noise: true         # Add noise but keep bandwidth same
  noise_scale: 5.0
```

### Example 3: Increase Bandwidth (Alternative)
```yaml
set_level_matched_sampling:
  use_propensity: true
  propensity_kernel: "tricubic"
  propensity_bandwidth: 75.0  # Increase from 50 to 75 (50% increase)
  # This achieves similar effect to relaxed_matching but manually
```

## How to Diagnose Over-Matching

Check your matching quality summary output. Signs of over-matching:

1. **CDS_length_pct differences < 1%**: Too perfect
2. **All matched sets have nearly identical means**: No variance
3. **P-values become non-significant**: Statistical power lost

### Expected Matching Quality

Good matching (preserves power):
- CDS_length_pct: 2-5% difference
- WB_pct: 3-8% difference
- Some variance across matched sets

Over-matching (kills power):
- CDS_length_pct: < 1% difference
- WB_pct: < 2% difference
- All matched sets nearly identical

## Technical Details

### How Noise is Applied

When `add_noise=True` or `relaxed_matching=True`:
```python
target_mean_adjusted = target_mean + np.random.normal(0, noise_scale)
target_mean_adjusted = np.clip(target_mean_adjusted, 0, 100)  # Clamp to [0, 100]
```

### How Relaxed Matching Works

1. **Noise addition**: Adds `noise_scale` percentile units of Gaussian noise to target mean
2. **Bandwidth increase**: Multiplies bandwidth by 1.5x (50% increase)
3. **Result**: More diverse gene selection with preserved variance

### Why This Preserves Statistical Power

- **Variance preservation**: Noise and wider bandwidth allow matched sets to have natural variation
- **Prevents perfect matching**: Null distribution won't be identical to real distribution
- **Maintains confounder control**: Still matches on CDS_length and WB, just not perfectly
- **Biological signal preserved**: Real differences between disease and control genes remain detectable

## Recommendations

1. **Start with `relaxed_matching: true`** if you're seeing over-matching issues
2. **Adjust `noise_scale`** based on your matching quality:
   - If still too tight: increase to 7-10
   - If too loose: decrease to 3-5
3. **Monitor matching quality** in your output logs
4. **Compare p-values** before and after enabling relaxed matching

## Migration Guide

If you're currently experiencing over-matching:

1. **Enable relaxed matching**:
   ```yaml
   relaxed_matching: true
   noise_scale: 5.0
   ```

2. **Run your analysis** and check matching quality

3. **Adjust if needed**:
   - Still too tight? Increase `noise_scale` to 7-10
   - Too loose? Decrease `noise_scale` to 3-5

4. **Compare results**: Check if p-values improve (become more significant)

## References

- Propensity score matching in causal inference
- Variance preservation in statistical matching
- Over-matching problem in observational studies

