# LOEUF Matching Bug Fix

## Problem Identified

LOEUF matching was performing worse than CDS_length and WB matching. Investigation revealed:

1. **All variables were weighted equally** (weight = 1.0) in the distance calculation
2. **LOEUF often has larger differences** between real and matched genes (200%+ in some cases)
3. **LOEUF was dominating the distance calculation**, making it harder to find good matches on CDS_length and WB

## Root Cause

In `PropensityWeightedGenes`, the distance calculation was:
```python
distances += weight * var_distances  # All weights = 1.0
```

When LOEUF has larger differences, it dominates the Euclidean distance, causing:
- Genes that match well on CDS_length and WB but poorly on LOEUF → High total distance
- Genes that match poorly on CDS_length and WB but well on LOEUF → Lower total distance (if LOEUF differences are smaller)

This led to suboptimal matching where LOEUF differences were prioritized over CDS_length and WB differences.

## Solution

Added **variable-specific weighting** to reduce LOEUF's influence:

1. **Default weights**:
   - `CDS_length`: 1.0 (full weight)
   - `WB`: 1.0 (full weight)
   - `LOEUF`: 0.5 (half weight) ⭐

2. **New parameter**: `--loeuf_weight` (default: 0.5)

3. **Distance calculation now uses**:
   ```python
   weight = variable_weights.get(var, 1.0)
   distances += weight * var_distances
   ```

## Why This Works

- **LOEUF is less informative** for expression bias matching than CDS_length and WB
- **LOEUF often has extreme values** (disease genes have very low LOEUF = high constraint)
- **Reducing LOEUF weight** allows better matching on CDS_length and WB while still considering LOEUF

## Usage

### Command Line
```bash
python scripts/script_generate_geneweights.py \
    --use_propensity \
    --loeuf_weight 0.5 \  # Default: 0.5 (half weight)
    ...
```

### Config File (via Snakefile)
The `loeuf_weight` parameter can be added to `config.yaml`:
```yaml
set_level_matched_sampling:
  loeuf_weight: 0.5  # Weight for LOEUF (default: 0.5)
```

## Expected Improvements

After this fix:
- ✅ Better matching on CDS_length and WB
- ✅ LOEUF still considered but doesn't dominate
- ✅ More balanced distance calculation
- ✅ Better overall matching quality

## Testing Recommendations

1. **Compare matching quality** before and after fix:
   - Check CDS_length_pct differences (should improve)
   - Check WB_pct differences (should improve)
   - Check LOEUF_pct differences (may be slightly worse, but that's OK)

2. **Adjust `loeuf_weight`** if needed:
   - If LOEUF still dominating: decrease to 0.3-0.4
   - If LOEUF not considered enough: increase to 0.6-0.7
   - Default 0.5 is a good starting point

3. **Monitor p-values**: Should improve as matching quality improves

## Technical Details

The fix modifies the distance calculation in `PropensityWeightedGenes`:
- Before: `distances += 1.0 * var_distances` (all equal)
- After: `distances += weight * var_distances` (variable-specific)

This is applied to the squared distance before taking the square root, so it properly scales each variable's contribution to the total Euclidean distance.

