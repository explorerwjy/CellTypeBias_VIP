# SCZ Gene Set Size Stability Analysis — Design

**Date**: 2026-03-10
**Reviewer Point**: #2 (Gene set size and stability of estimates)
**Notebook**: `notebooks/Number_Gene_Effect.py`

## Problem

Reviewer asks: "In Figure S4, mutation bias correlations remain stable even when the SCZ gene set is expanded to 200 genes. This is surprising: including genes without true association should add noise."

Current reply argues weighting protects against dilution. But data shows SCZ weights are flat (mean 6.8 for added genes vs 8.9 for top-61), so at N=200 the top-61 carry only 37% of total weight. Weighting alone doesn't explain stability.

## Three-Part Answer

### Argument 1: Added SCZ genes carry concordant cell-type signal

**Added SCZ genes (62-200) have correlated bias profiles** with the top-61 (Spearman r=0.60, all CTs; r=0.48, neurons). Random genes show r≈-0.08. MGE interneurons are a top supercluster for both subsets. This means added genes aren't noise — they're weaker but concordant signal.

### Argument 2: Weighting further stabilizes (secondary)

Unweighted analysis drops faster than weighted, confirming mutation-count weighting helps. But this is a secondary factor since SCZ weights are relatively flat.

### Argument 3: Effect size diverges across disorders — justifies N=61

Even though correlation (pattern) is preserved, **effect size magnitude** decreases asymmetrically:

| N genes | SCZ CGE | ASD(ID) CGE | DDD CGE |
|---------|---------|-------------|---------|
| 61      | 0.150   | 0.156       | 0.167   |
| 100     | 0.084   | 0.144       | 0.181   |
| 200     | 0.088   | 0.124       | 0.160   |

At N=61, disorders have comparable effect sizes enabling fair cross-disorder comparison. Beyond 61, SCZ dilutes faster (case-control noise) while ASD/DDD retain signal (strong de novo). Discrimination (std of neuronal effects) also degrades asymmetrically: SCZ std drops from 0.060 to 0.050, while ASD/DDD stay ~1.6× higher.

This justifies selecting N=61: it's the regime where cross-disorder comparisons are fair.

## New Analyses (all in `notebooks/Number_Gene_Effect.py`)

### Analysis 1: Split-half bias comparison (SCZ only)
- Compute bias separately for: top-61, added (62-200), random (139 genes)
- Two-panel scatter: top-61 bias (x) vs added-genes bias (y); top-61 vs random
- Color by neuronal/non-neuronal, annotate Spearman r
- Print top 5 superclusters for each subset

### Analysis 2: Sliding window correlation decay (SCZ only)
- For windows of 61 genes starting at rank k (step=10, k=1..140)
- Y-axis: Spearman r of window-bias vs top-61 reference bias
- Secondary y-axis: mean P-value of genes in window
- Shows signal persists well beyond rank 61

### Analysis 3: Cross-disorder effect size comparison
- For SCZ, ASD(ID), DDD at N = 20, 40, 61, 80, 100, 150, 200
- Plot CGE interneuron mean effect vs N genes (3 colored lines)
- Plot neuronal effect std (discrimination) vs N genes
- Vertical line at N=61 to highlight the chosen threshold
- Shows: comparable at 61, SCZ diverges beyond

## Reply Text Revision
- Lead with concordant signal (Argument 1: r=0.60 vs random r≈0)
- Weighting as amplifying factor (Argument 2)
- Effect size divergence justifies N=61 (Argument 3)
