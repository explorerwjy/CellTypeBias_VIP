# Specificity Inflation Simulation: Technical Specification

This document describes the negative binomial (NB) simulation used to demonstrate that
extreme specificity values in the Siletti et al. brain atlas arise from sampling noise,
not genuine expression enrichment. The simulation is implemented in
`notebooks_rebuttal/FigS4_Specificity_Cap.py` (Section 3) and originally in
`notebooks_rebuttal/Specificity_ZINB_Simulation.py`.

---

## 1. Overview

**Goal**: Show that genes with uniform true expression across all 461 cell types
(true specificity = 1.0 everywhere) produce artifactually inflated specificity in
low-UMI cell types, purely due to sampling noise.

**Approach**:
1. Fit an NB noise model from real single-cell data
2. Simulate gene counts using that noise model + real atlas parameters
3. Compute specificity using the same formula as the real pipeline
4. Show that low-UMI cell types produce specificity >> 2 for lowly expressed genes

---

## 2. Noise Model Fitting

### 2.1 Data source

Raw single-cell count matrices per cluster from Siletti et al. (2023), stored as
`/mnt/data0/HumanBrainCellType/cluster_GeneXCell/{cluster_id}.csv.gz`.

Each file: genes (rows) × cells (columns), raw UMI counts.

### 2.2 Representative clusters

Three clusters spanning the UMI depth range:

| Cluster ID | Name | N cells | Total UMI/cell | Fitted θ | Median Var/Mean |
|-----------|------|---------|----------------|----------|-----------------|
| 0 | Miscellaneous (B-cell) | 105 | 2,259 | 1.186 | 1.040 |
| 44 | Oligodendrocyte | 101,039 | 5,994 | 1.938 | 1.179 |
| 200 | DG granule neuron | 23,293 | 17,663 | 0.789 | 1.234 |

### 2.3 Fitting procedure

For each cluster:

1. **Load raw counts**: Read up to 2,000 cells (columns) per cluster for efficiency.
2. **Compute per-gene mean and variance**: `gene_mean = df.mean(axis=1)`, `gene_var = df.var(axis=1, ddof=1)`. Filter to genes with mean > 0.
3. **Bin by mean expression**: Divide genes into 30 log-spaced bins by `log10(mean)`. Compute the median of mean and variance within each bin (requires ≥10 genes per bin).
4. **Fit NB dispersion θ**: Fit the NB variance model to the binned medians:

   ```
   Var(X) = μ + μ² / θ
   ```

   using `scipy.optimize.curve_fit` with initial guess θ = 1.0 and bounds [0.01, 1000].

   When θ → ∞, the NB reduces to Poisson (Var = μ). All clusters show Var/Mean > 1, confirming NB (not Poisson) is the appropriate model.

### 2.4 Global θ

The global θ used for simulation is the **median** of the three fitted θ values:

```
θ_global = median(1.186, 1.938, 0.789) = 1.186
```

**Interpretation**: θ ≈ 1.2 means substantial overdispersion — the variance is roughly
μ + μ²/1.2 ≈ μ + 0.83μ², so for a gene with mean count 10, variance ≈ 10 + 83 = 93,
compared to Poisson variance of 10.

---

## 3. Expression Level Sweep

### 3.1 Determining expression levels to simulate

Expression levels are drawn from the **real** gene expression distribution:

1. **Load cluster-average expression matrix**: `Human.CT.Exp.Entrez.csv` (genes × 461 cell types, raw expression values — not TPM, not specificity).
2. **Compute per-gene mean expression**: Average across all 461 cell types.
3. **Compute mean cluster total**: `mean_cluster_total = mean(sum of expression per cell type)`.
4. **Convert to expression fraction**: For each percentile of the gene-mean distribution:

   ```
   true_frac = percentile_value(gene_means_nonzero, P) / mean_cluster_total
   ```

   This gives the fraction of a cell's total UMI expected to come from this gene.

### 3.2 Percentiles simulated

| Percentile | true_frac | Interpretation |
|-----------|-----------|----------------|
| P10 | 7.3e-08 | Very lowly expressed gene (bottom 10%) |
| P25 | 1.3e-06 | Low expression |
| P50 | 1.2e-05 | Median expression |
| P75 | 4.3e-05 | Moderately high expression |
| P90 | 1.1e-04 | Highly expressed |
| P95 | 2.0e-04 | Very highly expressed (top 5%) |

---

## 4. Simulation Procedure

### 4.1 Function: `simulate_gene_specificity()`

**Parameters**:
- `anno_df`: Cell type annotation DataFrame (461 rows, with "Total UMI" and "Number of cells")
- `true_frac`: Gene expression as fraction of total UMI (**uniform across all cell types** — this is the key: true specificity = 1.0)
- `n_sims`: Number of simulated genes (default: 500)
- `theta`: NB overdispersion parameter (default: 1.186)
- `n_cells`: Number of cells to simulate per cell type (default: 500)
- `min_tpm`: TPM floor — values below this are set to 0 (default: 0.1)
- `seed`: Random seed (default: 42)

### 4.2 Step-by-step algorithm

For each of the 500 simulated genes, independently and in parallel across all 461 cell types:

#### Step 1: Sample per-cell library sizes

For each cell type `ct` with observed total UMI = `T_ct`:

```
cell_total_UMI[i] ~ LogNormal(mean=log(T_ct), sigma=0.3)    for i = 1..500 cells
```

This models cell-to-cell variation in library size within a cluster. The lognormal
is parameterized so the **median** cell UMI equals the cluster's observed total UMI.
σ = 0.3 gives a coefficient of variation of ~30%, typical for 10x Chromium data.

**Note**: `log_total_umis` uses the natural log (`np.log`), matching numpy's `lognormal`
parameterization where `mean` is the mean of the underlying normal distribution.

#### Step 2: Sample gene counts per cell

For each cell `i` in cell type `ct`:

```
μ_i = true_frac × cell_total_UMI[i]
```

Then sample gene count:

```
count[i] ~ NegativeBinomial(n=θ, p=θ/(θ + μ_i))
```

This uses numpy's parameterization where:
- `n` = number of successes (= θ, the dispersion parameter)
- `p` = success probability = θ / (θ + μ)
- E[count] = n(1-p)/p = θ × μ/θ = μ
- Var[count] = n(1-p)/p² = μ + μ²/θ

If θ is None or ≥ 1000, Poisson sampling is used instead.

#### Step 3: Compute cluster-average count

```
cluster_avg_count[gene, ct] = mean(count[i] for i in 1..500)
cluster_avg_total_UMI[gene, ct] = mean(cell_total_UMI[i] for i in 1..500)
```

Both are averaged across the 500 simulated cells. This mirrors the real preprocessing
where cluster-level expression is the mean across cells in that cluster.

#### Step 4: TPM normalization

```
TPM[gene, ct] = (cluster_avg_count[gene, ct] / cluster_avg_total_UMI[gene, ct]) × 10^6
```

Values below `min_tpm` (0.1) are set to 0. This matches the real pipeline's low-expression filter.

#### Step 5: Compute specificity (fold-enrichment)

```
specificity[gene, ct] = TPM[gene, ct] / mean(TPM[gene, :])
```

where the mean is across all 461 cell types. This is the same formula used in the
real preprocessing pipeline (`preprocessing.py` / `CellType_PSY.py`).

**Expected result**: Since `true_frac` is identical for all cell types, the true
specificity is 1.0 everywhere. Any deviation from 1.0 is purely noise-driven.

### 4.3 Output

For each cell type, the function returns:
- `mean_spec`: Mean specificity across 500 simulated genes
- `max_spec`: Maximum specificity across 500 simulated genes
- `p95_spec`: 95th percentile of specificity across 500 simulated genes
- `total_umi`: The cell type's total UMI (from annotation)
- `is_neuronal`: Whether the cell type is neuronal

### 4.4 Memory management

The simulation processes cells in chunks of 500 to avoid allocating a
`(500 genes × 461 cell types × 500 cells)` array at once. Counts and UMIs are
accumulated via running sums, then divided by `n_cells` at the end.

---

## 5. Why Specificity Inflates

The inflation mechanism has two interacting components:

### 5.1 Low UMI → noisy TPM

For a cell type with total UMI = 2,000 (e.g., B-cells, cluster 0):
- A gene at P10 expression (true_frac = 7.3e-08) has expected count per cell: 2000 × 7.3e-08 = 0.000146
- Most cells will have 0 counts for this gene
- Occasionally, by chance, a cell gets 1 count → that count represents 1/2000 = 5e-04 of total UMI → TPM = 500
- The cluster average is dominated by these rare nonzero cells

For a cell type with total UMI = 50,000 (e.g., large neurons):
- Same gene: expected count per cell = 50000 × 7.3e-08 = 0.00365
- Still mostly zeros, but the 1-count events represent only 1/50000 = 2e-05 → TPM = 20
- 25× lower TPM from the same stochastic event

### 5.2 Fold-enrichment amplifies the asymmetry

Specificity = TPM(ct) / mean(TPM across all cell types).

If one low-UMI cell type happens to detect the gene while most others don't:
- Numerator: high TPM (from the 1-count event in a low-UMI background)
- Denominator: low mean TPM (most cell types have TPM = 0)
- Result: specificity >> 1, potentially 10–30×

This is **not** biological enrichment — it's a sampling artifact amplified by normalization.

### 5.3 NB overdispersion makes it worse

With θ ≈ 1.2, the variance is μ + 0.83μ², much larger than Poisson (μ).
This means the occasional nonzero counts can be larger (2, 3, ...) rather than just 0 or 1,
further inflating the TPM in low-UMI cell types.

---

## 6. Key Results

### 6.1 Expression level is the primary driver

| Expression level | Max simulated specificity (low-UMI) | Max simulated specificity (high-UMI) |
|-----------------|-------------------------------------|--------------------------------------|
| P10 | ~35× | ~2× |
| P25 | ~12× | ~1.5× |
| P50 | ~3× | ~1.2× |
| P75–P95 | ~1.5× | ~1.0× |

At P10 expression (bottom 10% of genes), low-UMI cell types show simulated specificity
up to 35×, far above the true value of 1.0 and well above the 2× cap.

### 6.2 UMI depth is the secondary driver

Within each expression level, specificity inflation is monotonically decreasing
with total UMI. The Low UMI group (<10k) always shows the most inflation.

### 6.3 Cap at 2× effectively neutralizes inflation

A cap at 2× clips the extreme tail while preserving the ~1.0 specificity for
well-sampled cell types and moderately/highly expressed genes.

---

## 7. Validation Against Empirical Data

### 7.1 Simulated vs empirical max specificity

Per-cell-type comparison (P10 expression level):
- Spearman(simulated max, empirical max) ≈ 0.7–0.8
- Top-20 overlap: ~65% of cell types with highest simulated inflation are also in the top-20 by empirical max specificity

### 7.2 Simulated noise vs empirical fraction clipped

The simulated specificity SD at P25 expression correlates with the empirical fraction
of genes exceeding the 2× cap (Spearman ρ ≈ 0.7, p < 1e-50), confirming that the
simulation captures the same cell types that are empirically problematic.

---

## 8. Parameters and Reproducibility

| Parameter | Value | Source/Rationale |
|-----------|-------|------------------|
| `n_sims` | 500 | Number of simulated genes per expression level |
| `n_cells` | 500 | Cells per cell type (fixed, not from real data) |
| `theta` | 1.186 | Median of 3 fitted θ values from real data |
| `sigma` (lognormal) | 0.3 | Cell-to-cell UMI variation (~30% CV) |
| `min_tpm` | 0.1 | TPM floor matching real pipeline |
| `seed` | 42 | Reproducibility |
| `clip_threshold` | 1.942 | 2 × mean(empirical specificity matrix) |

### 8.1 Data files

| File | Description |
|------|-------------|
| `/mnt/data0/HumanBrainCellType/cluster_GeneXCell/{id}.csv.gz` | Raw single-cell counts per cluster (for noise model fitting) |
| `dat/ExpMats/HumanCT.TPM.0.1.Filt.Spec.clip100.0.lowexp.cut1e4.csv` | Unclipped specificity matrix (461 cell types) |
| `~/Work/CellType_Psy/dat/HumanCTExpressionMats/Human.CT.Exp.Entrez.csv` | Cluster-average expression (for expression level calibration) |
| `dat/annotation.xlsx` | Cell type annotation (Total UMI, Number of cells, Supercluster) |
| `results/figures/plot_data/zinb_data.pkl` | Pre-computed simulation results (pickle) |
| `results/figures/specificity_cap/noise_model_fit.csv` | Fitted θ values per cluster |

### 8.2 Code files

| File | Role |
|------|------|
| `notebooks_rebuttal/FigS4_Specificity_Cap.py` | Consolidated figure notebook (Section 3 = simulation) |
| `notebooks_rebuttal/Specificity_ZINB_Simulation.py` | Original standalone simulation notebook |

---

## 9. Limitations and Assumptions

1. **Uniform true expression**: The simulation assumes all cell types express the gene at the same fraction of total UMI. Real genes have genuine differential expression, so the simulation isolates the noise component only.

2. **Fixed n_cells = 500**: All cell types are simulated with 500 cells, regardless of their actual cell count (which ranges from 105 to 101,039). This is conservative for small clusters (which have fewer real cells and thus even noisier estimates) and slightly generous for large clusters.

3. **Single global θ**: The simulation uses one θ for all cell types, but the fitted values range from 0.79 to 1.94. Using a single value is a simplification; per-cluster θ would likely increase inflation in the smallest clusters (which tend to have lower θ, i.e., more overdispersion).

4. **Lognormal UMI variation**: Cell-to-cell library size variation is modeled as lognormal with σ = 0.3. The actual distribution may differ, but the key finding (low-UMI cell types show inflation) is robust to this choice because it is driven by the cluster-level UMI, not the within-cluster variation.

5. **No gene-level covariates**: The simulation doesn't model gene length, GC content, or other technical covariates that might affect detection probability. These would add additional noise sources that further support the need for capping.
