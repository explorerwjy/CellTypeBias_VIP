# Specificity Cap Analysis

## Why Clipping at 2x Mean Is Necessary

Expression specificity is computed as fold-enrichment:

```
S(g, ct) = TPM(g, ct) / mean_TPM(g)
```

A value of 1.0 means the gene is expressed at the average level across cell types; values > 1 indicate enrichment. Without clipping, small cell types with low sequencing depth develop inflated specificity scores (up to ~97x) because sparse expression + TPM normalization concentrates signal in the few detected genes.

The cap at **2x mean specificity (~2.0)** is a conservative guard against these technical artifacts.

## Key Findings

### 1. Low-UMI cell types are disproportionately affected

| Metric | Neuronal (n=378) | Non-neuronal (n=83) |
|--------|-----------------|---------------------|
| Mean Total UMI | ~24,000 | ~8,000 |
| Mean fraction of genes exceeding cap | 3.86% | 13.65% |
| Median max specificity | ~30 | ~55 |
| Fraction of specificity sum from tail | 15.4% | 56.7% |

- **Spearman correlation** (Total UMI vs fraction clipped): rho = -0.497, p = 3.6e-30
- **Mann-Whitney U** (neuronal < non-neuronal fraction clipped): p = 2.7e-45
- Non-neuronal cell types have **3.5x** more genes exceeding the cap than neuronal types

### 2. Heavy tail dominates bias without cap

Although only ~5.6% of all specificity values exceed the cap, these extreme values carry disproportionate weight in the bias calculation (weighted mean specificity):

- **Non-neuronal**: 56.7% of total specificity sum comes from genes exceeding the cap
- **Neuronal**: 15.4% of total specificity sum comes from genes exceeding the cap
- Without capping, a single gene with specificity 50x would outweigh 25 genes at specificity 2x

### 3. Non-neuronal superclusters have the most extreme specificity

Superclusters ranked by median max specificity (highest first):
1. Vascular
2. Upper rhombic lip
3. Fibroblast
4. Committed oligodendrocyte precursor
5. Bergmann glia
6. Microglia

All non-neuronal. Neuronal superclusters (CGE interneuron, MGE interneuron, Deep-layer IT, etc.) have lower max specificity values, closer to or below the cap.

### 4. SCZ bias is robust across cap values

Using SCZ as an example (61 risk genes), we tested cap = 1x, 2x, 3x, and no cap:

| Cap level | Effect on results |
|-----------|-------------------|
| **Cap = 1x** | Over-clips: compresses effect sizes (~0.25 max vs 0.31 at 2x), reduces discriminability between cell types |
| **Cap = 2x** (default) | Clean biological signal: MGE interneurons and deep-layer IT neurons at top |
| **Cap = 3x** | Very similar to 2x, confirming robustness |
| **No cap** | Rankings shift: CT 314 (Miscellaneous, small neuronal type) jumps to #1 with EFFECT=0.56 (2x the next), driven by inflated specificity in a low-UMI cluster |

All top 20 cell types remain neuronal across all cap levels for SCZ, because SCZ risk genes are heavily enriched for neuronal expression. The main effect of removing the cap is that small, poorly-sampled cell types generate spurious extreme rankings — even within neuronal populations.

## Mechanism

The inflation arises from the interaction of:
1. **Low cell counts** -> fewer UMIs per cell type -> noisier TPM estimates
2. **TPM normalization** -> in low-UMI types, a few detected transcripts get very high TPM
3. **Fold-enrichment** -> TPM / mean_TPM amplifies the noise into extreme specificity values

This is not biological signal — it is a technical artifact of sparse sampling.

## Notebook

Full analysis with figures: `notebooks_rebuttal/Specificity_Cap_Analysis.py` (paired `.ipynb`)

### Figures produced

All saved to `results/figures/specificity_cap/`:

| Figure | Description |
|--------|-------------|
| `specificity_cap_analysis.{pdf,png}` | 4-panel: (A) UMI vs fraction clipped, (B) neuronal vs non-neuronal boxplot, (C) max specificity by supercluster, (D) example cell type distributions |
| `specificity_tail_contribution.{pdf,png}` | 2-panel: overall specificity distribution with tail highlighted, tail contribution to bias by UMI |
| `scz_bias_cap_comparison.{pdf,png}` | 4-panel: SCZ top-20 cell types at cap = 1x, 2x, 3x, no cap |

## Conclusion

The specificity cap at 2x mean:
- Affects only 5.6% of values overall
- Has minimal impact on neuronal cell types (3.86% of genes affected vs 13.65% for non-neuronal)
- Prevents technical artifacts from dominating bias calculations in small, low-UMI cell types
- Results are robust: cap at 3x produces nearly identical rankings to 2x
