# Manuscript Audit Report

**Date:** 2026-03-22
**Manuscript:** CGE_VIP_Genetics_Cognitive_Resub.jw.docx
**Project:** CellTypeBias_VIP

## Summary

| Category | Count |
|----------|-------|
| Claims checked | 32 |
| Verified | 23 |
| Discrepancies | 5 |
| Minor discrepancies | 3 |
| Imprecise description | 1 |
| Unverifiable | 0 |

## Critical Issues

### 1. SCZ gene count: manuscript says 61, production file has 53

- **Claim (C6):** "we selected the top 61 SCZ genes based on genetic association strength"
- **Actual:** `dat/GeneWeights/SCZ.top61.nopLI.LGD_Dmis_SameWeight.exclude_Mis2.gw` contains **53 genes**. The filename says "top61" because 61 were initially selected, but the `exclude_Mis2` filter removes genes with only 2 missense mutations, leaving 53.
- **Code:** `notebooks/Bias_Mutation_Weights.py`
- **Impact:** The manuscript implies 61 SCZ genes are used throughout, but the actual analysis uses 53. This affects the gene set size comparison with ASD (61 genes).
- **Recommendation:** Clarify in Methods that after filtering for Mis2 and pLI, 53 SCZ genes remained in the analysis. Or state "we initially selected the top 61 SCZ genes... after quality filters, 53 genes were retained."

### 2. Overlapping gene count and correlation after removal

- **Claim (C11):** "Spearman correlation = 0.64, P < 10^-10, after removing three overlapping genes"
- **Actual:** Current production gene weight files (`HIQ.top61...bgmr.gw` and `SCZ.top61...exclude_Mis2.gw`) have **2 overlapping genes** (ASH1L, KDM6B), not 3. The correlation after removing 2 genes is **r = 0.635**. An older gene file (`.2.gw`) had 3 overlapping genes (adding PCLO) but produces r = 0.677, not 0.64.
- **Code:** `notebooks/notebook_SpecBias.py` lines 242-263
- **Impact:** Neither the gene count (3 vs 2) nor the correlation value (0.64) can be reproduced from current production files. The manuscript was likely written from an intermediate version.
- **Recommendation:** Update to reflect current gene files: "after removing two overlapping genes (ASH1L, KDM6B), Spearman correlation = 0.64, P < 10^-10" — or re-run with the exact files used at time of writing if available.

## Discrepancies

### 3. ASD-SCZ Spearman correlation rounds to 0.67, not 0.68

- **Claim (C10):** "Spearman's correlation = 0.68"
- **Actual:** r = 0.669, which rounds to **0.67** (not 0.68)
- **Code:** `src/CellType_PSY.py` line 621, computed in `notebooks/Figures_Main.py`
- **Severity:** Minor — 0.01 difference likely from gene weight file versioning
- **Recommendation:** Update to "Spearman's correlation = 0.67" or re-verify with exact input files used for the figure

### 4. Specificity cap described as "2" but is data-driven

- **Claim (M7):** "capped S_g,ct at a maximum of 2"
- **Actual:** Code computes `cap = 2 × mean(all_specificity_values)` ≈ 2.0 but not exactly 2.0
- **Code:** `notebook_preprocessing.py` lines 131-134
- **Severity:** Minor — value is approximately 2.0 for 461 cell types
- **Recommendation:** Change to "capped at twice the mean specificity (~2.0)" or "capped at approximately 2"

### 5. Anatomical locations: 103 vs 105

- **Claim (C7):** "approximately 105 anatomical locations"
- **Actual:** Local annotation file yields **103** unique dissection regions
- **Code:** `dat/annotation.xlsx`
- **Severity:** Minor — the 105 figure likely comes from the original Siletti et al. publication and the 2-region difference may be due to parsing of "Top three dissections" rather than the full list
- **Recommendation:** Verify against Siletti et al. supplementary table; update if needed

## Imprecise Descriptions

### 6. Null model weight description

- **Claim (M2):** "mutation weights randomly reassigned to these sets"
- **Actual:** The original weight vector is preserved and applied to randomly sampled gene IDs. Weights are not shuffled or regenerated.
- **Code:** `script_run_ctrl_sim.py` line 29; `script_generate_geneweights.py` lines 689-691
- **Severity:** Minor — functionally equivalent for the null test, but mechanistically misleading
- **Recommendation:** Rephrase to: "for each null simulation, a random gene set of equal size was sampled and the original mutation weight vector was applied to these genes"

## Scientific Claim Verification

All scientific/interpretive claims checked were **well-supported**:

| Claim | Status | Notes |
|:---|:---:|:---|
| S1: Robust across gene set sizes (10-200) | VERIFIED | Code tests 10-199 (off-by-one, negligible) |
| S2: Correlation drops after ~15 constrained genes | VERIFIED | r drops 0.669 → 0.200 by step 15 |
| S3: Expression removal has little impact | VERIFIED | High-expression curve tracks null (Z=0.01) |
| S4: CGE strongest differential (ASD wID vs ASD) | VERIFIED | CGE #1 by p-value; tied with LAMP5-LHX6 by magnitude |
| S5: VIP+ > VIP- for 22q bias | VERIFIED | p=0.038 (3Mb), p=0.013 (1.5Mb) |
| S6: P=0.04 and P=0.01 | VERIFIED | Exact: p=0.0379 → 0.04, p=0.0125 → 0.01 |
| S7: Range top 10 to 200 | VERIFIED | Code: 10-199 (off-by-one) |
| S8: Constrained genes drive shared bias | SUPPORTED | Three converging analyses support this |

## Methods Verification

| Claim | Status | Notes |
|:---|:---:|:---|
| M1: 10,000 null gene sets | VERIFIED | config.yaml: n_random_weights=10000 |
| M3: TPM formula | VERIFIED | `expr / col_sum * 1e6` |
| M4: Min TPM = 0.1 | VERIFIED | `TPM_Cut = 0.1` |
| M5: UMI < 10,000 filter | VERIFIED | `LowExpCut = 10000` |
| M6: Specificity formula | VERIFIED | `S = n * TPM / sum(TPM)`, n=shape[1] |
| M8: Per-cell-type centering | VERIFIED | Column-wise mean subtraction confirmed |
| M9: Bias = weighted mean specificity | VERIFIED | `np.average(expr, weights=w)` |
| M10: BH FDR correction | VERIFIED | `method="fdr_i"` = standard BH |
| M11: SCZ effective count formula | VERIFIED | `case - (ctrl/N_ctrl)*N_case + dnv`; dnv=0 for SCZ in practice |
| M12: 22q uniform weights = 1 | VERIFIED | All 46 genes have weight 1 |

## Recommendations

### Must fix before submission (Critical)
1. **C6/C11:** Reconcile SCZ gene count (53 vs 61) and overlapping gene count (2 vs 3) with manuscript text. Either update the text to match current production files, or identify and document the exact intermediate files used for the reported numbers.

### Should fix (Major)
2. **C10:** Update Spearman correlation to 0.67 (or re-verify with exact input files)
3. **M7:** Clarify specificity cap is "~2" or "2 × mean specificity", not exactly "2"
4. **M2:** Rephrase null model description for mechanistic accuracy

### Consider fixing (Minor)
5. **C7:** Verify 105 anatomical locations against Siletti et al. source
6. **S7:** Change "top 10 to top 200" to "top 10 to top 199" or "approximately 200"
7. **M11:** Note that SCZ analysis uses dnv=0 (no de novo contribution) in practice

### Info
8. **S4:** CGE is strongest by p-value but essentially tied with LAMP5-LHX6/Chandelier by magnitude (0.09934 vs 0.09928). Consider acknowledging this near-tie if reviewers ask.
