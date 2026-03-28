# Figure Changes Plan

Detailed figure-by-figure plan for the revised manuscript, based on the rebuttal letter commitments and current figure inventory.

---

## Main Figures

After revision, Figure 1 (study schematic) moves to Graphical Abstract. All subsequent figures shift down by one.

### Graphical Abstract (was Fig 1)

**Source**: `docs/MS/Figures/Fig1 - Study schematic.png`
**Changes**: None to the figure itself. Remove from main figure numbering; use as Graphical Abstract.

---

### NEW Fig 1: Shared and Distinct Cell Type Biases (was Fig 2)

**Source**: `docs/MS/Figures/Fig2 - Share and Distinct CT.png`
**Redesign**: Two-column parallel layout contrasting cross-disorder (ASD vs SCZ) with within-disorder cognitive contrast (ASD w/o ID vs ASD with ID).

**Design rationale**: The ASD with/without ID comparison is the strongest internal control — same disorder, same study design (de novo), same cohort — the only variable is cognitive function (IQ). Placing it side-by-side with ASD vs SCZ lets the reader immediately see that CGE interneurons uniquely track the cognitive axis: MGE differentiates disorders but NOT IQ groups, while CGE differentiates both.

#### Left column: ASD vs SCZ (cross-disorder)

| Panel | Content | Notes |
|-------|---------|-------|
| **A** | ASD vs SCZ bias scatter (all cell types, R=0.68) | Soften "strong correlation" in legend |
| **C** | Supercluster bar chart (ASD blue, SCZ orange) | Add legend note: MSN divergence is biologically interesting but not primary focus |
| **E** | MGE interneuron scatter (ASD vs SCZ, ΔBias=-0.096, P=2e-5) | Shows MGE differentiates ASD from SCZ |
| **G** | CGE interneuron scatter (ASD vs SCZ, ΔBias=-0.088, P=5e-5) | Shows CGE also differentiates ASD from SCZ |

#### Right column: ASD w/o ID vs ASD with ID (cognitive contrast)

| Panel | Content | Notes |
|-------|---------|-------|
| **B** | ASD vs ASD-with-ID bias scatter (all cell types) | **NEW.** Parallel to panel A. Shows overall correlation between ASD subgroups. |
| **D** | Supercluster bar chart (ASD blue, ASD-with-ID orange) | **NEW.** Parallel to panel C. Key visual: CGE bar should be the main differentiator; most other superclusters similar. |
| **F** | MGE interneuron scatter (ASD vs ASD with ID) | **NEW.** Should show MGE does NOT strongly differentiate (small ΔBias, NS or weak P). This is the critical negative control — MGE tracks disorder identity, not IQ. |
| **H** | CGE interneuron scatter (ASD vs ASD with ID, ΔBias=-0.093, P=4e-5) | Currently panel G. Shows CGE is the key differentiator for cognitive impairment. |

#### Standalone panels (below or between columns)

| Panel | Content | Notes |
|-------|---------|-------|
| **I** | Bias correlation vs gene removal (LOEUF/expression) | Currently panel B. Keep as-is; can go below the two-column grid. |

**Visual punchline**: Left column shows both MGE and CGE differentiate ASD from SCZ. Right column shows ONLY CGE differentiates ASD w/o ID from ASD with ID — MGE does not. This directly demonstrates CGE interneurons track cognitive impairment specifically.

**Moved to FigS7**: The ASD-with-ID vs SCZ comparisons (R2.1-2: CGE comparable P=0.9, MGE less biased P=7e-3) go to supplementary rather than crowding the main figure.

---

### NEW Fig 2: CGE Interneuron Convergence (was Fig 3)

**Source**: `docs/MS/Figures/Fig3 - CGE_VIP biases.png`
**Current panels**: A-B

| Panel | Current content | Change needed |
|-------|----------------|---------------|
| A | CGE interneuron bias boxplot across cohorts (VNR+, ASD, VNR-, SCZ, ASD with ID, DD/ID) | Update labels: "ASD with ID", "ASD without ID" |
| B | Bias difference (Group1 - Group2) for CGE/MGE/Lamp5 across comparisons | **Add comparison column**: "SCZ → ASD with ID" showing CGE comparable (P=0.9) and MGE less biased (P=7e-3). This is the key R2.1-2 addition. |

**Panel changes**:
- Panel A: Relabel only. SCZ protective (OR<1) and SCZ rm NDD stay in supplement.
- Panel B: Add one new comparison column. Explicitly highlight the CGE vs MGE differential pattern across all comparisons.

---

### NEW Fig 3: 22q11.2 VIP+ Bias (was Fig 4)

**Source**: `docs/MS/Figures/Fig4 - 22q11.2 VIP+.png`
**Current panels**: A-B

| Panel | Current content | Change needed |
|-------|----------------|---------------|
| A | Interneuron type boxplot (CGE, Lamp5/LHX6, MGE) for 22q11.2 | Keep as-is |
| B | VIP+ vs VIP- within CGE for human 22q11.2 and mouse 16qA13 | Keep as-is |

**No panel changes needed.** Legend should note FDR caveat (no individual clusters pass FDR in full 22q analysis, Fig S9).

---

### Fig 4: VIP+ Calcium Dynamics (was Fig 5, experimental)

**Source**: Not in `docs/MS/Figures/` (experimental, from collaborators)
**Current panels**: A-H (imaging setup, velocity correlation, transient frequency/amplitude by subtype)

| Change | Details |
|--------|---------|
| Add "n" of mice | To legend for all subtype comparisons (R2 Minor 4) |
| Clarify narrative | Phenotypes = intrinsic excitability + encoding, not motor (R1) |
| **NEW panel?** | Consider adding ephys characterization panel if data from Erica arrives. Otherwise new supplementary figure. |

---

### Fig 5: Spatial Coding and Population Activity (was Fig 6, experimental)

**Source**: Not in `docs/MS/Figures/` (experimental, from collaborators)
**Current panels**: A-E (SVM decoding, place fields, pairwise correlations)

| Change | Details |
|--------|---------|
| Add "n" of mice | To legend (R2 Minor 4) |
| Transition text | Clarify rationale connecting activity (Fig 4) → spatial coding (Fig 5) in manuscript text (R2 Minor 5) |

---

## Supplementary Figures

### Existing supplementary figures (no renumbering needed for S1-S9)

| Figure | Content | Changes needed |
|--------|---------|---------------|
| **FigS1** | FSIQ Distribution | None |
| **FigS2** | Temporal Expression | None |
| **FigS3** | Bias Calculation Workflow | None |
| **FigS4** | Effect of Number of Genes | **Update legend/text**: Reframe as showing stability across broad range, not just N=61. Reference new gene set expansion figure (R3.2). |
| **FigS5** | Complete Biases (ASD, ASD with ID, SCZ, DD/ID) | **Add panel E**: VNR- complete bias (if not already present). Update labels to "ASD with ID". Currently has panels A-D. |
| **FigS6** | Bias Corr by Expression Level | None |
| **FigS7** | Shared and Distinct CT (detailed superclusters) | **Major update**: (1) Add ASD vs ASD-with-ID supercluster scatter panels (Hipp CA1-3, Upper-layer IT, Deep-layer IT, etc.) to mirror existing ASD vs SCZ panels. (2) Add ASD-with-ID vs SCZ comparison panels (CGE comparable P=0.9, MGE less biased P=7e-3) — moved here from main figure per R2.1-2. (3) Update all labels. |
| **FigS8** | Cross-disorder MGE/Lamp5 bias | None |
| **FigS9** | Complete Biases for 22q11.2 | **Mark FDR**: Add explicit FDR annotation showing no individual clusters pass multiple-testing correction. Add note to legend. Currently shows -log10(P) boxplot without FDR marking. |
| **FigS10** | VIP+ Cell imaging setup (experimental) | Add "n" of mice |

### New supplementary figures (S11 onwards, or interleave as needed)

Proposed numbering for new supplementary figures:

| New Fig | Content | Source file | Panels | Status |
|---------|---------|-------------|--------|--------|
| **FigS11** | QQ plots: random vs matched null + negative controls | `systematic_comparison/FigS_QQ_plot_with_negctrl.png` | 8 panels: SCZ, ASD(all), ASD(highIQ), ASD(lowIQ), DDD, 22q, HDL, Alanine. Each shows random (gray) and Best-of-N (blue) with neuron/non-neuron split. | Generated |
| **FigS12** | Significance across matching methods (heatmap) | `systematic_comparison/Fig_significance_heatmap.png` | 3 panels: SCZ, ASD_LIQ, DDD_61. Heatmap of -log10(q-value) across Random, Gene-by-gene, Rejection, Best-of-N, PropWeight, SIS. Key superclusters highlighted. | Generated |
| **FigS13** | Gene set expansion: real vs random additions | — | Need to generate. Show bias-profile correlation as gene set expands: real ranked genes vs weight-transferred random additions for each disorder. | **NOT YET GENERATED** |
| **FigS14** | Downsampling stability across disorders | `FigSX_downsampling_stability_with_overlap.pdf` | 3 panels: ASD, SCZ, DDD. Dual-axis: bias correlation (colored) + gene overlap (gray) vs sample fraction. Shows de novo (ASD/DDD) stable at 10%, SCZ unstable below 50%. | Generated |
| **FigS15** | Negative control traits bias profiles | `FigSX_negative_controls_boxplot_full.pdf` | 10 panels: ASD w/o ID, ASD with ID, HDL, IBD, SCZ, DDD, HbA1c, Parkinson's, Alanine AT, Alzheimer's. Full supercluster boxplots. Brain disorders show neuronal enrichment; non-brain show no pattern. | Generated |
| **FigS16** | SCZ protective-direction (OR<1) gene bias | `FigSX_SCZ_protective_bias.pdf` | 2 panels: (Left) mutation bias by supercluster boxplot, CGE at bottom = most negative. (Right) -log10(P) showing CGE depleted. | Generated |
| **FigS17** | SCZ with NDD genes removed | `FigSX_SCZ_rmNDD_combined.pdf` | 2 panels: (A) Scatter SCZ vs SCZ-rm-NDD285 (rho=0.923). (B) CGE boxplot: ASD w/o ID vs SCZ vs SCZ rm NDD61 vs SCZ rm NDD285. CGE signal persists (P<1e-6). | Generated |
| **FigS18** | Cross-species CCKBC mapping and 22q bias | `FigSX_CCKBC_22q_bias.pdf` | 2 panels: (Left) Three-way boxplot: VIP- CCKBC vs VIP+ CCKBC vs VIP+ ISI. (Right) CCKBC fraction vs 22q bias scatter (rho=0.05, NS). | Generated |
| **FigS19** | DRD1-5 expression across superclusters | `drd_expression/DRD_per_cluster_all5.png` | 5 sub-panels (DRD1-5) showing TPM by supercluster. CGE interneurons highlighted in red — shows DRD expression is present but not unique to VIP/CGE. | Generated |
| **FigS20a** | Specificity cap: empirical inflation by UMI depth | `specificity_cap/specificity_tail_contribution.png` | 2 panels: (A) Specificity distribution with cap=2 marked. (B) Fraction of specificity sum from genes exceeding cap vs UMI depth — non-neuronal (blue) dominates tail. | Generated |
| **FigS20b** | Specificity cap: NB simulation | `specificity_cap/zinb_simulation_main.png` | 3 panels: (A) Heatmap of median max simulated specificity by expression level x UMI. (B) Max simulated specificity by UMI tertile. (C) P95 simulated specificity vs UMI. | Generated |
| **FigS20c** | Specificity cap: sensitivity sweep | `specificity_cap/cap_sensitivity_figure.png` | 3 panels: (A) Spearman rho vs cap=2 across disorders. (B) Mean supercluster bias vs cap level. (C) Supercluster rank vs cap level. Neuronal stable 1.5-3x; non-neuronal spike at 5-10x. | Generated |
| **FigS20d** | Specificity cap: TDEP comparison + LOO stability | `specificity_cap/tdep_loo_combined.png` | 3 panels: (A) Neuronal fraction in top-N: capped vs TDEP uncapped. (B) LOO stability violins (capped vs uncapped, SCZ/ASD). (C) Rank scatter for most influential gene removal. | Generated |

### Experimental supplementary figures (waiting on data)

| New Fig | Content | Status |
|---------|---------|--------|
| **FigSXX** | VIP+ ephys characterization (Erica & Steph) | Waiting on raw data |
| **FigSXX** | CCKBC waveform analysis limitation evidence | Waiting on wetlab |

---

## Panel-Level Change Summary

### Panels to CREATE for main Fig 1 (two-column redesign):

1. **Fig 1B**: Overall scatter, ASD vs ASD-with-ID (parallel to panel A)
2. **Fig 1D**: Supercluster bar chart, ASD vs ASD-with-ID (parallel to panel C)
3. **Fig 1F**: MGE scatter, ASD vs ASD-with-ID — expected: no strong differential (critical negative control)
4. **Fig 1H**: CGE scatter, ASD vs ASD-with-ID (current panel G, re-lettered) — the key result

### Panels to ADD to existing supplementary figures:

5. **FigS7**: ASD vs ASD-with-ID supercluster scatter panels (parallel to existing ASD vs SCZ panels A-D)
6. **FigS7**: ASD-with-ID vs SCZ comparison panels (CGE P=0.9, MGE P=7e-3) — R2.1-2 request moved here
7. **Fig 2B** (was Fig 3): Add "SCZ → ASD with ID" comparison column

### Panels to MODIFY:

1. **Fig 1A**: Soften "strong correlation" in legend (R=0.68)
2. **Fig 1C**: Current panel C; add legend note about MSN
3. **Fig 2A**: Relabel to "ASD with ID", "ASD without ID"
4. **FigS5**: Add VNR- panel if missing; update labels
5. **FigS9**: Add FDR annotation (none significant after correction)

### Panels with NO changes:

- Fig 3A, B (was Fig 4)
- FigS1, S2, S3, S6, S8

---

## Fig 1 Panel Layout (Visual Reference)

```
┌─────────────────────────────────┬─────────────────────────────────┐
│   ASD vs SCZ (cross-disorder)   │  ASD vs ASD+ID (cognitive only) │
├─────────────────────────────────┼─────────────────────────────────┤
│ A. All cell types scatter       │ B. All cell types scatter       │
│    (R=0.68)                     │    (ASD vs ASD-with-ID)         │
├─────────────────────────────────┼─────────────────────────────────┤
│ C. Supercluster bar chart       │ D. Supercluster bar chart       │
│    (ASD blue, SCZ orange)       │    (ASD blue, ASD+ID orange)    │
├─────────────────────────────────┼─────────────────────────────────┤
│ E. MGE scatter                  │ F. MGE scatter                  │
│    ASD vs SCZ                   │    ASD vs ASD+ID                │
│    ΔBias=-0.096, P=2e-5         │    ΔBias=?, P=NS?  ← KEY       │
│    (MGE differentiates)         │    (MGE does NOT differentiate) │
├─────────────────────────────────┼─────────────────────────────────┤
│ G. CGE scatter                  │ H. CGE scatter                  │
│    ASD vs SCZ                   │    ASD vs ASD+ID                │
│    ΔBias=-0.088, P=5e-5         │    ΔBias=-0.093, P=4e-5         │
│    (CGE differentiates)         │    (CGE differentiates) ← KEY   │
└─────────────────────────────────┴─────────────────────────────────┘
                    ┌─────────────────────┐
                    │ I. Gene removal      │
                    │    (LOEUF/expr)      │
                    └─────────────────────┘

Punchline: Left column → both MGE & CGE differentiate ASD from SCZ
           Right column → ONLY CGE differentiates → CGE tracks cognition
```

---

## Label Consistency Pass

The following labels must be updated throughout ALL figures:
- "ASD" → "ASD without ID" or "ASD (w/o ID)" (when referring to high-IQ subset)
- "ASD with ID" or "ASD (with ID)" (when referring to low-IQ subset)
- Check axis labels, legends, and figure titles

Figures affected: Fig 1, Fig 2, FigS5, FigS7, FigS9, and all new supplementary figures.

---

## Numbering Cross-Reference

| Old number | New number | Title |
|------------|-----------|-------|
| Fig 1 | Graphical Abstract | Study schematic |
| Fig 2 | Fig 1 | Shared and Distinct Cell Type Biases (two-column redesign) |
| Fig 3 | Fig 2 | CGE Interneuron Convergence |
| Fig 4 | Fig 3 | 22q11.2 VIP+ Bias |
| Fig 5 | Fig 4 | VIP+ Calcium Dynamics |
| Fig 6 | Fig 5 | Spatial Coding / Population Activity |
| FigS1-S10 | FigS1-S10 | No change |
| — | FigS11-S20d | New supplementary figures (see table above) |
