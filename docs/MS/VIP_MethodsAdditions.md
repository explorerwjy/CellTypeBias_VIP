# Methods Additions Plan

Additions to the Methods section, organized by where they insert into the current flow.

## Current Methods Structure (paragraph numbers from docx)

1. **Rare mutation genetic data** (104-107) — gene sets, cohort sources, N=61
2. **Human whole brain single-cell data** (109-110) — Siletti atlas
3. **Gene expression specificity score** (112-123) — TPM → specificity → cap at 2 → centering
4. **Disorder mutation bias** (125-146) — weighted bias, Eqs 1-6
5. **Statistical significance** (148-149) — permutation framework (random null, 10k sims)
6. **Sensitivity to gene properties** (151-152) — LOEUF/expression removal analysis
7. **Difference in cell type mutation biases** (154-161) — differential bias, Mann-Whitney, FDR
8. Viruses → ... → Transient Analysis (experimental methods, 163-238)
9. **Statistical analyses** (241-244) — ANOVA, t-test, Mann-Whitney
10. **Code and Data availability** (246-249)

---

## Additions, in order of insertion

### A1. Expand "Rare mutation genetic data" (after para 107)

**Add**: ASD subgroup definition + expected-count normalization

- [ ] Define "ASD with ID" (IQ ≤ 70) and "ASD without ID" (IQ > 70) explicitly. State this is an analytical categorization, not inherent to source datasets. (R3.7)
- [ ] Document that expected-count normalization (already applied to SCZ per Singh et al.) is now extended to ASD and NDD gene sets to account for gene length and sequencing depth confounds. (R3.1)

---

### A2. Expand "Gene expression specificity score" — specificity cap (after para 119)

**Current text** (para 119): "To prevent extreme values from a small number of genes from dominating downstream bias computations, we capped the specificity score at 2× the mean..."

**Expand with**:
- [ ] Rationale for 2x cap: extreme specificity values (up to ~97x) arise from sparse sampling in cell types with low UMI depth, not genuine enrichment. Without capping, a single gene at 50x would outweigh 25 genes at 2x, allowing technical artifacts to dominate bias estimates for low-UMI cell types. (R3 Minor 2)
- [ ] Brief statement: robustness validated via NB simulation, cap sensitivity sweep (1x-10x), comparison with TDEP-sLDSC specificity, and leave-one-out stability analysis (see Supplementary Note 1 and Fig S20a-d). (R3 Minor 2)

*Keep it to ~3-4 sentences. All detail goes to Supp Note 1.*

---

### A3. Expand "Statistical significance" (para 148-149) — add matched null model

**Current text**: Describes only random null (size-matched).

**Expand with**:
- [ ] Describe both null models and their rationale: (R3.1)
  - **Random null** (primary): 10,000 random gene sets matched on set size. Serves as baseline.
  - **Confounder-matched null** (Best-of-N): For each gene set, generate 1,000 candidate null sets matched on gene length (CDS bases), whole-brain expression level (WB), and evolutionary conservation (mean phastCons), then select the candidate whose joint distribution best matches the real gene set (minimizing summed KS statistics). 10,000 matched null sets per trait.
- [ ] State that for SCZ, confounders (gene length, mutation rate) were already accounted for in the expected mutation counts from the burden analysis (Singh et al.). Extension to ASD/NDD described above (A1).
- [ ] Explain why LOEUF was NOT included as a matching variable: mutation-intolerant genes are genuinely enriched among disorder genes due to negative selection, and matching on LOEUF would remove the biologically meaningful signal. LOEUF is a consequence of the biology, not a technical confounder. (R3.1)
- [ ] Report that FDR values under both null models are provided for all traits in Supplementary Table X. (R3.1)

---

### A4. NEW section: "Gene set size and robustness" (after Sensitivity section, ~para 152)

- [ ] State that mutation bias results are robust to gene set size across the range N=10 to N=200 (Fig S4). Explain why: in the testable range, additional genes carry genuine (though weaker) association signal, so they refine rather than dilute cell-type rankings. (R3.2)
- [ ] Note that the default N=61 corresponds to the ASD gene set size and serves as a practical reference point; results are not specific to this cutoff. (R3.2)
- [ ] Reference new gene set expansion analysis (real ranked genes vs weight-transferred random additions; Fig S13) to confirm that added genes carry signal. (R3.2)

---

### A5. NEW section: "Impact of genetic architecture on inference stability" (after A4)

- [ ] Brief paragraph on downsampling analysis: cohorts subsampled at 10-100% in 10% increments, 100 bootstrap iterations per fraction. At each fraction, new gene weights derived from the subsampled cohort, and cell-type bias profiles recomputed. Stability measured as Spearman correlation with full-cohort bias and gene overlap with full gene set. (R3.3)
- [ ] State key finding: de novo designs (ASD, DDD) achieve r > 0.9 at 10% sampling (~4,300 probands for ASD); case-control (SCZ) requires >50% (~12,000 cases) due to weaker per-gene signals and polygenic architecture. (R3.3)
- [ ] Point to Supplementary Note 3 for detailed analysis of per-gene effect size differences and statistical power considerations. (R3.3)

---

### A6. NEW section: "Control analyses" (after A5 or after "Difference in cell type mutation biases")

- [ ] **Negative control traits**: Applied the full pipeline to non-brain traits from UK Biobank rare-variant burden analyses (HDL cholesterol, IBD, HbA1c, Alanine aminotransferase) and neurological traits (Parkinson's disease, Alzheimer's disease). Top 61 genes per trait selected by LoF burden p-value. Used as specificity controls — these should not show neuronal enrichment. (R3.3)
- [ ] **SCZ protective-direction genes**: Selected top 61 genes with OR < 1 from SCHEMA. If CGE enrichment in SCZ reflects risk biology, protective-direction genes should show the opposite pattern (CGE depletion). (R3.4)
- [ ] **SCZ with NDD genes removed**: To test whether CGE signal in SCZ is driven by overlap with NDD genes, removed genes overlapping with DDD top 61 (1 gene) and top 285 (8 genes) and recomputed bias. CGE enrichment persists (P < 1e-6 for both). (R3.7)
- [ ] **VNR(+) as internal control**: Clarify that VNR(+) gene set was included primarily as a control analysis — these genes show opposite bias direction, consistent with negative selection against cognition-reducing variants. (R3.4)

---

### A7. NEW section: "Cross-species CCKBC mapping" (after A6, brief)

- [ ] 2-3 sentences: Integrated mouse M1 patch-seq data (Allen Institute; 5,764 cells including 333 CCKBCs) with human Siletti atlas using Harmony and scVI. Five human CGE clusters received majority CCKBC assignments; these split by VIP expression into VIP-negative (low 22q bias) and VIP-positive (high 22q bias, comparable to ISI clusters). See Supplementary Note 2 for full methods and results. (R2.6, R3.6)

---

## Checklist summary

| # | Where it goes | What | Reviewer |
|---|---------------|------|----------|
| A1 | After "Rare mutation genetic data" | ASD subgroup definition + expected-count normalization | R3.7, R3.1 |
| A2 | Expand specificity cap paragraph | Cap rationale + point to Supp Note 1 | R3 Minor 2 |
| A3 | Expand "Statistical significance" | Both null models + LOEUF rationale + Supp Table ref | R3.1 |
| A4 | New section after Sensitivity | Gene set size robustness | R3.2 |
| A5 | New section after A4 | Genetic architecture / downsampling | R3.3 |
| A6 | New section after A5 | Control analyses (neg ctrl, SCZ protective, SCZ rm NDD, VNR+) | R3.3, R3.4, R3.7 |
| A7 | New section after A6 | Cross-species CCKBC mapping (brief) | R2.6, R3.6 |
