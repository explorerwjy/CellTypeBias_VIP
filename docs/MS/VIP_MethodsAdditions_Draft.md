# Methods Additions — Draft Text

Each section below corresponds to an insertion point in the current Methods. Text is ready for copy-editing into the manuscript docx.

---

## A1. ASD subgroup definition and expected-count normalization

*Insert after current paragraph 107 (end of "Rare mutation genetic data" section), before "Human whole brain single-cell data".*

### ASD subgroup definition

ASD probands were stratified by cognitive function into two subgroups based on full-scale IQ (FSIQ) scores available in the SPARK cohort. "ASD with ID" was defined as probands with FSIQ ≤ 70, and "ASD without ID" as probands with FSIQ > 70. This categorization is defined by our analysis and does not correspond to inherent clinical labels in the source datasets; the underlying de novo mutation data are drawn from the same cohort and sequencing pipeline, ensuring that the only systematic difference between subgroups is cognitive ability. De novo mutation counts and gene weights were recomputed separately for each subgroup using the same statistical framework applied to the full ASD cohort.

### Expected-count normalization

For the schizophrenia gene set, the expected mutation counts used in the case-control burden analysis already account for gene-level confounders including gene length (coding sequence size) and regional mutation rate (Singh et al., 2022). To ensure comparable confounder adjustment across all gene sets, we extended this expected-count normalization to the ASD and NDD (DDD) de novo gene sets. For each gene, we computed the expected number of de novo loss-of-function and damaging missense mutations given the gene's coding sequence length and trinucleotide-context-specific mutation rates, and used the ratio of observed to expected counts as the mutation weight. This normalization ensures that differences in gene length do not systematically bias cell-type enrichment estimates.

---

## A2. Specificity cap rationale (expanded)

*Expand the existing text in paragraph 119. Replace the brief mention of capping with the following.*

To prevent extreme values from a small number of genes from dominating downstream bias computations, we capped the specificity score at 2× the mean specificity (i.e., S_cap = min(S, 2)). This cap guards against technical inflation of specificity in cell types with low sequencing depth. In our atlas of 461 cell types, uncapped specificity values can reach up to ~97×, driven primarily by cell populations with low total UMI counts (e.g., vascular, fibroblast, microglia) where the fold-enrichment calculation amplifies sampling noise. Without capping, a single gene with specificity of 50× would outweigh 25 genes at 2× in the bias calculation, allowing technical artifacts in a handful of genes to dominate the result for poorly sampled cell types. The 2× threshold was chosen to sit in the optimal range that preserves discriminability among neuronal subtypes while preventing non-neuronal inflation from distorting rankings: cell-type-level Spearman rank correlations between cap = 2× and alternative caps in the 1.5–3× range exceed 0.96 for all disorders tested. The robustness of this choice was validated through negative binomial simulation of sampling noise, a cap sensitivity sweep across seven levels (1×–10×), comparison with the TDEP-sLDSC specificity metric (Yao et al., 2025), and leave-one-out gene stability analysis (Supplementary Note 1; Supplementary Figures S20a–d).

---

## A3. Statistical significance — both null models

*Replace/expand paragraph 148-149 ("Statistical significance of the cell type mutation biases").*

### Statistical significance of the cell type mutation biases

To determine the statistical significance of the calculated cell type mutation biases, we employed a permutation-based framework comparing observed bias against null distributions derived from random gene sets.

**Random null model.** For each disorder or phenotype, we generated 10,000 random gene sets matched on set size, drawn from genes with valid expression specificity scores in the atlas. For each random gene set, we computed the cell type bias profile using the same weighted-average procedure (Equation 6). The empirical P-value for each cell type was calculated as the fraction of null gene sets producing an equal or greater bias. P-values were corrected for multiple testing across cell types using the Benjamini-Hochberg false discovery rate (FDR) procedure.

**Confounder-matched null model (Best-of-N).** To evaluate whether cell-type enrichment results could be explained by gene-level confounders, we implemented an alternative permutation framework in which null gene sets are matched to the disorder gene set on three properties: gene length (number of CDS bases), whole-brain expression level (mean expression across BrainSpan developmental stages), and evolutionary conservation (mean phastCons score). All three variables were converted to percentile ranks (0–100 scale) to ensure comparability. For each null simulation, we generated 1,000 candidate gene sets of the same size, sampled uniformly from the full gene universe (excluding input genes), and scored each candidate by the weighted mean absolute difference in percentile-space means between the candidate and target gene sets across the three matching variables, with gene length weighted 3× relative to whole-brain expression and conservation. The candidate with the smallest distance was selected. This procedure was repeated 10,000 times per trait to construct the matched null distribution. As expected, confounder matching reduced the nominal significance of several traits, because many genes with matched constraint and brain expression are themselves enriched for CNS-related functions; however, the core findings—CGE enrichment in disorders with cognitive impairment—remained significant under the matched null (Supplementary Table X; Supplementary Figure S11).

**Rationale for not matching on LOEUF.** We intentionally did not include LOEUF (loss-of-function observed/expected upper bound fraction) as a matching variable. Mutation-intolerant genes (low LOEUF) are enriched among disorder-associated genes because negative selection preferentially removes damaging variants in functionally important genes—this enrichment reflects the biology of disease risk, not a technical confounder. Matching on LOEUF would therefore remove biologically meaningful signal rather than correcting for a nuisance variable.

FDR-corrected significance values under both the random and confounder-matched null models are reported for all traits in Supplementary Table X.

---

## A4. Gene set size and robustness

*Insert as a new subsection after "Sensitivity of the correlation between ASD and SCZ cell type mutation biases to gene properties" (paragraph 152).*

### Robustness to gene set size

To assess the sensitivity of mutation-bias estimates to gene set size, we recalculated cell type bias profiles using gene sets ranging from the top 10 to the top 200 genes for each disorder, ranked by statistical significance of association. Across this range, cell-type-level Spearman correlations with the primary gene set (N = 61) remained high (r > 0.8 for SCZ, ASD, and DDD), indicating that the bias profiles are not sensitive to the specific choice of N = 61 (Figure S4). The default of 61 genes corresponds to the number of genes reaching genome-wide significance in the ASD exome sequencing study (Zhou et al., 2022) and serves as a practical reference point; results are qualitatively similar across a broad range of gene set sizes.

The stability of bias profiles as additional genes are included reflects the fact that, within the testable range, lower-ranked genes still carry genuine association signal—albeit with weaker individual effects—and therefore refine rather than dilute cell-type rankings. To directly confirm this, we compared bias-profile stability when appending real ranked disease genes (ranks 62–N) versus appending weight-transferred random genes to the fixed core set of 61 genes. Real gene additions maintained high correlation with the core-set profile, whereas random additions introduced noise and reduced correlation, demonstrating that the additional ranked genes carry cell-type-relevant signal (Supplementary Figure S13).

---

## A5. Impact of genetic architecture on inference stability

*Insert as a new subsection after A4.*

### Impact of genetic architecture on inference stability

The statistical power to detect cell-type-specific mutation bias depends on the genetic architecture of the disorder and the study design used to identify risk genes. To empirically characterize this dependence, we performed a downsampling analysis: for each disorder (ASD, SCZ, DDD), we subsampled the original cohort at fractions ranging from 10% to 100% in 10% increments, with 100 bootstrap iterations per fraction. At each fraction, de novo mutation counts (ASD, DDD) or case-control burden statistics (SCZ) were recomputed from the subsampled data, new gene weights derived, and cell-type bias profiles calculated across all 461 cell types. Stability was measured as the Spearman correlation between the subsampled and full-cohort bias profiles, along with the gene overlap between the subsampled and full gene sets.

De novo-driven disorders (ASD, DDD) showed high stability even at small sample fractions: ASD achieved r > 0.9 at 10% sampling (~4,300 probands), and DDD maintained r > 0.95 across all fractions. Schizophrenia, which relies on case-control burden testing with weaker per-gene effect sizes, required substantially larger sample fractions to achieve comparable stability (r > 0.9 only above 50%, corresponding to ~12,000 cases). These differences reflect the distinct genetic architectures: de novo studies yield high per-gene signal (relative risk ~24–88 per variant) enabling stable gene rankings even in small cohorts, whereas case-control burden studies produce modest aggregate effects (OR ≈ 1.26 for constrained gene PTVs in SCZ), resulting in noisier per-gene estimates that require larger samples to stabilize (Supplementary Note 3; Supplementary Figure S14).

---

## A6. Control analyses

*Insert as a new subsection after A5.*

### Control analyses

To further validate the specificity of CGE interneuron enrichment, we performed several control analyses targeting distinct alternative explanations.

**Non-brain and neurological negative control traits.** To confirm that neuronal cell-type enrichment is specific to psychiatric and cognitive phenotypes rather than an artifact of the analytical framework, we applied the full mutation-bias pipeline to non-brain traits derived from UK Biobank exome-wide rare-variant burden analyses: HDL cholesterol, inflammatory bowel disease (IBD), glycated hemoglobin (HbA1c), and alanine aminotransferase (ALT). For each trait, we selected the top 61 genes ranked by loss-of-function burden P-value and assigned uniform weights (weight = 1 for all genes). We additionally tested two neurological diseases — Parkinson's disease and Alzheimer's disease — using the same gene selection procedure (top 61 and top 60 genes, respectively, by LoF burden P-value from UK Biobank ICD-based phenotypes). Non-brain traits showed no systematic enrichment in neuronal superclusters, with bias profiles dominated by non-neuronal cell types such as glia, vascular, and ependymal populations. The neurological diseases showed modest and diffuse neuronal enrichment that did not recapitulate the selective CGE interneuron pattern observed for psychiatric disorders with cognitive impairment (Supplementary Figure S15).

**SCZ protective-direction genes.** If CGE interneuron enrichment in schizophrenia reflects genuine risk biology, genes with protective direction of effect should show the opposite pattern. We selected the top 61 genes from the SCHEMA exome sequencing study for which rare damaging variants were associated with reduced schizophrenia risk (odds ratio < 1). Consistent with prediction, these protective-direction genes showed CGE interneuron depletion—the lowest bias among all neuronal superclusters—providing a directional control for the risk-gene finding (Supplementary Figure S16).

**SCZ with NDD genes removed.** To test whether CGE interneuron enrichment in schizophrenia is driven by overlap with neurodevelopmental disorder (NDD) genes, we quantified the gene overlap between the SCZ gene set and two NDD gene sets (DDD top 61: 1 overlapping gene; DDD top 285: 8 overlapping genes) and recomputed SCZ mutation bias after removing overlapping genes. CGE interneuron enrichment persisted in both analyses (P < 1 × 10⁻⁶), with cell-type-level Spearman correlation between the full and reduced SCZ profiles exceeding 0.92, confirming that the SCZ CGE signal is not driven by NDD gene overlap (Supplementary Figure S17).

**VNR(+) as internal control.** The verbal-numerical reasoning positive-direction gene set (VNR+), comprising genes for which rare coding mutations are associated with higher cognitive scores, was included primarily as an internal control analysis. Under strong negative selection on cognition, these genes are expected to show low or opposite mutation bias compared to risk gene sets. Consistent with this, VNR(+) genes showed the lowest CGE bias among all traits tested (Figure 2A), supporting the interpretation that CGE enrichment tracks the direction of cognitive impairment.

---

## A7. Cross-species CCKBC mapping (brief)

*Insert after "Assignment of subtype identity based on immunostaining" (para 208, after the VIPo definition), before "Pre-processing of Ca2+ imaging data" (para 210). This placement follows the subtype definitions and naturally bridges to: "we also asked whether these mouse-defined subtypes map onto human transcriptomic clusters."*

### Cross-species interneuron subtype mapping

To investigate whether the 22q11.2 mutation bias within CGE interneurons differs between cholecystokinin-expressing basket cell (CCKBC) and interneuron-selective (ISI) subtypes, we performed cross-species transcriptomic integration. Mouse patch-seq transcriptomic data from the Allen Institute (5,764 cells from visual and motor cortex, including 333 transcriptomically defined CCKBCs; Gouwens et al., 2020; Scala et al., 2021) were projected into the human Siletti et al. atlas space using two independent methods (Harmony batch correction and scVI latent-space integration via scArches). Human CGE clusters were classified as CCKBC-enriched based on the fraction of mapped mouse CCKBC cells received via k-nearest-neighbor voting (k = 30). The five imputed CCKBC clusters segregated by VIP expression into VIP-negative (low 22q bias) and VIP-positive (high 22q bias comparable to ISI clusters), indicating that VIP expression status, rather than CCKBC identity, is the primary predictor of 22q11.2 mutation bias within CGE interneurons. Full details of the integration procedure, electrophysiology signature transfer analysis, and VIP+/VIP− bias comparison are provided in Supplementary Note 2 (Supplementary Figure S18).

---

## Results text additions (brief mentions for main text)

*These are short sentences to insert into the Results narrative, not the Methods section.*

### Negative control traits (insert in Results, CGE convergence section, after presenting psychiatric results)

> To confirm that the neuronal cell-type enrichment observed for psychiatric disorders is not an artifact of the analytical framework, we applied the same pipeline to non-brain traits from UK Biobank rare-variant burden analyses (HDL cholesterol, IBD, HbA1c, alanine aminotransferase). None of these traits showed systematic enrichment in neuronal cell types, and their bias profiles were dominated by non-neuronal populations (Supplementary Figure S15). We additionally tested two neurological diseases (Parkinson's and Alzheimer's); while these showed modest and diffuse neuronal enrichment, they did not recapitulate the selective CGE interneuron pattern observed for psychiatric phenotypes with cognitive impairment (Supplementary Figure S15).

### SCZ protective genes (insert in Results, after SCZ CGE enrichment finding)

> As a directional control, we analyzed schizophrenia genes with protective direction of effect (OR < 1, i.e., genes for which rare damaging variants are associated with reduced risk). Protective-direction genes showed CGE interneuron depletion — the lowest bias among all neuronal superclusters — consistent with the prediction that risk and protective genes should show opposite cell-type patterns (Supplementary Figure S16).

### SCZ rm NDD (insert in Results, after CGE convergence across disorders)

> To rule out the possibility that CGE enrichment in schizophrenia is driven by shared neurodevelopmental disorder (NDD) genes, we removed overlapping genes (1 gene from DDD top 61; 8 genes from DDD top 285) and recomputed SCZ mutation bias. CGE interneuron enrichment persisted in both analyses (Mann-Whitney P < 1 × 10⁻⁶), confirming that the SCZ CGE signal reflects disorder-specific biology rather than NDD gene overlap (Supplementary Figure S17).

### ★ Cross-species CCKBC/ISI subtype mapping (insert in Results, mouse section, after para 84 — the subtype correlation analysis ending with "...perturbations in specific subtypes of VIP+ interneurons can disrupt network dynamics and information processing important for cognition.")

> Given the subtype-specific differences observed between CCKBC and ISI interneurons in the Df(16)A+/- mouse, we asked whether this functional distinction is reflected in differential 22q11.2 mutation bias in the human brain atlas. Using cross-species transcriptomic integration, we mapped mouse CCKBCs onto human CGE clusters and found that the five imputed CCKBC-enriched clusters segregated by VIP expression: VIP-negative CCKBCs showed low 22q bias, while VIP-positive CCKBCs showed bias comparable to VIP+ ISI clusters (Supplementary Figure S18). CCKBC identity itself did not predict differential 22q bias (Spearman ρ = 0.05, P = 0.86); rather, VIP expression status was the primary determinant of mutation vulnerability within CGE interneurons. These results are consistent with the known difficulties of cross-species Sncg subtype classification (AUROC = 0.50; Bakken et al., 2021) and suggest that the CGE/VIP supercluster — rather than individual subtypes — is the most robust unit for cross-species integration of mutation-bias signals (Supplementary Note 2).

### ★ 22q11.2 FDR caveat (insert in Results, 22q section, after presenting CGE/VIP supercluster enrichment)

> We note that while CGE and VIP+ interneurons showed the strongest mutation bias among interneuron superclusters for the 22q11.2 deletion interval (Figure 3A), no individual cell type cluster reached significance after FDR correction (Figure S9). This is expected given the composition of the 22q11.2 gene set: unlike ASD, SCZ, or DDD gene sets in which every gene carries strong individual statistical evidence of disease association, the 22q11.2 deletion interval contains 46 genes with uniform weights, only a subset of which are likely to contribute to disease risk. The inclusion of functionally neutral "passenger" genes dilutes the cell-type-specific signal, reducing power at the individual-cluster level. Accordingly, the 22q11.2 analysis is not designed for de novo cell-type discovery but rather to test the hypothesis — motivated by convergent evidence from ASD, SCZ, and NDD analyses — that CGE/VIP interneurons are preferentially affected within this locus.

### ★ FigS9 legend addition (append to existing Figure S9 legend)

> Note: no individual cell type cluster reached FDR-corrected significance (q < 0.1) for the 22q11.2 gene set. This reflects the heterogeneous composition of the deletion interval, in which all 46 genes are weighted equally regardless of their contribution to disease risk, unlike the statistically prioritized gene sets used for ASD, SCZ, and DDD. The 22q11.2 analysis is therefore interpreted at the supercluster level as a hypothesis-driven test rather than a genome-wide screen.

---

## Completeness Check


| Reviewer point | Topic                                  | Covered by | Status                                                               |
| -------------- | -------------------------------------- | ---------- | -------------------------------------------------------------------- |
| R3.1           | Both null models described             | A3         | ✓                                                                    |
| R3.1           | FDR under each model in Supp Table     | A3         | ✓                                                                    |
| R3.1           | LOEUF rationale                        | A3         | ✓                                                                    |
| R3.1           | Expected-count normalization extended  | A1         | ✓                                                                    |
| R3.2           | Gene set size robustness               | A4         | ✓                                                                    |
| R3.2           | ASD as practical reference for N=61    | A4         | ✓                                                                    |
| R3.3           | Genetic architecture / downsampling    | A5         | ✓                                                                    |
| R3.3           | Negative control traits                | A6         | ✓ (moved from A5 to A6)                                              |
| R3.4           | SCZ protective (OR<1)                  | A6         | ✓                                                                    |
| R3.4           | VNR(+) as control                      | A6         | ✓                                                                    |
| R3.7           | ASD subgroup definition                | A1         | ✓                                                                    |
| R3.7           | SCZ rm NDD                             | A6         | ✓                                                                    |
| R3 Minor 2     | Specificity cap rationale + validation | A2         | ✓                                                                    |
| R2.6/R3.6      | Cross-species CCKBC (Methods)          | A7         | ✓ (placed after immunostaining subtype definitions)                  |
| R2.6/R3.6      | Cross-species CCKBC (Results)          | Results ★  | ✓ (after mouse subtype correlation para 84)                         |
| R3.6            | 22q FDR caveat (Results)               | Results ★  | ✓                                                                    |
| R3.6            | FigS9 legend FDR note                  | Legend ★   | ✓                                                                    |
| R3.1           | Positive control gene sets             | —          | **NOT YET** — need to add if we include positive controls in Methods |


