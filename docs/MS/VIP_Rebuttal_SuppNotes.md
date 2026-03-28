# Supplementary Notes

These supplementary notes accompany the reviewer response and contain detailed methods, simulations, and analyses that support the main rebuttal but are too lengthy for inclusion in the point-by-point response.

---

## Supplementary Note 1: Specificity Cap Validation

*Supporting response to Reviewer 3, Minor Point 2*

### 1.1 Rationale for capping

Expression specificity is computed as fold-enrichment: $S(g, \text{ct}) = \text{TPM}(g, \text{ct}) / \overline{\text{TPM}(g)}$, where values greater than 1 indicate above-average expression in a given cell type. We cap this score at 2x the global mean specificity as a conservative guard against technical inflation in cell types with low sequencing depth.

In our atlas of 461 cell types, uncapped specificity values can reach up to ~97x. These extreme scores are artifacts of sparse sampling, not genuine expression enrichment, and are driven primarily by small cell populations (e.g., vascular, fibroblast, microglia, and cerebellum neurons) with low total UMI counts. The inflation arises because cell types with low total UMI per cell produce noisier TPM estimates, and the fold-enrichment calculation amplifies this noise into extreme specificity values. Quantitatively, without the cap, 56.7% of the total specificity sum in non-neuronal cell types comes from genes exceeding the threshold, compared to only 15.4% for neuronal types (Spearman R between UMI depth and fraction of genes exceeding cap = -0.497, p = 3.6 x 10^-30). This means that in the bias calculation---a weighted mean of specificities---a single uncapped gene at 50x would outweigh 25 genes at 2x, allowing technical artifacts in a handful of genes to dominate the result for poorly-sampled cell types (Supplementary Figure MC2a).

### 1.2 Negative binomial simulation

To demonstrate the inflation mechanism, we performed a simulation using an empirically calibrated negative binomial (NB) noise model. We first fit the mean-variance relationship ($\text{Var} = \mu + \mu^2 / \theta$) from raw single-cell counts across representative clusters spanning the UMI range and confirmed overdispersion (Var/Mean > 1) in all clusters. We then simulated genes with "uniform" true expression across all 461 cell types---meaning the true specificity is exactly 1.0 everywhere---using each cluster's actual total UMI depth and the fitted NB noise parameters. Any deviation from 1.0 in the simulated specificity is therefore purely due to sampling noise.

The results (Supplementary Figure MC2b) reveal that specificity inflation is jointly driven by two factors: gene expression level and cell-type UMI depth. For lowly expressed genes (bottom 5th-10th percentile of the real expression distribution), simulated specificity in low-UMI cell types can reach 20x, orders of magnitude above the true value of 1.0. In contrast, moderately or highly expressed genes remain near 1.0. Even at the 95th percentile of simulated specificity, low-expression genes in low-UMI cell types substantially exceed the 2x cap, while high-expression genes stay flat near 1.0. These results demonstrate that the extreme specificity values we observe are an expected and unavoidable consequence of NB sampling noise in scRNA-seq data, particularly for lowly expressed genes in cell types with low total UMI.

### 1.3 Cap sensitivity sweep

To confirm the robustness of our choice, we recomputed mutation bias for SCZ, ASD (without ID), and DDD across seven cap levels (1x, 1.5x, 2x, 2.5x, 3x, 5x, 10x). Cell-type-level Spearman rank correlations between cap = 2x and caps in the 1.5-3x range exceed 0.96 for all three disorders (Supplementary Figure MC2c, Panel A), confirming high stability. At cap = 1x, the correlation drops sharply (ASD: 0.41) because over-clipping compresses specificity variation and makes neuronal superclusters indistinguishable. At higher caps (5-10x), ASD rankings diverge substantially (rho < 0.90) as non-neuronal artifacts begin to dominate. The mechanism is visible in Panel B-C: non-neuronal cell types (Vascular, Microglia, Astrocyte) with inflated specificity due to low UMI depth gain bias steeply as the cap is raised, while neuronal types (CGE, MGE, MSN) plateau by 2x. In terms of supercluster ranking, at cap = 2x non-neuronal types rank 20th or lower, safely below all key neuronal populations. At cap = 10x, Microglia climbs from rank 26 to 3 and Vascular from 20 to 5, overtaking MSN (which drops from rank 10 to 25). The 2x cap thus sits in the optimal range: preserving discriminability among neuronal subtypes while preventing technical inflation in low-UMI populations from distorting rankings.

### 1.4 Comparison with TDEP-sLDSC specificity and leave-one-out stability

The most directly comparable approach is the cell-type specificity used in TDEP-sLDSC (Naqvi et al., *Nature Communications* 2024), which computes specificity as the proportion of a gene's total expression attributed to each cell type, applied to the same Siletti et al. brain atlas used in our study. Critically, TDEP-sLDSC uses this specificity for rank-based thresholding---selecting the top 10% most specific genes per cell type for LDSC analysis, rather than as a continuous multiplier. Because each cell type's gene set is defined independently, extreme specificity values do not bias comparisons across cell types. Our mutation-bias framework, by contrast, uses specificity as a continuous weight in a gene-level weighted sum, where the absolute magnitude of the specificity score directly determines the bias estimate.

To test whether the TDEP specificity metric is compatible with our framework, we computed mutation bias using their specificity matrix (scaled to fold-enrichment and mean-centered per gene, matching our preprocessing) for all three disorders. The results (Supplementary Figure MC2d, Panel A) show that without capping, non-neuronal cell types dominated the top-ranked positions, particularly for ASD without ID, where the neuronal fraction in the top-5 dropped to 0% under TDEP, compared to 40% with our capped metric. SCZ and DDD maintained 100% neuronal enrichment with capping but showed reduced enrichment under TDEP. Across all three disorders, our capped metric maintained consistently higher neuronal enrichment in top-ranked cell types, while the uncapped metric showed substantially reduced enrichment.

Beyond overall cap-level sensitivity, capping also prevents single-gene dominance in the bias calculation. Without capping, a single gene with extreme specificity (e.g., 21x or higher) can single-handedly drive a cell type to the top of the rankings, making results unstable. To quantify this, we performed a leave-one-out (LOO) stability analysis: for each gene in the SCZ and ASD gene sets, we removed it and recomputed mutation bias across all 461 cell types, then measured the Spearman rank correlation with the full-set result (Supplementary Figure MC2d, Panels B-C). With capping at 2x, removing any single gene produces minimal change: the minimum Spearman rho across all genes is 0.955 for SCZ and 0.978 for ASD. Without capping, however, LOO stability drops substantially---minimum rho = 0.909 (SCZ) and 0.911 (ASD)---and the distributions shift markedly lower. The most influential genes are those with extreme uncapped specificity in specific cell types: *CACNA1G* (uncapped specificity of 21x in Purkinje neurons) and *TBR1* (11.7x in layer 6 IT neurons). Removing these genes under the uncapped metric causes large rank rearrangements, whereas the same removal under the capped metric has negligible effect.

### Figures

- **Supplementary Figure MC2a**: Empirical specificity inflation is driven by low UMI depth. (A) Fraction of genes exceeding the specificity cap (2x) as a function of total UMI per cell type. (B) Specificity distribution for representative cell types spanning the UMI range.

- **Supplementary Figure MC2b**: Negative binomial simulation demonstrates that specificity inflation arises from sampling noise in lowly expressed genes. (A) Heatmap of median maximum simulated specificity. (B) Maximum simulated specificity across UMI tertiles. (C) 95th-percentile simulated specificity vs. total UMI, stratified by expression level.

- **Supplementary Figure MC2c**: Mutation bias rankings are robust to specificity cap level in the 1.5-3x range. (A) Spearman correlation between cap = 2x and alternative cap levels. (B) Mean supercluster bias across cap levels for SCZ. (C) Supercluster ranking across cap levels.

- **Supplementary Figure MC2d**: TDEP-sLDSC uncapped specificity produces non-neuronal dominance, and capping ensures leave-one-out stability. (A) Neuronal fraction in top-N ranked cell types. (B) LOO stability distributions for SCZ and ASD. (C) Rank scatter for most influential gene removal.

---

## Supplementary Note 2: Cross-Species CCKBC Mapping

*Supporting responses to Reviewer 2, Point 6 and Reviewer 3, Point 6*

### 2.1 Cross-species transcriptomic integration

We integrated mouse patch-seq data from the Allen Institute (5,764 cells from visual and motor cortex (Gouwens et al., 2020; Scala et al., 2021), including 333 transcriptomically defined CCKBCs from the Vip and Sncg subclass) with the human Siletti et al. atlas (175,000 CGE interneurons, 21 clusters) using two independent integration methods: Harmony batch correction and scVI latent-space integration (scArches). In both methods, we projected mouse cells into the human reference space and assigned each to its nearest human cluster via k-nearest-neighbor voting (k = 30). Five human CGE clusters (277, 278, 279, 280, 281) received a majority of CCKBC assignments (Harmony CCKBC fraction > 50%), with moderate concordance between methods (Spearman rho = 0.46, P = 0.03). Among these, clusters 279, 280, and 281 showed the highest CCKBC fractions under both Harmony (67-96%) and scVI (86-90%) and are annotated as VIP-expressing (INT-VIP) in the Siletti atlas, while clusters 277 and 278 are VIP-negative.

### 2.2 Electrophysiology signature transfer

As an independent check, we tested whether the mouse CCKBC electrophysiological signature transfers to the imputed human clusters. We identified 5 shared electrophysiological feature pairs between mouse M1 patch-seq data and human patch-seq data (Lee et al., 2023): AP width, AP threshold, ISI CV, input resistance, and membrane time constant. We defined a mouse CCKBC direction vector from these features and scored human CGE cells for CCKBC-likeness in this shared feature space. k-nearest-neighbor classification (k = 5) showed that the fraction of mouse CCKBC neighbors did not differ between imputed CCKBC and VIP+ ISI human cells (Mann-Whitney P = 0.79; Kruskal-Wallis three-group P = 0.96), and logistic regression trained on mouse cells yielded low CCKBC probability across all human groups. This result indicates that the mouse CCKBC electrophysiological signature does not distinguish human CGE subtypes, consistent with the known evolutionary divergence in Sncg interneuron identity (Bakken et al., 2021), and further supports the conclusion that VIP lineage — rather than subtype-level identity — is the most informative predictor of mutation vulnerability.

### 2.3 VIP+/VIP- CCKBC split and 22q mutation bias

The five imputed CCKBC clusters split by VIP expression into two groups with divergent properties:

- **VIP-negative CCKBCs** (clusters 277, 278): Express RELN, SP8, and NDNF; consistently exhibit low 22q11.2 mutation bias.

- **VIP-positive CCKBCs** (clusters 279, 280, 281): Express VIP alongside CCK and CNR1; exhibit high 22q mutation bias comparable to VIP+ ISI clusters.

We tested whether imputed CCKBCs (all five clusters combined) show differential 22q bias compared to VIP+ ISI clusters. Four 22q-related gene sets were tested: 22q deletion interval, 22q mouse model deletion, highly expressed 22q genes, and 22q DEGs (day 75). None reached significance (one-tailed Mann-Whitney P = 0.37, 0.63, 0.50, and 0.11, respectively). The three-way comparison revealed that VIP expression status, rather than CCKBC identity, is the primary predictor of 22q mutation bias.

### 2.4 Limitations of cross-species subtype mapping

Cross-species classification of Sncg subtypes (Scala et al., 2021) yields near-chance accuracy (AUROC = 0.50; Bakken et al., 2021), the worst of any interneuron subclass. CCK mRNA is expressed ubiquitously across human neurons (Darmanis et al., 2015), which is also what we observed in Siletti et al. (2023), precluding its use as a discriminative marker. These results demonstrate that:

1. Cross-species subtype mapping at the CCKBC/ISI level is limited by evolutionary divergence in Sncg interneuron identity, precluding direct one-to-one correspondence.
2. VIP-negative CCKBC clusters show low 22q bias, while VIP-positive CCKBC clusters show bias comparable to VIP+ ISI cells, indicating that VIP expression status rather than CCKBC identity drives 22q vulnerability.
3. The broader CGE/VIP supercluster---rather than individual subtypes---remains the most robust unit for cross-species integration of mutation-bias signals.

### Figure

- **Supplementary Figure R2.6/R3.6**: Imputed human CCKBC identity does not predict differential 22q11.2 mutation bias within CGE interneurons. (Left) Three-way comparison of 22q deletion mutation bias across VIP-negative CCKBCs, VIP-positive CCKBCs, and VIP+ ISI clusters. (Right) CCKBC fraction (Harmony mapping) vs 22q deletion bias across all 16 VIP/CCKBC clusters (Spearman rho = 0.05, P = 0.86).

---

## Supplementary Note 3: Genetic Architecture and Downsampling Stability

*Supporting response to Reviewer 3, Point 3*

### 3.1 De novo vs. case-control power differences

The reduced stability of the SCZ cell-type bias estimates at low sample fractions (compared to ASD and DDD) reflects two distinct but interrelated factors:

**Study design differences**: De novo studies use Poisson rate tests that directly compare observed mutation counts against a gene-specific background mutation rate, yielding relatively precise per-gene statistics even at moderate sample sizes. Case-control burden tests, by contrast, compare mutation frequencies between two groups using Fisher's exact test, which requires substantially larger sample sizes to detect the same effect.

**Per-gene effect size differences**: Cross-disorder comparisons using the same statistical framework (Nguyen et al., 2017) estimated that the mean relative risk per de novo loss-of-function variant is ~88 for developmental disorders, and ~24 for ASD, but only ~12 for schizophrenia, a 2-8-fold difference. This pattern is consistent with the observation that ASD and NDD de novo genes show high individual penetrance: top ASD genes such as CHD8 and SCN2A have no loss-of-function carriers among >213,000 unaffected individuals (Wigdor et al., 2023), and the aggregate de novo protein-truncating variant (PTV) enrichment in constrained genes is 3.5-fold for ASD (Satterstrom et al., 2020) and ~10-30-fold for the top DDD genes (Kaplanis et al., 2020). These large effects arise because trio-based de novo designs are naturally suited to early-onset developmental disorders in which highly penetrant mutations arise spontaneously and are subject to strong negative selection.

By contrast, schizophrenia is a later-onset disorder for which large-scale trio sequencing is logistically more difficult to conduct; consequently, gene discovery has relied primarily on case-control designs, and the aggregate PTV burden across constrained genes is only OR = 1.26 (Singh et al., 2022), reflecting a more distributed, polygenic rare-variant architecture in which each individual gene contributes a modest effect.

### 3.2 Implications for mutation-bias stability

Because the mutation-bias metric weights genes by their mutation burden, the weaker and noisier per-gene signals in SCZ propagate more uncertainty into the cell-type bias estimates, and the very low gene overlap at small sample fractions (only 12% at 10% sampling) indicates that the SCZ gene rankings themselves are not yet stable at current sample sizes. These results suggest that the current SCZ cohort size (~24,000 cases) is near the minimum required for stable cell-type inference, and that future larger case-control sequencing studies will be needed to refine SCZ cell-type conclusions.

### Figure

- **Supplementary Figure R3.3A**: Cell-type mutation bias estimates as a function of cohort sample size. Downsampling analysis showing stability of cell-type bias across ASD (de novo), SCZ (case-control), and DDD (de novo). Lines and ribbons indicate mean +/- s.d. across 100 iterations.

---

## References

Bakken, T. E., et al. (2021). Comparative cellular analysis of motor cortex in human, marmoset and mouse. *Nature*, 598, 111-119.

Darmanis, S., et al. (2015). A survey of human brain transcriptome diversity at the single cell level. *PNAS*, 112, 7285-7290.

Gouwens, N. W., et al. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. *Cell*, 183, 935-953.

Kaplanis, J., et al. (2020). Evidence for 28 genetic disorders discovered by combining healthcare and research data. *Nature*, 586, 757-762.

Naqvi, S., et al. (2024). TDEP-sLDSC for cell-type-specific enrichment. *Nature Communications*.

Nguyen, H. T., et al. (2017). Integrated Bayesian analysis of rare exonic variants. *Genome Medicine*, 9, 114.

Satterstrom, F. K., et al. (2020). Large-scale exome sequencing study implicates both developmental and functional changes in the neurobiology of autism. *Cell*, 180, 568-584.

Scala, F., et al. (2021). Phenotypic variation of transcriptomic cell types in mouse motor cortex. *Nature*, 598, 144-150.

Siletti, K., et al. (2023). Transcriptomic diversity of cell types across the adult human brain. *Science*, 382, eadd7046.

Singh, T., et al. (2022). Rare coding variants in ten genes confer substantial risk for schizophrenia. *Nature*, 604, 509-516.

Wigdor, E. M., et al. (2023). The female protective effect against autism spectrum disorder. *Cell Genomics*, 2, 100134.
