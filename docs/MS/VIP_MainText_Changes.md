# Main-Text Changes Checklist

Changes to the manuscript implied by reviewer responses. Organized by section.
Items marked **[NEW]** were identified during the 2025-03-21 cross-check of the rebuttal letter against this checklist.

---

## Introduction

- Reorganize to streamline narrative and pivot earlier toward the central finding that CGE-derived interneurons represent a shared locus of vulnerability for cognitive dysfunction (R1)
- Reduce emphasis on atlas-level mapping; foreground CGE convergence result (R1)
- Distinguish similarities and differences in neuronal impact across disorders (R2.M3)
- **[NEW]** Reposition the mouse experimental narrative: clarify that in vivo work is a circuit-level test of how CGE-lineage interneuron dysfunction translates rare mutation burden into altered cognitive representations, not a motor study (R1)

## Results

### Shared & distinct biases (current Fig 2)

- Add full ASD with ID vs. ASD comparison for both MGE and CGE interneurons (R2.1-2, Figure 2G equivalent for MGE)
- Add SCZ vs. ASD with ID comparison (R2.1-2, Figure 3)
- Explicitly highlight CGE vs. MGE differential biases when comparing ASD with ID to SCZ (R2.1-2)
- **[NEW]** Add text contextualizing MSN/excitatory divergence (Fig 2C/D) as biologically interesting but not the primary focus of the manuscript (R2.M3)

### CGE convergence across phenotypes (current Fig 3)

- Clarify that within-disorder bias and between-disorder contrast are distinct analytical frameworks (R3.5)
- Note that CGE interneurons are not significantly enriched in ASD without ID in within-disorder analysis (R3.5)
- Include VNR(+) as internal control analysis; clarify it is primarily a control (R3.4)
- Add SCZ protective-direction gene analysis (OR < 1) showing CGE depletion (R3.4)
- Add SCZ with NDD genes removed analysis confirming CGE signal persists (R3.7)
- **[NEW]** Reframe Figure S4 discussion: emphasize stability across a broad range of gene-set sizes, not specific to N=61; reference new real-vs-random expansion analysis (R3.2)

### 22q11.2 and VIP (current Fig 4)

- Present cross-species CCKBC mapping in Supplementary Information; discuss limitations of cross-species VIP subtype mapping (R2.6, R3.6)
- **[NEW]** State main-text conclusion from CCKBC analysis: VIP expression status, not CCKBC identity, is the primary predictor of 22q mutation bias within CGE interneurons (R2.6, R3.6)
- Explicitly mark FDR-adjusted significance in 22q11.2 figure (Fig. S9); note that no individual clusters pass FDR (R3.6)
- Clarify rationale for focusing on CGE/VIP rather than excitatory neurons in functional experiments: CGE/VIP show selective pattern tracking cognitive impairment, unlike broadly enriched excitatory neurons (R3.6)
- **[NEW]** Add caveat that 22q11.2 analysis is inherently underpowered for genome-wide inference; goal is hypothesis-driven cross-dataset test, not de novo discovery (R3.6)

### Mouse model results (current Figs 5-6)

- Clarify that our data demonstrate VIP interneuron dysfunction consistent with, but not directly establishing, a disinhibitory mechanism (R2.3)
- Expand mechanistic framework for reward-related VIP modulation; note dopamine receptor expression evidence but temper interpretation (R2.5)
- Clarify transition rationale between Figure 5 (activity) and Figure 6 (spatial coding) (R2 Minor 5)
- **[NEW]** Clarify that the in vivo phenotypes reflect alterations in intrinsic excitability and encoding of spatial/reward information, not motor function (R1)

### Throughout

- Consistently label subgroups as "ASD with ID" and "ASD without ID" throughout (R3.7)
- Soften "strong correlation" language (e.g., R = 0.68) and similar overstatements throughout (R2 Minor 2)

## Methods

- Describe both null models (random and confounder-matched Best-of-N) and their rationale (R3.1)
- Report FDR values under each null model for all traits in new supplementary table (R3.1)
- **[NEW]** Document extension of expected-count normalization from SCZ to ASD and NDD gene sets (R3.1)
- **[NEW]** Explain rationale for not treating LOEUF/constraint as a nuisance variable in matching — it captures biologically meaningful signal, not a technical confounder (R3.1)
- Discuss how gene set size and genetic architecture jointly influence robustness of mutation-bias estimates (R3.2)
- Explain ASD gene set as practical reference point for default gene count (R3.2)
- Add discussion of inference stability dependence on genetic architecture, sample size, and study design (R3.3)
- Document specificity cap rationale (2x), NB simulation, cap sensitivity, TDEP-sLDSC comparison, and leave-one-out stability in Methods (R3 Minor 2)
- Clarify that ASD subgroup categorization is defined by analysis, not inherent to source datasets (R3.7)

## Discussion

- Add paragraph on non-VIP effects in the 22q11 mouse: what would happen recording from SST, PV, SNCG, or LAMP5 interneurons; could VIP changes relate to altered excitatory inputs (R3 General)
- **[NEW]** Address whether VIP changes in 22q mouse are intrinsic (as predicted by computational bias) or reflect altered network inputs; discuss as limitation or cite relevant evidence (R3.8)
- Reference prior work on PV/SOM dynamics as convergent evidence for circuit interpretation (R2.3)
- Identify optogenetic rescue experiments as important future direction; discuss limitations in context of multi-gene 22q11 deletion (R2.3)
- Discuss technical limitations preventing calcium-transient analysis of other CGE subtypes (R2.1-2)

## Figures

### Main figure renumbering

- Move Figure 1 (Study schematic) to Graphical Abstract; remove from main text (R2 Minor 1)
- **[NEW]** Renumber all subsequent main figures: current Fig 2 → Fig 1, Fig 3 → Fig 2, Fig 4 → Fig 3, Fig 5 → Fig 4, Fig 6 → Fig 5. Update all in-text references accordingly
- Retain Figure 2D (now Fig 1D) in main text but sharpen surrounding narrative (R2 Minor 3)

### Figure updates

- Add "n" of mice to all figure legends (R2 Minor 4)
- Add positive control gene sets with well-established cell-type preferences (R3.1)
- Add non-brain negative control traits figure (R3.3, new Supplementary Figure)
- **[NEW]** Update FigS9 to explicitly mark FDR-adjusted significance (currently unclear if updated)
- **[NEW]** Copy generated rebuttal figures from `results/figures/` into `docs/MS/Figures/` with final numbering

### New experimental figures (wetlab)

- **[NEW]** New main/supplementary figure: VIP+ ephys characterization study (Erica & Steph; waiting on data)
- **[NEW]** New supplementary figure: evidence that CCKBC activity cannot be analyzed by waveform analysis (wetlab)

## Supplementary Information

### New tables

- New supplementary table: FDR values under random and matched null models for all traits (R3.1)
- New supplementary table: complete mutation-bias estimates for all traits in long format (trait x cell type), with gene lists and mutation counts (R3 Minor 1)

### New supplementary figures

- New supplementary figure: QQ plots under random vs. matched null (R3.1, Figure R3.1)
- New supplementary figure: gene set expansion analysis (real vs. random additions) (R3.2, Figure R3.2) — **figure not yet generated**
- New supplementary figure: downsampling stability across disorders (R3.3, Figure R3.3A)
- New supplementary figure: negative control traits bias profiles, including Parkinson's and Alzheimer's alongside HDL/IBD/ALT (R3.3, Figure R3.3B)
- New supplementary figure: SCZ protective-direction gene bias (R3.4, Figure R3.4a)
- New supplementary figure: SCZ with NDD genes removed (R3.7, Figure R3.7)
- New supplementary figure: cross-species CCKBC mapping and 22q bias (R2.6/R3.6, Figure R2.6/R3.6)
- New supplementary figure: DRD expression in superclusters (R2.5, Figure R2.5)
- New supplementary figures: specificity cap validation (MC2a-d) (R3 Minor 2)

### New supplementary notes

- Cross-species CCKBC mapping detailed methods in Supplementary Note (R2.6/R3.6)
- Specificity cap validation detailed methods in Supplementary Note (R3 Minor 2)
- Genetic architecture / downsampling detailed analysis in Supplementary Note (R3.3)

---

## Production Status

Generated figures in `results/figures/` awaiting transfer to `docs/MS/Figures/`:

- `FigSX_downsampling_stability*.pdf` → R3.3A
- `FigSX_negative_controls_*.pdf` → R3.3B
- `FigSX_SCZ_protective_bias.pdf` → R3.4
- `FigSX_SCZ_rmNDD_combined.pdf` → R3.7
- `FigSX_CCKBC_22q_bias.pdf` → R2.6/R3.6
- `drd_expression/DRD_per_cluster_all5.pdf` → R2.5
- `specificity_cap/*.pdf` → MC2a-d
- `systematic_comparison/FigS_QQ_plot_with_negctrl.pdf` → R3.1

Figures NOT yet generated:

- Gene set expansion (real vs. random additions) — R3.2
- Positive control gene sets — R3.1
- Ephys characterization (waiting on raw data from Erica)
- CCKBC waveform analysis evidence (wetlab)

Rebuttal response NOT yet written:

- R3.8 (intrinsic vs. network VIP dysfunction) — paragraphs empty in rebuttal docx

