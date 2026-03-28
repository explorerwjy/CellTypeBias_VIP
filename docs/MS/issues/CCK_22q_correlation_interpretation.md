# Issue: CCK TPM Correlation with 22q Bias and Consistency with Mouse CCKBC Findings

**Date**: 2026-03-23
**Status**: Open — needs discussion with co-authors
**Related analyses**: `notebooks_rebuttal/VIP_Expression_22q_Correlation.py`

---

## Summary

When correlating CGE marker gene expression with 22q11.2 mutation bias across 21 CGE interneuron clusters, the metric used (TPM vs specificity) changes the results for CCK substantially:

| Gene | TPM rho | Specificity rho | N clusters clipped at 2x |
|------|---------|-----------------|--------------------------|
| VIP  | 0.462 (P=0.035)* | 0.469 (P=0.032)* | 14/21 |
| CCK  | 0.412 (P=0.064) | 0.170 (P=0.460) | 17/21 |
| CALB2 | 0.242 (P=0.291) | 0.400 (P=0.073) | 12/21 |

VIP is the top positive correlate under both metrics. The discrepancy for CCK arises because CCK is broadly expressed (433/461 cell types, mean TPM 161 globally) and 17/21 CGE clusters exceed the 2x specificity cap, collapsing their rank variation. The TPM correlation (0.41) is real in a statistical sense but may be misleading biologically.

## The concern

The mouse experimental data (Bhatt et al., hippocampal CA1) shows:
- **ISI cells (ISI2, ISI3)**: reduced Ca2+ amplitudes, disrupted locomotion correlation, altered reward responses in Df(16)A+/- mice
- **CCKBCs**: activity preserved in these specific behavioral paradigms

If we report that CCK TPM positively correlates with 22q bias in human data, this **appears to contradict** the mouse finding that CCKBCs are unaffected. A reviewer could flag this inconsistency.

## Why specificity is the appropriate metric

1. **Matches the bias pipeline**: Our mutation bias calculation uses specificity (fold-enrichment, capped at 2x), not TPM. Reporting correlations on the same scale the pipeline uses is internally consistent.

2. **Specificity answers the right question**: "Which gene's *enrichment pattern* across CGE subtypes predicts 22q bias?" For a ubiquitously expressed gene like CCK, high TPM everywhere means it doesn't discriminate subtypes — the specificity metric correctly captures this.

3. **Clipping is not an artifact for CCK**: CCK hits the ceiling because it genuinely is not specific to any CGE subtype. The clipping accurately reflects that CCK expression doesn't distinguish high-bias from low-bias CGE clusters.

4. **VIP survives clipping**: Despite 14/21 clusters being clipped, VIP's correlation is essentially unchanged (0.462 -> 0.469) because the 7 unclipped clusters (VIP-low types: LAMP5, Pax6, Sncg) carry meaningful rank variation.

## Key counterexample: Cluster 278

Cluster 278 (Pax6-labeled, VIP-negative) has CCK TPM = 3391 (4th highest of all CGE clusters) but the **lowest** 22q bias of all 21 clusters (0.014). If CCK expression drove 22q vulnerability, this cluster should rank near the top. It doesn't.

## Partial correlation analysis

Controlling for VIP, CCK retains partial rho = 0.34 (P = 0.127, n.s.). Controlling for CCK, VIP retains partial rho = 0.41 (P = 0.067). Neither is significant with N=21, but VIP retains a larger fraction of its signal.

## Mouse data limitations to consider

1. **Region**: Mouse data is hippocampal CA1 only; human data is whole-brain atlas. CCK basket cells may be differentially vulnerable across brain regions.

2. **Task specificity**: "CCKBC activity preserved" refers to specific behavioral correlates (locomotion, reward anticipation). Molecular or synaptic deficits in CCKBCs may exist but not manifest in these particular Ca2+ imaging paradigms.

3. **Developmental timepoint**: The mouse model captures adult functional phenotypes. The human genetic data reflects developmental vulnerability patterns that may differ from adult circuit function.

## Literature: CCKBCs ARE implicated in psychiatric disorders

Contrary to the "CCKBCs are fine" interpretation from the mouse data, substantial literature shows CCKBC involvement in SCZ and cognitive deficits:

### Schizophrenia (strongest evidence)
- **Curley & Lewis (2012)** *J Physiol* — Reduced CCK mRNA and CB1R mRNA/protein in DLPFC of SCZ subjects. Shifts perisomatic inhibition balance between CCK and PV basket cells, disrupting gamma oscillations.
- **Eggan et al. (2008)** *Arch Gen Psychiatry* — CB1R mRNA reduced ~15%, protein reduced ~12-14% in DLPFC area 9. Interpreted as compensatory downregulation.
- **Dienel & Lewis (2018)** *Neurobiol Dis* — CCK mRNA downregulated in layers 2-superficial 3 of PFC in SCZ.
- **Chou et al. (2023)** *Neurobiol Dis* — Cell-type-specific CB1R alterations: decreased in inhibitory (CCKBC) boutons, increased in excitatory boutons.

### Intellectual disability
- **Jahncke et al. (2023)** *eLife* — Dystroglycanopathy (causes ID + seizures): dystroglycan required for CCK basket cell axon targeting; loss causes synaptic defects.
- **del Pino et al. (2017)** *Nat Neurosci* — ErbB4 knockout in CCKBCs disrupts spatial coding, impairs spatial memory and novelty recognition.

### Cognitive circuits
- **Nguyen et al. (2020)** *J Neurosci* — CCK interneurons required for working memory retrieval in mPFC (not PV).
- **Liu et al. (2020)** *eLife* — CCK interneurons gate hippocampal-to-PFC communication via endocannabinoid signaling.
- **Whissell et al. (2019)** *eNeuro* — Chemogenetic activation of CCK-GABA neurons enhances cognition broadly.

### SCZ risk genes directly affect CCKBCs
- **Kotzadimitriou et al. (2018)** *eNeuro* — NRG1 (SCZ risk gene) overexpression reduces NMDAR currents ~50% specifically in CCK interneurons.
- **del Pino et al. (2017)** — ErbB4 (NRG1 receptor, SCZ risk gene) is functionally critical for CCKBC wiring.

### 22q11.2 specifically
- No direct studies of CCKBCs in 22q models exist. This is a gap in the literature. However, Dgcr8 haploinsufficiency disrupts interneuron migration broadly, and the hippocampal circuit dysfunction in 22q models is consistent with potential CCKBC involvement.

## Options for proceeding

### Option A: Report specificity correlation only (current approach)
- VIP specificity is the only significant positive correlate (rho=0.47, P=0.032)
- CCK specificity is non-significant (rho=0.17, P=0.46)
- Consistent with bias pipeline; avoids the TPM discrepancy
- Risk: a careful reviewer who checks TPM will notice the difference

### Option B: Report both metrics with explanation
- Show specificity as primary (matches pipeline)
- Acknowledge TPM correlation for CCK in supplementary, explain that CCK's broad expression (94% of cell types) means TPM variation within CGE does not reflect subtype-specific enrichment
- More transparent, pre-empts reviewer concerns

### Option C: Reframe the analysis
- Instead of "which marker best predicts 22q bias," frame as "VIP expression gradient explains 22q bias variation across CGE subtypes"
- Show VIP scatter only (with cluster annotations showing VIP+/VIP- types)
- Move multi-marker comparison to supplement
- Avoids direct CCK comparison

### Option D: Nuance the CCKBC interpretation in the manuscript
- Soften "CCKBC activity was preserved" to note this is task-specific and region-specific
- Add a sentence acknowledging the literature on CCKBC involvement in SCZ (Lewis lab findings)
- Frame the mouse result as: "Within the hippocampal CA1 behavioral paradigms tested, ISI cells showed the most prominent functional deficits, while CCKBC activity patterns for locomotion and reward were relatively preserved. This does not exclude CCKBC involvement through other mechanisms or in other brain regions."

## Recommendation

**Combine Options B + D**: Use specificity as the primary metric with a brief note about TPM. Soften the CCKBC interpretation in the manuscript to acknowledge it's region- and task-specific, and that the broader literature implicates CCKBCs in SCZ via the CB1R/endocannabinoid axis (Lewis lab). The story remains that VIP lineage is the dominant axis of 22q vulnerability within CGE, while leaving room for CCKBC involvement through distinct mechanisms.

## Decision needed

- [ ] Which metric to use as primary for Figure: specificity vs TPM vs both?
- [ ] How to frame the CCKBC mouse finding: "preserved" vs "relatively preserved in tested paradigms"
- [ ] Whether to cite Lewis lab CCKBC-SCZ findings in Discussion
- [ ] Whether this analysis belongs in main text, supplement, or is dropped entirely
