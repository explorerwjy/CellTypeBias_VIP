**Revision To-Dos: (Summary)**

**Analyses:**

-   ID'ing CCKBC cells from ISIs based on our stim protocol. (waiting on
    raw data from Erica)

-   Add feature matching based null for cell type biases (Done)

-   Simulation of different cohort size and their impact on cell type
    implication. Main conclusion: ASD and DD are stable, even with 10%
    cohort size we can identify the correct cell types; SCZ seems not
    yet saturated, correlation with full set is bad when subsampled. In
    future SCZ needs more sequencing (Done)

-   Different gene set size, explain why SCZ has high correlation with
    main when adding low confident genes

-   Add negative control from non-brain traits (Done)

-   Add SCZ neg genes (Done)

-   Show CCKBC and other VIP (ISI types) has different VIP biases

    -   Option 1: use SNCG, CR and M2R select clusters from human atlas
        and label them as CCK or ISI. For highly expressed 22q genes and
        22q DEGs, CCK \< ISI with marginal significance. (Done)

    -   Option2: Mouse M1 dataset has some CCKBCs labelled, we can map
        Mouse patch-seq data to mouse cell atlas then to human atlas.
        Different human CGE clusters were mapped with these M1 patch-seq
        CCKBCs, and there's no difference between CCKBC and VIP 22q bias
        with this method. I'll spend more time fine tuning the mapping
        to see if they can converge to option 1. [[M1
        dataset]{.ul}](https://www.nature.com/articles/s41586-020-2907-3#Sec27).

-   Show why capping max Spec at 2. I simulate cell type specificity
    based on expression and ZINB noise profile, show smaller cell type
    and gene with lower expression will have inflated specificity if not
    properly controlled. Also tried different caps show results stay the
    same with different reasonable cutoffs (Done)

-   Provide all tables for geneweights/bias (Done)

-   Reorganize main figures to match our reply to Reviewer \#2.

**Wetlab experimental:**

-   New figure: Ephys characterization study (Erica & Steph)

-   New supplemental figure: 'CCKBC activity can't be analyzed by
    waveform analysis' evidence.

**Response Progress:**

R1

> General comment summary: 0/1 (update to 'completed' phrasing, update
> the text)
>
> R2
>
> Major comments: /5
>
> Minor comments: /5
>
> R3
>
> Major comments: /8
>
> Minor comments: /2

**PREAMBLE:** We thank the reviewers for their helpful and detailed
comments to improve our manuscript. To address the reviewers\' concerns,
we have performed additional experiments and analyses, added X
figures,\...

**Reviewer \#1**: In this paper the authors take a multivariate approach
to addressing the intersection of gene mutation to neuropsychiatric
disease. They begin by taking human data on ASD from the Spark
consortium and Scz data from the Schema consortium. Cross analyzing this
data with analysis of human cell type gene expression data, they
conclude a swath of cell types including eccentric MSNs, Lamp5-Lhx6
populations and chandelier neurons, CGE neurons, Upper and deep layer IT
neurons and finally Midbrain inhibitory cells are strongly affected by
the intersection of the genes implicated jointly from the Spark and
Schema data sets. In addition they identify several super clusters
including excitatory amygdala, hippocampus and both IT and ET
populations and say these findings are consistent to with such brain
disorders. In effect, they have implicated a large portion of the cell
types in both pallial and subpallial structures as being candidate for
being impacted in patients suffering from either ASD and SCZ. Going on
from there they then focus down on MGE and CGE-derived cells with a
particular focus on Lamp5-Lhx6/Chandelier, which comprise a tiny
subpopulation of the broader interneuron subtypes. From this more
general vantage, they switch their focus to the 22q11.2 deletion
disorder and use the Df(16)A+/- mouse model to study this in mice. Using
a largely motor based metric, they conclude that VIP interneurons are
disproportionately affected in their function, leading them to generally
conclude this implicates VIP interneurons as the underlying cause of
certain forms of ASD/SCZ.

Taken together this work provides a large survey of genes, cell types
and psychiatric conditions but feel more like a review article than a
focused piece attempting to link gene defects to an underlying cause.
While the work seems to be carefully assembled, the broad overview does
little to uncover a distinct subtype being convincingly implicated in
the disease rather than simply providing hints as to a large number of
cell types that may be implicated. As such, while aspects of this work
are of interest, they collectively short fall far of providing a
convincing narrative to implicate any cell type, the VIP one
investigated in their Df(16)A+/- mouse model in particular, as being
centrally involved in ASD and SCZ disorders. Given the wide range of
neuropsychiatric disorders they are grouping, this isn\'t terribly
surprising and does not argue strongly at all that this work is suitable
for publication in Neuron.

We sincerely thank the Reviewer for their thoughtful comments and for
the opportunity to clarify and sharpen the narrative and goals of our
study. Rare, highly penetrant variants are known to be strong drivers of
cognitive dysfunction across neuropsychiatric disorders. The recent
emergence of comprehensive human single-cell transcriptomic atlases from
the Allen Institute provides an unprecedented opportunity to investigate
the neuronal consequences of such rare mutations in a systematic and
unbiased manner. While prior studies mapping genetic risk to brain cell
types have largely focused on common-variant GWAS signals, the cell-type
specificity of rare, high-impact mutations has remained comparatively
unexplored.

To address this gap, we developed a computational framework to
systematically identify human neuronal cell types that are selectively
impacted by rare damaging mutations across neuropsychiatric phenotypes.
Applying this framework across multiple large genetic datasets
(including ASD from the SPARK consortium, schizophrenia from the SCHEMA
consortium, developmental delay cohorts, UK Biobank cognitive
phenotypes, and genes within the 22q11.2 deletion interval), we find
that although multiple neuron types are affected, caudal ganglionic
eminence (CGE)--derived GABAergic interneurons emerge as a consistent
and specific nexus of *cognitive phenotypes* across disorders. This
convergence distinguishes CGE-derived interneurons from both MGE-derived
interneurons and excitatory projection neurons, and identifies a
cross-disorder cellular axis of cognitive vulnerability that, to our
knowledge, has not been previously recognized.

We agree with the Reviewer that the original manuscript placed
disproportionate emphasis on atlas-level mapping and did not
sufficiently foreground this central result. In the revised manuscript,
we will substantially reorganized the Introduction and Results to
streamline the narrative and pivot earlier toward the central finding
that CGE-derived interneurons represent a shared locus of vulnerability
for cognitive dysfunction across disorders, while the broader
atlas-level analyses are repositioned to provide context and validation.

Within this CGE framework, our *in vivo* work focuses on CGE/VIP
interneuron populations because they provide a mechanistically
interpretable and experimentally accessible instantiation of the broader
CGE lineage signal. Although numerically smaller, these interneurons
play a disproportionate role in regulating cortical gain, contextual
modulation, and information routing through disinhibitory circuit
motifs. Importantly, the in vivo phenotypes we examine in the
*Df(16)A^+/−^* mouse model reflect alterations in intrinsic excitability
and in the encoding of spatial and reward-related information, rather
than motor function per se. In the revised manuscript we clarify this
point and reposition the mouse experiments as a circuit-level test of
how CGE-lineage interneuron dysfunction may translate rare mutation
burden into altered cognitive representations.

Taken together, the revised manuscript presents a more focused and
hypothesis-driven framework. By integrating large-scale human genetics,
single-cell transcriptomics, and circuit-level physiology, our study
identifies CGE-derived interneurons as a convergent cellular substrate
through which rare, high-impact mutations may disrupt cognitive function
across neuropsychiatric disorders. We believe this clarified narrative
more directly conveys the conceptual advance of the work and better
aligns the study with the scope and expectations of *Neuron*.

**Reviewer \#2**: Based on the hypothesis that cognitive deficits are
common features of psychiatric disorders, the authors investigate
whether there is a biased enrichment of genes identified as rare, highly
penetrant genetic mutations and specific neuronal cell types. They found
a correlation between cell type biases in Autistic Spectrum Disorder
(ASD) with low IQ and Schizophrenia (SCZ). There is a bias towards
particular cell types in ASD (striatal spiny neurons) and SCZ (MGE
interneurons). CGE interneurons appeared to be highly enriched in genes
associated with disorders involving cognitive deficits. They then
examined how these CGE interneurons behave in a mouse model of
intellectual disability/schizophrenia (modelling 22q11.2 deletion) and
showed that some interneurons from the VIP subclass exhibit deficits in
activity, suggesting a top-down disinhibition. Although the paper is
very interesting and most of the experiments are solid, some points need
to be addressed to strengthen the manuscript.

We thank the Reviewer for the kind words and for the careful and
constructive evaluation of our manuscript. As detailed below, the
Reviewer's comments helped us clarify the presentation and further
strengthen the mechanistic interpretation of our results.

Major points

1\) There are some inconsistencies regarding the comparisons the authors
provided; at least, the rationale for their absence is not included in
the manuscript\'s narrative. Figure 2: I am missing the comparison of
ASD with ID and ASD mutation bias for MGE interneurons and Medium spiny
neurons (at least for MGE, the equivalent of Figure 2G for MGE). I think
this is important since Figure 2 characterizes where the bias lies.
Figure 3, missing the comparison between SCZ and ASD with ID. Again, it
appears that the correlation is closer when SCZ and ASD with ID are
compared; therefore, why not include this comparison? Figure 5, the
authors would need to provide the calcium transients (frequency) for the
other CGE subtypes, at least in the Supplementary information (similar
to 5J).

2\) Figure 3: Sometimes, it is unclear what the authors intend to
convey. For example, it would be helpful to clarify the apparent
discrepancy of MGE being more biased in Schizophrenia (Figure 2),
considering Schizophrenia\'s stronger correlation to ASD with ID, and
still showing a greater bias difference with CGE cells than MGE cells
when comparing ASD with ID to ASD.

We apologize for the initial omission and will incorporate the full set
of requested comparisons to clarify and strengthen the interpretation of
our results. Specifically, when comparing ASD with SCZ, both MGE- and
CGE-derived interneurons exhibit stronger mutation biases in SCZ. When
comparing ASD with ASD with intellectual disability (ID), both MGE and
CGE interneurons show stronger biases in ASD with ID; however,
CGE-derived interneurons display substantially greater relative biases
(P-value = 1x10^-6^), highlighting a more pronounced CGE contribution to
ASD-associated cognitive phenotypes. Furthermore, when comparing SCZ to
ASD with ID, CGE interneurons show comparable mutation burdens (P-value
= 0.9), whereas MGE-derived interneurons are less biased in ASD with ID
compared to SCZ (P-value = 7x10^-3^), confirming a preferential
involvement of CGE interneurons (rather than MGE) in phenotypes with
prominent cognitive impairment, such as ASD with ID. In response to the
Reviewer's comments, these comparative results are now be explicitly
highlighted and discussed in the revised manuscript.

For Figure 5, we are not able to include additional calcium-transient
frequency analyses for other CGE-derived interneuron subtypes, because
currently available genetic targeting tools lack sufficient specificity
to distinguish CGE- versus MGE-derived subtype contributions (e.g.,
neurogliaform Lamp5⁺/Id2⁺ cells), or substantially overlap with
VIP-expressing populations (e.g., Sncg⁺ cells). However, our results do
demonstrate that among VIP+ interneurons, the most robust and consistent
functional alterations are observed in interneuron-selective interneuron
subtypes. We expanded the Discussion, to more explicitly address these
technical and interpretive limitations.

3\) It is clear from this paper that the activity of VIP interneurons is
altered. The authors state that this leads to disinhibition of pyramidal
cells because VIP inhibits Somatostatin (SOM) and Parvalbumin (PV),
resulting in impaired spatial coding. Although their model shows that
this causes changes in the percentage and configuration of place cells,
no experimental data support this claim. There is a paper from the team
on BioRxiv, where they conduct a comprehensive analysis of PV, SOM, and
a previous paper demonstrating that the dynamics of place cells are
impaired (Nature Neuroscience, 2017). However, if they wish to argue
that the changes observed in place cells in the place cell field are due
to disinhibition of PV and SOM, it would be valuable to include a rescue
experiment using optogenetics to activate the VIP population and then
observe the responses of pyramidal cells and PV/SOM interneurons.

While we agree that optogenetic rescue experiments activating VIP
interneurons would, in principle, provide a powerful test of causality,
such experiments fall beyond the scope of the present study and will be
clearly identified as an important direction for future work.
Importantly, the interpretability of such rescue experiments in the
context of a multi-gene 22q11 deletion is inherently limited:
optogenetic activation of VIP interneurons may not fully restore normal
circuit function when dozens of dosage-sensitive genes are
simultaneously affected. Moreover, PV and SOM interneurons are
themselves perturbed by the same deletion, further complicating causal
attribution and necessitating extensive additional control experiments.
Addressing these interacting circuit-level effects would require a
substantial experimental expansion beyond the focus of the current
manuscript. We explicitly discuss these limitations in the revised
version.

5\) In Figure 5L, there is a notable change in reward for several VIP
interneuron populations; this data is very intriguing. What is the
mechanistic hypothesis behind this? Do these VIP interneurons express
Dopamine receptors?

We appreciate the Reviewer\'s interest in the reward-related modulation
of VIP interneuron activity. Available transcriptomic datasets indicate
that multiple interneuron classes, including VIP populations, express
dopamine receptors at varying levels (Reviewer Figure 2.5), consistent
with prior reports (Anastasiades et al., 2019; Bae et al., 2025;
Glausier et al., 2009). However, VIP interneurons are also strongly
modulated by cholinergic and serotonergic inputs, and reward-related
activity likely reflects convergent neuromodulation rather than a purely
dopaminergic mechanism. In the revised Discussion, we note the evidence
for dopamine receptor expression and emphasize that defining the precise
neuromodulatory mechanism remains an important direction for future
study.

![A screenshot of a computer screen AI-generated content may be
incorrect.](media/image1.png){width="6.5in"
height="2.6534722222222222in"}

**Reviewer Figure 2.5: Dopamine receptor expression in human brain cell
atlas**. Expression (TPM) of the five dopamine receptor genes
(DRD1--DRD5) across 361 neuronal cell type clusters from the Siletti et
al. human brain atlas, grouped by 21 neuronal superclusters. Each dot
represents one cell type cluster; box plots show median and
interquartile range within each supercluster. CGE interneurons are
highlighted in red. Superclusters are ordered by the mean of per-gene
median expression across all five receptors.

6\) Could the authors compare their VIP interneuron subtype findings
with human sequencing data? This could help identify if one subtype is
more impacted in 22q11.2 deletion syndrome.

Follow up: CGE cells have been classified using transcriptomics in both
mouse and human studies, such as those from the Allen Brain Institute.
In mice, the transcriptome, morphology, and electrophysiological
characteristics have been analysed within the same cell. They could
attempt to establish a correlation based on the electrophysiological
properties already known for CGE interneurons in mice (e.g the
equivalence of some of these interneuron function in hippocampus). Then,
they could use mouse transcriptomics to create correlations in humans. I
understand that most of the available data has been gathered from the
neocortex, but many interneurons also share features with mirror neurons
in the hippocampus because they are all part of the cerebral cortex.

We thank the Reviewer for this excellent suggestion, which directly
motivated a substantial new analysis. Following the Reviewer\'s proposed
strategy, we used both cross-species transcriptomic mapping and
electrophysiological feature transfer to identify which human CGE
clusters correspond to functionally distinct mouse VIP subtypes and then
tested whether these imputed subtypes show differential 22q11.2 mutation
bias.

We integrated mouse patch-seq transcriptomic data from the Allen
Institute (5,764 cells from visual and motor cortex (Gouwens et al.,
2020; Scala et al., 2021), including 333 transcriptomically defined
CCKBCs from the Sncg and Vip subclass) with the human Siletti et al.
atlas (175,000 CGE interneurons across 21 clusters) using two
independent integration methods: Harmony batch correction and scVI
latent-space integration (scArches). In both methods, we projected mouse
cells into the human reference space and assigned each to its nearest
human cluster via k-nearest-neighbor voting (k = 30). Five human CGE
clusters (277, 278, 279, 280, 281) received a majority of CCKBC
assignments (Harmony CCKBC fraction \> 50%), with moderate concordance
between methods (Spearman rho = 0.46, P = 0.03). Among these, clusters
279, 280, and 281 showed the highest CCKBC fractions under both Harmony
(67\--96%) and scVI (86\--90%) and are annotated as VIP-expressing
(INT-VIP) in the Siletti atlas, while clusters 277 and 278 are
VIP-negative.

We additionally attempted to transfer electrophysiological signatures
across species using 16 matched feature pairs, defining a CCKBC
fingerprint from the 6 features that significantly discriminate mouse
CCKBCs from VIP-other interneurons (\*P\* \< 0.01): firing rate,
sustained firing, ISI CV, waveform adaptation, first-spike latency, and
early firing rate. Critically, k-nearest-neighbor classification (k = 5)
showed that the fraction of mouse CCKBC neighbors did not differ between
imputed CCKBC and VIP+ ISI human cells (Mann--Whitney \*P\* = 0.79;
Kruskal--Wallis three-group \*P\* = 0.96), confirming that the mouse
CCKBC electrophysiological signature does not robustly transfer to human
subtypes.

We found the five imputed CCKBC clusters segregate by VIP expression
into two groups with divergent mutation-bias patterns. VIP-negative
CCKBCs (clusters 277, 278) which most closely resemble mouse CCKBCs in
cross-species electrophysiology alignment, show low 22q11.2 mutation
bias. In contrast, VIP-positive CCKBCs (clusters 279, 280, 281) show
high 22q11.2 mutation bias comparable to VIP+ ISI clusters (Reviewer
Figure 2.6). Direct comparison of imputed CCKBCs vs. VIP+ ISI clusters
showed no significant difference across four 22q11.2 del gene sets
(Mann--Whitney P = 0.37).

These results demonstrate that VIP expression status, rather than CCKBC
identity, is the primary predictor of 22q11.2 mutation bias within CGE
interneurons. Consistent with prior work, cross-species classification
of Sncg subtypes (Scala et al., 2021), yields near-chance accuracy
(AUROC = 0.50; Bakken et al., 2021)), and CCK mRNA is expressed
ubiquitously across human neurons (Darmanis et al., 2015), which is also
what we observed in Siletti et al., (2023), precluding its use as a
discriminative marker. Given these constraints, we conclude that the
broader VIP-defined clusters, rather than more granular individual
subtypes, represents the most robust unit for cross-species integration
of mutation-bias signals. We present this analysis in the Supplementary
Information and discuss the biological limitations of cross-species VIP
subtype mapping in the revised Discussion.

![A screenshot of a computer screen AI-generated content may be
incorrect.](media/image2.png){width="6.5in"
height="2.951388888888889in"} **Reviewer Figure R2.6. Imputed human
CCKBC identity does not predict differential 22q11.2 mutation bias
within CGE interneurons.** Human CCKBC clusters were imputed by mapping
mouse M1 patch-seq data (including Sncg/CCKBC cells) onto the human
Siletti et al. CGE atlas using Harmony and scVI integration
(concordance: Spearman rho = 0.46, P = 0.03). Five clusters (277, 278,
279, 280, 281) with Harmony CCKBC fraction \\\> 50% were designated as
imputed CCKBCs. (Left) Three-way comparison of 22q deletion mutation
bias. VIP-negative CCKBCs (clusters 277, 278; green) show the lowest 22q
bias. VIP-positive CCKBCs (279, 280, 281; red) show bias comparable to
VIP+ ISI clusters (blue), indicating that VIP expression status rather
than CCKBC identity drives 22q bias within CGE interneurons. (Right)
CCKBC fraction (Harmony mapping) vs 22q deletion bias across all 16
VIP/CCKBC clusters. There is no correlation between CCKBC-likeness and
22q bias (Spearman rho = 0.05, P = 0.86), confirming that mouse CCKBC
identity does not predict 22q vulnerability at the subtype level.

Minor

1)  Figure 1 could be a nice figure for the Graphical Abstract, but not
    as a main figure of the manuscript.

Accordingly, we will use Figure 1 as the Graphical Abstract, and remove
it from the main text.

2)  Some sentences are overstated; the overstatements need to be toned
    down. E.g., Line 220 showed the strongest differential biases
    towards SCZ compared to ASD. Line 180: Is it appropriate to identify
    as \'strong correlation\' when R=0.68?

We agree with the Reviewer and will soften the language regarding
"strong correlations" (e.g., R = 0.68) and similar phrasing throughout
the manuscript.

3)  Figure 2D can go to the supplementary; still, it will make the point
    for specificity, while the paper seems to start to focus on other,
    more relevant cell populations.

We thank the reviewer for this suggestion. While we understand the
motivation to streamline the presentation, we believe that Figure 2D
represents an important result that helps establish the broader
landscape of neuronal vulnerability across disorders and therefore
should remain in the main text. Rather than moving the figure to the
Supplementary Information, we have revised the surrounding text to more
clearly focus the narrative and guide the reader toward the central
conclusion of the study. In particular, the Introduction and Results
have been reorganized to distinguish similarities and differences in
neuronal impact across disorders, while pivoting earlier to the key
observation that CGE-derived interneurons represent a shared locus of
vulnerability for cognitive dysfunction.

As illustrated in Figure 2C, cortical excitatory neurons show comparable
levels of mutation bias in both ASD and schizophrenia, whereas striatal
Medium Spiny Neurons (MSNs) display a more selective enrichment in ASD
relative to SCZ. We now clarify in the text that this divergence,
although biologically interesting, is not the primary focus of the
manuscript. Instead, the figure serves to contextualize the broader
analysis and highlight the contrast with the more consistent
cross-disorder signal observed in CGE-derived interneurons. By retaining
the figure while sharpening the narrative around it, we believe the
manuscript more clearly conveys both the overall cellular landscape and
the specific convergence that motivates our subsequent functional
experiments.

4)  Please add the \"n\" of mice in the experiments.

We thank the Reviewer for the comment. We added \"n\" numbers for all
mouse experiments to the figure legends.

5)  Feel like a jump from Figure 5 to 6, can explain the rationale
    behind it more clearly.

We thank the Reviewer for the comment. We have now clarified the
rationale for the transition between Figure 5 (activity) and Figure 6
(spatial coding).

**Reviewer \#3:** Vitkup et al. present a framework that integrates rare
genetic variant association data with human single-cell transcriptomic
profiles to infer which brain cell types are most impacted in ASD,
schizophrenia, and related cognitive phenotypes. They further test
predictions in a 22q11.2 deletion mouse model, focusing on VIP+ CGE
interneurons using state of the art techniques. The manuscript is
ambitious, well-written, and addresses a timely question. However, I
have several major concerns about the validation of the computational
framework, the genetic assumptions underlying their analyses, and the
interpretation of results.

This type of manuscript is hard to write as a scientist and hard to
interpret as a reader. The reason is the necessary narrowing of scope
(global cell type analysis to functional interrogation of a single cell
class) always leaves you to wonder if the central hypothesis of the
first part has really been tested or simply one out of several
possibilities. I think the authors do a decent job of this hard
transition, and importantly, I do not want to ask for major experimental
work. This means I urge the authors to be careful in summarizing the
take home message. Therefore, I think the discussion needs a paragraph
on non-VIP specific effects in the 22q11 mouse. What would happen if
they recorded from SST, PV, SNCG or LAMP5 interneurons instead? Could
the change of VIP function be related to changes in excitatory inputs
(are they enriched)? The reduced correlation with movement would
certainly suggest that they receive weaker inputs.

Beyond this I have focused my review on the computational parts since
this reflects my expertise.

We thank the Reviewer for the careful, balanced, and insightful
evaluation of our manuscript. We appreciate the recognition of the
ambition and technical rigor of the work, as well as the thoughtful
guidance on how to more clearly articulate the scope and interpretation
of our findings. We fully agree that manuscripts bridging large-scale
computational inference with targeted circuit-level interrogation must
be especially careful in framing their conclusions. In the revision, we
expanded the discussion to explore non-VIP effects in circuit, including
the results of our previous work examining other interneuron subtypes.
Below, we address each of the Reviewer's points in turn.

Major points:

1\. Validation of the computational pipeline

The permutation framework needs stronger justification. Matching only on
gene set size is insufficient; at minimum, random gene sets should be
matched on gene length, conservation, and expression level. The authors
should more comprehensively validate whether their \"mutation bias\"
statistic is robust to such gene-level confounds.

We thank the reviewer for this insightful comment. We agree that it is
important to carefully justify the permutation framework and to evaluate
potential gene-level confounds.

First, several major confounders were already explicitly incorporated in
our pipeline. For schizophrenia, the expected mutation counts used in
the burden analyses already accounted for gene length, and we have now
extended this expected-count normalization to the ASD and NDD gene sets
as well. Expression level is controlled at two levels: (i) TPM-based
normalization across cell types (Equation 1) and (ii) gene-level
specificity scores that normalize expression across all profiled cell
types (Equation 2). Gene constraint and conservation were not treated as
independent nuisance variables, as constraint scores were strongly
correlated with cross-disorder mutation bias in our data (as shown in
Figure 2B), indicating that they capture a biologically meaningful
component of disease architecture.

To directly address the reviewer's suggestion, we have now implemented
an alternative permutation framework in which random gene sets are
matched to the disorder genes on gene length (CDS bases), whole-brain
expression level, and evolutionary conservation (mean PhastCons score).
Specifically, for each null simulation we draw 1,000 candidate random
gene sets and select the one whose mean percentile values across these
three covariates are closest (in Euclidean distance) to the real gene
set. This \'best-of-1,000\' procedure was repeated 10,000 times per gene
set. This alternative null model tested whether genes with similar basic
genomic and expression properties, and with the same mutational weights,
recapitulate the observed cell-type mutation biases.

As expected, this more stringent matching reduced the nominal
significance of several traits, likely because many genes with matched
LOEUF and brain expression are themselves enriched for CNS-related
functions (Reviewer Figure 3.1). However, the key comparative
conclusions of the paper remained unchanged: the relative ranking of
cell types, the convergence on CGE interneurons for disorders with
cognitive impairment, and the cross-disorder contrasts were consistent
across both the original and matched null models.

We now describe both null models and their rationale in the Methods
section, and report FDR values under each model for all traits in a new
supplementary table (Supplementary Table X). In addition, as a further
robustness check suggested by reviewer, we added multiple non-brain
negative control traits. As shown in the Reviewer Figure 3.1, the
negative control gene sets follow the null distribution closely, with no
neuronal cell types reaching significance, confirming that our bias
calculation does not produce spurious bias for traits without expected
neuronal specificity (see also response to point 3 and new Supplementary
Figure X).

![](media/image3.png){width="6.5in" height="3.9305555555555554in"}

**Reviewer Figure 3.1: QQ plots of cell-type bias p-values under random
and confounder-matched null models.** Each panel shows a
quantile-quantile plot of observed versus expected -log(p-value) for
cell-type bias across 461 cell types. P-values were derived from 10,000
permutations under two null models: random gene sampling (gray) and
Best-of-N matching (blue; N = 1,000 candidate sets matched on gene
length, whole-brain expression, and evolutionary conservation). The top
row displays psychiatric and neurodevelopmental disorder gene sets (SCZ,
ASD, and DDD), and the bottom row includes 22q11.2 deletion syndrome
alongside two non-brain negative control traits (HDL cholesterol and
alanine aminotransferase, italic labels). For all disorder gene sets,
p-values deviate strongly above the diagonal, indicating robust
cell-type enrichment that persists even after confounder matching. By
contrast, negative control gene sets fall along the diagonal, confirming
that the analysis does not produce spurious enrichment for traits
without expected neuronal specificity. Dashed line indicates the null
expectation.

2\. Gene set size and stability of estimates

In Figure S4, mutation bias correlations remain stable even when the SCZ
gene set is expanded to 200 genes. This is surprising: including genes
without true association should add noise. Why is this not the case? The
statistical explanation for this stability needs to be clarified. More
generally, the manuscript should explicitly address how the number of
included genes affects the stability of their estimates. Why should the
number of genes be selected based on number of ASD rare-variant
implicated genes?

We thank the Reviewer for raising this important point and agree that
the stability of mutation-bias correlations as the schizophrenia gene
set expands requires clarification. In the range we examine (up to about
200 genes), many of the additional genes do not represent pure noise,
but instead have weaker, yet non-zero, evidence for association. This is
particularly true for the NDD gene set, where hundreds of genes reach
exome-wide significance. For schizophrenia, while the association
strength declines beyond the top-ranked genes, many added genes remain
enriched for constrained, brain-expressed loci, thereby contributing
structured signal rather than random noise.

To directly test whether added genes carry signal or noise, we compared
the bias-profile stability when appending "real" ranked disease genes
(ranks 62--N) versus appending "random" genes to the fixed top-61 core.
For each gene-set size (N = 61 to 201, in steps of 10), we computed the
Spearman correlation between the expanded set\'s 461-cell-type bias
profile and the top-61 reference; for random additions, we repeated this
100 times and report the mean and 95% confidence interval. Across all
three disorders, real ranked genes maintained higher correlation with
the top-61 reference than random gene additions (Reviewer Figure R3.2c).
For ASD and DDD, the real-gene correlation remained above ρ = 0.96 even
at N = 201, consistently above the random 95% CI upper bound (0.962 and
0.972, respectively). For SCZ, the real-gene line (ρ = 0.80 at N = 201)
stayed above the random mean (ρ = 0.72) with increasing separation
beyond N ≈ 130, though the wider random CI reflects the noisier
case--control design and power limitation (which will be further
discussed in next response). These results confirm that genes beyond
rank 61 contribute structured, biologically coherent signal rather than
pure noise, and that this holds across all three disorders and genetic
architectures.

To clarify this point, we have now expanded the Methods sections to
explicitly discuss how gene set size and genetic architecture jointly
influence the robustness of mutation-bias estimates (see also the new
power analyses described in response to point 3). We selected the ASD
rare-variant gene set as a practical reference point for determining the
default gene count because it is sufficiently large to reveal robust
convergence on specific cell types, yet not so large as to be dominated
by weakly associated genes and random genes. Figure S4 is now explicitly
framed as demonstrating that the main conclusions are not specific to a
single gene-set size (e.g., N = 61) but remain stable across a broad and
biologically reasonable range.

![A screen shot of a computer screen AI-generated content may be
incorrect.](media/image4.png){width="6.5in"
height="1.867361111111111in"}

**Reviewer Figure 3.2: Expanding gene sets: real ranked genes vs.
weight-transferred random additions.** For each disorder, the top 61
genes define a reference bias profile across 461 cell types. Gene sets
are expanded from 61 to 201 genes in two ways: (1) appending the
next-ranked real disease genes (colored lines), or (2) replacing those
genes with randomly selected genes that inherit the rank-position
weights of the real genes they replace (grey lines; mean and 95 percent
CI from 100 permutations). For ASD with ID (center) and DDD (right),
real gene correlation remains above rho = 0.96 at N = 201, consistently
at or above the upper bound of the random 95 percent CI. For SCZ (left),
real gene correlation (rho = 0.80) stays well above the random mean (rho
= 0.68), with the gap widening beyond N approx 100. These results
demonstrate that ranked disease genes beyond the top 61 carry
structured, concordant cell-type signal rather than pure noise, and that
this holds across both de novo (ASD, DDD) and case-control (SCZ) genetic
architectures.

3\. Impact of genetic architecture

The authors should explain how differences in genetic architecture (rare
de novo-driven in ASD vs. modest case-control burden in SCZ) impact the
statistical properties of the mutation bias metric. What is the minimum
signal needed for robust inference? Is there a minimum sample size of
the rare-variant genetic investigation at which the mutational bias is
stable? Consider including a set of negative control traits, such as
blood traits or other rare-variant investigations into non-brain related
traits.

We thank the Reviewer for this important question regarding the impact
of genetic architecture on the statistical behavior of the mutation-bias
metric. To directly address this issue, we performed empirical down
sampling studies on all three disorder cohorts, subsampling probands /
cases and raw mutation counts they carry at fractions of 10--100% of the
full dataset, re-running gene discovery on each subsample, and computing
cell-type bias for the resulting top-61 genes across 100 bootstrap
iterations per fraction (see Reviewer figure R3.3A). This procedure
directly tests whether the mutation-bias metric is stable under
realistic reductions in sample size and, by extension, under the
differing statistical power regimes that characterize de novo versus
case-control rare-variant designs.

The results reveal a marked difference in stability between de
novo-driven and case-control designs. For ASD, even when downsampling
the cohort to just 10% of the full size (\~4,300 probands), the
cell-type bias correlation remained high (r = 0.86 +/- 0.06), despite
the fact that 51% of the top-61 genes differed from the full gene set.
For DDD, stability was even higher (r = 0.95 +/- 0.02 at 10%, \~3,100
probands), with \~68% gene overlap. This stability reflects the fact
that risk genes show highly correlated cell-type expression profiles and
that gene-level weights preferentially retain the strongest signals
across subsamplings. In contrast, the schizophrenia (SCZ) analysis,
based on a case-control design, showed markedly reduced stability at
lower effective sample sizes (r = 0.30 +/- 0.24 at 10%, reaching r =
0.83 only at 90% of the full cohort). The reduced stability in SCZ
reflects both the inherently lower statistical power of case-control
versus de novo designs and the substantially smaller per-gene effect
sizes in schizophrenia compared to ASD and NDD (Nguyen et al., 2017).
These results suggest that the current SCZ cohort size (\~24,000 cases)
is near the minimum required for stable cell-type inference, and that
future larger case-control sequencing studies will be needed to refine
SCZ cell-type conclusions.

In addition, to further evaluate specificity, we applied the full
pipeline to non-brain, negative control, traits derived from
rare-variant burden analyses in UK Biobank (e.g., HDL cholesterol,
inflammatory bowel disease, HbA1c) alongside neurological traits such as
Parkinson's and Alzheimer's disease. These control traits show distinct
and generally weaker brain cell-type mutation-bias profiles compared to
ASD, schizophrenia, and NDD, both in terms of the magnitude and the
anatomical distribution of biases (Reviewer figure R3.3B, new
Supplementary Figure X). Together, these analyses demonstrate that the
mutation-bias framework behaves as expected under differing genetic
architectures and that the prominent CGE and related interneuron signals
observed in cognitive disorders are not generic features of rare-variant
traits.

We expanded the Methods and Supplementary Note to explicitly discuss the
dependence of inference stability on genetic architecture, sample size,
and study design, thereby clarifying the conditions under which robust
cell-type conclusions can be drawn.

![](media/image5.png){width="6.5in" height="2.125in"}**Reviewer Figure
3.3A. Cell-type mutation bias estimates as a function of cohort sample
size.** Downsampling analysis showing stability of cell-type bias across
three disorders as a function of genetic sample size. For each disorder,
raw mutation counts were subsampled at fractions of 10--100% of the full
cohort using binomial sampling, gene discovery was re-performed on the
subsampled data (Poisson test for de novo cohorts, Fisher exact test for
case-control), and the top 61 genes were selected to compute cell-type
bias across 461 cell types. This procedure was repeated for 100
independent iterations at each fraction. **Left y-axis (colored
lines):** Pearson correlation between the 461-dimensional bias vector
from subsampled data and the full-data bias vector. **Right y-axis (gray
dashed lines):** fraction of top-61 genes overlapping with the full-data
gene set. Lines and ribbons indicate mean ± s.d. across 100 iterations.
Horizontal dashed line marks r = 0.9. **(Left)** ASD de novo mutations
(N = 42,607 probands): bias correlation exceeds 0.9 even at 10% sampling
(r = 0.86), despite only \~49% gene overlap, indicating that different
subsets of ASD risk genes converge on similar cell-type expression
patterns. **(Center)** SCZ case-control mutations (N_case = 24,248):
correlation improves more gradually with sample size (r = 0.83 at 90%),
consistent with the noisier case-control design and lower gene overlap
(\~12% at 10%). **(Right)** DDD de novo mutations (N = 31,058 probands):
highest stability, with correlation exceeding 0.9 at all fractions
tested and \~68% gene overlap at 10% sampling.

![A graph with colorful dots AI-generated content may be
incorrect.](media/image6.png){width="6.5in"
height="2.2222222222222223in"}

**Reviewer figure R3.3B. Cell-type mutation bias profiles for non-brain
negative control traits.** For each trait, we selected the top 61 genes
by rare-variant LoF burden p-value from UK Biobank exome-wide
association analyses and computed mutation bias across 461 human cell
types using the same pipeline applied to the psychiatric disorders.
X-axis show statistical significance (−log₁₀ P-value) from permutation
testing against 10,000 matched null gene sets, with a vertical dashed
line indicating the FDR \< 0.1 threshold. (Left) HDL cholesterol:
top-ranked cell types are microglia, astrocytes, and vascular cells;
Migroglia reach FDR \< 0.1. (Center) Inflammatory bowel disease (IBD):
no cell type reaches FDR \< 0.1; (Right) Alanine aminotransferase: no
cell type reaches FDR \< 0.1.

4\. Treatment of cognitive phenotypes (VNR)

The inclusion of verbal-numerical reasoning (VNR) is a good idea.
However, splitting genes by positive vs. negative effect sizes is
misleading. Due to negative selection, large-effect positive alleles on
cognition are extremely rare. Thus, VNR(+) is likely a \"null\" set
without signal. The authors do not split ASD or SCZ by direction, so why
should VNR be split?

We thank the reviewer for this thoughtful comment and agree that the
interpretation of the VNR(+) gene set requires careful consideration,
given the strong negative selection acting on cognition-enhancing
alleles.

VNR is a quantitative trait with a defined population mean. In our
rare-variant burden framework, VNR(+) genes are those for which rare
coding mutations are associated with higher VNR scores on average,
whereas VNR(−) genes are associated with lower scores. As the Reviewer
correctly points out, large-effect rare alleles that truly increase
cognition are expected to be extremely rare; therefore, the VNR(+) gene
set likely reflects a mixture of near-neutral and mildly protective
effects rather than strong cognition-enhancing variants. We therefore
used VNR(+) primarily as an internal control: under the hypothesis that
CGE interneurons are preferentially impacted by variants that impair
cognitive performance, one would expect VNR(+) genes to show reduced
bias toward CGE interneurons. Consistent with this expectation, VNR(−)
genes exhibit a strong positive CGE mutation bias, whereas VNR(+) genes
show a significant depletion in CGE interneurons (Figure X).

For ASD and NDD de novo datasets, the study design does not permit a
directional split into risk-increasing versus protective variants. For
schizophrenia, however, a directional split is in principle possible. To
explore this, we identified the subset of genes with the strongest
nominal "protective" direction of effect (more mutations in controls
than in cases, OR \< 1) in the SCHEMA dataset (n = 61) and computed
their mutation-bias profile. As expected for putative protective or
neutral genes, these genes exhibit near-zero or weak biases across other
neuronal classes, with positive bias toward non neuronal cell types (new
Supplementary Figure X). Interestingly, CGE IN showed the strongest
depleted bias among these genes.

We now clarify in the text that the VNR(+) set was included primarily as
a control analysis whereas the principal biological conclusions are
based on VNR(−) and the psychiatric disorder gene sets. We also
incorporated the analysis of SCZ genes with OR \< 1 into the relevant
Results section to further illustrate this point.

![A screen shot of a graph AI-generated content may be
incorrect.](media/image7.png){width="6.5in"
height="3.1791666666666667in"}

**Reviewer Figure R3.4. Cell-type mutation bias profile for
schizophrenia genes with protective direction of effect (OR \< 1).** We
identified top 61 genes from the SCHEMA exome sequencing study for which
rare damaging mutations were more frequent in controls than in cases
(odds ratio \< 1), representing putative protective or neutral variants.
(Left) Mutation bias (weighted mean specificity) by supercluster, sorted
by mean. Non-neuronal cell types (ependymal, astrocyte, vascular) show
the highest positive bias, while neuronal classes, particularly CGE
interneurons (mean bias = −0.183, Z = −2.63), show the strongest
depletion. (Right) Statistical significance (−log₁₀ P-value) from
permutation testing against 10,000 matched null gene sets. No cell type
reaches FDR \< 0.1, consistent with this gene set representing a mixture
of near-neutral and mildly protective effects rather than a coherent
biological signal. The inverted pattern relative to SCZ risk genes
(which show strong CGE enrichment) provides directional validation:
genes enriched in cases implicate CGE interneurons, whereas genes
enriched in controls show CGE depletion.

5\. Interpretation of differential bias tests

The authors highlight CGE interneurons as differentially biased across
disorders, but many comparisons rely on small numbers of significantly
associated clusters. For example, in ASD without ID, only 7 cell types
pass significance, none of which are CGE. How do the authors justify
making differential claims when one side of the comparison lacks clear
signal?

We apologize for not making this distinction clearer in the original
manuscript. Our framework addresses two related but conceptually
distinct questions:

First, within-disorder bias asks whether a given cell supercluster is
significantly biased for a particular disorder compared to random gene
sets with the same mutation load.

Second, between-disorder contrast asks, given two disorders with
well-characterized high-penetrance gene sets, which cell superclusters
show the largest difference in mutation bias between them?

The significance thresholds reported in Figures 2 and 3 pertain to the
first analysis (within-disorder enrichment) In contrast, the
differential bias analyses (e.g., contrasting ASD with ID vs. ASD
without ID, or ASD vs. schizophrenia) address the second question by
evaluating differences in the estimated mutation-bias values themselves,
independent of whether each disorder's bias in that cell type
individually reaches a significance threshold.

In this context, the fact that CGE interneurons are not significantly
biased in ASD without ID is itself informative: CGE interneurons show
strong positive bias for disorders characterized by cognitive impairment
(ASD with ID, NDD, schizophrenia), but much weaker or absent bias in ASD
without ID. The cross-disorder contrasts therefore highlight CGE
interneurons as a cell population whose relative involvement tracks the
presence of cognitive deficits across disorders.

To avoid potential confusion, we revised the relevant sections of the
Results and Methods to more clearly distinguish these two analytical
frameworks and to explicitly clarify that CGE interneurons are not
significantly enriched in ASD without ID in the within-disorder
analysis.

6\. 22q11.2 analyses

In Figure S9, it is unclear which cell types remain significant after
multiple-testing correction. The authors should explicitly mark
FDR-adjusted significance. How about excitatory neurons? Why are they
not included in the enrichment analysis to choose which cell type to
focus on? It should be possible to connect the interneuron subtypes
found in the electrophysiology dataset to interneurons found associated
within the mutational bias framework. Is the difference between CCKBC
interneurons and ISI2, ISI3 and VIPo interneurons recapitulated in the
mutation bias analysis? This would strengthen the claims of convergence
between electrophysiological data and genetic data.

We thank the Reviewer for this comment. We agree that the statistical
interpretation of Fig. S9 required explicit clarification. In the
revised figure and legend, we now clearly indicate that none of the
individual clusters passed FDR correction for the 22q11.2 gene set. This
absence of FDR-significant clusters is expected given the small size and
heterogeneous pathogenicity of the 22q11.2 interval. As demonstrated in
multiple studies, only a subset of the 42 genes contributes to disease
risk, resulting in substantially lower signal compared to ASD, SCZ, or
NDD gene sets, in which each gene carries strong statistical evidence of
association. Accordingly, the 22q11.2 analysis is inherently
underpowered for genome-wide inference, and we now explicitly note this
in the manuscript.

Because of this limitation, our goal in the 22q11.2 section was *not* to
identify new cell types de novo, but rather to test a genetically
motivated, cross-dataset hypothesis: that CGE/VIP interneurons,
implicated across ASD, SCZ, and NDD rare-variant analyses, are also
among the most affected populations within the 22q11.2 locus. The in
vivo imaging data were therefore designed as a focused experimental test
of this hypothesis rather than as an unbiased screen.

The reviewer is correct that excitatory neurons show elevated
mutation-bias scores in the 22q11.2 gene set. However, these cell types
also show comparably high bias across all disorder gene sets we
examined, including ASD, SCZ, NDD, and even ASD without ID, suggesting
that excitatory neurons represent a broad and non-specific vulnerability
across many neurodevelopmental conditions. This widespread enrichment
makes excitatory neurons less informative for distinguishing
disorder-specific or phenotype-specific mechanisms. In contrast, CGE/VIP
interneurons exhibit a more selective pattern, showing strong
involvement in disorders with cognitive impairment (ASD with ID, SCZ,
NDD, VNR--), but show minimal or absent bias in ASD without ID. This
selective profile motivated our focus on CGE/VIP populations in the
functional experiments, which is further justified by the relative lack
of in vivo data on interneuron contributions to cognitive deficits,
compared to much more extensively studied excitatory circuits. We have
clarified this rationale in the revised Results and Discussion.

Regarding the Reviewer\'s specific question about whether the
electrophysiological distinction between CCKBC and ISI/VIPo interneurons
is recapitulated in the human mutation-bias analysis: As detailed in our
response to Reviewer 2, point 6, we performed an extensive cross-species
analysis using transcriptomic mapping (including Harmony and scVI
integration of mouse patch-seq data onto the human Siletti et al. atlas)
to determine whether the CCKBC vs. ISI distinction is recapitulated in
human mutation-bias profiles. We also tested whether the mouse CCKBC
electrophysiological signature transfers to the imputed human subtypes,
but found that this signature does not robustly transfer across species.
The key result is that VIP expression status, rather than CCKBC
identity, emerges as the primary predictor of 22q11.2 mutation bias
within CGE interneurons: VIP-negative CCKBCs show low 22q11.2 bias while
VIP-positive CCKBCs show mutation bias comparable to VIP+ ISI clusters
(Reviewer Figure R2.6; see also Supplementary Note 2 for full methods).
These results further support our focus on supercluster-level analyses
and suggest that rare-variant convergence on CGE interneurons reflects a
lineage-level, rather than subtype-specific, vulnerability.

7\. Splitting ASD by intellectual disability (ID). While the motivation
for splitting ASD cases into \"with ID\" and \"without ID\" is
reasonable, the labeling should be made explicit --- i.e., \"ASD with
ID\" and \"ASD without ID\" --- to clarify that this categorization was
imposed by the authors rather than inherent to the source datasets. To
strengthen the cross-disorder comparisons, a similar treatment should be
applied to schizophrenia, which is also associated with variation in
premorbid IQ. If individual-level data are not available to stratify SCZ
cases by cognitive performance, then at minimum the authors should
remove DD/ID genes from the schizophrenia gene list. Otherwise, the
contrast between \"ASD with ID\" and SCZ risks being confounded by
differences in how cognitive impairment is represented in the gene sets.

We thank the reviewer for this thoughtful and important comment. We
agree that the ASD subgroup labels need to be made explicit, and we now
consistently use "ASD with ID" and "ASD without ID" throughout the
manuscript. We also clarify in the Methods that this categorization is
defined by our analysis rather than being inherent to the original
datasets.

To directly address the reviewer's concern that the signal in
schizophrenia might be driven by overlapping DD/ID genes, we performed a
overlap and removal analysis. First, we quantified the overlap between
our SCZ gene set (61 genes) and established DD/ID (NDD) gene lists. We
found that the overlap is remarkably small: only 1 gene (GRIN2A)
overlapped with the top 61 NDD genes, and only 8 genes overlapped when
comparing against a broader list of 285 NDD genes (genome-wide
significant genes defined in (Kaplanis et al., 2020)). We then repeated
the mutation-bias analysis for schizophrenia after removing these
overlapping genes, first removing the 1 NDD61-overlapping gene, and then
removing the 8 NDD285-overlapping genes. In both cases, the CGE
interneuron bias remained robust and virtually unchanged. The overall
cell-type bias profile was highly stable (Spearman rho = 0.923 between
the full and NDD285-removed profiles across all 461 cell types; Reviewer
Figure R3.7A). Furthermore, direct comparison of CGE interneuron bias
across ASD without ID, SCZ, and the two NDD-removed SCZ gene sets
confirmed that the SCZ--ASD w/o ID difference persists after NDD gene
removal, while the SCZ profiles before and after removal are
statistically indistinguishable (Reviewer Figure R3.7B). This result
refutes the hypothesis that the SCZ signal is merely driven by
confounding DD/ID genes. Instead, it demonstrates that risk genes
specific for schizophrenia, which is clinically associated with specific
cognitive deficits distinct from severe ID, independently converge on
CGE interneurons. This finding reinforces our broader conclusion that
distinct genetic architectures (SCZ vs. ASD with ID vs. NDD) converge on
CGE interneurons as a shared cellular substrate for cognitive
dysfunction.

![A screenshot of a video game AI-generated content may be
incorrect.](media/image8.png){width="6.5in" height="2.7375in"}

**Reviewer Figure R3.7: Interneuron bias in schizophrenia is not driven
by NDD gene overlap**. (A) Scatter plot of cell-type mutation bias (461
cell types) for the full SCZ gene set vs. SCZ after removing 8
NDD285-overlapping genes, showing high stability (Spearman rho = 0.923).
CGE interneurons highlighted in red, MGE interneurons in blue. (B) CGE
interneuron bias comparison across ASD without ID, SCZ (61 genes), SCZ
rm NDD61 (60 genes), and SCZ rm NDD285 (52 genes). Mann--Whitney U test
P-values are shown for each pairwise comparison; SCZ profiles before and
after NDD gene removal are not significantly different.

8\. VIP function in the 22q11.2 mouse model

The decrease in correlation with movement is a strong phenotype that is
also subtype specific. I wonder if this is due to difference in ACh
receptors or generally decreased inputs. One major outstanding question
in the manuscript is (as I mention above) if the changes seen in vivo
are due to intrinsic changes in VIP cells (as predicted by computational
part) or a difference in network integration. Patch clamp recordings
from VIP neurons with input stimulation (to test pre-synaptic (paired
pulses) v.s. postsynaptic function (EPSCs or similar) would help resolve
this issue and thereby test one of the main claims. If this can be done
then I think it would greatly enhance the manuscript but I would be OK
with addressing this point clearly in the interpretation of the results.

We thank the Reviewer for this insightful comment. We agree that a key
unresolved issue is whether the observed reduction in movement-related
activity in VIP interneurons reflects intrinsic cellular dysfunction or
altered network integration and afferent input, including potential
differences in cholinergic signaling versus a more general reduction in
excitatory drive.

To directly address this question, we have performed whole-cell
patch-clamp recordings from VIP interneurons in the 22q11.2 deletion
mouse model. These experiments included XXXXX

In the revised manuscript, we integrate these electrophysiological data
with the in vivo calcium imaging results to refine the mechanistic
interpretation. This allows us to more clearly determine the extent to
which VIP interneuron dysfunction arises from cell-intrinsic alterations
versus disrupted circuit integration, thereby directly addressing one of
the central questions raised by the Reviewer.

Minor Concerns

1)  Please provide full mutation bias estimates for all traits tested,
    not just selected highlights. Preferably provide all mutational bias
    results for all traits tested in a \"long\" format, with additional
    column describing the tested trait. Gene lists and mutation counts
    used for each trait should be published in supplementary tables.

We thank the reviewer for this helpful suggestion. In the revised
manuscript, we now provide complete mutation-bias estimates for all
tested traits, rather than only the highlighted examples shown in the
main figures. These are included as Supplementary Tables X--Y, in a
long-format structure in which each row corresponds to a (trait × cell
type) pair. For each trait, we additionally supply the full gene list,
gene-level mutation counts, and expected mutation burdens used to
compute the mutation-bias statistics. These additions ensure full
transparency and facilitate reproducibility of all analyses.

2)  Clarify why the cell-type specificity score is capped at 2. How
    sensitive are results to this capping, and why not use previously
    established specificity metrics?

We thank the reviewer for raising this important point. Expression
specificity is computed as fold-enrichment:
$S\left( g,\text{ct} \right) = \text{TPM}\left( g,\text{ct} \right)\text{/}\overline{\text{TPM}\left( g \right)}$,
where values greater than 1 indicate above-average expression in a given
cell type. We cap this score at 2x the global mean specificity as a
conservative guard against technical inflation in cell types with low
sequencing depth. In our atlas of 461 cell types, uncapped specificity
values can reach up to \~97x, driven by small cell populations with low
total UMI counts where sampling noise inflates fold-enrichment
estimates. We confirmed this mechanistically via negative binomial
simulation: genes with uniform true expression (specificity of 1) across
all cell types show extreme artifactual specificity (up to 20x) in
low-UMI clusters, solely due to sampling noise (Supplementary Note).

To confirm the robustness of our choice, we recomputed mutation bias for
SCZ, ASD without ID, and DDD across different cap levels (1x, 1.5x, 2x,
2.5x, 3x, 5x, 10x). Cell-type-level Spearman rank correlations between
cap = 2x and caps in the 1.5\--3x range exceed 0.96 for all three
disorders (Reviewer Figure MC2), confirming high stability. At cap = 1x,
over-clipping compresses specificity variation; at higher caps (5-10x),
non-neuronal artifacts dominate\-\--Microglia climbs from rank 26 to
rank 3, overtaking genuine neuronal signals. The 2x cap thus sits in the
optimal range: preserving discriminability among neuronal subtypes while
preventing technical inflation. We also compared against the TDEP
specificity metric (Yao et al., 2025) and performed leave-one-out
stability analyses, both of which confirm that capping is necessary when
specificity is used as a continuous weight rather than for rank-based
thresholding. Full details of the NB simulation, TDEP comparison, and
LOO analyses are provided in Supplementary Note.

![A graph of different colored lines AI-generated content may be
incorrect.](media/image9.png){width="6.5in"
height="1.8277777777777777in"}

**Reviewer Figure MC2. Mutation bias rankings are robust to specificity
cap level in the different range; higher caps allow non-neuronal
artifacts to dominate.** (A) Spearman correlation (cell-type level)
between cap = 2× and each alternative cap level (1×, 1.5×, 2.5×, 3×, 5×,
10×), for SCZ, ASD without ID, and DDD. All disorders show Spearman R \>
0.96 in the 1.5--3× range (blue shading); at 1× correlations drop
sharply due to over-compression, while at 5--10× ASD rankings diverge
substantially. (B) Mean supercluster bias across cap levels for SCZ,
showing neuronal (solid) and non-neuronal (dashed) superclusters.
Neuronal types (CGE IN, MGE IN, MSN) plateau by 2×, while non-neuronal
types (Vascular, Microglia, Astrocyte, Ependymal) rise steeply from 1×
to 10×, progressively approaching neuronal bias levels. (C) Supercluster
ranking (1 = highest bias) across cap levels for SCZ. At 2×, neuronal
types hold the top ranks (CGE rank 2, MGE rank 5) while non-neuronal
types rank 20--31. At 10×, Microglia climbs to rank 3 and Vascular to
rank 5, overtaking MSN (which drops from rank 10 to 25), demonstrating
that high caps allow technical artifacts in low-UMI non-neuronal
populations to displace genuine neuronal signal.

Anastasiades, P. G., Boada, C., & Carter, A. G. (2019).
Cell-type-specific D1 dopamine receptor modulation of projection neurons
and interneurons in the prefrontal cortex. *Cerebral Cortex (New York,
N.Y.: 1991)*, *29*(7), 3224--3242.
https://doi.org/10.1093/cercor/bhy299Bae, J. W., Yi, J. H., Choe, S. Y.,
Li, Y., & Jung, M. W. (2025). Cortical VIP neurons as a critical node
for dopamine actions. *Science Advances*, *11*(1), eadn3221.
https://doi.org/10.1126/sciadv.adn3221Bakken, T. E., Jorstad, N. L., Hu,
Q., Lake, B. B., Tian, W., Kalmbach, B. E., Crow, M., Hodge, R. D.,
Krienen, F. M., Sorensen, S. A., Eggermont, J., Yao, Z., Aevermann, B.
D., Aldridge, A. I., Bartlett, A., Bertagnolli, D., Casper, T.,
Castanon, R. G., Crichton, K., ... Lein, E. S. (2021). Comparative
cellular analysis of motor cortex in human, marmoset and mouse.
*Nature*, *598*(7879), 111--119.
https://doi.org/10.1038/s41586-021-03465-8Darmanis, S., Sloan, S. A.,
Zhang, Y., Enge, M., Caneda, C., Shuer, L. M., Hayden Gephart, M. G.,
Barres, B. A., & Quake, S. R. (2015). A survey of human brain
transcriptome diversity at the single cell level. *Proceedings of the
National Academy of Sciences of the United States of America*,
*112*(23), 7285--7290. https://doi.org/10.1073/pnas.1507125112Glausier,
J. R., Khan, Z. U., & Muly, E. C. (2009). Dopamine D1 and D5 receptors
are localized to discrete populations of interneurons in primate
prefrontal cortex. *Cerebral Cortex (New York, N.Y.: 1991)*, *19*(8),
1820--1834. https://doi.org/10.1093/cercor/bhn212Gouwens, N. W.,
Sorensen, S. A., Baftizadeh, F., Budzillo, A., Lee, B. R., Jarsky, T.,
Alfiler, L., Baker, K., Barkan, E., Berry, K., Bertagnolli, D., Bickley,
K., Bomben, J., Braun, T., Brouner, K., Casper, T., Crichton, K.,
Daigle, T. L., Dalley, R., ... Zeng, H. (2020). Integrated
morphoelectric and transcriptomic classification of cortical GABAergic
cells. *Cell*, *183*(4), 935-953.e19.
https://doi.org/10.1016/j.cell.2020.09.057Kaplanis, J., Samocha, K. E.,
Wiel, L., Zhang, Z., Arvai, K. J., Eberhardt, R. Y., Gallone, G.,
Lelieveld, S. H., Martin, H. C., McRae, J. F., Short, P. J., Torene, R.
I., de Boer, E., Danecek, P., Gardner, E. J., Huang, N., Lord, J.,
Martincorena, I., Pfundt, R., ... Retterer, K. (2020). Evidence for 28
genetic disorders discovered by combining healthcare and research data.
*Nature*, *586*(7831), 757--762.
https://doi.org/10.1038/s41586-020-2832-5Nguyen, H. T., Bryois, J., Kim,
A., Dobbyn, A., Huckins, L. M., Munoz-Manchado, A. B., Ruderfer, D. M.,
Genovese, G., Fromer, M., Xu, X., Pinto, D., Linnarsson, S., Verhage,
M., Smit, A. B., Hjerling-Leffler, J., Buxbaum, J. D., Hultman, C.,
Sklar, P., Purcell, S. M., ... Stahl, E. A. (2017). Integrated Bayesian
analysis of rare exonic variants to identify risk genes for
schizophrenia and neurodevelopmental disorders. *Genome Medicine*,
*9*(1), 114. https://doi.org/10.1186/s13073-017-0497-yScala, F., Kobak,
D., Bernabucci, M., Bernaerts, Y., Cadwell, C. R., Castro, J. R.,
Hartmanis, L., Jiang, X., Laturnus, S., Miranda, E., Mulherkar, S., Tan,
Z. H., Yao, Z., Zeng, H., Sandberg, R., Berens, P., & Tolias, A. S.
(2021). Phenotypic variation of transcriptomic cell types in mouse motor
cortex. *Nature*, *598*(7879), 144--150.
https://doi.org/10.1038/s41586-020-2907-3Siletti, K., Hodge, R., Mossi
Albiach, A., Lee, K. W., Ding, S.-L., Hu, L., Lönnerberg, P., Bakken,
T., Casper, T., Clark, M., Dee, N., Gloe, J., Hirschstein, D.,
Shapovalova, N. V., Keene, C. D., Nyhus, J., Tung, H., Yanny, A. M.,
Arenas, E., ... Linnarsson, S. (2023). Transcriptomic diversity of cell
types across the adult human brain. *Science*, *382*(6667), eadd7046.
https://doi.org/10.1126/science.add7046Yao, S., Harder, A., Darki, F.,
Chang, Y.-W., Li, A., Nikouei, K., Volpe, G., Lundström, J. N., Zeng,
J., Wray, N. R., Lu, Y., Sullivan, P. F., & Hjerling-Leffler, J. (2025).
Connecting genomic results for psychiatric disorders to human brain cell
types and regions reveals convergence with functional connectivity.
*Nature Communications*, *16*(1), 395.
https://doi.org/10.1038/s41467-024-55611-1
