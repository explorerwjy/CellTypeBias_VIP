**Rare mutations implicate CGE interneurons as a vulnerable axis of
cognitive deficits across psychiatric disorders**

Stephanie A. Herrlinger^1,\ 2^\*, Jiayao Wang^3^\*, Bovey Y Rao^1,2,4^,
Jonathan Chang^3^, Joseph A. Gogos ^1,2,5,6,7^†, Attila Losonczy^1,2^†
and Dennis Vitkup^3^†

^1^ Department of Neuroscience, Columbia University, New York, NY,
United States

^2^ Mortimer B. Zuckerman Mind Brain Behavior Institute, Columbia
University, New York, NY, United States

^3^ Department of Systems Biology and Department Biomedical Informatics,
Columbia University, New York, NY, United States

^4^ Doctoral Program in Neurobiology and Behavior, Columbia University,
New York, NY, United States

^5^ Department of Physiology, Columbia University, New York, NY, United
States

^6^ Department of Psychiatry, Columbia University, New York, NY, United
States

^7^ Stavros Niarchos Foundation Center for Precision Psychiatry and
Mental Health, Columbia University, New York, NY, United States

\*These authors contributed equally to this work

† Correspondence should be addressed to JG
([jag90\@columbia.edu](http://jag90@columbia.edu)), AL
([al2856\@columbia.edu](http://al2856@columbia.edu)) and DV
(<dv2121@cumc.columbia.edu>)

**Abstract**

**Neuropsychiatric disorders such as autism spectrum disorder (ASD) and
schizophrenia (SCZ) share genetic risk factors, including genes affected
by rare high-penetrance single nucleotide variants (SNVs) and copy
number variants (CNVs). ASD and SCZ exhibit both overlapping and
distinct clinical phenotypes. Cognitive deficits and intellectual
disability---critical predictors of long-term outcomes---are common to
both conditions. To investigate shared and disorder-specific
neurobiological impact of highly penetrant rare mutations in ASD and
SCZ, we analyzed human single-nucleus whole-brain sequencing data to
identify strongly affected brain cell types. Our analysis revealed
caudal ganglionic eminence (CGE)-derived GABAergic interneurons as a key
nexus for cognitive deficits across these disorders. Notably, genes
within 22q11.2 deletions, known to confer a high risk for SCZ, ASD, and
cognitive impairment, showed a strong expression bias toward vasoactive
intestinal peptide-expressing cells (VIP+) among CGE subtypes. To
explore perturbations of VIP+ GABAergic interneurons in the 22q11.2
deletion syndrome *in vivo*, we examined their activity in the
*Df(16)A^+/-^* mouse model during a spatial navigation task and observed
reduced activity along with altered responses to random rewards. At the
population level, VIP+ interneurons exhibited impaired spatial encoding
and diminished subtype-specific activity suggesting deficient
disinhibition in CA1 microcircuits in the hippocampus, a region
essential for learning and memory. Overall, these results demonstrate
the crucial role of CGE-derived interneurons in mediating cognitive
processes that are disrupted across a range of psychiatric and
neurodevelopmental disorders.**

**\
**

**Introduction**

Psychiatric disorders exhibit distinct clinical symptoms, yet frequently
share genetic risk factors and underlying cellular and molecular
mechanisms^1--6^. For example, schizophrenia (SCZ) is characterized by
psychosis, including hallucinations and delusions, as well as cognitive
control impairments leading to disorganization, impulsivity, and
difficulties with working memory^7--11^. In contrast, autism spectrum
disorder (ASD) is characterized by restricted interests, difficulties
with social communication, and cognitive rigidity along with repetitive
behaviors^8,12^. Several copy number variations (CNVs), such as the
22q11.2 deletion, have been implicated in both SCZ and ASD, as well as
in other psychiatric disorders^13,14^. Among their diverse symptoms,
cognitive impairments and intellectual disability (ID) are frequently
comorbid with neuropsychiatric phenotypes^15^ and are especially
pronounced in cases with high genetic risk^15--18^. Cognitive deficits
in psychiatric disorders are particularly debilitating, and span a broad
range of functions, including memory, attention, reasoning,
decision-making and social cognition^19^. The severity of these deficits
is a strong predictor of patients' long-term health outcomes and overall
well-being^20--22^.

Despite advances in the genetics of disorders such as ASD and
SCZ^1,2,23^, therapeutic progress in treating their cognitive symptoms
remains limited, largely due to gaps in understanding the
neurobiological mechanisms through which genetic risk factors impact
cognition. To better understand the relationships between genetics and
psychiatric phenotypes, researchers are increasingly integrating genetic
analyses with single-cell and single-nucleus sequencing data to identify
cell types affected across psychiatric disorders^24--26^. Animal models,
particularly mouse models that replicate rare, high-risk genetic
insults, play a critical role in this integrative approach. The
evolutionary conservation of brain regions, cell types, and core
biological processes between mice and humans^27--30^ facilitates the
detailed examination of genetic risk factors at both the cellular and
circuit levels^31--36^.

In this study, we investigated the neurobiological impact of rare,
highly penetrant mutations in ASD, SCZ, and 22q11.2 deletion syndrome
across diverse brain cell types, using recently released human
single-nucleus whole-brain RNA sequencing data^37^ (Figure 1, left). Our
analyses revealed distinct disorder-specific mutation biases affecting
different cell types, with caudal ganglionic eminence (CGE)-derived
interneurons consistently implicated across genetic insults associated
with cognitive dysfunction, including ASD with intellectual disability,
SCZ, and Developmental Disorders and Intellectual Disability (DD/ID)
(Figure 1, middle). Furthermore, our analysis of genes within the
22q11.2 deletion revealed that VIP+ CGE-derived interneurons were
strongly affected in this condition. Using *in vivo* two-photon calcium
imaging of VIP+ GABAergic interneuron activity in a mouse model of
22q11.2 deletion, we observed specific alterations in
interneuron-selective interneuron (ISI) cell types, which primarily
inhibit other interneurons (Figure 1, right). These findings implicate
ISI dysfunction in circuit-level abnormalities contributing to cognitive
and behavioral phenotypes of this syndrome, and likely many other
psychiatric conditions.

![](media/image1.png){width="6.5in"
height="3.0569444444444445in"}**Figure 1. Overview of study.** **Left:**
Integration of genetic and single-cell transcriptomic data was used to
identify cell type-specific disorder biases. Genes implicated in ASD and
SCZ^1,2^ and genes in the 22q11.2 deletion were used to calculate cell
type mutation biases based on human brain snRNA-seq data from 461 cell
types^37^. **Middle:** Disorder cell type bias analysis revealed shared
cell type risk patterns between distinct disorders with overlapping
cognitive phenotypes. **Right:** Experimental analysis of the 22q11.2
deletion mouse model was used to identify cell type specific deficits
*in vivo* using 2-photon calcium imaging to capture neuronal activity
during a spatial navigation task.

**Results**

**Identifying disorder-specific cell type mutation biases by leveraging
genetic and single-nucleus expression data**

Individually, highly penetrant rare mutations typically have
substantially greater functional impacts on psychiatric phenotypes
compared to common variants^1,38^. The availability of
well-characterized rare variant datasets from ASD and SCZ cohorts has
allowed us to perform direct cross-disorder comparisons. For ASD, we
utilized a recent meta-analysis performed by the SPARK consortium^2^,
which comprised *de novo* mutations from 42,607 ASD probands. This
analysis identified 61 genes with exome-wide significance by de novo
mutations. For schizophrenia, we used a case-control study involving
24,248 schizophrenia patients and 97,322 control individuals from the
SCHEMA consortium^1^. In our analysis, we selected the top 61 SCZ genes
based on genetic association strength to match the number of considered
ASD-associated genes.

A significant fraction of probands in the ASD cohort exhibit low IQs (IQ
≤ 70, which we refer to as "ASD with ID")^39,40^ (Figure S1). To
facilitate a more appropriate comparison with SCZ, where the average
premorbid IQ is approximately 94^41,42^, we selected higher IQ (IQ \>70)
ASD cases from the SPARK cohort (hereafter referred to as "ASD", with an
average IQ approximately 95). Notably, individuals with SCZ who carry
rare mutations have IQ scores that are only marginally lower than those
in the broader schizophrenia population^43^. Furthermore, genes
implicated in higher-IQ ASD also share more similar developmental
expression profiles with SCZ genes compared to those associated with
lower-IQ ASD (Figure S2). It is important to emphasize that even though
SCZ and ASD cases have comparable overall IQs, the nature of cognitive
dysfunction in these conditions is often distinct. While SCZ is usually
characterized by broad deficits across multiple cognitive domains, such
as working memory and executive function^44,45^, high-functioning ASD
cases tend to show more selective, domain-specific cognitive
impairments, particularly in social cognition and communication^46,47^.

To identify the brain cell types most impacted by rare coding mutations
in ASD and SCZ, we used a recent single-nucleus sequencing dataset
collected from the entire human brain^37^. This dataset encompassed
approximately 105 anatomical locations dissected from three different
donors across the forebrain, midbrain, and hindbrain, with 3,369,219
individual nuclei classified into 31 cell type superclusters,
representing broad categories of major brain cell classes, and 461 cell
type clusters, which provide finer-grained subdivisions based on robust
transcriptomic similarity. (Figure 1). To evaluate the impact of
mutations on brain cell types in ASD and SCZ, we calculated mutation
biases across brain cell types by integrating gene expression with
disorder-specific mutation burden (Figure S3, see Methods). We first
computed a cell type specificity score for each gene, quantifying how
selectively a gene is expressed in a given cell type relative to
others^26^. To that end, we first quantified gene expression levels
using transcripts per million (TPM) for each cell type, and then used
TPM values to compute specificity scores capturing the enrichment of
each gene's expression in different cell types. Next, to assess how
strongly each cell type is implicated in the disorders, we integrated
the expression specificity scores with each gene's mutation burden,
using *de novo* mutation counts for ASD and effective case mutation
counts for SCZ. The resulting cell-type mutation bias captured how
selectively disease genes are expressed in each cell type, weighted by
their mutational burden in diagnosed patient cohorts. To evaluate
statistical significance of our mutation biases, we compared the
observed bias for each cell type cluster to a null distribution
generated from 10,000 random gene sets matched in size to the
disorder-specific gene set, with the original mutation weights randomly
reassigned to these sets. Empirical P‑values were calculated from these
simulations and corrected for multiple testing across 461 clusters.
These mutation bias results were robust across a range of gene set sizes
(Figure S4).

**Shared cell type mutation biases in ASD and SCZ**

The cell types that showed strong mutation biases in ASD included
amygdala excitatory neurons, medium spiny neurons, hippocampal dentate
gyrus neurons, and deep-layer intratelencephalic (IT) neurons (Figure
S5A). The same cell type superclusters also showed strong mutational
impacts in ASD with ID, with additional significant cell types including
excitatory neurons from the cortex and hippocampus, as well as cortical
interneurons (Figure S5B). For SCZ, cells with strong mutational impacts
included several cell types implicated in recent GWAS studies of
SCZ^25,26^. Specifically, MGE interneurons were the most strongly
impacted cell type, with eccentric medium spiny neurons, Lamp5-LHX6 and
Chandelier neurons, CGE interneurons, upper- and deep-layer IT neurons,
and midbrain-derived inhibitory neurons also strongly affected. Other
significantly impacted cell types included amygdala excitatory neurons,
mammillary body neurons, and hippocampal dentate gyrus and CA1-3 neurons
(Figure S5C).

Several cell type superclusters were strongly impacted by mutations in
both ASD and SCZ, including excitatory neurons in the amygdala and
hippocampus, as well as cortical excitatory neurons (such as upper-layer
IT, deep-layer IT, deep-layer corticothalamic, and layer 6b neurons)
(Figure 2C, Figure S5A-S5C, Figure S7A-S7D). These findings are
consistent with several previous studies demonstrating frequent
perturbations of these neurons across brain disorders. For example,
amygdala neurons and their cortical and hippocampal projections have
been implicated in a broad range of behavioral disturbances in SCZ and
ASD, including deficits in emotion, motivation, and social
interactions^48,49^. Moreover, the hippocampus is known to play a key
role in social and cognitive deficits in ASD^50^, as well as in learning
deficits and positive symptoms associated with SCZ^51,52^. Deep-layer IT
neurons, which are central to sensory processing, decision-making, and
learning, are also known to be strongly affected in both disorders^53^.
Overall, these shared neuronal disruptions may underlie the overlapping
dysfunction of brain circuits mediating psychiatric phenotypes^54^.

Our analysis revealed a strong correlation between ASD and SCZ mutation
biases across neuronal cell type clusters (Figure 2A), with cell types
showing strong mutation biases in ASD usually also showing strong
mutation biases in SCZ (Spearman's correlation = 0.68, P \< 10^-10^).
Notably, the correlation between cell type biases was not due to the
overlap between ASD and SCZ in implicated genes (Spearman correlation =
0.64, P \< 10^-10^, after removing three overlapping genes). To
understand what influences the shared cell type biases between the
disorders, we investigated the role of genes with strong constraint
against loss-of-function mutations, as constrained genes often have a
substantial mutation burden in psychiatric disorders^1,2,39,55,56^. When
we progressively removed genes from the considered disorder-specific
sets (one gene at a time, ranked by their loss-of-function constraint
LOEUF score^57^, see Methods), the correlation between ASD and SCZ
substantially decreased after the removal of approximately 15 genes
(Figure 2B, red); in contrast, removing the least constrained genes
preserved the substantial correlation between mutation biases (Figure
2B, blue). We next asked whether gene expression levels across cells
might also contribute to the shared correlation of mutation biases
between ASD and SCZ, as both disorders involve genes highly expressed in
the brain^3,58^. Interestingly, in contrast to gene mutation constraint,
removing highly expressed genes had little impact on the cell type bias
correlation (Figure S6). These findings suggest that highly constrained
genes harboring disease-associated mutations tend to be preferentially
expressed in cell types that are primarily targeted by mutations in both
disorders.

**Distinct cell type mutation biases in ASD and SCZ**

Given the overall similarity in cell type mutation biases between ASD
and SCZ, we next sought to identify cell types preferentially affected
in one disorder relative to the other. To that end, we calculated, for
each cell type, the differential mutation bias between ASD and SCZ
across cell type clusters within each cell type supercluster. This
differential bias measure quantifies how strongly each cell type
supercluster is preferentially affected in one disorder compared to the
other, potentially underlying the unique behavioral phenotypes
associated with each disorder. We observed that mutations in ASD and SCZ
are differentially biased toward several specific neuron types (Figure
2C). Specifically, medium spiny neurons (MSNs) showed substantially
stronger mutation biases in ASD compared to SCZ (Figure 2D). This agrees
with previous findings implicating MSNs in ASD in computational^58^ and
experimental studies^59,60^. In contrast, GABAergic interneurons,
including MGE-derived, CGE-derived, and Lamp5-LHX6/Chandelier
interneurons showed the strongest differential biases towards SCZ
compared to ASD (Figure 2C, 2E and 2F), which is consistent with
previous studies implicating interneurons in SCZ^25,26,61,62^. As
mentioned above, cognitive impairments are widespread in schizophrenia,
while cognitive symptoms in high-functioning ASD are often more
localized and specific^44--47^. Considering this, we investigated which
neuron types were more affected in ASD probands with ID compared to ASD
probands without ID. This analysis showed that all three interneuron
superclusters had substantially stronger biases toward ASD with ID
compared to high-functioning ASD, with CGE interneurons showing the
strongest differential bias (Figure S7E, Figure 2G). Overall, these
results suggest a potential role of interneurons, particularly
CGE-derived interneurons, in cognitive dysfunctions observed in both ASD
with ID and SCZ.

![](media/image3.png){width="6.03284886264217in"
height="5.774389763779528in"}

**Figure 2. ASD and SCZ Mutation bias towards different brain cell
types**. **A**. The correlation of mutation biases between ASD and SCZ
across all neuronal cell types (Spearman R = 0.68, P \< $10^{- 10}$).
Each dot represents a neuron cell type cluster. **B**. Spearman
correlation of cell type bias between ASD and SCZ genes, assessed by
sequential gene removal. Genes were removed one at a time from each
disorder's gene set in parallel, either starting from the most
constrained (blue) or least constrained (red) based on LOEUF score. At
each step, ASD--SCZ bias correlations were recomputed using mutation
biases obtained from the remaining genes. The gray line shows the
average correlation for random gene removals, with error bars indicating
±1 standard deviation. **C**. Bar plot comparing mutation biases across
cell type superclusters between ASD (blue) and SCZ (orange). All
neuronal cell types at supercluster level are displayed and sorted by
average differential bias between ASD and SCZ. Error bars represent
standard error across cluster biases from corresponding supercluster.
Asterisks indicate significance level between groups \*: P \< 0.05.
\*\*: P \< 0.01. \*\*\*: P \< 0.001, \*\*\*\*: P \< 0.0001, Mann-Whitney
U test with FDR correction. **D-G.** Detailed comparison of mutation
biases in specific neuron superclusters, each point representing a cell
type cluster from the corresponding supercluster, P-values were
calculated based on Mann--Whitney U test with FDR correction. **(D)**
Medium spiny neurons show stronger bias in ASD. **(E)** MGE interneurons
exhibit stronger bias in SCZ, **(F)** CGE interneurons display bias
toward SCZ and **(G)** CGE interneurons show stronger bias when
comparing ASD with ID to ASD.

**CGE interneurons are strongly affected in multiple psychiatric
phenotypes with cognitive impairments**

To further explore the involvement of CGE interneurons in cognitive
function, we analyzed several additional cohorts with related phenotypes
(Figure 3A, see Methods). Consistent with our observation of stronger
CGE interneuron bias in SCZ and ASD with ID, we also observed a strong
bias toward CGE interneurons (Figure 3A; Figure S7H) in the
developmental disorders and intellectual disability (DD/ID) cohort^56^,
which is characterized by severe cognitive impairments. Furthermore, we
analyzed rare protein-coding variants available in the UK Biobank
associated with cognitive functions, including verbal-numerical
reasoning (VNR), reaction time (RT), and educational attainment
(EDU)^63^. We found that genes negatively associated with VNR (VNR-),
i.e., genes with mutations leading to worse VNR performance, an
assessment of fluid intelligence and best proxy to general intelligence
among the traits^64^ (see Methods), exhibited a pronounced bias toward
CGE interneurons compared to genes positively associated with VNR (VNR+,
Figure 3A and S7G). When we compare bias contrast of all three major
interneuron groups, we find in contrast to CGE interneurons, MGE and
Lamp5-LHX6/Chandelier interneurons showed relatively weaker mutation
biases across the additional cohorts we examined (Figure 3B, Figure S8).
The consistently strong CGE interneuron mutation bias across different
cohorts suggests that these neurons play a key role in cognitive
processes.

![](media/image5.png){width="6.5in" height="2.748611111111111in"}

**Figure 3. CGE interneurons are impacted in multiple cohorts with
cognitive abnormalities.** **A.** CGE interneuron mutation biases across
different cohorts shows stronger bias in disorders with cognitive
deficits. Each dot represents a CGE interneuron subtype. Comparison
includes ASD (without ID), verbal-numerical reasoning negative and
positively associated genes (VNR- and VNR+), schizophrenia (SCZ), ASD
with ID, and developmental disorders/intellectual disability (DD/ID).
P-values indicate the significance of differences in CGE interneuron
mutation biases between groups (e.g., SCZ, DD/ID and ASD with ID vs.
ASD, VNR-- vs. VNR+) and were calculated using the Mann--Whitney U test
with FDR correction. **B.** Pairwise mutation bias differences between
disorders across 3 major interneuron superclusters. For each disorder
pair (Group1 -- Group2), the y-axis shows the difference in mutation
bias per cell type cluster. Each dot represents one cluster from the
three interneuron superclusters (CGE, MGE, LAMP5-LHX6/Chandelier).
Colored boxplots indicate the supercluster of origin. Statistical
comparisons were conducted between CGE clusters and each of the other
two superclusters (MGE, LAMP5-LHX6/Chandelier) using the Mann--Whitney U
test. Significant (0.05/8 comparisons) differences are annotated with
the corresponding p-values.

**VIP+ interneurons are strongly affected by genes in the 22q11.2
deletion**

Copy number variants (CNVs) consistently stand out as highly penetrant
genetic risk factors for SCZ and ASD^15,36,65,66^. 22q11.2 deletion
syndrome is a neurodevelopmental disorder characterized by a
significantly elevated risk of schizophrenia in adulthood, along with
comorbidities such as ASD and intellectual disability^14,67^. Cognitive
impairment is a key feature of 22q11.2 deletion syndrome, with nearly
half of affected individuals displaying an IQ below 70^68,69^. Similar
to most disease-associated CNVs, the full range of phenotypic symptoms
results from the synergistic effects of multiple genes within the
22q11.2 deletion^14^. Given the significant SCZ risk and cognitive
disability phenotypes associated with the 22q11.2 deletion, we next
investigated the bias toward CGE and other interneurons among the
deletion genes. The 22q11.2 deletions occur in two main forms: a 3 Mb
deletion (containing 46 genes) and a smaller 1.5 Mb deletion (containing
28 genes). Analysis of cell type biases for genes in the 3 Mb deletion
(Figure S9) revealed a substantial bias toward CGE-derived interneurons
compared to other major interneuron types (Figure 4A). Notably, among
CGE interneurons, VIP+ cell types exhibited significantly stronger
biases compared to VIP- CGE cell types (Figure 4B). This pattern
remained consistent and significant when examining genes within the
1.5 Mb deletion, which corresponds to the 16qA13 deletion in mouse
(Figure 4B). The selective impact on VIP+ interneurons compared to other
CGE interneuron types suggests a particularly important role in
mediating the severe cognitive and behavioral phenotypes associated with
22q11.2 deletion syndrome. CGE-derived VIP+ interneurons modulate neural
circuits by reducing inhibition and enhancing excitatory
activity^70,71^, an essential mechanism for cognitive functions such as
perception, attention, learning and memory^72^. They achieve this by
suppressing Sst+, Cck+, and possibly Pvalb+ GABAergic
interneurons^73,74^, thereby indirectly amplifying pyramidal neuron
activity. Moreover, VIP+ interneurons interact with cholinergic,
serotonergic, and dopaminergic systems^74--79^ to maintain the
excitation-inhibition balance critical for network stability and
effective cognitive processing^72^.

![](media/image7.png){width="6.5in" height="3.063888888888889in"}

**Figure 4. VIP+ CGE interneurons are preferentially affected by 22q11.2
deletion genes. (A)** Boxplot showing expression biases for three major
interneuron types (CGE interneurons, Lamp5/LHX6 and Chandelier cells,
and MGE interneurons) in the 22q11.2 deletion region. CGE interneurons
exhibit the strongest bias compared to other interneuron populations
(P-values shown for pairwise comparisons, Mann-Whitney U-test). **(B)**
The CGE interneuron cell type bias across different 22q11.2 gene sets,
including complete deletion genes (3MB deletion, 46 genes) and 28 genes
included in 1.5MB deletion which is equivalent to mouse 16qA13 deletion.
VIP+ CGE types shown in red and VIP- CGE types shown in blue. Difference
between VIP+ and VIP- bias is significant for all 22q11.2 genes and
1.5MB deletion genes (P=0.04 and 0.01, respectively, Mann-Whitney
U-test.)

**Altered VIP+ interneuron activity and population specificity in a
mouse model of 22q11.2 deletion syndrome**

Having identified a robust genetic association between mutation biases
toward CGE interneurons and cognitive phenotypes across multiple genetic
disorders, and a prominent bias of VIP+ cell types compared to other
interneurons in 22q11.2 deletion syndrome, we next experimentally
investigated whether and how this bias is reflected in the altered
function of these interneurons within hippocampal circuitry. To this
end, we used the *Df(16)A*^+/-^ mouse model, which faithfully
recapitulates the human 1.5 Mb 22q11.2 deletion^13,33,35,80^. This
model, known for deficits in working and episodic memory, synaptic
plasticity, and cortical and hippocampal development, has proven
invaluable for uncovering novel biological pathways relevant to
cognitive dysfunction and psychiatric disorders^14,35,36,67,80--83^. We
therefore examined the physiological impact of the 22q11.2 deletion on
specific CGE/VIP+ cell types, thus advancing our understanding of the
cellular mechanisms by which CGE-derived VIP+ interneurons mediate
cognitive processes are disrupted in 22q11.2 deletion syndrome.

To examine VIP+ interneuron activity *in vivo*, we crossed
*Df(16)A^+/-^* mice with VIP-Cre mice and locally injected an adenovirus
carrying a floxed GCaMP8m or GCaMP7f calcium indicator into the CA1
region to selectively express GCaMP in VIP+ interneurons (Figure S10A).
Using two-photon three-dimensional acousto-optic deflector
(3D-AOD)^84--86^ microscopy, we then recorded VIP+ interneuron activity
in a 3D volume within the dorsal hippocampal CA1 region (Figure
S10B-S10D). Tissue registration and *post-hoc*
immunohistochemistry^32,84,87^ (see Methods, Figure S10E, S10F) allowed
us to identify four major VIP+ interneuron subtypes:
interneuron-selective interneurons type 2 (ISI2), interneuron-selective
interneurons type 3 (ISI3), VIP+ cholecystokinin basket cells (CCKBCs),
and a group of VIP-expressing cells, referred to as \"VIP+ only\" (VIPo)
cells, that were immunonegative for all other tested markers^70,88,89^.
The distribution of these subtypes across hippocampal layers did not
significantly differ between wild type and *Df(16)A^+/-^* mice (Figures
S10G and S10H), and comparable numbers of each interneuron subtype were
recorded from both wild type and mutant genotypes (Figures S10I). Thus,
our analysis did not find evidence that the 22q11.2 deletion affects the
number or proportion of VIP+ interneuron subtypes in the hippocampal CA1
region.

As previously established, global interneuron activity is strongly
driven by ambulatory state^90--93^, and in treadmill-based experiments
we consistently observed that locomotion is the greatest predictor of
global interneuron activity^32,84,94^. To investigate locomotion-related
activity of hippocampal VIP+ interneurons, mice were next trained to run
on a head-fixed, cue-rich 2-meter treadmill for randomly distributed
water rewards while undergoing calcium imaging (Figure 5A)^32,84,87^.
Behavioral metrics, including laps run, licking, and region occupancy on
the belt during the 10-minute trials, were similar across genotypes
(Figures S10J-S10M). Expected locomotion-related activity patterns were
preserved across genotypes, with ISI3, ISI2, and VIPo subtypes
positively correlated with movement, and VIP+ CCKBCs anticorrelated
(Figure 5B)^84,95^. However, in *Df(16)A^+/-^* mice, ISI3 and ISI2 cell
activity exhibited a significantly reduced correlation with locomotion
(Figure 5C), and ISI3, ISI2, and VIPo interneurons showed visibly
reduced calcium transient amplitudes (Figure 5B). Due to the spike-like
activity of ISI3, ISI2, and VIPo cells we were able to perform a
transient waveform analysis (see Methods^96^). This analysis confirmed
reduced amplitudes in ISI3, ISI2, and VIPo cells (Figures 5D-J), while
VIPo cell transients uniquely occurred at a higher frequency (Figure
5J). Reward responses, accessed by Peri-event Time Averages (PETA, see
Methods), exhibited genotype-dependent differences, including reduced
VIPo response amplitudes and enhanced ISI2 and ISI3 response amplitudes
(Figures 5K, 5L), suggesting that ISIs in mutant mice exhibit heightened
sensitivity to random salient features^94^. Overall, these results show
that VIP+ interneurons in the hippocampal CA1 region exhibit reduced *in
vivo* activity dynamics in the 22q11.2 deletion model and that reward
responses are significantly altered among specific VIP+ subtypes,
demonstrating physiological consequences of the genetic perturbation
predicted by our computational analysis.

![A collage of graphs and diagrams AI-generated content may be
incorrect.](media/image9.png){width="5.047619203849519in"
height="6.2636854768153984in"}

**Figure 5: Baseline changes in VIP+ cell type calcium dynamics in
mutant mice. A.** Imaging and behavior schematic: mice are trained to
run on a cue-rich belt to expect random water rewards to be deposited
across the belt (bottom) while head-fixed and imaged with 2-photon
3D-AOD microscopy (top). **B.** Example traces from all four VIP+ cell
types detected in a single session from a Wildtype littermate (left) and
a *Df(16)A^+/-^* (right) mouse. **C.** Pearson correlation coefficients
of mouse velocity against average subtype calcium activity. **D.**
Waveform comparison of average detected calcium transient from wildtype
littermates (black) and *Df(16)A^+/-^* (brown) ISI3 cells. **E.**
Average amplitude of detected ISI3 waveforms compared between genotypes.
**F.** Waveform comparison of average detected calcium transients from
wildtype littermates (black) and Df(16)A^+/-^ (brown) ISI2 cells. **G.**
Average amplitude of detected ISI2 waveforms compared between genotypes.
**H.** Waveform comparison of average detected calcium transients from
wildtype littermates (black) and *Df(16)A^+/-^* (brown) VIPo cells.
**I.** Average amplitude of detected VIPo waveforms compared between
genotypes. **J.** Number of transients from VIPo cells per cell between
genotypes. **K.** Averaged z-scored reward-onset peri-event time
averages (PETAs) between -1.0 and 1.0 z-scores by cell type across all
cells within cell type during the Random Foraging task for 3 seconds
pre- and 3 seconds post-event. **L.** Post - Pre PETA responses per cell
type before and after reward onset comparing the 3 seconds pre- to 3
seconds post-event. *CCKBC:* VIP+ and cholecystokinin + basket cell;
*ISI3:* Interneuron-selective interneuron type 3; *ISI2:*
Interneuron-selective interneuron type 2; *VIPo*: VIP+/Cal-/M2R-/CCK-
cells in layer oriens or layer pyramidale. \*: P \< 0.05. \*\*: P \<
0.01. \*\*\*\*: P \< 0.0001. See Table S1 for detailed statistics.

Hippocampal place cells, a subset of pyramidal neurons, are active in
specific environmental locations and provide cellular substrates for
episodic memory and goal-directed spatial learning^97,98^. In
*Df(16)A^+/−^* mice, our previous study revealed significant impairments
in learning performance when both contextual cues and reward locations
were altered during a spatial navigation task^82^. These learning
deficits were associated with reduced long-term stability of place
fields, impaired context-related remapping and plasticity, and a lack of
reward-related reorganization. Further investigation showed that
parvalbumin- and somatostatin-expressing GABAergic interneurons, which
are directly inhibited by VIP+ cell activity, displayed diminished place
fields and heightened reward-related activity, likely contributing to
the observed place cell deficits^32^. Changes in the activity dynamics
of VIP+ cells in dorsal CA1 may also reflect altered spatial information
encoding by these cells. To explore this, we next evaluated how well
mouse position could be decoded based on VIP+ subtype activity during
random foraging (Figure 5A) by training a support vector machine using
all recorded interneuron activity (Figure 6A). Decoding error did not
differ between genotypes when the total number of imaged cells per
session was not normalized. However, when controlling for the number of
cells tested (Figure 6B) or the number of laps run per session (Figure
6C), decoding error was significantly higher in *Df(16)A^+/-^* mice,
indicating that individual mutant cells carry less spatial information
than cells in wild type mice. Furthermore, analysis of place preference
in VIP+ cell activity revealed that fewer VIP+ cells in *Df(16)A^+/-^*
mice exhibited significant place fields, with a reduced proportion
displaying at least one place field (Figure 6D). These findings
demonstrate that VIP+ interneurons in the 22q11.2 deletion model exhibit
markedly impaired spatial encoding within the CA1 microcircuit.

Local brain circuits balance coordinated and unique activity depending
on input properties, connectivity, and brain state. Taking advantage of
our simultaneous recordings across multiple VIP+ subtypes in our
experiments, we computed pairwise Pearson's correlation coefficients
between calcium activity traces within and between different VIP+
interneuron subtypes (see Methods). We found that VIP+ cells in
*Df(16)A^+/-^* mice exhibited unusually high correlations between
subtypes during random foraging. This suggests a loss of
subtype-specific heterogeneity of *in vivo* activity dynamics,
particularly among disinhibitory VIP+ subtypes, indicating a potential
microcircuit-level disruption and impaired fine-scale coordination
between VIP+ interneuron subtypes (Figures 6E and 6F). These findings
reveal that VIP+ interneurons in the CA1 region of *Df(16)A^+/-^* mice
exhibit impaired spatial encoding --- reflected in increased decoding
errors and reduced place preference --- and that this disruption is also
accompanied by a loss of subtype-specific heterogeneity, evident from
unusually high inter-subtype correlations during the foraging task.
These physiological consequences align with our computational analyses
and illuminate how perturbations in specific subtypes of VIP+
interneurons can disrupt network dynamics and information processing
important for cognition.

![A close-up of a graph AI-generated content may be
incorrect.](media/image10.png){width="6.375in"
height="3.104407261592301in"}

**Figure 6. Reduced specificity in *Df(16)A^+/-^* VIP+ interneuron
population activity. A.** All recorded VIP+ interneuron activity was
used to train a support vector machine (*SVM*) to decode discrete
position bins (*D.E.*: decoding error) across the belt. Each bin is
equivalent to 2 cm. **B.** SVM results bootstrapped 100 times by number
of cells used per genotype per bin. **C.** SVM results bootstrapped 100
times by the number of laps used per genotype per bin. **D.** Proportion
of positive place fields detected using the Peak Place method. **E.**
Pairwise Pearson's correlation coefficient value matrix between cell
types of (∆*F*/*F*) traces across the trial during the Random Foraging
paradigm from WT (left) and *Df(16)A^+/-^* (right). **F.** Comparative
differences between cell type correlated activity between genotypes.
*CCKBC:* VIP+ and cholecystokinin + basket cell; *ISI3:*
Interneuron-selective interneuron type 3; *ISI2:* Interneuron-selective
interneuron type 2; *VIPo*: VIP+/Cal-/M2R-/CCK- cells in layer oriens or
layer pyramidale. \*: P \< 0.05. \*\*: P \< 0.01. \*\*\*: P \< 0.001.
See Table S1 for detailed statistics.

**\
**

**Discussion**

Cognitive dysfunction is central to psychopathology, as many psychiatric
conditions share overlapping yet challenging-to-treat cognitive
impairments^19--22,40,68^. For instance, schizophrenia often presents
with deficits in executive function and working memory^8,12^. Similarly,
ASD cases driven by rare damaging mutations^58,99^ are frequently
associated with pronounced cognitive deficits. Cognitive disabilities
are also common in genetic conditions such as 22q11.2 deletion syndrome,
where affected individuals frequently display difficulties with
attention, working memory, and executive processes^14,67,68^. In this
study, we focused on analyzing neuronal impact of rare highly penetrant
genetic mutations. The availability of high-quality genetic datasets,
each containing dozens of genes with significant disorder
associations^1,2,56^, now enables direct quantitative comparisons of
cell type mutation biases between the disorders. We observed a strong
correlation between cell type biases in ASD and SCZ, which was primarily
driven by mutations in highly constrained genes. We also identified
several cell types that were preferentially affected in each disorder,
possibly contributing to disorder-specific symptoms. For example, medium
spiny neurons of the striatum showed significantly stronger mutation
biases in ASD compared to SCZ, consistent with ASD\'s core symptoms of
repetitive behaviors and social impairments, which have frequently been
associated with striatal dysfunction^59,60^. In contrast, multiple
interneuron subtypes exhibited stronger mutation biases in SCZ,
consistent with recent GWAS-based studies^25,26^, particularly
parvalbumin-expressing interneurons crucial for gamma oscillations and
cognitive processing^100,101^, as well as somatostatin-expressing
interneurons important for information processing^61^. Supporting these
findings, transcriptomic studies have observed alterations in these
interneuron populations in SCZ patients^62^.

Our genetic analysis revealed dysfunction in CGE-derived GABAergic
interneurons as a central vulnerability contributing to cognitive
deficits across multiple psychiatric conditions (Figure 3).
Specifically, we found that CGE interneurons are more affected in SCZ
compared to high-functioning ASD, consistent with cognitive impairments
being a core and defining symptom of SCZ^21,22^. We also observed a
strong genetic impact on CGE interneurons in low-functioning compared to
high-functioning ASD, as well as in intellectual disability and 22q11.2
deletion syndrome. CGE interneurons also showed strong biases in UK
Biobank gene sets associated with lower cognitive performance (e.g.,
verbal-numerical reasoning)^63^.

Within the broader class of CGE interneurons, VIP+ subtypes were
specifically targeted by genes in the 22q11.2 deletion region, a highly
penetrant CNV conferring risk for schizophrenia, ASD, and intellectual
disability. To probe the functional consequences, we evaluated VIP+ CGE
derived interneuron activity *in vivo* using the *Df(16)A^+/−^* mouse
model. We observed significant defects at both the cellular and circuit
levels, likely contributing to broader hippocampal dysfunction
previously reported in this mouse model^32,82^. Specifically, VIP+
interneurons exhibited dampened calcium activity, suggesting a reduction
in the top-down disinhibition^72,102^ normally exerted by these cells in
cortical circuits. In addition, we found a reduction in the spatial
information carried by their activity. Because synaptic plasticity is
critical for fine-tuning synaptic connections and forming dynamic
spatial representations and episodic memories, this reduction likely
reflects impaired plasticity mechanisms. Such impairment may hinder
MGE-derived interneurons and pyramidal neurons from effectively encoding
and updating spatial information^32,82^, ultimately disrupting the
integration of spatial signals within the hippocampal network.
Furthermore, VIP+ interneuron subtypes exhibited increased correlation
in their *in vivo* activity, particularly among interneuron-selective
interneuron subtypes. This suggests a loss of their normally distinct
activity patterns, leading to reduced informational diversity in the
circuit. Impaired disinhibition within the CA1 microcircuit in patients
with 22q11.2 deletion syndrome likely disrupts the processing and
integration of salient inputs, resulting in cognitive deficits
associated with hippocampal dysfunction. Given that VIP+ disinhibitory
motifs are highly conserved across cortical circuits^70,71,103^, the
deficits we uncovered in the CA1 circuit are likely pervasive across
neocortical brain regions.

Because VIP+ interneurons inhibit other major interneuron subtypes,
their activation increases pyramidal neuron excitability, thereby
improving the detection of low-contrast signals^89,104^. Thus, these
interneurons are essential for learning and attention, as they respond
to salient stimuli by adjusting overall neural excitability. For
instance, during fear conditioning, blocking VIP+ interneurons disrupts
projection neuron activity and impairs fear learning^71^. In the
hippocampal CA1 region, proper disinhibitory VIP+ interneuron function
is vital for spatial learning as excessive inhibition of pyramidal
neurons hinders place cell remapping and performance, while their
activation enhances learning outcomes^94^. Similarly, in the prefrontal
cortex, activation of VIP+ interneurons improves working memory and
attention, whereas activation of Pvalb+, Sst+, or pyramidal neurons
often leads to detrimental effects^79,94,104,105^. Of note, impaired
disinhibition and CGE-derived VIP+ interneuron dysfunction have also
been previously described in other neurodevelopmental
disorders^106--108^. In Rett syndrome, for example, selective deletion
of Mecp2 in CGE/VIP+ interneurons, but not in Pvalb+ or Sst+
interneurons, leads to deficits in cortical state transitions and social
memory^107^. In Dravet syndrome, SCN1A mutations lead to persistent VIP+
interneuron dysfunction in juvenile and adult mutant mice, contributing
to social and cognitive deficits^106,108^.

While further research is needed to link specific disease-related
behaviors to CGE/VIP+ interneuron circuit dysfunction, our combined
computational and experimental results suggest that disruption of VIP+
interneuron activity may be a key driver of cognitive deficits in
22q11.2 deletion syndrome and potentially in other
cognitively-debilitating neurodevelopmental and psychiatric
disorders^109^. These findings position CGE/VIP+ interneurons as a
critical point of dysfunction. Given the persistent cognitive
comorbidities across multiple psychiatric disorders, pharmacologically
targeting CGE interneuron dysfunction is a potential avenue for
therapeutic intervention.

**METHODS**

**Rare mutation genetic data for neuropsychiatric disorders and
phenotypes**

For the analysis on cell type mutation biases of ASD and SCZ, we
primarily used the aggregated mutation data obtained from the SFARI
Simons Foundation Powering Autism Research for Knowledge (SPARK) cohort
and the Schizophrenia Exome Sequencing Meta-Analysis (SCHEMA) cohort. In
the SPARK study^2^, de novo mutations were identified through familial
sequencing of 42,607 ASD probands and 8,267 unaffected siblings. We
focused our analysis on de novo likely gene-disrupting (LGD) mutations
and damaging missense (Dmis) mutations. Dmis mutations were identified
based on the REVEL score^110^, which estimates the probability of a
missense mutation is pathogenic; mutations with REVEL scores higher than
0.5 were defined as deleterious. For our main analysis, we used 61 genes
that reached exome-wide significance (P-value \< 1.3x10^-6^) identified
in the SPARK study, using the DeNovoWEST statistical model^56^. In the
SCHEMA study^1^, which analyzed 24,248 schizophrenia cases and 97,322
controls, to match the number of analyzed ASD genes, we selected the top
61 genes with the strongest association for our analysis (p-value \<
1.0x10^-3^ or FDR \< 0.3). Dmis mutations in SCZ were identified based
on the MPC score^111^, and we included in our analysis LGD and Dmis
mutations with MPC scores greater than 3. We also investigated cell type
biases across different cohorts with cognitive abnormalities.
Specifically, we investigated genes within the 22q11.2 deletion region
as defined by
[Karayiorgou](https://www.nature.com/articles/nrn2841#auth-Maria-Karayiorgou-Aff1-Aff2)
*et al.*^14^. We also analyzed genes associated with verbal-numerical
reasoning (VNR) from the UK Biobank study by Chen *et al.*^63^. We
separately analyzed genes with positive and negative impacts (indicated
by the beta coefficients of gene-phenotype associations), focusing on
European population data for optimal statistical power. Additionally, we
analyzed genes associated with developmental disorders and intellectual
disability (DD/ID) from Kaplanis *et al.*^56^, ranking the top genes
according to DeNovoWEST p-values.

To validate the robustness of gene set selection, we tested our bias
calculations with varying numbers of genes. Our mutation biases were not
sensitive to the selection of 61 genes, when different numbers of genes
were tested, ranged from top 10 to top 200, ranked by association
significance to the disorder, bias correlations with our main selection
set remained very strong across a broad range of gene set sizes (Figure
S4). We chose 61 genes to balance statistical power from a larger gene
set with the individual gene\'s confidence for association with each
disorder.

**Human whole brain single-cell data**

For our single-cell analysis, we used the recent human brain
single-nucleus sequencing dataset ^37^, comprising 3,369,219 nuclei
sampled from 105 regions across the entire brain of three donors. This
comprehensive single-cell brain atlas was classified, based on cellular
transcriptional similarities and developmental origins, into 461 cell
type clusters organized into 31 cell type superclusters.

**Gene expression specificity score for individual cell types**

To identify neuronal populations preferentially impacted by rare
mutations, we calculated a gene expression specificity score for each
cell type cluster (Figure S3). To that end, for each gene $g$ in each
cell type $\text{ct}$, we first quantified gene expression using
Transcripts Per Million (TPM), calculated from raw UMI counts as:

$\text{TPM}_{g,ct} = \frac{\text{UMI}_{g,ct}}{\sum_{g\mathbf{'}\ }^{}\text{UMI}_{g\mathbf{'},ct}} \times 10^{6}$
(1)

where $\text{UMI}_{g,ct}$ is the expression count of gene $g$ in a cell
type $\text{ct}$, and the sum in the denominator is across all genes
$g\mathbf{'}$ in that cell type. In this way, TPM values normalize
expression across cell types by accounting for the differences in the
total transcript abundance. To filter out genes with very low expression
in individual cell types and reduce noise in the specificity
calculation, we applied a minimum TPM threshold of 0.1 for each gene in
each cell types. Additionally, genes with total UMI counts below 10,000
across all 3 million cells of the human brain cell atlas were excluded
from the analysis, removing genes with very low global expression that
leads to bad specificity estimate.

Next, we computed a cell type specificity score to quantify how
selectively a gene is expressed in each cell type. For each gene $g$,
specificity toward cell type $\text{ct}$ was calculated as:

$S_{g,ct} = n \cdot \frac{\text{TPM}g,ct}{\sum_{\text{ct}'}^{}\text{TPM}_{g,ct'}}$
(2)

where $n$ is the total number of cell types. In case of the uniform
expression of a gene across all cell types, the expected value of
$S_{g,ct}$ is 1*,* and higher values indicate stronger cell type
specific expression. To prevent extreme values from a small number of
genes from dominating downstream bias computations, we capped $S_{g,ct}$
at a maximum of 2, which represents the double of the expected value.

To further ensure comparability across cell types with different
expression distributions, we applied a per-cell type centering step. For
each cell type, $\text{ct}$, we subtracted the average specificity
across all human genes:

$S'_{g,ct} = S_{g,ct} - \mu_{\text{ct}},\ \text{where}\ \mu_{\text{ct}} = \frac{1}{m}\sum_{g = 1}^{m}S_{g,ct}$
(3)

with $m$ denoting the total number of genes. The resulting centered
specificity score $S'_{g,ct}$ reflects how strongly a particular gene is
enriched in a given cell type relative to all genes in that cell type.

**Disorder mutation bias for individual cell types**

To calculate cell type-specific mutation biases for each disorder
cohort, we weighted gene expression specificities ($S'_{g,ct}$, Equation
3) with their corresponding genetic burdens. This weighted approach
accounts for the relative genetic contributions (proportional to the
genetic burden) of individual genes to the disorder. Specifically, for
the ASD cohort, we quantified each gene\'s contribution with the
following mutation weight:

$w_{g}$ $= n_{g,\ \ LGD} + \ n_{g,Dmis}$ (4)

where $n_{g,LGD}$ and $n_{g,Dmis}$ are the numbers of observed *de novo*
LGD and Dmis mutations in gene $g$.

For the SCZ cohort, we calculated mutation biases by integrating both
*de novo* and inherited mutations identified in case-control studies.
Since the genetic evidence comes from these two distinct sources, we
computed a mutation counts that combine their contributions. For
case-control mutations, the effective count for each gene represents the
number of excess mutations in cases after subtracting the expected
counts, which were calculated by scaling control mutations by the ratio
of the number of cases to controls. We applied this approach separately
for each mutation class (LGD and Dmis) and then combined these values
with observed *de novo* mutations:

$n_{g,\ mutType}^{\text{eff}}$ ~=~ $n_{g,\ mutType}^{\text{case}}$ ~-~
$\frac{n_{g,\ mutType}^{\text{ctrl}}}{N_{\text{ctrl}}} \times N_{\text{case}}\ \  + n_{g,\ mutType}^{\text{de\ novo}}$
(5)

$w_{g}$ $= n_{g,\ LGD}^{\text{eff}} + \ n_{g,\ Mis3}^{\text{eff}}\ $

where $n_{g,\ mutType}^{\text{eff}}$ is the effective count of mutation
in gene $g$, $n_{g,\ mutType}^{\text{case}}$ is number of mutations
observed in case cohort in gene $g$, $n_{g,\ mutType}^{\text{ctrl}}$is
number of mutations observed in controls in gene $g$.
$n_{g,\ mutType}^{\text{de\ novo}}\ $ is number of de novo mutations in
gene $g$. $\text{mutType}$ includes LGD and Mis3 mutations.
$N_{\text{case}}$ and $N_{\text{ctrl}}$ are total number of cases and
controls in SCZ cohort, respectively.

For genes associated with the UK Biobank cognitive phenotypes, we
utilized each gene\'s effect size, i.e., beta coefficient from
gene-phenotype association tests, as the gene weight. This approach
assumes that genes with larger phenotypic effects contribute more
strongly to cell type biases relevant to cognition. For genes in the
22q11.2 deletion, we applied uniform weights for all genes ($w_{g}$ =
1), as each individual gene effectively contributes one LGD mutation to
the CNV.

Finally, based on the calculated mutation weights for each gene
(Equations 3-5), the mutation bias for all disorder-associated genes
toward a cell type $\text{ct}$ was calculated as follows:

$\text{Mutation}{\_ Bias}_{\text{ct}} = \ \frac{\sum_{g \in G}^{\ }{w_{\text{g\ }} \cdot \ S'_{g,ct}}}{\sum_{g \in G}^{\ }w_{\text{g\ }}}$
(6)

where the summation is over genes $g$ in the set $G$ containing all
implicated genes for a disorder or a phenotype, and $w_{g}$ is the
mutation weight of each gene $g$.

**Statistical significance of the cell type mutation biases**

To determine the statistical significance of the calculated cell type
mutation biases, we employed a permutation-based framework using matched
random gene sets. For each disorder or phenotype, we generated 10,000
random gene sets (selected from genes with valid expression specificity
score) matched in size to the target gene set (e.g., ASD, SCZ, DD/ID)
and transferred the mutation weights from the original genes to the
randomly selected ones. Mutation biases were then computed for each
random gene set using the same method as for the original data. At the
cluster level (461 cell types), we compared the observed mutation bias
of each cluster to the distribution of biases from the 10,000 matched
random gene sets. Empirical P-values were calculated per cluster and
corrected for multiple comparisons across all 461 clusters using the
Benjamini-Hochberg false discovery rate (FDR) procedure.

**Sensitivity of the correlation between ASD and SCZ cell type mutation
biases to gene properties**

To assess how gene-level properties---specifically mutation constraint
and expression level---influence the correlation of cell type mutation
biases between ASD and SCZ, we evaluated changes in this correlation as
a function of these properties. To that end, for each disorder, we
independently ranked its associated genes by either mutation constraint
or expression level. Gene constraint was quantified using LOEUF scores
from the gnomAD v4.0 release ^57^; if a gene lacked a LOEUF score in
v4.0, we used the corresponding value from gnomAD v2.1. Gene expression
levels were obtained from the BrainSpan developmental transcriptome
dataset^112^, calculated as the average log10 RPKM across all samples.
We then iteratively removed genes, one at a time, from each ranked gene
list for ASD and SCZ, removing the gene with the highest or the lowest
value in each set at each step. This ensured that the number of removed
genes remained equal between the two considered disorders. After each
step, we recalculated the cell type mutation biases for both ASD and SCZ
using the remaining genes and computed the Spearman correlation between
their mutation biases across all neuronal cell type clusters. As a
control, we also performed 1,000 iterations of random gene removal
(matched by number) and computed the corresponding distribution of
correlations. This approach allowed us to determine whether the observed
similarity between ASD and SCZ cell type biases was primarily driven by
genes with high constraint or high expression levels.

> **Difference in cell type mutation biases between ASD and SCZ**

To identify disorder-specific mutational biases at the supercluster
level, we assessed the difference in the mutation bias between disorders
(ASD vs. SCZ), quantified by the average difference in mutation biases
across cell type clusters forming each supercluster. This differential
bias measure describing how strongly each particular supercluster is
preferentially affected in one disorder compared to another, was
calculated as follows:

$\Delta\text{Bias}_{s}^{\text{Phen}o_{\text{AB}}} = \frac{\sum_{ct\_ 1}^{ct\_ n}{(\text{Bias}_{\text{ct}}^{\text{Phen}o_{A}} - \text{Bias}_{\text{ct}}^{\text{Phen}o_{B}})}}{n}$
(7)

where $\Delta\text{Bias}_{s}^{\text{Phen}o_{\text{AB}}}$ is differential
bias of supercluster $s$ between phenotype A and B.
$\text{Bias}_{\text{ct}}^{\text{Phen}o_{A}}$ is mutation bias for a cell
type cluster in supercluster $s$ for phenotype A, and
$\text{Bias}_{\text{ct}}^{\text{Phen}o_{B}}$ is mutation bias for a cell
type cluster in supercluster $s$ for phenotype B. The summation is over
all *n* cell type clusters within the supercluster.

To evaluate whether these differential biases were statistically
significant, we applied the Mann-Whitney U test to compare mutation bias
distributions between disorders for each supercluster (e.g., CGE
interneurons, MGE interneurons, medium spiny neurons). Multiple
hypothesis correction was performed using the Benjamini-Hochberg
procedure to control the false discovery rate across 31 superclusters.

**Viruses**

Cre-dependent recombinant adeno-associated virus (rAAV) for GCaMP7f
(rAAV1-syn-FLEX-jGCaMP7f-WPRE, Addgene \#1944820-AAV1, titer: ≥1:10\^13
vg/mL) and GCaMP8m ((rAAV1-syn-FLEX-jGCaMP8m-WPRE, Addgene \#
162378-AAV1, titer: ≥1:10\^13 vg/mL) were used to express GCaMP7f or
GCaMP8m in the right, dorsal hippocampus of VIP-Cre mice crossed with
Df(16)A^+/-^ mice. The Cre-dependent recombinase enabled expression of
calcium sensor GCaMP7f or GCaMP8m in VIP+ interneurons (INs).

**Viral injections and hippocampal window/headpost implant**

During isoflurane anesthesia, Df(16)A^+/-^ x VIP-cre mice were placed
onto a stereotaxic surgery instrument. Meloxicam and bupivacaine were
administered subcutaneously to reduce discomfort and provide analgesia.
The skin was cut along the midline to expose the skull, and a craniotomy
was performed over the right hemisphere with a dental drill. A Nanoject
injector was loaded with a sterile glass capillary filled with mineral
oil and rAAV1-syn-FLEX-jGCaMP7f-WPRE or rAAV1-syn-FLEX-jGCaMP8m-WPRE and
lowered into dorsal CA1 of the right hemisphere. Injections were made at
AP -2.2, ML -1.75, DV - 1.8, -1.6, -1.4, -1.2, -1.1, -1.0, -0.8 relative
to bregma with 64 nl of virus injected at each DV level. The pipette was
held within the brain for 5 minutes after the last injection and slowly
raised from the brain. The skin was sutured, and a topical antibiotic
was applied to the suture site. 1.0 mL of saline was subcutaneously
given to each mouse, while the mice recovered in their home cage with a
heat pad. Four days after viral injection, a circular cannula with a
glass window was inserted into the cavity produced by craniotomy as
previously described^113^. Similarly as described above, mice were
anesthetized with isoflurane and placed on a stereotaxic surgery
instrument. Meloxicam and bupivacaine were subcutaneously given to the
mice. The skull was exposed, and a 3 mm craniotomy was performed over
the right hippocampus at AP -2.2, ML -1.75 relative to bregma, centered
around the injection site. After removal of the dura, the cortex was
slowly aspirated with negative pressure and the consistent application
of ice-cold cortex buffer. When the horizontal fibers over CA1 were
exposed, aspiration was stopped and a sterile collagen membrane was
applied to stop bleeding. A 3 mm circular cannula with a glass window
was placed into the craniotomy and glued into the site with Vetbond. A
headpost was then attached to the posterior of the skull with layers of
dental cement. 1.0 mL of saline was subcutaneously given to each mouse,
while the mice recovered in their home cage with a heat pad. All mice
were monitored for 3 days of post-operative care until behavior training
began.

**Behavioral training and paradigms**

Following the post-op recovery period, mice were acclimated to handling
and head-fixation for a few days, while being water restricted to 85-90%
of their original weight. Mice were trained to run on a fabric, cue-free
burlap belt with random rewards that were automatically dispensed at a
decreasing number of random locations. Mice were trained for \~1 week on
the cue-free belt and then were transitioned to a cue rich belt during
the imaging experiments.

**Random Foraging Task:**

Water restricted mice were imaged on a cue rich familiar belt for 2 days
with 3 random rewards distributed across the length of the belt,
randomized every lap.

**AOD-based two-photon calcium imaging**

After the mice were acclimated and trained in the respective tasks, a
high-resolution structural Z-stack was captured 24 hours prior to
imaging, while the mouse was anesthetized with isoflurane. The reference
Z-stack was used to identify the 3D position of the GCaMP-positive
neurons for access by the custom-modified AOD (acousto-optic deflector)
microscope (Femto3D-ATLAS, Femtomics Ltd). The AOD microscope was used
as previously described in Geiller et al. 2020, and Vancura and Geiller
et al. 2023.

To provide stable transmission parameters during chronic imaging in the
entire 3D scanning volume, the AOD microscope was extended with a
high-speed and precision beam stabilization unit, which was directly
attached to the AOD scan head, sensitive to input beam misalignment. The
beam stabilization unit consisted of two quadrant detectors (PDQ80A and
TPA101, Thorlabs) and two broadband dielectric mirrors (Thorlabs)
mounted on motorized mirror mounts (Femtonics). Beam alignment was
performed by the LaserControl software (Femtonics).

A water immersion objective (x16 Nikon CF175) was used to focus at the
dorsal CA1 pyramidal layer and then fixed in position. A tunable
two-photon laser (Coherent Ultra II) was adjusted to wavelength = 920
nm. The high-resolution structural Z-stack was captured in CA1 from
\~100 um below to glass to \~500 um to span the entirety of CA1 from the
*stratum oriens* to the *stratum lacunosum-moleculare* with a 800x800
pixel capture with a resolution of 1.25 um per pixel with a 4 um step
size. Laser power and photomultiplier tube (PMT) detectors (GaAsP,
H10770PA-40, Hamamatsu) were compensated in power and gain throughout
the Z-stack capture. After the Z-stack was completed, the mouse was
returned to its home cage and given 24 hours to recover before
functional imaging.

Using the integrated software from Femtonics (MES), the Z-stack was
examined, and 100-150 interneurons were manually selected to generate a
spreadsheet of the x-y-z positions at the center of each respective
cell. These x-y-z positions served as the regions of interest (ROIs)
that were then used on subsequent days by the AOD microscope to find the
position of the corresponding interneurons. For each day of functional
imaging in each mouse, the identical field of view was identified within
the CA1 pyramidal layer and then the ROIs were loaded in using the
software for AOD 3D imaging. Each imaging session was taken for 10
minutes, with approximately 5-10 Hz for each recording experiment (frame
rate depended on the number of ROIs, ROI size, imaging resolution, and
scanning speed).

**Immunohistochemistry**

The immunohistochemical (IHC) staining strategy is consistent with our
previously published methods^84,87^. In brief: mice were collected the
last day after imaging and perfused with 20 mL Phosphate Buffered Saline
(PBS) followed by 20 mL of 4% paraformaldehyde (PFA, Electron Microscopy
Sciences) and the brain extracted for overnight fixation in 4% PFA. The
next day, the brains were washed 3 x 5 minutes with 1x PBS and
subsequently vibratome sliced for horizontal sections from the dorsal
side of the brain into 75 µM sections. The relevant sections were then
taken through our IHC protocol. Sections were permeabilized for 2 x 20
minutes (0.03% Triton-X 100 (Sigma-Aldrich), 1x PBS), blocked for 45
minutes (0.03% Triton-X 100, 10% Normal Donkey Serum (NDS), 1x PBS), and
incubated with primary antibodies (Table 1) for 1 hour at room
temperature and overnight at 4 degrees with agitation. After 48 hours,
primary antibodies were washed out (3 x 5 minutes 1x PBS) and secondary
antibodies were applied for 2 hours (Table 1). Secondary antibodies were
removed and the sections were washed with 1x PBS 5 x 15 minutes, and the
sections were mounted for confocal imaging.

  Cell Marker                                Round   Primary Antibody (dilution)           Secondary Antibody (dilution)              Catalog \#
  ------------------------------------------ ------- ------------------------------------- ------------------------------------------ -------------------------------------
  Cholecystokinin (CCK)                      1       rabbit anti-proCCK (1:500)            donkey anti-rabbit DyLight 405 (1:300)     Frontier Institute CCK-pro-Rb-Af350
  Calretinin (CalR)                          1       Guinea pig anti-Calretinin (1:1000)   donkey anti-guinea pig Alexa 568 (1:300)   Swant CRgp7
  Muscarinic Cholinergic M2 Receptor (M2R)   1       rat anti-M2R (1:2000)                 donkey anti-rat Alexa 647 (1:300)          Millipore MAB367

**Table 1. Antibodies used in this study**

**Confocal imaging**

The NikonA1 confocal microscope used to take multi-channel fluorescence
images (405 nm, 488 nm, 561 nm, and 640 nm excitation for the four blue,
green, red and far-red colored channels) of labeled tissue sections.
Image capturing yielded 2048 x 2048 pixel z-stacks of the entire depth
of the section (with 3 microns between each image of a given stack).
Z-stacks were viewable via Fiji (Image J).

**Registration of confocal images to *in vivo* z-stacks and
identification of immunopositivity/negativity**

Unique groupings of cells from *in vivo* z-stacks were used as visual
fiducial markers to manually register *ex vivo* 75 μm z-stack images in
ImageJ, section by section. Care was taken to note the orientation of
the *ex vivo* sections regarding the *in vivo* z-stack images. The
coordinates taken from the AOD software when selecting ROIs for imaging
against the *in vivo* z-stack were used to identify ROI ID numbers in
the registered immunohistochemically stained images.

**Assignment of subtype identity based on immunostaining**

Interneuron subtypes were assigned based on our previously published
classification strategy^84^.

*Cholecystokinin-expressing basket cells (CCKBCs):*

GCaMP-expressing, immunopositive for proCCK, immunonegative for CalR and
M2R^114,115^. This population represents perisomatic-targeting basket
cells.

*Interneuron-selective interneuron type 2 (ISI2s):*

GCaMP-expressing, immunonegative for CalR, CCK and M2R. Present in
*stratum radiatum* or *stratum lacunosum moleculare* CA1
layers^114,115^.

*Interneuron-selective interneuron type 3 (ISI3s):*

GCaMP-expressing, immunopositive for CalR and immunonegative for CCK and
M2R. Present in all CA1 layers^104,114--117^.

*VIP only (VIPos):*

GCaMP-expressing, immunonegative for CalR, CCK and M2R. Present in
*stratum oriens* or *stratum pyamidale* CA1 layers. This category
excludes CCKBCs, ISI2s, and ISI3s, and could possibly include long-range
projecting interneurons^118^, but it may not be restricted to them;
therefore we called them VIP only for our purposes in this study.

**Pre-processing of Ca2+ imaging data**

Raw movies of each individual cell were independently motion corrected
using the whole-frame cross registration implemented in the SIMA
package^119^. ROIs were then hand-drawn over each imaged cell.
Fluorescence was extracted from each drawn ROI using SIMA. DeltaF/F was
then computed based on a 1st percentile baseline value over a 30-second
rolling window. The DeltaF/F traces were then smoothed through an
exponential filter. For analyses involving velocity regressed traces,
each animal's velocity was smoothed with a Hanning filter and then
covariate linear regression was applied to each fluorescence trace where
the smoothed velocity is regressed out.

**Velocity Correlation**

For each cell, velocity correlation was computed by taking Pearson\'s
correlation coefficient between the DeltaF/F trace and animal's smoothed
velocity with a -5 second to +5 second rolling lag. The highest absolute
value of the correlation coefficients in this time window was taken as
the correlation coefficient for the respective cell.

**Peri-Event Time Averages (PETA)**

For salient stimuli such as reward distribution events, a PETA was
constructed around each event respectively. These events were z-score
and then averaged for each respective cell, and the z-scored averages
were visualized using a heatmap or an averaged line-plot. A Post-Pre
metric was calculated based on the z-scored averages as a measure for
event modulation.

**Spatial Tuning Curve**

The 2 meter belt was divided into 100 bins (\~2cm per bin). Spatial
tuning curves for each cell during an experiment were generated by
averaging the velocity regressed DeltaF/F in each bin, when the animal
was running over 5 cm/s.

**Place field calculation**

The spatial tuning curves were shuffled 1000 times by circularly
rotating the velocity-regressed ∆F/F tuning curves around the 100
spatial bins to generate a shuffle distribution. Cells were classified
as being spatially selective (i.e., having a place field) if the cell's
tuning curve exceeded the 95th percentile (or was under the 5th
percentile for negatively selective cells) with at least 10 sequential
bins over/under the cutoffs from the shuffled distribution. The number
of bins over or under the 95th and 5th percentile curves were visualized
and a cutoff of 30 consecutive bins was used for the figures as it
represented the approximate 75th quartile for place field width.

**Position Decoding**

A linear classifier support vector machine (SVM) was trained on the
tuning curves of each respective cell and the decoding error was
visualized. The SVM was cross-validated by training the classifier on
(n-1) laps and evaluating on the held-out lap and this process was
repeated such that every lap was utilized as a test lap at least once.
Furthermore, the SVM was bootstrapped with different numbers of input
cells and input laps to compare the position decoding efficacy between
genotypes.

**Cross-correlation Analysis**

Pairwise Pearson's correlation coefficients were computed for every pair
of cells within a given experiment. The average coefficient value for
each subtype was then calculated per trial pair and averaged across
mice.

**Transient Analysis**

Online Active Set method to Infer Spikes (OASIS) was used to generate
deconvolved spike events based on transient events within the dF/F
traces^96^. We understand that these are not true spiking events, and we
are using OASIS to algorithmically identify transient events rather than
interpret the deconvolved spiking events.

Frequency: The deconvolved events were filtered to remove events that
were within 10 frames. The total number of filtered deconvolved events
by OASIS was normalized by the total duration of the recording.

For each cell, the average transient was segmented within a 15-second
time window and computed by averaging along aligned deconvolved spike
times. If multiple detected events were within 10 frames, the events
were treated as a single transient. A cell was only used if there were
at least 3 detected transients within the total trace to exclude cells
without obvious GCaMP-Ca^2+^ dynamics. For each average transient, we
computed the median value based on 5 seconds pre peak transient and the
range 5-10 seconds after the peak transient (given that the average
transient took significantly less than 5 seconds to resolve) to act as
the baseline.

Amplitude: The difference between the max value of the average transient
and the baseline value was computed

**QUANTIFICATION AND STATISTICAL ANALYSIS**

**Statistical analyses**

Subtype comparisons were conducted using one-way ANOVA tests with
Tukey's range tests with correction for multiple testing if appropriate.
For comparisons across genotypes, unpaired t tests were utilized if the
data followed a normal distribution. Otherwise, for non-normally
distributed comparisons, the Mann-Whitney U test was used. \*, P \<
0.05, \*\*, P \< 0.01, \*\*\*, P \< 0.001. All data analysis and
visualization were done with custom software built on Python version
2.7.15.

**Code and Data availability:**

All data used in this report are publicly available. ASD genes and
mutations counts can be found in Supplementary table from Zhou et al^2^.
SCZ genes and mutation counts can be found in Supplementary table from
Singh et al^1^. DD/ID genes and mutations counts can be found in
Supplementary table from Kaplanis et al^56^. Genes implicated in UKBB
adult cognitive functions can be found in Supplementary table from Chen
et al^63^. Genes in 22q11.2 3Mb deletion and 1.5 Mb deletion can be
found in Karayiorgou et al^14^. Human brain cell atlas can be found
(<https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443>).
Individual genetics and phenotypic data used in our study from the SPARK
project are available through the Simons Foundation
(<https://base.sfari.org>) following the access procedures described in
Zhou et al. Gene constraint score from gnomAD can be found at
(<https://gnomad.broadinstitute.org/news/2024-03-gnomad-v4-0-gene-constraint/>).
Brainspan expression data can be found at (https://www.brainspan.org/).

The code used for our implementation of the Cell type bias with analysis
pipeline, and code to generating all figures and tables, is available on
GitHub: (<https://github.com/explorerwjy/CellTypeBias_VIP.git>).

**Acknowledgements**

S.A.H. is supported by the Burroughs Welcome grant. This work has been
supported by the grants NIMHR01MH124923 to D.V. and J.A.G,
NIMHR01MH124047 to J.A.G, Stavros Niarchos Foundation to J.A.G., and
NIMHR01MH124047 NIMHR01MH124867 NINDSR01NS133381 to A.L. We thank
members of the Gogos, Losonczy, and Vitkup labs for fruitful scientific
discussions.

**Author Contributions**

S.H., J.W., J.A.G., A.L., and D.V. conceived the study. J.W. and J.C.
performed genetic bias and computational analyses. S.H. and B.R.
collected and analyzed experimental data. S.H., J.W., J.A.G., A.L., and
D.V. wrote the manuscript.

**\
**

**Supplemental Tables and Figures**

  Figure   Statistics
  -------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  5C       *CCKBC:* CCK basket cells student's t-test: T statistic: -1.313, P = 0.02. *ISI3*: Interneuron-Selective Interneuron Type 3: Student's T test Statistic: -3.284, P = 0.001. *ISI2:* Interneuron-Selective Interneuron Type 2: Student's T test Statistic: 2.191, P = 0.03. *VIP+ Alone:* Student's T-test Statistic: 1.105, P = 0.3
  5E       Mann-Whitney test, t = 104931.0, P = $4 \times 10^{- 5}$.
  5G       Mann-Whitney test, t =24290.0, P = 0.04.
  5I       Mann-Whitney test, t = 42011.0, P = $2 \times 10^{- 13}$.
  5J       Mann-Whitney test, t = 54521.5, P = $1 \times 10^{- 3}$
  5L       *CCKBC*: cholecystokinin-expressing basket cells Mann-Whitney U test: P = 0.3, *ISI3*: interneuron-selective interneuron type 3 cells Mann-Whitney U test: P = 0.0029, *ISI2*: interneuron-selective interneuron type 2 cells Mann-Whitney U test: P = 0.003, *VIP-Alone*: Mann-Whitney U test: P = 0.04. For Mann Whitney U tests: median ± inter-quartile range.
  6A       Averaged per animal: Mann-Whitney U test: statistic = 3.5, P = 0.144
  6B       10 cells: Mann-Whitney U test: P = $6 \times 10^{- 10}$, 30 cells: Mann-Whitney U test: P = $8 \times 10^{- 11}$, 50 cells: Mann-Whitney U test: P = $3 \times 10^{- 16}$, 70 cells: Mann-Whitney U test: P = $4 \times 10^{- 18}$`<!-- -->`{=html}8, 90 cells: Mann-Whitney U test: P = $3 \times 10^{- 21}$.
  6C       5 laps: Mann-Whitney U test: P = $3 \times 10^{- 6}$, 10 laps: Mann-Whitney U test: P =$2 \times 10^{- 21}$, 15 laps: Mann-Whitney U test: P =$3 \times 10^{- 15}$.
  6D       0 fields: Mann-Whitney U test: p = 0.056, 1 field: Mann-Whitney U test: P = 0.03, 2 fields: Mann-Whitney U test: P = 0.1, 3 fields: Mann-Whitney U test: P = 0.06.
  6F       Welch's T-test comparisons: ISI2 vs VIPo: P = 0.007. ISI3 vs ISI3: P = 0.043. ISI3 vs VIPo: P = 0.04. VIPo vs ISI2: P = 0.007. VIPo vs ISI3: P = 0.04. Welch's T-test.

**Table S1: Statistical detail summary.**

Boxplots show the 25th, 50th (median), and 75th quartile ranges, with
the whiskers extending to 1.5 interquartile ranges below or above the
25th or 75th quartiles, respectively. Outliers are defined as values
extending beyond. the whisker ranges. For Welch's t tests: mean ± sem.
For Mann Whitney U tests: median ± inter-quartile range. n.s.: not
significant. \*: P \< 0.05. \*\*: P \< 0.01. \*\*\*: P \< 0.001.

**Table S2-S7: Human brain Cell type biases at cluster level for
different psychiatric disorders and traits, including SCZ, ASD, ASD with
ID, DD/ID, 22q11.2 and VNR- genes.**

**Table S8-S11: Human brain Cell type biases Bias contrast at
supercluster level of different comparisons, including ASD to SCZ, ASD
to ASD with ID, ASD to DD/ID, and VNR- to VNR+.**

![](media/image11.png){width="6.5in" height="2.7180555555555554in"}

**Figure S1: Full-scale IQ distribution among ASD probands in the SPARK
and ASC cohorts. A).** The distribution of FSIQ scores among ASD
probands from the SPARK and ASC cohorts. The histogram (blue bars) shows
the FSIQ distribution with an overlaid density plot (red line). A
two-component Gaussian mixture model was fitted to the data (blue and
orange dashed curves), revealing distinct subpopulations within the FSIQ
distribution. The vertical dashed line at FSIQ = 70 indicates a commonly
used threshold for intellectual disability (IQ \< 70). FSIQ phenotype
data of SPARK and ASC probands were obtained from SFARI base and
Satterstrom *et al*.^120^ **B).** Gene-specific FSIQ patterns across the
top 61 high-confidence ASD genes, represented by vertical lines.
Individual probands are represented as dots, colored by mutation type,
with likely gene-disrupting (LGD) mutations in red and damaging missense
(Dmis) mutations in blue. The green diamonds indicate the average FSIQ
in the cohorts for each gene. Genes are arranged along the x-axis in
based on *de novo* west P-values from Zhou *et al*.^2^.

![A graph of different colored lines AI-generated content may be
incorrect.](media/image13.png){width="6.5in"
height="4.196527777777778in"}

**Figure S2. Temporal expression patterns of disorder-associated genes
across human brain development stages.** Average expression levels
(log10 RPKM) of genes associated with different disorders across human
brain developmental stages, from embryonic to adult periods. Gene sets
include ASD genes (green), ASD with ID genes (orange), SCZ genes (blue),
and all human genes (dashed yellow). Error bars represent standard error
of the mean. SCZ genes show similar expression levels to ASD genes,
while ASD with ID genes show higher expression levels. The vertical
dashed line marks the transition from prenatal to postnatal periods.
Expression data were obtained from the BrainSpan Atlas of Human Brain
Development^112^.

![](media/image14.png){width="6.5in" height="2.93125in"}

**Figure S3. Workflow for computing cell type mutation biases from
disorder-associated genes.**

This schematic illustrates the updated two-step procedure used to
compute disorder mutation bias for each brain cell type. **(Left)**:
Starting from a comprehensive single-cell gene expression atlas of 461
human brain cell types, we first compute gene expression levels as TPM
(transcripts per million) for each gene across all cell types.
Expression specificity scores are then derived by normalizing each
gene's TPM in a given cell type relative to its total TPM across all
cell types. **(Right)**: For each disorder, cell type--specific mutation
bias is calculated by aggregating expression specificity scores across
all disorder-associated genes, weighted by their mutational burden. This
yields a composite mutation bias score for each cell type, reflecting
its relative vulnerability to rare pathogenic mutations in the disorder.

![](media/image16.png){width="6.5in" height="4.0125in"}

**Figure S4. Impact of number of genes included in the mutation bias
analysis.** Mutation bias correlations were calculated between the
primary gene set (61 genes) and gene sets of varying sizes (10--200
genes). For each disorder, cell type--specific biases were recalculated
while progressively increasing the number of genes included, starting
from the most significantly associated and adding progressively less
significant genes. Panels show results for: **A)** SCZ; **B)** ASD;
**C)** ASD with ID and **D)** DDD. Blue lines indicate the Spearman
correlation between the mutation bias profile of the expanded set and
the primary 61-gene set; red lines indicate the --log10 (P-value) of the
least significant gene included in the set.

![](media/image18.png){width="6.5in" height="5.779166666666667in"}

**Figure S5. Comprehensive analysis of mutation biases across brain cell
types in psychiatric disorders.** Box plots show the distribution of
mutation bias significance (--log₁₀(P)) for major brain cell type
superclusters in: **A)** ASD, **B)** ASD with intellectual disability
(ID), **C)** schizophrenia (SCZ), and **D)** developmental disorders and
intellectual disability (DD/ID). Each point represents a cell type
(cluster) within a given supercluster. P-values were calculated by
comparing the observed mutation bias of each cluster to a null
distribution generated from 10,000 random gene sets (see Methods). Red
diamonds indicate the average --log₁₀(P) for each supercluster, and red
vertical lines indicate the median. The red dashed line corresponds to
an FDR threshold of 0.05; the lighter red dashed line corresponds to FDR
= 0.1.

![](media/image20.png){width="6.5in" height="3.890972222222222in"}

**Figure S6. Impact of gene expression levels on ASD-SCZ cell type bias
correlation.** Effect of sequentially removing genes based on their
expression levels (based on brain span expression^112^) on the cell type
mutation bias correlation between ASD and SCZ. Red line shows
correlation changes when removing genes with lowest expression levels
first, while blue dashed line shows removal of highest expressing genes
first. Grey line and error bars show mean ± standard deviation for
random gene removal.

![](media/image22.png){width="6.5in"
height="5.911111111111111in"}**Figure S7. ASD and SCZ Mutation bias
towards different cell type superclusters.** **A-D.** Detailed
comparison of mutation biases in specific neuron superclusters between
ASD and SCZ, P-value was calculated based on Mann--Whitney U test with
FDR correction: **(A)** Hippocampal CA1-3 neurons. **(B)** Upper-layer
IT. **(C)** Deep-layer IT **(D)** Deep-layer corticothalamic and 6b
neurons. **(E)** Bar plot comparing mutation biases across major neuron
types between ASD (blue) and ASD with ID (orange). Asterisks indicate
significance level between groups (\*: P \< 0.05. \*\*: P \< 0.01.
\*\*\*: P \< 0.001, \*\*\*\*: P \< 0.0001, Mann--Whitney U test with FDR
correction). All neuronal cell types at supercluster level is displayed.
Error bars represent standard error across clusters from corresponding
supercluster. **(F)** CGE interneurons show similar bias in SCZ genes
compared to ASD with ID genes. **(G)** CGE interneurons show stronger
bias in VNR- genes compared to VNR+ genes (verbal-numerical reasoning
negative and positively associated genes). **(H)** CGE interneurons show
stronger bias in DD/ID genes compared to ASD genes.

![](media/image24.png){width="6.5in" height="2.2354166666666666in"}

**Figure S8. Mutation bias comparison across psychiatric disorders for
MGE interneurons and LAMP5-LHX6/Chandelier interneurons**. **A**:
Mutation bias distribution for MGE interneurons across different
cohorts. **B**: Mutation bias distribution for LAMP5-LHX6 and Chandelier
cells across different cohorts. P-values for pairwise comparison were
shown, Mann--Whitney U test with FDR correction.

![](media/image26.png){width="6.5in" height="6.410416666666666in"}

**Figure S9. Comprehensive analysis of mutation biases across brain cell
types in 22q11.2 deletion.** Box plots showing significance of mutation
biases for all major brain cell types (superclusters). Each bar
represents the bias distribution of the human adult brain neuron type
belonging to a supercluster, X-axis shows the significance of biases as
-log10(P-value). Red diamond shows average of each supercluster, red
vertical line shows medium of each supercluster.

![](media/image28.png){width="5.602083333333334in"
height="6.731944444444444in"}

**Figure S10. VIP+ Cell imaging and targeting strategy.**

**A.** Schematic for mouse crossing and viral injection and calcium
imaging setup. **B.** 3D reconstruction of Z-stack from VIP-cre mice
**C.** Schematic of 3D regions of interest (ROIs) imaged with AOD
2-photon microscopy. **D.** Example chessboard scan view of imaged ROIs.
**E.** Manual registration of *in vivo* Z-stack section with *ex vivo*
75 m transverse section for *post hoc* subtype classification with
immunohistochemistry. **F.** Histological strategy for defining
VIP-expressing cell subtypes. **G.** 3D distribution of all cells imaged
from both genotypes. **H.** Distribution of VIP+ subtypes across CA1
hippocampal layers. **I.** Number of cells recorded per mouse (n=4,
ISI3: ISI2: VIP+ Alone: CCKBC: ). **J.** Average velocity in cm per
second plotted by location on the belt of WT (black) and *Df(16)A^+/-^*
(brown) mice. **K.** Average percent occupancy on the belt by position
per genotype plotted by location on the belt of WT (black) and
*Df(16)A^+/-^* (brown) mice. **L.** Total laps across mice within
genotype (WT: 14.4 ± 1.1 laps, *Df(16)A^+/-^*: 16.2 ± 2.0 laps, Mann
Whitney U = 51.5, P = 0.2). **M**. Licking proportion of WT (black) and
*Df(16)A^+/-^* (brown). There is no significant difference in total
licks across the belt during Random Foraging in WT (*n* = 4) or
*Df(16)A^+/-^* (*n* = 4) depicted by box-and-whisker plot, WT: 1916 ±
304 licks, *Df(16)A^+/-^*: 2200 ± 489 licks, Mann-Whitney U Test, P =
0.4.

**References**

1\. Singh, T., Poterba, T., Curtis, D., Akil, H., Al Eissa, M., Barchas,
J.D., Bass, N., Bigdeli, T.B., Breen, G., Bromet, E.J., et al. (2022).
Rare coding variants in ten genes confer substantial risk for
schizophrenia. Nature *604*, 509--516.2. Zhou, X., Feliciano, P., Shu,
C., Wang, T., Astrovskaya, I., Hall, J.B., Obiajulu, J.U., Wright, J.R.,
Murali, S.C., Xu, S.X., et al. (2022). Integrating de novo and inherited
variants in 42,607 autism cases identifies mutations in new
moderate-risk genes. Nat. Genet. *54*, 1305--1319.3. Gilman, S.R.,
Chang, J., Xu, B., Bawa, T.S., Gogos, J.A., Karayiorgou, M., and Vitkup,
D. (2012). Diverse types of genetic variation converge on functional
gene networks involved in schizophrenia. Nat. Neurosci. *15*,
1723--1728.4. Gilman, S.R., Iossifov, I., Levy, D., Ronemus, M., Wigler,
M., and Vitkup, D. (2011). Rare de novo variants associated with autism
implicate a large functional network of genes involved in formation and
function of synapses. Neuron *70*, 898--907.5. Gandal, M.J., Haney,
J.R., Parikshak, N.N., Leppa, V., Ramaswami, G., Hartl, C., Schork,
A.J., Appadurai, V., Buil, A., Werge, T.M., et al. (2018). Shared
molecular neuropathology across major psychiatric disorders parallels
polygenic overlap. Science *359*, 693--697.6. Rees, E., Creeth, H.D.J.,
Hwu, H.-G., Chen, W.J., Tsuang, M., Glatt, S.J., Rey, R., Kirov, G.,
Walters, J.T.R., Holmans, P., et al. (2021). Schizophrenia, autism
spectrum disorders and developmental disorders share specific disruptive
coding mutations. Nat. Commun. *12*, 5353.7. Kaladjian, A., Jeanningros,
R., Azorin, J.-M., Anton, J.-L., and Mazzola-Pomietto, P. (2011).
Impulsivity and neural correlates of response inhibition in
schizophrenia. Psychol. Med. *41*, 291--299.8. Trevisan, D.A.,
Foss-Feig, J.H., Naples, A.J., Srihari, V., Anticevic, A., and
McPartland, J.C. (2020). Autism spectrum disorder and schizophrenia are
better differentiated by positive symptoms than negative symptoms.
Front. Psychiatry *11*, 548.9. Andreasen, N.C., and Olsen, S. (1982).
Negative v positive schizophrenia. Definition and validation. Arch. Gen.
Psychiatry *39*, 789--794.10. Forbes, N.F., Carrick, L.A., McIntosh,
A.M., and Lawrie, S.M. (2009). Working memory in schizophrenia: a
meta-analysis. Psychol. Med. *39*, 889--905.11. Green, M., and Walker,
E. (1985). Neuropsychological performance and positive and negative
symptoms in schizophrenia. J. Abnorm. Psychol. *94*, 460--469.12. Spek,
A.A., and Wouters, S.G.M. (2010). Autism and schizophrenia in high
functioning adults: Behavioral differences and overlap. Res. Autism
Spectr. Disord. *4*, 709--717.13. Schneider, M., Debbané, M., Bassett,
A.S., Chow, E.W.C., Fung, W.L.A., van den Bree, M., Owen, M., Murphy,
K.C., Niarchou, M., Kates, W.R., et al. (2014). Psychiatric disorders
from childhood to adulthood in 22q11.2 deletion syndrome: results from
the International Consortium on Brain and Behavior in 22q11.2 Deletion
Syndrome. Am. J. Psychiatry *171*, 627--639.14. Karayiorgou, M., Simon,
T.J., and Gogos, J.A. (2010). 22q11.2 microdeletions: linking DNA
structural variation to brain dysfunction and schizophrenia. Nature
Reviews Neuroscience 2010 11:6 *11*, 402--416.15. Rees, E., Kendall, K.,
Pardiñas, A.F., Legge, S.E., Pocklington, A., Escott-Price, V., MacCabe,
J.H., Collier, D.A., Holmans, P., O'Donovan, M.C., et al. (2016).
Analysis of intellectual disability copy number variants for association
with schizophrenia. JAMA Psychiatry *73*, 963--969.16. Miller, D.T.,
Adam, M.P., Aradhya, S., Biesecker, L.G., Brothman, A.R., Carter, N.P.,
Church, D.M., Crolla, J.A., Eichler, E.E., Epstein, C.J., et al. (2010).
Consensus statement: chromosomal microarray is a first-tier clinical
diagnostic test for individuals with developmental disabilities or
congenital anomalies. Am. J. Hum. Genet. *86*, 749--764.17. Ropers, H.H.
(2008). Genetics of intellectual disability. Curr. Opin. Genet. Dev.
*18*, 241--250.18. Vissers, L.E.L.M., Gilissen, C., and Veltman, J.A.
(2016). Genetic studies in intellectual disability and related
disorders. Nat. Rev. Genet. *17*, 9--18.19. Millan, M.J., Agid, Y.,
Brüne, M., Bullmore, E.T., Carter, C.S., Clayton, N.S., Connor, R.,
Davis, S., Deakin, B., DeRubeis, R.J., et al. (2012). Cognitive
dysfunction in psychiatric disorders: characteristics, causes and the
quest for improved therapy. Nat. Rev. Drug Discov. *11*, 141--168.20.
Harvey, P.D., and Strassnig, M. (2012). Predicting the severity of
everyday functional disability in people with schizophrenia: cognitive
deficits, functional capacity, symptoms, and health status. World
Psychiatry *11*, 73--79.21. Bowie, C.R., and Harvey, P.D. (2006).
Cognitive deficits and functional outcome in schizophrenia.
Neuropsychiatr. Dis. Treat. https://doi.org/10.2147/nedt.s24531.22.
Bowie, C.R., Reichenberg, A., Patterson, T.L., Heaton, R.K., and Harvey,
P.D. (2006). Determinants of real-world functional performance in
schizophrenia subjects: Correlations with cognition, functional
capacity, and symptoms. Am. J. Psychiatry *163*, 418--425.23.
Trubetskoy, V., Pardiñas, A.F., Qi, T., Panagiotaropoulou, G., Awasthi,
S., Bigdeli, T.B., Bryois, J., Chen, C.Y., Dennison, C.A., Hall, L.S.,
et al. (2022). Mapping genomic loci implicates genes and synaptic
biology in schizophrenia. Nature *604*, 502--508.24. Bryois, J., Skene,
N.G., Hansen, T.F., Kogelman, L.J.A., Watson, H.J., Liu, Z., Eating
Disorders Working Group of the Psychiatric Genomics Consortium,
International Headache Genetics Consortium, 23andMe Research Team,
Brueggeman, L., et al. (2020). Genetic identification of cell types
underlying brain complex traits yields insights into the etiology of
Parkinson's disease. Nat. Genet. *52*, 482--493.25. Duncan, L.E., Li,
T., Salem, M., Li, W., Mortazavi, L., Senturk, H., Shahverdizadeh, N.,
Vesuna, S., Shen, H., Yoon, J., et al. (2025). Mapping the cellular
etiology of schizophrenia and complex brain phenotypes. Nat. Neurosci.,
1--11.26. Yao, S., Harder, A., Darki, F., Chang, Y.-W., Li, A., Nikouei,
K., Volpe, G., Lundström, J.N., Zeng, J., Wray, N., et al. (2025).
Connecting genomic results for psychiatric disorders to human brain cell
types and regions reveals convergence with functional connectivity.
Nature Communications. https://doi.org/10.1101/2024.01.18.24301478.27.
Broglio, C., Martín-Monzón, I., Ocaña, F.M., Gómez, A., Durán, E.,
Salas, C., and Rodríguez, F. (2015). Hippocampal pallium and map-like
memories through vertebrate evolution. J. Behav. Brain Sci. *05*,
109--120.28. Rodríguez, F., López, J.C., Vargas, J.P., Broglio, C.,
Gómez, Y., and Salas, C. (2002). Spatial memory and hippocampal pallium
through vertebrate evolution: insights from reptiles and teleost fish.
Brain Res. Bull. *57*, 499--503.29. Buzsáki, G. (2002). Theta
oscillations in the hippocampus. Neuron *33*, 325--340.30. Pelkey, K.A.,
Vargish, G.A., Pellegrini, L.V., Calvigioni, D., Chapeton, J., Yuan, X.,
Hunt, S., Cummins, A.C., Eldridge, M.A.G., Pickel, J., et al. (2023).
Evolutionary conservation of hippocampal mossy fiber synapse properties.
Neuron *111*, 3802-3818.e5.31. Hamm, J.P., Peterka, D.S., Gogos, J.A.,
and Yuste, R. (2017). Altered Cortical Ensembles in Mouse Models of
Schizophrenia. Neuron *94*, 153-167.e8.32. Herrlinger, S.A., Rao, B.Y.,
Conde Paredes, M.E., Tuttman, A.L., Arain, H., Varol, E., Gogos, J.A.,
and Losonczy, A. (2024). Disorganized Inhibitory Dynamics and Functional
Connectivity in Hippocampal area CA1 of 22q11.2 Deletion Mutant Mice.
bioRxiv, 2024.04.28.591464.
https://doi.org/10.1101/2024.04.28.591464.33. Gogos, J.A., Crabtree, G.,
and Diamantopoulou, A. (2020). The abiding relevance of mouse models of
rare mutations to psychiatric neuroscience and therapeutics. Schizophr.
Res. *217*, 37--51.34. Diamantopoulou, A., and Gogos, J.A. (2019).
Neurocognitive and Perceptual Processing in Genetic Mouse Models of
Schizophrenia: Emerging Lessons. Neuroscientist *25*, 597--619.35.
Arguello, P.A., and Gogos, J.A. (2006). Modeling madness in mice: one
piece at a time. Neuron *52*, 179--196.36. Gordon, A., Forsingdal, A.,
Klewe, I.V., Nielsen, J., Didriksen, M., Werge, T., and Geschwind, D.H.
(2019). Transcriptomic networks implicate neuronal energetic
abnormalities in three mouse models harboring autism and
schizophrenia-associated mutations. Molecular Psychiatry 2019 26:5 *26*,
1520--1534.37. Siletti, K., Hodge, R., Mossi Albiach, A., Lee, K.W.,
Ding, S.-L., Hu, L., Lönnerberg, P., Bakken, T., Casper, T., Clark, M.,
et al. (2023). Transcriptomic diversity of cell types across the adult
human brain. Science *382*, eadd7046.38. Geschwind, D.H., and Flint, J.
(2015). Genetics and genomics of psychiatric disease. Science *349*,
1489--1494.39. Iossifov, I., O'Roak, B.J., Sanders, S.J., Ronemus, M.,
Krumm, N., Levy, D., Stessman, H.A., Witherspoon, K.T., Vives, L.,
Patterson, K.E., et al. (2014). The contribution of de novo coding
mutations to autism spectrum disorder. Nature *515*, 216--221.40.
Matson, J.L., and Shoemaker, M. (2009). Intellectual disability and its
relationship to autism spectrum disorders. Res. Dev. Disabil. *30*,
1107--1114.41. Woodberry, K.A., Giuliano, A.J., and Seidman, L.J.
(2008). Premorbid IQ in schizophrenia: a meta-analytic review. Am. J.
Psychiatry *165*, 579--587.42. Khandaker, G.M., Barnett, J.H., White,
I.R., and Jones, P.B. (2011). A quantitative meta-analysis of
population-based studies of premorbid intelligence and schizophrenia.
Schizophr. Res. *132*, 220--227.43. Creeth, H.D.J., Rees, E., Legge,
S.E., Dennison, C.A., Holmans, P., Walters, J.T.R., O'Donovan, M.C., and
Owen, M.J. (2022). Ultrarare coding variants and cognitive function in
schizophrenia. JAMA Psychiatry *79*, 963--970.44. Nuechterlein, K.H.,
Green, M.F., Kern, R.S., Baade, L.E., Barch, D.M., Cohen, J.D., Essock,
S., Fenton, W.S., Frese, F.J., 3rd, Gold, J.M., et al. (2008). The
MATRICS Consensus Cognitive Battery, part 1: test selection,
reliability, and validity. Am. J. Psychiatry *165*, 203--213.45.
Mesholam-Gately, R.I., Giuliano, A.J., Goff, K.P., Faraone, S.V., and
Seidman, L.J. (2009). Neurocognition in first-episode schizophrenia: a
meta-analytic review. Neuropsychology *23*, 315--336.46. Minshew, N.J.,
Goldstein, G., and Siegel, D.J. (1997). Neuropsychologic functioning in
autism: profile of a complex information processing disorder. J. Int.
Neuropsychol. Soc. *3*, 303--316.47. Hill, E.L. (2004). Executive
dysfunction in autism. Trends Cogn. Sci. *8*, 26--32.48. Aleman, A., and
Kahn, R. (2005). Strange feelings: Do amygdala abnormalities dysregulate
the emotional brain in schizophrenia? Prog. Neurobiol. *77*,
283--298.49. Baron-Cohen, S., Ring, H.A., Bullmore, E.T., Wheelwright,
S., Ashwin, C., and Williams, S.C. (2000). The amygdala theory of
autism. Neurosci. Biobehav. Rev. *24*, 355--364.50. Banker, S.M., Gu,
X., Schiller, D., and Foss-Feig, J.H. (2021). Hippocampal contributions
to social and cognitive deficits in autism spectrum disorder. Trends
Neurosci. *44*, 793--807.51. Lieberman, J.A., Girgis, R.R., Brucato, G.,
Moore, H., Provenzano, F., Kegeles, L., Javitt, D., Kantrowitz, J.,
Wall, M.M., Corcoran, C.M., et al. (2018). Hippocampal dysfunction in
the pathophysiology of schizophrenia: a selective review and hypothesis
for early detection and intervention. Mol. Psychiatry *23*,
1764--1772.52. Zierhut, K., Bogerts, B., Schott, B., Fenker, D., Walter,
M., Albrecht, D., Steiner, J., Schütze, H., Northoff, G., Düzel, E., et
al. (2010). The role of hippocampus dysfunction in deficient memory
encoding and positive symptoms in schizophrenia. Psychiatry Res. *183*,
187--194.53. Moberg, S., and Takahashi, N. (2022). Neocortical layer 5
subclasses: From cellular properties to roles in behavior. Front.
Synaptic Neurosci. *14*, 1006773.54. Jutla, A., Foss-Feig, J., and
Veenstra-VanderWeele, J. (2022). Autism spectrum disorder and
schizophrenia: An updated conceptual review. Autism Res. *15*,
384--412.55. Palmer, D.S., Howrigan, D.P., Chapman, S.B., Adolfsson, R.,
Bass, N., Blackwood, D., Boks, M.P.M., Chen, C.-Y., Churchhouse, C.,
Corvin, A.P., et al. (2022). Exome sequencing in bipolar disorder
identifies AKAP11 as a risk gene shared with schizophrenia. Nat. Genet.
*54*, 541--547.56. Kaplanis, J., Samocha, K.E., Wiel, L., Zhang, Z.,
Arvai, K.J., Eberhardt, R.Y., Gallone, G., Lelieveld, S.H., Martin,
H.C., McRae, J.F., et al. (2020). Evidence for 28 genetic disorders
discovered by combining healthcare and research data. Nature *586*,
757--762.57. Karczewski, K.J., Francioli, L.C., Tiao, G., Cummings,
B.B., Alföldi, J., Wang, Q., Collins, R.L., Laricchia, K.M., Ganna, A.,
Birnbaum, D.P., et al. (2020). The mutational constraint spectrum
quantified from variation in 141,456 humans. Nature *581*, 434--443.58.
Chang, J., Gilman, S.R., Chiang, A.H., Sanders, S.J., and Vitkup, D.
(2015). Genotype to phenotype relationships in autism spectrum
disorders. Nat. Neurosci. *18*, 191--198.59. Peça, J., Feliciano, C.,
Ting, J.T., Wang, W., Wells, M.F., Venkatraman, T.N., Lascola, C.D., Fu,
Z., and Feng, G. (2011). Shank3 mutant mice display autistic-like
behaviours and striatal dysfunction. Nature *472*, 437--442.60.
Rothwell, P.E., Fuccillo, M.V., Maxeiner, S., Hayton, S.J., Gokce, O.,
Lim, B.K., Fowler, S.C., Malenka, R.C., and Südhof, T.C. (2014).
Autism-associated neuroligin-3 mutations commonly impair striatal
circuits to boost repetitive behaviors. Cell *158*, 198--212.61. Van
Derveer, A.B., Bastos, G., Ferrell, A.D., Gallimore, C.G., Greene, M.L.,
Holmes, J.T., Kubricka, V., Ross, J.M., and Hamm, J.P. (2021). A role
for somatostatin-positive interneurons in neuro-oscillatory and
information processing deficits in schizophrenia. Schizophr. Bull. *47*,
1385--1398.62. Batiuk, M.Y., Tyler, T., Dragicevic, K., Mei, S.,
Rydbirk, R., Petukhov, V., Deviatiiarov, R., Sedmak, D., Frank, E.,
Feher, V., et al. (2022). Upper cortical layer-driven network impairment
in schizophrenia. Sci. Adv. *8*, eabn8367.63. Chen, C.-Y., Tian, R., Ge,
T., Lam, M., Sanchez-Andrade, G., Singh, T., Urpa, L., Liu, J.Z.,
Sanderson, M., Rowley, C., et al. (2023). The impact of rare protein
coding genetic variation on adult cognitive function. Nat. Genet. *55*,
927--938.64. Fawns-Ritchie, C., and Deary, I.J. (2020). Reliability and
validity of the UK Biobank cognitive tests. PLoS One *15*, e0231627.65.
Marshall, C.R., Howrigan, D.P., Merico, D., Thiruvahindrapuram, B., Wu,
W., Greer, D.S., Antaki, D., Shetty, A., Holmans, P.A., Pinto, D., et
al. (2017). Contribution of copy number variants to schizophrenia from a
genome-wide study of 41,321 subjects. Nat. Genet. *49*, 27--35.66.
Matsunami, N., Hadley, D., Hensel, C.H., Christensen, G.B., Kim, C.,
Frackelton, E., Thomas, K., da Silva, R.P., Stevens, J., Baird, L., et
al. (2013). Identification of rare recurrent copy number variants in
high-risk autism families and their prevalence in a large ASD
population. PLoS One *8*, e52239.67. Robin, N.H., and Shprintzen, R.J.
(2005). Defining the clinical spectrum of deletion 22q11.2. J. Pediatr.
*147*, 90--96.68. Tang, K.L., Antshel, K.M., Fremont, W.P., and Kates,
W.R. (2015). Behavioral and psychiatric phenotypes in 22q11.2 deletion
syndrome. J. Dev. Behav. Pediatr. *36*, 639--650.69. Swillen, A.,
Devriendt, K., Legius, E., Eyskens, B., Dumoulin, M., Gewillig, M., and
Fryns, J.P. (1997). Intelligence and psychosocial adjustment in
velocardiofacial syndrome: a study of 37 children and adolescents with
VCFS. J. Med. Genet. *34*, 453--458.70. Kepecs, A., and Fishell, G.
(2014). Interneuron cell types are fit to function. Nature *505*,
318--326.71. Chamberland, S., and Topolnik, L. (2012). Inhibitory
control of hippocampal inhibitory neurons. Front. Neurosci. *0*, 165.72.
Letzkus, J.J., Wolff, S.B.E., and Lüthi, A. (2015). Disinhibition, a
circuit mechanism for associative learning and memory. Neuron *88*,
264--276.73. Pfeffer, C.K., Xue, M., He, M., Huang, Z.J., and Scanziani,
M. (2013). Inhibition of inhibition in visual cortex: the logic of
connections between molecularly distinct interneurons. Nat. Neurosci.
*16*, 1068--1076.74. Fu, Y., Tucciarone, J.M., Espinosa, J.S., Sheng,
N., Darcy, D.P., Nicoll, R.A., and Stryker M (2014). A cortical circuit
for gain control by behavioral state. Cell *156*, 1139--1152.75.
Bennett, C., Arroyo, S., Berns, D., and Hestrin, S. (2012). Mechanisms
generating dual-component nicotinic EPSCs in cortical interneurons. J.
Neurosci. *32*, 17287--17296.76. Letzkus, J.J., Wolff, S.B.E., Meyer,
E.M.M., Tovote, P., Courtin, J., Herry, C., and Lüthi, A. (2011). A
disinhibitory microcircuit for associative fear learning in the auditory
cortex. Nature *480*, 331--335.77. Poorthuis, R.B., Enke, L., and
Letzkus, J.J. (2014). Cholinergic circuit modulation through
differential recruitment of neocortical interneuron types during
behaviour: Cholinergic circuit modulation via interneurons. J. Physiol.
*592*, 4155--4164.78. Froemke, R.C., Merzenich, M.M., and Schreiner,
C.E. (2007). A synaptic memory trace for cortical receptive field
plasticity. Nature *450*, 425--429.79. Lee, A.T., Cunniff, M.M., See,
J.Z., Wilke, S.A., Luongo, F.J., Ellwood, I.T., Ponnavolu, S., and
Sohal, V.S. (2019). VIP interneurons contribute to avoidance behavior by
regulating information flow across hippocampal-prefrontal networks.
Neuron *102*, 1223-1234.e4.80. Stark, K.L., Xu, B., Bagchi, A., Lai,
W.S., Liu, H., Hsu, R., Wan, X., Pavlidis, P., Mills, A.A., Karayiorgou,
M., et al. (2008). Altered brain microRNA biogenesis contributes to
phenotypic deficits in a 22q11-deletion mouse model. Nat. Genet. *40*,
751--760.81. Drew, L.J., Stark, K.L., Fénelon, K., Karayiorgou, M.,
Macdermott, A.B., and Gogos, J.A. (2011). Evidence for altered
hippocampal function in a mouse model of the human 22q11.2
microdeletion. Mol. Cell. Neurosci. *47*, 293--305.82. Zaremba, J.D.,
Diamantopoulou, A., Danielson, N.B., Grosmark, A.D., Kaifosh, P.W.,
Bowler, J.C., Liao, Z., Sparks, F.T., Gogos, J.A., and Losonczy, A.
(2017). Impaired hippocampal place cell dynamics in a mouse model of the
22q11.2 deletion. Nat. Neurosci. *20*, 1612--1623.83. Sun, Z., Williams,
D.J., Xu, B., and Gogos, J.A. (2018). Altered function and maturation of
primary cortical neurons from a 22q11.2 deletion mouse model of
schizophrenia. Transl. Psychiatry *8*.
https://doi.org/10.1038/s41398-018-0132-8.84. Geiller, T., Vancura, B.,
Terada, S., Troullinou, E., Chavlis, S., Tsagkatakis, G., Tsakalides,
P., Ócsai, K., Poirazi, P., Rózsa, B.J., et al. (2020). Large-Scale 3D
Two-Photon Imaging of Molecularly Identified CA1 Interneuron Dynamics in
Behaving Mice. Neuron *108*, 968-983.e9.85. Szalay, G., Judák, L.,
Katona, G., Ócsai, K., Juhász, G., Veress, M., Szadai, Z., Fehér, A.,
Tompa, T., Chiovini, B., et al. (2016). Fast 3D imaging of spine,
dendritic, and neuronal assemblies in behaving animals. Neuron *92*,
723--738.86. Katona, G., Szalay, G., Maák, P., Kaszás, A., Veress, M.,
Hillier, D., Chiovini, B., Vizi, E.S., Roska, B., and Rózsa, B. (2012).
Fast two-photon in vivo imaging with three-dimensional random-access
scanning in large tissue volumes. Nat. Methods *9*, 201--208.87.
Vancura, B., Geiller, T., Grosmark, A., Zhao, V., and Losonczy, A.
(2023). Inhibitory control of sharp-wave ripple duration during learning
in hippocampal recurrent networks. Nat. Neurosci. *26*, 788--797.88.
Somogyi, P., and Klausberger, T. (2005). Defined types of cortical
interneurone structure space and spike timing in the hippocampus. J.
Physiol. *562*, 9--26.89. Guet-McCreight, A., Skinner, F.K., and
Topolnik, L. (2020). Common principles in functional organization of
VIP/calretinin cell-driven disinhibitory circuits across cortical areas.
Front. Neural Circuits *14*, 32.90. Varga, C., Golshani, P., and
Soltesz, I. (2012). Frequency-invariant temporal ordering of
interneuronal discharges during hippocampal oscillations in awake mice.
Proc. Natl. Acad. Sci. U. S. A. *109*, E2726-34.91. Lapray, D.,
Lasztoczi, B., Lagler, M., Viney, T.J., Katona, L., Valenti, O.,
Hartwich, K., Borhegyi, Z., Somogyi, P., and Klausberger, T. (2012).
Behavior-dependent specialization of identified hippocampal
interneurons. Nat. Neurosci. *15*, 1265--1271.92. Katona, L., Lapray,
D., Viney, T.J., Oulhaj, A., Borhegyi, Z., Micklem, B.R., Klausberger,
T., and Somogyi, P. (2014). Sleep and movement differentiates actions of
two types of somatostatin-expressing GABAergic interneuron in rat
hippocampus. Neuron *82*, 872--886.93. Somogyi, P., Katona, L.,
Klausberger, T., Lasztóczi, B., and Viney, T.J. (2014). Temporal
redistribution of inhibition over neuronal subcellular domains underlies
state-dependent rhythmic change of excitability in the hippocampus.
Philos. Trans. R. Soc. Lond. B Biol. Sci. *369*, 20120518.94. Turi,
G.F., Li, W.-K., Chavlis, S., Pandi, I., Hare, J.O., Priestley, J.B.,
Grosmark, A.D., Liao, Z., Ladow, M., Zhang, J.F., et al. (2019).
Vasoactive Intestinal Polypeptide-Expressing Interneurons in the
Hippocampus Support Goal-Oriented Spatial Learning. Neuron, 1--16.95.
Dudok, B., Klein, P.M., Hwaun, E., Lee, B.R., Yao, Z., Fong, O., Bowler,
J.C., Terada, S., Sparks, F.T., Szabo, G.G., et al. (2021). Alternating
sources of perisomatic inhibition during behavior. Neuron *109*,
997-1012.e9.96. Friedrich, J., Zhou, P., and Paninski, L. (2017). Fast
online deconvolution of calcium imaging data. PLoS Comput. Biol. *13*,
e1005423.97. Buzsáki, G., and Moser, E.I. (2013). Memory, navigation and
theta rhythm in the hippocampal-entorhinal system. Nature Neuroscience
2013 16:2 *16*, 130--138.98. O'Keefe, J., and Dostrovsky, J. (1971). The
hippocampus as a spatial map. Preliminary evidence from unit activity in
the freely-moving rat. Brain Res. *34*, 171--175.99. Robinson, E.B.,
Samocha, K.E., Kosmicki, J.A., McGrath, L., Neale, B.M., Perlis, R.H.,
and Daly, M.J. (2014). Autism spectrum disorder severity reflects the
average contribution of de novo and familial influences. Proc. Natl.
Acad. Sci. U. S. A. *111*, 15161--15165.100. Gonzalez-Burgos, G., and
Lewis, D.A. (2012). NMDA receptor hypofunction, parvalbumin-positive
neurons, and cortical gamma oscillations in schizophrenia. Schizophr.
Bull. *38*, 950--957.101. Buzsáki, G., and Draguhn, A. (2004). Neuronal
oscillations in cortical networks. Science *304*, 1926--1929.102.
Möhler, H., and Rudolph, U. (2017). Disinhibition, an emerging
pharmacology of learning and memory. F1000Res. *6*, 101.103.
Lovett-Barron, M., and Losonczy, A. (2014). Behavioral consequences of
GABAergic neuronal diversity. Curr. Opin. Neurobiol. *26*, 27--33.104.
Luo, X., Guet-McCreight, A., Villette, V., Francavilla, R., Marino, B.,
Chamberland, S., Skinner, F.K., and Topolnik, L. (2020). Synaptic
mechanisms underlying the network state-dependent recruitment of
VIP-expressing interneurons in the CA1 hippocampus. Cereb. Cortex *30*,
3667--3685.105. Kamigaki, T., and Dan, Y. (2017). Delay activity of
specific prefrontal interneuron subtypes modulates memory-guided
behavior. Nat. Neurosci. *20*, 854--863.106. Goff, K.M., and Goldberg,
E.M. (2021). A role for vasoactive intestinal peptide interneurons in
neurodevelopmental disorders. Dev. Neurosci. *43*, 168--180.107.
Mossner, J.M., Batista-Brito, R., Pant, R., and Cardin, J.A. (2020).
Developmental loss of MeCP2 from VIP interneurons impairs cortical
function and behavior. Elife *9*.
https://doi.org/10.7554/eLife.55639.108. Goff, K.M., Liebergall, S.R.,
Jiang, E., Somarowthu, A., and Goldberg, E.M. (2023). VIP interneuron
impairment promotes in vivo circuit dysfunction and autism-related
behaviors in Dravet syndrome. Cell Rep. *42*, 112628.109. Murray, J.D.,
Anticevic, A., Gancsos, M., Ichinose, M., Corlett, P.R., Krystal, J.H.,
and Wang, X.-J. (2014). Linking microcircuit dysfunction to cognitive
impairment: effects of disinhibition associated with schizophrenia in a
cortical working memory model. Cereb. Cortex *24*, 859--872.110.
Ioannidis, N.M., Rothstein, J.H., Pejaver, V., Middha, S., McDonnell,
S.K., Baheti, S., Musolf, A., Li, Q., Holzinger, E., Karyadi, D., et al.
(2016). REVEL: An Ensemble Method for Predicting the Pathogenicity of
Rare Missense Variants. Am. J. Hum. Genet. *99*, 877--885.111. Samocha,
K.E., Kosmicki, J.A., Karczewski, K.J., O'Donnell-Luria, A.H.,
Pierce-Hoffman, E., MacArthur, D.G., Neale, B.M., and Daly, M.J. (2017).
Regional missense constraint improves variant deleteriousness
prediction. bioRxiv, 148353. https://doi.org/10.1101/148353.112. Kang,
H.J., Kawasawa, Y.I., Cheng, F., Zhu, Y., Xu, X., Li, M., Sousa, A.M.M.,
Pletikos, M., Meyer, K.A., Sedmak, G., et al. (2011). Spatio-temporal
transcriptome of the human brain. Nature *478*, 483--489.113.
Lovett-Barron, M., Kaifosh, P., Kheirbek, M.A., Danielson, N., Zaremba,
J.D., Reardon, T.R., Turi, G.F., Hen, R., Zemelman, B.V., and Losonczy,
A. (2014). Dendritic inhibition in the hippocampus supports fear
learning. Science *343*, 857--863.114. Acsády, L., Görcs, T.J., and
Freund, T.F. (1996). Different populations of vasoactive intestinal
polypeptide-immunoreactive interneurons are specialized to control
pyramidal cells or interneurons in the hippocampus. Neuroscience *73*,
317--334.115. Acsády, L., Arabadzisz, D., and Freund, T.F. (1996).
Correlated morphological and neurochemical features identify different
subsets of vasoactive intestinal polypeptide-immunoreactive interneurons
in rat hippocampus. Neuroscience *73*, 299--315.116. Tyan, L.,
Chamberland, S., Magnin, E., Camiré, O., Francavilla, R., David, L.S.,
Deisseroth, K., and Topolnik, L. (2014). Dendritic inhibition provided
by interneuron-specific cells controls the firing rate and timing of the
hippocampal feedback inhibitory circuitry. J. Neurosci. *34*,
4534--4547.117. Pelkey, K.A., Chittajallu, R., Craig, M.T., Tricoire,
L., Wester, J.C., and Mcbain, X.C.J. (2017). Hippocampal GABAergic
inhibitory interneurons. Physiol. Rev. *97*, 1619--1747.118.
Francavilla, R., Villette, V., Luo, X., Chamberland, S., Muñoz-Pino, E.,
Camiré, O., Wagner, K., Kis, V., Somogyi, P., and Topolnik, L. (2018).
Connectivity and network state-dependent recruitment of long-range
VIP-GABAergic neurons in the mouse hippocampus. Nat. Commun. *9*,
5043.119. Kaifosh, P., Zaremba, J.D., Danielson, N.B., and Losonczy, A.
(2014). SIMA: Python software for analysis of dynamic fluorescence
imaging data. Front. Neuroinform. *8*, 80.120. Satterstrom, F.K.,
Kosmicki, J.A., Wang, J., Breen, M.S., De Rubeis, S., An, J.-Y., Peng,
M., Collins, R., Grove, J., Klei, L., et al. (2020). Large-Scale Exome
Sequencing Study Implicates Both Developmental and Functional Changes in
the Neurobiology of Autism. Cell *180*, 568-584.e23.
