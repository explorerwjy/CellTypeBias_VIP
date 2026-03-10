**Revision To-Dos: (Summary)**



**Analyses:**



ID’ing CCKBC cells from ISIs based on our stim protocol. (waiting on raw data from Erica)

Add feature matching based null for cell type biases (Done)

Simulation of different cohort size and their impact on cell type implication. Main conclusion: ASD and DD are stable, even with 10% cohort size we can identify the correct cell types; SCZ seems not yet saturated, correlation with full set is bad when subsampled. In future SCZ needs more sequencing (Done)

Different gene set size, explain why SCZ has high correlation with main when adding low confident genes

Add negative control from non-brain traits (Done)

Add SCZ neg genes (Done)

Show CCKBC and other VIP (ISI types) has different VIP biases

Option 1: use SNCG, CR and M2R select clusters from human atlas and label them as CCK or ISI. For highly expressed 22q genes and 22q DEGs, CCK < ISI with marginal significance. (Done)

Option2: Mouse M1 dataset has some CCKBCs labelled, we can map Mouse patch-seq data to mouse cell atlas then to human atlas. Different human CGE clusters were mapped with these M1 patch-seq CCKBCs, and there’s no difference between CCKBC and VIP 22q bias with this method. I’ll spend more time fine tuning the mapping to see if they can converge to option 1. M1 dataset.

Show why capping max Spec at 2. I simulate cell type specificity based on expression and ZINB noise profile, show smaller cell type and gene with lower expression will have inflated specificity if not properly controlled. Also tried different caps show results stay the same with different reasonable cutoffs  (Done)

Provide all tables for geneweights/bias (Done)

Reorganize main figures to match our reply to Reviewer #2.



**Wetlab experimental:**

New figure: Ephys characterization study (Erica & Steph)

New supplemental figure: ‘CCKBC activity can’t be analyzed by waveform analysis’ evidence.



**Response Progress:**

R1

General comment summary: 0/1 (update to ‘completed’ phrasing, update the text)



R2

Major comments: /5

Minor comments: /5



R3

Major comments: /8

Minor comments: /2

PREAMBLE: We thank the reviewers for their helpful and detailed comments to improve our manuscript. To address the reviewers' concerns, we have performed additional experiments and analyses, added X figures,...



Reviewer #1: In this paper the authors take a multivariate approach to addressing the intersection of gene mutation to neuropsychiatric disease. They begin by taking human data on ASD from the Spark consortium and Scz data from the Schema consortium. Cross analyzing this data with analysis of human cell type gene expression data, they conclude a swath of cell types including eccentric MSNs, Lamp5-Lhx6 populations and chandelier neurons, CGE neurons, Upper and deep layer IT neurons and finally Midbrain inhibitory cells are strongly affected by the intersection of the genes implicated jointly from the Spark and Schema data sets. In addition they identify several super clusters including excitatory amygdala, hippocampus and both IT and ET populations and say these findings are consistent to with such brain disorders. In effect, they have implicated a large portion of the cell types in both pallial and subpallial structures as being candidate for being impacted in patients suffering from either ASD and SCZ. Going on from there they then focus down on MGE and CGE-derived cells with a particular focus on Lamp5-Lhx6/Chandelier, which comprise a tiny subpopulation of the broader interneuron subtypes. From this more general vantage, they switch their focus to the 22q11.2 deletion disorder and use the Df(16)A+/- mouse model to study this in mice. Using a largely motor based metric, they conclude that VIP interneurons are disproportionately affected in their function, leading them to generally conclude this implicates VIP interneurons as the underlying cause of certain forms of ASD/SCZ.



Taken together this work provides a large survey of genes, cell types and psychiatric conditions but feel more like a review article than a focused piece attempting to link gene defects to an underlying cause. While the work seems to be carefully assembled, the broad overview does little to uncover a distinct subtype being convincingly implicated in the disease rather than simply providing hints as to a large number of cell types that may be implicated. As such, while aspects of this work are of interest, they collectively short fall far of providing a convincing narrative to implicate any cell type, the VIP one investigated in their Df(16)A+/- mouse model in particular, as being centrally involved in ASD and SCZ disorders. Given the wide range of neuropsychiatric disorders they are grouping, this isn't terribly surprising and does not argue strongly at all that this work is suitable for publication in Neuron.



We sincerely thank the Reviewer for their thoughtful comments and for the opportunity to clarify and sharpen the narrative and goals of our study. Rare, highly penetrant variants are known to be strong drivers of cognitive dysfunction across neuropsychiatric disorders. The recent emergence of comprehensive human single-cell transcriptomic atlases from the Allen Institute provides an unprecedented opportunity to investigate the neuronal consequences of such rare mutations in a systematic and unbiased manner. While prior studies mapping genetic risk to brain cell types have largely focused on common-variant GWAS signals, the cell-type specificity of rare, high-impact mutations has remained comparatively unexplored.



To address this gap, we developed a computational framework to systematically identify human neuronal cell types that are selectively impacted by rare damaging mutations across neuropsychiatric phenotypes. Applying this framework across multiple large genetic datasets (including ASD from the SPARK consortium, schizophrenia from the SCHEMA consortium, developmental delay cohorts, UK Biobank cognitive phenotypes, and genes within the 22q11.2 deletion interval), we find that although multiple neuron types are affected, caudal ganglionic eminence (CGE)–derived GABAergic interneurons emerge as a consistent and specific nexus of cognitive phenotypes across disorders. This convergence distinguishes CGE-derived interneurons from both MGE-derived interneurons and excitatory projection neurons, and identifies a cross-disorder cellular axis of cognitive vulnerability that, to our knowledge, has not been previously recognized.



We agree with the Reviewer that the original manuscript placed disproportionate emphasis on atlas-level mapping and did not sufficiently foreground this central result. In the revised manuscript, we will substantially reorganized the Introduction and Results to streamline the narrative and pivot earlier toward the central finding that CGE-derived interneurons represent a shared locus of vulnerability for cognitive dysfunction across disorders, while the broader atlas-level analyses are repositioned to provide context and validation.



Within this CGE framework, our in vivo work focuses on CGE/VIP interneuron populations because they provide a mechanistically interpretable and experimentally accessible instantiation of the broader CGE lineage signal. Although numerically smaller, these interneurons play a disproportionate role in regulating cortical gain, contextual modulation, and information routing through disinhibitory circuit motifs. Importantly, the in vivo phenotypes we examine in the Df(16)A+/− mouse model reflect alterations in intrinsic excitability and in the encoding of spatial and reward-related information, rather than motor function per se. In the revised manuscript we clarify this point and reposition the mouse experiments as a circuit-level test of how CGE-lineage interneuron dysfunction may translate rare mutation burden into altered cognitive representations.



Taken together, the revised manuscript presents a more focused and hypothesis-driven framework. By integrating large-scale human genetics, single-cell transcriptomics, and circuit-level physiology, our study identifies CGE-derived interneurons as a convergent cellular substrate through which rare, high-impact mutations may disrupt cognitive function across neuropsychiatric disorders. We believe this clarified narrative more directly conveys the conceptual advance of the work and better aligns the study with the scope and expectations of Neuron.



Reviewer #2: Based on the hypothesis that cognitive deficits are common features of psychiatric disorders, the authors investigate whether there is a biased enrichment of genes identified as rare, highly penetrant genetic mutations and specific neuronal cell types. They found a correlation between cell type biases in Autistic Spectrum Disorder (ASD) with low IQ and Schizophrenia (SCZ). There is a bias towards particular cell types in ASD (striatal spiny neurons) and SCZ (MGE interneurons). CGE interneurons appeared to be highly enriched in genes associated with disorders involving cognitive deficits. They then examined how these CGE interneurons behave in a mouse model of intellectual disability/schizophrenia (modelling 22q11.2 deletion) and showed that some interneurons from the VIP subclass exhibit deficits in activity, suggesting a top-down disinhibition. Although the paper is very interesting and most of the experiments are solid, some points need to be addressed to strengthen the manuscript.



We thank the Reviewer for the kind words and for the careful and constructive evaluation of our manuscript. As detailed below, the Reviewer’s comments helped us clarify the presentation and further strengthen the mechanistic interpretation of our results.



Major points



1) There are some inconsistencies regarding the comparisons the authors provided; at least, the rationale for their absence is not included in the manuscript's narrative. Figure 2: I am missing the comparison of ASD with ID and ASD mutation bias for MGE interneurons and Medium spiny neurons (at least for MGE, the equivalent of Figure 2G for MGE). I think this is important since Figure 2 characterizes where the bias lies. Figure 3, missing the comparison between SCZ and ASD with ID. Again, it appears that the correlation is closer when SCZ and ASD with ID are compared; therefore, why not include this comparison? Figure 5, the authors would need to provide the calcium transients (frequency) for the other CGE subtypes, at least in the Supplementary information (similar to 5J).



2) Figure 3: Sometimes, it is unclear what the authors intend to convey. For example, it would be helpful to clarify the apparent discrepancy of MGE being more biased in Schizophrenia (Figure 2), considering Schizophrenia's stronger correlation to ASD with ID, and still showing a greater bias difference with CGE cells than MGE cells when comparing ASD with ID to ASD.



We sincerely apologize for the initial omission and will incorporate the full set of requested comparisons to clarify and strengthen the interpretation of our results. Specifically, when comparing ASD with SCZ, both MGE- and CGE-derived interneurons exhibit stronger mutation biases in SCZ. When comparing ASD with ASD with intellectual disability (ID), both MGE and CGE interneurons show stronger biases in ASD with ID; however, CGE-derived interneurons display substantially greater relative biases (P-value = 1x10-6), highlighting a more pronounced CGE contribution to ASD-associated cognitive phenotypes. Furthermore, when comparing SCZ to ASD with ID, CGE interneurons show comparable mutation burdens (P-value = 0.9), whereas MGE-derived interneurons are less biased in ASD with ID compared to SCZ (P-value = 7x10-3), confirming a preferential involvement  of CGE interneurons (rather than MGE) in phenotypes with prominent cognitive impairment, such as ASD with ID. In response to the Reviewer’s comments, these comparative results are now be explicitly highlighted and discussed in the revised manuscript.



For Figure 5, we are not able to calcium-transient frequency analyses for CGE-derived interneuron subtypes, currently available genetic targeting tools  (e.g., neurogliaform Lamp5⁺/Id2⁺ cells), or substantially overlap with VIP-expressing populations (e.g., Sncg⁺ cells). However, our results do demonstrate that among VIP+ interneurons, the most robust and consistent functional alterations are observed in interneuron-selective interneuron subtypes. We expanded iscussion these technical and interpretive limitations.



3) It is clear from this paper that the activity of VIP interneurons is altered. The authors state that this leads to disinhibition of pyramidal cells because VIP inhibits Somatostatin (SOM) and Parvalbumin (PV), resulting in impaired spatial coding. Although their model shows that this causes changes in the percentage and configuration of place cells, no experimental data support this claim. There is a paper from the team on BioRxiv, where they conduct a comprehensive analysis of PV, SOM, and a previous paper demonstrating that the dynamics of place cells are impaired (Nature Neuroscience, 2017). However, if they wish to argue that the changes observed in place cells in the place cell field are due to disinhibition of PV and SOM, it would be valuable to include a rescue experiment using optogenetics to activate the VIP population and then observe the responses of pyramidal cells and PV/SOM interneurons.



Reply Plan: We agree with the Reviewer that the causal link between altered VIP interneuron activity and disinhibition of pyramidal neurons should be carefully framed. In the current manuscript, we do not claim direct experimental of PV/SOM disinhibition; rather, we interpret our findings within a well-established circuit framework in which VIP interneurons preferentially inhibit SOM and PV interneurons. In the revision, we :

·   	Explicitly acknowledge this limitation in the Results and Discussion, and that our data demonstrate VIP interneuron dysfunction and impaired spatial codingare consistent with, but do not directly , a disinhibitory mechanism.

·   	Cite prior work (including our previous Nature Neuroscience study and related analyses of PV/SOM dynamics) as converg evidence supporting this circuit interpretation.



While we agree that optogenetic rescue experiments activating VIP interneurons would, in principle, provide a powerful test of causality, such experiments fall beyond the scope of the present study and will be clearly identified as an important direction for future work. Importantly, the interpretability of such rescue experiments in the context of a multi-gene 22q11 deletion is inherently limited: optogenetic activation of VIP interneurons may not fully restore normal circuit function when dozens of dosage-sensitive genes are simultaneously affected. Moreover, PV and SOM interneurons are themselves perturbed by the same deletion, further complicating causal attribution and necessitating extensive additional control experiments. Addressing these interacting circuit-level effects would require a substantial experimental expansion beyond the focus of the current manuscript. We explicitly discuss these limitations in the revised version.



5) In Figure 5L, there is a notable change in reward for several VIP interneuron populations; this data is very intriguing. What is the mechanistic hypothesis behind this? Do these VIP interneurons express Dopamine receptors?



We appreciate the Reviewer’s interest in the reward-related modulation of VIP interneuron activity and agree that these observations are intriguing. In the revised Discussion, we expanded the mechanistic framework to propose that these VIP interneurons integrate neuromodulatory and contextual signals associated with reward processing.



Available transcriptomic datasets indicate that multiple interneuron classes, including VIP populations, express dopamine receptors at varying levels, consistent with prior reports. However, the presence of dopamine receptor transcripts does not, by itself, establish a primary or exclusive dopaminergic mechanism underlying the reward-related effects observed here. VIP interneurons are known to be strongly modulated by diverse long-range inputs and neuromodulators (including cholinergic and serotonergic systems), and it is therefore likely that reward-related activity reflects convergent modulation rather than a purely dopamine-driven process.

Accordingly, in the revised manuscript we tempered the interpretation, noted the evidence for dopamine receptor expression in interneurons, and emphasized that defining the precise neuromodulatory mechanism underlying the observed reward sensitivity of VIP interneurons remains an important direction for future study.





6) Could the authors compare their VIP interneuron subtype findings with human sequencing data? This could help identify if one subtype is more impacted in 22q11.2 deletion syndrome.

Follow up: CGE cells have been classified using transcriptomics in both mouse and human studies, such as those from the Allen Brain Institute. In mice, the transcriptome, morphology, and electrophysiological characteristics have been analysed within the same cell. They could attempt to establish a correlation based on the electrophysiological properties already known for CGE interneurons in mice (e.g the equivalence of some of these interneuron function in hippocampus). Then, they could use mouse transcriptomics to create correlations in humans. I understand that most of the available data has been gathered from the neocortex, but many interneurons also share features with mirror neurons in the hippocampus because they are all part of the cerebral cortex.



We thank the Reviewer for this excellent suggestion, which directly motivated a substantial new analysis. Following the Reviewer's proposed strategy, we used both cross-species transcriptomic mapping and electrophysiological feature transfer to identify which human CGE clusters correspond to functionally distinct mouse VIP subtypes, and then tested whether these imputed subtypes show differential mutation bias.

**Cross-species transcriptomic mapping.** We integrated mouse patch-seq data from the Allen Institute (5,764 cells from visual and motor cortex, including 333 transcriptomically defined CCKBCs from the Sncg subclass) with the human Siletti et al. atlas (175,000 CGE interneurons, 21 clusters) using two independent integration methods: Harmony batch correction and scVI latent-space integration. In both methods, we projected mouse cells into the human reference space and assigned each to its nearest human cluster via k-nearest-neighbor voting (k = 30). Both methods converged on the same set of human clusters — 279, 280, and 281 — as the primary targets of mouse CCKBCs (Harmony: 67–96% CCKBC fraction; scVI: 86–90%), with significant concordance between methods (Spearman rho = 0.46, P = 0.03). These three clusters are annotated as VIP-expressing (INT-VIP) in the Siletti atlas.

**Electrophysiology signature transfer.** As the Reviewer suggested, we leveraged mouse patch-seq data in which transcriptomic identity and electrophysiology are measured in the same cell. We defined the mouse CCKBC electrophysiological fingerprint — higher ISI coefficient of variation, higher input resistance, and longer membrane time constant relative to other VIP types — and then scored 188 human CGE patch-seq cells (from Lee & Bhatt, mapped to Siletti clusters via scANVI) for similarity to this signature. The electrophysiology-based ranking of human clusters correlated positively with the transcriptomic CCKBC fraction from the Harmony mapping (Spearman rho = 0.49), providing independent cross-modal support for the cluster assignments.

**Mutation-bias comparison.** Based on these convergent results, we designated five human CGE clusters (277, 278, 279, 280, 281) as imputed "CCKBC-like" (Harmony CCKBC fraction > 50%) and compared their mutation-bias scores against the remaining 16 non-CCKBC CGE clusters across multiple disorders. For the 22q11.2 gene set, imputed CCKBCs showed slightly lower bias than non-CCKBC clusters, consistent in direction with the mouse electrophysiology finding, but this difference was not significant (Mann-Whitney P = 0.97) — likely reflecting the low statistical power inherent to the small 22q11.2 gene set. For other disorders, the pattern was reversed: imputed CCKBCs showed significantly *higher* bias for ASD (P = 0.015) and schizophrenia (P = 0.0002).

However, we urge caution in interpreting these subtype-level results, because the cross-species mapping of VIP subtypes faces fundamental biological limitations. CCKBCs are defined in mouse as a discrete subtype within the Sncg subclass, but no corresponding discrete cluster exists in human transcriptomic atlases. Cross-species classification of Sncg subtypes yields near-chance accuracy (AUROC = 0.50; Bakken et al., Nature 2021), the worst of any interneuron subclass, and CCK mRNA is expressed ubiquitously across human cortical neurons (Darmanis et al., PNAS 2015), precluding its use as a discriminative marker. Given these constraints, we conclude that the broader CGE/VIP supercluster — rather than individual subtypes — represents the most robust unit for cross-species integration of mutation-bias signals. We present this analysis in the Supplementary Information and discuss the biological limitations of cross-species VIP subtype mapping in the revised Discussion.











Minor



Figure 1 could be a nice figure for the Graphical Abstract, but not as a main figure of the manuscript.



Accordingly, we will use Figure 1 as the Graphical Abstract, and remove it from the main text.



Some sentences are overstated; the overstatements need to be toned down. E.g., Line 220 showed the strongest differential biases towards SCZ compared to ASD. Line 180: Is it appropriate to identify as 'strong correlation' when R=0.68?



We agree with the Reviewer and will soften the language regarding “strong correlations” (e.g., R = 0.68) and similar phrasing throughout the manuscript.



Figure 2D can go to the supplementary; still, it will make the point for specificity, while the paper seems to start to focus on other, more relevant cell populations.



Reply Plan: In the resubmission, we will make a clear distinction between the similarities and differences in neuronal impact across disorders. We will substantially reorganize the Introduction and Results to streamline the narrative and pivot earlier to the central finding that CGE-derived interneurons represent a shared locus of vulnerability for cognitive dysfunction across disorders.





We agree with the reviewer’s suggestion. We find it intriguing that while cortical excitatory neurons are targeted at similar levels in both disorders, striatal Medium Spiny Neurons (MSNs) show a distinct specificity, being significantly more biased in ASD compared to SCZ. However, we recognize that this divergence is a separate story from our central narrative, which focuses on the convergence of cognitive deficits onto CGE interneurons. To maintain the narrative focus, we have moved this analysis to the Supplementary Information (now Figure S7) as suggested.



We thank the reviewer for this suggestion. While we understand the motivation to streamline the presentation, we believe that Figure 2D represents an important result that helps establish the broader landscape of neuronal vulnerability across disorders and therefore should remain in the main text. Rather than moving the figure to the Supplementary Information, we have revised the surrounding text to more clearly focus the narrative and guide the reader toward the central conclusion of the study. In particular, the Introduction and Results have been reorganized to distinguish similarities and differences in neuronal impact across disorders, while pivoting earlier to the key observation that CGE-derived interneurons represent a shared locus of vulnerability for cognitive dysfunction.

As illustrated in Figure 2D, cortical excitatory neurons show comparable levels of mutation bias in both ASD and schizophrenia, whereas striatal Medium Spiny Neurons (MSNs) display a more selective enrichment in ASD relative to SCZ. We now clarify in the text that this divergence, although biologically interesting, is not the primary focus of the manuscript. Instead, the figure serves to contextualize the broader analysis and highlight the contrast with the more consistent cross-disorder signal observed in CGE-derived interneurons. By retaining the figure while sharpening the narrative around it, we believe the manuscript more clearly conveys both the overall cellular landscape and the specific convergence that motivates our subsequent functional experiments.





Please add the "n" of mice in the experiments.





Feel like a jump from Figure 5 to 6, can explain the rationale behind it more clearly.



Reply Plan: We thank the Reviewer for the comment. We will add specific "n" numbers for all mouse experiments to the figure legends and clarify the rationale for the transition between Figure 5 (activity) and Figure 6 (spatial coding).



Reviewer #3: Vitkup et al. present a framework that integrates rare genetic variant association data with human single-cell transcriptomic profiles to infer which brain cell types are most impacted in ASD, schizophrenia, and related cognitive phenotypes. They further test predictions in a 22q11.2 deletion mouse model, focusing on VIP+ CGE interneurons using state of the art techniques. The manuscript is ambitious, well-written, and addresses a timely question. However, I have several major concerns about the validation of the computational framework, the genetic assumptions underlying their analyses, and the interpretation of results.

This type of manuscript is hard to write as a scientist and hard to interpret as a reader. The reason is the necessary narrowing of scope (global cell type analysis to functional interrogation of a single cell class) always leaves you to wonder if the central hypothesis of the first part has really been tested or simply one out of several possibilities. I think the authors do a decent job of this hard transition, and importantly, I do not want to ask for major experimental work. This means I urge the authors to be careful in summarizing the take home message. Therefore, I think the discussion needs a paragraph on non-VIP specific effects in the 22q11 mouse. What would happen if they recorded from SST, PV, SNCG or LAMP5 interneurons instead? Could the change of VIP function be related to changes in excitatory inputs (are they enriched)? The reduced correlation with movement would certainly suggest that they receive weaker inputs.



Beyond this I have focused my review on the computational parts since this reflects my expertise.



We thank the Reviewer for the careful, balanced, and insightful evaluation of our manuscript. We appreciate the recognition of the ambition and technical rigor of the work, as well as the thoughtful guidance on how to more clearly articulate the scope and interpretation of our findings. We fully agree that manuscripts bridging large-scale computational inference with targeted circuit-level interrogation must be especially careful in framing their conclusions. In the revision, we expanded the discussion to explore non-VIP effects in circuit, including the results of our previous work examining other interneuron subtypes. Below, we address each of the Reviewer’s points in turn.



Major points:

1. Validation of the computational pipeline

The permutation framework needs stronger justification. Matching only on gene set size is insufficient; at minimum, random gene sets should be matched on gene length, conservation, and expression level. The authors should more comprehensively validate whether their "mutation bias" statistic is robust to such gene-level confounds.



We thank the reviewer for this insightful comment. We agree that it is important to carefully justify the permutation framework and to evaluate potential gene-level confounds.



First, several major confounders were already explicitly incorporated in our pipeline. For schizophrenia, the expected mutation counts used in the burden analyses already accounted for gene length, and we have now extended this expected-count normalization to the ASD and NDD gene sets as well. Expression level is controlled at two levels: (i) TPM-based normalization across cell types (Equation 1) and (ii) gene-level specificity scores that normalize expression across all profiled cell types (Equation 2). Gene constraint and conservation were not treated as independent nuisance variables, as constraint scores were strongly correlated with cross-disorder mutation bias in our data, indicating that they capture a biologically meaningful component of disease architecture.



To directly address the reviewer’s suggestion, we have now implemented an alternative permutation framework in which random gene sets are matched to the disorder genes on gene length, LOEUF, and brain expression level, in addition to overall mutation load. This alternative null model tested whether genes with similar basic genomic and expression properties, and with the same total number of mutations, recapitulate the observed cell-type mutation biases. As expected, this more stringent matching reduced the nominal significance of several traits, likely because many genes with matched LOEUF and brain expression are themselves enriched for CNS-related functions. Importantly, however, the key comparative conclusions of the paper remained unchanged: the relative ranking of cell types, the convergence on CGE interneurons for disorders with cognitive impairment, and the cross-disorder contrasts were consistent across both the original and matched null models.



We now describe both null models and their rationale in the Methods section, and report FDR values under each model for all traits in a new supplementary table (Supplementary Table X). In addition, as a further robustness check, we added positive control gene sets with well-established cell-type preferences and non-brain negative control traits which support the specificity of the framework (see response to point 3 and new Supplementary Figure X).





Reviewer figure X: QQ plots of cell-type bias p-values under random and confounder-matched null models. Each panel shows a quantile-quantile plot of observed versus expected -log(p-value) for cell-type bias across 461 cell types. P-values were derived from 10,000 permutations under two null models: random gene sampling (gray) and Best-of-N matching (blue; N = 1,000 candidate sets matched on gene length, whole-brain expression, and evolutionary conservation). The top row displays psychiatric and neurodevelopmental disorder gene sets (SCZ, ASD, and DDD), and the bottom row includes 22q11.2 deletion syndrome alongside two non-brain negative control traits (HDL cholesterol and  alanine aminotransferase, italic labels). For all disorder gene sets, p-values deviate strongly above the diagonal, indicating robust cell-type enrichment that persists even after confounder matching. By contrast, negative control gene  sets fall along the diagonal, confirming that the analysis does not produce spurious enrichment for traits without expected neuronal specificity. Dashed line indicates the null expectation.





2. Gene set size and stability of estimates

In Figure S4, mutation bias correlations remain stable even when the SCZ gene set is expanded to 200 genes. This is surprising: including genes without true association should add noise. Why is this not the case? The statistical explanation for this stability needs to be clarified. More generally, the manuscript should explicitly address how the number of included genes affects the stability of their estimates. Why should the number of genes be selected based on number of ASD rare-variant implicated genes?



We thank the Reviewer for raising this important point and agree that the stability of mutation-bias correlations as the schizophrenia gene set expands requires clarification. To address this directly, we performed three new analyses that together explain the observed stability and justify the choice of N = 61 genes.

First, the stability primarily reflects the fact that added SCZ genes (ranks 62–200) carry concordant cell-type signal, not random noise. When we computed mutation bias separately for the top-61 genes and the added genes (62–200), the two bias profiles were strongly correlated across all 461 cell types (Spearman ρ = 0.60, P = 4.4 × 10⁻⁴⁷; ρ = 0.48 for neuronal cell types alone). In contrast, a size-matched set of 139 random genes showed near-zero correlation with the top-61 profile (ρ = −0.08, P = 0.08). Both the top-61 and added gene subsets implicated MGE interneurons and amygdala excitatory neurons among their top-ranked superclusters, confirming that genes beyond rank 61 contribute structured, biologically coherent signal rather than noise. A sliding-window analysis further demonstrated that 61-gene windows starting at progressively lower ranks maintain ρ > 0.30 with the top-61 reference, even at the lowest-ranked genes (Reviewer Figures R3.2a,b).

Second, our mutation-bias statistic uses gene-level mutation counts as weights (Equations 4–6), which further stabilizes estimates. When we performed an unweighted analysis in which all genes contributed equally, we observed a substantially larger drop in correlation as the gene set expanded, consistent with the Reviewer’s expectation that weakly associated genes introduce noise (Figure S4). However, because SCZ gene weights are relatively flat (mean weight 6.8 for added genes vs 8.9 for the top-61), weighting alone cannot explain the stability—the concordant cell-type signal described above is the primary factor.

Third, although the correlation (i.e., the cell-type pattern) is preserved when expanding the gene set, the effect size magnitude decreases asymmetrically across disorders. At N = 61, SCZ, ASD, and DDD show comparable CGE interneuron effects (0.150, 0.156, and 0.167, respectively), enabling fair cross-disorder comparison. Beyond N = 61, the SCZ effect dilutes markedly (to 0.088 at N = 200, a 41% drop), while ASD and DDD retain substantially more signal (0.124 and 0.160, respectively). This asymmetry reflects the inherently noisier case–control design underlying SCZ gene discovery compared to de novo-based designs for ASD and DDD. Neuronal effect discrimination (standard deviation across neuronal cell types) shows the same pattern: SCZ discrimination degrades faster than ASD/DDD as the gene set expands, reducing the ability to distinguish disorder-specific cell-type profiles (Reviewer Figure R3.2c).

We therefore selected N = 61—the number of exome-wide significant ASD genes—as the default gene-set size because it represents the regime where (i) all three disorders have comparable effect sizes, enabling meaningful cross-disorder comparison, and (ii) the cell-type signal is maximally discriminative. Figure S4 is now explicitly framed as demonstrating that the main conclusions are not specific to this particular gene-set size but remain qualitatively stable across a broad range, while the choice of N = 61 is justified as the most informative operating point for cross-disorder analyses.





3. Impact of genetic architecture

The authors should explain how differences in genetic architecture (rare de novo-driven in ASD vs. modest case-control burden in SCZ) impact the statistical properties of the mutation bias metric. What is the minimum signal needed for robust inference? Is there a minimum sample size of the rare-variant genetic investigation at which the mutational bias is stable? Consider including a set of negative control traits, such as blood traits or other rare-variant investigations into non-brain related traits.



We thank the Reviewer for this important question regarding the impact of genetic architecture on the statistical behavior of the mutation-bias metric. To directly address this issue, we performed simulation-based power analyses and empirical downsampling studies designed to model realistic rare-variant burden scenarios.



In these simulations, we independently varied (i) the fraction of truly associated genes, (ii) sample size and resulting mutation counts, (iii) case–control vs. trio-based designs, and (iv) the degree of cell-type-specific convergence. These analyses quantify how each of these factors influences the stability of the inferred cell-type biases and the minimum signal needed for robust inference (new Supplementary Note and Supplementary Figure X).



The results reveal a marked stability for the de novo-driven disorders (ASD and NDD). For ASD, even when downsampling the cohort to just 10% of the full size, the cell-type bias correlation remained high (e.g., r > 0.8), despite the fact that 52% of the genes selected in the subsample differed from the full gene set. This stability reflects  three factors: first, risk genes show highly correlated cell-type expression profiles, such that different subsets of risk genes tend to be enriched in the same cell types (e.g., medium spiny neurons, cortical projection neurons); second, a consistent "core" signal, from the top genes with the highest mutation burden. persists  across downsamplings, driving the core biology; and finally, our use of gene-level weights, such that  that genes with more mutations (the strong signal) are more likely to be retained, further stabilizing the overall cell-type correlation. In contrast, the schizophrenia (SCZ) analysis, based on a case-control design, showed reduced stability at lower effective sample sizes (e.g., r = 0.67 at 10% downsampling). This reduced stability is expected because the case-control studies inherently have more noise than de novo studies, and the lower gene overlap (only 12% at 10% sampling) indicates genuine signal instability at lower sample sizes.



In addition, to further evaluate specificity, we applied the full pipeline to non-brain, negative control, traits derived from rare-variant burden analyses in UK Biobank (e.g., HDL cholesterol, inflammatory bowel disease, HbA1c) alongside neurological traits such as Parkinson’s and Alzheimer’s disease. These control traits show distinct and generally weaker brain cell-type mutation-bias profiles compared to ASD, schizophrenia, and NDD, both in terms of the magnitude and the anatomical distribution of biases (new Supplementary Figure X). Together, these analyses demonstrate that the mutation-bias framework behaves as expected under differing genetic architectures and that the prominent CGE and related interneuron signals observed in cognitive disorders are not generic features of rare-variant traits.

We expanded the Results and Methods to explicitly discuss the dependence of inference stability on genetic architecture, sample size, and study design, thereby clarifying the conditions under which robust cell-type conclusions can be drawn.

Reviewer figure X. Cell-type mutation bias estimates as a function of cohort sample size. Downsampling analysis showing stability of cell-type bias across three disorders as a function of genetic sample size. For each disorder, raw mutation counts were subsampled at fractions of 10–100% of the full cohort using binomial sampling, gene discovery was re-performed on the subsampled data (Poisson test for de novo cohorts, Fisher exact test for case-control), and the top 61 genes were selected to compute cell-type bias across 461 cell types. This procedure was repeated for 100 independent iterations at each fraction. Left y-axis (colored lines): Pearson correlation between the 461-dimensional bias vector from subsampled data and the full-data bias vector. Right y-axis (gray dashed lines): fraction of top-61 genes overlapping with the full-data gene set. Lines and ribbons indicate mean ± s.d. across 100 iterations. Horizontal dashed line marks r = 0.9. (Left) ASD de novo mutations (N = 42,607 probands): bias correlation exceeds 0.9 even at 10% sampling (r = 0.86), despite only ~49% gene overlap, indicating that different subsets of ASD risk genes converge on similar cell-type expression patterns. (Center) SCZ case-control mutations (N_case = 24,248): correlation improves more gradually with sample size (r = 0.83 at 90%), consistent with the noisier case-control design and lower gene overlap (~12% at 10%). (Right) DDD de novo mutations (N = 31,058 probands): highest stability, with correlation exceeding 0.9 at all fractions tested and ~68% gene overlap at 10% sampling.







4. Treatment of cognitive phenotypes (VNR)

The inclusion of verbal-numerical reasoning (VNR) is a good idea. However, splitting genes by positive vs. negative effect sizes is misleading. Due to negative selection, large-effect positive alleles on cognition are extremely rare. Thus, VNR(+) is likely a "null" set without signal. The authors do not split ASD or SCZ by direction, so why should VNR be split?



We thank the reviewer for this thoughtful comment and agree that the interpretation of the VNR(+) gene set requires careful consideration, given the strong negative selection acting on cognition-enhancing alleles.

VNR is a quantitative trait with a defined population mean. In our rare-variant burden framework, VNR(+) genes are those for which rare coding mutations are associated with higher VNR scores on average, whereas VNR(−) genes are associated with lower scores. As the Reviewer correctly points out, large-effect rare alleles that truly increase cognition are expected to be extremely rare; therefore, the VNR(+) gene set likely reflects a mixture of near-neutral and mildly protective effects rather than strong cognition-enhancing variants. We therefore used VNR(+) primarily as an internal control: under the hypothesis that CGE interneurons are preferentially impacted by variants that impair cognitive performance, one would expect VNR(+) genes to show reduced bias toward CGE interneurons. Consistent with this expectation,  VNR(−) genes exhibit a strong positive CGE mutation bias, whereas VNR(+) genes show a significant depletion in CGE interneurons (Figure X).

For ASD and NDD de novo datasets, the study design does not permit a directional split into risk-increasing versus protective variants. For schizophrenia, however, a directional split is in principle possible. To explore this, we identified the subset of genes with the strongest nominal “protective” direction of effect (more mutations in controls than in cases, OR < 1) in the SCHEMA dataset (n = 61) and computed their mutation-bias profile. As expected for putative protective or neutral genes, these genes exhibit near-zero or weak biases across other neuronal classes, with positive bias toward non neuronal cell types (new Supplementary Figure X). Interestingly, CGE IN showed the strongest depleted bias among these genes.

We now clarify in the text that the VNR(+) set was included primarily as a control analysis whereas the principal  biological conclusions are based on VNR(−) and the psychiatric disorder gene sets. We also incorporated the analysis of SCZ genes with OR < 1 into the relevant Results section to further illustrate this point.





5. Interpretation of differential bias tests

The authors highlight CGE interneurons as differentially biased across disorders, but many comparisons rely on small numbers of significantly associated clusters. For example, in ASD without ID, only 7 cell types pass significance, none of which are CGE. How do the authors justify making differential claims when one side of the comparison lacks clear signal?



We apologize for not making this distinction clearer in the original manuscript. Our framework addresses two related but conceptually distinct questions:

First, within-disorder bias asks whether a given cell supercluster significantly biased for a particular disorder compared to random gene sets with the same mutation load.

Second, between-disorder contrast asks, given two disorders with well-characterized high-penetrance gene sets, which cell superclusters show the largest difference in mutation bias between them?

The significance thresholds reported in Figures 2 and 3 pertain to the first analysis (within-disorder enrichment) In contrast, the differential bias analyses (e.g., contrasting ASD with ID vs. ASD without ID, or ASD vs. schizophrenia) address the second question by evaluating differences in the estimated mutation-bias values themselves, independent of whether each disorder’s bias in that cell type individually reaches a significance threshold.

In this context, the fact that CGE interneurons are not significantly biased in ASD without ID is itself informative: CGE interneurons show strong positive bias for disorders chacterized by cognitive impairment (ASD with ID, NDD, schizophrenia), but much weaker or absent bias in ASD without ID. The cross-disorder contrasts therefore highlight CGE interneurons as a cell population whose relative involvement tracks the presence of cognitive deficits across disorders.

To avoid potential confusion, we revised the relevant sections of the Results and Methods to more clearly distinguish these two analytical frameworks and to explicitly clarify that CGE interneurons are not significantly enriched in ASD without ID in the within-disorder analysis..



6. 22q11.2 analyses

In Figure S9, it is unclear which cell types remain significant after multiple-testing correction. The authors should explicitly mark FDR-adjusted significance. How about excitatory neurons? Why are they not included in the enrichment analysis to choose which cell type to focus on? It should be possible to connect the interneuron subtypes found in the electrophysiology dataset to interneurons found associated within the mutational bias framework. Is the difference between CCKBC interneurons and ISI2, ISI3 and VIPo interneurons recapitulated in the mutation bias analysis? This would strengthen the claims of convergence between electrophysiological data and genetic data.





We thank the Reviewer for this comment and agree that the presentation of the 22q11.2 analyses required clarification.



We agree that the statistical interpretation of Fig. S9 required explicit clarification. In the revised figure and legend, we now clearly indicate that none of the individual clusters passed FDR correction for the 22q11.2 gene set. This absence of FDR-significant clusters is expected given the small size and heterogeneous pathogenicity of the 22q11.2 interval. As demonstrated in multiple studies, only a subset of the 42 genes contributes to disease risk, resulting in substantially lower signal compared to ASD, SCZ, or NDD gene sets, in which each gene carries strong statistical evidence of association. Accordingly, the 22q11.2 analysis is inherently underpowered for genome-wide inference, and we now explicitly note this in the manuscript.



Because of this limitation, our goal in the 22q11.2 section was not to identify new cell types de novo, but rather to test a genetically motivated, cross-dataset hypothesis: that CGE/VIP interneurons, implicated across ASD, SCZ, and NDD rare-variant analyses, are also among the most affected populations within the 22q11.2 locus. The in vivo imaging data were therefore designed as a focused experimental test of this hypothesis rather than as an unbiased screen.



The reviewer is correct that excitatory neurons show elevated mutation-bias scores in the 22q11.2 gene set. However, these cell types also show comparably high bias across all disorder gene sets we examined—including ASD, SCZ, NDD, and even ASD without ID—suggesting that excitatory neurons represent a broad and non-specific vulnerability across many neurodevelopmental conditions. This widespread enrichment makes excitatory neurons less informative for distinguishing disorder-specific or phenotype-specific mechanisms.



In contrast, CGE/VIP interneurons exhibit a more selective pattern, showing strong involvement in disorders with cognitive impairment (ASD with ID, SCZ, NDD, VNR–), but show minimal or absent bias in ASD without ID. This selective profile motivated our focus on CGE/VIP populations in the functional experiments, which is further justified by the relative lack of in vivo data on interneuron contributions to cognitive deficits, compared to much more extensively studied excitatory circuits. We have clarified this rationale in the revised Results and Discussion.



Regarding the Reviewer's specific question about whether the electrophysiological distinction between CCKBC and ISI/VIPo interneurons is recapitulated in the human mutation-bias analysis: we performed an extensive cross-species analysis to address this point, but found that the mouse VIP subtype distinction does not cleanly translate to human transcriptomic data. CCK basket cells (CCKBCs) are defined in mouse as a discrete electrophysiological and transcriptomic subtype within the Sncg subclass (VIP–, CCK+, fast-spiking). However, no discrete CCKBC cluster exists in current human transcriptomic atlases. Cross-species classification of Sncg subtypes yields near-chance performance (AUROC = 0.50; Bakken et al., Nature 2021), the worst of any interneuron subclass, and CCK mRNA is expressed ubiquitously across human neurons rather than being interneuron-specific (Darmanis et al., PNAS 2015), precluding its use as a discriminative marker.

Despite this fundamental limitation, we attempted to impute human CCKBC identity using two independent approaches. First, we performed cross-species transcriptomic mapping by integrating mouse patch-seq data (5,764 cells, including 333 CCKBCs from the Sncg subclass) with the human Siletti et al. atlas (175,000 CGE interneurons) using both Harmony batch correction and scVI latent-space integration. Both methods identified the same set of human CGE clusters (279, 280, 281) as receiving the highest fraction of mouse CCKBCs (Harmony: 67–96% CCKBC fraction; scVI: 86–90%), with significant concordance between the two methods (Spearman rho = 0.46, P = 0.03). Second, we performed an electrophysiology signature transfer: we defined the mouse CCKBC electrophysiological fingerprint (higher ISI CV, higher input resistance, longer time constant) relative to other VIP interneurons within mouse, then scored human CGE cells from the Lee & Bhatt human patch-seq dataset for similarity to this signature. The ephys-based ranking of human clusters correlated positively with the transcriptomic CCKBC fraction (Spearman rho = 0.49), providing independent support from a second modality.

We then tested whether these imputed "human CCKBC" clusters (n = 5 clusters: 277, 278, 279, 280, 281; defined as Harmony CCKBC fraction > 50%) show differential mutation bias compared to the remaining 16 non-CCKBC CGE clusters. Contrary to the prediction from the mouse electrophysiology data, imputed CCKBCs showed *higher* — not lower — mutation bias than non-CCKBC CGE clusters for ASD (Mann-Whitney P = 0.015), ASD without ID (P = 0.04), and schizophrenia (P = 0.0002). For the 22q11.2 gene set, the direction was consistent with the mouse prediction (imputed CCKBCs slightly lower) but did not approach significance (P = 0.97), likely reflecting the low statistical power of the 22q11.2 analysis noted above. We note that these results should be interpreted cautiously: the imputed CCKBC classification is inherently noisy due to the poor cross-species conservation of the Sncg subclass, and the sample sizes are small (5 vs. 16 clusters).

Taken together, these analyses demonstrate that: (1) a direct one-to-one mapping between mouse CCKBC/ISI subtypes and human CGE clusters is not currently feasible due to evolutionary divergence in Sncg interneuron identity; (2) even when we impute human CCKBC clusters through convergent transcriptomic and electrophysiological evidence, the mouse subtype-specific bias pattern does not replicate; and (3) the broader CGE/VIP supercluster — rather than individual subtypes within it — remains the most robust and reproducible unit for cross-species integration of mutation-bias signals. These results reinforce our decision to focus the main analyses at the supercluster level and are consistent with the interpretation that rare-variant convergence on CGE interneurons reflects a lineage-level, rather than subtype-specific, vulnerability.



7. Splitting ASD by intellectual disability (ID). While the motivation for splitting ASD cases into "with ID" and "without ID" is reasonable, the labeling should be made explicit — i.e., "ASD with ID" and "ASD without ID" — to clarify that this categorization was imposed by the authors rather than inherent to the source datasets. To strengthen the cross-disorder comparisons, a similar treatment should be applied to schizophrenia, which is also associated with variation in premorbid IQ. If individual-level data are not available to stratify SCZ cases by cognitive performance, then at minimum the authors should remove DD/ID genes from the schizophrenia gene list. Otherwise, the contrast between "ASD with ID" and SCZ risks being confounded by differences in how cognitive impairment is represented in the gene sets.



We thank the reviewer for this thoughtful and important comment. We agree that the ASD subgroup labels need to be made explicit, and we now consistently use “ASD with ID” and “ASD without ID” throughout the manuscript. We also clarify in the Methods that this categorization is defined by our analysis rather than being inherent to the original datasets.



To directly address the reviewer’s concern that the signal in schizophrenia might be driven by overlapping DD/ID genes, we performed a comprehensive overlap and removal analysis. First, we quantified the overlap between our high-confidence SCZ gene set and established DD/ID (NDD) gene lists. We found that the overlap is remarkably small: among the SCZ genes included in our main analysis (61 genes), only 1 gene overlapped with the top 61 DD/ID genes, and only 9 genes overlapped when comparing against a broader list of 297 DD/ID genes (genome-wide significant genes). We then repeated the mutation-bias analyses for schizophrenia after removing these overlapping genes, and found that the CGE interneuron bias in SCZ remained robust and virtually unchanged following this removal. The overall cell-type bias profile was also highly stable (Spearman correlation > 0.92 with the original profile). This result refutes the hypothesis that the SCZ signal is merely driven by confounding DD/ID genes. Instead, it demonstrates that risk genes specific for schizophrenia, which is clinically associated with specific cognitive deficits distinct from severe ID, independently converge on CGE interneurons. This finding reinforces our broader conclusion that distinct genetic architectures (SCZ vs. ASD with ID vs. NDD) converge on CGE interneurons as a shared cellular substrate for cognitive dysfunction. We have included this sensitivity analysis in the Supplementary Information (Figure SXX) and clarified the interpretation in the text.





8. VIP function in the 22q11.2 mouse model

The decrease in correlation with movement is a strong phenotype that is also subtype specific. I wonder if this is due to difference in ACh receptors or generally decreased inputs. One major outstanding question in the manuscript is (as I mention above) if the changes seen in vivo are due to intrinsic changes in VIP cells (as predicted by computational part) or a difference in network integration. Patch clamp recordings from VIP neurons with input stimulation (to test pre-synaptic (paired pulses) v.s. postsynaptic function (EPSCs or similar) would help resolve this issue and thereby test one of the main claims. If this can be done then I think it would greatly enhance the manuscript but I would be OK with addressing this point clearly in the interpretation of the results.



Reply Plan: We thank the Reviewer for this comment. We agree that a key unresolved question is whether the observed in vivo alterations in VIP interneuron activity reflect intrinsic cellular dysfunction or changes in network integration and afferent input. To directly address this issue, we will include whole-cell patch-clamp recordings from VIP interneurons in the 22q11.2 deletion mouse model. In the revised manuscript, we will integrate these electrophysiological results with the in vivo calcium-imaging data to refine the mechanistic interpretation, clarifying the extent to which VIP interneuron dysfunction arises from cell-intrinsic alterations versus altered circuit integration.



Minor Concerns

Please provide full mutation bias estimates for all traits tested, not just selected highlights. Preferably provide all mutational bias results for all traits tested in a "long" format, with additional column describing the tested trait. Gene lists and mutation counts used for each trait should be published in supplementary tables.



We thank the reviewer for this helpful suggestion. In the revised manuscript, we now provide complete mutation-bias estimates for all tested traits, rather than only the highlighted examples shown in the main figures. These are included as Supplementary Tables X–Y, in a long-format structure in which each row corresponds to a (trait × cell type) pair. For each trait, we additionally supply the full gene list, gene-level mutation counts, and expected mutation burdens used to compute the mutation-bias statistics. These additions ensure full transparency and facilitate reproducibility of all analyses.



Clarify why the cell-type specificity score is capped at 2. How sensitive are results to this capping, and why not use previously established specificity metrics?





We thank the reviewer for raising this point. The specificity cap, set conservatively at 2x mean fold-enrichment, is a necessary guard against technical inflation, particularly in cell types with low sequencing depth.



In our atlas of 461 cell types, uncapped specificity values can reach up to $\sim$97x. These extreme scores are artifacts of sparse sampling, not genuine expression enrichment, and are driven primarily by small non-neuronal populations (e.g., vascular, fibroblast, microglia) with low total UMI. Quantitatively, without the cap, $56.7%$ of the total specificity sum in non-neuronal cell types comes from genes exceeding the threshold, compared to only $15.4%$ for neuronal types (Spearman $\rho$ between UMI depth and fraction of genes exceeding cap $ = -0.497, p = 3.6e-30$).



To confirm the robustness of our choice, we recomputed SCZ mutation bias under four cap levels ($1\times, 2\times, 3\times$, and uncapped). A cap of $3\times$ produced rankings nearly identical to the default $2\times$, confirming its suitability. Conversely, a cap of $1\times$ over-clips and reduces discriminability, while no cap yields spurious extreme bias scores for poorly-sampled cell types. These analyses are presented in the new Supplementary Figure Sx.



Finally, regarding previously established metrics (e.g., $\tau$, pSI), our fold-enrichment approach was chosen because it is directly interpretable in terms of expression magnitude and is compatible with the weighted-sum bias calculations required by our analytical framework. Rank-based metrics, in contrast, discard the quantitative information needed for this weighted analysis.







We appreciate the reviewer’s insightful question regarding the capping of the cell-type specificity score. Our specificity metric is otherwise identical to previously established measures (e.g., normalized expression specificity across all cell types; Methods Equation 2-3), with the sole exception that we introduce an upper cap of 2.0. This cap is required because of the interaction between specificity and the mutation-bias metric (Methods Equation 6), where specificity is used as a multiplicative weight, rather than as a ranking metric to define cell type specific expression gene sets.



In single-cell datasets, genes with very low absolute expression, particularly in cell types with limited sequencing depth, or low cell volume—can exhibit artificially inflated specificity scores due to noise or dropout artifacts. Although these genes are biologically lowly expressed, their specificity values can become extremely large (≫ 2), causing them to exert disproportionate influence on the mutation-bias statistic despite low or noisy expression.

[insert a figure show cell type expression specify, distribution cross cell types; and a figure show correlation between cell total UMI vs largest specificity]



Previous studies using this specificity metric typically did not use specificity directly as a continuous multiplier, but instead relied on rank-based thresholds (e.g., “top 10% most specific genes”). Such rank-based methods are robust to extreme specificity values. However, because our mutation-bias framework multiplies expression specificity with mutation counts, extreme specificity outliers can dominate the score unless appropriately bounded.



The specificity metric is normalized so that the average gene-specificity across all cell types equals 1. A cap at 2.0 therefore represents a symmetric and biologically interpretable upper bound: allowing weights up to 2× the average specificity, while preventing unrealistic inflation from noisy low-expression genes. Empirically, we observed that highly expressed, biologically meaningful marker genes almost never exceed specificity = 2 (Reviewer Figure XXX), while inflated values are almost exclusively associated with low-expression genes or clusters with limited sequencing depth. Thus, the cap selectively removes artifacts without affecting real signal.



To address the reviewer’s question about robustness, we now include additional supplementary analyses (Supplementary Figure X) showing: (i) the distribution of specificity scores before and after capping; (ii) the effect of alternative caps (e.g., 1.5, 2.5, 3.0) on mutation-bias estimates; (iii) comparisons of mutation-bias results using uncapped specificity, rank-based specificity, and the capped metric. Across these analyses, the qualitative mutation-bias patterns are preserved, while the cap reduces variance and prevents outlier-driven instability. We describe these findings in the revised Methods and Supplementary Note.









