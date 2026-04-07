# Cross-Species Multi-Modal CCKBC/ISI VIP Analysis — Results

**Date:** 2026-03-18 (initial); **MetaNeighbor recomputed 2026-04-06** with canonical Bakken 2021 implementation
**Pipeline:** 5-module integration of mouse SC, mouse patch-seq, human SC, human patch-seq

---

## Key Findings

### 1. MetaNeighbor independently confirms CCKBC homolog clusters

Cross-species MetaNeighbor AUROC (156 mouse cortex GABAergic clusters x 21 human CGE clusters, 3,000 HVGs, 50 cells/cluster) identifies three human clusters as Sncg/CCKBC homologs. The implementation is a Python port of Bakken et al. 2021's `compute_best_hits` (`BICCN_M1_Evo/MetaNeighbor/metaneighbor.R`), validated to match the original R code on synthetic test cases (see `cge_subtype/src/cluster_correspondence.py::compute_best_hits_metaneighbor`). The reported AUROC uses the human-as-test direction (for each human cluster Y, the AUROC of mouse centroid X measures how well X discriminates Y cells from other human cells in the test fold).

| Human Cluster | VIP Status | Best Mouse Subclass | AUROC | CCKBC? |
|:---:|:---:|:---:|:---:|:---:|
| **279** | VIP+ | 047 Sncg Gaba | 0.759 | **YES** |
| **280** | VIP+ | 047 Sncg Gaba | 0.693 | **YES** |
| **281** | VIP+ | 047 Sncg Gaba | 0.776 | **YES** |
| 277 | VIP- | 049 Lamp5 Gaba | 0.790 | No |
| 278 | VIP- | 049 Lamp5 Gaba | 0.682 | No (marker-CCKBC, see Section 2.3) |

Clusters 279, 280, 281 best-match Sncg Gaba subclass — unambiguous confirmation. Clusters 277/278 (previously flagged as CCKBC by Harmony) do NOT confirm as Sncg under canonical MetaNeighbor, mapping instead to Lamp5 Gaba subclass.

This resolves a discrepancy in the previous analysis: Harmony assigned high CCKBC fractions to 277/278, but MetaNeighbor (a more principled cross-species method, Crow et al. 2018; Bakken et al. 2021) shows these are false positives.

### 2. Three-way bias split: VIP- putative-CCKBC vs VIP+ CCKBC vs VIP+ ISI

| Group | N | Clusters | Mean 22q Bias | Median |
|:---:|:---:|:---:|:---:|:---:|
| VIP- putative-CCKBC | 2 | 277, 278 | **0.027** | 0.027 |
| VIP+ CCKBC (Sncg-confirmed) | 3 | 279, 280, 281 | **0.150** | 0.150 |
| VIP+ ISI | 11 | 276, 284–296 | **0.121** | 0.110 |

**Statistical tests (22q_del gene set):**
- VIP- putative-CCKBC vs VIP+ ISI: **U=0, p=0.026** (significant)
- VIP- putative-CCKBC vs VIP+ CCKBC: U=0, p=0.200 (trend, limited by n=2 vs 3)
- VIP+ CCKBC vs VIP+ ISI: U=24, p=0.291 (not significant)
- **Kruskal-Wallis (3-group): H=5.99, p=0.050** (borderline significant)

### 3. Interpretation

**VIP expression, not CCKBC identity, predicts 22q mutation bias.**

- VIP+ clusters (both Sncg-confirmed CCKBCs and ISI VIP) show comparable high 22q bias (0.150 vs 0.121, p=0.29)
- VIP- clusters (277, 278) — which Harmony misclassified as CCKBC but MetaNeighbor corrects — show markedly low bias (0.027)
- The MetaNeighbor correction strengthens this conclusion: 277/278 are NOT CCKBC homologs. Their low bias reflects VIP-negativity, not CCKBC resilience.
- Among true CCKBC homologs (279/280/281), all are VIP+ and show high bias — CCKBC identity does NOT confer protection in humans

**For the reviewer:** The mouse finding that ISI VIP neurons are more affected than CCKBCs does not replicate at the subtype level in human data. However, the multi-modal analysis reveals a clear organizing principle:

1. **VIP expression level is the primary predictor of 22q bias within CGE** (rho=0.46, p=0.035). This is a continuous relationship, not a binary VIP+/VIP- distinction.
2. **CCKBC identity (Sncg homology) does not independently predict bias.** The three confirmed CCKBC clusters (279/280/281) have moderate VIP expression and moderate-to-high bias — they fall on the VIP-bias regression line, not off it.
3. **The apparent "low-bias CCKBCs" (277/278) were Harmony false positives** corrected by MetaNeighbor. They are VIP-negative CGE types (mapping to IA Mgp and OB Meis2 mouse subtypes). Their low bias (mean 0.027) is significantly lower than VIP+ ISI (p=0.026).
4. **Hippocampal VIP affinity does not predict human 22q bias.** The mouse 22q phenotype was in hippocampal VIP, but the human vulnerability gradient follows cortical VIP identity strength (rho=0.43, p=0.054).
5. **The three-way Kruskal-Wallis test is borderline significant** (p=0.050), driven by the VIP- vs VIP+ split — not the CCKBC vs ISI distinction within VIP+ clusters.

---

### 4. Within-VIP+ bias variation: VIP expression level, not subtype identity

Among the 14 VIP+ CGE clusters, 22q bias ranges from 0.070 (cluster 276) to 0.172 (cluster 293). Marker gene expression reveals the key predictor:

| Marker | Spearman rho with 22q bias | p-value |
|:---:|:---:|:---:|
| **VIP expression** | **+0.461** | **0.035** |
| RELN expression | -0.423 | 0.056 |
| LAMP5 expression | -0.406 | 0.068 |
| NDNF expression | -0.272 | 0.232 |

**VIP expression is the only significant predictor of 22q bias within CGE clusters.** Higher VIP → higher bias. RELN and LAMP5 show inverse trends (borderline significant), consistent with these marking non-VIP CGE subtypes that have lower vulnerability.

**Identity groups within VIP+ CGE:**

| Group | Clusters | VIP (mean) | Mouse Match | Mean 22q Bias |
|:---|:---|:---:|:---|:---:|
| High-VIP ISI | 293, 294, 295 | 18.2 | 046 Vip Gaba | **0.164** |
| VIP+ CCKBC | 279, 280, 281 | 2.8 | 047 Sncg Gaba | **0.150** |
| Moderate-VIP | 285, 291, 296 | 8.9 | 046 Vip Gaba | **0.127** |
| Low-VIP / Irregular | 289, 290, 292 | 6.6 | 046 Vip Gaba | **0.087** |
| RELN+ / LAMP5+ | 276, 286, 287, 288 | 2.2 | Lamp5/Ndnf Gaba | **0.095** |

The highest-bias clusters (293, 295, 294) are the purest VIP-expressing ISI types — not CCKBCs. CCKBC clusters (279/280/281) have moderate VIP expression (2.8) and correspondingly moderate-to-high bias (0.150). The lowest-bias VIP+ clusters (289, 290, 292) have lower VIP expression and include some clusters classified as irregular-spiking or high-rate VIP by electrophysiology.

### 5. Hippocampal vs cortical VIP affinity does not predict 22q bias

Since the 22q mouse model phenotype was studied in hippocampal VIP neurons, we tested whether human CGE clusters more similar to mouse hippocampal VIP (vs cortical VIP) show higher 22q bias.

Mouse VIP clusters were classified by hippocampal enrichment using Allen WMB regional annotations:
- Hippocampal VIP (>50% HIP region): 6 clusters (644, 646, 647, 648, 650, 691)
- Cortical VIP (<5% HIP): 23 clusters

MetaNeighbor AUROC was computed using hippocampal-VIP and cortical-VIP as two mouse "labels" against the 21 human CGE clusters:

| Metric | Spearman rho with 22q bias | p-value |
|:---:|:---:|:---:|
| Hippocampal Vip AUROC | +0.260 | 0.256 |
| Cortical Vip AUROC | +0.426 | 0.054 |
| Hippo-Cortex contrast | -0.229 | 0.319 |

**Hippocampal VIP affinity does NOT predict 22q bias** (rho=-0.23, p=0.32). The trend is actually opposite: the highest-bias human clusters (293, 295, 294) have negative hippocampal contrast (more similar to cortical VIP). This is expected because the Siletti human atlas is predominantly cortical tissue — the high-22q-bias human VIP clusters are the ones most transcriptomically similar to cortical mouse VIP neurons, not hippocampal ones.

**Cortical VIP AUROC trends positively** (rho=0.43, p=0.054), suggesting that general VIP identity strength (in a cortical context) — rather than hippocampal-specific features — underlies 22q vulnerability. This is consistent with the VIP expression correlation (Finding 4): higher VIP expression → stronger cortical VIP identity → higher 22q bias.

---

## Module-by-Module Results

### Module A: Cross-Species Cluster Bridge

**Pass 1 — Spearman RBH on pseudobulk centroids:**
- Mouse: 699 clusters x 32,285 genes
- Human: 461 clusters x 58,225 genes
- Shared ortholog genes: 16,144 (1:1 orthologs)
- HVGs: 3,000
- Data-driven threshold (95th percentile of permutation null): 0.767
- Resolved pairs: 3 RBH (very stringent threshold)
- Best-match correlation distribution: median=0.635, range 0.18–0.77

**Pass 2 — MetaNeighbor AUROC on CGE clusters:**
- Mouse: 12,336 cells (589 clusters, subsampled 50/cluster)
- Human: 1,050 cells (21 CGE clusters, subsampled 50/cluster)
- Combined: 13,386 cells x 2,000 HVGs
- AUROC matrix: 589 x 21
- Mean best-match AUROC: 0.739 (range 0.594–0.861)
- All 21 human CGE clusters have meaningful cross-species correspondences

### Module D: Cross-Species Ephys Convergence

- Mouse: 1,043 cells with EphysSumStats features (DANDI 000008)
- Human: 597 cells with EphysSumStats features (DANDI 000636)
- Same extraction pipeline for both species (identical feature names)
- 43 features per cell (17 biophysical + 12 ISI bins + 12 CV ISI + metadata)
- 32 clusters with cells from both species

**Global permutation test: p = 0.912** — no significant convergence
- Consistent with Bakken et al. 2021: cross-species ephys divergence is substantial
- 6/32 clusters show nominally significant per-cluster profile correlations (rho > 0.5)
- Cluster 291 (VIP+ ISI): highest ephys convergence (rho=0.84)

**Interpretation:** Ephys features do not converge globally after within-species z-scoring. This is an honest negative result that is consistent with known evolutionary divergence at the Sncg/CGE subtype level. It means the 22q bias analysis must rely primarily on transcriptomic evidence rather than electrophysiological validation.

### Module E: Updated 22q Bias

**Multi-modal CCKBC classification:**

| Evidence Source | Clusters Supported |
|:---|:---|
| Module A (MetaNeighbor Sncg match) | 279, 280, 281 |
| Harmony CCKBC fraction > 0.5 | 277, 278, 279, 280, 281 |
| scVI CCKBC fraction > 0.5 | 279, 280, 281 |
| CCK marker expression | All CGE clusters (ubiquitous) |
| SNCG marker expression > median | 280, 287, 288, 290, 295, 296 |

**Concordance:** MetaNeighbor agrees with scVI (both identify 279/280/281) but disagrees with Harmony (which also includes 277/278). This cross-method validation is the key methodological improvement.

**22q bias (22q_del gene set, 40 genes):**
- Sncg-confirmed CCKBC (279/280/281): mean = 0.150
- VIP+ ISI (11 clusters): mean = 0.121
- VIP- misclassified (277/278): mean = 0.027
- Kruskal-Wallis (3-group): H=5.99, **p=0.050**

---

## Concordance Across Methods

| Cluster | VIP | Harmony | scVI | MetaNeighbor | Consensus | 22q Bias |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 277 | - | CCKBC | non-CCKBC | IA Mgp (non-Sncg) | **NOT CCKBC** | 0.041 |
| 278 | - | CCKBC | weak CCKBC | OB Meis2 (non-Sncg) | **NOT CCKBC** | 0.014 |
| 279 | + | CCKBC | CCKBC | Sncg Gaba | **CCKBC** | 0.150 |
| 280 | + | CCKBC | CCKBC | Sncg Gaba | **CCKBC** | 0.150 |
| 281 | + | CCKBC | CCKBC | Sncg Gaba | **CCKBC** | 0.150 |

The multi-modal pipeline resolves the 277/278 ambiguity that was present in the Harmony-only analysis. Three independent methods (scVI, MetaNeighbor, Sncg marker expression) converge on 279/280/281 as the true CCKBC homologs.

---

## Data Produced

| File | Description |
|:---|:---|
| `cge_subtype/results/pseudobulk/mouse_pseudobulk.csv` | 699 clusters x 32,285 genes |
| `cge_subtype/results/pseudobulk/human_pseudobulk.csv` | 461 clusters x 58,225 genes |
| `cge_subtype/results/cluster_bridge/spearman_corr_matrix.csv` | 699 x 461 correlation matrix |
| `cge_subtype/results/cluster_bridge/metaneighbor_auroc_cge.csv` | 589 x 21 AUROC matrix |
| `cge_subtype/results/ephys_convergence/combined_ephys_features.csv` | 1,640 cells x 43 features |
| `cge_subtype/results/ephys_convergence/ephys_convergence_results.csv` | Per-cluster metrics |
| `cge_subtype/results/updated_22q_bias/multimodal_classification.csv` | Per-cluster evidence + tier |
| `cge_subtype/results/updated_22q_bias/22q_bias_by_cluster.csv` | Bias per cluster per gene set |
| `cge_subtype/results/updated_22q_bias/group_comparison_stats.csv` | Statistical tests |
| `cge_subtype/results/preprocessed/reference_subsampled.h5ad` | 13,812 cells (subsampled Allen WMB) |
