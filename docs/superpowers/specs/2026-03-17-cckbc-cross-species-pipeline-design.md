# Cross-Species Multi-Modal CCKBC/ISI VIP Pipeline Design

**Date:** 2026-03-17
**Status:** Draft
**Context:** Reviewer 2 Q6 / Reviewer 3 Q6 — Robustify cross-species CCKBC vs ISI VIP 22q bias analysis

## Problem

The current analysis directly maps mouse patch-seq cells to human single-cell clusters (cross-species + cross-modality in one jump). While the results are reasonable (VIP- CCKBC weakly biased, VIP+ CCKBC biased like ISI VIP), the methodology is vulnerable to reviewer criticism because it relies on a single mapping path.

Four datasets are available but underutilized:

| Dataset | Species | Modality | Source | Cells |
|---------|---------|----------|--------|-------|
| Mouse single-cell | Mouse | snRNA-seq (10x) | Allen WMB-10Xv3 | ~456K isocortex |
| Mouse patch-seq | Mouse | Smart-seq2 + ephys | DANDI 000008 (M1) | 1,330 |
| Human single-cell | Human | snRNA-seq | Siletti et al. 2023 | 461 clusters |
| Human patch-seq | Human | Smart-seq2 + ephys | DANDI 000636 (Lee & Dalley) | 711 sessions |

Current edges: mouse patch-seq -> human SC (direct), human patch-seq -> human SC (scANVI), mouse patch-seq -> mouse SC (scVI, noisy). Missing: mouse SC -> human SC (cluster-level bridge).

## Goal

Build a robust multi-modal pipeline that:
1. Leverages all 4 datasets as stepping stones rather than one cross-species+cross-modality jump
2. Validates mapping consistency via direct vs indirect path concordance
3. Uses electrophysiology (from same EphysSumStats pipeline on both species) as independent functional validation
4. Recalculates 22q bias on CCKBC vs ISI groupings supported by multi-modal evidence
5. If the refined pipeline shows significant ISI > CCKBC 22q bias, that is a bonus but not required

## Prerequisites

### Conda Environment
All code runs in the `gencic` conda environment. Additional packages to install:
- `harmonypy` — Harmony batch correction (Module B primary method)
- `neuroCombat` — ComBat batch correction (Module D Tier 2 sensitivity analysis, optional)

### Packages Already Available
- `scanpy`, `anndata` — single-cell analysis (already in gencic)
- `scipy`, `sklearn` — statistics and ML (already in gencic)

### Packages NOT Required
- `scvi-tools` — existing scVI results from atlas_matching are used as a comparator baseline; no new scVI training is needed unless results are missing. If scVI hyperparameter sweep is desired, install separately.
- `MetaNeighbor` (R package) — NOT used. Module A Pass 2 implements the MetaNeighbor algorithm natively in Python using scipy/sklearn (~50 lines: Spearman cell-cell correlation -> k-NN graph -> AUROC per cluster pair).

### Required Existing Outputs
- Ortholog mapping: /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/cross_species/orthologs/ortholog_mapping.csv (22,563 rows)
- Existing scVI mapping results: /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/scvi_mapped/
- Human patch-seq -> Siletti scANVI results: cge_subtype/results/human_scvi_mapping_results.csv
- Mouse patch-seq -> Human SC (direct): cge_subtype/results/patchseq_mapping_results.csv
- EphysSumStats processed data: /mnt/data0/DANDI/Processed/000008/ and /mnt/data0/DANDI/Processed/000636/

### Estimated Resources
- Module A: ~2-4 hours (pseudobulk computation RAM-bound, MetaNeighbor fast on subsampled data)
- Module B: ~8-12 hours (hyperparameter sweep, parallelizable on 10 cores)
- Module C: <1 hour (table joins and statistics)
- Module D: ~2-3 hours (aggregation + statistical tests)
- Module E: <1 hour (bias calculation + plotting)
- Disk: ~5 GB for intermediate results (pseudobulk matrices, subsampled references)

---

## Pipeline Architecture

```
                    Module A (cluster-level bridge)
Mouse SC (Allen WMB) ──────────────────────> Human SC (Siletti 461)
   ^                                              ^
   | Module B (Harmony,                           | (existing scANVI,
   |  improved + hyperparameter sweep)            |  already done)
   |                                              |
Mouse Patch-seq ──────────────────────────> Human SC (Siletti 461)
                  (existing Harmony/scVI,          ^
                   direct path)                    |
                                            Human Patch-seq
                                            (Lee & Dalley)

Module C: Compare direct vs indirect paths for mouse patch-seq
Module D: Ephys convergence at shared human clusters
Module E: 22q bias on CCKBC vs ISI groupings
```

Existing work kept as-is:
- Human patch-seq -> Human SC (scANVI, 96% purity)
- Mouse patch-seq -> Human SC (direct Harmony/scVI)
- 22q bias calculation machinery (HumanCT_AvgZ_Weighted())

---

## Module A: Mouse SC -> Human SC (Cluster-Level Bridge)

### Goal
Establish which mouse clusters correspond to which human clusters (461) at the cluster level (no individual cell-to-cell mapping required). Mouse cluster identity is defined by the `cluster_alias` column in the Allen WMB-10Xv3 metadata; the exact count (~300-400 cortical clusters) should be confirmed during data exploration after filtering to the cortical region subset defined in cross_species_config.yaml.

### Two-Pass Strategy

**Pass 1 — Spearman RBH on pseudobulk centroids (resolve easy clusters):**
1. Compute pseudobulk (mean expression) per cluster for each species. Load one atlas at a time to control RAM (backed mode, chunk by cluster).
2. Convert mouse genes to human orthologs using existing ortholog table (cge_subtype/scripts/00_prepare_orthologs.py). Restrict to 1:1 orthologs (~15-16K genes).
3. Subset to top 3,000 HVGs by variance across cluster centroids in combined matrix.
4. Spearman correlation matrix (mouse ~360 x human ~461).
5. Assign reciprocal best hits (RBH). The correlation threshold for "resolved" pairs is determined empirically: compute the full RBH correlation distribution, inspect the histogram for a natural gap between high-confidence matches and ambiguous ones, and set the threshold at that gap. If no clear gap exists, use a data-driven cutoff such as the inflection point of the sorted RBH correlation curve or a permutation-based null (shuffle cluster labels, compute RBH correlations, set threshold at the 95th percentile of the null).
6. Expected: most major classes (Pvalb, Sst, L2/3 IT, L5 ET, etc.) resolve here.

**Pass 2 — MetaNeighbor-style AUROC on unresolved clusters (resolve hard clusters):**
1. Take clusters that failed RBH or had ambiguous matches (multiple near-ties).
2. Subsample cells from these clusters only (~100-200 cells/cluster).
3. Compute MetaNeighbor unsupervised AUROCs (Crow et al. 2018) using a Python-native implementation: (a) Spearman cell-cell correlation matrix on shared ortholog HVGs, (b) for each cell, rank neighbors by correlation, (c) compute AUROC for same-cluster vs different-cluster neighbors using scipy.stats or sklearn.metrics. This is ~50 lines of Python, no R dependency needed.
4. Expected: CGE subtypes (VIP, Sncg, Lamp5 boundary) land here.
5. Bakken et al. 2021 showed Sncg subclass has near-chance AUROC (~0.5) across species — this is a known limitation to report honestly.

### Output
- Full cross-species cluster correspondence table with method tag (Spearman-RBH vs MetaNeighbor) and confidence score per pair
- CGE-specific AUROC matrix
- Which mouse VIP/Sncg clusters map to which human CGE clusters (276-296)

### RAM Management
- Peak RAM during pseudobulk: load atlas in backed mode (sc.read_h5ad(..., backed='r')), compute pseudobulk in chunks by cluster, never hold full dense matrix
- MetaNeighbor pass: only subsampled cells from unresolved clusters (~4-8K cells for CGE), negligible RAM

### Key Data Paths
All paths are also documented in cge_subtype/configs/cross_species_config.yaml (single source of truth).
- Mouse SC: /mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-{1,2}-raw.h5ad
- Mouse SC metadata: /mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv
- Human SC: /mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_*.h5ad (per-supercluster files; see cross_species_config.yaml for full list)
- Ortholog mapping: /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/cross_species/orthologs/ortholog_mapping.csv
- Ortholog generation script: cge_subtype/scripts/00_prepare_orthologs.py

---

## Module B: Mouse Patch-seq -> Mouse SC (Improved)

### Goal
Improve the existing noisy mapping (~84% accuracy, neuron->glia and Glut/GABA cross-mapping errors) using Harmony as an alternative, evaluated against Cre-line ground truth.

### Steps
1. **Subsample reference** — Allen WMB-10Xv3, 100-200 cells/cluster (~50K total). Avoids RAM crashes.
2. **Harmony integration** — PCA (50 dims) + Harmony with batch key = technology (10x vs Smart-seq2).
3. **k-NN label transfer** — k=30 nearest neighbors in Harmony-corrected space, majority vote for cluster assignment, confidence = fraction agreeing.
4. **Neuronal-only filtering** — Second pass with non-neuronal cells removed (script 06 exists in atlas_matching). Compare neuron->glia rate.
5. **Evaluate against Cre-line ground truth:**
   - V1 (Gouwens): corresponding_AIT2.3.1_alias -> subclass
   - M1 (Tasic): RNA family + RNA type
   - Metrics: subclass accuracy, confusion matrix, neuron->glia rate
6. **Pick best method** — Highest subclass accuracy with neuron->glia < 5%.

### Hyperparameter Sweep
Harmony:
- n_pcs: [30, 50, 75]
- n_neighbors for k-NN: [15, 30, 50]
- theta (diversity penalty): [1.0, 2.0]
- HVG count: [2000, 3000, 5000]

scVI (comparator, optional — uses existing trained models from atlas_matching as baseline):
- Use pre-trained models from /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/models/ as comparators
- If new scVI runs are desired, requires installing scvi-tools in gencic env
- n_latent: [30, 50], n_hidden: [128, 256], Likelihood: [NB, Poisson] (8 combos)

Total: ~24 Harmony combinations (primary) + existing scVI baselines (comparator). Parallelizable, run overnight on 10 cores.

**Selection criterion:** Highest CGE subclass accuracy (Vip, Sncg, Lamp5) with neuron->glia rate < 5%.

### Output
- Mouse patch-seq -> mouse cluster assignment table with confidence scores
- Accuracy benchmarks (confusion matrices against Cre-line labels)
- Method comparison table (Harmony vs scVI vs scVI-neuronal-only, across hyperparameters)

### Key Data Paths
- Mouse patch-seq V1 counts: /home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_counts.csv
- Mouse patch-seq M1 counts: /home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_counts.csv.gz
- V1 metadata: /home/jw3514/Work/NeurSim/TransEphys/dat/expression/V1_patchseq_metadata.csv (ground truth: corresponding_AIT2.3.1_alias)
- M1 metadata: /home/jw3514/Work/NeurSim/hh_sbi/data/m1_patchseq_meta_data.tsv (ground truth: RNA family, RNA type)
- Mouse SC reference (subsampled): see Module A paths
- Existing scVI models: /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/models/
- Existing mapping results: /home/jw3514/Work/NeurSim/TransEphys/atlas_matching/results/scvi_mapped/

---

## Module C: Path Concordance Validation

### Goal
Show that direct (mouse patch-seq -> human SC) and indirect (mouse patch-seq -> mouse SC -> human SC) paths give consistent answers.

### Steps
1. **Build indirect path assignments:**
   - Mouse patch-seq -> mouse cluster (Module B)
   - Mouse cluster -> human cluster (Module A)
   - Chain = indirect human cluster assignment per mouse patch-seq cell

2. **Compare direct vs indirect** for each mouse patch-seq cell.

3. **Concordance metrics at 3 levels:**
   - Supercluster (MGE, CGE, excitatory, etc.): expect >90% agreement
   - Subclass (Pvalb, Sst, Vip, Sncg, Lamp5): expect moderate-high
   - Cluster (Siletti 276-296 for CGE): expect lower but informative
   - Cohen's kappa at each level (chance-corrected agreement)

4. **CCKBC focus:**
   - For M1 CCKBCs (~19 cells from Sncg subclass + Vip Sncg type in M1 patch-seq, per CCKBC definition in atlas_matching/docs/cckbc_definition.md): do both paths agree on which human CGE clusters?
   - Agreement = strong evidence for CCKBC correspondence
   - Disagreement = honest reporting, reinforces poor Sncg cross-species conservation

5. **Handling ambiguous mappings in the chain:**
   - When Module A produces many-to-many correspondences (likely for CGE subtypes), use the top-ranked human cluster for the primary concordance analysis
   - Additionally report a "relaxed" concordance: indirect assignment agrees with direct if the direct cluster is among the top-3 Module A candidates for that mouse cluster
   - Flag clusters with ambiguous Module A mappings and report coverage (fraction of cells with unambiguous indirect assignments) alongside concordance

6. **Visualization:**
   - Sankey/alluvial: mouse subclass -> human cluster (direct vs indirect side by side)
   - Concordance heatmap at each resolution
   - Scatter: direct confidence vs indirect confidence, colored by agreement

### Output
- Concordance table with agreement rates at 3 levels
- Cohen's kappa statistics
- CCKBC-focused concordance summary
- Figures for reviewer response

---

## Module D: Cross-Species Ephys Convergence

### Goal
Validate that mouse and human patch-seq cells mapping to the same human cluster share electrophysiological features — functional evidence independent of transcriptomics.

### Key Advantage: Same Extraction Pipeline
Both mouse (DANDI 000008) and human (DANDI 000636) were processed by the same in-house EphysSumStats pipeline with identical feature extraction logic, producing ~77 columns per cell in analysis.csv. This eliminates the feature-naming mismatch that limited the previous analysis to 5 manually mapped features.

- Mouse: /mnt/data0/DANDI/Processed/000008/ (1,330 cells)
- Human: /mnt/data0/DANDI/Processed/000636/ (711 sessions)

### Ephys Harmonization (3 tiers)

**Tier 1 — Standard (what the field does, primary analysis):**
- Aggregate per-cell summaries from analysis.csv (select representative spiking sweep: target 5-60 Hz, fallback rheobase)
- Match columns directly (same EphysSumStats names). Drop ISI bins beyond bin 11 from human data.
- Remove QC/metadata columns, keep ~30-40 biologically meaningful features
- Log-transform skewed features (firing rates, ISI values) before z-scoring
- Within-species z-scoring (standard in Lee et al. 2023, Bakken et al. 2021, Cadwell/Scala 2023)
- Rank-based non-parametric tests (Mann-Whitney, Kruskal-Wallis)

**Tier 2 — Sensitivity analysis:**
- Apply ComBat with species as batch, cell-type label as covariate (neuroCombat Python package)
- Compare results with Tier 1 — if conclusions stable, result is robust to harmonization method

**Tier 3 — Exploratory (if time permits):**
- Sparse PCA on z-scored features (following Gouwens et al. 2020) to ~20-30 components
- Compare in sPCA space for multivariate ephys signatures

### Analysis Steps
1. **Identify shared clusters** — Human CGE clusters (276-296) with both mouse and human patch-seq cells mapped to them.
2. **Per-cluster comparison** — For each shared cluster, compare z-scored ephys profiles (mouse vs human). Metrics: Spearman correlation of mean profiles, Euclidean distance, Mann-Whitney per feature.
3. **Global permutation test** — Across all shared clusters: are same-cluster cross-species pairs more similar than different-cluster pairs? Shuffle cluster labels, recompute, get p-value.
4. **CCKBC-specific** — Do mouse CCKBC cells mapping to putative CCKBC clusters show narrow AP / fast-spiking? Do VIP+ ISI clusters show irregular-spiking?

### Output
- Table of shared clusters with cell counts per species
- Per-cluster ephys profile comparison (correlation + distance)
- Global permutation test result
- CCKBC vs ISI ephys profile comparison
- Visualization: paired ephys radar plots or heatmap

### Realistic Expectations
Recording conditions differ (M1 cortex mouse vs human cortex, different labs). A trend-level result is still useful. Flat result is honest and publishable (consistent with Bakken et al.'s Sncg conservation findings).

### Key References
- Lee et al. 2023 Science (human GABA patch-seq, no formal batch correction)
- Bakken et al. 2021 Nature (BICCN cross-species, qualitative ephys comparison)
- Cadwell/Scala et al. 2023 Science (L1 interneurons, Kruskal-Wallis + Dunn's)
- Gouwens et al. 2020 Cell (IPFX -> sPCA pipeline)
- Schwider & Ramezani 2026 arXiv (transfer learning for cross-species ephys)
- Existing cross-species comparison: TransEphys/jaxley_fit/notebooks/cross_species_ephys_comparison.py

---

## Module E: Updated 22q Bias on Refined Groupings

### Goal
Recalculate 22q bias using CCKBC/ISI groupings supported by multi-modal evidence from Modules A-D.

### CCKBC Classification (Multi-Modal Confidence Score)
For each human CGE cluster (276-296), assign CCKBC confidence based on converging evidence:
- Module A: Does the cluster's best-matching mouse cluster belong to Sncg/CCKBC subclass?
- Module B->C: Do mouse CCKBC patch-seq cells map here via both direct and indirect paths?
- Module D: Does the cluster's ephys profile resemble mouse CCKBC (fast-spiking, narrow AP)?
- Existing: Marker expression (CCK+, SNCG+, VIP status)

Classification tiers:
- 3+ lines of evidence = "high-confidence CCKBC"
- 1-2 lines = "tentative CCKBC"
- 0 lines = "ISI VIP" (or other CGE)

### 22q Bias Calculation
Same machinery: HumanCT_AvgZ_Weighted() on Siletti specificity matrix.
4 gene sets: 22q_del, 22q_mouse, 22q_HighExp, 22q_DEG_d75.

### Group Comparisons
- **Primary:** High-confidence CCKBC vs ISI VIP (Mann-Whitney U)
- **Secondary:** 3-way split by VIP expression (VIP- CCKBC vs VIP+ CCKBC vs VIP+ ISI)
- **Correlation:** CCKBC confidence score vs 22q bias across all CGE clusters (Spearman)

### Sensitivity Analysis
- Classification based on transcriptomic evidence only vs ephys evidence only
- Across all 4 gene sets
- Bootstrap confidence intervals on the group difference

### Output
- Updated bias comparison table with multi-modal CCKBC classification
- Statistical tests (p-values, effect sizes)
- Sensitivity analysis results
- Figures: updated boxplot, scatter (confidence vs bias), evidence summary panel

### Expected Outcome
VIP expression split likely persists (VIP- CCKBC weak, VIP+ CCKBC strong). Now supported by robust multi-modal classification. Significant ISI vs CCKBC difference would be a bonus.

---

## Summary

| Module | Input | Output | Status |
|--------|-------|--------|--------|
| A | Mouse SC + Human SC atlases | Cluster correspondence table | New |
| B | Mouse patch-seq + Mouse SC ref | Improved mapping + benchmarks | Improved |
| C | Modules A + B + existing direct mapping | Concordance metrics | New |
| D | Mouse + human patch-seq ephys (EphysSumStats) | Ephys convergence validation | New (5 -> ~40 features) |
| E | Modules A-D + 22q gene sets | Updated bias comparison | Updated |

### Dependencies
- Module A: independent (can start immediately)
- Module B: independent (can start immediately, parallelize with A)
- Module C: depends on A + B
- Module D: depends on existing mappings (can start in parallel with A/B for data prep)
- Module E: depends on A + B + C + D (final step)

### Key Constraints
- RAM: Never load full atlas into memory. Use backed mode, subsample, pseudobulk.
- Parallelization: Default n_processes=10 to avoid hogging resources.
- The analysis is supplementary validation. If individual modules produce null/negative results, report honestly — this strengthens the paper's credibility.
