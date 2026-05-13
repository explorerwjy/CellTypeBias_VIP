# CGE Subtype Analysis — Data Sources

## Mouse M1 Patch-seq

- **Raw NWB files:** `/mnt/data0/DANDI/000008/`
- **Processed ephys (EphysSumStats):** `/mnt/data0/DANDI/Processed/000008/`
- **Original source repo (with cell typing):** `/home/jw3514/Work/NeurSim/mini-atlas/`
- **M1 ephys features (custom extraction):** `/home/jw3514/Work/NeurSim/TransEphys/dat/expression/M1_patchseq_ephys_features.csv` (1,033 cells)
- **Cell type reference:** Table 1 from Scala et al. 2021 (Nature 598:144-150), saved at `docs/41586_2020_2907_Tab1_ESM.jpg`
  - CCKBC types are in the **Sncg subclass** (e.g., Sncg Gaba types in Allen WMB)
  - Also some VIP subclass cells express CCK (Vip/CCK basket cells)

### CCKBC identification
- Sncg subclass from Allen taxonomy = transcriptomically defined CCKBCs
- `harmony_patchseq_mapping_results.csv` column `is_cckbc` flags these cells
- 333 CCKBCs total in the 5,764 M1+V1 patchseq cells (Sncg + Vip/Sncg overlap)

## Mouse V1 Patch-seq

- **Raw NWB files:** `/mnt/data0/DANDI/000020/`

## Human GABAergic Patch-seq (Lee & Dalley / Berg et al.)

- **Raw NWB files:** `/mnt/data0/DANDI/000636/`
- **Processed ephys (EphysSumStats):** `/mnt/data0/DANDI/Processed/000636/`
- **Ground truth metadata:** `/home/jw3514/Work/NeurSim/human_patchseq_gaba/data/`
  - `LeeDalley_manuscript_metadata.csv` — 778 cells with TRUE transcriptomic labels
  - `LeeDalley_ephys_fx.csv` — 704 cells with ephys features (Allen IPFX extraction)
  - Key columns: `Revised_subclass_label` (PVALB/SST/VIP/LAMP5), `Transcriptomic_type` (44 types)
  - **Brain region: Middle Temporal Gyrus (MTG)**, NOT motor cortex

### Human VIP cells
- 122 VIP cells (by `Revised_subclass_label`), 111 with ephys
- 21 transcriptomic types including VIP SERPINF1 (possible CCKBC homolog) and various ISI types
- NO explicit "CCKBC" label in human data — must be inferred

### WARNING
- `cge_subtype/results/human_scvi_mapping_results.csv` is from scVI inference, NOT ground truth
  - Contains obvious errors (132 "Microglia" from interneuron patchseq data)
  - Do NOT use as cell type labels; use `LeeDalley_manuscript_metadata.csv` instead
