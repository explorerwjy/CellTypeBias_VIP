"""
Phase 1 of NegCtrl β-vs-uniform sensitivity analysis.

Builds 17 dev gene-weight files in dat/GeneWeights/_dev_betaweighted/:
  - 8 NegCtrl_<trait>_beta.gw.csv : same Entrez IDs as production NegCtrl files,
                                    weights replaced with abs(BETA Burden) from raw genebass.
  - 9 <set>_uniform.gw.csv : same Entrez IDs as production disorder/cognition files,
                             weights replaced with 1.0.

Format: 2-column CSV (Entrez_ID, weight), no header.

Plan: docs/plans/2026-05-13-negctrl-beta-vs-uniform.md
"""

import sys
from pathlib import Path
import pandas as pd
import yaml

PROJ_DIR = Path("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP")
GW_DIR = PROJ_DIR / "dat" / "GeneWeights"
DEV_DIR = GW_DIR / "_dev_betaweighted"
DEV_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import LoadGeneINFO

HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

with open(PROJ_DIR / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Part A: NegCtrl β-weighted files (8)
# ---------------------------------------------------------------------------
CTRL_DIR = PROJ_DIR / "dat" / "CTRL"
NEGCTRL_MAP = {
    "NegCtrl_HDL":       CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30760-both_sexes--irnt_2025_11_25_15_49_42.csv",
    "NegCtrl_Alanine":   CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30620-both_sexes--irnt_2025_12_16_20_18_07.csv",
    "NegCtrl_HbA1c":     CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30750-both_sexes--irnt_2025_11_25_16_32_28.csv",
    "NegCtrl_T2D":       CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-T2D_custom-both_sexes--custom_2025_11_25_16_14_20.csv",
    "NegCtrl_IBD":       CTRL_DIR / "gene-burden-results-exomes_pLoF_categorical-IBD_custom2-both_sexes--custom_2025_11_25_16_13_49.csv",
    "NegCtrl_RBC":       CTRL_DIR / "gene-burden-results-exomes_pLoF_continuous-30070-both_sexes--irnt_2025_12_16_20_41_11.csv",
    "NegCtrl_Parkinson": CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131022-both_sexes--_2025_11_25_16_56_06.csv",
    "NegCtrl_Alzheimer": CTRL_DIR / "gene-burden-results-exomes_pLoF_icd_first_occurrence-131036-both_sexes--_2025_11_25_17_00_04.csv",
}

print("=" * 72)
print("Part A: NegCtrl β-weighted files")
print("=" * 72)
print(f"{'name':<22} {'n_prod':>7} {'n_matched':>10} {'β_min':>10} {'β_max':>10}")
print("-" * 72)

for setname, raw_path in NEGCTRL_MAP.items():
    prod_path = Path(cfg["gene_sets"][setname]["geneweights"])
    prod_df = pd.read_csv(prod_path, header=None, names=["Entrez_ID", "weight"])
    n_prod = len(prod_df)

    raw = pd.read_csv(raw_path)
    pcol = next(c for c in raw.columns if "P" in c and "SKATO" in c)
    bcol = next(c for c in raw.columns if "BETA" in c and "Burden" in c)

    raw["EntrezID"] = raw["Gene Name"].map(
        lambda x: int(GeneSymbol2Entrez[x]) if x in GeneSymbol2Entrez and pd.notnull(GeneSymbol2Entrez.get(x)) else None
    )

    sym2beta = dict(zip(raw["EntrezID"].dropna().astype(int), raw[bcol]))
    out_rows = []
    missing = 0
    for eid in prod_df["Entrez_ID"].astype(int):
        beta = sym2beta.get(eid)
        if beta is None or pd.isna(beta):
            missing += 1
            continue
        out_rows.append((eid, abs(float(beta))))

    out_df = pd.DataFrame(out_rows, columns=["Entrez_ID", "weight"])
    out_path = DEV_DIR / f"{setname}_beta.gw.csv"
    out_df.to_csv(out_path, header=False, index=False)

    print(f"{setname:<22} {n_prod:>7} {len(out_df):>10} {out_df['weight'].min():>10.4f} {out_df['weight'].max():>10.4f}"
          + (f"  ({missing} missing)" if missing else ""))


# ---------------------------------------------------------------------------
# Part B: Uniform-weighted disorder + cognition files (9)
# ---------------------------------------------------------------------------
UNIFORM_SETS = ["ASD_All", "SCZ", "DDD_61",
                "UKBB_VNR_Pos", "UKBB_VNR_Neg",
                "UKBB_EDU_Pos", "UKBB_EDU_Neg",
                "UKBB_RT_Pos",  "UKBB_RT_Neg"]

print()
print("=" * 72)
print("Part B: Uniform-weighted disorder + cognition files")
print("=" * 72)
print(f"{'name':<22} {'n':>5} {'orig_w_min':>11} {'orig_w_max':>11}")
print("-" * 72)

for setname in UNIFORM_SETS:
    prod_path = Path(cfg["gene_sets"][setname]["geneweights"])
    prod_df = pd.read_csv(prod_path, header=None, names=["Entrez_ID", "weight"])
    out_df = prod_df.copy()
    out_df["weight"] = 1.0

    # Strip the UKBB_ prefix for filename compactness and append _uniform suffix
    short = setname.replace("UKBB_", "")
    out_path = DEV_DIR / f"{short}_uniform.gw.csv"
    out_df.to_csv(out_path, header=False, index=False)

    print(f"{setname:<22} {len(out_df):>5} {prod_df['weight'].min():>11.4f} {prod_df['weight'].max():>11.4f}")


# ---------------------------------------------------------------------------
# Verify
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print(f"All files written to {DEV_DIR}")
print("=" * 72)
written = sorted(DEV_DIR.glob("*.gw.csv"))
for p in written:
    df = pd.read_csv(p, header=None)
    print(f"  {p.name:<40} {len(df):>4} rows, weight range [{df[1].min():.3f}, {df[1].max():.3f}]")
print(f"Total: {len(written)} files")
