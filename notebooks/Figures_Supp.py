# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: gencic
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Supplementary Figures and Tables
#
# Figures S2, S5, S7-S9 and Tables S2-S11.
#
# **This notebook only loads pre-computed data and plots.**
# Run `Compute_BiasContrasts.ipynb` first to generate the required CSV files.
#
# Each multi-panel figure section saves individual panels, then assembles a
# composite PDF at the end.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os
import io
from contextlib import contextmanager

ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/"
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *
from matplotlib.image import imread

os.chdir(f"{ProjDIR}/notebooks/")

import matplotlib.font_manager as fm
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)
fm._load_fontmanager(try_read_cache=False)
plt.style.use('seaborn-v0_8-whitegrid')

FIG_DIR = os.path.join(ProjDIR, "results/figures/supp/")
os.makedirs(FIG_DIR, exist_ok=True)


# %%
# -- Panel capture utility --
@contextmanager
def save_panel(filepath, dpi=300):
    """Intercept plt.show() to save the figure before displaying inline."""
    _orig = plt.show
    plt.show = lambda *a, **kw: None
    try:
        yield
    finally:
        plt.show = _orig
        fig = plt.gcf()
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', transparent=True,
                    facecolor='none')
        plt.close(fig)
        from IPython.display import Image as IPImage, display as ipdisplay
        ipdisplay(IPImage(filename=filepath))


def assemble_panels(panel_paths, panel_labels, layout, figsize, out_path, dpi=300):
    """Assemble saved panel PNGs into a multi-panel composite figure."""
    nrows = max(r for r, _ in layout) + 1
    ncols = max(c.stop for _, c in layout)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='white')
    gs = fig.add_gridspec(nrows, ncols, hspace=0.08, wspace=0.08)

    for path, label, (row, col) in zip(panel_paths, panel_labels, layout):
        ax = fig.add_subplot(gs[row, col])
        img = imread(path)
        ax.imshow(img)
        ax.axis('off')
        ax.text(-0.02, 1.05, label, transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='bottom', ha='right')

    fig.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.show()
    print(f"Saved composite: {out_path}")


# %% [markdown]
# ## Load data

# %%
# Bias DataFrames
Bias_Save_Dir = os.path.join(ProjDIR, "results/main_results/random/Centering/")

ASD_All_Bias = pd.read_csv(Bias_Save_Dir + "ASD_All_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)
VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)

# %%
# Pre-computed contrasts
CONTRAST_DIR = os.path.join(ProjDIR, "results/main_results/contrasts/")

ASD_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_SCZ_contrast.csv", index_col=0)
ASD_wID_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_wID_vs_SCZ_contrast.csv", index_col=0)
HIQ_LIQ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_ASD_wID_contrast.csv", index_col=0)
ASD_DDD_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_DDD_contrast.csv", index_col=0)
SCZ_ASD_wID_Contrast = pd.read_csv(CONTRAST_DIR + "SCZ_vs_ASD_wID_contrast.csv", index_col=0)
VNR_Contrast = pd.read_csv(CONTRAST_DIR + "VNR_neg_vs_pos_contrast.csv", index_col=0)
EDU_Contrast = pd.read_csv(CONTRAST_DIR + "EDU_neg_vs_pos_contrast.csv", index_col=0)
all_contrasts_df = pd.read_csv(CONTRAST_DIR + "all_contrasts_fdr.csv")
brainspan_df = pd.read_csv(CONTRAST_DIR + "brainspan_expression.csv")

# Neuron-filtered contrasts
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
HIQ_LIQ_Contrast_Neurons = HIQ_LIQ_Contrast[HIQ_LIQ_Contrast.index.isin(Neurons)]
ASD_DDD_Contrast_Neurons = ASD_DDD_Contrast[ASD_DDD_Contrast.index.isin(Neurons)]
VNR_Contrast_Neurons = VNR_Contrast[VNR_Contrast.index.isin(Neurons)]
EDU_Contrast_Neurons = EDU_Contrast[EDU_Contrast.index.isin(Neurons)]

EffLabel = "EFFECT"

# %% [markdown]
# ---
# ## Figure S2 — BrainSpan developmental expression (single panel)

# %%
palette = sns.color_palette("Set2", 5)

Time = ['mean_2A', 'mean_2B', 'mean_3A', 'mean_3B', 'mean_4', 'mean_5',
        'mean_6', 'mean_7', 'mean_8', 'mean_9', 'mean_10', 'mean_11']

labels_time = [
    "Embryonic", "Early fetal", "Early mid-fetal", "Late mid-fetal", "Late fetal",
    "Early infancy", "Late infancy", "Early Childhood", "Late childhood", "Adolescence",
    "Young adulthood", "Adulthood"
]

gene_set_styles = {
    "ASD w/o ID": dict(color=palette[0], marker='o', label="ASD genes"),
    "ASD with ID": dict(color=palette[1], marker='^', label="ASD with ID genes"),
    "SCZ": dict(color=palette[2], marker='s', label="SCZ genes"),
    "ALL": dict(color=palette[4], marker='*', label="ALL Genes", ls="--", markersize=12, alpha=0.8, zorder=2),
}

fig, ax = plt.subplots(dpi=300, figsize=(12, 8), constrained_layout=True)

xvals = np.arange(len(Time))
for gene_set, style in gene_set_styles.items():
    subset = brainspan_df[brainspan_df["gene_set"] == gene_set].sort_values("stage", key=lambda s: s.map({t: i for i, t in enumerate(Time)}))
    shift = {"ASD w/o ID": 0.05, "SCZ": -0.05}.get(gene_set, 0.0)
    ax.errorbar(
        x=xvals + shift,
        y=subset["mean"].values,
        yerr=subset["sem"].values,
        linewidth=2, markersize=style.get("markersize", 8),
        alpha=style.get("alpha", 0.95), zorder=style.get("zorder", 3),
        capsize=5, elinewidth=2, capthick=2,
        color=style["color"], marker=style["marker"], label=style["label"],
        ls=style.get("ls", "-"),
    )

ax.axvline(x=5, color='grey', ls="--", linewidth=1.5, zorder=1)

ax.set_xticks(xvals)
ax.set_xticklabels(labels_time, rotation=30, ha='right', fontsize=13, weight='bold')
ax.set_ylabel("Brain expression level (log10(RPKM))", fontsize=18, weight='bold')
ax.set_xlabel("Developmental stage", fontsize=16, weight='bold', labelpad=10)

sns.despine(ax=ax)
ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.5, zorder=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.7, 1), fancybox=True, shadow=True,
          ncol=2, fontsize=14, frameon=True, borderaxespad=0.5, title="Gene sets", title_fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=13, width=1.5)

fig.savefig(FIG_DIR + "FigureS2.pdf", dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.show()
print(f"Saved: {FIG_DIR}FigureS2.pdf")

# %% [markdown]
# ---
# ## Figure S3 — BrainSpan gene removal impact on ASD-SCZ bias correlation
#
# Pre-generated by Similarity_ASD_SCZ.spec notebook.

# %%
import shutil
from PIL import Image as PILImage

BRAINSPAN_REMOVAL_SRC = os.path.join(ProjDIR, "results/figures/BrainSpan_gene_removal_ASD_SCZ_correlation.png")
shutil.copy2(BRAINSPAN_REMOVAL_SRC, FIG_DIR + "FigureS3.png")

img = PILImage.open(FIG_DIR + "FigureS3.png")
img.save(FIG_DIR + "FigureS3.pdf", "PDF", resolution=300)
print(f"Saved: {FIG_DIR}FigureS3.pdf")

from IPython.display import Image as IPImage, display as ipdisplay
ipdisplay(IPImage(filename=FIG_DIR + "FigureS3.png"))

# %% [markdown]
# ---
# ## Figure S5 — SuperCluster Bias BoxPlots
#
# | Panel | Content |
# |-------|---------|
# | A | ASD (all) |
# | B | ASD w/o ID |
# | C | ASD with ID |
# | D | SCZ |
# | E | DD/ID |

# %% [markdown]
# ### Panel A — ASD

# %%
with save_panel(FIG_DIR + "FigS5_A.png"):
    SuperClusterBias_BoxPlot(ASD_All_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.05)

# %% [markdown]
# ### Panel B — ASD w/o ID

# %%
with save_panel(FIG_DIR + "FigS5_B.png"):
    SuperClusterBias_BoxPlot(HighIQ_ASD_Bias, "ASD w/o ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.1)

# %% [markdown]
# ### Panel C — ASD with ID

# %%
with save_panel(FIG_DIR + "FigS5_C.png"):
    SuperClusterBias_BoxPlot(LowIQ_ASD_Bias, "ASD with ID", NeuroOnly=False, sortby="mean", EffectCol="-logP")

# %% [markdown]
# ### Panel D — SCZ

# %%
with save_panel(FIG_DIR + "FigS5_D.png"):
    SuperClusterBias_BoxPlot(SCZ_Bias, "SCZ", NeuroOnly=False, sortby="mean", EffectCol="-logP")

# %% [markdown]
# ### Panel E — DD/ID

# %%
with save_panel(FIG_DIR + "FigS5_E.png"):
    SuperClusterBias_BoxPlot(DDD_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP")

# %% [markdown]
# ### Composite Figure S5

# %%
assemble_panels(
    panel_paths=[
        FIG_DIR + "FigS5_A.png",
        FIG_DIR + "FigS5_B.png",
        FIG_DIR + "FigS5_C.png",
        FIG_DIR + "FigS5_D.png",
        FIG_DIR + "FigS5_E.png",
    ],
    panel_labels=["A", "B", "C", "D", "E"],
    layout=[
        (0, slice(0, 2)),   # A — left
        (0, slice(2, 4)),   # B — right
        (1, slice(0, 2)),   # C — left
        (1, slice(2, 4)),   # D — right
        (2, slice(0, 2)),   # E — left (last row, half width)
    ],
    figsize=(20, 26),
    out_path=FIG_DIR + "FigureS5.pdf",
)

# %% [markdown]
# ---
# ## Figure S7 — Individual cell type comparisons
#
# | Panel | Content |
# |-------|---------|
# | A | Hippo CA1-3 (ASD w/o ID vs SCZ) |
# | B | Upper IT (ASD w/o ID vs SCZ) |
# | C | Deep IT (ASD w/o ID vs SCZ) |
# | D | Deep CT6b (ASD w/o ID vs SCZ) |
# | E | ASD w/o ID vs ASD with ID bar plot |
# | F | CGE: ASD with ID vs SCZ |
# | G | CGE: VNR+ vs VNR- |
# | H | CGE: ASD w/o ID vs DD/ID |

# %% [markdown]
# ### Panel A — Hippocampal CA1-3

# %%
with save_panel(FIG_DIR + "FigS7_A.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Hippocampal CA1-3", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.15, 0.12))

# %% [markdown]
# ### Panel B — Upper-layer IT

# %%
with save_panel(FIG_DIR + "FigS7_B.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Upper-layer intratelencephalic", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, pval="Mann_Whitney_FDR", loc=(0.15, 0.12))

# %% [markdown]
# ### Panel C — Deep-layer IT

# %%
with save_panel(FIG_DIR + "FigS7_C.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer intratelencephalic", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.3), pval="Mann_Whitney_FDR")

# %% [markdown]
# ### Panel D — Deep-layer CT6b

# %%
with save_panel(FIG_DIR + "FigS7_D.png"):
    CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer corticothalamic and 6b", ASD_SCZ_Contrast,
                    "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.1, 0.05))

# %% [markdown]
# ### Panel E — ASD w/o ID vs ASD with ID bar plot

# %%
with save_panel(FIG_DIR + "FigS7_E.png"):
    plot_bias_comparison(HIQ_LIQ_Contrast_Neurons, "ASD w/o ID", "ASD with ID",
                         p_test="Mann_Whitney_FDR", legend_anchor=(0.9, 1.0))

# %% [markdown]
# ### Panel F — CGE: ASD with ID vs SCZ

# %%
with save_panel(FIG_DIR + "FigS7_F.png"):
    CompareSingleCT(LowIQ_ASD_Bias, SCZ_Bias, "CGE interneuron", ASD_wID_SCZ_Contrast,
                    "ASD with ID Mutation Bias", "SCZ Mutation Bias", loc=(0.08, 0.21))

# %% [markdown]
# ### Panel G — CGE: VNR+ vs VNR-

# %%
with save_panel(FIG_DIR + "FigS7_G.png"):
    CompareSingleCT(VNR_Pos_Bias, VNR_Neg_Bias, "CGE interneuron", VNR_Contrast,
                    "VNR + Mutation Bias", "VNR - Mutation Bias", loc=(0.0, -0.05))

# %% [markdown]
# ### Panel H — CGE: ASD w/o ID vs DD/ID

# %%
with save_panel(FIG_DIR + "FigS7_H.png"):
    CompareSingleCT(HighIQ_ASD_Bias, DDD_Bias, "CGE interneuron", ASD_DDD_Contrast,
                    "ASD w/o ID Mutation Bias", "DD/ID Mutation Bias", loc=(0.1, 0.05))

# %% [markdown]
# ### Composite Figure S7

# %%
assemble_panels(
    panel_paths=[
        FIG_DIR + "FigS7_A.png", FIG_DIR + "FigS7_B.png",
        FIG_DIR + "FigS7_C.png", FIG_DIR + "FigS7_D.png",
        FIG_DIR + "FigS7_E.png", FIG_DIR + "FigS7_F.png",
        FIG_DIR + "FigS7_G.png", FIG_DIR + "FigS7_H.png",
    ],
    panel_labels=["A", "B", "C", "D", "E", "F", "G", "H"],
    layout=[
        (0, slice(0, 2)), (0, slice(2, 4)),   # A, B
        (1, slice(0, 2)), (1, slice(2, 4)),   # C, D
        (2, slice(0, 4)),                      # E (full width)
        (3, slice(0, 2)),                      # F
        (3, slice(2, 4)),                      # G (was col 2-3)
        (4, slice(0, 2)),                      # H
    ],
    figsize=(20, 36),
    out_path=FIG_DIR + "FigureS7.pdf",
)

# %% [markdown]
# ---
# ## Figure S8 — Multi-disorder comparison (MGE, LAMP5)
#
# | Panel | Content |
# |-------|---------|
# | A | MGE interneuron |
# | B | LAMP5-LHX6 and Chandelier |

# %%
datasets = {
    'ASD w/o ID': HighIQ_ASD_Bias,
    'ASD with ID': LowIQ_ASD_Bias,
    'VNR+': VNR_Pos_Bias,
    'VNR-': VNR_Neg_Bias,
    'DD/ID': DDD_Bias,
    'SCZ': SCZ_Bias
}
TestPairs = [("VNR+", "VNR-"), ("ASD w/o ID", "SCZ"), ("ASD w/o ID", "ASD with ID"),
             ("SCZ", "ASD with ID"), ("ASD w/o ID", "DD/ID")]

# %% [markdown]
# ### Panel A — MGE interneuron

# %%
with save_panel(FIG_DIR + "FigS8_A.png"):
    plot_mutation_bias_comparison_V2("MGE interneuron", datasets, Anno, all_contrasts_df, TestPairs=TestPairs)

# %% [markdown]
# ### Panel B — LAMP5-LHX6 and Chandelier

# %%
with save_panel(FIG_DIR + "FigS8_B.png"):
    plot_mutation_bias_comparison_V2("LAMP5-LHX6 and Chandelier", datasets, Anno, all_contrasts_df, TestPairs=TestPairs)

# %% [markdown]
# ### Composite Figure S8

# %%
assemble_panels(
    panel_paths=[FIG_DIR + "FigS8_A.png", FIG_DIR + "FigS8_B.png"],
    panel_labels=["A", "B"],
    layout=[(0, slice(0, 1)), (0, slice(1, 2))],
    figsize=(14, 7),
    out_path=FIG_DIR + "FigureS8.pdf",
)

# %% [markdown]
# ---
# ## Figure S9 — 22q11.2 SuperCluster Bias BoxPlot (single panel)

# %%
with save_panel(FIG_DIR + "FigureS9.png"):
    SuperClusterBias_BoxPlot(X22q_Bias, "22q11.2", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=0.1)

# Save single-panel figure as PDF too
from PIL import Image as PILImage
img = PILImage.open(FIG_DIR + "FigureS9.png")
img.save(FIG_DIR + "FigureS9.pdf", "PDF", resolution=300)
print(f"Saved: {FIG_DIR}FigureS9.pdf")

# %% [markdown]
# ---
# ## Supplementary Tables S2-S7: Cluster-level biases

# %%
import openpyxl

def prepare_bias_table_v1(df):
    """Prepare a bias DataFrame for supplementary table output."""
    df_out = df.copy(deep=True)
    df_out.index.name = "Cluster"
    columns_to_keep = [
        "EFFECT", "P-value", "q-value", "Class", "Supercluster", "Subtype", "Neurotransmitter",
        "Top three regions", "Top three dissections", "Number of cells"]
    df_out = df_out[columns_to_keep]
    df_out.rename(columns={"EFFECT": "Bias"}, inplace=True)
    df_out = df_out.sort_values(by=["P-value", "Bias"], ascending=[True, False])
    return df_out

def prepare_bias_table_v2(BiasPos, BiasNeg, Name1, Name2):
    df_out = pd.DataFrame()
    df_out["Bias {}".format(Name1)] = BiasPos["EFFECT"]
    df_out["P-value {}".format(Name1)] = BiasPos["P-value"]
    df_out["q-value {}".format(Name1)] = BiasPos["q-value"]
    df_out["Bias {}".format(Name2)] = BiasNeg["EFFECT"]
    df_out["P-value {}".format(Name2)] = BiasNeg["P-value"]
    df_out["q-value {}".format(Name2)] = BiasNeg["q-value"]
    df_out["Class"] = BiasPos["Class"]
    df_out["Supercluster"] = BiasPos["Supercluster"]
    df_out["Subtype"] = BiasPos["Subtype"]
    df_out["Neurotransmitter"] = BiasPos["Neurotransmitter"]
    df_out["Top three regions"] = BiasPos["Top three regions"]
    df_out["Top three dissections"] = BiasPos["Top three dissections"]
    df_out["Number of cells"] = BiasPos["Number of cells"]
    df_out.index.name = "Cluster"
    df_out = df_out.sort_values(by=["P-value {}".format(Name2), "Bias {}".format(Name2)], ascending=[True, False])
    return df_out

SCZ_Bias_toST = prepare_bias_table_v1(SCZ_Bias)
HighIQ_ASD_Bias_toST = prepare_bias_table_v1(HighIQ_ASD_Bias)
LowIQ_ASD_Bias_toST = prepare_bias_table_v1(LowIQ_ASD_Bias)
X22q_Bias_toST = prepare_bias_table_v1(X22q_Bias)
DDD_Bias_toST = prepare_bias_table_v1(DDD_Bias)
UKBB_VNR_Bias_toST = prepare_bias_table_v2(VNR_Pos_Bias, VNR_Neg_Bias, "VNR+", "VNR-")

# %%
SuppTabOutDir = os.path.join(ProjDIR, "dat/suppl.data/")
excel_path = SuppTabOutDir + "SupTab_VIP.xlsx"

# Remove existing data sheets (keep Table_of_contents and experiment sheets)
wb = openpyxl.load_workbook(excel_path)
sheets_to_remove = [s for s in wb.sheetnames if s not in [
    "Table_of_contents",
    "Mouse Experiment Statistical detail  Statistical detail summary.",
    "Experiement Statistical detail "
]]
for s in sheets_to_remove:
    wb.remove(wb[s])
wb.save(excel_path)
wb.close()

# %%
with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    Name_Dict = {
        "SCZ_Bias": "Table_S2_Cluster_Bias_SCZ",
        "HighIQ_ASD_Bias": "Table_S3_Cluster_Bias_ASD_woID",
        "LowIQ_ASD_Bias": "Table_S4_Cluster_Bias_ASD_ID",
        "X22q_Bias": "Table_S5_Cluster_Bias_22q11.2",
        "DDD_Bias": "Table_S6_Cluster_Bias_DD",
        "VNR_Bias": "Table_S7_Cluster_Bias_VNR",
    }
    DF_list = [
        ("SCZ_Bias", SCZ_Bias_toST),
        ("HighIQ_ASD_Bias", HighIQ_ASD_Bias_toST),
        ("LowIQ_ASD_Bias", LowIQ_ASD_Bias_toST),
        ("X22q_Bias", X22q_Bias_toST),
        ("DDD_Bias", DDD_Bias_toST),
        ("VNR_Bias", UKBB_VNR_Bias_toST),
    ]
    for df_name, DF in DF_list:
        sheet_name = Name_Dict.get(df_name, df_name)
        DF.to_excel(writer, sheet_name=sheet_name)

print("Tables S2-S7 written.")

# %% [markdown]
# ---
# ## Supplementary Tables S8-S11: SuperCluster bias contrasts

# %%
def process_BiasContrast_df(df):
    columns_to_drop = ["Wilcoxon_P", "Wilcoxon_FDR", "Bonferroni_P"]
    df = df.drop(columns=[col for col in columns_to_drop if col in df.columns])
    disorder_name_1 = df.columns[0].replace("Bias_", "")
    disorder_name_2 = df.columns[1].replace("Bias_", "")
    df.rename(columns={"Bias_Diff": "Bias_Diff_{}_{}".format(disorder_name_1, disorder_name_2)}, inplace=True)
    return df

ASD_SCZ_Contrast_toST = process_BiasContrast_df(ASD_SCZ_Contrast_Neurons)
HIQ_LIQ_Contrast_toST = process_BiasContrast_df(HIQ_LIQ_Contrast_Neurons)
DDD_ASD_Contrast_toST = process_BiasContrast_df(ASD_DDD_Contrast_Neurons)
VNR_Contrast_toST = process_BiasContrast_df(VNR_Contrast_Neurons)

# %%
with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    Name_Dict = {
        "BiasContrast_ASD_SCZ": "Table_S8_BiasContrast_ASD_SCZ",
        "BiasContrast_HIQ_LIQ": "Table_S9_BiasContrast_HIQ_LIQ",
        "BiasContrast_HIQ_DDD": "Table_S10_BiasContrast_HIQ_DDD",
        "BiasContrast_VNR": "Table_S11_BiasContrast_VNR",
    }
    DF_list = [
        ("BiasContrast_ASD_SCZ", ASD_SCZ_Contrast_toST),
        ("BiasContrast_HIQ_LIQ", HIQ_LIQ_Contrast_toST),
        ("BiasContrast_HIQ_DDD", DDD_ASD_Contrast_toST),
        ("BiasContrast_VNR", VNR_Contrast_toST),
    ]
    for df_name, DF in DF_list:
        sheet_name = Name_Dict.get(df_name, df_name)
        DF.to_excel(writer, sheet_name=sheet_name)

print("Tables S8-S11 written.")
