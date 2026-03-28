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
#     display_name: Python (gencic)
#     language: python
#     name: gencic
# ---

# %% [markdown]
# # Exploratory: Bias Correlations and Contrasts
#
# Scatter plots, overall contrasts, individual cell type comparisons, VNR/EDU contrasts.
# This is an exploratory notebook for interactive analysis — not for final figures.

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os

from pathlib import Path
import yaml
with open("/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/config/config.yaml") as f:
    _cfg = yaml.safe_load(f)
PROJ_DIR = Path(_cfg["ProjDIR"])
sys.path.insert(0, str(PROJ_DIR / "src"))
from CellType_PSY import *

plt.style.use('seaborn-v0_8-whitegrid')

# %% [markdown]
# ## Load data

# %%
Bias_Save_Dir = str(PROJ_DIR / "results/main_results/random/Centering/") + "/"

HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)
VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)

# Pre-computed contrasts
CONTRAST_DIR = str(PROJ_DIR / "results/main_results/contrasts/") + "/"
ASD_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_SCZ_contrast.csv", index_col=0)
ASD_wID_SCZ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_wID_vs_SCZ_contrast.csv", index_col=0)
HIQ_LIQ_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_ASD_wID_contrast.csv", index_col=0)
ASD_DDD_Contrast = pd.read_csv(CONTRAST_DIR + "ASD_woID_vs_DDD_contrast.csv", index_col=0)
VNR_Contrast = pd.read_csv(CONTRAST_DIR + "VNR_neg_vs_pos_contrast.csv", index_col=0)
EDU_Contrast = pd.read_csv(CONTRAST_DIR + "EDU_neg_vs_pos_contrast.csv", index_col=0)

# %% [markdown]
# ## Bias Correlations (CompareCT scatter plots)

# %%
HIQ_ASD_SCZ, fig = CompareCT(HighIQ_ASD_Bias, SCZ_Bias, "ASD w/o ID", "SCZ", effectlabel="EFFECT", SuperClusters=ALL_CTs)

# %%
LIQ_ASD_SCZ, fig = CompareCT(LowIQ_ASD_Bias, SCZ_Bias, "ASD with ID", "SCZ", effectlabel="EFFECT", SuperClusters=ALL_CTs)

# %%
HIQ_LIQ_ASD, fig = CompareCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "ASD w/o ID", "ASD with ID", effectlabel="EFFECT", SuperClusters=ALL_CTs)

# %% [markdown]
# ## Bias Contrast Plots (PlotBiasContrast)

# %%
eff_label = "EFFECT"
PlotBiasContrast(HIQ_ASD_SCZ, label1="{}_ASD w/o ID".format(eff_label), label2="{}_SCZ".format(eff_label),
                 name1="ASD w/o ID Mutation Bias", name2="SCZ Mutation Bias",
                 title="", neur_only=True)

# %%
PlotBiasContrast(LIQ_ASD_SCZ, label1="{}_ASD with ID".format(eff_label), label2="{}_SCZ".format(eff_label),
                 name1="ASD with ID Mutation Bias", name2="SCZ Mutation Bias",
                 title="")

# %% [markdown]
# ## Overall Contrasts (bar plots)

# %%
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
plot_bias_comparison(ASD_SCZ_Contrast_Neurons, "ASD w/o ID", "SCZ",
                     p_test="Mann_Whitney_FDR", legend_anchor=(0.15, 0.9))

# %%
ASD_wID_SCZ_Contrast_Neurons = ASD_wID_SCZ_Contrast[ASD_wID_SCZ_Contrast.index.isin(Neurons)]
plot_bias_comparison(ASD_wID_SCZ_Contrast_Neurons, "ASD with ID", "SCZ",
                     p_test="Mann_Whitney_FDR", legend_anchor=(0.15, 1.0))

# %%
HIQ_LIQ_Contrast_Neurons = HIQ_LIQ_Contrast[HIQ_LIQ_Contrast.index.isin(Neurons)]
plot_bias_comparison(HIQ_LIQ_Contrast_Neurons, "ASD w/o ID", "ASD with ID",
                     p_test="Mann_Whitney_FDR", legend_anchor=(0.9, 1.0))

# %% [markdown]
# ## Individual CT Contrasts

# %%
EffLabel = "EFFECT"

# ASD w/o ID vs SCZ — various cell types
CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Hippocampal CA1-3", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.15, 0.12))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Upper-layer intratelencephalic", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, pval="Mann_Whitney_FDR", loc=(0.15, 0.12))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer intratelencephalic", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.3), pval="Mann_Whitney_FDR")

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Amygdala excitatory", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.15, 0.29))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Medium spiny neuron", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, pval="Mann_Whitney_FDR", loc=(0.05, 0.23))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "Deep-layer corticothalamic and 6b", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", efflabel=EffLabel, loc=(0.1, 0.05))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "MGE interneuron", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.05))

CompareSingleCT(HighIQ_ASD_Bias, SCZ_Bias, "CGE interneuron", ASD_SCZ_Contrast,
                "ASD w/o ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.05))

# %%
# ASD with ID vs SCZ
CompareSingleCT(LowIQ_ASD_Bias, SCZ_Bias, "MGE interneuron", ASD_wID_SCZ_Contrast,
                "ASD with ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.05))

CompareSingleCT(LowIQ_ASD_Bias, SCZ_Bias, "CGE interneuron", ASD_wID_SCZ_Contrast,
                "ASD with ID Mutation Bias", "SCZ Mutation Bias", loc=(0.1, 0.25))

# %% [markdown]
# ## VNR / EDU contrasts

# %%
CompareSingleCT(VNR_Pos_Bias, VNR_Neg_Bias, "CGE interneuron", VNR_Contrast,
                "VNR + Mutation Bias", "VNR - Mutation Bias", loc=(0.0, -0.05))

CompareSingleCT(EDU_Pos_Bias, EDU_Neg_Bias, "CGE interneuron", EDU_Contrast,
                "EDU + Mutation Bias", "EDU - Mutation Bias", loc=(0.05, 0.02))

# %% [markdown]
# ## ASD w/o ID vs DD/ID — CGE

# %%
CompareSingleCT(HighIQ_ASD_Bias, DDD_Bias, "CGE interneuron", ASD_DDD_Contrast,
                "ASD w/o ID Mutation Bias", "DD/ID Mutation Bias", loc=(0.1, 0.05))

# %% [markdown]
# ## ASD w/o ID vs ASD with ID — various cell types

# %%
CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "CGE interneuron",
                HIQ_LIQ_Contrast_Neurons, "ASD w/o ID Mutation Bias",
                "ASD with ID\nMutation Bias", loc=(0.12, 0.05), efflabel=EffLabel)

CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "Deep-layer corticothalamic and 6b",
                HIQ_LIQ_Contrast_Neurons, "ASD w/o ID Mutation Bias",
                "ASD with ID\nMutation Bias", loc=2, efflabel=EffLabel)

CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "Amygdala excitatory",
                HIQ_LIQ_Contrast_Neurons, "ASD w/o ID Mutation Bias",
                "ASD with ID\nMutation Bias", efflabel=EffLabel, loc=(0.25, 0.17))

CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "Upper-layer intratelencephalic",
                HIQ_LIQ_Contrast_Neurons, "ASD w/o ID Mutation Bias",
                "ASD with ID\nMutation Bias", loc=2, efflabel=EffLabel)

CompareSingleCT(HighIQ_ASD_Bias, LowIQ_ASD_Bias, "Deep-layer intratelencephalic",
                HIQ_LIQ_Contrast_Neurons, "ASD w/o ID Mutation Bias",
                "ASD with ID\nMutation Bias", loc=2, efflabel=EffLabel)
