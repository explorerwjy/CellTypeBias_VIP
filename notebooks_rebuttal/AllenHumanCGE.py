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

# %%
# %load_ext autoreload
# %autoreload 2
import sys
import os

ProjDIR = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/" # Change to your project directory
sys.path.insert(1, f'{ProjDIR}/src/')
from CellType_PSY import *
import scanpy as sc
HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()

try:
    #os.chdir(f"{ProjDIR}/notebooks_rebuttal/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:
    print(f"Unexpected error: {e}")

import anndata as ad
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection, multipletests

# %%
CGE_raw = "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_CGE_interneuron.h5ad"
#CGE_raw = "/mnt/data0/HumanBrainCellType/SuperCluster_CPM/Supercluster_CGE_interneuron.logCPM.h5ad"
CGE_backed = ad.read_h5ad(CGE_raw, backed='r')

# %%
# Check var columns and set feature_name as index
print("Current var columns:", list(CGE_backed.var.columns))
print("Current var index name:", CGE_backed.var.index.name)
print("Shape:", CGE_backed.var.shape)

# Check if feature_name column exists
if 'feature_name' in CGE_backed.var.columns:
    # Set feature_name as index
    CGE_backed.var.set_index('feature_name', inplace=True)
    print("\nSet 'feature_name' as var index")
    print("New var index name:", CGE_backed.var.index.name)
    print("Sample var index:", list(CGE_backed.var.index[:5]))
else:
    print("\nWarning: 'feature_name' column not found in var")
    print("Available columns:", list(CGE_backed.var.columns))

# %%
CGE_backed

# %%
CGE_backed.var

# %%
CGE_backed.obs.head(2)

# %%

# %%
Makders = ["VIP", "CCK", "SNCG", "CHRM2", "CALB2"]
for marker in Makders:
    print(marker in CGE_backed.var.index)

# %%
import matplotlib.pyplot as plt
import numpy as np

# Fix typo in variable name - use uppercase markers for human
Markers = Makders if 'Makders' in locals() else ["VIP", "CCK", "SNCG", "CHRM2", "CALB2"]

# Create figure with subplots
n_markers = len(Markers)
n_cols = 3
n_rows = (n_markers + n_cols - 1) // n_cols  # Ceiling division

fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
axes = axes.flatten() if n_markers > 1 else [axes]

for idx, marker in enumerate(Markers):
    ax = axes[idx]
    
    # Check if marker exists in the index (CGE_backed for human data)
    if marker in CGE_backed.var.index:
        # Get expression values for this marker
        marker_idx = CGE_backed.var.index.get_loc(marker)
        
        # Extract expression values (handle sparse matrix and backed mode)
        # For backed objects, X is accessed directly and may be sparse
        expr_slice = CGE_backed.X[:, marker_idx]
        if hasattr(expr_slice, 'toarray'):
            expression = expr_slice.toarray().flatten()
        elif hasattr(expr_slice, 'A1'):  # scipy sparse matrix
            expression = expr_slice.A1  # A1 is flattened array
        else:
            expression = np.array(expr_slice).flatten()
        
        # Plot histogram
        ax.hist(expression, bins=50, edgecolor='black', alpha=0.7)
        ax.set_title(f'{marker} Expression Distribution', fontsize=12, fontweight='bold')
        ax.set_xlabel('Expression Level', fontsize=10)
        ax.set_ylabel('Number of Cells', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        mean_expr = np.mean(expression)
        median_expr = np.median(expression)
        ax.axvline(mean_expr, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_expr:.2f}')
        ax.axvline(median_expr, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_expr:.2f}')
        ax.legend(fontsize=8)
        
        print(f"{marker}: Mean={mean_expr:.2f}, Median={median_expr:.2f}, Max={np.max(expression):.2f}, Non-zero cells={np.sum(expression > 0)}/{len(expression)} ({np.sum(expression > 0)/len(expression)*100:.1f}%)")
    else:
        ax.text(0.5, 0.5, f'{marker}\nnot found', 
                ha='center', va='center', fontsize=14, color='red')
        ax.set_title(f'{marker} - Not Found', fontsize=12, fontweight='bold')

# Hide unused subplots
for idx in range(n_markers, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.show()


# %%

# %%

# %%

# %%
# Extract expression values for each marker (adapted for human CGE_backed)
def get_expression(adata, gene_name):
    """Extract expression values for a gene (works with backed mode)"""
    if gene_name not in adata.var.index:
        print(f"Warning: {gene_name} not found in var.index")
        return None
    marker_idx = adata.var.index.get_loc(gene_name)
    
    # Handle backed mode and sparse matrices
    expr_slice = adata.X[:, marker_idx]
    if hasattr(expr_slice, 'toarray'):
        return expr_slice.toarray().flatten()
    elif hasattr(expr_slice, 'A1'):  # scipy sparse matrix
        return expr_slice.A1  # A1 is flattened array
    else:
        return np.array(expr_slice).flatten()

# Use uppercase gene names for human (VIP, CCK, SNCG, CALB2)
# Get expression for all markers from CGE_backed
Vip_expr = get_expression(CGE_backed, "VIP")
Cck_expr = get_expression(CGE_backed, "CCK")
Sncg_expr = get_expression(CGE_backed, "SNCG")
Calb2_expr = get_expression(CGE_backed, "CALB2")

# Define cutoff
CPM_Cut_off = 0.1

# Create boolean masks for each subset
# CCK_BC: CCK > CPM_Cut_off; SNCG > CPM_Cut_off; CALB2 < CPM_Cut_off
CCK_BC_mask = (Cck_expr > CPM_Cut_off) & (Sncg_expr > CPM_Cut_off) & (Calb2_expr < CPM_Cut_off)

# ISI2: VIP > CPM_Cut_off; CCK < CPM_Cut_off; SNCG < CPM_Cut_off; CALB2 < CPM_Cut_off
ISI2_mask = (Vip_expr > CPM_Cut_off) & (Sncg_expr < CPM_Cut_off) & (Calb2_expr < CPM_Cut_off)

# ISI3: VIP > CPM_Cut_off; CCK < CPM_Cut_off; SNCG < CPM_Cut_off; CALB2 > CPM_Cut_off
ISI3_mask = (Vip_expr > CPM_Cut_off) & (Sncg_expr < CPM_Cut_off) & (Calb2_expr > CPM_Cut_off)

# Make subsets mutually exclusive (prioritize CCK_BC, then ISI2, then ISI3)
# Start with ISI2 and ISI3 (they don't overlap with CCK_BC)
ISI2_mask_exclusive = ISI2_mask & ~CCK_BC_mask
ISI3_mask_exclusive = ISI3_mask & ~CCK_BC_mask

# Create subset AnnData objects (use to_memory() for backed objects)
CGE_CCK_BC = CGE_backed[CCK_BC_mask].to_memory()
CGE_ISI2 = CGE_backed[ISI2_mask_exclusive].to_memory()
CGE_ISI3 = CGE_backed[ISI3_mask_exclusive].to_memory()

# Print summary
print("=" * 60)
print("CGE Subset Classification Results (Human)")
print("=" * 60)
print(f"CPM Cutoff: {CPM_Cut_off}")
print()
print(f"CCK_BC: {CCK_BC_mask.sum()} cells")
print(f"  Criteria: CCK > {CPM_Cut_off}, SNCG > {CPM_Cut_off}, CALB2 < {CPM_Cut_off}")
print()
print(f"ISI2: {ISI2_mask_exclusive.sum()} cells")
print(f"  Criteria: VIP > {CPM_Cut_off}, CCK < {CPM_Cut_off}, SNCG < {CPM_Cut_off}, CALB2 < {CPM_Cut_off}")
print()
print(f"ISI3: {ISI3_mask_exclusive.sum()} cells")
print(f"  Criteria: VIP > {CPM_Cut_off}, CCK < {CPM_Cut_off}, SNCG < {CPM_Cut_off}, CALB2 > {CPM_Cut_off}")
print()
print(f"Total classified: {CCK_BC_mask.sum() + ISI2_mask_exclusive.sum() + ISI3_mask_exclusive.sum()} cells")
print(f"Unclassified: {CGE_backed.n_obs - (CCK_BC_mask.sum() + ISI2_mask_exclusive.sum() + ISI3_mask_exclusive.sum())} cells")
print("=" * 60)

# Store masks in a dictionary for easy access
CGE_subset_masks = {
    'CCK_BC': CCK_BC_mask,
    'ISI2': ISI2_mask_exclusive,
    'ISI3': ISI3_mask_exclusive
}

# %%
# Report ROI (tissue) distribution for each subset (Human CGE)
import pandas as pd

# Create a summary DataFrame
roi_summary = []

# Get ROI counts for each subset
subsets = {
    'CCK_BC': CGE_CCK_BC,
    'ISI2': CGE_ISI2,
    'ISI3': CGE_ISI3
}

for subset_name, subset_adata in subsets.items():
    roi_counts = subset_adata.obs['tissue'].value_counts()
    total = subset_adata.n_obs
    
    for roi, count in roi_counts.items():
        roi_summary.append({
            'Subset': subset_name,
            'Tissue': roi,
            'Count': count,
            'Percentage': (count / total * 100) if total > 0 else 0,
            'Total': total
        })

# Get unclassified cells
unclassified_mask = ~(CCK_BC_mask | ISI2_mask_exclusive | ISI3_mask_exclusive)
unclassified_adata = CGE_backed[unclassified_mask].to_memory()
unclassified_roi_counts = unclassified_adata.obs['tissue'].value_counts()
unclassified_total = unclassified_adata.n_obs

for roi, count in unclassified_roi_counts.items():
    roi_summary.append({
        'Subset': 'Unclassified',
        'Tissue': roi,
        'Count': count,
        'Percentage': (count / unclassified_total * 100) if unclassified_total > 0 else 0,
        'Total': unclassified_total
    })

# Create DataFrame and display
roi_df = pd.DataFrame(roi_summary)

# Pivot table for better visualization
roi_pivot = roi_df.pivot_table(
    index='Subset',
    columns='Tissue',
    values='Count',
    fill_value=0
)

# Add total column
roi_pivot['Total'] = roi_pivot.sum(axis=1)

print("=" * 80)
print("Tissue Distribution by CGE Subset (Human)")
print("=" * 80)
print()
print(roi_pivot)
print()

# Also show percentages
roi_pct = roi_df.pivot_table(
    index='Subset',
    columns='Tissue',
    values='Percentage',
    fill_value=0
)
roi_pct['Total_Cells'] = roi_pivot['Total']

print("=" * 80)
print("Tissue Distribution (Percentages)")
print("=" * 80)
print()
for subset in roi_pct.index:
    print(f"\n{subset}:")
    total = roi_pct.loc[subset, 'Total_Cells']
    for tissue in roi_pct.columns:
        if tissue != 'Total_Cells':
            count = roi_pivot.loc[subset, tissue]
            pct = roi_pct.loc[subset, tissue]
            print(f"  {tissue}: {count} cells ({pct:.1f}%)")
    print(f"  Total: {int(total)} cells")

print()
print("=" * 80)
print("Overall Tissue Distribution in CGE_backed")
print("=" * 80)
overall_tissue = CGE_backed.obs['tissue'].value_counts()
for tissue, count in overall_tissue.items():
    pct = count / CGE_backed.n_obs * 100
    print(f"{tissue}: {count} cells ({pct:.1f}%)")
print(f"Total: {CGE_backed.n_obs} cells")
print("=" * 80)

# %%
# Define important 22q11.2 genes in human genome
# These are human gene symbols (uppercase)
important_22q11_genes = [
    'COMT',      # Catechol-O-methyltransferase - critical for dopamine metabolism
    'DGCR8',     # DiGeorge syndrome critical region gene 8 - microRNA processing
    'PRODH',     # Proline dehydrogenase - associated with schizophrenia
    'TBX1',      # T-box transcription factor 1 - important for development
    'GP1BB',     # Glycoprotein Ib beta - platelet function
    'UFD1L',     # Ubiquitin fusion degradation 1 like
    'HIRA',      # Histone cell cycle regulator
    'TXNRD2',    # Thioredoxin reductase 2
    'SLC25A1',   # Solute carrier family 25 member 1
    'DGCR2',     # DiGeorge syndrome critical region gene 2
    'DGCR14',    # DiGeorge syndrome critical region gene 14
    'SEPT5',     # Septin 5
    'SNAP29',    # Synaptosome associated protein 29
    'ZDHHC8',    # Zinc finger DHHC-type containing 8
    'RTN4R',     # Reticulon 4 receptor
    'ARVCF',     # Armadillo repeat containing, catenin delta
    'GSCL',      # Goosecoid like
    'CLDN5',     # Claudin 5
    'HIC2',      # HIC ZBTB transcriptional repressor 2
    'UBA1',      # Ubiquitin like modifier activating enzyme 1
    'P2RX6'      # Purinergic receptor P2X 6
]

# Check which genes are available in the dataset
available_genes = []
missing_genes = []

for gene in important_22q11_genes:
    if gene in CGE_backed.var.index:
        available_genes.append(gene)
    else:
        missing_genes.append(gene)

print("=" * 80)
print("22q11.2 Genes Availability Check (Human)")
print("=" * 80)
print(f"\nAvailable genes ({len(available_genes)}): {available_genes}")
print(f"\nMissing genes ({len(missing_genes)}): {missing_genes}")

# Select top 5 most important genes that are available
# Priority order: COMT, DGCR8, PRODH, TBX1, then others
priority_order = ['COMT', 'DGCR8', 'PRODH', 'TBX1', 'GP1BB', 'UFD1L', 'HIRA', 
                  'TXNRD2', 'SLC25A1', 'DGCR2', 'ZDHHC8', 'RTN4R', 'ARVCF']

selected_genes = []
for gene in priority_order:
    if gene in available_genes and len(selected_genes) < 5:
        selected_genes.append(gene)

# If we don't have 5 yet, add more from available_genes
if len(selected_genes) < 5:
    for gene in available_genes:
        if gene not in selected_genes:
            selected_genes.append(gene)
            if len(selected_genes) >= 5:
                break

print(f"\nSelected top 5 genes for plotting: {selected_genes}")
print("=" * 80)

# %%
# Plot expression levels for selected 22q11.2 genes across subsets (Human CGE)
import warnings
warnings.filterwarnings('ignore')

import seaborn as sns
import matplotlib.pyplot as plt

# Prepare data for plotting
plot_data = []

for gene in selected_genes:
    # Extract expression for each subset directly from subset AnnData objects
    for subset_name, subset_adata in [('CCK_BC', CGE_CCK_BC), ('ISI2', CGE_ISI2), ('ISI3', CGE_ISI3)]:
        if subset_adata.n_obs > 0 and gene in subset_adata.var.index:
            gene_idx = subset_adata.var.index.get_loc(gene)
            
            # Extract expression values
            if hasattr(subset_adata.X, 'toarray'):
                expressions = subset_adata.X[:, gene_idx].toarray().flatten()
            elif hasattr(subset_adata.X, 'A1'):  # scipy sparse matrix
                expressions = subset_adata.X[:, gene_idx].A1
            else:
                expressions = np.array(subset_adata.X[:, gene_idx]).flatten()
            
            # Add to plot data
            for expr in expressions:
                plot_data.append({
                    'Gene': gene,
                    'Subset': subset_name,
                    'Expression': expr
                })

plot_df = pd.DataFrame(plot_data)

# Create the plot with 3 columns per row
n_cols = 3
n_rows = (len(selected_genes) + n_cols - 1) // n_cols  # Ceiling division
fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.5*n_cols, 6*n_rows))
axes = axes.flatten() if len(selected_genes) > 1 else [axes]

for idx, gene in enumerate(selected_genes):
    ax = axes[idx]
    gene_data = plot_df[plot_df['Gene'] == gene]
    
    # Create violin plot with box plot overlay
    sns.violinplot(data=gene_data, x='Subset', y='Expression', ax=ax, 
                   palette=['#1f77b4', '#ff7f0e', '#2ca02c'], inner='box')
    
    # Increase font sizes by 50%
    ax.set_title(f'{gene}', fontsize=21, fontweight='bold')
    ax.set_xlabel('Subset', fontsize=18)
    ax.set_ylabel('Expression Level (log2(CPM+1))', fontsize=18)
    ax.tick_params(labelsize=13.5)  # 50% increase from default ~9pt
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add statistics text
    for subset in ['CCK_BC', 'ISI2', 'ISI3']:
        subset_data = gene_data[gene_data['Subset'] == subset]['Expression']
        if len(subset_data) > 0:
            mean_expr = subset_data.mean()
            median_expr = subset_data.median()
            # Position text at top of plot
            ax.text(list(gene_data['Subset'].unique()).index(subset), 
                   ax.get_ylim()[1] * 0.95,
                   f'μ={mean_expr:.2f}\nmed={median_expr:.2f}',
                   ha='center', va='top', fontsize=12, 
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Hide unused subplots
for idx in range(len(selected_genes), len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.suptitle('22q11.2 Genes Expression Across CGE Subsets (Human)', 
             fontsize=24, fontweight='bold', y=1.02)
plt.show()

# Print summary statistics
print("\n" + "=" * 80)
print("Summary Statistics for 22q11.2 Genes (Human)")
print("=" * 80)
for gene in selected_genes:
    print(f"\n{gene}:")
    gene_data = plot_df[plot_df['Gene'] == gene]
    for subset in ['CCK_BC', 'ISI2', 'ISI3']:
        subset_data = gene_data[gene_data['Subset'] == subset]['Expression']
        if len(subset_data) > 0:
            print(f"  {subset}:")
            print(f"    Mean: {subset_data.mean():.3f}")
            print(f"    Median: {subset_data.median():.3f}")
            print(f"    Std: {subset_data.std():.3f}")
            print(f"    Cells expressing (>0): {np.sum(subset_data > 0)}/{len(subset_data)} ({np.sum(subset_data > 0)/len(subset_data)*100:.1f}%)")
print("=" * 80)

# %%
CGE_CCK_BC.obs.columns

# %%
CGE_CCK_BC.obs['cluster_id'].value_counts()

# %%
CGE_ISI2.obs['cluster_id'].value_counts()

# %%
CGE_ISI3.obs['cluster_id'].value_counts()

# %%
Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/results/random/Centering/"
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
X22q_mouse_Bias = pd.read_csv(Bias_Save_Dir + "22q_small_del_bias_addP.csv", index_col=0)
DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)

# %%
CGE_ISI3_bias = []
for i in CGE_ISI3.obs['cluster_id'].values:
    CGE_ISI3_bias.append(X22q_mouse_Bias.loc[i, "EFFECT"])
CGE_ISI2_bias = []
for i in CGE_ISI2.obs['cluster_id'].values:
    CGE_ISI2_bias.append(X22q_mouse_Bias.loc[i, "EFFECT"])
CGE_CCK_BC_bias = []
for i in CGE_CCK_BC.obs['cluster_id'].values:
    CGE_CCK_BC_bias.append(X22q_mouse_Bias.loc[i, "EFFECT"])

# Prepare for annotation
CGE_ISI3_bias = np.array(CGE_ISI3_bias)
CGE_ISI2_bias = np.array(CGE_ISI2_bias)
CGE_CCK_BC_bias = np.array(CGE_CCK_BC_bias)

plt.figure(figsize=(8, 6))
n1, bins1, patches1 = plt.hist(CGE_ISI3_bias, bins=100, alpha=0.6, label="ISI3")
n2, bins2, patches2 = plt.hist(CGE_ISI2_bias, bins=100, alpha=0.6, label="ISI2")
n3, bins3, patches3 = plt.hist(CGE_CCK_BC_bias, bins=100, alpha=0.6, label="CCK_BC")

plt.xlabel("Bias ('EFFECT')")
plt.ylabel("Number of Cells")
plt.title("Distribution of X22q_mouse_Bias EFFECT for CGE interneuron subtypes")
legend = plt.legend(title="Neuron type")

# DEBUGGED: Use get_legend_handles_labels() instead of non-existent legendHandles attribute
handles, labels = plt.gca().get_legend_handles_labels()
label_color_map = {}
for handle, label in zip(handles, labels):
    facecolor = None
    if hasattr(handle, "get_facecolor"):
        fc = handle.get_facecolor()
        # fc is usually (N,4) array of RGBA per handle, or tuple
        if hasattr(fc, "__len__") and len(fc) == 1:
            color = fc[0]
        else:
            color = fc
        facecolor = color
    label_color_map[label] = facecolor if facecolor is not None else "black"

# Annotate with mean for each neurontype, using matching legend color
means = [
    (CGE_ISI3_bias.mean(), "ISI3"),
    (CGE_ISI2_bias.mean(), "ISI2"),
    (CGE_CCK_BC_bias.mean(), "CCK_BC")
]
for mean_val, label in means:
    color = label_color_map.get(label, "black")
    plt.axvline(mean_val, linestyle="--", linewidth=2, label=f"{label} mean", alpha=0.9, color=color)
    plt.text(
        mean_val,
        plt.ylim()[1]*0.85,
        f"{label} mean\n{mean_val:.2f}",
        rotation=90,
        va='center',
        ha='right' if mean_val > plt.xlim()[1]/2 else 'left',
        color=color,
        fontsize=9
    )

plt.show()

# %%

# %%
