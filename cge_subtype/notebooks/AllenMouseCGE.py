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
    os.chdir(f"{ProjDIR}/notebooks_rebuttal/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:
    print(f"Unexpected error: {e}")


import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection, multipletests

# %%
Cortex1_raw = "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-log2.h5ad"
Cortex2_raw = "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-2-log2.h5ad"
HPF_raw = "/mnt/data0/AllenMouseSC/abc_download_root/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-HPF-log2.h5ad"

# %%
import anndata as ad

# First, read in backed mode to inspect metadata without loading full data into memory
print("Reading metadata in backed mode...")
cortex1_backed = ad.read_h5ad(Cortex1_raw, backed='r')
cortex2_backed = ad.read_h5ad(Cortex2_raw, backed='r')
hpf_backed = ad.read_h5ad(HPF_raw, backed='r')

print(f"\nCortex1 shape: {cortex1_backed.shape} (cells x genes)")
print(f"Cortex1 obs columns: {list(cortex1_backed.obs.columns)}")
print(f"\nCortex2 shape: {cortex2_backed.shape} (cells x genes)")
print(f"Cortex2 obs columns: {list(cortex2_backed.obs.columns)}")
print(f"\nHPF shape: {hpf_backed.shape} (cells x genes)")
print(f"HPF obs columns: {list(hpf_backed.obs.columns)}")

# %%
hpf_backed.var["gene_symbol"]

# %%
#CellMeta = pd.read_csv("/mnt/data0/AllenMouseSC/abc_download_root/metadata/WMB-10X/20230830/cell_metadata.csv")
CellMeta = pd.read_csv("/home/jw3514/Work/NeurSim/TransEphys/dat/expression/reference_atlas/cell_metadata.csv")
CellMeta_subset = CellMeta[CellMeta["feature_matrix_label"].isin(["WMB-10Xv3-Isocortex-1", "WMB-10Xv3-Isocortex-2", "WMB-10Xv3-HPF"])]
CellMeta_subset.set_index("cell_label", inplace=True)

# %%
CellMeta_CGE = CellMeta_subset[CellMeta_subset["class"] == "06 CTX-CGE GABA"]
CellMeta_CGE.shape


# %%
cortex2_backed.obs

# %%
# Get the cell labels we want to keep from CellMeta_CGE
cells_to_keep = set(CellMeta_CGE.index.values)
print(f"Total CGE cells to extract: {len(cells_to_keep)}")

# Function to filter and load cells matching CellMeta_CGE
def filter_by_cell_labels(adata_backed, cell_labels_set, library_name):
    """Filter h5ad file to only include cells in cell_labels_set and set gene_symbol as index"""
    # Match cell_label (index) in h5ad to cell_label in CellMeta
    mask = adata_backed.obs.index.isin(cell_labels_set)
    n_matched = mask.sum()
    print(f"  {library_name}: Found {n_matched} matching cells (out of {adata_backed.n_obs})")
    
    if n_matched == 0:
        print(f"    Warning: No cells matched for {library_name}")
        return None
    
    # Subset and load into memory (use to_memory() for backed objects)
    adata_filtered = adata_backed[mask].to_memory()
    
    # Keep Ensembl ID as index during filtering to avoid duplicate issues during merge
    # We'll set gene_symbol as index after merging
    
    return adata_filtered

# Filter each dataset
print("\nFiltering Cortex1...")
cortex1_filtered = filter_by_cell_labels(cortex1_backed, cells_to_keep, "Cortex1")

print("Filtering Cortex2...")
cortex2_filtered = filter_by_cell_labels(cortex2_backed, cells_to_keep, "Cortex2")

print("Filtering HPF...")
hpf_filtered = filter_by_cell_labels(hpf_backed, cells_to_keep, "HPF")

# Combine all filtered datasets
datasets_to_merge = []
if cortex1_filtered is not None:
    datasets_to_merge.append(cortex1_filtered)
if cortex2_filtered is not None:
    datasets_to_merge.append(cortex2_filtered)
if hpf_filtered is not None:
    datasets_to_merge.append(hpf_filtered)

if len(datasets_to_merge) > 0:
    print(f"\nMerging {len(datasets_to_merge)} datasets...")
    CGE_IN_CPM = ad.concat(datasets_to_merge, join='inner', index_unique=None)
    
    # Set gene_symbol as the gene index after merging
    # Get gene_symbol from one of the original datasets
    if 'gene_symbol' in datasets_to_merge[0].var.columns:
        # Map Ensembl IDs to gene_symbol
        gene_symbol_map = datasets_to_merge[0].var['gene_symbol'].to_dict()
        # Add gene_symbol to var if not already present
        if 'gene_symbol' not in CGE_IN_CPM.var.columns:
            CGE_IN_CPM.var['gene_symbol'] = CGE_IN_CPM.var.index.map(gene_symbol_map)
        
        # Check for duplicates and make unique if needed
        if CGE_IN_CPM.var['gene_symbol'].duplicated().any():
            print(f"    Warning: Found duplicate gene symbols. Making index unique...")
            n_duplicates = CGE_IN_CPM.var['gene_symbol'].duplicated().sum()
            print(f"    Found {n_duplicates} duplicate gene symbols")
            # Make duplicates unique by appending number suffix using groupby
            counts = CGE_IN_CPM.var.groupby('gene_symbol').cumcount()
            CGE_IN_CPM.var['gene_symbol_unique'] = CGE_IN_CPM.var['gene_symbol']
            # Add suffix only to duplicates (where count > 0)
            mask = counts > 0
            CGE_IN_CPM.var.loc[mask, 'gene_symbol_unique'] = (
                CGE_IN_CPM.var.loc[mask, 'gene_symbol'] + '_' + counts[mask].astype(str)
            )
            CGE_IN_CPM.var.set_index('gene_symbol_unique', inplace=True)
        else:
            # No duplicates, safe to set as index
            CGE_IN_CPM.var.set_index('gene_symbol', inplace=True)
        
        print(f"    Set gene_symbol as gene index")
    
    print(f"Final merged dataset CGE_IN_CPM: {CGE_IN_CPM.shape} (cells x genes)")
    print(f"Total CGE cells loaded: {CGE_IN_CPM.n_obs}")
else:
    print("Error: No cells were matched. Check that the AnnData index (cell_label) matches CellMeta_CGE index.")

# %%
CGE_IN_CPM.var


# %%
"Cnr1" in CGE_IN_CPM.var.index

# %%
Makders = ["Vip", "Cck", "Sncg", "Chrm2", "Calb2", "Cnr1"]

# %%
import matplotlib.pyplot as plt
import numpy as np

# Fix typo in variable name
Markers = Makders if 'Makders' in locals() else ["Vip", "Cck", "Sncg", "Chrm2", "Calb2"]

# Create figure with subplots
n_markers = len(Markers)
n_cols = 3
n_rows = (n_markers + n_cols - 1) // n_cols  # Ceiling division

fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
axes = axes.flatten() if n_markers > 1 else [axes]

for idx, marker in enumerate(Markers):
    ax = axes[idx]
    
    # Check if marker exists in the index
    if marker in CGE_IN_CPM.var.index:
        # Get expression values for this marker
        marker_idx = CGE_IN_CPM.var.index.get_loc(marker)
        
        # Extract expression values (handle sparse matrix)
        if hasattr(CGE_IN_CPM.X, 'toarray'):
            expression = CGE_IN_CPM.X[:, marker_idx].toarray().flatten()
        else:
            expression = CGE_IN_CPM.X[:, marker_idx].flatten()
        
        # Plot histogram
        ax.hist(expression, bins=50, edgecolor='black', alpha=0.7)
        ax.set_title(f'{marker} Expression Distribution', fontsize=12, fontweight='bold')
        ax.set_xlabel('Expression Level (log2)', fontsize=10)
        ax.set_ylabel('Number of Cells', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        mean_expr = np.mean(expression)
        median_expr = np.median(expression)
        ax.axvline(mean_expr, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_expr:.2f}')
        ax.axvline(median_expr, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_expr:.2f}')
        ax.legend(fontsize=8)
        
        print(f"{marker}: Mean={mean_expr:.2f}, Median={median_expr:.2f}, Max={np.max(expression):.2f}, Non-zero cells={np.sum(expression > 0)}/{len(expression)}")
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
# CPM_Cut_off = 2
# CCK_BC: Cck > CPM_Cut_off; Sncg > CPM_Cut_off; Calb2 < CPM_Cut_off
# ISI2: Vip > CPM_Cut_off; Cck < CPM_Cut_off; Sncg < CPM_Cut_off; Calb2 < CPM_Cut_off
# ISI3: Vip > CPM_Cut_off; Cck < CPM_Cut_off; Sncg < CPM_Cut_off; Calb2 > CPM_Cut_off




# %%
# Extract expression values for each marker
def get_expression(adata, gene_name):
    """Extract expression values for a gene"""
    if gene_name not in adata.var.index:
        print(f"Warning: {gene_name} not found in var.index")
        return None
    marker_idx = adata.var.index.get_loc(gene_name)
    if hasattr(adata.X, 'toarray'):
        return adata.X[:, marker_idx].toarray().flatten()
    else:
        return adata.X[:, marker_idx].flatten()

# Get expression for all markers
Vip_expr = get_expression(CGE_IN_CPM, "Vip")
Cck_expr = get_expression(CGE_IN_CPM, "Cck")
Sncg_expr = get_expression(CGE_IN_CPM, "Sncg")
Calb2_expr = get_expression(CGE_IN_CPM, "Calb2")

# Define cutoff
CPM_Cut_off = 2

# Create boolean masks for each subset
# CCK_BC: Cck > CPM_Cut_off; Sncg > CPM_Cut_off; Calb2 < CPM_Cut_off
CCK_BC_mask = (Cck_expr > CPM_Cut_off) & (Sncg_expr > CPM_Cut_off) & (Calb2_expr < CPM_Cut_off)

# ISI2: Vip > CPM_Cut_off; Cck < CPM_Cut_off; Sncg < CPM_Cut_off; Calb2 < CPM_Cut_off
ISI2_mask = (Vip_expr > CPM_Cut_off) & (Cck_expr < CPM_Cut_off) & (Sncg_expr < CPM_Cut_off) & (Calb2_expr < CPM_Cut_off)

# ISI3: Vip > CPM_Cut_off; Cck < CPM_Cut_off; Sncg < CPM_Cut_off; Calb2 > CPM_Cut_off
ISI3_mask = (Vip_expr > CPM_Cut_off) & (Cck_expr < CPM_Cut_off) & (Sncg_expr < CPM_Cut_off) & (Calb2_expr > CPM_Cut_off)

# Make subsets mutually exclusive (prioritize CCK_BC, then ISI2, then ISI3)
# Start with ISI2 and ISI3 (they don't overlap with CCK_BC)
ISI2_mask_exclusive = ISI2_mask & ~CCK_BC_mask
ISI3_mask_exclusive = ISI3_mask & ~CCK_BC_mask

# Create subset AnnData objects
CGE_CCK_BC = CGE_IN_CPM[CCK_BC_mask].copy()
CGE_ISI2 = CGE_IN_CPM[ISI2_mask_exclusive].copy()
CGE_ISI3 = CGE_IN_CPM[ISI3_mask_exclusive].copy()

# Print summary
print("=" * 60)
print("CGE Subset Classification Results")
print("=" * 60)
print(f"CPM Cutoff: {CPM_Cut_off}")
print()
print(f"CCK_BC: {CCK_BC_mask.sum()} cells")
print(f"  Criteria: Cck > {CPM_Cut_off}, Sncg > {CPM_Cut_off}, Calb2 < {CPM_Cut_off}")
print()
print(f"ISI2: {ISI2_mask_exclusive.sum()} cells")
print(f"  Criteria: Vip > {CPM_Cut_off}, Cck < {CPM_Cut_off}, Sncg < {CPM_Cut_off}, Calb2 < {CPM_Cut_off}")
print()
print(f"ISI3: {ISI3_mask_exclusive.sum()} cells")
print(f"  Criteria: Vip > {CPM_Cut_off}, Cck < {CPM_Cut_off}, Sncg < {CPM_Cut_off}, Calb2 > {CPM_Cut_off}")
print()
print(f"Total classified: {CCK_BC_mask.sum() + ISI2_mask_exclusive.sum() + ISI3_mask_exclusive.sum()} cells")
print(f"Unclassified: {CGE_IN_CPM.n_obs - (CCK_BC_mask.sum() + ISI2_mask_exclusive.sum() + ISI3_mask_exclusive.sum())} cells")
print("=" * 60)

# Store masks in a dictionary for easy access
CGE_subset_masks = {
    'CCK_BC': CCK_BC_mask,
    'ISI2': ISI2_mask_exclusive,
    'ISI3': ISI3_mask_exclusive
}

# %%
CGE_CCK_BC

# %%
import matplotlib.pyplot as plt

# Get VIP expression vector for all CGE interneurons
vip_expr_values = CGE_CCK_BC.obs_vector("Vip")

plt.figure(figsize=(8, 5))
plt.hist(vip_expr_values, bins=50, color='skyblue', edgecolor='k')
plt.xlabel('VIP expression (CPM)')
plt.ylabel('Number of cells')
plt.title('Histogram of VIP expression in CGE interneurons')
plt.tight_layout()
plt.show()

# %%
# Report ROI distribution for each subset
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
    roi_counts = subset_adata.obs['anatomical_division_label'].value_counts()
    total = subset_adata.n_obs
    
    for roi, count in roi_counts.items():
        roi_summary.append({
            'Subset': subset_name,
            'ROI': roi,
            'Count': count,
            'Percentage': (count / total * 100) if total > 0 else 0,
            'Total': total
        })

# Get unclassified cells
unclassified_mask = ~(CCK_BC_mask | ISI2_mask_exclusive | ISI3_mask_exclusive)
unclassified_adata = CGE_IN_CPM[unclassified_mask]
unclassified_roi_counts = unclassified_adata.obs['anatomical_division_label'].value_counts()
unclassified_total = unclassified_adata.n_obs

for roi, count in unclassified_roi_counts.items():
    roi_summary.append({
        'Subset': 'Unclassified',
        'ROI': roi,
        'Count': count,
        'Percentage': (count / unclassified_total * 100) if unclassified_total > 0 else 0,
        'Total': unclassified_total
    })

# Create DataFrame and display
roi_df = pd.DataFrame(roi_summary)

# Pivot table for better visualization
roi_pivot = roi_df.pivot_table(
    index='Subset',
    columns='ROI',
    values='Count',
    fill_value=0
)

# Add total column
roi_pivot['Total'] = roi_pivot.sum(axis=1)

print("=" * 80)
print("ROI Distribution by CGE Subset")
print("=" * 80)
print()
print(roi_pivot)
print()

# Also show percentages
roi_pct = roi_df.pivot_table(
    index='Subset',
    columns='ROI',
    values='Percentage',
    fill_value=0
)
roi_pct['Total_Cells'] = roi_pivot['Total']

print("=" * 80)
print("ROI Distribution (Percentages)")
print("=" * 80)
print()
for subset in roi_pct.index:
    print(f"\n{subset}:")
    total = roi_pct.loc[subset, 'Total_Cells']
    for roi in roi_pct.columns:
        if roi != 'Total_Cells':
            count = roi_pivot.loc[subset, roi]
            pct = roi_pct.loc[subset, roi]
            print(f"  {roi}: {count} cells ({pct:.1f}%)")
    print(f"  Total: {int(total)} cells")

print()
print("=" * 80)
print("Overall ROI Distribution in CGE_IN_CPM")
print("=" * 80)
overall_roi = CGE_IN_CPM.obs['anatomical_division_label'].value_counts()
for roi, count in overall_roi.items():
    pct = count / CGE_IN_CPM.n_obs * 100
    print(f"{roi}: {count} cells ({pct:.1f}%)")
print(f"Total: {CGE_IN_CPM.n_obs} cells")
print("=" * 80)

# %%
# Define important 22q11.2 genes in mouse genome
# These are orthologs of human 22q11.2 genes
important_22q11_genes = [
    'Comt',      # Catechol-O-methyltransferase - critical for dopamine metabolism
    'Dgcr8',     # DiGeorge syndrome critical region gene 8 - microRNA processing
    'Prodh',     # Proline dehydrogenase - associated with schizophrenia
    'Tbx1',      # T-box transcription factor 1 - important for development
    'Gp1bb',     # Glycoprotein Ib beta - platelet function
    'Ufd1l',     # Ubiquitin fusion degradation 1 like
    'Hira',      # Histone cell cycle regulator
    'Txnrd2',    # Thioredoxin reductase 2
    'Slc25a1',   # Solute carrier family 25 member 1
    'Dgcr2',     # DiGeorge syndrome critical region gene 2
    'Dgcr14',    # DiGeorge syndrome critical region gene 14
    'Sept5',     # Septin 5
    'Snap29',    # Synaptosome associated protein 29
    'Zdhhc8',    # Zinc finger DHHC-type containing 8
    'Rtn4r',     # Reticulon 4 receptor
    'Arvcf',     # Armadillo repeat containing, catenin delta
    'Gscl',      # Goosecoid like
    'Cldn5',     # Claudin 5
    'Hic2',      # HIC ZBTB transcriptional repressor 2
    'Uba1',      # Ubiquitin like modifier activating enzyme 1
    'T10',       # T10 protein
    'P2rx6'      # Purinergic receptor P2X 6
]

# Check which genes are available in the dataset
available_genes = []
missing_genes = []

for gene in important_22q11_genes:
    if gene in CGE_IN_CPM.var.index:
        available_genes.append(gene)
    else:
        missing_genes.append(gene)

print("=" * 80)
print("22q11.2 Genes Availability Check")
print("=" * 80)
print(f"\nAvailable genes ({len(available_genes)}): {available_genes}")
print(f"\nMissing genes ({len(missing_genes)}): {missing_genes}")

# Select top 5 most important genes that are available
# Priority order: Comt, Dgcr8, Prodh, Tbx1, then others
priority_order = ['Comt', 'Dgcr8', 'Prodh', 'Tbx1', 'Gp1bb', 'Ufd1l', 'Hira', 
                  'Txnrd2', 'Slc25a1', 'Dgcr2', 'Zdhhc8', 'Rtn4r', 'Arvcf']

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
# Plot expression levels for selected 22q11.2 genes across subsets
import warnings
warnings.filterwarnings('ignore')

import seaborn as sns

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
            else:
                expressions = subset_adata.X[:, gene_idx].flatten()
            
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
    ax.set_ylabel('Expression Level (log2)', fontsize=18)
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
plt.suptitle('22q11.2 Genes Expression Across CGE Subsets', 
             fontsize=24, fontweight='bold', y=1.02)
plt.show()

# Print summary statistics
print("\n" + "=" * 80)
print("Summary Statistics for 22q11.2 Genes")
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

# %%
