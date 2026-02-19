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


try:
    os.chdir(f"{ProjDIR}/notebooks/")
    print(f"Current working directory: {os.getcwd()}")
except FileNotFoundError as e:
    print(f"Error: Could not change directory - {e}")
except Exception as e:  
    print(f"Unexpected error: {e}")    


import yaml
with open(ProjDIR + '/config/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# %%
expression_matrix = config['analysis_types']['Centering']
print(expression_matrix)
HCT_Spec_MAT = pd.read_csv(ProjDIR + expression_matrix, index_col=0)
HCT_Spec_MAT.columns = HCT_Spec_MAT.columns.astype(int)

match_table = pd.read_csv("../dat/Variable_2_Match_master_table_pct.csv", index_col=0)
match_table.head()

# %%
import matplotlib.font_manager as fm
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
fm.fontManager.addfont(font_path)  # Only if you're adding a new font file
fm._load_fontmanager(try_read_cache=False)
#plt.style.use('seaborn-v0_8-paper')
plt.style.use('seaborn-v0_8-whitegrid')

# %%
Bias_Save_Dir = "/home/jw3514/Work/CellType_Psy/CellTypeBias_VIP/results/main_results/matched_WB_mean_phastCons_n_CDS_bases_Best1000/Centering/"

ASD_All_Bias = pd.read_csv(Bias_Save_Dir + "ASD_All_bias_addP.csv", index_col=0)
SCZ_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_bias_addP.csv", index_col=0)
HighIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_bias_addP.csv", index_col=0)
LowIQ_ASD_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_bias_addP.csv", index_col=0)
X22q_Bias = pd.read_csv(Bias_Save_Dir + "22q_del_bias_addP.csv", index_col=0)
try:
    DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_61_bias_addP.csv", index_col=0)
except:
    DDD_Bias = pd.read_csv(Bias_Save_Dir + "DDD_bias_addP.csv", index_col=0)

DDD_297_Bias = pd.read_csv(Bias_Save_Dir + "DDD_297_bias_addP.csv", index_col=0)

try:
    SCZ_100_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_100_bias_addP.csv", index_col=0)
    SCZ_200_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_200_bias_addP.csv", index_col=0)
    SCZ_500_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_500_bias_addP.csv", index_col=0)
    ASD_HIQ_100_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_100_bias_addP.csv", index_col=0)
    ASD_LIQ_100_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_100_bias_addP.csv", index_col=0)
except:
    pass

# VNR_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Pos_bias_addP.csv", index_col=0)
# VNR_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_VNR_Neg_bias_addP.csv", index_col=0)
# EDU_Pos_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Pos_bias_addP.csv", index_col=0)
# EDU_Neg_Bias = pd.read_csv(Bias_Save_Dir + "UKBB_EDU_Neg_bias_addP.csv", index_col=0)


# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ### Cell type bias and rank for each disorders

# %% [markdown]
# #### Plot Pvalue BoxPlot (Figure S4)

# %%
fdr_cut = 0.1

# Create a 2x3 subplot grid
fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120, facecolor='none')
fig.patch.set_alpha(0.0)
plt.style.use('seaborn-v0_8-whitegrid')

# Flatten axes array for easier indexing
axes = axes.flatten()

# Plot each dataset on its own subplot
SuperClusterBias_BoxPlot(ASD_All_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[0])
SuperClusterBias_BoxPlot(HighIQ_ASD_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[1])
SuperClusterBias_BoxPlot(LowIQ_ASD_Bias, "ASD with ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[2])
#SuperClusterBias_BoxPlot(X22q_Bias, "22q11.2", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])
SuperClusterBias_BoxPlot(SCZ_Bias, "SCZ", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[4])
SuperClusterBias_BoxPlot(DDD_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[5])
SuperClusterBias_BoxPlot(DDD_297_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])

plt.tight_layout()
plt.show()

# %%
# fdr_cut = 0.1

# # Arrange panels for a nice logical flow: left->right, top to bottom
# # [0] ASD_HIQ_100, [1] ASD_LIQ_100, [2] SCZ_100, [3] SCZ_200, [4] SCZ_500, [5] DDD_297

# # Create a 2x3 subplot grid
# fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=120, facecolor='none')
# fig.patch.set_alpha(0.0)
# plt.style.use('seaborn-v0_8-whitegrid')
# axes = axes.flatten()

# # Fill panels in a logical order, one dataset per panel
# SuperClusterBias_BoxPlot(ASD_HIQ_100_Bias, "ASD (High IQ Top 100)", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[0])
# SuperClusterBias_BoxPlot(ASD_LIQ_100_Bias, "ASD (Low IQ Top 100)", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[1])
# SuperClusterBias_BoxPlot(SCZ_100_Bias, "SCZ Top 100", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[2])
# SuperClusterBias_BoxPlot(SCZ_200_Bias, "SCZ Top 200", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])
# SuperClusterBias_BoxPlot(SCZ_500_Bias, "SCZ Top 500", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[4])
# SuperClusterBias_BoxPlot(DDD_297_Bias, "DD/ID Top 297", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[5])

# # Make spacing nice and add supertitle
# plt.tight_layout()
# plt.subplots_adjust(top=0.88)
# fig.suptitle("Cell-type bias (-logP) distributions for select matched gene sets", fontsize=18, y=1.02)
# plt.show()

# %%
fdr_cut = 0.05

# Create a 2x3 subplot grid
fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120, facecolor='none')
fig.patch.set_alpha(0.0)
plt.style.use('seaborn-v0_8-whitegrid')

# Flatten axes array for easier indexing
axes = axes.flatten()

# Plot each dataset on its own subplot
SuperClusterBias_BoxPlot(ASD_All_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[0])
SuperClusterBias_BoxPlot(HighIQ_ASD_Bias, "ASD", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[1])
SuperClusterBias_BoxPlot(LowIQ_ASD_Bias, "ASD with ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[2])
#SuperClusterBias_BoxPlot(X22q_Bias, "22q11.2", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])
SuperClusterBias_BoxPlot(SCZ_Bias, "SCZ", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[4])
SuperClusterBias_BoxPlot(DDD_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[5])
SuperClusterBias_BoxPlot(DDD_297_Bias, "DD/ID", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])

plt.tight_layout()
plt.show()

# %%
# fdr_cut = 0.05

# # Arrange panels for a nice logical flow: left->right, top to bottom
# # [0] ASD_HIQ_100, [1] ASD_LIQ_100, [2] SCZ_100, [3] SCZ_200, [4] SCZ_500, [5] DDD_297

# # Create a 2x3 subplot grid
# fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=120, facecolor='none')
# fig.patch.set_alpha(0.0)
# plt.style.use('seaborn-v0_8-whitegrid')
# axes = axes.flatten()

# # Fill panels in a logical order, one dataset per panel
# SuperClusterBias_BoxPlot(ASD_HIQ_100_Bias, "ASD (High IQ Top 100)", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[0])
# SuperClusterBias_BoxPlot(ASD_LIQ_100_Bias, "ASD (Low IQ Top 100)", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[1])
# SuperClusterBias_BoxPlot(SCZ_100_Bias, "SCZ Top 100", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[2])
# SuperClusterBias_BoxPlot(SCZ_200_Bias, "SCZ Top 200", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[3])
# SuperClusterBias_BoxPlot(SCZ_500_Bias, "SCZ Top 500", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[4])
# SuperClusterBias_BoxPlot(DDD_297_Bias, "DD/ID Top 297", NeuroOnly=False, sortby="mean", EffectCol="-logP", fdr_cut=fdr_cut, ax=axes[5])

# # Make spacing nice and add supertitle
# plt.tight_layout()
# plt.subplots_adjust(top=0.88)
# fig.suptitle("Cell-type bias (-logP) distributions for select matched gene sets", fontsize=18, y=1.02)
# plt.show()

# %% [markdown]
# #### [END] Plot Pvalue BoxPlot

# %%
#### Compare Real vs Matched Genes Properties Across Multiple Gene Sets

# Define gene sets to analyze
gene_sets_config = {
    'SCZ': {'file': 'SCZ_random_geneweights.csv', 'display_name': 'SCZ'},
    'ASD_All': {'file': 'ASD_All_random_geneweights.csv', 'display_name': 'ASD'},
    'ASD_HIQ': {'file': 'ASD_HIQ_random_geneweights.csv', 'display_name': 'ASD (High IQ)'},
    'ASD_LIQ': {'file': 'ASD_LIQ_random_geneweights.csv', 'display_name': 'ASD (Low IQ)'},
    'X22q': {'file': '22q_del_random_geneweights.csv', 'display_name': '22q11.2'},
    'DDD': {'file': ['DDD_61_random_geneweights.csv', 'DDD_random_geneweights.csv'], 'display_name': 'DD/ID'},
}

# Properties to compare - using percentile versions for better comparison
properties_to_compare = ['WB_pct', 'mean_phastCons_pct', 'n_CDS_bases_pct']
n_matches_to_use = 10  # Number of matched sets to analyze

# %%
match_table


# %%
# Function to analyze a single gene set
def analyze_gene_set(gene_set_name, file_name, display_name, match_table, properties, n_matches=1000):
    """
    Analyze real vs matched genes for a single gene set.
    
    Args:
        file_name: Can be a string or list of strings (to try multiple file names)
    
    Returns:
        dict with 'real_genes', 'real_stats', 'matched_stats', 'matched_properties'
    """
    # Handle multiple possible file names
    if isinstance(file_name, list):
        file_names_to_try = file_name
    else:
        file_names_to_try = [file_name]
    
    null_genes_df = None
    used_file = None
    
    for fn in file_names_to_try:
        try:
            null_genes_df = pd.read_csv(f"{Bias_Save_Dir}/null_weights/{fn}", index_col=0)
            used_file = fn
            break
        except FileNotFoundError:
            continue
    
    if null_genes_df is None:
        raise FileNotFoundError(f"Could not find any of the files: {file_names_to_try}")
    
    try:
        
        # Get real genes (from index)
        real_genes_all = null_genes_df.index.values
        
        # Filter to only genes that exist in match_table using intersection
        # This handles type mismatches automatically
        real_genes_index = null_genes_df.index.intersection(match_table.index)
        real_genes = real_genes_index.values
        
        missing_genes = null_genes_df.index.difference(match_table.index).values
        
        if len(real_genes) == 0:
            raise ValueError(f"No genes from {gene_set_name} found in match_table")
        
        if len(missing_genes) > 0:
            print(f"  Warning: {len(missing_genes)} genes not found in match_table (using {len(real_genes)}/{len(real_genes_all)})")
            if len(missing_genes) <= 5:
                print(f"    Missing genes: {missing_genes.tolist()}")
            else:
                print(f"    Missing genes (first 5): {missing_genes[:5].tolist()}...")
        
        # Get matched gene columns (skip first column which is GeneWeight)
        matched_genes_cols = null_genes_df.columns[1:].values[:n_matches]
        
        # Extract properties for real genes
        available_props = [prop for prop in properties if prop in match_table.columns]
        real_stats = {}
        real_properties = {}
        
        for prop in available_props:
            # Use only genes that exist in match_table (use index for safer lookup)
            values = match_table.loc[real_genes_index, prop].dropna().values
            real_properties[prop] = values
            if len(values) > 0:
                real_stats[prop] = {
                    'mean': values.mean(),
                    'median': np.median(values),
                    'std': values.std(),
                    'min': values.min(),
                    'max': values.max(),
                    'n': len(values)
                }
        
        # Extract properties for matched genes
        matched_properties = {prop: [] for prop in available_props}
        matched_means_per_set = {prop: [] for prop in available_props}
        
        for match_col in matched_genes_cols:
            matched_genes = null_genes_df[match_col].dropna().values
            
            for prop in available_props:
                matched_genes_in_table = [g for g in matched_genes if g in match_table.index]
                if len(matched_genes_in_table) > 0:
                    prop_values = match_table.loc[matched_genes_in_table, prop].dropna().values
                    if len(prop_values) > 0:
                        matched_properties[prop].append(prop_values)
                        matched_means_per_set[prop].append(prop_values.mean())
        
        # Calculate summary statistics for matched genes
        matched_stats = {}
        for prop in available_props:
            all_matched_values = np.concatenate(matched_properties[prop]) if len(matched_properties[prop]) > 0 else np.array([])
            if len(all_matched_values) > 0:
                matched_stats[prop] = {
                    'mean': all_matched_values.mean(),
                    'median': np.median(all_matched_values),
                    'std': all_matched_values.std(),
                    'min': all_matched_values.min(),
                    'max': all_matched_values.max(),
                    'n': len(all_matched_values),
                    'mean_per_set': matched_means_per_set[prop]
                }
        
        return {
            'name': gene_set_name,
            'display_name': display_name,
            'real_genes': real_genes,
            'real_stats': real_stats,
            'real_properties': real_properties,
            'matched_stats': matched_stats,
            'matched_properties': matched_properties,
            'n_real_genes': len(real_genes),
            'n_matched_sets': len(matched_genes_cols),
            'available_props': available_props,
            'file_used': used_file
        }
    except FileNotFoundError as e:
        print(f"Warning: File not found for {gene_set_name}: {e}")
        return None
    except Exception as e:
        print(f"Error analyzing {gene_set_name}: {e}")
        return None

# Analyze all gene sets
print("Loading and analyzing gene sets...")
gene_set_results = {}

for gene_set_name, config in gene_sets_config.items():
    print(f"\nAnalyzing {config['display_name']}...")
    result = analyze_gene_set(
        gene_set_name, 
        config['file'], 
        config['display_name'],
        match_table,
        properties_to_compare,
        n_matches_to_use
    )
    if result is not None:
        gene_set_results[gene_set_name] = result
        file_info = f" (file: {result['file_used']})" if 'file_used' in result else ""
        print(f"  ✓ {result['n_real_genes']} real genes, {result['n_matched_sets']} matched sets{file_info}")

print(f"\n✓ Successfully analyzed {len(gene_set_results)} gene sets")


# %%
#### Summary Statistics Table

# Create summary table comparing real vs matched for all gene sets
summary_data = []

for gene_set_name, result in gene_set_results.items():
    for prop in result['available_props']:
        if prop in result['real_stats'] and prop in result['matched_stats']:
            real_mean = result['real_stats'][prop]['mean']
            matched_mean = result['matched_stats'][prop]['mean']
            mean_diff = matched_mean - real_mean
            mean_diff_pct = (mean_diff / real_mean * 100) if real_mean != 0 else np.nan
            
            summary_data.append({
                'Gene Set': result['display_name'],
                'Property': prop,
                'Real Mean': real_mean,
                'Matched Mean': matched_mean,
                'Difference': mean_diff,
                'Difference %': mean_diff_pct,
                'N Real Genes': result['n_real_genes'],
                'N Matched Sets': result['n_matched_sets']
            })

summary_df = pd.DataFrame(summary_data)
print("\n" + "="*80)
print("SUMMARY: Real vs Matched Gene Properties")
print("="*80)
print(summary_df.to_string(index=False))
print("\n")


# %%
#### Comparison Plots: Real vs Matched Across All Gene Sets

import os
from scipy.stats import mannwhitneyu

# Ensure output directory exists
fig_outdir = "../results/figures"
os.makedirs(fig_outdir, exist_ok=True)

# Get all available properties (should be same for all gene sets)
if len(gene_set_results) > 0:
    available_props = list(gene_set_results.values())[0]['available_props']
    
    # Create comparison plots for each property
    for prop in available_props:
        n_sets = len(gene_set_results)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120)
        fig.patch.set_alpha(0.0)
        axes = axes.flatten()
        
        plot_idx = 0
        for gene_set_name, result in gene_set_results.items():
            if plot_idx >= len(axes):
                break
                
            ax = axes[plot_idx]
            
            if prop in result['real_properties'] and prop in result['matched_properties']:
                real_values = result['real_properties'][prop]
                all_matched_values = np.concatenate(result['matched_properties'][prop]) if len(result['matched_properties'][prop]) > 0 else np.array([])
                
                if len(real_values) > 0 and len(all_matched_values) > 0:
                    # Create violin plot
                    data_to_plot = [real_values, all_matched_values]
                    labels = ['Real\nGenes', 'Matched\nGenes']
                    
                    parts = ax.violinplot(data_to_plot, positions=[1, 2], widths=0.6, showmeans=True, showmedians=True)
                    
                    # Customize violin plot
                    for pc in parts['bodies']:
                        pc.set_facecolor('lightblue')
                        pc.set_alpha(0.7)
                    
                    # Add box plot overlay
                    bp = ax.boxplot(data_to_plot, positions=[1, 2], widths=0.3, patch_artist=False,
                                   boxprops=dict(color='black', linewidth=1.5),
                                   medianprops=dict(color='red', linewidth=2),
                                   showmeans=False)
                    
                    ax.set_xticks([1, 2])
                    ax.set_xticklabels(labels, fontsize=12)
                    ax.set_ylabel(prop, fontsize=13, fontweight='bold')
                    ax.set_title(f'{result["display_name"]}\n(n={result["n_real_genes"]})', fontsize=14, fontweight='bold')
                    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
                    
                    # Calculate p-value using Mann-Whitney U test
                    try:
                        stat, pval = mannwhitneyu(all_matched_values, real_values, alternative='two-sided')
                        # Format p-value for display (scientific notation for small values)
                        if pval < 0.001:
                            pval_str = f'P = {pval:.2e}'.replace('e-0', ' x 10^-').replace('e-', ' x 10^-')
                        elif pval < 0.01:
                            pval_str = f'P = {pval:.2e}'.replace('e-0', ' x 10^-').replace('e-', ' x 10^-')
                        else:
                            pval_str = f'P = {pval:.3f}'
                        
                        # Color code p-value: red if significant, black if not
                        pval_color = 'red' if pval < 0.05 else 'black'
                    except Exception as e:
                        pval_str = 'P = N/A'
                        pval_color = 'black'
                        pval = 1.0
                    
                    # Get y-axis limits to position the significance bar
                    y_min, y_max = ax.get_ylim()
                    y_range = y_max - y_min
                    # Position bar near the top of the plot
                    bar_y = y_max - 0.05 * y_range
                    
                    # Draw horizontal bar connecting the two violins for p-value annotation
                    bar_x_start = 1
                    bar_x_end = 2
                    bar_height = 0.01 * y_range  # Small vertical offset for the bar
                    
                    # Draw the horizontal line (bar)
                    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    # Draw vertical lines at the ends
                    ax.plot([bar_x_start, bar_x_start], [bar_y - bar_height, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    ax.plot([bar_x_end, bar_x_end], [bar_y - bar_height, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    
                    # Add p-value text above the bar
                    ax.text((bar_x_start + bar_x_end) / 2, bar_y + 0.02 * y_range, 
                           pval_str, ha='center', va='bottom', fontsize=11,
                           color=pval_color, fontweight='bold' if pval < 0.05 else 'normal',
                           zorder=10)
                    
                    # Add mean difference in text box (keep this for reference)
                    real_mean = real_values.mean()
                    matched_mean = all_matched_values.mean()
                    mean_diff = matched_mean - real_mean
                    
                    # Create text box with mean difference (smaller, bottom left)
                    textstr = f'Δ={mean_diff:.2f}'
                    ax.text(0.05, 0.05, textstr, 
                            transform=ax.transAxes, ha='left', fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
                            zorder=10)
            
            plot_idx += 1
        
        # Hide unused subplots
        for idx in range(plot_idx, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle(f'Real vs Matched: {prop}', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        fig_name = os.path.join(fig_outdir, f"compare_violin_real_vs_matched_{prop}.png")
        plt.savefig(fig_name, bbox_inches='tight')
        plt.show()
else:
    print("No gene sets were successfully analyzed.")


# %%
#### Distribution of Matched Set Means: Comparison Across Gene Sets

# Plot distribution of matched set means for each property
if len(gene_set_results) > 0:
    available_props = list(gene_set_results.values())[0]['available_props']
    
    for prop in available_props:
        n_sets = len(gene_set_results)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120)
        fig.patch.set_alpha(0.0)
        axes = axes.flatten()
        
        plot_idx = 0
        for gene_set_name, result in gene_set_results.items():
            if plot_idx >= len(axes):
                break
                
            ax = axes[plot_idx]
            
            if prop in result['real_stats'] and prop in result['matched_stats']:
                real_mean = result['real_stats'][prop]['mean']
                matched_means = result['matched_stats'][prop]['mean_per_set']
                
                if len(matched_means) > 0:
                    # Plot histogram of matched means
                    ax.hist(matched_means, bins=30, alpha=0.7, color='lightblue', 
                           edgecolor='black', label='Matched sets')
                    
                    # Add vertical line for real mean
                    ax.axvline(real_mean, color='red', linestyle='--', linewidth=2.5, 
                              label=f'Real mean ({real_mean:.2f})')
                    
                    # Calculate how many matched sets are within reasonable range
                    matched_mean_array = np.array(matched_means)
                    std_real = result['real_stats'][prop]['std']
                    within_1std = np.sum(np.abs(matched_mean_array - real_mean) <= std_real)
                    within_2std = np.sum(np.abs(matched_mean_array - real_mean) <= 2*std_real)
                    
                    ax.set_xlabel(f'{prop}', fontsize=11, fontweight='bold')
                    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
                    ax.set_title(f'{result["display_name"]}\n(n={len(matched_means)} sets)', 
                               fontsize=12, fontweight='bold')
                    ax.legend(fontsize=9)
                    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
                    
                    # Add statistics text
                    textstr = f'Within 1σ: {within_1std} ({100*within_1std/len(matched_means):.1f}%)\n'
                    textstr += f'Within 2σ: {within_2std} ({100*within_2std/len(matched_means):.1f}%)'
                    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                           verticalalignment='top', 
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            
            plot_idx += 1
        
        # Hide unused subplots
        for idx in range(plot_idx, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle(f'Distribution of Matched Set Means: {prop}', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.show()


# %%
#### Side-by-Side Comparison: Mean Differences Across Gene Sets

# Create a heatmap-style comparison showing mean differences
if len(gene_set_results) > 0:
    available_props = list(gene_set_results.values())[0]['available_props']
    
    # Prepare data for heatmap
    heatmap_data = []
    for gene_set_name, result in gene_set_results.items():
        row = {'Gene Set': result['display_name']}
        for prop in available_props:
            if prop in result['real_stats'] and prop in result['matched_stats']:
                real_mean = result['real_stats'][prop]['mean']
                matched_mean = result['matched_stats'][prop]['mean']
                mean_diff_pct = ((matched_mean - real_mean) / real_mean * 100) if real_mean != 0 else np.nan
                row[prop] = mean_diff_pct
        heatmap_data.append(row)
    
    heatmap_df = pd.DataFrame(heatmap_data).set_index('Gene Set')
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    sns.heatmap(heatmap_df, annot=True, fmt='.1f', cmap='RdBu_r', center=0,
                cbar_kws={'label': 'Difference (%)'}, ax=ax, linewidths=0.5)
    ax.set_title('Mean Difference: (Matched - Real) / Real × 100%', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Property', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene Set', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    print("\nMean Difference Percentage (Matched - Real) / Real × 100%:")
    print(heatmap_df.to_string())


# %%
#### Statistical Tests: Real vs Matched Comparison

# Perform statistical tests for each gene set and property
from scipy.stats import mannwhitneyu

stat_test_results = []

for gene_set_name, result in gene_set_results.items():
    for prop in result['available_props']:
        if prop in result['real_properties'] and prop in result['matched_properties']:
            real_values = result['real_properties'][prop]
            all_matched_values = np.concatenate(result['matched_properties'][prop]) if len(result['matched_properties'][prop]) > 0 else np.array([])
            
            if len(real_values) > 0 and len(all_matched_values) > 0:
                try:
                    stat, pval = mannwhitneyu(all_matched_values, real_values, alternative='two-sided')
                    
                    real_mean = real_values.mean()
                    matched_mean = all_matched_values.mean()
                    mean_diff = matched_mean - real_mean
                    
                    stat_test_results.append({
                        'Gene Set': result['display_name'],
                        'Property': prop,
                        'Real Mean': real_mean,
                        'Matched Mean': matched_mean,
                        'Mean Difference': mean_diff,
                        'Mann-Whitney U': stat,
                        'P-value': pval,
                        'Significant (p<0.05)': 'Yes' if pval < 0.05 else 'No'
                    })
                except Exception as e:
                    print(f"Error in statistical test for {result['display_name']}, {prop}: {e}")

if len(stat_test_results) > 0:
    stat_test_df = pd.DataFrame(stat_test_results)
    print("\n" + "="*100)
    print("STATISTICAL TESTS: Mann-Whitney U Test (Real vs Matched)")
    print("="*100)
    print(stat_test_df.to_string(index=False))
    
    # Create visualization of p-values
    fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    
    # Prepare data for grouped bar plot
    gene_sets = stat_test_df['Gene Set'].unique()
    props = stat_test_df['Property'].unique()
    x = np.arange(len(gene_sets))
    width = 0.25
    
    for i, prop in enumerate(props):
        prop_data = stat_test_df[stat_test_df['Property'] == prop]
        pvals = []
        for gs in gene_sets:
            gs_data = prop_data[prop_data['Gene Set'] == gs]
            if len(gs_data) > 0:
                pvals.append(-np.log10(gs_data['P-value'].values[0]))
            else:
                pvals.append(0)
        
        ax.bar(x + i*width, pvals, width, label=prop, alpha=0.8)
    
    ax.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax.set_title('Statistical Significance: Real vs Matched Genes', fontsize=14, fontweight='bold')
    ax.set_xticks(x + width * (len(props) - 1) / 2)
    ax.set_xticklabels(gene_sets, rotation=45, ha='right')
    Pcut = 0.01
    ax.axhline(-np.log10(Pcut), color='red', linestyle='--', linewidth=2, label='p='+str(Pcut))
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    plt.tight_layout()
    plt.show()
else:
    print("No statistical test results available.")


# %%
#### Variance Tests: Real vs Matched Comparison

# Perform variance tests for each gene set and property
from scipy.stats import levene, bartlett

variance_test_results = []

for gene_set_name, result in gene_set_results.items():
    for prop in result['available_props']:
        if prop in result['real_properties'] and prop in result['matched_properties']:
            real_values = result['real_properties'][prop]
            all_matched_values = np.concatenate(result['matched_properties'][prop]) if len(result['matched_properties'][prop]) > 0 else np.array([])
            
            if len(real_values) > 0 and len(all_matched_values) > 0:
                try:
                    # Calculate variances
                    real_var = real_values.var(ddof=1)  # Sample variance
                    matched_var = all_matched_values.var(ddof=1)
                    variance_ratio = matched_var / real_var if real_var > 0 else np.nan
                    
                    # Levene's test for equal variances (more robust to non-normality)
                    # Tests null hypothesis: variances are equal
                    levene_stat, levene_pval = levene(real_values, all_matched_values)
                    
                    # Also calculate Bartlett's test (assumes normality, but good for comparison)
                    try:
                        bartlett_stat, bartlett_pval = bartlett(real_values, all_matched_values)
                    except:
                        bartlett_stat, bartlett_pval = np.nan, np.nan
                    
                    # Calculate standard deviations
                    real_std = real_values.std(ddof=1)
                    matched_std = all_matched_values.std(ddof=1)
                    std_ratio = matched_std / real_std if real_std > 0 else np.nan
                    
                    variance_test_results.append({
                        'Gene Set': result['display_name'],
                        'Property': prop,
                        'Real Variance': real_var,
                        'Matched Variance': matched_var,
                        'Variance Ratio': variance_ratio,
                        'Real Std': real_std,
                        'Matched Std': matched_std,
                        'Std Ratio': std_ratio,
                        'Levene Statistic': levene_stat,
                        'Levene P-value': levene_pval,
                        'Bartlett Statistic': bartlett_stat,
                        'Bartlett P-value': bartlett_pval,
                        'Significant (p<0.05)': 'Yes' if levene_pval < 0.05 else 'No'
                    })
                except Exception as e:
                    print(f"Error in variance test for {result['display_name']}, {prop}: {e}")

if len(variance_test_results) > 0:
    variance_test_df = pd.DataFrame(variance_test_results)
    print("\n" + "="*120)
    print("VARIANCE TESTS: Levene's Test for Equal Variances (Real vs Matched)")
    print("="*120)
    print(variance_test_df.to_string(index=False))
    
    # Create visualization of variance ratios
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
    
    # Left plot: Variance ratios
    ax1 = axes[0]
    gene_sets = variance_test_df['Gene Set'].unique()
    props = variance_test_df['Property'].unique()
    x = np.arange(len(gene_sets))
    width = 0.25
    
    for i, prop in enumerate(props):
        prop_data = variance_test_df[variance_test_df['Property'] == prop]
        ratios = []
        for gs in gene_sets:
            gs_data = prop_data[prop_data['Gene Set'] == gs]
            if len(gs_data) > 0:
                ratios.append(gs_data['Variance Ratio'].values[0])
            else:
                ratios.append(np.nan)
        
        bars = ax1.bar(x + i*width, ratios, width, label=prop, alpha=0.8)
        # Color bars: green if ratio close to 1 (good match), red if far from 1
        for bar, ratio in zip(bars, ratios):
            if not np.isnan(ratio):
                if 0.8 <= ratio <= 1.2:
                    bar.set_facecolor('green')
                elif 0.6 <= ratio <= 1.5:
                    bar.set_facecolor('lightgreen')
                elif 0.4 <= ratio <= 2.0:
                    bar.set_facecolor('yellow')
                else:
                    bar.set_facecolor('orange')
    
    ax1.axhline(1.0, color='red', linestyle='--', linewidth=2, label='Equal variance (ratio=1)')
    ax1.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Variance Ratio (Matched / Real)', fontsize=12, fontweight='bold')
    ax1.set_title('Variance Comparison: Matched vs Real Genes', fontsize=14, fontweight='bold')
    ax1.set_xticks(x + width * (len(props) - 1) / 2)
    ax1.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # Right plot: P-values from Levene's test
    ax2 = axes[1]
    for i, prop in enumerate(props):
        prop_data = variance_test_df[variance_test_df['Property'] == prop]
        pvals = []
        for gs in gene_sets:
            gs_data = prop_data[prop_data['Gene Set'] == gs]
            if len(gs_data) > 0:
                pval = gs_data['Levene P-value'].values[0]
                pvals.append(-np.log10(pval) if pval > 0 else 0)
            else:
                pvals.append(0)
        
        ax2.bar(x + i*width, pvals, width, label=prop, alpha=0.8)
    
    Pcut = 0.05
    ax2.axhline(-np.log10(Pcut), color='red', linestyle='--', linewidth=2, label=f'p={Pcut}')
    ax2.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax2.set_ylabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax2.set_title('Statistical Significance: Variance Differences (Levene\'s Test)', fontsize=14, fontweight='bold')
    ax2.set_xticks(x + width * (len(props) - 1) / 2)
    ax2.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.show()
    
    # Create heatmap of variance ratios
    heatmap_data = []
    for gene_set_name, result in gene_set_results.items():
        row = {'Gene Set': result['display_name']}
        for prop in result['available_props']:
            prop_data = variance_test_df[
                (variance_test_df['Gene Set'] == result['display_name']) & 
                (variance_test_df['Property'] == prop)
            ]
            if len(prop_data) > 0:
                row[prop] = prop_data['Variance Ratio'].values[0]
            else:
                row[prop] = np.nan
        heatmap_data.append(row)
    
    variance_heatmap_df = pd.DataFrame(heatmap_data).set_index('Gene Set')
    
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    sns.heatmap(variance_heatmap_df, annot=True, fmt='.2f', cmap='RdYlGn', center=1.0,
                cbar_kws={'label': 'Variance Ratio (Matched/Real)'}, ax=ax, linewidths=0.5,
                vmin=0.5, vmax=2.0)
    ax.set_title('Variance Ratio Heatmap: Matched / Real', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Property', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene Set', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    print("\nVariance Ratio Summary (Matched / Real):")
    print("  Ratio = 1.0: Perfect variance match")
    print("  Ratio < 1.0: Matched genes have lower variance")
    print("  Ratio > 1.0: Matched genes have higher variance")
    print("\n" + variance_heatmap_df.to_string())
else:
    print("No variance test results available.")


# %% [markdown]
# ## Testing on log10 CDS instead of CDS
#
# This section repeats the analysis using log10(CDS_length) instead of CDS_length_pct to see if p-values change.
#

# %%
# Create log10 CDS column from CDS_length
# First check if CDS_length column exists
if 'CDS_length' in match_table.columns:
    # Create log10 CDS, handling zeros and negative values
    match_table_log10 = match_table.copy()
    match_table_log10['log10_CDS_length'] = np.log10(match_table_log10['CDS_length'].clip(lower=1))
    print(f"Created log10_CDS_length column")
    print(f"CDS_length range: {match_table['CDS_length'].min():.2f} to {match_table['CDS_length'].max():.2f}")
    print(f"log10_CDS_length range: {match_table_log10['log10_CDS_length'].min():.2f} to {match_table_log10['log10_CDS_length'].max():.2f}")
    print(f"Number of valid log10_CDS_length values: {match_table_log10['log10_CDS_length'].notna().sum()}")
else:
    print("Error: CDS_length column not found in match_table")
    print(f"Available columns: {match_table.columns.tolist()}")


# %%
#### Re-analyze with log10 CDS: Real vs Matched Comparison

# Properties to compare using log10 CDS instead of CDS_length_pct
properties_to_compare_log10 = ['log10_CDS_length', 'WB_pct', 'LOEUF_pct']

# Re-analyze all gene sets with log10 CDS
print("Loading and analyzing gene sets with log10 CDS...")
gene_set_results_log10 = {}

for gene_set_name, config in gene_sets_config.items():
    print(f"\nAnalyzing {config['display_name']} with log10 CDS...")
    result = analyze_gene_set(
        gene_set_name, 
        config['file'], 
        config['display_name'],
        match_table_log10,  # Use the table with log10 CDS
        properties_to_compare_log10,
        n_matches_to_use
    )
    if result is not None:
        gene_set_results_log10[gene_set_name] = result
        file_info = f" (file: {result['file_used']})" if 'file_used' in result else ""
        print(f"  ✓ {result['n_real_genes']} real genes, {result['n_matched_sets']} matched sets{file_info}")

print(f"\n✓ Successfully analyzed {len(gene_set_results_log10)} gene sets with log10 CDS")


# %%
#### Statistical Tests: Real vs Matched Comparison (log10 CDS)

# Perform statistical tests for each gene set using log10 CDS
from scipy.stats import mannwhitneyu

stat_test_results_log10 = []

for gene_set_name, result in gene_set_results_log10.items():
    for prop in result['available_props']:
        if prop in result['real_properties'] and prop in result['matched_properties']:
            real_values = result['real_properties'][prop]
            all_matched_values = np.concatenate(result['matched_properties'][prop]) if len(result['matched_properties'][prop]) > 0 else np.array([])
            
            if len(real_values) > 0 and len(all_matched_values) > 0:
                try:
                    stat, pval = mannwhitneyu(all_matched_values, real_values, alternative='two-sided')
                    
                    real_mean = real_values.mean()
                    matched_mean = all_matched_values.mean()
                    mean_diff = matched_mean - real_mean
                    
                    stat_test_results_log10.append({
                        'Gene Set': result['display_name'],
                        'Property': prop,
                        'Real Mean': real_mean,
                        'Matched Mean': matched_mean,
                        'Mean Difference': mean_diff,
                        'Mann-Whitney U': stat,
                        'P-value': pval,
                        'Significant (p<0.05)': 'Yes' if pval < 0.05 else 'No'
                    })
                except Exception as e:
                    print(f"Error in statistical test for {result['display_name']}, {prop}: {e}")

if len(stat_test_results_log10) > 0:
    stat_test_df_log10 = pd.DataFrame(stat_test_results_log10)
    print("\n" + "="*100)
    print("STATISTICAL TESTS: Mann-Whitney U Test (Real vs Matched) - Using log10 CDS")
    print("="*100)
    print(stat_test_df_log10.to_string(index=False))
    
    # Create visualization of p-values
    fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    
    # Prepare data for grouped bar plot
    gene_sets = stat_test_df_log10['Gene Set'].unique()
    props = stat_test_df_log10['Property'].unique()
    x = np.arange(len(gene_sets))
    width = 0.25
    
    for i, prop in enumerate(props):
        prop_data = stat_test_df_log10[stat_test_df_log10['Property'] == prop]
        pvals = []
        for gs in gene_sets:
            gs_data = prop_data[prop_data['Gene Set'] == gs]
            if len(gs_data) > 0:
                pvals.append(-np.log10(gs_data['P-value'].values[0]))
            else:
                pvals.append(0)
        
        ax.bar(x + i*width, pvals, width, label=prop, alpha=0.8)
    
    ax.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax.set_title('Statistical Significance: Real vs Matched Genes (log10 CDS)', fontsize=14, fontweight='bold')
    ax.set_xticks(x + width * (len(props) - 1) / 2)
    ax.set_xticklabels(gene_sets, rotation=45, ha='right')
    Pcut = 0.01
    ax.axhline(-np.log10(Pcut), color='red', linestyle='--', linewidth=2, label='p='+str(Pcut))
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    plt.tight_layout()
    plt.show()
else:
    print("No statistical test results available.")


# %%
#### P-value Comparison: CDS_length_pct vs log10_CDS_length

# Compare p-values between the two analyses
print("\n" + "="*100)
print("P-VALUE COMPARISON: CDS_length_pct vs log10_CDS_length")
print("="*100)

# Get p-values from original analysis (CDS_length_pct)
pvals_original = {}
if len(stat_test_results) > 0:
    stat_test_df_original = pd.DataFrame(stat_test_results)
    for _, row in stat_test_df_original.iterrows():
        if row['Property'] == 'CDS_length_pct':
            key = row['Gene Set']
            pvals_original[key] = row['P-value']

# Get p-values from log10 CDS analysis
pvals_log10 = {}
if len(stat_test_results_log10) > 0:
    for _, row in stat_test_df_log10.iterrows():
        if row['Property'] == 'log10_CDS_length':
            key = row['Gene Set']
            pvals_log10[key] = row['P-value']

# Create comparison table
comparison_data = []
for gene_set in pvals_original.keys():
    if gene_set in pvals_log10:
        p_original = pvals_original[gene_set]
        p_log10 = pvals_log10[gene_set]
        p_ratio = p_log10 / p_original if p_original > 0 else np.nan
        log10_ratio = np.log10(p_log10) - np.log10(p_original) if p_original > 0 and p_log10 > 0 else np.nan
        
        comparison_data.append({
            'Gene Set': gene_set,
            'P-value (CDS_length_pct)': p_original,
            'P-value (log10_CDS_length)': p_log10,
            'P-value Ratio (log10/original)': p_ratio,
            'log10(P) Difference': log10_ratio,
            'Change': 'Increased' if p_log10 > p_original else 'Decreased' if p_log10 < p_original else 'Same'
        })

if len(comparison_data) > 0:
    comparison_df = pd.DataFrame(comparison_data)
    print("\n" + comparison_df.to_string(index=False))
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
    
    # Left plot: Side-by-side p-values
    ax1 = axes[0]
    gene_sets = comparison_df['Gene Set'].values
    x = np.arange(len(gene_sets))
    width = 0.35
    
    p_orig = comparison_df['P-value (CDS_length_pct)'].values
    p_log10 = comparison_df['P-value (log10_CDS_length)'].values
    
    bars1 = ax1.bar(x - width/2, -np.log10(p_orig), width, label='CDS_length_pct', alpha=0.8)
    bars2 = ax1.bar(x + width/2, -np.log10(p_log10), width, label='log10_CDS_length', alpha=0.8)
    
    ax1.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax1.set_ylabel('-log10(P-value)', fontsize=12, fontweight='bold')
    ax1.set_title('P-value Comparison: CDS_length_pct vs log10_CDS_length', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.3, axis='y')
    Pcut = 0.05
    ax1.axhline(-np.log10(Pcut), color='red', linestyle='--', linewidth=2, label=f'p={Pcut}')
    
    # Right plot: Ratio of p-values
    ax2 = axes[1]
    p_ratio = comparison_df['P-value Ratio (log10/original)'].values
    colors = ['green' if r < 1 else 'red' if r > 1 else 'gray' for r in p_ratio]
    bars = ax2.bar(gene_sets, p_ratio, alpha=0.8, color=colors, edgecolor='black')
    ax2.axhline(1.0, color='black', linestyle='--', linewidth=2, label='No change (ratio=1)')
    ax2.set_xlabel('Gene Set', fontsize=12, fontweight='bold')
    ax2.set_ylabel('P-value Ratio (log10 / original)', fontsize=12, fontweight='bold')
    ax2.set_title('P-value Ratio: log10 CDS / CDS_length_pct', fontsize=14, fontweight='bold')
    ax2.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, ratio in zip(bars, p_ratio):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{ratio:.2f}x',
                ha='center', va='bottom' if ratio > 1 else 'top', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY:")
    print("="*80)
    print(f"Number of gene sets analyzed: {len(comparison_data)}")
    print(f"P-values increased (less significant): {sum(comparison_df['Change'] == 'Increased')}")
    print(f"P-values decreased (more significant): {sum(comparison_df['Change'] == 'Decreased')}")
    print(f"Mean p-value ratio: {comparison_df['P-value Ratio (log10/original)'].mean():.3f}")
    print(f"Median p-value ratio: {comparison_df['P-value Ratio (log10/original)'].median():.3f}")
    print("\nInterpretation:")
    print("  Ratio < 1: log10 CDS gives more significant (smaller) p-value")
    print("  Ratio > 1: log10 CDS gives less significant (larger) p-value")
    print("  Ratio = 1: No change in p-value")
else:
    print("No comparison data available.")


# %%
#### Comparison Plots: Real vs Matched (log10 CDS)

# Create comparison plots for log10 CDS similar to the original analysis
if len(gene_set_results_log10) > 0:
    available_props = list(gene_set_results_log10.values())[0]['available_props']
    
    # Focus on log10_CDS_length for comparison
    prop = 'log10_CDS_length'
    if prop in available_props:
        n_sets = len(gene_set_results_log10)
        fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120)
        fig.patch.set_alpha(0.0)
        axes = axes.flatten()
        
        plot_idx = 0
        for gene_set_name, result in gene_set_results_log10.items():
            if plot_idx >= len(axes):
                break
                
            ax = axes[plot_idx]
            
            if prop in result['real_properties'] and prop in result['matched_properties']:
                real_values = result['real_properties'][prop]
                all_matched_values = np.concatenate(result['matched_properties'][prop]) if len(result['matched_properties'][prop]) > 0 else np.array([])
                
                if len(real_values) > 0 and len(all_matched_values) > 0:
                    # Create violin plot
                    data_to_plot = [real_values, all_matched_values]
                    labels = ['Real\nGenes', 'Matched\nGenes']
                    
                    parts = ax.violinplot(data_to_plot, positions=[1, 2], widths=0.6, showmeans=True, showmedians=True)
                    
                    # Customize violin plot
                    for pc in parts['bodies']:
                        pc.set_facecolor('lightgreen')
                        pc.set_alpha(0.7)
                    
                    # Add box plot overlay
                    bp = ax.boxplot(data_to_plot, positions=[1, 2], widths=0.3, patch_artist=False,
                                   boxprops=dict(color='black', linewidth=1.5),
                                   medianprops=dict(color='red', linewidth=2),
                                   showmeans=False)
                    
                    ax.set_xticks([1, 2])
                    ax.set_xticklabels(labels, fontsize=12)
                    ax.set_ylabel('log10(CDS_length)', fontsize=13, fontweight='bold')
                    ax.set_title(f'{result["display_name"]}\n(n={result["n_real_genes"]})', fontsize=14, fontweight='bold')
                    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
                    
                    # Calculate p-value using Mann-Whitney U test
                    try:
                        stat, pval = mannwhitneyu(all_matched_values, real_values, alternative='two-sided')
                        # Format p-value for display
                        if pval < 0.001:
                            pval_str = f'P = {pval:.2e}'.replace('e-0', ' x 10^-').replace('e-', ' x 10^-')
                        elif pval < 0.01:
                            pval_str = f'P = {pval:.2e}'.replace('e-0', ' x 10^-').replace('e-', ' x 10^-')
                        else:
                            pval_str = f'P = {pval:.3f}'
                        
                        pval_color = 'red' if pval < 0.05 else 'black'
                    except Exception as e:
                        pval_str = 'P = N/A'
                        pval_color = 'black'
                        pval = 1.0
                    
                    # Get y-axis limits to position the significance bar
                    y_min, y_max = ax.get_ylim()
                    y_range = y_max - y_min
                    bar_y = y_max - 0.05 * y_range
                    
                    # Draw horizontal bar connecting the two violins for p-value annotation
                    bar_x_start = 1
                    bar_x_end = 2
                    bar_height = 0.01 * y_range
                    
                    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    ax.plot([bar_x_start, bar_x_start], [bar_y - bar_height, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    ax.plot([bar_x_end, bar_x_end], [bar_y - bar_height, bar_y], 
                           color='black', linewidth=1.5, zorder=10)
                    
                    # Add p-value text above the bar
                    ax.text((bar_x_start + bar_x_end) / 2, bar_y + 0.02 * y_range, 
                           pval_str, ha='center', va='bottom', fontsize=11,
                           color=pval_color, fontweight='bold' if pval < 0.05 else 'normal',
                           zorder=10)
                    
                    # Add mean difference
                    real_mean = real_values.mean()
                    matched_mean = all_matched_values.mean()
                    mean_diff = matched_mean - real_mean
                    
                    textstr = f'Δ={mean_diff:.3f}'
                    ax.text(0.05, 0.05, textstr, 
                            transform=ax.transAxes, ha='left', fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
                            zorder=10)
            
            plot_idx += 1
        
        # Hide unused subplots
        for idx in range(plot_idx, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle(f'Real vs Matched: log10(CDS_length)', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.show()
    else:
        print(f"Property {prop} not found in available properties")
else:
    print("No gene sets were successfully analyzed with log10 CDS.")


# %% [markdown]
# ## Diversity Analysis

# %%
#### Gene Diversity Analysis: Fraction of Genes Appearing in 50% of Null Sets

# This analysis tests the "diversity" of matched genes
# Low diversity = many genes appear in most null sets (same genes repeatedly selected)
# High diversity = genes are more evenly distributed across null sets

def calculate_gene_diversity(gene_set_name, file_name, n_matches=1000, threshold=0.10):
    """
    Calculate what fraction of unique genes appear in at least threshold% of null sets.
    
    Args:
        gene_set_name: Name of the gene set
        file_name: File name(s) for the null weights
        n_matches: Number of matched sets to analyze
        threshold: Threshold fraction (default 0.90 = 90%)
    
    Returns:
        dict with diversity metrics
    """
    # Handle multiple possible file names
    if isinstance(file_name, list):
        file_names_to_try = file_name
    else:
        file_names_to_try = [file_name]
    
    null_genes_df = None
    used_file = None
    
    for fn in file_names_to_try:
        try:
            null_genes_df = pd.read_csv(f"{Bias_Save_Dir}/null_weights/{fn}", index_col=0)
            used_file = fn
            break
        except FileNotFoundError:
            continue
    
    if null_genes_df is None:
        return None
    
    # Get matched gene columns (skip first column which is GeneWeight)
    matched_genes_cols = null_genes_df.columns[1:].values[:n_matches]
    n_sets = len(matched_genes_cols)
    
    # Collect all unique genes that appear in matched sets
    all_matched_genes = {}
    
    for match_col in matched_genes_cols:
        matched_genes = null_genes_df[match_col].dropna().values
        for gene in matched_genes:
            if gene not in all_matched_genes:
                all_matched_genes[gene] = 0
            all_matched_genes[gene] += 1
    
    # Calculate fraction of genes appearing in at least threshold% of sets
    threshold_count = int(np.ceil(n_sets * threshold))
    genes_above_threshold = [g for g, count in all_matched_genes.items() if count >= threshold_count]
    # Also store genes with their counts for detailed analysis
    genes_above_threshold_with_counts = [(g, count) for g, count in all_matched_genes.items() if count >= threshold_count]
    # Sort by count (descending)
    genes_above_threshold_with_counts.sort(key=lambda x: x[1], reverse=True)
    total_unique_genes = len(all_matched_genes)
    
    fraction_above_threshold = len(genes_above_threshold) / total_unique_genes if total_unique_genes > 0 else 0
    
    # Calculate distribution statistics
    gene_counts = np.array(list(all_matched_genes.values()))
    
    return {
        'gene_set': gene_set_name,
        'file_used': used_file,
        'n_null_sets': n_sets,
        'threshold': threshold,
        'threshold_count': threshold_count,
        'total_unique_genes': total_unique_genes,
        'genes_above_threshold': len(genes_above_threshold),
        'fraction_above_threshold': fraction_above_threshold,
        'mean_appearances': gene_counts.mean(),
        'median_appearances': np.median(gene_counts),
        'max_appearances': gene_counts.max(),
        'min_appearances': gene_counts.min(),
        'gene_counts': gene_counts,
        'all_matched_genes': all_matched_genes,
        'genes_above_threshold_list': genes_above_threshold,
        'genes_above_threshold_with_counts': genes_above_threshold_with_counts
    }

# Calculate diversity for all gene sets
print("="*80)
print("GENE DIVERSITY ANALYSIS: Fraction of Genes Appearing in ≥50% of Null Sets")
print("="*80)
print("\nThis metric tests diversity:")
print("  - Low diversity (high fraction) = same genes repeatedly selected")
print("  - High diversity (low fraction) = genes more evenly distributed\n")

diversity_results = []
diversity_dict = {}

for gene_set_name, config in gene_sets_config.items():
    result = calculate_gene_diversity(
        gene_set_name,
        config['file'],
        n_matches=n_matches_to_use,
        threshold=0.50
    )
    
    if result is not None:
        diversity_dict[gene_set_name] = result
        diversity_results.append({
            'Gene Set': config['display_name'],
            'N Null Sets': result['n_null_sets'],
            'Total Unique Genes': result['total_unique_genes'],
            'Genes in ≥50% Sets': result['genes_above_threshold'],
            'Fraction in ≥50%': f"{result['fraction_above_threshold']*100:.2f}%",
            'Mean Appearances': f"{result['mean_appearances']:.1f}",
            'Median Appearances': f"{result['median_appearances']:.1f}",
            'Max Appearances': result['max_appearances'],
            'Min Appearances': result['min_appearances']
        })

diversity_df = pd.DataFrame(diversity_results)
print(diversity_df.to_string(index=False))
print("\n")


# %%
#### Genes Appearing in ≥50% of Null Sets

# Display the specific genes that appear in ≥50% of null sets for each disorder
print("="*80)
print("GENES APPEARING IN ≥50% OF NULL SETS")
print("="*80)
print("\n")

# Try to get gene symbols if available
try:
    HGNC, ENSID2Entrez, GeneSymbol2Entrez, Entrez2Symbol = LoadGeneINFO()
    has_gene_symbols = True
except:
    has_gene_symbols = False
    print("Note: Gene symbol mapping not available, showing Entrez IDs only\n")

for gene_set_name, result in diversity_dict.items():
    display_name = gene_sets_config[gene_set_name]['display_name']
    genes_with_counts = result['genes_above_threshold_with_counts']
    
    print(f"\n{display_name}:")
    print(f"  Total genes in ≥50% of sets: {len(genes_with_counts)}")
    print(f"  Threshold: ≥{result['threshold_count']} appearances out of {result['n_null_sets']} sets")
    
    if len(genes_with_counts) > 0:
        print(f"\n  Top genes (showing up to 50):")
        print(f"  {'Gene ID':<12} {'Symbol':<15} {'Appearances':<12} {'% of Sets':<10}")
        print(f"  {'-'*12} {'-'*15} {'-'*12} {'-'*10}")
        
        for i, (gene_id, count) in enumerate(genes_with_counts[:50]):
            pct = (count / result['n_null_sets']) * 100
            symbol = "N/A"
            if has_gene_symbols:
                try:
                    symbol = Entrez2Symbol.get(int(gene_id), "N/A")
                except:
                    symbol = "N/A"
            print(f"  {str(gene_id):<12} {symbol:<15} {count:<12} {pct:>6.1f}%")
        
        if len(genes_with_counts) > 50:
            print(f"\n  ... and {len(genes_with_counts) - 50} more genes")
    else:
        print("  No genes appear in ≥50% of null sets")
    
    print()

print("="*80)


# %%
#### Visualization: Gene Diversity Comparison (50% Threshold)

# Create visualizations comparing diversity across gene sets
if len(diversity_dict) > 0:
    # 1. Bar plot of fraction of genes in ≥50% of sets
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=120)
    
    # Left plot: Fraction in ≥50%
    ax1 = axes[0]
    gene_sets = [gene_sets_config[gs]['display_name'] for gs in diversity_dict.keys()]
    fractions = [diversity_dict[gs]['fraction_above_threshold'] * 100 for gs in diversity_dict.keys()]
    
    bars = ax1.bar(gene_sets, fractions, alpha=0.7, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Fraction of Genes in ≥50% of Sets (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Gene Diversity: Fraction Appearing in ≥50% of Null Sets', fontsize=13, fontweight='bold')
    ax1.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax1.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # Color bars based on diversity (lower = more diverse)
    for i, (bar, frac) in enumerate(zip(bars, fractions)):
        if frac < 10:
            bar.set_facecolor('green')  # High diversity
        elif frac < 25:
            bar.set_facecolor('lightgreen')  # Medium-high diversity
        elif frac < 50:
            bar.set_facecolor('yellow')  # Medium diversity
        else:
            bar.set_facecolor('orange')  # Low diversity
    
    # Add value labels on bars
    for bar, frac in zip(bars, fractions):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{frac:.1f}%',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Right plot: Distribution of gene appearances
    ax2 = axes[1]
    
    # Create box plot of gene appearance counts
    data_to_plot = [diversity_dict[gs]['gene_counts'] for gs in diversity_dict.keys()]
    labels = [gene_sets_config[gs]['display_name'] for gs in diversity_dict.keys()]
    
    bp = ax2.boxplot(data_to_plot, labels=labels, patch_artist=True,
                     boxprops=dict(linewidth=1.5),
                     medianprops=dict(color='red', linewidth=2),
                     whiskerprops=dict(linewidth=1.5),
                     capprops=dict(linewidth=1.5))
    
    # Color boxes
    colors = []
    for frac in fractions:
        if frac < 10:
            colors.append('green')
        elif frac < 25:
            colors.append('lightgreen')
        elif frac < 50:
            colors.append('yellow')
        else:
            colors.append('orange')
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax2.set_ylabel('Number of Appearances Across Null Sets', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution of Gene Appearances', fontsize=13, fontweight='bold')
    ax2.set_xticklabels(labels, rotation=45, ha='right')
    ax2.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # Add threshold line
    if len(diversity_dict) > 0:
        threshold_count = list(diversity_dict.values())[0]['threshold_count']
        ax2.axhline(threshold_count, color='red', linestyle='--', linewidth=2, 
                   label=f'50% threshold ({threshold_count} sets)')
        ax2.legend()
    
    plt.tight_layout()
    plt.show()
    
    # 2. Histogram showing distribution of appearance frequencies for each gene set
    n_sets = len(diversity_dict)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=120)
    fig.patch.set_alpha(0.0)
    axes = axes.flatten()
    
    plot_idx = 0
    for gene_set_name, result in diversity_dict.items():
        if plot_idx >= len(axes):
            break
        
        ax = axes[plot_idx]
        gene_counts = result['gene_counts']
        threshold_count = result['threshold_count']
        
        # Create histogram
        ax.hist(gene_counts, bins=30, alpha=0.7, color='lightblue', edgecolor='black')
        ax.axvline(threshold_count, color='red', linestyle='--', linewidth=2.5,
                  label=f'50% threshold ({threshold_count})')
        ax.axvline(gene_counts.mean(), color='green', linestyle='--', linewidth=2,
                  label=f'Mean ({gene_counts.mean():.1f})')
        
        ax.set_xlabel('Number of Appearances', fontsize=11, fontweight='bold')
        ax.set_ylabel('Number of Genes', fontsize=11, fontweight='bold')
        ax.set_title(f"{gene_sets_config[gene_set_name]['display_name']}\n"
                    f"Fraction in ≥50%: {result['fraction_above_threshold']*100:.1f}%",
                    fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, linestyle='--', alpha=0.3, axis='y')
        
        plot_idx += 1
    
    # Hide unused subplots
    for idx in range(plot_idx, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Distribution of Gene Appearance Frequencies Across Null Sets', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.show()
    
else:
    print("No diversity results available.")


# %% [markdown]
# ## Bias contrast

# %%
# DDD_297_Bias = pd.read_csv(Bias_Save_Dir + "DDD_297_bias_addP.csv", index_col=0)
# SCZ_200_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_200_bias_addP.csv", index_col=0)
# SCZ_500_Bias = pd.read_csv(Bias_Save_Dir + "SCZ_500_bias_addP.csv", index_col=0)
# ASD_HIQ_100_Bias = pd.read_csv(Bias_Save_Dir + "ASD_HIQ_100_bias_addP.csv", index_col=0)
# ASD_LIQ_100_Bias = pd.read_csv(Bias_Save_Dir + "ASD_LIQ_100_bias_addP.csv", index_col=0)

# %%
ASD_HIQ_DF = ASD_HIQ_100_Bias
ASD_LIQ_DF = ASD_LIQ_100_Bias
SCZ_DF = SCZ_Bias

# %%
name1="ASD"
name2="SCZ"
EffLabel = "EFFECT"
ASD_SCZ_Contrast = compare_biases(ASD_HIQ_DF, SCZ_DF, name1=name1, name2=name2, efflabel=EffLabel)
ASD_SCZ_Contrast_Neurons = ASD_SCZ_Contrast[ASD_SCZ_Contrast.index.isin(Neurons)]
plot_bias_comparison(ASD_SCZ_Contrast_Neurons, name1, name2, p_test="Mann_Whitney_FDR", legend_anchor=(0.15, 0.9))

# %%
name1="LIQ_ASD"
name2="SCZ"
EffLabel = "EFFECT"
LIQ_ASD_SCZ_Contrast = compare_biases(ASD_LIQ_DF, SCZ_DF, name1=name1, name2=name2, efflabel=EffLabel,neurons=Neurons)
LIQ_ASD_SCZ_Contrast_Neurons = LIQ_ASD_SCZ_Contrast[LIQ_ASD_SCZ_Contrast.index.isin(Neurons)]
#plot_bias_comparison(LIQ_ASD_SCZ_Contrast_Neurons, name1, name2, p_test="Bonferroni_P", legend_anchor=(0.15, 1.0))
plot_bias_comparison(LIQ_ASD_SCZ_Contrast_Neurons, name1, name2, p_test="Mann_Whitney_FDR", legend_anchor=(0.15, 1.0))   


# %%
name1="ASD"
name2="ASD with ID"
EffLabel = "EFFECT"
#HIQ_LIQ_ASD_Contrast  = compare_biases(HighIQ_ASD_Bias, LowIQ_ASD_Bias, name1="ASD", name2="ASD with ID", efflabel="EFFECT", neurons=ALL_CTs)
HIQ_LIQ_ASD_Contrast  = compare_biases(ASD_HIQ_DF, ASD_LIQ_DF, name1="ASD", name2="ASD with ID", efflabel=EffLabel, neurons=ALL_CTs)
# Set index to CT column before filtering to avoid SettingWithCopyWarning
HIQ_LIQ_ASD_Contrast_Neurons = HIQ_LIQ_ASD_Contrast[HIQ_LIQ_ASD_Contrast.index.isin(Neurons)]
plot_bias_comparison(HIQ_LIQ_ASD_Contrast_Neurons, name1, name2, p_test="Mann_Whitney_FDR", legend_anchor=(0.9, 1.0)) 
#plot_bias_comparison(HIQ_LIQ_ASD_Contrast_Neurons, name1, name2, p_test="Bonferroni_P", legend_anchor=(0.9, 1.0))
#plot_bias_comparison(HIQ_LIQ_ASD_Contrast_Neurons, name1, name2, p_test="Wilcoxon_FDR", legend_anchor=(0.9, 0))

# %%

# %%
