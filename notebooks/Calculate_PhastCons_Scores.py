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
# # Calculate phastCons 100-way Scores for Protein-Coding Genes
#
# This notebook calculates phastCons conservation scores for all protein-coding genes in GENCODE v19 (hg19).
#
# **Input files:**
# - GENCODE v19 GTF: `/home/jw3514/Work/Resources/gencode.v19.annotation.gtf.gz`
# - phastCons 100-way BigWig: `/home/jw3514/Work/Resources/hg19.100way.phastCons.bw`
#
# **Strategy:**
# 1. Parse GTF to extract all CDS regions for protein-coding genes
# 2. Merge CDS coordinates (union across all transcripts per gene)
# 3. Calculate mean, median, and max phastCons scores across merged CDS
# 4. Map Ensembl gene IDs to Entrez IDs and gene symbols
# 5. Export results
#
# **Output:**
# - CSV with columns: EntrezID, GeneSymbol, EnsemblID, mean_phastCons, median_phastCons, max_phastCons, n_CDS_bases

# %%
import gzip
import pandas as pd
import numpy as np
import pyBigWig
from collections import defaultdict
from tqdm.auto import tqdm
import mygene
import warnings
warnings.filterwarnings('ignore')

print(f"pyBigWig version: {pyBigWig.__version__}")
print(f"pandas version: {pd.__version__}")


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

# %% [markdown]
# ## Step 1: Parse GENCODE v19 GTF and Extract CDS Regions

# %%
# File paths
GTF_FILE = "/home/jw3514/Work/Resources/gencode.v19.annotation.gtf.gz"
PHASTCONS_FILE = "/home/jw3514/Work/Resources/hg19.100way.phastCons.bw"
OUTPUT_FILE = str(PROJ_DIR / "results/phastCons_scores_all_genes.csv") + "/"

print(f"GTF file: {GTF_FILE}")
print(f"phastCons file: {PHASTCONS_FILE}")
print(f"Output file: {OUTPUT_FILE}")


# %%
def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into a dictionary."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            key, value = item.split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

def extract_cds_regions(gtf_file):
    """
    Extract all CDS regions for protein-coding genes from GENCODE GTF.
    Returns a dictionary: {gene_id: {'symbol': str, 'chr': str, 'cds': [(start, end), ...]}}
    """
    gene_data = defaultdict(lambda: {'symbol': None, 'chr': None, 'cds': []})
    
    print("Parsing GTF file...")
    with gzip.open(gtf_file, 'rt') as f:
        for line in tqdm(f, desc="Reading GTF"):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            # Parse attributes
            attrs = parse_gtf_attributes(attributes)
            
            # Only process protein-coding genes
            if attrs.get('gene_type') != 'protein_coding':
                continue
            
            gene_id = attrs.get('gene_id', '').split('.')[0]  # Remove version
            
            # Extract CDS regions
            if feature == 'CDS':
                gene_data[gene_id]['chr'] = chrom
                gene_data[gene_id]['symbol'] = attrs.get('gene_name', gene_id)
                gene_data[gene_id]['cds'].append((int(start), int(end)))
    
    print(f"\nFound {len(gene_data)} protein-coding genes with CDS annotations")
    return dict(gene_data)

# Extract CDS regions
gene_cds_data = extract_cds_regions(GTF_FILE)

# %%
# Check the first few genes
print("Sample of extracted genes:")
for i, (gene_id, data) in enumerate(list(gene_cds_data.items())[:3]):
    print(f"\n{gene_id} ({data['symbol']}):")
    print(f"  Chromosome: {data['chr']}")
    print(f"  Number of CDS regions: {len(data['cds'])}")
    print(f"  First 3 CDS: {data['cds'][:3]}")


# %% [markdown]
# ## Step 2: Merge CDS Regions (Union Across Transcripts)

# %%
def merge_intervals(intervals):
    """
    Merge overlapping intervals.
    Input: list of (start, end) tuples
    Output: list of merged (start, end) tuples
    """
    if not intervals:
        return []
    
    # Sort by start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:  # Overlapping
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    
    return merged

# Merge CDS regions for each gene
print("Merging CDS regions...")
for gene_id in tqdm(gene_cds_data.keys(), desc="Merging CDS"):
    gene_cds_data[gene_id]['merged_cds'] = merge_intervals(gene_cds_data[gene_id]['cds'])

print("\nMerging complete!")
print("\nExample - before and after merging:")
example_gene = list(gene_cds_data.keys())[0]
print(f"Gene: {example_gene} ({gene_cds_data[example_gene]['symbol']})")
print(f"  Original CDS count: {len(gene_cds_data[example_gene]['cds'])}")
print(f"  Merged CDS count: {len(gene_cds_data[example_gene]['merged_cds'])}")

# %% [markdown]
# ## Step 3: Calculate phastCons Scores
#
# For each gene's merged CDS regions, we'll extract phastCons scores and calculate:
# - **Mean**: Average conservation across all CDS bases
# - **Median**: Median conservation score
# - **Max**: Maximum conservation score
# - **Total CDS length**: Number of bases in merged CDS

# %%
# Open phastCons BigWig file
print(f"Opening phastCons file: {PHASTCONS_FILE}")
bw = pyBigWig.open(PHASTCONS_FILE)
print(f"Chromosomes in BigWig: {list(bw.chroms().keys())[:10]}...")  # Show first 10


# %%
def calculate_phastcons_for_gene(gene_id, gene_info, bw_file):
    """
    Calculate phastCons scores for a gene's merged CDS regions.
    Returns: dict with mean, median, max phastCons and total CDS length
    """
    chrom = gene_info['chr']
    merged_cds = gene_info['merged_cds']
    
    if not merged_cds:
        return {'mean': np.nan, 'median': np.nan, 'max': np.nan, 'n_bases': 0}
    
    # Collect all phastCons values across merged CDS
    all_scores = []
    total_bases = 0
    
    for start, end in merged_cds:
        try:
            # pyBigWig uses 0-based coordinates, GTF is 1-based
            # Convert GTF 1-based to 0-based
            scores = bw_file.values(chrom, start - 1, end)
            
            # Filter out None values (regions without data)
            valid_scores = [s for s in scores if s is not None]
            all_scores.extend(valid_scores)
            total_bases += (end - start + 1)
        except:
            # Handle chromosomes not in BigWig or out-of-bounds regions
            continue
    
    if not all_scores:
        return {'mean': np.nan, 'median': np.nan, 'max': np.nan, 'n_bases': total_bases}
    
    return {
        'mean': np.mean(all_scores),
        'median': np.median(all_scores),
        'max': np.max(all_scores),
        'n_bases': total_bases
    }

# Calculate phastCons for all genes
print("Calculating phastCons scores for all genes...")
phastcons_results = {}

for gene_id in tqdm(gene_cds_data.keys(), desc="Computing phastCons"):
    phastcons_results[gene_id] = calculate_phastcons_for_gene(gene_id, gene_cds_data[gene_id], bw)

# Close BigWig file
bw.close()
print("\nphastCons calculation complete!")

# %%
# Check results
print("Sample phastCons results:")
for i, (gene_id, scores) in enumerate(list(phastcons_results.items())[:5]):
    symbol = gene_cds_data[gene_id]['symbol']
    print(f"{gene_id} ({symbol}):")
    print(f"  Mean: {scores['mean']:.4f}, Median: {scores['median']:.4f}, Max: {scores['max']:.4f}")
    print(f"  CDS bases: {scores['n_bases']}")

# Summary statistics
n_genes_with_scores = sum(1 for s in phastcons_results.values() if not np.isnan(s['mean']))
print(f"\n{n_genes_with_scores} / {len(phastcons_results)} genes have phastCons scores")

# %% [markdown]
# ## Step 4: Map Ensembl Gene IDs to Entrez IDs and Gene Symbols
#
# We'll use the `mygene.info` API to map Ensembl gene IDs to Entrez IDs.

# %%
# Initialize mygene client
mg = mygene.MyGeneInfo()

# Get all Ensembl gene IDs
ensembl_ids = list(gene_cds_data.keys())
print(f"Mapping {len(ensembl_ids)} Ensembl IDs to Entrez IDs...")

# Query in batches (mygene has a limit of 1000 per query)
batch_size = 1000
entrez_mapping = {}

for i in tqdm(range(0, len(ensembl_ids), batch_size), desc="Querying mygene"):
    batch = ensembl_ids[i:i+batch_size]
    
    # Query mygene
    results = mg.querymany(
        batch, 
        scopes='ensembl.gene',
        fields='entrezgene,symbol',
        species='human',
        returnall=True
    )
    
    # Parse results
    for result in results['out']:
        ensembl_id = result['query']
        entrez_id = result.get('entrezgene', None)
        symbol = result.get('symbol', None)
        
        entrez_mapping[ensembl_id] = {
            'entrez_id': entrez_id,
            'symbol': symbol if symbol else gene_cds_data[ensembl_id]['symbol']
        }

# Fill in missing mappings
for ensembl_id in ensembl_ids:
    if ensembl_id not in entrez_mapping:
        entrez_mapping[ensembl_id] = {
            'entrez_id': None,
            'symbol': gene_cds_data[ensembl_id]['symbol']
        }

n_with_entrez = sum(1 for v in entrez_mapping.values() if v['entrez_id'] is not None)
print(f"\nMapped {n_with_entrez} / {len(ensembl_ids)} genes to Entrez IDs")

# %%
# Check some mappings
print("Sample ID mappings:")
for i, (ensembl_id, mapping) in enumerate(list(entrez_mapping.items())[:10]):
    print(f"{ensembl_id} -> Entrez: {mapping['entrez_id']}, Symbol: {mapping['symbol']}")

# %% [markdown]
# ## Step 5: Combine Results and Export

# %%
# Build final dataframe
print("Building final dataframe...")
results_list = []

for ensembl_id in gene_cds_data.keys():
    phast_scores = phastcons_results[ensembl_id]
    mapping = entrez_mapping[ensembl_id]
    
    results_list.append({
        'EnsemblID': ensembl_id,
        'EntrezID': mapping['entrez_id'],
        'GeneSymbol': mapping['symbol'],
        'Chromosome': gene_cds_data[ensembl_id]['chr'],
        'mean_phastCons': phast_scores['mean'],
        'median_phastCons': phast_scores['median'],
        'max_phastCons': phast_scores['max'],
        'n_CDS_bases': phast_scores['n_bases']
    })

df_results = pd.DataFrame(results_list)

# Sort by Entrez ID
df_results = df_results.sort_values('EntrezID')

print(f"\nFinal dataframe shape: {df_results.shape}")
print(f"Columns: {list(df_results.columns)}")
df_results.head(10)

# %%
# Summary statistics
print("Summary Statistics:")
print("=" * 60)
print(f"Total genes: {len(df_results)}")
print(f"Genes with Entrez ID: {df_results['EntrezID'].notna().sum()}")
print(f"Genes with phastCons scores: {df_results['mean_phastCons'].notna().sum()}")
print(f"\nphastCons score distribution:")
print(df_results[['mean_phastCons', 'median_phastCons', 'max_phastCons']].describe())

# %%
# Export to CSV
print(f"Exporting results to: {OUTPUT_FILE}")
df_results.to_csv(OUTPUT_FILE, index=False)
print("Export complete!")

# Also create a version with only genes that have Entrez IDs
OUTPUT_FILE_ENTREZ_ONLY = OUTPUT_FILE.replace('.csv', '_EntrezOnly.csv')
df_entrez = df_results[df_results['EntrezID'].notna()].copy()
df_entrez['EntrezID'] = df_entrez['EntrezID'].astype(int)
df_entrez.to_csv(OUTPUT_FILE_ENTREZ_ONLY, index=False)
print(f"\nAlso exported Entrez-only version ({len(df_entrez)} genes): {OUTPUT_FILE_ENTREZ_ONLY}")

# %% [markdown]
# ## Validation and Quality Control

# %%
# Plot distribution of phastCons scores
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Mean phastCons
axes[0].hist(df_results['mean_phastCons'].dropna(), bins=50, edgecolor='black')
axes[0].set_xlabel('Mean phastCons')
axes[0].set_ylabel('Number of genes')
axes[0].set_title('Distribution of Mean phastCons Scores')

# Median phastCons
axes[1].hist(df_results['median_phastCons'].dropna(), bins=50, edgecolor='black', color='orange')
axes[1].set_xlabel('Median phastCons')
axes[1].set_ylabel('Number of genes')
axes[1].set_title('Distribution of Median phastCons Scores')

# Max phastCons
axes[2].hist(df_results['max_phastCons'].dropna(), bins=50, edgecolor='black', color='green')
axes[2].set_xlabel('Max phastCons')
axes[2].set_ylabel('Number of genes')
axes[2].set_title('Distribution of Max phastCons Scores')

plt.tight_layout()
plt.show()

# %%
# Check some well-known genes
test_genes = ['BRCA1', 'TP53', 'MECP2', 'SCN2A', 'CACNA1C']
print("\nphastCons scores for well-known genes:")
print("=" * 80)
for gene in test_genes:
    gene_data = df_results[df_results['GeneSymbol'] == gene]
    if not gene_data.empty:
        row = gene_data.iloc[0]
        print(f"{gene}:")
        print(f"  Entrez ID: {row['EntrezID']}")
        print(f"  Mean: {row['mean_phastCons']:.4f}, Median: {row['median_phastCons']:.4f}, Max: {row['max_phastCons']:.4f}")
        print(f"  CDS bases: {row['n_CDS_bases']}")
    else:
        print(f"{gene}: Not found")

# %%
# Top 20 most conserved genes by mean phastCons
print("\nTop 20 most conserved genes (by mean phastCons):")
print("=" * 80)
top_conserved = df_results.nlargest(20, 'mean_phastCons')[['GeneSymbol', 'EntrezID', 'mean_phastCons', 'median_phastCons', 'max_phastCons', 'n_CDS_bases']]
print(top_conserved.to_string(index=False))

# %% [markdown]
# ## Summary
#
# **Outputs:**
# 1. `phastCons_scores_all_genes.csv` - All protein-coding genes
# 2. `phastCons_scores_all_genes_EntrezOnly.csv` - Only genes with Entrez IDs
#
# **Columns:**
# - `EnsemblID`: Ensembl gene ID
# - `EntrezID`: Entrez gene ID (if available)
# - `GeneSymbol`: Gene symbol
# - `Chromosome`: Chromosome location
# - `mean_phastCons`: Mean phastCons score across all CDS bases
# - `median_phastCons`: Median phastCons score
# - `max_phastCons`: Maximum phastCons score
# - `n_CDS_bases`: Total number of CDS bases (merged)
