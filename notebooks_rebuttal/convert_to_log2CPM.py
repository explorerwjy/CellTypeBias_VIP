#!/usr/bin/env python3
"""
Convert h5ad file from raw UMI counts to log2(CPM+1) transformation.
Processes data in chunks to avoid out-of-memory errors.

Usage:
    python convert_to_log2CPM.py [input_file] [output_file]
    
    Or with memory limit:
    ulimit -v 104857600  # 100GB in KB
    python convert_to_log2CPM.py
    
    Or using systemd-run (Linux):
    systemd-run --scope -p MemoryMax=100G python convert_to_log2CPM.py
"""

import numpy as np
import anndata as ad
import h5py
import sys
import os
import gc
from scipy import sparse

def get_memory_usage_mb():
    """Get current memory usage in MB"""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        return None

def convert_to_log2cpm_chunked(input_file, output_file, chunk_size=10000, max_memory_gb=100):
    """
    Convert AnnData file from raw counts to log2(CPM+1) in chunks.
    
    Parameters:
    -----------
    input_file : str
        Path to input h5ad file
    output_file : str
        Path to output h5ad file
    chunk_size : int
        Number of cells to process at a time (default: 10000)
    max_memory_gb : float
        Maximum memory to use in GB (default: 100)
    """
    
    print(f"Reading input file: {input_file}")
    print(f"Output file: {output_file}")
    print(f"Chunk size: {chunk_size} cells")
    print(f"Max memory: {max_memory_gb} GB")
    
    # Read in backed mode to inspect structure
    print("\nStep 1: Reading metadata...")
    adata_backed = ad.read_h5ad(input_file, backed='r')
    n_cells, n_genes = adata_backed.shape
    print(f"Dataset shape: {n_cells} cells × {n_genes} genes")
    
    # Set feature_name as index if needed
    if 'feature_name' in adata_backed.var.columns and adata_backed.var.index.name != 'feature_name':
        print("Setting feature_name as var index...")
        adata_backed.var.set_index('feature_name', inplace=True)
    
    # Calculate number of chunks needed
    n_chunks = (n_cells + chunk_size - 1) // chunk_size
    print(f"Will process in {n_chunks} chunks")
    
    # Initialize output arrays
    print("\nStep 2: Initializing output arrays...")
    # We'll write directly to output file using h5py for better memory control
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create output AnnData structure
    print("\nStep 3: Processing chunks and writing output...")
    
    # Process in chunks and accumulate results
    processed_chunks = []
    total_counts_per_cell = np.zeros(n_cells, dtype=np.float64)
    
    # First pass: calculate total counts per cell
    print("First pass: Calculating total counts per cell...")
    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, n_cells)
        chunk_indices = slice(start_idx, end_idx)
        
        # Load chunk
        chunk = adata_backed[chunk_indices].to_memory()
        
        # Get counts matrix
        if sparse.issparse(chunk.X):
            chunk_counts = chunk.X.toarray()
        else:
            chunk_counts = chunk.X
        
        # Calculate total counts per cell for this chunk
        chunk_totals = chunk_counts.sum(axis=1).astype(np.float64)
        total_counts_per_cell[start_idx:end_idx] = chunk_totals
        
        if (i + 1) % 10 == 0:
            print(f"  Processed {i + 1}/{n_chunks} chunks...")
    
    # Avoid division by zero
    total_counts_per_cell = np.maximum(total_counts_per_cell, 1.0)
    print(f"Total counts per cell range: [{total_counts_per_cell.min():.0f}, {total_counts_per_cell.max():.0f}]")
    
    # Second pass: calculate CPM and log2 transform
    print("\nSecond pass: Calculating CPM and log2(CPM+1)...")
    
    # Collect transformed chunks (using sparse matrices to save memory)
    transformed_chunks = []
    mem_usage = get_memory_usage_mb()
    if mem_usage:
        print(f"Memory usage before second pass: {mem_usage:.1f} MB")
    
    # Process chunks and transform data
    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, n_cells)
        chunk_indices = slice(start_idx, end_idx)
        
        # Load chunk
        chunk = adata_backed[chunk_indices].to_memory()
        
        # Get counts matrix
        if sparse.issparse(chunk.X):
            chunk_counts = chunk.X.toarray()
        else:
            chunk_counts = chunk.X
        
        # Calculate CPM for this chunk
        chunk_totals = total_counts_per_cell[start_idx:end_idx]
        CPM = (chunk_counts / chunk_totals[:, np.newaxis]) * 1_000_000
        
        # Apply log2(CPM + 1)
        log2cpm = np.log2(CPM + 1).astype(np.float32)
        
        # Store as sparse matrix to save memory
        transformed_chunks.append(sparse.csr_matrix(log2cpm))
        
        # Clean up
        del chunk, chunk_counts, CPM, log2cpm
        
        if (i + 1) % 10 == 0:
            mem_usage = get_memory_usage_mb()
            mem_str = f" ({mem_usage:.1f} MB)" if mem_usage else ""
            print(f"  Processed {i + 1}/{n_chunks} chunks{mem_str}...")
            # Force garbage collection periodically
            gc.collect()
    
    # Concatenate all chunks
    print("\nStep 4: Concatenating chunks...")
    output_X = sparse.vstack(transformed_chunks)
    del transformed_chunks  # Free memory
    gc.collect()
    
    # Create output AnnData object
    print("Step 5: Creating output AnnData object...")
    output_adata = ad.AnnData(
        X=output_X,
        obs=adata_backed.obs.copy(),
        var=adata_backed.var.copy()
    )
    
    # Copy additional metadata
    if hasattr(adata_backed, 'obsm'):
        output_adata.obsm = adata_backed.obsm.copy()
    if hasattr(adata_backed, 'uns'):
        output_adata.uns = adata_backed.uns.copy()
    
    print("\nStep 6: Writing output file...")
    output_adata.write(output_file)
    
    mem_usage = get_memory_usage_mb()
    mem_str = f" ({mem_usage:.1f} MB)" if mem_usage else ""
    
    print(f"\n{'='*80}")
    print(f"Conversion complete!{mem_str}")
    print(f"{'='*80}")
    print(f"Output file: {output_file}")
    print(f"Final shape: {output_adata.shape}")
    print(f"Expression range: [{output_adata.X.data.min():.2f}, {output_adata.X.data.max():.2f}]")
    print(f"{'='*80}")
    
    return output_adata


if __name__ == "__main__":
    # Default paths
    input_file = "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_CGE_interneuron.h5ad"
    output_file = "/mnt/data0/HumanBrainCellType/SuperTypeRawDat/Supercluster_CGE_interneuron_log2CPM.h5ad"
    chunk_size = 10000  # Process 10k cells at a time
    
    # Allow command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    if len(sys.argv) > 3:
        chunk_size = int(sys.argv[3])
        print(f"Using custom chunk size: {chunk_size}")
    
    print("=" * 80)
    print("Converting to log2(CPM+1)")
    print("=" * 80)
    print(f"If you encounter memory issues, reduce chunk_size (current: {chunk_size})")
    print("=" * 80)
    
    # Convert
    convert_to_log2cpm_chunked(
        input_file=input_file,
        output_file=output_file,
        chunk_size=chunk_size,
        max_memory_gb=100
    )
