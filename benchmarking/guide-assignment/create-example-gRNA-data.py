#!/usr/bin/env python3
"""
Generate synthetic gRNA count data for benchmarking pipeline testing.

Creates a dataset with:
- 1000 cells and 10 gRNAs
- Exactly one gRNA per cell (100 cells per gRNA)
- gRNA-specific mean counts from N(10,1)
- Cells with gRNA: counts from Poisson(mean count)
- Cells without gRNA: counts from Poisson(0.1)
"""

import numpy as np
import pandas as pd
import anndata as adata
from pathlib import Path
import os
from scipy.sparse import csr_matrix

def generate_grna_data():
    """Generate synthetic gRNA count data according to specifications."""
    
    # Configuration parameters
    n_grnas = 3
    cells_per_grna = 100
    background_rate = 0.1
    grna_mean_count_mean = 10
    grna_mean_count_std = 1
    random_seed = 42

    n_cells = n_grnas * cells_per_grna
    
    # Set random seed for reproducibility
    np.random.seed(random_seed)
    
    # Generate gRNA names
    grna_names = [f"gRNA_{i+1}" for i in range(n_grnas)]
    
    # Generate cell barcodes
    cell_barcodes = [f"CELL_{i+1:04d}" for i in range(n_cells)]
    
    # Generate mean counts for each gRNA from N(mean, std)
    grna_mean_counts = np.random.normal(grna_mean_count_mean, grna_mean_count_std, n_grnas)
    # Ensure means are positive
    grna_mean_counts = np.maximum(grna_mean_counts, background_rate)
    
    print(f"Generated gRNA mean counts: {grna_mean_counts}")
    
    # Initialize count matrix (gRNAs x cells)
    count_matrix = np.zeros((n_grnas, n_cells))
    
    # Assign exactly one gRNA per cell (100 cells per gRNA)
    cell_assignments = np.repeat(range(n_grnas), cells_per_grna)
    
    # Generate counts for each cell
    for cell_idx in range(n_cells):
        assigned_grna = cell_assignments[cell_idx]
        
        for grna_idx in range(n_grnas):
            if grna_idx == assigned_grna:
                # Cell has this gRNA: sample from Poisson(mean count)
                count_matrix[grna_idx, cell_idx] = np.random.poisson(grna_mean_counts[grna_idx])
            else:
                # Cell doesn't have this gRNA: sample from Poisson(0.1)
                count_matrix[grna_idx, cell_idx] = np.random.poisson(background_rate)
    
    # Create AnnData object (transpose so cells are obs, gRNAs are vars)
    adata_obj = adata.AnnData(
        X=csr_matrix(count_matrix.T),
        obs=pd.DataFrame(index=cell_barcodes),
        var=pd.DataFrame(index=grna_names)
    )
    
    # Add metadata
    adata_obj.var['mean_count'] = grna_mean_counts
    adata_obj.obs['assigned_grna'] = [grna_names[cell_assignments[i]] for i in range(n_cells)]
    adata_obj.obs['batch'] = 1
    
    # Add some basic statistics
    adata_obj.obs['total_count'] = np.array(adata_obj.X.sum(axis=1)).flatten()
    adata_obj.var['total_grna_count'] = np.array(adata_obj.X.sum(axis=0)).flatten()
    
    return adata_obj

def main():
    """Generate data and save to specified location."""
    
    # Generate the data
    print("Generating synthetic gRNA count data...")
    adata_obj = generate_grna_data()
    
    # Create output directory
    output_path = Path("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/example").expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save as H5AD file
    output_file = output_path / "gRNA_counts.h5ad"
    print(f"Saving data to {output_file}")
    adata_obj.write_h5ad(output_file)
    
    # Print summary statistics
    print(f"\nDataset summary:")
    print(f"- Cells: {adata_obj.n_obs}")
    print(f"- gRNAs: {adata_obj.n_vars}")
    print(f"- Total counts: {adata_obj.X.sum():.0f}")
    print(f"- Mean counts per gRNA: {adata_obj.X.sum(axis=1).mean():.2f}")
    print(f"- Mean counts per cell: {adata_obj.X.sum(axis=0).mean():.2f}")
    
    print(f"\nFirst 5 gRNAs x first 5 cells:")
    print(adata_obj.X[:5, :5])
    
    print(f"\nData saved successfully to {output_file}")

if __name__ == "__main__":
    main()