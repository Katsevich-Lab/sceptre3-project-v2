#!/usr/bin/env python3
import sys
import scipy.io
import anndata as ad
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_h5ad_to_mtx.py <input.h5ad> <output.mtx>")
        sys.exit(1)
    
    input_h5ad = sys.argv[1]
    output_mtx = sys.argv[2]
    
    # Load h5ad file
    print(f"Loading {input_h5ad}...")
    adata = ad.read_h5ad(input_h5ad)
    
    # CLEANSER expects genes x cells format, so we transpose if needed
    # Assuming input is cells x genes (standard AnnData format)
    matrix = adata.X.T  # Transpose to genes x cells
    
    # Convert to integer counts (CLEANSER expects integer count data)
    import numpy as np
    from scipy import sparse
    
    # Ensure we have integer data
    if hasattr(matrix, 'toarray'):
        # It's already sparse, just convert data to int
        matrix.data = np.round(matrix.data).astype(int)
    else:
        # Convert dense to sparse coordinate format with integer data
        matrix = np.round(matrix).astype(int)
        matrix = sparse.coo_matrix(matrix)
    
    # Write to Matrix Market coordinate format
    print(f"Writing Matrix Market format to {output_mtx}...")
    scipy.io.mmwrite(output_mtx, matrix)
    
    print(f"Successfully converted {input_h5ad} to {output_mtx}")
    print(f"Matrix shape (genes x cells): {matrix.shape}")

if __name__ == "__main__":
    main()