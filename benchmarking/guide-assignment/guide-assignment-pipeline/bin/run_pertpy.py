#!/usr/bin/env python3
import sys
import pandas as pd
import anndata as ad
import pertpy as pt
from pathlib import Path

input_h5ad = sys.argv[1]

# Load the h5ad data
print(f"Loading {input_h5ad}...")
adata = ad.read_h5ad(input_h5ad)

# Create pertpy object for guide assignment
print("Running pertpy guide assignment...")
pertpy_obj = pt.pp.GuideAssignment()

# Use mixture model assignment method (Poisson-Gaussian mixture model)
# This matches the probabilistic sophistication of crispat and CLEANSER
pertpy_obj.assign_mixture_model(adata, assigned_guides_key="assigned_guide")

# Extract assignments from the results
# pertpy adds assignments to adata.obs with the specified key
assignments = adata.obs['assigned_guide'].copy()

# Convert to standardized format (cell_id, grna_id)
standardized_df = pd.DataFrame({
    'cell_id': adata.obs.index,
    'grna_id': assignments
})

# Remove any unassigned cells (if pertpy marks them as NaN or specific value)
standardized_df = standardized_df.dropna()

print(f"Assigned {len(standardized_df)} cells to guides")

# Write standardized output
standardized_df.to_csv("assignments_pertpy.csv", index=False)