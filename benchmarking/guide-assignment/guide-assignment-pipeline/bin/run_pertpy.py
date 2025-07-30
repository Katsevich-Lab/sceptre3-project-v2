#!/usr/bin/env python3
import sys
import pandas as pd
import anndata as ad
import pertpy as pt

input_h5ad = sys.argv[1]

# Load data and run pertpy guide assignment
adata = ad.read_h5ad(input_h5ad)
pertpy_obj = pt.pp.GuideAssignment()
pertpy_obj.assign_mixture_model(adata, assigned_guides_key="assigned_guide")

# Convert to standardized format and write output
standardized_df = pd.DataFrame({
    'cell_id': adata.obs.index,
    'grna_id': adata.obs['assigned_guide']
})

standardized_df.to_csv("assignments_pertpy.csv", index=False)