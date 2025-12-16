#!/usr/bin/env python3

import jax
print("JAX devices (inside run_pertpy.py):", jax.devices(), flush=True)

import sys
import pandas as pd
import anndata as ad
import pertpy as pt

input_h5ad = sys.argv[1]
dataset_id = sys.argv[2]

# Load data
adata = ad.read_h5ad(input_h5ad)

# Determine max_assignments_per_cell based on dataset_id
if 'simulated' in dataset_id:
    # For simulated data, gRNAs are iid so a cell could be assigned to all of them
    max_assignments_per_cell = adata.n_vars
elif 'replogle-rd7' in dataset_id:
    max_assignments_per_cell = 10
elif 'gasperini' in dataset_id:
    max_assignments_per_cell = 50
else:
    raise ValueError(f"Unrecognized dataset_id: {dataset_id}")

print(f"max_assignments_per_cell set to {max_assignments_per_cell} for dataset {dataset_id}", flush=True)

# Run pertpy guide assignment
pertpy_obj = pt.pp.GuideAssignment()
pertpy_obj.assign_mixture_model(adata, assigned_guides_key="assigned_guide", max_assignments_per_cell=max_assignments_per_cell)

# Convert to standardized format and write output
standardized_df = pd.DataFrame({
    'cell_id': adata.obs.index,
    'grna_id': adata.obs['assigned_guide']
})

standardized_df.to_csv("assignments_pertpy.csv", index=False)
