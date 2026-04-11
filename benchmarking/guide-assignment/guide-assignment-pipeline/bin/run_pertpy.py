#!/usr/bin/env python3

import jax
print("JAX devices (inside run_pertpy.py):", jax.devices(), flush=True)

import sys
import time
import pandas as pd
import anndata as ad
import pertpy as pt

input_h5ad = sys.argv[1]
dataset_id = sys.argv[2]

# Load data
adata = ad.read_h5ad(input_h5ad)

# this doesn't need to be set intelligently. We want to see every gRNA for each cell
max_assignments_per_cell = 1000  

# setting priors and MCMC params
# actually for the paper we are just using all-defaults right now
#mixture_params = {
#    'gaussian_mean_prior': (8,3),
#    # 'poisson_rate_prior': 1,  # default ok
#    'fraction_positive_expected': 0.05,
#    'num_warmup': 100,
#    'num_samples': 200
#}
#print("Mixture params:", ', '.join([f'{k}={v}' for k,v in mixture_params.items()]), flush=True)


# Run pertpy guide assignment
print("Running pertpy guide assignment...", flush=True)
# print("Note: JAX will compile on first iterations (watch for 'Compiling...' messages)", flush=True)

pertpy_obj = pt.pp.GuideAssignment()

# Time the assignment
start_time = time.time()
pertpy_obj.assign_mixture_model(
    adata,
    assigned_guides_key="assigned_guide",
    max_assignments_per_cell=max_assignments_per_cell
    # **mixture_params
)
elapsed_time = time.time() - start_time

print(f"Assignment completed in {elapsed_time/60:.2f} minutes ({elapsed_time:.1f} seconds)", flush=True)

# Convert to standardized format and write output
standardized_df = pd.DataFrame({
    'cell_id': adata.obs.index,
    'grna_id': adata.obs['assigned_guide']
})

standardized_df.to_csv("assignments_pertpy.csv", index=False)
