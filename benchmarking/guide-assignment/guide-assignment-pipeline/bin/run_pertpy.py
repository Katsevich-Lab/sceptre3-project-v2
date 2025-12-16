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

# this doesn't need to be set intelligently. We want to see every gRNA for each cell
max_assignments_per_cell = 1000  

# setting priors and MCMC params
mixture_params = {
    'gaussian_mean_prior': (8,3),
    # 'poisson_rate_prior': 1,  # default ok
    'fraction_positive_expected': 0.05,
    'num_warmup': 100,
    'num_samples': 200
}
print("Mixture params:", ', '.join([f'{k}={v}' for k,v in mixture_params.items()]), flush=True)


# Run pertpy guide assignment
pertpy_obj = pt.pp.GuideAssignment()
pertpy_obj.assign_mixture_model(
    adata,
    assigned_guides_key="assigned_guide",
    max_assignments_per_cell=max_assignments_per_cell,
    **mixture_params
)

# Convert to standardized format and write output
standardized_df = pd.DataFrame({
    'cell_id': adata.obs.index,
    'grna_id': adata.obs['assigned_guide']
})

standardized_df.to_csv("assignments_pertpy.csv", index=False)
