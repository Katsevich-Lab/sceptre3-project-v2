#!/usr/bin/env python3

# <--- checking that everything seems right for the GPU ---> #
import os
# Friendlier VRAM behavior + compile logging (also set in nextflow.config)
# os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
# os.environ.setdefault("JAX_LOG_COMPILES", "1")

# Initialize JAX persistent compilation cache
try:
    from jax.experimental.compilation_cache import initialize_cache
    cache_dir = os.environ.get("JAX_COMPILATION_CACHE_DIR", ".jax_cache")
    os.makedirs(cache_dir, exist_ok=True)
    initialize_cache(cache_dir)
    print("JAX persistent cache:", cache_dir, flush=True)
except Exception as e:
    print("JAX cache init skipped:", e, flush=True)

import jax
print("JAX devices (inside run_pertpy.py):", jax.devices(), flush=True)

# <--- end GPU checks ---> #


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
