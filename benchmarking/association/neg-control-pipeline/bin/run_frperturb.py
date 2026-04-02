#!/usr/bin/env python3
"""
Wrapper script for FR-Perturb negative control analysis

Key differences from positive control:
- No --control-perturbation-name specified (compare each pseudo-target against average)
- Uses _full_log1p covariates (already log1p-transformed in h5ad)

Inputs (from dataset_dir):
- response_matrix.h5ad: Gene expression matrix (AnnData format) with perturbation
  assignments in obs["perturbation"] column (pseudo-target names like "non-targeting1")

Outputs (created by FR-Perturb in current directory):
- frperturb_results.log
- frperturb_results_LFCs.txt
- frperturb_results_pvals.txt
- frperturb_results_qvals.txt
"""

import sys
import os
import subprocess

dataset_dir = sys.argv[1]
dataset_id = sys.argv[2]

print(f"Running FR-Perturb negative control analysis", flush=True)
print(f"Dataset directory: {dataset_dir}", flush=True)
print(f"Dataset ID: {dataset_id}", flush=True)

# Determine which dataset this is
DATASET_NAMES = ["gasperini", "replogle"]
dataset_name = [name for name in DATASET_NAMES if name in dataset_id.lower()]
if len(dataset_name) != 1:
    raise ValueError(f"Could not determine dataset from dataset_id: {dataset_id}")
dataset_name = dataset_name[0]
print(f"Detected dataset: {dataset_name}", flush=True)

# Path to FR-Perturb repo (relative to bin directory)
# bin/ -> neg-control-pipeline/ -> association/ -> external/
script_dir = os.path.dirname(os.path.abspath(__file__))
frperturb_repo = os.path.join(script_dir, "../../external/FR-Perturb")
frperturb_script = os.path.join(frperturb_repo, "run_FR_Perturb.py")

# Input h5ad file (must have perturbation column in obs)
input_h5ad = os.path.join(dataset_dir, "response_matrix.h5ad")

# Set covariates based on dataset
# Note: Using _full_log1p versions (already log1p-transformed in h5ad)
covariates_lookup = {
    "replogle": "response_n_nonzero_full_log1p,response_n_umis_full_log1p,grna_n_nonzero_full_log1p,grna_n_umis_full_log1p",
    "gasperini": "response_n_nonzero_full_log1p,response_n_umis_full_log1p,grna_n_nonzero_full_log1p,grna_n_umis_full_log1p,prep_batch"
}
covariates = covariates_lookup[dataset_name]
print(f"Covariates: {covariates}", flush=True)

# P-value computation settings
use_skew_t = True  # Set to True to enable --fit-zero-pval
num_perms = 1000

if use_skew_t:
    print(f"Use skew_t: {use_skew_t}", flush=True)
else:
    print(f"Number of permutations: {num_perms}", flush=True)

delimiter_lookup = {
    "replogle": None,   # low MOI: single perturbation per cell, no delimiter needed
    "gasperini": "@"    # high MOI: multiple targets concatenated with "@"
}
delimiter = delimiter_lookup[dataset_name]

# Build command list
# NOTE: No --control-perturbation-name specified - each pseudo-target compared against average
cmd = [
    frperturb_script,
    "--input-h5ad", input_h5ad,
    "--perturbation-column-name", "perturbation",
    "--covariates", covariates,
    "--compute-pval",
    "--out", "frperturb_results"
]

if delimiter is not None:
    cmd.extend(["--perturbation-delimiter", delimiter])

# Conditionally add skew_t flag at the end
if use_skew_t:
    cmd.append("--fit-zero-pval")
else:
    cmd.extend(["--num-perms", str(num_perms)])

print(cmd, flush=True)

# Run FR-Perturb
print("Running FR-Perturb...", flush=True)
subprocess.run(cmd, check=True)


print("FR-Perturb negative control analysis complete!", flush=True)
print("Output files: frperturb_results.log, frperturb_results_LFCs.txt, "
      "frperturb_results_pvals.txt, frperturb_results_qvals.txt", flush=True)
