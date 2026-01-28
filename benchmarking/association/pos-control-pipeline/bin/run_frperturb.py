#!/usr/bin/env python3
"""
Wrapper script for FR-Perturb positive control analysis

Inputs (from dataset_dir):
- response_matrix.h5ad: Gene expression matrix (AnnData format) with perturbation
  assignments in obs["perturbation"] column

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

print(f"Running FR-Perturb positive control analysis", flush=True)
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
# bin/ -> pos-control-pipeline/ -> association/ -> external/
script_dir = os.path.dirname(os.path.abspath(__file__))
frperturb_repo = os.path.join(script_dir, "../../external/FR-Perturb")
frperturb_script = os.path.join(frperturb_repo, "run_FR_Perturb.py")

# Input h5ad file (must have perturbation column in obs)
input_h5ad = os.path.join(dataset_dir, "response_matrix.h5ad")

# Set covariates based on dataset
covariates_lookup = {
    "replogle": "grna_n_nonzero_subset,grna_n_umis_subset",
    "gasperini": "grna_n_nonzero_subset,grna_n_umis_subset,prep_batch"
}
covariates = covariates_lookup[dataset_name]
print(f"Covariates: {covariates}", flush=True)

# Run FR-Perturb
# Output will be created as frperturb_results.*, frperturb_results_*.txt
print("Running FR-Perturb...", flush=True)
subprocess.run([
    frperturb_script,
    "--input-h5ad", input_h5ad,
    "--perturbation-column-name", "perturbation",
    "--control-perturbation-name", "non-targeting",
    "--covariates", covariates,
    "--compute-pval",
    "--fit-zero-pval",
    "--perturbation-delimiter", ":",
    # "--num-perms", 100000,
    "--out", "frperturb_results"  # Creates frperturb_results.log, etc.
], check=True)

print("FR-Perturb analysis complete!", flush=True)
print("Output files: frperturb_results.log, frperturb_results_LFCs.txt, "
      "frperturb_results_pvals.txt, frperturb_results_qvals.txt", flush=True)
