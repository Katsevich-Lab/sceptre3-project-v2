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

# Path to FR-Perturb repo (relative to bin directory)
# bin/ -> pos-control-pipeline/ -> association/ -> external/
script_dir = os.path.dirname(os.path.abspath(__file__))
frperturb_repo = os.path.join(script_dir, "../../external/FR-Perturb")
frperturb_script = os.path.join(frperturb_repo, "run_FR_Perturb.py")

# Input h5ad file (must have perturbation column in obs)
input_h5ad = os.path.join(dataset_dir, "response_matrix.h5ad")

# Run FR-Perturb
# Output will be created as frperturb_results.*, frperturb_results_*.txt
print("Running FR-Perturb...", flush=True)
subprocess.run([
    frperturb_script,
    "--input-h5ad", input_h5ad,
    "--perturbation-column-name", "perturbation",
    "--control-perturbation-name", "non-targeting",
    "--covariates", "grna_n_nonzero_subset,grna_n_umis_subset",
    "--compute-pval",
    "--fit-zero-pval",
    "--perturbation-delimiter", ":",
    # "--num-perms", 100000,
    "--out", "frperturb_results"  # Creates frperturb_results.log, etc.
], check=True)

print("FR-Perturb analysis complete!", flush=True)
print("Output files: frperturb_results.log, frperturb_results_LFCs.txt, "
      "frperturb_results_pvals.txt, frperturb_results_qvals.txt", flush=True)
