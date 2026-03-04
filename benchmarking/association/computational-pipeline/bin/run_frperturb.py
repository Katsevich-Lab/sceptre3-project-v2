#!/usr/bin/env python3
"""
Wrapper script for FR-Perturb computational benchmarking

Inputs (from dataset_dir):
- response_matrix.h5ad: Gene expression matrix (AnnData format) with perturbation
  assignments in obs["perturbation"] column (all log1p'd already if needed)

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

print(f"Running FR-Perturb computational benchmarking", flush=True)
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
# bin/ -> computational-pipeline/ -> association/ -> external/
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

moi_lookup = {
    "replogle": "low",
    "gasperini": "high"
}
moi = moi_lookup[dataset_name]

# P-value computation settings
use_skew_t = True  # Set to True to enable --fit-zero-pval
num_perms = 5000

if use_skew_t:
    print(f"Use skew_t: {use_skew_t}", flush=True)
else:
    print(f"Number of permutations: {num_perms}", flush=True)

# Build command list
cmd = [
    frperturb_script,
    "--input-h5ad", input_h5ad,
    "--perturbation-column-name", "perturbation",
    "--covariates", covariates,
    "--compute-pval",
    "--perturbation-delimiter", ":",
    "--out", "frperturb_results"
]

if moi == "low":
    cmd.extend(["--control-perturbation-name", "non-targeting"])

# Conditionally add skew_t flag at the end
if use_skew_t:
    cmd.append("--fit-zero-pval")
else:
    cmd.extend(["--num-perms", str(num_perms)])

print(cmd, flush=True)

# Run FR-Perturb
print("Running FR-Perturb...", flush=True)
subprocess.run(cmd, check=True)


print("FR-Perturb computational benchmarking complete!", flush=True)
print("Output files: frperturb_results.log, frperturb_results_LFCs.txt, "
      "frperturb_results_pvals.txt, frperturb_results_qvals.txt", flush=True)
