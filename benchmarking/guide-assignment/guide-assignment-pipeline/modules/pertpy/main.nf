// Nextflow module for pertpy guide assignment (simple + strict)

// NOTE: pertpy ignores the 'cpus' and 'memory' fields of the run config files.
// The 'gpu' block of nextflow.config is where its resources are set

process PERTPY_ASSIGN {
  label 'pertpy' // you can also add 'gpu' if you like
  tag "${dataset_id}"
  stageInMode 'symlink'
  conda "${moduleDir}/environment.yml"

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_pertpy.csv"), emit: assignments

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignments_pertpy_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  # --- GPU visibility (logs land in nf-logs/*.out) ---
  #  export LD_LIBRARY_PATH="\$CONDA_PREFIX/lib:\${LD_LIBRARY_PATH:-}"  # just unsetting this is better??
  
  unset LD_LIBRARY_PATH
  echo "LD_LIBRARY_PATH is unset"
  nvidia-smi || true



  ## confirming that GPU is present
  python - <<'PY'
  import os, jax
  print("JAX backend:", jax.default_backend()) 
  print("JAX devices:", jax.devices()) 
  PY

  # Run pertpy guide assignment
  python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad"
  """
}
