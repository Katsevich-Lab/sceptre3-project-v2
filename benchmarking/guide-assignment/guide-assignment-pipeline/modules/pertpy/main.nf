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

  # --- GPU visibility ---
  unset LD_LIBRARY_PATH
  echo "LD_LIBRARY_PATH is unset"
  nvidia-smi || true

  ## confirming that GPU is present
  python - <<'PY'
  import os, jax
  print("JAX backend:", jax.default_backend())
  print("JAX devices:", jax.devices())
  PY

  # Export JAX/GPU settings for conda environment
  export JAX_LOG_COMPILES=1
  export JAX_PLATFORMS=cuda
  export XLA_PYTHON_CLIENT_PREALLOCATE=false
  export XLA_PYTHON_CLIENT_MEM_FRACTION=0.85

  # JAX persistent compilation cache
  export JAX_COMPILATION_CACHE_DIR=${projectDir}/.jax_cache
  export JAX_COMPILATION_CACHE_MAX_SIZE=${25L * 1024 * 1024 * 1024}
  export JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS=0
  export JAX_PERSISTENT_CACHE_MIN_ENTRY_SIZE_BYTES=-1
  export JAX_PERSISTENT_CACHE_ENABLE_XLA_CACHES=all

  # Other caches and settings
  export NUMBA_CACHE_DIR=${projectDir}/.numba_cache
  export MPLCONFIGDIR=${projectDir}/.mplconfig
  export PYTHONNOUSERSITE=1

  # Thread limits (not critical for GPU-bound pertpy, uncomment if needed)
  # export OMP_NUM_THREADS=1
  # export OPENBLAS_NUM_THREADS=1
  # export MKL_NUM_THREADS=1

  # Run pertpy guide assignment
  python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad" ${dataset_id}
  """
}
