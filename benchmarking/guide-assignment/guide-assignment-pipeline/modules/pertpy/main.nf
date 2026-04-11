// Nextflow module for pertpy guide assignment (simple + strict)

// NOTE: pertpy ignores the 'cpus' field of the run config files.
// The 'memory' field is used from the config CSV (like other methods).
// The GPU queue and time are set in nextflow.config

process PERTPY_ASSIGN {
  label 'pertpy'
  tag "${dataset_id}"
  stageInMode 'symlink'
  conda "${moduleDir}/environment.yml"

  memory { resources.memory }

  // Store gpu settings in task.ext for access in nextflow.config
  ext.gpu_queue = { resources.gpu_queue }
  ext.gpu_time = { resources.gpu_time }

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

  export JAX_PLATFORMS=cpu
  export JAX_ENABLE_X64=0
  export NUMBA_CACHE_DIR=${projectDir}/.numba_cache
  export MPLCONFIGDIR=${projectDir}/.mplconfig
  export PYTHONNOUSERSITE=1

  # Run pertpy guide assignment
  python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad" ${dataset_id}
  """
}
