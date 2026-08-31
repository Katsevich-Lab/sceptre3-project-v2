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
  time   { resources.time }   // -> SGE -l h_rt AND queue routing (>=4h -> hpc3.q)

  // Store gpu settings in task.ext for access in nextflow.config
  ext.gpu_queue = { resources.gpu_queue }
  ext.gpu_time = { resources.gpu_time }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_pertpy.csv"), emit: assignments
  path("pertpy_${dataset_id}.time.txt"), optional: true, emit: timing

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignments_pertpy_${dataset_id}.csv" }

  publishDir "${outdir}/monitoring",
             mode: 'copy',
             pattern: '*.time.txt'

  script:
  """
  set -euo pipefail

  export JAX_PLATFORMS=cpu
  export JAX_ENABLE_X64=0
  export NUMBA_CACHE_DIR=${projectDir}/.numba_cache
  export MPLCONFIGDIR=${projectDir}/.mplconfig
  export PYTHONNOUSERSITE=1

  # Run pertpy guide assignment, measuring peak memory & elapsed time
  # (parity with the cleanser module)
  /usr/bin/time -v -o pertpy_${dataset_id}.time.txt \\
    python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad" ${dataset_id}

  # Print a one-line summary into .command.out for convenience
  awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
       /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' pertpy_${dataset_id}.time.txt

  """
}
