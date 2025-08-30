// Nextflow module for pertpy guide assignment (simple + strict)

process PERTPY_ASSIGN {
  label 'pertpy'
  tag "${dataset_id}"
  stageInMode 'copy'                 // copy inputs into task dir (no symlinks)

  container "${moduleDir}/pertpy.sif"

  cpus { resources.cpus }
  memory { resources.memory }

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

  # Writable caches for matplotlib/numba
  export MPLCONFIGDIR="\$PWD/.mplconfig";      mkdir -p "\$MPLCONFIGDIR"
  export NUMBA_CACHE_DIR="\$PWD/.numba_cache"; mkdir -p "\$NUMBA_CACHE_DIR"

  # Run from the dataset subdir so the script writes output there
  pushd ${dataset_dir} >/dev/null
  python ${projectDir}/bin/run_pertpy.py "\$PWD/grna_counts_pertpy.h5ad"
  # Move the output up to the task work dir with the expected name
  mv assignments_pertpy.csv ..
  popd >/dev/null
  """
}
