// modules/crispat/main.nf
// Nextflow DSL2 process: CRISPAT_ASSIGN

process CRISPAT_ASSIGN {
  label 'crispat'
  tag "${dataset_id}"

  // Use the env file that lives next to this module
  conda "${moduleDir}/environment.yml"

  // Resources set dynamically from configs/{run_id}_config.csv
  cpus { resources.cpus }
  memory { resources.memory }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_crispat.csv"), emit: assignments

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignments_crispat_${dataset_id}.csv" }

  // Use triple *single* quotes so $VARS are passed to bash untouched.
  // Inject Nextflow vars with !{...}.
  script:
  '''
  set -euo pipefail

  export PYTHONNOUSERSITE=1
  [ -n "${PYTHONPATH:-}" ] && unset PYTHONPATH
  export MPLCONFIGDIR="$PWD/.mplconfig";      mkdir -p "$MPLCONFIGDIR"
  export NUMBA_CACHE_DIR="$PWD/.numba_cache"; mkdir -p "$NUMBA_CACHE_DIR"
  export XDG_CACHE_HOME="$PWD/.cache";        mkdir -p "$XDG_CACHE_HOME"

  # Run crispat guide assignment
  python !{projectDir}/bin/run_crispat.py !{dataset_dir}/grna_counts_crispat.h5ad
  '''
}

