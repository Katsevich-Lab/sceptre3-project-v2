// modules/mixscale/main.nf
process MIXSCALE_COMPUTATIONAL {
  tag "${dataset_id}"

  container "${params.mixscale_sif}"

  cpus   { num_cpus }
  memory { mem_str }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(mem_str), val(num_cpus)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_computational_mixscale.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_computational_mixscale_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # Fair benchmarking: disable thread-level parallelism, use process-level only
  # export OMP_NUM_THREADS=1

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # Run mixscale computational benchmarking; --vanilla avoids reading host/user profiles or writing history
  Rscript --vanilla ${projectDir}/bin/run_mixscale.R ${dataset_dir} ${dataset_id}

  """
}
