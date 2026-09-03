// modules/sceptre_v030/main.nf
process SCEPTRE_V030_POSCTRL {
  tag "${dataset_id}"

  container "${params.sceptre_v030_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_on_target_sceptre_v030.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_on_target_sceptre_v030_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # v0.3.0 has no internal parallelism; pin threading so it cannot grab extra cores.
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  Rscript --vanilla ${projectDir}/bin/run_sceptre_v030.R ${dataset_dir} ${dataset_id}

  """
}
