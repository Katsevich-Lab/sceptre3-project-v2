// modules/sceptre_v030/main.nf
process SCEPTRE_V030_NEGCTRL {
  tag "${dataset_id}"

  container "${params.sceptre_v030_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_neg_control_sceptre_v030.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_neg_control_sceptre_v030_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # v0.3.0 has no internal parallelism, so this runs on one core regardless of the
  # cpus allocation; pin the threading layers so it cannot grab extra cores via BLAS.
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # The image installs sceptre pinned to tag v0.3.0 as the sole `sceptre`.
  Rscript --vanilla ${projectDir}/bin/run_sceptre_v030.R ${dataset_dir} ${dataset_id}

  """
}
