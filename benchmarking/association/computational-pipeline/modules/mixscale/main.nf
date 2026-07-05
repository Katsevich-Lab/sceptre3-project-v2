// modules/mixscale/main.nf
process MIXSCALE_COMPUTATIONAL {
  tag "${dataset_id}"

  container "${params.mixscale_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
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

  # Fair benchmarking: pin EVERY threading layer to 1 core so no method sneaks a
  # second core via BLAS. Set BEFORE Rscript loads BLAS. Verify via trace %cpu ~100.
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # Run mixscale computational benchmarking; --vanilla avoids reading host/user profiles or writing history
  Rscript --vanilla ${projectDir}/bin/run_mixscale.R ${dataset_dir} ${dataset_id}

  """
}
