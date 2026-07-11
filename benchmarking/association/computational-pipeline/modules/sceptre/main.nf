// modules/sceptre/main.nf
process SCEPTRE_COMPUTATIONAL {
  tag "${dataset_id}"

  container "${params.sceptre_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_computational_sceptre.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_computational_sceptre_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # Fair benchmarking: pin EVERY threading layer to 1 core so no method sneaks a
  # second core via BLAS/OpenMP/etc. Set BEFORE R/Python starts and loads libraries.
  export OMP_NUM_THREADS=1
  export OMP_THREAD_LIMIT=1
  export OMP_DYNAMIC=FALSE
  export OMP_MAX_ACTIVE_LEVELS=1

  export OPENBLAS_NUM_THREADS=1
  export OPENBLAS_DEFAULT_NUM_THREADS=1
  export GOTO_NUM_THREADS=1

  export MKL_NUM_THREADS=1
  export MKL_DYNAMIC=FALSE

  export NUMEXPR_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export BLIS_NUM_THREADS=1
  export RAYON_NUM_THREADS=1
  export RCPP_PARALLEL_NUM_THREADS=1
  export R_DATATABLE_NUM_THREADS=1

  echo "Thread-control environment:"
  env | grep -E 'OMP|OPENBLAS|GOTO|MKL|NUMEXPR|VECLIB|BLIS|RAYON|RCPP|DATATABLE' | sort || true

  echo "Nextflow task cpus: ${task.cpus}"
  echo "NSLOTS: \${NSLOTS:-unset}"
  echo "Host: \$(hostname)"

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # Run sceptre computational benchmarking; --vanilla avoids reading host/user profiles or writing history
  Rscript --vanilla ${projectDir}/bin/run_sceptre.R ${dataset_dir} ${dataset_id}

  """
}
