// Nextflow module for FR-Perturb computational benchmarking

process FRPERTURB_COMPUTATIONAL {
  tag "${dataset_id}"

  // Use conda environment from FR-Perturb repo (cached after first run)
  conda "${params.frperturb_conda_env}"

  cpus { resources.cpus }
  memory { resources.memory }
  time { resources.time }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  // Collect all 4 output files created by FR-Perturb
  tuple val(dataset_id), val(method), path("frperturb_results*"), emit: results

  // Publish to frperturb subdirectory, preserving filenames but adding dataset suffix
  publishDir "${outdir}/frperturb",
             mode: 'copy',
             saveAs: { filename ->
               // frperturb_results.log -> frperturb_results_replogle-rd7_neg_control_100genes.log
               filename.replace("frperturb_results", "frperturb_results_${dataset_id}")
             }

  script:
  """
  set -euo pipefail

  echo "Dataset: ${dataset_id}"
  echo "Dataset directory: ${dataset_dir}"
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

  # Ensure writable temp directory (prevents issues on some cluster nodes)
  export TMPDIR="\$PWD/tmp"
  mkdir -p "\$TMPDIR"

  # Run FR-Perturb wrapper (outputs frperturb_results.* files to pwd)
  python ${projectDir}/bin/run_frperturb.py ${dataset_dir} ${dataset_id}

  # Verify outputs exist
  ls -lh frperturb_results*
  """
}
