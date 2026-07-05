// Nextflow module for FR-Perturb computational benchmarking

process FRPERTURB_COMPUTATIONAL {
  tag "${dataset_id}"

  // Use conda environment from FR-Perturb repo (cached after first run)
  conda "${params.frperturb_conda_env}"

  cpus { resources.cpus }
  memory { resources.memory }

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
  # second core via BLAS. Set BEFORE python loads numpy/BLAS. Verify via trace %cpu ~100.
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1

  # Ensure writable temp directory (prevents issues on some cluster nodes)
  export TMPDIR="\$PWD/tmp"
  mkdir -p "\$TMPDIR"

  # Run FR-Perturb wrapper (outputs frperturb_results.* files to pwd)
  python ${projectDir}/bin/run_frperturb.py ${dataset_dir} ${dataset_id}

  # Verify outputs exist
  ls -lh frperturb_results*
  """
}
