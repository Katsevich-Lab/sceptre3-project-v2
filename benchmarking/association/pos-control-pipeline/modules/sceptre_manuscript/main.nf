// modules/sceptre_manuscript/main.nf
process SCEPTRE_MANUSCRIPT_POSCTRL {
  tag "${dataset_id}"

  container "${params.sceptre_manuscript_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_on_target_sceptre_manuscript.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_on_target_sceptre_manuscript_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # Number of workers for mclapply (gene/gRNA precompute + pair loop).
  export NCPUS="${task.cpus}"

  # One thread per worker: parallelism here is process-level via mclapply.
  export OMP_NUM_THREADS=1

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  Rscript --vanilla ${projectDir}/bin/run_sceptre_manuscript.R ${dataset_dir} ${dataset_id}

  """
}
