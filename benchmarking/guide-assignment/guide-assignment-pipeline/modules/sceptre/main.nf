// modules/sceptre/main.nf
process SCEPTRE_ASSIGN {
  tag "${dataset_id}"

  container params.sceptre_sif

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> Slurm --time

  stageInMode 'symlink' 

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignment_matrix_sceptre.rds"), emit: assignments

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignment_matrix_sceptre_${dataset_id}.rds" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # Run sceptre; --vanilla avoids reading host/user profiles or writing history.
  # pin_cores.sh: single-threaded, on exactly 1 core (see nextflow.config).
  ${projectDir}/bin/pin_cores.sh 1 \\
    Rscript --vanilla ${projectDir}/bin/run_sceptre.R ${dataset_dir} ${dataset_id}

  """
}

