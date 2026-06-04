// modules/script/main.nf
//
// Runs a "script_*" guide-assignment variant by dispatching to bin/run_script.R,
// which sources bin/script/{variant}.R based on the method name.
// Reuses the sceptre.sif container (R + Matrix + sceptre's transitive deps).

process SCRIPT_ASSIGN {
  tag "${method}/${dataset_id}"

  container "${projectDir}/modules/sceptre/sceptre.sif"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignment_matrix_script.rds"), emit: assignments
  path "assignment_extras_script.rds", optional: true, emit: extras

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { fname ->
               fname.replaceAll(/_script\.rds$/, "_${method}_${dataset_id}.rds")
             }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  echo "method:      ${method}"
  echo "cpus:        ${task.cpus}"
  ls -l "${dataset_dir}" || true

  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  Rscript --vanilla ${projectDir}/bin/run_script.R \\
      ${dataset_dir} ${dataset_id} ${method} ${task.cpus}
  """
}
