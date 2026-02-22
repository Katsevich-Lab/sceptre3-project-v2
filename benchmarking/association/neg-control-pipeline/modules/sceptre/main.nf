// modules/sceptre/main.nf
process SCEPTRE_NEGCTRL {
  tag "${dataset_id}"

  container "${params.sceptre_sif}"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_neg_control_sceptre.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_neg_control_sceptre_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  echo "dataset_dir: ${dataset_dir}"
  ls -l "${dataset_dir}" || true

  # Set number of processors for R (from Nextflow task allocation)
  export NCPUS="${task.cpus}"

  # prevents too many resources being used
  # NOTE: this is only for neg-control since i will always be 
  # running sceptre in parallel here
  export OMP_NUM_THREADS=1


  # R needs writable temp and (if any package tries) a user lib dir
  export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
  export R_TMPDIR="\$TMPDIR"
  export R_USER="\$PWD"
  export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

  # Run sceptre negative control analysis; --vanilla avoids reading host/user profiles or writing history
  Rscript --vanilla ${projectDir}/bin/run_sceptre.R ${dataset_dir} ${dataset_id}

  """
}
