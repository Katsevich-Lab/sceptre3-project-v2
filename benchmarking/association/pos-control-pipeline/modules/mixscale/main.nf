// modules/mixscale/main.nf
process MIXSCALE_POSCTRL {
  tag "${dataset_id}"

  container "${moduleDir}/mixscale.sif"

  cpus { resources.cpus }
  memory { resources.memory }

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("association_on_target_mixscale.csv"), emit: results

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "association_on_target_mixscale_${dataset_id}.csv" }

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

  # Run mixscale positive control analysis; --vanilla avoids reading host/user profiles or writing history
  Rscript --vanilla ${projectDir}/bin/run_mixscale.R ${dataset_dir} ${dataset_id}

  """
}
