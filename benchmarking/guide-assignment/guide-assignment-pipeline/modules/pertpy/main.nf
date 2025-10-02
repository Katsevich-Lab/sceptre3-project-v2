// Nextflow module for pertpy guide assignment (simple + strict)

process PERTPY_ASSIGN {
  label 'pertpy'
  tag "${dataset_id}"
  stageInMode 'symlink' 

  container "${moduleDir}/pertpy.sif"

  cpus { resources.cpus }
  memory { resources.memory }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_pertpy.csv"), emit: assignments

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignments_pertpy_${dataset_id}.csv" }

  script:
  """
  set -euo pipefail

  # Run pertpy guide assignment
  python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad"
  """
}
