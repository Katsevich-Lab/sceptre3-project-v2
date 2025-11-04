// Nextflow module for pertpy guide assignment (simple + strict)

process PERTPY_ASSIGN {
  label 'pertpy','gpu' // so the process for gpu in nextflow.config also applies here
  tag "${dataset_id}"
  stageInMode 'symlink' 

  // container "${moduleDir}/pertpy.sif"
  conda "${moduleDir}/environment.yml"

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

# --- GPU visibility (logs land in nf-logs/*.out) ---
 nvidia-smi || true
  python - <<'PY'
import os, jax
print("JAX version:", jax.__version__)
print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))
print("JAX devices:", jax.devices())
PY


  # Run pertpy guide assignment
  python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad"
  """
}
