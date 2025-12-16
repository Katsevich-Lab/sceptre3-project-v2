process CRISPAT_ASSIGN {
  label 'crispat'
  tag "${dataset_id}"

  // uses modules/crispat/environment.yml (next to this file)
  conda "${moduleDir}/environment.yml"

  cpus { resources.cpus }
  memory { resources.memory }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_crispat.csv"), emit: assignments

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { "assignments_crispat_${dataset_id}.csv" }

  // Triple *single* quotes: shell $VARS pass through unchanged.
  // Inject Nextflow vars with !{...}.
script:
"""
set -euo pipefail

echo "NF projectDir: ${projectDir}"
echo "dataset_dir: ${dataset_dir}"
ls -l "${dataset_dir}" || true

# Python isolation - prevent interference from user site-packages
export PYTHONNOUSERSITE=1
[ -n "\${PYTHONPATH:-}" ] && unset PYTHONPATH

# Cache directory for Python packages
export XDG_CACHE_HOME="\$PWD/.cache"
mkdir -p "\$XDG_CACHE_HOME"

python "${projectDir}/bin/run_crispat.py" "${dataset_dir}/grna_matrix.h5ad"
"""

}
