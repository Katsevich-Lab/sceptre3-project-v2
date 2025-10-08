process CLEANSER_ASSIGN {
  label 'cleanser'
  tag "${dataset_id}"

  conda "${moduleDir}/environment.yml"

  cpus { resources.cpus }
  memory { resources.memory }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_cleanser.csv"), emit: assignments
  path("cleanser_${dataset_id}.time.txt"), emit: timing

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { filename ->
               filename.endsWith('.csv') ? "assignments_cleanser_${dataset_id}.csv" : null
             }

  publishDir "${outdir}/monitoring",
             mode: 'copy',
             pattern: '*.time.txt'

script:
"""
set -euo pipefail

echo "NF projectDir: ${projectDir}"
echo "dataset_dir: ${dataset_dir}"
ls -l "${dataset_dir}" || true

export PYTHONNOUSERSITE=1
[ -n "\${PYTHONPATH:-}" ] && unset PYTHONPATH
export XDG_CACHE_HOME="\$PWD/.cache";        mkdir -p "\$XDG_CACHE_HOME"

# Measure peak memory & elapsed time for cleanser
/usr/bin/time -v -o cleanser_${dataset_id}.time.txt \\
  python "${projectDir}/bin/run_cleanser.py" "${dataset_dir}/grna_matrix.mtx" "${dataset_id}"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' cleanser_${dataset_id}.time.txt
"""

}
