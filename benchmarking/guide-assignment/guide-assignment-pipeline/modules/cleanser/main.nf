process CLEANSER_ASSIGN {
  label 'cleanser'
  tag "${dataset_id}"

  conda "${moduleDir}/environment.yml"

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> SGE -l h_rt AND queue routing (>=4h -> hpc3.q)

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_cleanser.csv"), emit: assignments
  path("cleanser_${dataset_id}.time.txt"), optional: true, emit: timing

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

# Python isolation - prevent interference from user site-packages
export PYTHONNOUSERSITE=1
[ -n "\${PYTHONPATH:-}" ] && unset PYTHONPATH

# Cache directories
export XDG_CACHE_HOME="\$PWD/.cache"
export CMDSTANPY_CACHE_DIR="\$PWD/.cmdstanpy_cache"
mkdir -p "\$XDG_CACHE_HOME" "\$CMDSTANPY_CACHE_DIR"

# Measure peak memory & elapsed time for cleanser
/usr/bin/time -v -o cleanser_${dataset_id}.time.txt \\
  python "${projectDir}/bin/run_cleanser.py" "${dataset_dir}/grna_matrix.mtx" "${dataset_id}"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' cleanser_${dataset_id}.time.txt
"""

}
