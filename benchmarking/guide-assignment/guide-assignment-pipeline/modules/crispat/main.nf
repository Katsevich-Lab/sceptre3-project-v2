process CRISPAT_ASSIGN {
  tag "${dataset_id}"

  // Exact lock of the env verified on Betty; environment.yml (next to this file)
  // holds the requirements and why. The only place this env is chosen.
  conda "${moduleDir}/environment.lock.yml"

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> Slurm --time

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_crispat.csv"), emit: assignments
  path("crispat_${dataset_id}.time.txt"), optional: true, emit: timing

  publishDir "${outdir}",
             mode: 'copy',
             // Only the CSV: without this guard the .time.txt was ALSO renamed to
             // this name, and whichever copy landed last won (the timing report
             // silently replaced the assignments on Betty, 2026-09-26).
             saveAs: { filename ->
               filename.endsWith('.csv') ? "assignments_crispat_${dataset_id}.csv" : null
             }

  publishDir "${outdir}/monitoring",
             mode: 'copy',
             pattern: '*.time.txt'

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

# TIME LIMIT ENFORCED IN-BAND. Slurm enforces --time (= task.time) by killing
# the whole job, which shows up in the trace only as an ambiguous scheduler kill.
# So `timeout` fires 5 minutes earlier and exits 124 -- an unambiguous "hit the
# time limit". /usr/bin/time stays OUTSIDE the timeout so peak-RSS telemetry is
# still written when the limit fires.
# Measure peak memory & elapsed time (parity with the cleanser module)
# pin_cores.sh: single-threaded, on exactly 1 core (see nextflow.config).

/usr/bin/time -v -o crispat_${dataset_id}.time.txt \\
  timeout -k 60s ${Math.max(60, task.time.toSeconds() - 300)}s \\
  ${projectDir}/bin/pin_cores.sh 1 \\
  python "${projectDir}/bin/run_crispat.py" "${dataset_dir}/grna_matrix.h5ad"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' crispat_${dataset_id}.time.txt

"""

}
