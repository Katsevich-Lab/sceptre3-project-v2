process CLEANSER_ASSIGN {
  label 'cleanser'
  tag "${dataset_id}"

  conda "${moduleDir}/environment.yml"

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> Slurm --time

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

# Cache directory. NOTE: there is deliberately no CMDSTANPY_CACHE_DIR here --
# cmdstanpy does not read that variable (it appears nowhere in the package), so
# setting it only created an empty unused dir per task and implied, wrongly, that
# compiled Stan models were task-local. They are not: CmdStan lives in the conda
# env (CMDSTAN is set by its activate.d script) and a compiled model is written
# next to its .stan file in site-packages/cleanser/, so it PERSISTS in the env.
export XDG_CACHE_HOME="\$PWD/.cache"
mkdir -p "\$XDG_CACHE_HOME"

# Stan scratch goes in the task work dir, NOT the node-local default TMPDIR.
# cmdstanpy writes one CSV per chain holding 3 length-N arrays per draw, so a
# guide with N cells costs ~1000*3*N numbers per chain, x4 chains. The largest
# replogle warm-up guide (N=112,007) is ~27 GB on its own, which overruns
# /mnt/tmp on short.q and kills all 4 chains with exit 1.
export TMPDIR="\$PWD/stan_tmp"
mkdir -p "\$TMPDIR"

# TIME LIMIT ENFORCED IN-BAND. Slurm enforces --time (= task.time) by killing
# the whole job, which shows up in the trace only as an ambiguous scheduler kill.
# So `timeout` fires 5 minutes earlier and exits 124 -- an unambiguous "hit the
# time limit". /usr/bin/time stays OUTSIDE the timeout so peak-RSS telemetry is
# still written when the limit fires.
# Measure peak memory & elapsed time for cleanser
# pin_cores.sh: EXACTLY 4 cores, one per MCMC chain, however many the task holds
# (cmdstanpy runs min(node cores, 4) chains at once; each chain is single-threaded
# via the env scope in nextflow.config). Fails if the task has fewer than 4 cores.
/usr/bin/time -v -o cleanser_${dataset_id}.time.txt \\
  timeout -k 60s ${Math.max(60, task.time.toSeconds() - 300)}s \\
  ${projectDir}/bin/pin_cores.sh 4 \\
  python "${projectDir}/bin/run_cleanser.py" "${dataset_dir}/grna_matrix.mtx" "${dataset_id}"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' cleanser_${dataset_id}.time.txt
"""

}
