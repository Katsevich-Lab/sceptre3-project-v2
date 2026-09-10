// modules/fishash/main.nf
//
// Runs fishash via bin/run_fishash.R. Reads the cleanser/ input subdirectory --
// fishash() takes a bare count matrix, which is what that Matrix Market file
// already holds (see the input_subdir mapping in main.nf).

process FISHASH_ASSIGN {
  tag "${dataset_id}"

  container params.fishash_sif

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> SGE -l h_rt AND queue routing (>=4h -> hpc3.q)

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_fishash.csv"), emit: assignments
  path("fishash_${dataset_id}.time.txt"), optional: true, emit: timing

  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { filename ->
               filename.endsWith('.csv') ? "assignments_fishash_${dataset_id}.csv" : null
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

# R needs writable temp and (if any package tries) a user lib dir
export TMPDIR="\$PWD/tmp";           mkdir -p "\$TMPDIR"
export R_TMPDIR="\$TMPDIR"
export R_USER="\$PWD"
export R_LIBS_USER="\$PWD/.Rlibs";   mkdir -p "\$R_LIBS_USER"

# TIMEOUT ENFORCED IN-BAND, not by the scheduler. Measured on HPC3 2026-08-31:
# a task requesting `-l h_rt=00:03:00` ran for 10m and exited 0, so SGE does NOT
# enforce h_rt here (h_rt is requestable and short.q caps at 04:05:00, yet the
# limit never fired). On mem.q -- where every full-dataset task lands -- the
# ceiling is effectively a year, so an unbounded method would run until the
# nextflow driver dies and qdel's it. `timeout` exits 124, which reaches the
# trace as an unambiguous "hit the time limit". /usr/bin/time stays OUTSIDE the
# timeout so peak-RSS telemetry is still written when the limit fires.
# GNU time comes from the `time` apt package in fishash.def: unlike the
# conda-backed methods this task runs entirely inside the container, so the
# host's /usr/bin/time is not reachable.
/usr/bin/time -v -o fishash_${dataset_id}.time.txt \\
  timeout -k 60s ${task.time.toSeconds()}s \\
  Rscript --vanilla "${projectDir}/bin/run_fishash.R" "${dataset_dir}/grna_matrix.mtx" "${dataset_id}"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' fishash_${dataset_id}.time.txt
"""

}
