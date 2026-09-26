// modules/fishashplus/main.nf
//
// Runs fishash+ via bin/run_fishashplus.R. Reads the sceptre/ input subdirectory
// (grna_matrix.rds, the same matrix sceptre gets; see the input_subdir mapping in
// main.nf) and writes assignments by guide/cell NAME.

process FISHASHPLUS_ASSIGN {
  tag "${dataset_id}"

  container params.fishashplus_sif

  cpus { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> Slurm --time

  stageInMode 'symlink'

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_fishashplus.csv"), emit: assignments
  path("fishashplus_${dataset_id}.time.txt"), optional: true, emit: timing

  // Only the CSV under the assignments name (the .time.txt goes to monitoring/).
  publishDir "${outdir}",
             mode: 'copy',
             saveAs: { filename ->
               filename.endsWith('.csv') ? "assignments_fishashplus_${dataset_id}.csv" : null
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

# TIME LIMIT ENFORCED IN-BAND. Slurm enforces --time (= task.time) by killing
# the whole job, which shows up in the trace only as an ambiguous scheduler kill.
# So `timeout` fires 5 minutes earlier and exits 124 -- an unambiguous "hit the
# time limit". /usr/bin/time stays OUTSIDE the timeout so peak-RSS telemetry is
# still written when the limit fires.
# GNU time comes from the `time` apt package in fishashplus.def: this task runs
# entirely inside the container, so the host's /usr/bin/time is not reachable.
# pin_cores.sh: single-threaded, on exactly 1 core (see nextflow.config).
# fishashplus has no threading of its own; this is a guarantee, not a fix.

/usr/bin/time -v -o fishashplus_${dataset_id}.time.txt \\
  timeout -k 60s ${Math.max(60, task.time.toSeconds() - 300)}s \\
  ${projectDir}/bin/pin_cores.sh 1 \\
  Rscript --vanilla "${projectDir}/bin/run_fishashplus.R" "${dataset_dir}/grna_matrix.rds" "${dataset_id}"

# Print a one-line summary into .command.out for convenience
awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
     /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' fishashplus_${dataset_id}.time.txt
"""

}
