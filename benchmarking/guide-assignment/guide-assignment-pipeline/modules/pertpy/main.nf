// Nextflow module for pertpy guide assignment (simple + strict)

// NOTE: pertpy runs single-threaded (pin_cores.sh 1), but `cpus` still comes
// from the run config: on Betty, memory is allocated per CPU, so a large-memory
// task may need several CPUs even though it uses only one.

process PERTPY_ASSIGN {
  label 'pertpy'
  tag "${dataset_id}"
  stageInMode 'symlink'
  conda "${moduleDir}/environment.yml"

  cpus   { resources.cpus }
  memory { resources.memory }
  time   { resources.time }   // -> Slurm --time

  // Store gpu settings in task.ext for access in nextflow.config
  ext.gpu_queue = { resources.gpu_queue }
  ext.gpu_time = { resources.gpu_time }

  input:
  tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
  val outdir

  output:
  tuple val(dataset_id), val(method), path("assignments_pertpy.csv"), emit: assignments
  path("pertpy_${dataset_id}.time.txt"), optional: true, emit: timing

  publishDir "${outdir}",
             mode: 'copy',
             // Only the CSV: without this guard the .time.txt was ALSO renamed to
             // this name, and whichever copy landed last won (see crispat module).
             saveAs: { filename ->
               filename.endsWith('.csv') ? "assignments_pertpy_${dataset_id}.csv" : null
             }

  publishDir "${outdir}/monitoring",
             mode: 'copy',
             pattern: '*.time.txt'

  script:
  """
  set -euo pipefail

  export JAX_PLATFORMS=cpu
  export JAX_ENABLE_X64=0
  export NUMBA_CACHE_DIR=${projectDir}/.numba_cache
  export MPLCONFIGDIR=${projectDir}/.mplconfig
  export PYTHONNOUSERSITE=1

  # TIME LIMIT ENFORCED IN-BAND. Slurm enforces --time (= task.time) by killing
  # the whole job, which shows up in the trace only as an ambiguous scheduler kill.
  # So `timeout` fires 5 minutes earlier and exits 124 -- an unambiguous "hit the
  # time limit". /usr/bin/time stays OUTSIDE the timeout so peak-RSS telemetry is
  # still written when the limit fires.
  # Run pertpy guide assignment, measuring peak memory & elapsed time
  # (parity with the cleanser module)
  # pin_cores.sh: single-threaded, on exactly 1 core (see nextflow.config).
  /usr/bin/time -v -o pertpy_${dataset_id}.time.txt \\
    timeout -k 60s ${Math.max(60, task.time.toSeconds() - 300)}s \\
    ${projectDir}/bin/pin_cores.sh 1 \\
    python ${projectDir}/bin/run_pertpy.py "${dataset_dir}/grna_matrix.h5ad" ${dataset_id}

  # Print a one-line summary into .command.out for convenience
  awk '/Maximum resident set size/ {printf "Peak RAM: %.2f GiB\\n", \$NF/1024/1024} \\
       /Elapsed \\(wall clock\\) time/ {print "Elapsed:", \$0}' pertpy_${dataset_id}.time.txt

  """
}
