#!/usr/bin/env nextflow
// SGE graceful-failure smoke test.
//
// Run this ON WHARTON HPC3 (see tests/run_sge_smoke.sh). It deliberately kills
// some tasks and verifies the survivors still finish -- i.e. that one method
// hitting a memory or walltime limit cannot take the whole benchmark down.
//
// It loads the REAL pipeline nextflow.config, so it exercises the actual
// errorStrategy / clusterOptions / time-directive settings, not a copy.
//
// NOTE ON MEMORY: we do NOT allocate a real memory bomb. m_mem_free is a soft
// scheduling reservation on HPC3, so over-allocating is not killed by SGE -- it
// would just destabilise a shared node. Exit 137 is what a kernel OOM kill looks
// like to the job wrapper, so we simulate it.

nextflow.enable.dsl = 2

// Completes normally. MUST survive its dying siblings.
process SURVIVOR {
  tag "${name}"
  memory '2.GB'
  time   '30m'
  publishDir "${params.testout}", mode: 'copy'
  input:  val name
  output: path "out_${name}.txt"
  script: "sleep 60; echo survived > out_${name}.txt"
}

// Real SGE walltime kill: asks for 3 minutes, then sleeps 10.
process DIE_WALLTIME {
  tag "walltime"
  memory '2.GB'
  time   '3m'
  output: path "never_written.txt"
  script: "echo 'sleeping past my h_rt...'; sleep 600; touch never_written.txt"
}

// What a kernel OOM kill looks like to the wrapper (SIGKILL -> 137).
process DIE_OOM_SIM {
  tag "oom-sim"
  memory '2.GB'
  time   '10m'
  output: path "never_written.txt"
  script: "echo 'simulating OOM kill'; exit 137"
}

// An ordinary crash inside the method itself.
process DIE_CRASH {
  tag "crash"
  memory '2.GB'
  time   '10m'
  output: path "never_written.txt"
  script: "echo 'simulating a method crash'; exit 1"
}

// Verifies the >=14GB / >=4h routing path emits -q mem.q (body is instant; this
// may sit in the queue longer than the others because of the memory request).
process ROUTING_MEM {
  tag "routing"
  memory '20.GB'
  time   '6h'
  publishDir "${params.testout}", mode: 'copy'
  output: path "out_routing.txt"
  script: "echo routed > out_routing.txt"
}

workflow {
  SURVIVOR(Channel.of('survivor1', 'survivor2'))
  DIE_WALLTIME()
  DIE_OOM_SIM()
  DIE_CRASH()
  ROUTING_MEM()
}

workflow.onComplete {
  def s = workflow.stats
  println ""
  println "=" * 62
  println "SGE SMOKE TEST -- workflow finished"
  println "  succeeded : ${s.succeededCount}   (expect 3: survivor1, survivor2, routing)"
  println "  failed    : ${s.failedCount}  (${s.ignoredCount} ignored)   (expect 3 / 3)"
  println "=" * 62
  println "Workflow-level exit is 0 by design -- errorStrategy='ignore'."
  println "run_sge_smoke.sh checks the real pass/fail conditions."
}
