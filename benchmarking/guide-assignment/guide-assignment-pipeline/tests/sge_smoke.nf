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

// Hits the IN-BAND timeout, which is how the real modules now bound runtime.
//
// This used to ask SGE for h_rt=3m and sleep 600, expecting a scheduler kill.
// Measured on HPC3 2026-08-31, that task ran the full 10m and exited 0: SGE does
// NOT enforce h_rt here (it is requestable, and short.q caps at 04:05:00, yet
// the limit never fired). So the methods enforce their own limit with coreutils
// `timeout`, and this reproduces that mechanism. Expect exit 124.
process DIE_TIMEOUT {
  tag "timeout"
  memory '2.GB'
  time   '10m'
  output: path "never_written.txt"
  script: "echo 'about to exceed my in-band timeout...'; timeout -k 10s 5s sleep 600"
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
  memory '16.GB'   // just over the 14GB mem.q routing threshold -- no more
  time   '6h'      // than needed, since this reserves a slice of a shared node
  publishDir "${params.testout}", mode: 'copy'
  output: path "out_routing.txt"
  script: "echo routed > out_routing.txt"
}

workflow {
  SURVIVOR(Channel.of('survivor1', 'survivor2'))
  DIE_TIMEOUT()
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
