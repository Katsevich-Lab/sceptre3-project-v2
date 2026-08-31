#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import method modules directly
include { CRISPAT_ASSIGN } from './modules/crispat'
include { CLEANSER_ASSIGN } from './modules/cleanser'
include { PERTPY_ASSIGN } from './modules/pertpy'
include { SCEPTRE_ASSIGN } from './modules/sceptre'
include { SCRIPT_ASSIGN } from './modules/script'
// TODO: Add more methods as needed

workflow {
    // Ensure run_id is provided
    if (!params.run_id) {
        error "Please specify --run_id parameter"
    }
    
    // Compute output directory for this run
    outdir = "${params.out_base_dir}/${params.run_id}"
    
    // Derive config file from run_id
    config_file = "${projectDir}/configs/${params.run_id}_config.csv"
    
    // Check that config file exists
    if (!file(config_file).exists()) {
        error "Config file not found: ${config_file}"
    }
    
    // Create method-dataset combinations from run-specific config
    dataset_method_ch = Channel.fromPath(config_file)
        .splitCsv(header: true)
        .map { resource_row ->
            def dataset_id = resource_row.dataset
            def method = resource_row.method
            // All script_* variants consume the same R input format as sceptre,
            // so they share the sceptre/ input subdirectory.
            def input_subdir = method.startsWith('script_') ? 'sceptre' : method
            def dataset_dir = file("${params.dataset_base_dir}/${dataset_id}/${input_subdir}")
            // `time` drives BOTH the SGE -l h_rt AND the queue routing in
            // ~/.nextflow/config (which sends time >= 4h to hpc3.q). Without it a
            // sub-14GB task resource-matches onto short.q and is killed at 4h --
            // which would look like "the method can't handle this dataset" when
            // it is really a scheduling artifact. Optional column; defaults to
            // params.default_time for older config CSVs that lack it.
            def resources = [
                cpus: resource_row.cpus,
                memory: resource_row.memory,
                time: resource_row.time ?: params.default_time,
                gpu_queue: resource_row.gpu_queue ?: '',
                gpu_time: resource_row.gpu_time ?: ''
            ]

            [dataset_id, dataset_dir, method, resources]
        }
    
    // Route to appropriate method based on method name
    branched_ch = dataset_method_ch.branch {
        crispat: it[2] == 'crispat'
        cleanser: it[2] == 'cleanser'
        pertpy: it[2] == 'pertpy'
        sceptre: it[2] == 'sceptre'
        script: it[2].startsWith('script_')
        // TODO: Add more methods here
    }
        
    // Run crispat with explicit output directory
    crispat_results = CRISPAT_ASSIGN(branched_ch.crispat, outdir)
    
    // Run cleanser with explicit output directory
    cleanser_results = CLEANSER_ASSIGN(branched_ch.cleanser, outdir)
    
    // Run pertpy with explicit output directory
    pertpy_results = PERTPY_ASSIGN(branched_ch.pertpy, outdir)
    
    // Run sceptre with explicit output directory
    sceptre_results = SCEPTRE_ASSIGN(branched_ch.sceptre, outdir)

    // Run script_* variants (all share the SCRIPT_ASSIGN process; the wrapper
    // dispatches to bin/script/{variant}.R based on the method name)
    script_results = SCRIPT_ASSIGN(branched_ch.script, outdir)
    // TODO: Add more methods here
}

workflow.onComplete {
    // errorStrategy='ignore' means this workflow reports success even if every
    // method died, so summarise task outcomes explicitly rather than letting the
    // exit code imply everything worked.
    def s = workflow.stats
    println ""
    println "=" * 62
    println "Benchmarking pipeline completed."
    // NB: ignoredCount is a SUBSET of failedCount -- an ignored task is still a
    // failed task -- so report them together rather than as two separate totals.
    println "  succeeded : ${s.succeededCount}"
    println "  failed    : ${s.failedCount}  (${s.ignoredCount} ignored, run continued)"
    println "=" * 62
    if (s.failedCount > 0) {
        println "One or more methods did not finish. Per-task status, exit code and"
        println "peak memory are in the trace file; for a task the scheduler killed,"
        println "look up its native_id with: qacct -j <native_id> | grep -iE 'maxvmem|exit_status|failed'"
    }
    println "Results in: ${params.out_base_dir}/${params.run_id}"
}