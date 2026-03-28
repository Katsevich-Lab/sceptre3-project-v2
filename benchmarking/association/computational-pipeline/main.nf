#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import method modules
include { SCEPTRE_COMPUTATIONAL } from './modules/sceptre'
include { SCEPTREGAMPOI_COMPUTATIONAL } from './modules/sceptregampoi'
include { MIXSCALE_COMPUTATIONAL } from './modules/mixscale'
include { FRPERTURB_COMPUTATIONAL } from './modules/frperturb'

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
            def data_method = (method == 'sceptregampoi') ? 'sceptre' : method
            def dataset_dir = file("${params.dataset_base_dir}/${dataset_id}/${data_method}")
            def resources = [
                cpus: resource_row.cpus,
                memory: resource_row.memory
            ]

            [dataset_id, dataset_dir, method, resources]
        }

    // Route to appropriate method based on method name
    branched_ch = dataset_method_ch.branch {
        sceptre: it[2] == 'sceptre'
        sceptregampoi: it[2] == 'sceptregampoi'
        mixscale: it[2] == 'mixscale'
        frperturb: it[2] == 'frperturb'
    }

    // Run sceptre computational benchmarking
    sceptre_results = SCEPTRE_COMPUTATIONAL(branched_ch.sceptre, outdir)

    // Run sceptregampoi computational benchmarking
    sceptregampoi_results = SCEPTREGAMPOI_COMPUTATIONAL(branched_ch.sceptregampoi, outdir)

    // Run mixscale computational benchmarking
    mixscale_results = MIXSCALE_COMPUTATIONAL(branched_ch.mixscale, outdir)

    // Run FR-Perturb computational benchmarking
    frperturb_results = FRPERTURB_COMPUTATIONAL(branched_ch.frperturb, outdir)
}

workflow.onComplete {
    println "Computational benchmarking pipeline completed!"
    println "Results in: ${params.out_base_dir}/${params.run_id}"
    println "Performance metrics in trace.txt, timeline.html, report.html"
}
