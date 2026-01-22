#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import method modules
include { SCEPTRE_POSCTRL } from './modules/sceptre'
include { MIXSCALE_POSCTRL } from './modules/mixscale'
include { FRPERTURB_POSCTRL } from './modules/frperturb'

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
            def dataset_dir = file("${params.dataset_base_dir}/${dataset_id}/${method}")
            def resources = [
                cpus: resource_row.cpus,
                memory: resource_row.memory
            ]

            [dataset_id, dataset_dir, method, resources]
        }

    // Route to appropriate method based on method name
    branched_ch = dataset_method_ch.branch {
        sceptre: it[2] == 'sceptre'
        mixscale: it[2] == 'mixscale'
        frperturb: it[2] == 'frperturb'
    }

    // Run sceptre positive control analysis
    sceptre_results = SCEPTRE_POSCTRL(branched_ch.sceptre, outdir)

    // Run mixscale positive control analysis
    mixscale_results = MIXSCALE_POSCTRL(branched_ch.mixscale, outdir)

    // Run FR-Perturb positive control analysis
    frperturb_results = FRPERTURB_POSCTRL(branched_ch.frperturb, outdir)
}

workflow.onComplete {
    println "Positive control benchmarking pipeline completed!"
    println "Results in: ${params.out_base_dir}/${params.run_id}"
}
