#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import method modules directly
include { CRISPAT_ASSIGN } from './modules/crispat'
include { CLEANSER_ASSIGN } from './modules/cleanser'
include { PERTPY_ASSIGN } from './modules/pertpy'
include { SCEPTRE_ASSIGN } from './modules/sceptre'
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
            def dataset_dir = file("${params.dataset_base_dir}/${dataset_id}/${method}")
            def resources = [
                cpus: resource_row.cpus,
                memory: resource_row.memory,
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
    // TODO: Add more methods here
}

workflow.onComplete {
    println "Benchmarking pipeline completed!"
    println "Results in: ${params.out_base_dir}/${params.run_id}"
}