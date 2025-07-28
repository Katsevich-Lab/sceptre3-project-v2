#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import method modules directly
include { CRISPAT_ASSIGN } from './modules/crispat'
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
            def dataset_dir = file("${params.dataset_base_dir}/${dataset_id}")
            def resources = [cpus: resource_row.cpus, memory: resource_row.memory]
            
            [dataset_id, dataset_dir, method, resources]
        }
    
    // Route to appropriate method based on method name
    branched_ch = dataset_method_ch.branch {
        crispat: it[2] == 'crispat'
        // TODO: Add more methods here
    }
        
    // Run crispat with explicit output directory
    crispat_results = CRISPAT_ASSIGN(branched_ch.crispat, outdir)
    // TODO: Add more methods here
}

workflow.onComplete {
    println "Benchmarking pipeline completed!"
    println "Results in: ${params.out_base_dir}/${params.run_id}"
}