// Nextflow module for CLEANSER guide assignment

process CLEANSER_ASSIGN {
    tag "${dataset_id}"
    
    conda "${moduleDir}/environment.yml"
    
    // Resources set dynamically from configs/{run_id}_config.csv
    cpus { resources.cpus }
    memory { resources.memory }
    
    input:
    tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
    val outdir
    
    output:
    tuple val(dataset_id), val(method), path("assignments_cleanser.csv"), emit: assignments
        
    publishDir "${outdir}", 
               mode: 'copy',
               saveAs: { filename -> "assignments_cleanser_${dataset_id}.csv" }
    
    script:
    """
    # Run CLEANSER guide assignment
    python ${projectDir}/bin/run_cleanser.py ${dataset_dir}/grna_counts_cleanser.mtx
    """
}