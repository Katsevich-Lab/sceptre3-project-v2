// Nextflow module for crispat guide assignment

process CRISPAT_ASSIGN {
    tag "${dataset_id}"
    
    conda "${moduleDir}/environment.yml"
    
    // Resources set dynamically from configs/config.csv
    cpus { resources.cpus }
    memory { resources.memory }
    
    input:
    tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
    val outdir
    
    output:
    tuple val(dataset_id), val(method), path("assignments_crispat.csv"), emit: assignments
        
    publishDir "${outdir}", 
               mode: 'copy',
               saveAs: { filename -> "assignments_crispat_${dataset_id}.csv" }
    
    script:
    """
    # Run crispat guide assignment
    python ${projectDir}/bin/run_crispat.py ${dataset_dir}/grna_counts_crispat.h5ad
    """
}