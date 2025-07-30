// Nextflow module for pertpy guide assignment

process PERTPY_ASSIGN {
    tag "${dataset_id}"
    
    conda "${moduleDir}/environment.yml"
    
    // Resources set dynamically from configs/{run_id}_config.csv
    cpus { resources.cpus }
    memory { resources.memory }
    
    input:
    tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
    val outdir
    
    output:
    tuple val(dataset_id), val(method), path("assignments_pertpy.csv"), emit: assignments
        
    publishDir "${outdir}", 
               mode: 'copy',
               saveAs: { filename -> "assignments_pertpy_${dataset_id}.csv" }
    
    script:
    """
    # Run pertpy guide assignment
    python ${projectDir}/bin/run_pertpy.py ${dataset_dir}/grna_counts_pertpy.h5ad
    """
}