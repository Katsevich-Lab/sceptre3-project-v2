// Nextflow module for sceptre guide assignment

process SCEPTRE_ASSIGN {
    tag "${dataset_id}"
    
    container "${moduleDir}/sceptre.sif"
    
    // Resources set dynamically from configs/{run_id}_config.csv
    cpus { resources.cpus }
    memory { resources.memory }
    
    input:
    tuple val(dataset_id), path(dataset_dir), val(method), val(resources)
    val outdir
    
    output:
    tuple val(dataset_id), val(method), path("assignments_sceptre.csv"), emit: assignments
        
    publishDir "${outdir}", 
               mode: 'copy',
               saveAs: { filename -> "assignments_sceptre_${dataset_id}.csv" }
    
    script:
    """
    # Run sceptre guide assignment
    Rscript ${projectDir}/bin/run_sceptre.R ${dataset_dir}/sceptre_input
    """
}