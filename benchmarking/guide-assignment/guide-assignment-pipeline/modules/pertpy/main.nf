// Nextflow module for pertpy guide assignment

process PERTPY_ASSIGN {
    tag "${dataset_id}"
    stageInMode 'copy'                            // ensure inputs are copied, not symlinked

    container "${moduleDir}/pertpy.sif"

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
    set -euo pipefail

    # Writable caches (matplotlib/numba)
    export MPLCONFIGDIR="\$PWD/.mplconfig";      mkdir -p "\$MPLCONFIGDIR"
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"; mkdir -p "\$NUMBA_CACHE_DIR"

    echo "PWD: \$(pwd)"
    echo "Listing dataset_dir: ${dataset_dir}"
    ls -l "${dataset_dir}" || true

    # Avoid shell $dataset_dir; use Nextflow-substituted path + PWD
    cd ${dataset_dir}
    infile="\$PWD/grna_counts_pertpy.h5ad"
    [ -f "\$infile" ] || { echo "ERROR: Not found -> \$infile"; ls -l; exit 2; }

    # Run pertpy guide assignment
    python ${projectDir}/bin/run_pertpy.py "\$infile"
    """
}

