#!/bin/bash
#$ -N computational
#$ -cwd
#$ -j y
#$ -l m_mem_free=2G
#$ -pe openmp 1

set -euo pipefail

RUN_ID="run_repl_cells_50k_100k"

OUT_BASE="${LOCAL_BENCHMARKING_DIR}/association/computational/outputs"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR" nf-logs

# Copy config file to output directory for record keeping
cp "computational-pipeline/configs/${RUN_ID}_config.csv" "${OUT_DIR}/"

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow \
  -C ~/.nextflow/config \
  -C computational-pipeline/nextflow.config \
  run computational-pipeline/main.nf \
  --run_id "${RUN_ID}" \
  --out_base_dir "${OUT_BASE}" \
  -anew \
  -with-report   "${OUT_DIR}/report.html" \
  -with-trace    "${OUT_DIR}/trace.txt" \
  -with-timeline "${OUT_DIR}/timeline.html" \
  -with-dag      "${OUT_DIR}/dag.png"



echo "Artifacts: ${OUT_DIR}"
echo "SGE stdout/err logs: $(pwd)/nf-logs/"

# Move qsub cluster log to output directory (uses job name from #$ -N directive)
if [ -n "${JOB_ID:-}" ]; then
  CLUSTER_LOG="computational.o${JOB_ID}"
  if [ -f "${CLUSTER_LOG}" ]; then
    cp "${CLUSTER_LOG}" "${OUT_DIR}/"
    echo "Cluster log moved to: ${OUT_DIR}/${CLUSTER_LOG}"
  fi
fi
