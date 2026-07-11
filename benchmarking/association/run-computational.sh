#!/bin/bash
#$ -N computational
#$ -cwd
#$ -j y
#$ -l m_mem_free=2G
#$ -l h_rt=48:00:00      # the nf driver must OUTLIVE the whole pipeline; without this it
                         # resource-matches short.q (s_rt=4h) and, on dying, qdels its children
#$ -q hpc3.q             # long-job queue (h_rt max 8760h); driver is light, just orchestrates
#$ -pe openmp 1

set -euo pipefail

## these ended up being pretty small, so i'm trying bigger
# RUN_ID="run_all_max_t5"
RUN_ID="run_repl_max_upto_ng1000nt800"

OUT_BASE="${LOCAL_BENCHMARKING_DIR}/association/computational/outputs"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR" nf-logs

# Clean work directory for fresh benchmark (critical - prevents cached results)
rm -rf computational-pipeline/work

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
