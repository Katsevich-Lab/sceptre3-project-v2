#!/bin/bash
#$ -N grna-assignment
#$ -cwd
#$ -j y
#$ -l m_mem_free=2G
#$ -l h_rt=96:00:00      # the nf driver must OUTLIVE the whole pipeline; without this it
                         # resource-matches short.q (s_rt=4h) and, on dying, qdels its children.
                         # Must exceed the longest per-task `time` (default_time=48h) PLUS
                         # however long those tasks sit queued.
#$ -q hpc3.q             # long-job queue (h_rt max 8760h); driver is light, just orchestrates
#$ -pe openmp 1

set -euo pipefail

# RUN_ID="run_all_sims_medium"
# RUN_ID="run_all"
RUN_ID="test_pertpy_prealloc"
# RUN_ID="run_pertpy_gasp"
# RUN_ID="run_all_which_actually_finishes"
# RUN_ID="run_all_medium_pt2"

OUT_BASE="$(realpath -m "${LOCAL_BENCHMARKING_DIR}guide_assignment/outputs")"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR" nf-logs

# Copy config file to output directory for record keeping
cp "guide-assignment-pipeline/configs/${RUN_ID}_config.csv" "${OUT_DIR}/"

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow \
  -C ~/.nextflow/config \
  -C guide-assignment-pipeline/nextflow.config \
  run guide-assignment-pipeline/main.nf \
  --run_id "${RUN_ID}" \
  --out_base_dir "${OUT_BASE}" \
  -with-report   "${OUT_DIR}/report.html" \
  -with-trace    "${OUT_DIR}/trace.tsv" \
  -with-timeline "${OUT_DIR}/timeline.html" \
  -with-dag      "${OUT_DIR}/dag.png"



echo "Artifacts: ${OUT_DIR}"
echo "SGE stdout/err logs: $(pwd)/nf-logs/"

# Move qsub cluster log to output directory (uses job name from #$ -N directive)
if [ -n "${JOB_ID:-}" ]; then
  CLUSTER_LOG="grna-assignment.o${JOB_ID}"
  if [ -f "${CLUSTER_LOG}" ]; then
    cp "${CLUSTER_LOG}" "${OUT_DIR}/"
    echo "Cluster log moved to: ${OUT_DIR}/${CLUSTER_LOG}"
  fi
fi

