#!/bin/bash
set -euo pipefail

# Select which run to execute
# RUN_ID="test_100genes"
# RUN_ID="small_500genes"
RUN_ID="test_all_small"

OUT_BASE="$HOME/data/projects/sceptre3/benchmarking/association/computational/outputs"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR"

# Copy config file to output directory for record keeping
cp "computational-pipeline/configs/${RUN_ID}_config.csv" "${OUT_DIR}/"

# Nextflow environment setup
export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow \
  -C computational-pipeline/nextflow.config \
  run computational-pipeline/main.nf \
  --run_id "${RUN_ID}" \
  -profile local \
  -with-report   "${OUT_DIR}/report.html" \
  -with-trace    "${OUT_DIR}/trace.txt" \
  -with-timeline "${OUT_DIR}/timeline.html" \
  -with-dag      "${OUT_DIR}/dag.png"

echo "Computational benchmarking complete!"
echo "Results: ${OUT_DIR}"
echo "Trace file: ${OUT_DIR}/trace.txt"
