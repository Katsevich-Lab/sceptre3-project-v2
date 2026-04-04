#!/bin/bash
set -euo pipefail

RUN_ID="test_gasperini"
# RUN_ID="run_all_750genes"

OUT_BASE="$HOME/data/projects/sceptre3/benchmarking/association/neg-control/outputs"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR"

# Copy config file to output directory for record keeping
cp "neg-control-pipeline/configs/${RUN_ID}_config.csv" "${OUT_DIR}/"

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow \
  -C neg-control-pipeline/nextflow.config \
  run neg-control-pipeline/main.nf \
  --run_id "${RUN_ID}" \
  --out_base_dir "${OUT_BASE}" \
  -profile local \
  -with-report   "${OUT_DIR}/report.html" \
  -with-trace    "${OUT_DIR}/trace.tsv" \
  -with-timeline "${OUT_DIR}/timeline.html" \
  -with-dag      "${OUT_DIR}/dag.png"

echo "Artifacts: ${OUT_DIR}"
