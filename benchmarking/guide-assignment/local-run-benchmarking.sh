#!/bin/bash

## this script is for running locally, i.e. on a laptop
## as opposed to on the cluster

set -euo pipefail

# RUN_ID="replogle-rd7_small"
RUN_ID="gasperini_small"
# RUN_ID="gasperini_and_replogle-rd7"

OUT_BASE="$HOME/data/projects/sceptre3/benchmarking/guide_assignment/outputs"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR" nf-logs

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow run guide-assignment-pipeline/main.nf \
  --run_id "${RUN_ID}" \
  -profile local \
  --out_base_dir "${OUT_BASE}" \
  -with-report   "${OUT_DIR}/report.html" \
  -with-trace    "${OUT_DIR}/trace.tsv" \
  -with-timeline "${OUT_DIR}/timeline.html" \
  -with-dag      "${OUT_DIR}/dag.png"

echo "Artifacts: ${OUT_DIR}"
echo "SGE stdout/err logs: $(pwd)/nf-logs/"
