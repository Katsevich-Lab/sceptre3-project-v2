#!/bin/bash
# Launch the guide-assignment pipeline's Nextflow driver as a Slurm job.
#
# Usage (from benchmarking/guide-assignment/):
#   sbatch run-benchmarking.sh <RUN_ID>
# RUN_ID selects guide-assignment-pipeline/configs/<RUN_ID>_config.csv and names
# the output folder. Passing it here (rather than editing this file) keeps git
# clean; the record of each run lives with its results instead: the config CSV,
# nextflow.log, trace/report/timeline, the driver's Slurm log, and run-info.txt
# (git commit and any uncommitted changes).
#
# The driver is light (it only orchestrates; every method runs as its own Slurm
# job), but it must OUTLIVE the whole pipeline: if it dies, Nextflow cancels its
# children. So it asks for the partition maximum (7 days); tasks set their own
# --time from the config CSV. Account and partition are left to the cluster
# defaults; executor settings for the tasks come from ~/.nextflow/config.
#SBATCH -J grna-assignment
#SBATCH -c 1
#SBATCH -t 7-00:00:00
#SBATCH -o grna-assignment-%j.out

# Batch jobs don't get the interactive shell's setup, so load data paths here.
source ~/.research_config

# The driver builds the conda envs, so it needs conda. On Betty it's a module;
# elsewhere, whatever conda is already on PATH is used. (Before `set -u`:
# Betty's module scripts read unset variables.)
command -v conda >/dev/null || module load miniconda3/25.5.1

set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: sbatch run-benchmarking.sh <RUN_ID>" >&2
  exit 2
fi
RUN_ID=$1

CONFIG="guide-assignment-pipeline/configs/${RUN_ID}_config.csv"
if [ ! -f "$CONFIG" ]; then
  echo "No run config at $CONFIG" >&2
  exit 2
fi

OUT_BASE="$(realpath -m "${LOCAL_BENCHMARKING_DIR}guide_assignment/outputs")"
OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "$OUT_DIR"

# Record of what ran: the config, and the exact code version.
cp "$CONFIG" "${OUT_DIR}/"
{
  echo "run_id:     ${RUN_ID}"
  echo "date:       $(date)"
  echo "slurm_job:  ${SLURM_JOB_ID:-none}"
  echo "host:       $(hostname)"
  echo "git_branch: $(git rev-parse --abbrev-ref HEAD)"
  echo "git_commit: $(git rev-parse HEAD)"
  echo "uncommitted changes (full diff in run-info.diff if any):"
  git status --short | sed 's/^/  /'
} > "${OUT_DIR}/run-info.txt"
if [ -n "$(git status --short --untracked-files=no)" ]; then
  git diff HEAD > "${OUT_DIR}/run-info.diff"
fi

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="$PWD/.nextflow"
export NXF_APPTAINER=true
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export NXF_SINGULARITY_CMD=apptainer

nextflow \
  -log "${OUT_DIR}/nextflow.log" \
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

# Copy this driver's Slurm log next to the results
if [ -n "${SLURM_JOB_ID:-}" ]; then
  DRIVER_LOG="grna-assignment-${SLURM_JOB_ID}.out"
  if [ -f "${DRIVER_LOG}" ]; then
    cp "${DRIVER_LOG}" "${OUT_DIR}/"
    echo "Driver log copied to: ${OUT_DIR}/${DRIVER_LOG}"
  fi
fi
