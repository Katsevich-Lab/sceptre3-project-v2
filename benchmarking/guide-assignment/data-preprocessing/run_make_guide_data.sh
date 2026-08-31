#!/bin/bash
# run_make_guide_data.sh
# Build the guide-assignment benchmark input datasets on Wharton HPC3 (SGE).
#
#   qsub run_make_guide_data.sh gasperini
#   qsub run_make_guide_data.sh replogle
#
# Direct qsub jobs are NOT routed by the ~/.nextflow/config profile, so the
# queue/project/mem/time are set explicitly here. m_mem_free>=14G + the
# stat_ekatsevi_team project both route to mem.q, which allows long h_rt.
# Right-size m_mem_free after a first run with:  qacct -j <id> | grep -i maxvmem
#
#$ -N make-guide-data
#$ -cwd
#$ -j y
#$ -P stat_ekatsevi_team
#$ -q mem.q
#$ -l m_mem_free=32G
#$ -l h_rt=8:00:00
#$ -pe openmp 1

# ============================================================================
# DORMANT (paused 2026-08-30). Part of the NESTED GUIDE-SUBSET benchmark:
# datasets of 100/200/400/800 guides selected by perturbed sample size, built to
# measure how method memory and runtime SCALE with guide count.
#
# That is a different question from the one currently being run, which is "which
# methods can handle the ENTIRE grna matrix at all" -- see build_full_h5ad.R.
# Nothing outside this directory references these files.
#
# Paused, not abandoned: this code was unit- and integration-tested and works.
# ============================================================================

set -euo pipefail

DATASET="${1:-}"
case "$DATASET" in
  gasperini) BUILD_SCRIPT="build_gasperini.R" ;;
  replogle)  BUILD_SCRIPT="build_replogle.R"  ;;
  *) echo "Usage: qsub run_make_guide_data.sh {gasperini|replogle}" >&2; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# --- R environment --------------------------------------------------------
# Adjust to however R (with ondisc, Matrix, zellkonverter) is provided on HPC3
# -- e.g. `module load R/4.x` or activating a conda env. Left as a bare
# `Rscript` on PATH by default. NOTE: the datasets are normally built locally
# (like the association data-prep); this wrapper is just an optional HPC path.
# module load R/4.4.1

echo "Host: $(hostname)   Dataset: ${DATASET}   $(date)"
Rscript --vanilla "${SCRIPT_DIR}/${BUILD_SCRIPT}"
echo "Done: ${DATASET}   $(date)"
