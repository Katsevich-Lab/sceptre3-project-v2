#!/bin/bash
#$ -N scep-pipeline_gasperini_set-analysis
#$ -j y
#$ -cwd
#$ -pe openmp 1              # use the existing PE
#$ -l m_mem_free=20G         # scheduler reservation (per slot)

set -euo pipefail

SIF="../../images/sceptre/sceptre.sif"

apptainer exec \
  --bind "$HOME":"$HOME" \
  "$SIF" \
  R --vanilla -e 'source("set_analysis_parameters.R")'
