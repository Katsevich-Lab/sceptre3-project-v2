#!/bin/bash
#$ -N scep-pipeline_replogle-rd7_set-analysis
#$ -j y
#$ -cwd
#$ -q 16xl.q
#$ -pe openmp 1              # use the existing PE
#$ -l m_mem_free=20G         # scheduler reservation (per slot)

set -euo pipefail

R --vanilla -e 'source("set_analysis_parameters.R")'
