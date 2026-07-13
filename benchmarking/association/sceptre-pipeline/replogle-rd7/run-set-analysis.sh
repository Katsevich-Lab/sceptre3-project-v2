#!/bin/bash
#$ -N scep-pipeline_replogle-rd7_set-analysis
#$ -j y
#$ -cwd
#$ -pe openmp 1              # use the existing PE
#$ -l m_mem_free=20G         # scheduler reservation (per slot)
#$ -q 16xl.q                 # off mem.q (busy w/ mixscale); 20G proven OK here (same as the pipeline head job)
#$ -P stat_ekatsevi_team     # project with 16xl.q access (matches profile_16xl pods)
#$ -l h_rt=8:00:00           # builds the parameterized sceptre object; explicit walltime for safety

set -euo pipefail

R --vanilla -e 'source("set_analysis_parameters.R")'
