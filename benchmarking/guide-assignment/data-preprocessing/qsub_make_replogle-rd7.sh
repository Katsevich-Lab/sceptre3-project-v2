#!/bin/bash
#$ -N make_replogle
#$ -j y
#$ -cwd
#$ -V
#$ -pe openmp 1              # use the existing PE
#$ -l m_mem_free=50G         # scheduler reservation (per slot)
#$ -l h_vmem=55G             # hard per-slot cap; job is killed if exceeded

set -euo pipefail

# keep libraries single-threaded
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

R --vanilla -e 'source("make_replogle-rd7.R")'

