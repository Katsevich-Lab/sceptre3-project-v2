#!/bin/bash
#$ -l m_mem_free=50G
#$ -N make_gasperini # job name
#$ -j y              # merge stdout and stderr
#$ -cwd              # run in current working directory
#$ -V                # export environment variables
#$ -pe smp 1  # want to only use one slot so all memory is there

# avoids unintended copying
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1



R --vanilla -e 'source("make_gasperini.R")'
