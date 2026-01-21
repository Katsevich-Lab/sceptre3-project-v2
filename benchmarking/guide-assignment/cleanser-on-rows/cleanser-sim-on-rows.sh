#!/bin/bash
#$ -N cleanser-on-rows
#$ -cwd                     # run in current working directory
#$ -j y                     # join stdout/stderr
#$ -l m_mem_free=8G         # RAM per task
#$ -t 1-25                 # parallel tasks

source ~/.bashrc

path_to_rows="/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7_simulated/cleanser/separate_rows/"
conda activate cleanser-env
Rscript cleanser-run-on-row.R "${SGE_TASK_ID}" "${path_to_rows}"
