#!/bin/bash
#$ -N scep-pipe_assign-grnas_replogle-rd7
#$ -cwd
#$ -j y
#$ -pe openmp 2
#$ -l m_mem_free=4G
export NXF_OPTS="-Xms500M -Xmx4G"

source $HOME/.research_config
nextflow pull timothy-barry/sceptre-pipeline

##########################
# REQUIRED INPUT ARGUMENTS
##########################

dataset="replogle-rd7"

data_directory=$LOCAL_BENCHMARKING_DIR"guide_assignment/input_data/"$dataset"/sceptre-pipeline/"
# sceptre object
sceptre_object_fp=$data_directory"sceptre_object.rds"
# response ODM
response_odm_fp=$data_directory"response.odm"
# grna ODM
grna_odm_fp=$data_directory"grna.odm"
output_fp=$LOCAL_BENCHMARKING_DIR"guide_assignment/outputs/"$dataset"/sceptre-pipeline"

#################
# Invoke pipeline
#################
nextflow run timothy-barry/sceptre-pipeline -r main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_fp \
 --grna_assignment_method mixture \
 --pipeline_stop assign_grnas \
  --assign_grnas_memory "2GB" \
  -with-trace "$output_fp/tracing/trace.tsv" \
  -with-report "$output_fp/tracing/report.html" \
  -with-timeline "$output_fp/tracing/timeline.html" \
  -with-dag "$output_fp/tracing/flow.dot"
