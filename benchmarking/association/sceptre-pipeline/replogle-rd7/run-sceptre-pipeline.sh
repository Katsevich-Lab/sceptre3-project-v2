#!/bin/bash
#$ -N scep-pipe_assoc_replogle
#$ -cwd
#$ -j y
#$ -pe openmp 2
#$ -l m_mem_free=2G
export NXF_OPTS="-Xms500M -Xmx4G"

source $HOME/.research_config
nextflow pull timothy-barry/sceptre-pipeline

##########################
# REQUIRED INPUT ARGUMENTS
##########################

data_name="replogle-rd7"
assoc_dataset_name="replogle-rd7_comp_ngenes=100_ntargets=100_ncells=50k_n_nonzero_p=0.75"

data_dir=$LOCAL_BENCHMARKING_DIR"guide_assignment/input_data/"$data_name"/sceptre-pipeline/"
assoc_input_dir=$LOCAL_BENCHMARKING_DIR"association/computational/input_data/"$assoc_dataset_name"/sceptre-pipeline/"

sceptre_object_fp=$assoc_input_dir"sceptre_object.rds"

response_odm_fp=$data_dir"response.odm"
grna_odm_fp=$data_dir"grna.odm"


output_fp=$LOCAL_BENCHMARKING_DIR"association/computational/outputs/sceptre-pipeline/"$assoc_dataset_name
mkdir -p $output_fp


#################
# Invoke pipeline
#################
nextflow run timothy-barry/sceptre-pipeline -r main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_fp \
 --grna_assignment_method mixture \
  --assign_grnas_memory "2GB" \
  -with-trace "$output_fp/tracing/trace.tsv" \
  -with-report "$output_fp/tracing/report.html" \
  -with-timeline "$output_fp/tracing/timeline.html" \
  -with-dag "$output_fp/tracing/flow.dot"

