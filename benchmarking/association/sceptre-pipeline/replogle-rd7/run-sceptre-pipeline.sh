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
assoc_dataset_name="replogle-rd7_comp_ngenes=560_ntargets=225_ncells=90k_n_nonzero_p=0.75"
# CHANGE --pair_pod_size ??

data_dir=$LOCAL_BENCHMARKING_DIR"guide_assignment/input_data/"$data_name"/sceptre-pipeline/"
assoc_input_dir=$LOCAL_BENCHMARKING_DIR"association/computational/input_data/"$assoc_dataset_name"/sceptre-pipeline/"

sceptre_object_fp=$assoc_input_dir"sceptre_object.rds"

response_odm_fp=$data_dir"response.odm"
grna_odm_fp=$data_dir"grna.odm"

cells_to_remove_fp=$assoc_input_dir"cells_to_remove.csv"


output_fp=$LOCAL_BENCHMARKING_DIR"association/computational/outputs/sceptre-pipeline/"$assoc_dataset_name
mkdir -p $output_fp


#################
# Invoke pipeline
#################
nextflow run jdeu1023/sceptre-pipeline -r main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_fp \
 --grna_assignment_method mixture \
 --pair_pod_size 5000 \
 --additional_cells_to_remove $cells_to_remove_fp \
  --assign_grnas_memory "2GB" \
  -with-trace "$output_fp/tracing/trace.tsv" \
  -with-report "$output_fp/tracing/report.html" \
  -with-timeline "$output_fp/tracing/timeline.html" \
  -with-dag "$output_fp/tracing/flow.dot"

