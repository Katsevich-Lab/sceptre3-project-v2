#!/bin/bash
# Limit NF driver to 4 GB memory
source $HOME/.research_config
nextflow pull timothy-barry/sceptre-pipeline
export NXF_OPTS="-Xms500M -Xmx4G"

##########################
# REQUIRED INPUT ARGUMENTS
##########################
data_directory=$LOCAL_REPLOGLE_2022_DATA_DIR"/processed/rd7/"
project_directory=$LOCAL_SCEPTRE3_DATA_DIR"/replogle-2022/rd7/"
# sceptre object
sceptre_object_fp=$project_directory"sceptre_object.rds"
# response ODM
response_odm_fp=$data_directory"gene.odm"
# grna ODM
grna_odm_fp=$data_directory"grna.odm"

##################
# OUTPUT DIRECTORY
##################
output_directory=$HOME"/sceptre_outputs"

#################
# Invoke pipeline
#################
# nextflow run timothy-barry/sceptre-pipeline -r main \
nextflow run "/Users/timbarry/research_code/sceptre-pipeline/main.nf" \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_directory \
 --grna_assignment_method mixture \
 --pair_pod_size 500 \
 --grna_pod_size 25 \
 --response_n_nonzero_range_lower 0.07 \
 --trial \
 --pipeline_stop run_qc \
 -resume \