#$ -pe openmp 2
#$ -l m_mem_free=4G
export NXF_OPTS="-Xms500M -Xmx4G"

source $HOME/.research_config
nextflow pull timothy-barry/sceptre-pipeline

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
nextflow run timothy-barry/sceptre-pipeline -r main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_directory \
 --response_n_nonzero_range_lower 0.07 \
 --grna_assignment_method thresholding \
 --threshold 10 \
 --n_calibration_pairs 10000000 \
 --calibration_group_size 1 \
 --n_nonzero_trt_thresh 1 \
 --n_nonzero_cntrl_thresh 1 \
 --pipeline_stop run_calibration_check \
 --pair_pod_size 100000 
