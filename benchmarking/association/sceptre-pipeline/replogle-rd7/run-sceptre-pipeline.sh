#!/bin/bash
#$ -N scep-pipe_assoc_replogle
#$ -cwd
#$ -j y
#$ -l m_mem_free=20G
#$ -q 16xl.q
#$ -P stat_ekatsevi_team     # project with 16xl.q access (matches profile_16xl pods)
#$ -l h_rt=24:00:00          # driver must outlive all pods; explicit so it can't inherit a short default

set -euo pipefail

export NXF_OPTS="-Xms500M -Xmx4G"

source "$HOME/.research_config"
nextflow pull jdeu1023/sceptre-pipeline

##########################
# REQUIRED INPUT ARGUMENTS
##########################

data_name="${DATA_NAME:?DATA_NAME must be set}"
assoc_dataset_name="${ASSOC_DATASET_NAME:?ASSOC_DATASET_NAME must be set}"
npairs="${NPAIRS:?NPAIRS must be set}"
num_assoc_workers="${NUM_ASSOC_WORKERS:?NUM_ASSOC_WORKERS must be set}"
pair_pod_size=$(( (npairs + num_assoc_workers - 1) / num_assoc_workers ))

data_dir="${LOCAL_BENCHMARKING_DIR}guide_assignment/input_data/${data_name}/sceptre-pipeline/"
assoc_input_dir="${LOCAL_BENCHMARKING_DIR}association/computational/input_data/${assoc_dataset_name}/sceptre-pipeline/"

sceptre_object_fp="${assoc_input_dir}sceptre_object.rds"
response_odm_fp="${data_dir}response.odm"
grna_odm_fp="${data_dir}grna.odm"
cells_to_remove_fp="${assoc_input_dir}cells_to_remove.csv"

output_fp="${LOCAL_BENCHMARKING_DIR}association/computational/outputs/sceptre-pipeline/${assoc_dataset_name}"
mkdir -p "${output_fp}" "${output_fp}/tracing"

# Capture config path before cd-ing away from the script directory
config_fp="$(pwd)/nextflow.config"

# cd into output dir so .nextflow/ cache and .nextflow.log are isolated per dataset
cd "${output_fp}"

#################
# Invoke pipeline
#################
nextflow run jdeu1023/sceptre-pipeline -r main \
  -c "${config_fp}" \
  -profile profile_16xl \
  --sceptre_object_fp "${sceptre_object_fp}" \
  --response_odm_fp "${response_odm_fp}" \
  --grna_odm_fp "${grna_odm_fp}" \
  --output_directory "${output_fp}" \
  --grna_assignment_method maximum \
  --grna_pod_size 90 \
  --pair_pod_size "${pair_pod_size}" \
  --num_assoc_workers "${num_assoc_workers}" \
  --n_calibration_pairs 0 \
  --additional_cells_to_remove "${cells_to_remove_fp}" \
  --assign_grnas_memory "2GB" \
  --assign_grnas_time_per_grna "25s" \
  --run_association_analysis_time_per_pair "1.1s" \
  --combine_association_analysis_memory "20GB" \
  --combine_association_analysis_time "30m" \
  --n_nonzero_trt_thresh 0 \
  --n_nonzero_cntrl_thresh 0 \
  --response_n_umis_range_lower 0 \
  --response_n_umis_range_upper 1 \
  --response_n_nonzero_range_lower 0 \
  --response_n_nonzero_range_upper 1 \
  --p_mito_threshold 1 \
  -with-trace "${output_fp}/tracing/trace.tsv" \
  -with-report "${output_fp}/tracing/report.html" \
  -with-timeline "${output_fp}/tracing/timeline.html" \
  -with-dag "${output_fp}/tracing/flow.dot"
