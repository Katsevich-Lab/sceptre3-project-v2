source("~/.Rprofile")
library(sceptre)

dataset <- Sys.getenv("DATA_NAME")
if (dataset == "") stop("DATA_NAME environment variable not set.")

assoc_dataset_name <- Sys.getenv("ASSOC_DATASET_NAME")
if (assoc_dataset_name == "") stop("ASSOC_DATASET_NAME environment variable not set.")

if(!grepl(dataset, assoc_dataset_name)) {
  stop("`dataset` and `assoc_dataset_name` do not seem to match.")
}

# the full data is saved in guide_assignment for now
data_dir <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/", dataset, "/sceptre-pipeline/")

# loading the discovery pairs from regular sceptre for this
discovery_pairs <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/computational/input_data",
  assoc_dataset_name,
  "sceptre/discovery_pairs.rds"
) |>
  readRDS()

# making the directory just for this run
write_dir <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/computational/input_data",
  assoc_dataset_name,
  "sceptre-pipeline"
)
dir.create(write_dir, showWarnings = FALSE, recursive = TRUE)

sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp      = paste0(data_dir, "sceptre_object_initial.rds"),
  response_odm_file_fp   = paste0(data_dir, "response.odm"),
  grna_odm_file_fp       = paste0(data_dir, "grna.odm")
)

sceptre_object <- set_analysis_parameters(
  sceptre_object       = sceptre_object,
  discovery_pairs      = discovery_pairs,
  formula_object       = formula(~ log(response_n_nonzero + 1) + log(response_n_umis + 1) +
                                   log(grna_n_nonzero + 1) + log(grna_n_umis + 1) + prep_batch),
  resampling_mechanism = "crt"
)

write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = write_dir)
