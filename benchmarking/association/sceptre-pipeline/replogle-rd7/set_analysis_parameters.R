# source("~/.Rprofile")
# library(sceptre)
# dataset <- "replogle-rd7"
# import_dir <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/", dataset, "/sceptre-pipeline/")

# sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(import_dir, "sceptre_object_initial.rds"),
#                                                     response_odm_file_fp = paste0(import_dir, "response.odm"),
#                                                     grna_odm_file_fp = paste0(import_dir, "grna.odm"))

# sceptre_object <- set_analysis_parameters(
#   sceptre_object = sceptre_object,
#   resampling_mechanism = "permutations"
#   #formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
#   #                          log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
#   #                          response_p_mito + batch + response_s_score + response_g2m_score)
# )

# write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
#                                    directory_to_write = import_dir)



source("~/.Rprofile")
library(sceptre)

dataset <- "replogle-rd7"
# this would be pulled from a _config.csv
assoc_dataset_name <- "replogle-rd7_comp_ngenes=560_ntargets=225_ncells=90k_n_nonzero_p=0.75"

if(!grepl(dataset, assoc_dataset_name)) {
  stop("`dataset` and `assoc_dataset_name` do not seem to match.")
}

# the full data is saved in guide_assignment for now
data_dir <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/", dataset, "/sceptre-pipeline/")


# loading the discovery pairs from regular sceptre for this
discovery_pairs <-  file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/computational/input_data",
  assoc_dataset_name,
  "sceptre/discovery_pairs.rds"
) |>
readRDS()

# making the directory just for this run
write_dir = file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/computational/input_data",
  assoc_dataset_name,
  "sceptre-pipeline"
)
dir.create(write_dir, showWarnings = FALSE, recursive = TRUE)



sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(data_dir, "sceptre_object_initial.rds"),
                                                    response_odm_file_fp = paste0(data_dir, "response.odm"),
                                                    grna_odm_file_fp = paste0(data_dir, "grna.odm"))

# the empty one `set_analysis_parameters()` uses by default
positive_control_pairs = data.frame(grna_target = character(0), response_id =
    character(0))

sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  formula_object = formula(~ log(response_n_nonzero + 1) + log(response_n_umis + 1) +
                             log(grna_n_nonzero + 1) + log(grna_n_umis + 1)),
  resampling_mechanism = "permutations"
)

write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = write_dir)





