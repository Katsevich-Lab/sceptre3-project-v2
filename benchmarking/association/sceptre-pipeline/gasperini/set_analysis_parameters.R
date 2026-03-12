source("~/.Rprofile")
library(sceptre)
dataset <- "gasperini"
import_dir <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/", dataset, "/sceptre-pipeline/")

sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(import_dir, "sceptre_object_initial.rds"),
                                                    response_odm_file_fp = paste0(import_dir, "response.odm"),
                                                    grna_odm_file_fp = paste0(import_dir, "grna.odm"))

positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(sceptre_object, distance_threshold = 1e5)

sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  formula_object = formula(~ log(response_n_nonzero + 1) + log(response_n_umis + 1) +
                             log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                             response_p_mito + prep_batch),
  resampling_mechanism = "permutations"
)

write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = import_dir)




