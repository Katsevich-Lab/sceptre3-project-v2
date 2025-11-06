source("~/.Rprofile")
library(sceptre)
dataset <- "replogle-rd7"
import_dir <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/", dataset, "/sceptre-pipeline/")

sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(import_dir, "sceptre_object_initial.rds"),
                                                    response_odm_file_fp = paste0(import_dir, "response.odm"),
                                                    grna_odm_file_fp = paste0(import_dir, "grna.odm"))

sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  resampling_mechanism = "permutations"
  #formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
  #                          log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
  #                          response_p_mito + batch + response_s_score + response_g2m_score)
)

write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = import_dir)




