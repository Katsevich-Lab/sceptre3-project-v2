library(tidyverse)
library(sceptre)
conflicted::conflicts_prefer(dplyr::filter)
data_dir <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"), "/processed/rd7/")
project_dir <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "/replogle-2022/rd7/")
sceptre_object_fp <- paste0(project_dir, "sceptre_object_calib_check_run.rds")
response_odm_fp <- paste0(data_dir, "gene.odm")
grna_odm_fp <- paste0(data_dir, "grna.odm")
sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                    response_odm_file_fp = response_odm_fp,
                                                    grna_odm_file_fp = grna_odm_fp)

sceptre_object_cp <- sceptre_object
ex_pairs <- sceptre_object@calibration_result[1:10,] |>
  select(-p_value, -log_2_fold_change, -significant) |>
  mutate(grna_group = grna_target, grna_target = NULL)
sceptre_object_cp@negative_control_pairs <- ex_pairs
sceptre_object_cp@n_calibration_pairs <- nrow(ex_pairs)
sceptre_object_cp <- run_calibration_check_pt_2(sceptre_object_cp, output_amount = 3)
full_res <- sceptre_object_cp@calibration_result
null_dist_df <- full_res |> select(starts_with("z_null")) |> as.matrix()
hist(null_dist_df[6,], breaks = 40)

