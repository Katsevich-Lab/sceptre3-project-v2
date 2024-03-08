library(sceptre)
repl_offsite <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"))
sceptre3_rd7_offsite <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "replogle-2022/rd7/")
import_dir <- paste0(repl_offsite, "processed/rd7/")
sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(import_dir, "sceptre_object.rds"),
                                                    response_odm_file_fp = paste0(import_dir, "gene.odm"),
                                                    grna_odm_file_fp = paste0(import_dir, "grna.odm"))
# specify positive control pairs only, not discovery pairs
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
sceptre_object <- set_analysis_parameters(sceptre_object,
                                          positive_control_pairs = positive_control_pairs,
                                          formula_object = formula(~ log(grna_n_nonzero + 1) + log(grna_n_umis + 1) + log(response_n_nonzero) + log(response_n_umis)))
# save the sceptre_object with analysis parameters set
write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = sceptre3_rd7_offsite)
