library(sceptre)
repl_offsite <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"))
sceptre3_rd7_offsite <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "replogle-2022/rd7-remove-nts-add-cell-cycle/")
import_dir <- paste0(repl_offsite, "processed/rd7/remove-nts-add-cell-cycle/")
sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = paste0(import_dir, "sceptre_object.rds"),
                                                    response_odm_file_fp = paste0(import_dir, "gene.odm"),
                                                    grna_odm_file_fp = paste0(import_dir, "grna.odm"))

# Use crispat's analysis parameterss
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  side = "both", 
  grna_integration_strategy = "union", 
  control_group = 'nt_cells',
  formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
                            log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                            response_p_mito + batch + response_s_score + response_g2m_score),
  multiple_testing_alpha = 0.05
)

# save the sceptre_object with analysis parameters set
write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object,
                                   directory_to_write = sceptre3_rd7_offsite)