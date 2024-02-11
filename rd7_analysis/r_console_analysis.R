library(sceptre)

data_directory <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"), "/processed/rd7/")
project_directory <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "/replogle-2022/rd7/")
sceptre_object_fp <- paste0(project_directory, "sceptre_object.rds")
response_odm_fp <- paste0(data_directory, "gene.odm")
grna_odm_fp <- paste0(data_directory, "grna.odm")

# load the sceptre object
sceptre_object <- read_ondisc_backed_sceptre_object(sceptre_object_fp = sceptre_object_fp,
                                                    response_odm_file_fp = response_odm_fp,
                                                    grna_odm_file_fp = grna_odm_fp)
disc_pairs <- construct_trans_pairs(sceptre_object) |> dplyr::sample_n(100)

# set analysis params
sceptre_object <- set_analysis_parameters(sceptre_object = sceptre_object,
                                          discovery_pairs = disc_pairs,
                                          side = "both")
