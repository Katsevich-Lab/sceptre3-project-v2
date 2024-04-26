library(arrow)
library(dplyr)
library(katlabutils)

conflicts_prefer(dplyr::filter)
res_dir <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "nf_pipelines/gasp_trans/sceptre_outputs/trans_results")
ds <- open_dataset(res_dir)
result_1 <- read_parquet(paste0(res_dir, "/result_1.parquet"))
unique_nt_grnas <- grep(pattern = "^nt_", x = unique(result_1$grna_target), value = TRUE)

nt_res <- ds |>
    filter(grna_target == unique_nt_grnas[1]) |> 
    select(p_value, log_2_fold_change, response_id) |>
    arrange(p_value) |>
    collect()

head(nt_res)
