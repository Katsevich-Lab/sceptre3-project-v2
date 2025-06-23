library(sceptre)
library(dplyr)
library(ggplot2)

# Read configuration and set paths
sceptre3_rd7_offsite <- paste0(.get_config_path("LOCAL_SCEPTRE3_DATA_DIR"), "replogle-2022/rd7-cell-cycle/")

# Read calibration check results
calibration_results <- readRDS(paste0(sceptre3_rd7_offsite, "results_run_calibration_check.rds"))

# tabulate worst offender gRNAs
calibration_results |> 
  filter(significant) |> 
  count(grna_target) |>
  arrange(desc(n))
