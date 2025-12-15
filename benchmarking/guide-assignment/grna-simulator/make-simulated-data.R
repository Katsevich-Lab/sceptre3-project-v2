# makes simulated data based on Gasperini

rm(list=ls())
source("~/.Rprofile")
library(tidyverse)
library(sceptre)
library(ondisc)

source(file.path(
  .get_config_path("LOCAL_SCEPTRE3_REPO_DIR"),
  "benchmarking/guide-assignment/grna-simulator/grna-simulator.R"
))
# source(file.path(
#   .get_config_path("LOCAL_SCEPTRE3_REPO_DIR"),
#   "benchmarking/guide-assignment/grna-simulator/grna-simulator-poisson.R"
# ))

# do_transfer <- FALSE  # actually ssh this over? doesn't work in R due to OpenSSL version mismatch
num_cells <- 100000 # Inf to take all
source_data <- "replogle-rd7"
run_name <- paste0(source_data, "_simulated_medium")

stopifnot(source_data %in% c("gasperini", "replogle-rd7"))

# 1. loading the data and preparing everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
base_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data", source_data, "sceptre-pipeline/"
)


scep <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(base_fp, "sceptre_object.rds"),
  response_odm_file_fp = paste0(base_fp, "response.odm"),
  grna_odm_file_fp = paste0(base_fp, "grna.odm")
)
if(source_data == "gasperini") {
  cell_covariates <- scep@covariate_data_frame |>
    transmute(
      grna_n_nonzero, grna_n_umis,
      batch = ifelse(prep_batch == "prep_batch_2", 1, 0)
    ) 
} else if(source_data == "replogle-rd7") {
  cell_covariates <- scep@covariate_data_frame |>
    transmute(
      grna_n_nonzero, grna_n_umis
    ) 
}



# determining which gRNAs to run this on [no need to run again]
grna_odm <- ondisc::initialize_odm_from_backing_file(paste0(base_fp, "grna.odm"))

## only need to run this once ~~~~~~~~~~~~~~~~~~~~~
# get_active_grnas <- function(odm, row_summary_func, num_grnas = 5, seed = 12321) {
#   grna_summary_stats <- numeric(nrow(odm))
#   for(i in 1:nrow(odm)) {
#     grna_summary_stats[i] <- row_summary_func(odm[i,])
#   }
#   set.seed(seed)
#   highly_expressed <- (grna_summary_stats > quantile(grna_summary_stats, .9) &
#                          grna_summary_stats < quantile(grna_summary_stats, .99)) |>
#     which() |>
#     sample(num_grnas, replace = FALSE) |>
#     sort()
#   highly_expressed
# }
# get_active_grnas(grna_odm, function(x) sum(x > 0))
active_grnas = c(gasperini = 5034, `replogle-rd7` = 58)
grnas_to_use <- active_grnas[source_data]

grna_matrix_rows <- lapply(grnas_to_use, function(i) grna_odm[i,]) |>
  do.call(what = rbind) |>
  `dimnames<-`(list(
    rownames(grna_odm)[grnas_to_use],
    rownames(cell_covariates)
  ))

# grna_matrix_rows = grna_matrix_rows
# formula = ~ batch + log(grna_n_nonzero+1) + log(grna_n_umis+1)
# cell_covariates = cell_covariates
# num_grnas_per_perturbation = 100
# perturbation_effect_sizes = seq(from = 0, to = log(5), length.out = 5)
# prob_perturbed = 0.02


# 2. simulating the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
num_cells <- min(num_cells, nrow(cell_covariates))
grna_matrix_rows <- grna_matrix_rows[,1:num_cells, drop=FALSE]
cell_covariates <- cell_covariates[1:num_cells, ] 


sim_results <- simulate_grna_matrix(
  grna_matrix_rows = grna_matrix_rows,
  formula = ~ 1 + log(grna_n_nonzero+1) + log(grna_n_umis+1),
  cell_covariates = cell_covariates,
  num_grnas_per_perturbation = 25,
  perturbation_effect_sizes = data.frame(non_pert_effect = -4, pert_effect = 6, theta = 4.5),
  prob_perturbed = c(.001) 
)

plot_umi_hist_by_pert(
  grna_matrix_rows = grna_matrix_rows,
  sim_results = sim_results,
  grna_idx = 1, sim_idx = 1
)

# 3. saving to the appropriate places ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output_dir <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  paste0("guide_assignment/input_data/", run_name)
)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# 3.0: saving model info and true perturbation status
saveRDS(sim_results$perturbation_indicators, file = file.path(output_dir, "true_perturbation_status.rds"))
saveRDS(sim_results$model_info, file = file.path(output_dir, "model_info.rds"))
saveRDS(sim_results$sim_params, file = file.path(output_dir, "sim_params.rds"))

# 3.1: crispat
crispat_dir <- file.path(output_dir, "crispat")
dir.create(crispat_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing crispat to", crispat_dir))
R_to_h5ad(sim_results$simulated_matrix, file.path(crispat_dir, "grna_matrix.h5ad"))

# 3.2: pertpy
pertpy_dir <- file.path(output_dir, "pertpy")
dir.create(pertpy_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing pertpy to", pertpy_dir))
R_to_h5ad(sim_results$simulated_matrix, file.path(pertpy_dir, "grna_matrix.h5ad"))

# 3.3: cleanser
cleanser_dir <- file.path(output_dir, "cleanser")
dir.create(cleanser_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing cleanser to", cleanser_dir))
Matrix::writeMM(sim_results$simulated_matrix, file.path(cleanser_dir, "grna_matrix.mtx"))

# making smaller one just for CLEANSER
Matrix::writeMM(sim_results$simulated_matrix[,1:1e5], file.path(cleanser_dir, "grna_matrix.mtx"))

# 3.4: sceptre
sceptre_dir <- file.path(output_dir, "sceptre")
dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing sceptre gRNA matrix and other files to", sceptre_dir))
saveRDS(sim_results$simulated_matrix, file.path(sceptre_dir, "grna_matrix.rds"))

# 3.4.1: grna_target_data_frame
real_grna_target_df <- scep@grna_target_data_frame |>
  dplyr::filter(grna_id %in% rownames(grna_odm)[grnas_to_use])
sim_grna_target_df <- data.frame(
  grna_id = rownames(sim_results$simulated_matrix),
  grna_target = ""
)
for(i in 1:nrow(real_grna_target_df)) {
  matches_this_grna <- grepl(pattern = real_grna_target_df$grna_id[i], x = sim_grna_target_df$grna_id)
  sim_grna_target_df$grna_target[matches_this_grna] <- real_grna_target_df$grna_target[i]
}

# making NT gRNAs for sceptre
if(!any(sim_grna_target_df$grna_target == "non-targeting")) {
  sim_grna_target_df$grna_target[1:2] <- "non-targeting"
  warning("no NTs detected in `sim_grna_target_df` so first 2 are set to `non-targeting`.")
}

write.csv(sim_grna_target_df, file.path(sceptre_dir, "grna_target_data_frame.csv"), row.names=FALSE)

# 3.4.2: response_matrix
fake_response_mat <- Matrix::sparseMatrix(
  i = integer(0),  # No non-zero entries
  j = integer(0),
  x = numeric(0),
  dims = c(1L, ncol(sim_results$simulated_matrix)),
  dimnames = list(
    "FAKE_GENE_1",  # Single fake gene
    colnames(sim_results$simulated_matrix)  # Use same cell names as grna_matrix
  )
)
saveRDS(fake_response_mat, file.path(sceptre_dir, "response_matrix.rds"))


# 3.4.3: covariate_data_frame
cell_covariates |>
  dplyr::rename(true_grna_n_nonzero = grna_n_nonzero, true_grna_n_umis = grna_n_umis) |>
  write.csv(file.path(sceptre_dir, "covariate_data_frame.csv"), row.names=TRUE)

# 3.4.4 saving formula object for assign_grnas()

fmla <- readRDS(file.path(output_dir, "sim_params.rds"))$formula
saveRDS(fmla, file.path(sceptre_dir, "assign_grnas_formula.rds"))

cat("Command to transfer to HPC:\n")
cat(paste0("scp -r ", output_dir, " ", .get_config_path("REMOTE_BENCHMARKING_DIR"), "guide_assignment/input_data/", run_name, "\n"))
# to move to HPC
# if(do_transfer) {
#   cmd <- paste0("scp -r ", output_dir, " ", .get_config_path("REMOTE_BENCHMARKING_DIR"), "guide_assignment/input_data/", run_name)
#   system(cmd)
# }
