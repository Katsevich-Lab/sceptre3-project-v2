# makes simulated data based on Gasperini

source("~/.Rprofile")
library(tidyverse)
library(sceptre)
library(ondisc)

source(file.path(
  .get_config_path("LOCAL_SCEPTRE3_REPO_DIR"),
  "benchmarking/guide-assignment/grna-simulator/grna-simulator.R"
))

# do_transfer <- FALSE  # actually ssh this over? doesn't work in R due to OpenSSL version mismatch
run_name <- "gasperini_simulated_mid_rate_big_effects"

# 1. loading the data and preparing everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
base_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/gasperini/sceptre-pipeline/"
)


scep <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(base_fp, "sceptre_object.rds"),
  response_odm_file_fp = paste0(base_fp, "response.odm"),
  grna_odm_file_fp = paste0(base_fp, "grna.odm")
)
cell_covariates <- scep@covariate_data_frame |>
  transmute(
    grna_n_nonzero, grna_n_umis,
    batch = ifelse(prep_batch == "prep_batch_2", 1, 0),
  ) 


# determining which gRNAs to run this on [no need to run again]
grna_odm <- ondisc::initialize_odm_from_backing_file(paste0(base_fp, "grna.odm"))

# rowsums <- numeric(nrow(grna_odm))
# for(i in 1:nrow(grna_odm)) {
#   rowsums[i] <- sum(grna_odm[i,])
# }
# # getting a highly active but not maximally active gRNA
# grnas_to_consider <- rowsums > quantile(rowsums, .9) & rowsums < quantile(rowsums, .99)
# set.seed(12321)
# grnas_to_use <- which(grnas_to_consider) |> sample(5, replace = FALSE) |> sort()
# > grnas_to_use
# [1]  5034  7971  9796 10437 12669
grnas_to_use <- c(5034, 7971)#, 9796, 10437, 12669)
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
sim_results <- simulate_grna_matrix(
  grna_matrix_rows = grna_matrix_rows,
  # using pseudocount because there are 4 cells with exact 0 for each
  formula = ~ batch + log(grna_n_nonzero+1) + log(grna_n_umis+1),
  cell_covariates = cell_covariates,
  num_grnas_per_perturbation = 100,
  perturbation_effect_sizes = c(0, 1,  2, 3, 4, 5, 6), # seq(from = 0, to = log(5), length.out = 5),
  prob_perturbed = 0.008 #30 / 207324  # actual dataset had around 30 gRNAs per cell so this gives that
)



# do we get enough perturbed cells per perturbation?
# is binomial best? what about like Unif[15,45] for choosing the perturbed cells?

rr <- sim_results$perturbation_indicators |>
  rowSums()
# hist(rr, 50)
if(min(rr) < 10) {
  stop("Some gRNAs ended up with fewer than 10 cells expressing them.")
}

# ok i think this is good enough. 



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
# TODO this is hacky -- make more bulletproof later
string_for_nt <- "_effect_0.000_"  # this is used to represent exactly zero effect from the perturbation
is_nt <- grepl(pattern = string_for_nt, x = sim_grna_target_df$grna_id)
sim_grna_target_df$grna_target[is_nt] <- "non-targeting"
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



# to move to HPC
# if(do_transfer) {
#   cmd <- paste0("scp -r ", output_dir, " ", .get_config_path("REMOTE_BENCHMARKING_DIR"), "guide_assignment/input_data/", run_name)
#   system(cmd)
# }
