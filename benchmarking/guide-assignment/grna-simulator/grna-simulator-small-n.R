rm(list=ls())
source("~/.Rprofile")
library(tidyverse)
library(sceptre)
library(ondisc)
library(reticulate)
library(SingleCellExperiment)
library(zellkonverter)

write_to_h5ad <- function(mat, path_to_write) {
  env_name <- "r-anndata"
  curr_envs <- conda_list()$name
  if(!env_name %in% curr_envs) {
    conda_create("r-anndata", packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"),
               envname = "r-anndata", pip = TRUE)
  }
  use_condaenv("r-anndata", required = TRUE)
  py_config()
  
  sce <- SingleCellExperiment(assays = list(counts = mat))
  writeH5AD(sce, path_to_write)
}

plot_results <- function(grna_real, grna_sim, is_pert, sim_name) {
  
  logp1_trans <- function(base = 10) {
    trans <- function(x) log(x + 1, base)
    inv   <- function(x) base^x - 1  # correct inverse
    
    # Major breaks: at base^k - 1
    brks <- function(limits) {
      rng <- limits + 1
      scales::log_breaks(base = base)(rng) - 0
    }
    
    scales::trans_new(
      name = paste0("logp1-", base),
      transform = trans,
      inverse   = inv,
      breaks    = brks,
      domain = c(0, Inf)
    )
  }

  plt1 <- data.frame(
    umi = c(grna_real, grna_sim),
    group = rep(c("real", "sim"), each = length(grna_real))
  ) |>
    ggplot(aes(x = umi, fill = group)) +
    geom_histogram(bins=50) +
    facet_wrap(~group) +
    theme_bw()  +
    scale_y_continuous(trans = logp1_trans(10)) +
    scale_x_continuous(trans = logp1_trans(10)) +
    ggtitle(sim_name)
  
  plt2 <- data.frame(
    umi = grna_sim,
    pert_status = ifelse(is_pert == 1, "pert", "not pert")
  ) |>
    ggplot(aes(x = umi, fill = pert_status)) +
    geom_histogram(bins=50) +
    facet_wrap(~pert_status) +
    theme_bw()  +
    scale_y_continuous(trans = logp1_trans(10)) +
    scale_x_continuous(trans = logp1_trans(10)) +
    ggtitle("Simulated UMI counts by pert status")
  
  # printing some numeric comparisons
  comparison_df <- data.frame(
    real = c(max(grna_real),   mean(grna_real == 0),    var(grna_real),    mean(grna_real)),
    sim = c(max(grna_sim),     mean(grna_sim==0),      var(grna_sim),      mean(grna_sim))
  ) |>
    `rownames<-`(c("max", "frac_0", "var", "mean"))
  print(comparison_df)
  
  cowplot::plot_grid(plt1, plt2, nrow = 2)
}

# `sim_params` details: a list with 5 components, as described below:
# sim_params = list(
#   logit_prob_pert_increase = 2,  # intercept adjustment on link scale for P(perturbed)
#   non_pert_effect = -5, # link-scale intercept adjustment for non-perturbed NB dist mean
#   pert_effect = 5.25, # link-scale intercept adjustment for perturbed NB dist mean
#   theta_pert = 10, # theta to use for perturbed cells
#   theta_non_pert = .1 # theta to use for non-perturbed cells
# )
make_simulated_umi_counts <- function(
    umi_counts, # a vector of length n giving the actual observed UMI counts for the guide we are modeling
    num_guides_to_simulate, # int, the number of iid simulated copies to make (the simulated data will have this many rows)
    cell_covariates,  # a data.frame of length n giving cell covariates to regress on
    sceptre_assignments, # a boolean vector of length n indicating which cells are considered perturbed for the real guide in question
    formula, # a one-sided formula of the cell covariates
    num_cells, # the number of cells to simulate
    sim_params, # a list with the simulation parameters; see above
    guide_name, # the real guide name that was used
    clip_to_real_max = TRUE, # any simulated UMI counts above `max(umi_counts)` are set to that value
    seed=NULL
) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  mod_data <- cell_covariates |> mutate(umi_counts, sceptre_assignments)
  max_real_umi <- max(umi_counts)
  # 1. determine which cells we are using
  expressed_cells <- which(sceptre_assignments)
  non_expressed_cells <- setdiff(1:length(sceptre_assignments), expressed_cells)
  
  cell_idx_subsampled <- c(
    expressed_cells,
    sample(non_expressed_cells, num_cells - length(expressed_cells))
  ) |> sort()
  
  cat("Beginning to fit models: ")
  
  # 2. get the effect of the covariates on the mean of the UMI counts
  fmla_umi <- as.formula(paste("umi_counts ~", as.character(formula)[2]))
  mod_umi <- MASS::glm.nb(fmla_umi, data = mod_data)
  log_mean_umi <- fitted(mod_umi)[cell_idx_subsampled] |> log()
  
  # 3. get the effect of the covariates on the perturbation probability
  fmla_pert <- as.formula(paste("sceptre_assignments ~", as.character(formula)[2]))
  mod_pert <- glm(fmla_pert, data = mod_data, family = binomial())
  logit_prob_pert <- fitted(mod_pert)[cell_idx_subsampled] |> qlogis()
  
  cat("done.\n")
  
  # 4. for each guide we simulate, do the following:
  #   4a. simulate perturbation indicators for each cell
  #   4b. simulate UMI counts for each cell using those
  
  perturbation_indicators <- simulated_umi_counts <- matrix(0, num_guides_to_simulate, num_cells)
  
  cat("Beginning to simulate data:\n")
  
  for(sim_idx in 1:num_guides_to_simulate) {
    # 4a. simulate perturbation indicators for our subsampled cells, using `mod_pert`
    is_perturbed <- rbinom(
      num_cells, size = 1,
      plogis(logit_prob_pert + sim_params$logit_prob_pert_increase)
    ) |>
      as.logical()
    if(sim_idx == 1) {
      cat("The first simulated guide has", sum(is_perturbed), "perturbed cells.\n")
    }
    perturbation_indicators[sim_idx,] <- is_perturbed
    
    
    # 4b. simulate UMI counts, sampling from a 2 component NB mixture with means
    # coming from `mod_umi` 
    log_mean_umi_with_pert <- log_mean_umi + sim_params$pert_effect * is_perturbed +
      sim_params$non_pert_effect * (1 - is_perturbed)
    thetas <- ifelse(is_perturbed, sim_params$theta_pert, sim_params$theta_non_pert)
    simulated_counts <- MASS::rnegbin(
      num_cells,
      mu = exp(log_mean_umi_with_pert),
      theta = thetas
    )
    if(clip_to_real_max) simulated_counts <- pmin(simulated_counts, max_real_umi)
    simulated_umi_counts[sim_idx,] <- simulated_counts
  }
  cat("Done simulating.\n")
  # 5. collect and return results
  dimnames_ <- list(
    paste0(guide_name, "_sim_", 1:num_guides_to_simulate),
    paste0("cell_idx_", cell_idx_subsampled)
  )
  
  sim_results <- list(
    umi_counts = umi_counts,
    sceptre_assignments = sceptre_assignments,
    mod_umi = mod_umi,
    mod_pert = mod_pert,
    cell_idx_subsampled = cell_idx_subsampled,
    sim_params = sim_params,
    clip_to_real_max = clip_to_real_max,
    formula = formula,
    perturbation_indicators = perturbation_indicators |> Matrix::Matrix(sparse=TRUE) |> `dimnames<-`(dimnames_),
    simulated_umi_counts = simulated_umi_counts |> Matrix::Matrix(sparse=TRUE) |> `dimnames<-`(dimnames_)
  )
  return(sim_results)
}

source_data <- "replogle-rd7"

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
cell_covariates <- scep@covariate_data_frame |>
  transmute(
    grna_n_nonzero, grna_n_umis
  ) 


grna_odm <- ondisc::initialize_odm_from_backing_file(paste0(base_fp, "grna.odm"))

assignments_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs", source_data, "sceptre-pipeline/"
)
assn_mat <- readRDS(file.path(assignments_fp,  "grna_assignment_matrix.rds" ))

stopifnot(setequal(rownames(assn_mat), rownames(grna_odm)), dim(assn_mat) == dim(grna_odm))

guide <- "8645_TACC3_P1_ENSG00000013810"

run_name <- "replogle-rd7_simulated_100k_0.015-pert"
sim_results <- make_simulated_umi_counts(
  umi_counts = grna_odm[guide,],
  num_guides_to_simulate = 25,
  cell_covariates = cell_covariates,
  sceptre_assignments = assn_mat[guide,],
  formula = ~ 1 + log(grna_n_nonzero+1) + log(grna_n_umis+1),
  num_cells = 100000,
  sim_params = list(
    logit_prob_pert_increase = 2,
    non_pert_effect = -5,
    pert_effect = 5.25,
    theta_pert = 10,
    theta_non_pert = .1
  ),
  guide_name = guide,
  seed = 165756
)
# run_name <- "replogle-rd7_simulated_50k_0.015-pert"
# sim_results <- make_simulated_umi_counts(
#   umi_counts = grna_odm[guide,],
#   num_guides_to_simulate = 25,
#   cell_covariates = cell_covariates,
#   sceptre_assignments = assn_mat[guide,],
#   formula = ~ 1 + log(grna_n_nonzero+1) + log(grna_n_umis+1),
#   num_cells = 50000,
#   sim_params = list(
#     logit_prob_pert_increase = 1,
#     non_pert_effect = -5,
#     pert_effect = 5.25,
#     theta_pert = 10,
#     theta_non_pert = .1
#   ),
#   guide_name = guide,
#   seed = 165756
# )


idx <- 2
plot_results(
  grna_real = sim_results$umi_counts[sim_results$cell_idx_subsampled],
  grna_sim = sim_results$simulated_umi_counts[idx,],
  is_pert = sim_results$perturbation_indicators[idx,],
  sim_name = guide
)

# saving results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


output_dir <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  paste0("guide_assignment/input_data/", run_name)
)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# 3.0: saving model info and true perturbation status
saveRDS(sim_results$perturbation_indicators, file = file.path(output_dir, "true_perturbation_status.rds"))
saveRDS(sim_results, file = file.path(output_dir, "all_sim_results.rds"))


# 3.1: crispat
crispat_dir <- file.path(output_dir, "crispat")
dir.create(crispat_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing crispat to", crispat_dir))
write_to_h5ad(sim_results$simulated_umi_counts, file.path(crispat_dir, "grna_matrix.h5ad"))

# 3.2: pertpy
pertpy_dir <- file.path(output_dir, "pertpy")
dir.create(pertpy_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing pertpy to", pertpy_dir))
write_to_h5ad(sim_results$simulated_umi_counts, file.path(pertpy_dir, "grna_matrix.h5ad"))

# 3.3: cleanser
cleanser_dir <- file.path(output_dir, "cleanser")
dir.create(cleanser_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing cleanser to", cleanser_dir))
Matrix::writeMM(sim_results$simulated_umi_counts, file.path(cleanser_dir, "grna_matrix.mtx"))

# 3.4: sceptre
sceptre_dir <- file.path(output_dir, "sceptre")
dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
print(paste("  Writing sceptre gRNA matrix and other files to", sceptre_dir))
saveRDS(sim_results$simulated_umi_counts, file.path(sceptre_dir, "grna_matrix.rds"))

# 3.4.1: grna_target_data_frame
guide_target <- scep@grna_target_data_frame |>
  dplyr::filter(grna_id == guide) |> pull(grna_target)
sim_grna_target_df <- data.frame(
  grna_id = rownames(sim_results$simulated_umi_counts),
  grna_target = guide_target
)

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
  dims = c(1L, ncol(sim_results$simulated_umi_counts)),
  dimnames = list(
    "FAKE_GENE_1",  # Single fake gene
    colnames(sim_results$simulated_umi_counts)  # Use same cell names as grna_matrix
  )
)
saveRDS(fake_response_mat, file.path(sceptre_dir, "response_matrix.rds"))


# 3.4.3: covariate_data_frame
cell_covariates[sim_results$cell_idx_subsampled, ] |>
  dplyr::rename(true_grna_n_nonzero = grna_n_nonzero, true_grna_n_umis = grna_n_umis) |>
  write.csv(file.path(sceptre_dir, "covariate_data_frame.csv"), row.names=TRUE)

# 3.4.4 saving formula object for assign_grnas()
saveRDS(sim_results$formula, file.path(sceptre_dir, "assign_grnas_formula.rds"))

cat("Command to transfer to HPC:\n")
cat(paste0("scp -r ", output_dir, " ", .get_config_path("REMOTE_BENCHMARKING_DIR"), "guide_assignment/input_data/", run_name, "\n"))
