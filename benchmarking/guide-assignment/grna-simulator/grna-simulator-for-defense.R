# make simulated datasets


library(Matrix)
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)


rm(list=ls())
source("~/.Rprofile")

# Source utility functions
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "../../association/data-preparation/utils_data_prep.R"))

output_base_dir = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")

sim_data <- function(params, n = NULL, cell_scaling = NULL, seed=1) {
  set.seed(seed)
  if(is.null(n) && is.null(cell_scaling)) {
    stop("at least one of `n` and `cell_scaling` must be provided.")
  }
  if(!is.null(n) && !is.null(cell_scaling)) {
    if(length(n) != length(cell_scaling)) {
      stop("`n` and `cell_scaling` must have the same length.")
    }
  }
  if(is.null(cell_scaling)) {
    cell_scaling <- 1
  }
  if(is.null(n)) {
    n <- length(cell_scaling)
  }

  z <- rbinom(n = n, size = 1, prob = params["pi"])
  is_pert = z == 1
  num_pert = sum(is_pert)
  num_nonpert = n - num_pert
  cell_scaling_pert = cell_scaling[is_pert]
  
  # 1. simulate perturbed cells from NB; optionally conditioned on being positive
  ysim_pos = MASS::rnegbin(
    n = num_pert,
    mu = params["mu_pert"] * cell_scaling_pert,
    theta = params["theta_pert"]
  )

  # 2. simulate non-pert from NB; optionally, zero inflate
  ysim_neg <- MASS::rnegbin(
    n = num_nonpert,
    mu = params["mu_np"] * cell_scaling[z == 0],
    theta = params["theta_np"]
  )

  # collect
  ysim <- numeric(n)
  ysim[z == 0] <- ysim_neg
  ysim[z == 1] <- ysim_pos
  
  list(ysim=ysim, z=z)
}


path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data",
  "replogle-rd7",
  "sceptre-pipeline"
)

make_grna_simulation <- function(dataset_name, output_base_dir, params_list, grna_nonzero_normed, num_guides_per_sim, real_cell_locs, scep, seed) {

  # this should be a i,p,x list but i can get away with this for now
  num_cells = length(grna_nonzero_normed)
  grna_sims = true_pert_matrix = matrix(0, num_guides_per_sim * length(params_list), num_cells)
  k <- 1
  for(param_idx in seq_along(params_list)) {
    for(sim_idx in 1:num_guides_per_sim) {
      curr_sims = sim_data(params = params_list[[param_idx]], cell_scaling = grna_nonzero_normed, seed=k+seed)
      grna_sims[k,] <- curr_sims$ysim
      true_pert_matrix[k,] <- curr_sims$z
      k <- k + 1
    }
  }
  rownames(grna_sims) <- rownames(true_pert_matrix) <- lapply(
    names(params_list),
    function(param_name) paste0(param_name, "_grna_", 1:num_guides_per_sim)) |>
    unlist()
  cell_names = paste0("cell_idx_", real_cell_locs)
  colnames(grna_sims) <- colnames(true_pert_matrix) <- cell_names
  
  grna_sims = as(grna_sims, "RsparseMatrix")
  true_pert_matrix = as(true_pert_matrix, "RsparseMatrix")
  
  output_dir <- file.path(output_base_dir, dataset_name)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  setup_anndata_environment()
  
  
  # saving objects ~~~~~~~~~~~~~~~~
  
  cov_dat <- data.frame(
    real_cell_locs = real_cell_locs,
    grna_nonzero_normed =  grna_nonzero_normed
  )
  sim_params <- list(
    params_list = params_list, num_guides_per_sim = num_guides_per_sim
  )
  
  saveRDS(true_pert_matrix, file.path(output_dir, "true_pert_matrix.rds"))
  saveRDS(cov_dat, file.path(output_dir, "cell_scaling_and_locs.rds"))
  saveRDS(sim_params, file.path(output_dir, "sim_params.rds"))
  
  # 4a,b. crispat and pertpy ~~~~~~~~~~~~~~~~~~~~~~
  
  sce <- SingleCellExperiment(
    assays  = list(counts = grna_sims)
  )
  
  crispat_dir <- file.path(output_dir, "crispat")
  dir.create(crispat_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing crispat to", crispat_dir))
  
  writeH5AD(
    sce,
    file = file.path(crispat_dir, "grna_matrix.h5ad")
  )
  
  pertpy_dir <- file.path(output_dir, "pertpy")
  dir.create(pertpy_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing pertpy to", pertpy_dir))
  
  writeH5AD(
    sce,
    file = file.path(pertpy_dir, "grna_matrix.h5ad")
  )
  
  # 4c: cleanser  ~~~~~~~~~~~~~~~~~~~~~~
  
  cleanser_dir <- file.path(output_dir, "cleanser")
  dir.create(cleanser_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing cleanser to", cleanser_dir))
  Matrix::writeMM(grna_sims, file.path(cleanser_dir, "grna_matrix.mtx"))
  
  # 4d. sceptre  ~~~~~~~~~~~~~~~~~~~~~~
  
  sceptre_dir <- file.path(output_dir, "sceptre")
  dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing sceptre gRNA matrix to", sceptre_dir))
  
  covariates_subset <- scep@covariate_data_frame[real_cell_locs,] |>
    dplyr::select(response_n_umis_full = response_n_umis, response_n_nonzero_full = response_n_nonzero) |>
    as.data.frame()

  dummy_response_mat <- Matrix::sparseMatrix(
    i = integer(0),  # No non-zero entries
    j = integer(0),
    x = numeric(0),
    dims = c(1L, length(cell_names)),
    dimnames = list(
      "FAKE_GENE_1",  # Single fake gene
      cell_names  # Use same cell names as grna_matrix
    )
  )
  
  assign_fmla <- "~ log(response_n_umis_full + 1) + log(response_n_nonzero_full + 1) + log(grna_n_umis + 1) + log(grna_n_nonzero + 1)"

  grna_target_df <- data.frame(
    grna_id = rownames(grna_sims),
    grna_target = paste0("grna_target_", 1:nrow(grna_sims))
  )
  grna_target_df$grna_target[1] = "non-targeting"
  
  write_sceptre_output(
    response_matrix = dummy_response_mat,
    grna_matrix = grna_sims,
    cell_covariates = covariates_subset,
    grna_target_df,
    discovery_pairs = NULL,
    formula_object = assign_fmla,
    output_path = sceptre_dir
  )
  
  cat("scp -r", dataset_name, "hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
}

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

set.seed(111)
num_cells = 50000
locs = sample(nrow(scep@covariate_data_frame), num_cells)
grna_nonzero_normed = (scep@covariate_data_frame$grna_n_nonzero / mean(scep@covariate_data_frame$grna_n_nonzero))[locs]

params_cleanser <- c(
  pi = 0.004,
  mu_pert = 971.437,
  theta_pert = 0.126,
  mu_np = 0.066,
  theta_np = 2.26
)

mu_scaling_powers = (0:3)
mu_pert_seq <- params_cleanser["mu_pert"] / 2^mu_scaling_powers
params_list <- lapply(mu_pert_seq, function(mu) {
  curr_params = params_cleanser
  curr_params["mu_pert"] = mu
  curr_params
}) |>
  setNames(paste0("params_cleanser_mupow=", mu_scaling_powers))


make_grna_simulation(
  dataset_name = "sims_from_cleanser_50k_4_mu",
  output_base_dir = output_base_dir,
  params_list = params_list,
  grna_nonzero_normed = grna_nonzero_normed,
  num_guides_per_sim = 100,
  real_cell_locs = locs,
  scep = scep,
  seed=123
)





