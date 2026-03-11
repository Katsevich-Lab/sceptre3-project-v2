# utils_data_prep.R
# Shared utility functions for data preparation scripts
# Used by make_neg_control_replogle-rd7.R, make_pos_control_replogle-rd7.R, etc.

# Function 1: Convert ODM to sparse matrix ------------------------------
#' Convert OnDisk Matrix subset to sparse matrix
#'
#' @param odm OnDisk Matrix object
#' @param genes Character vector of gene names to extract
#' @param cell_idx Integer vector of cell indices to extract
#' @param set_rownames Logical, whether to set rownames to genes (default TRUE)
#' @return Sparse matrix with genes as rows, cells as columns
odm_to_sparse_matrix <- function(odm, genes, cell_idx, set_rownames = TRUE) {
  ilist <- jlist <- xlist <- vector("list", length(genes))
  for(i in seq_along(genes)) {
    curr_umis <- odm[genes[i], ][cell_idx]
    is_positive <- curr_umis > 0
    num_entries <- sum(is_positive)
    if(num_entries > 0) {
      ilist[[i]] <- rep(i, num_entries)
      jlist[[i]] <- which(is_positive)
      xlist[[i]] <- curr_umis[is_positive]
    }
  }

  result <- Matrix::sparseMatrix(
    i = unlist(ilist),
    j = unlist(jlist),
    x = unlist(xlist),
    dims = c(length(genes), length(cell_idx))
  )

  if(set_rownames) {
    rownames(result) <- genes
  }

  return(result)
}


# Function 2: Enforce single guide per cell (low MOI) -------------------
#' Enforce single guide per cell by random selection
#'
#' @param grna_indicator_matrix Sparse binary matrix (guides × cells)
#' @param random_seed Optional random seed for reproducibility
#' @return Modified grna_indicator_matrix with exactly 1 guide per cell
#' @details Randomly selects one guide per cell if multiple are expressed.
#'          Validates that all cells have exactly 1 guide and all guides have ≥1 cell.
enforce_single_guide_per_cell <- function(grna_indicator_matrix, random_seed = NULL) {
  if(!is.null(random_seed)) {
    set.seed(random_seed)
  }

  # Make a copy to avoid modifying the input
  result <- grna_indicator_matrix
  num_expressed_guides = Matrix::colSums(result)
  if(any(num_expressed_guides == 0)) {
    stop("Some cells have no expressed guides! Should not happen here.")
  }
  
  idx_to_enforce = which(num_expressed_guides > 1)
  for(cell_idx in idx_to_enforce) {
    expressed_guides <- which(result[, cell_idx] == 1)
    guide_to_keep <- sample(expressed_guides, 1)
    guides_to_remove <- setdiff(expressed_guides, guide_to_keep)
    result[guides_to_remove, cell_idx] <- 0
  }
  

  ## old slower way
  # for(cell_idx in seq_len(ncol(result))) {
  #   expressed_guides <- which(result[, cell_idx] == 1)
  # 
  #   if(length(expressed_guides) == 0) {
  #     stop("Cell ", cell_idx, " has no expressed guides! Should not happen.")
  #   }
  # 
  #   if(length(expressed_guides) > 1) {
  #     guide_to_keep <- sample(expressed_guides, 1)
  #     guides_to_remove <- setdiff(expressed_guides, guide_to_keep)
  #     result[guides_to_remove, cell_idx] <- 0
  #   }
  # }

  # Validate output
  stopifnot(all(Matrix::colSums(result) == 1))  # Every cell has exactly 1 guide

  return(result)
}


# Function 3: Write SCEPTRE format output --------------------------------
#' Write data in SCEPTRE format
#'
#' @param response_matrix Sparse matrix (genes × cells)
#' @param grna_matrix Sparse binary matrix (guides × cells)
#' @param cell_covariates Data frame of cell-level covariates
#' @param grna_target_df Data frame with columns: grna_id, grna_target
#' @param discovery_pairs Data frame with columns: grna_target, response_id
#' @param formula_object Formula string for SCEPTRE model
#' @param output_path Directory path to write outputs
write_sceptre_output <- function(
  response_matrix,
  grna_matrix,
  cell_covariates,
  grna_target_df,
  discovery_pairs,
  formula_object,
  output_path
) {
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  # Write files
  saveRDS(response_matrix, file.path(output_path, "response_matrix.rds"))
  saveRDS(grna_matrix, file.path(output_path, "grna_matrix.rds"))
  write.csv(cell_covariates, file.path(output_path, "cell_covariates.csv"),
            row.names = FALSE)
  write.csv(grna_target_df, file.path(output_path, "grna_target_data_frame.csv"),
            row.names = FALSE)
  saveRDS(formula_object, file.path(output_path, "formula_object.rds"))
  saveRDS(discovery_pairs, file.path(output_path, "discovery_pairs.rds"))

  cat("   SCEPTRE format written to:", output_path, "\n")
  cat("   Discovery pairs:", nrow(discovery_pairs), "\n")
}


# Function 4: Write Mixscale format output -------------------------------
#' Write data in Mixscale format
#'
#' @param response_matrix Sparse matrix (genes × cells)
#' @param cell_info Data frame with cell metadata (must have grna_target column)
#' @param output_path Directory path to write outputs
#' @param all_targets Character vector of all target names to test
write_mixscale_output <- function(
  response_matrix,
  cell_info,
  output_path,
  all_targets = NULL
) {
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  # Create cell names
  cell_names <- paste0("cell_idx_", cell_info$cell_idx)

  # Add cell names to response matrix
  response_mat_mixscale <- response_matrix |>
    `colnames<-`(cell_names)

  # Create assignment vector (cell names → targets)
  assign_vec <- cell_info$grna_target |>
    setNames(cell_names)

  # Write files
  saveRDS(response_mat_mixscale, file.path(output_path, "response_matrix.rds"))
  saveRDS(assign_vec, file.path(output_path, "assignments.rds"))
  if(!is.null(all_targets)) {
    saveRDS(all_targets, file.path(output_path, "mixscale_nt_targets.rds"))
  }


  cat("   Mixscale format written to:", output_path, "\n")
  if(!is.null(all_targets)) {
    cat("   Mixscale targets:", length(all_targets), "\n")
  }
}


# Function 5: Setup AnnData conda environment ----------------------------
#' Setup conda environment for AnnData/FR-Perturb
#'
#' @param env_name Name of conda environment (default "r-anndata")
#' @return TRUE if successful, FALSE otherwise
setup_anndata_environment <- function(env_name = "r-anndata") {
  library(reticulate)

  # Check if environment exists
  curr_envs <- conda_list()$name

  if(!env_name %in% curr_envs) {
    cat("Creating conda environment:", env_name, "\n")
    conda_create(env_name, packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"),
               envname = env_name, pip = TRUE)
  }

  use_condaenv(env_name, required = TRUE)
  py_config()

  return(TRUE)
}

# Function 6: Prepare FR-Perturb covariates ---------------------------------
#' Prepare covariates for FR-Perturb format
#'
#' @param cell_covariates Data frame of cell-level covariates
#' @param grna_targets Character vector of gRNA targets for each cell
#' @param covariates_to_log1p Character vector of covariate column names to log1p transform
#' @return Data frame with transformed covariates and perturbation column
prepare_frperturb_covariates <- function(
  cell_covariates,
  grna_targets,
  covariates_to_log1p
) {
  # Validate inputs
  stopifnot(all(covariates_to_log1p %in% names(cell_covariates)))
  stopifnot(nrow(cell_covariates) == length(grna_targets))

  # Apply log1p transformation to specified columns
  cell_covs_frpert <- cell_covariates
  for(col in covariates_to_log1p) {
    new_col_name <- paste0(col, "_log1p")
    cell_covs_frpert[[new_col_name]] <- log1p(cell_covariates[[col]])
  }

  # Remove original untransformed columns and add perturbation
  cell_covs_frpert <- cell_covs_frpert |>
    dplyr::select(-dplyr::all_of(covariates_to_log1p)) |>
    dplyr::mutate(perturbation = grna_targets)

  stopifnot(!any(is.na(cell_covs_frpert)))

  return(cell_covs_frpert)
}


# Function 7: Write FR-Perturb format output -----------------------------
#' Write data in FR-Perturb (AnnData/H5AD) format
#'
#' @param response_matrix Sparse matrix (genes × cells)
#' @param cell_names Character vector of cell names (must match ncol of response_matrix)
#' @param cell_covariates_frpert Data frame of prepared covariates (from prepare_frperturb_covariates)
#' @param output_path Directory path to write outputs
#' @param perturbation_in_covariates if TRUE, expect "perturbation" as a column in `cell_covariates_frpert`
write_frperturb_output <- function(
  response_matrix,
  cell_names,
  cell_covariates_frpert,
  output_path,
  perturbation_in_covariates = TRUE
) {
  if(perturbation_in_covariates && ! "perturbation" %in% colnames(cell_covariates_frpert)) {
    stop('"perturbation" is expected but not found as a colname of `cell_covariates_frpert`.')
  }
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  # Setup conda environment
  setup_anndata_environment()

  # Load required libraries
  library(SingleCellExperiment)
  library(zellkonverter)

  # Validate inputs
  stopifnot(length(cell_names) == ncol(response_matrix))
  stopifnot(nrow(cell_covariates_frpert) == ncol(response_matrix))

  # Add cell names to response matrix
  response_subset_frpert <- response_matrix |> `colnames<-`(cell_names)

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays  = list(counts = response_subset_frpert),
    colData = cell_covariates_frpert 
  )

  # Write H5AD file
  writeH5AD(
    sce,
    file = file.path(output_path, "response_matrix.h5ad"),
    X_name = "counts"
  )

  cat("   FR-Perturb format written to:", output_path, "\n")
}
