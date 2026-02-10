library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

dataset_name <- "gasperini_neg_control"

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/gasperini/sceptre-pipeline"
)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/gasperini/sceptre-pipeline"
)

# Helper function (MODIFIED from positive control to set rownames)
odm_to_sparse_matrix <- function(odm, genes, cell_idx) {
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
  Matrix::sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
                       dims = c(length(genes), length(cell_idx))) |>
    `rownames<-`(genes)
}

make_neg_control_gasperini <- function(
    num_nt_guides,
    num_genes,
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    num_mixscale_nts = 20,
    random_seed = 546
) {
  # Negative control data preparation:
  # 1. Select NT guides and downsample
  # 2. Select genes and downsample
  # 3. Create pseudo-targets for each NT guide
  # 4. Identify cells expressing these NT guides
  # 5. Subset data matrices
  # 6. Write outputs for SCEPTRE and FR-Perturb

  # Section 0: Random sampling (NEW) ------------------------------------
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  if(num_nt_guides == 0) {
    stop("`num_nt_guides` cannot be 0.")
  }
  
  # 0a. Identify all NT guides
  all_nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == "non-targeting"]
  cat("Found", length(all_nt_guides), "NT guides total.\n")

  # 0b. Identify all genes
  all_genes <- rownames(response_odm)
  cat("Found", length(all_genes), "genes total.\n")

  # Input validation
  if (length(all_nt_guides) == 0) {
    stop("No non-targeting guides found in `grna_target_df`")
  }

  # 0c. Downsample NT guides
  if (length(all_nt_guides) < num_nt_guides) {
    warning("Requested ", num_nt_guides, " NT guides but only ",
            length(all_nt_guides), " available. Using all available.")
    nt_guides <- all_nt_guides
    num_nt_guides <- length(all_nt_guides)
  } else {
    nt_guides <- sample(all_nt_guides, num_nt_guides)
  }
  cat("Selected", num_nt_guides, "NT guides.\n")

  # 0d. Downsample genes
  if (length(all_genes) < num_genes) {
    warning("Requested ", num_genes, " genes but only ",
            length(all_genes), " available. Using all available.")
    genes <- all_genes
    num_genes <- length(all_genes)
  } else {
    genes <- sample(all_genes, num_genes)
  }
  cat("Selected", num_genes, "genes.\n")

  # Section 1: Create pseudo-targets -------------------------
  # Create fake grna_target_df with unique pseudo-targets
  nt_grna_target_df <- data.frame(
    grna_id = nt_guides,
    grna_target = paste0("non-targeting", seq_len(num_nt_guides)),
    stringsAsFactors = FALSE
  )
  cat("Created", nrow(nt_grna_target_df), "pseudo-targets.\n")

  # Validation: ensure no duplicate pseudo-targets
  stopifnot(!any(duplicated(nt_grna_target_df$grna_target)))

  # Section 2: Identify cells (SIMILAR TO POSITIVE) ---------------------
  # 2a. Determine which cells express these NT guides
  cells_expressing_nt_guides <- vector("list", length(nt_guides)) |>
    setNames(nt_guides)

  for(guide in nt_guides) {
    cells_expressing_nt_guides[[guide]] <- which(scep_assn_mat[guide,])
  }

  total_cells_found <- unlist(cells_expressing_nt_guides) |> unique() |> length()
  cat("There are", total_cells_found, "cells expressing these NT guides.\n")

  # Validation: ensure cells found
  if (total_cells_found == 0) {
    stop("No cells found expressing the selected NT guides")
  }

  # 2b. Create cell_info data.frame (using pseudo-targets, not real targets)
  cell_info <- lapply(
    names(cells_expressing_nt_guides),
    function(guide) {
      data.frame(
        cell_id = cells_expressing_nt_guides[[guide]],
        grna_id = guide
      )
    }
  ) |>
    do.call(what = rbind) |>
    left_join(nt_grna_target_df |> dplyr::select(grna_id, grna_target), by = "grna_id") |>
    mutate(cell_name = rownames(cell_covariates)[cell_id])

  # Section 3: Prepare matrices (REUSE) ---------------------------------
  # 3a. Get unique cell indices
  all_cell_idx <- unique(cell_info$cell_id) |> sort()

  # 3b. Subset response matrix for selected genes (rownames now set by helper)
  response_subset <- odm_to_sparse_matrix(
    odm = response_odm,
    genes = genes,
    cell_idx = all_cell_idx
  )

  cat("Response matrix subset made with",
      nrow(response_subset), "genes and",
      ncol(response_subset), "cells.\n")

  # Validation: dimensions match
  stopifnot(nrow(response_subset) == num_genes & rownames(response_subset) == genes)
  stopifnot(ncol(response_subset) == length(all_cell_idx))

  # 3c. Subset cell covariates (rename to *_full for formula compatibility)
  cell_covariates_subset <- cell_covariates[all_cell_idx, ] |>
    transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full = response_n_umis,
      grna_n_nonzero_full = grna_n_nonzero,
      grna_n_umis_full = grna_n_umis,
      response_p_mito_full = response_p_mito,
      prep_batch
    )

  # Section 4: Create output directory (MODIFIED) -----------------------
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/neg-control/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)

  # Write shared files
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")

  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")

  ## Section 5a: SCEPTRE format (same as for pos-control except with pseudo-targets) -------------
  write_sceptre_fp <- file.path(write_fp, "sceptre")
  dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)

  # Use the pseudo-target grna_target_df (not the original)
  grna_target_df_kept <- nt_grna_target_df

  # Create binary assignment indicator matrix (SIMPLIFIED)
  # Convert ngRMatrix (logical sparse) to integer sparse by adding 0
  grna_indicator_matrix <- scep_assn_mat[nt_guides, all_cell_idx] + 0

  # Compute sceptre-specific covariates from the subsetted gRNA UMI data
  grna_matrix_subset <- odm_to_sparse_matrix(
    grna_odm,
    nt_guides,
    all_cell_idx
  )

  grna_n_nonzero_subset <- Matrix::colSums(grna_matrix_subset != 0)
  grna_n_umis_subset <- Matrix::colSums(grna_matrix_subset)

  cell_covariates_sceptre <- cbind(
    cell_covariates_subset,
    data.frame(
      grna_n_nonzero_subset = grna_n_nonzero_subset,
      grna_n_umis_subset = grna_n_umis_subset
    )
  )

  # Save formula as character string for sceptre
  formula_string <- "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1) + prep_batch"

  # Create discovery pairs (Cartesian product of all pseudo-targets × all genes)
  discovery_pairs <- expand.grid(
    grna_target = nt_grna_target_df$grna_target,
    response_id = rownames(response_subset),
    stringsAsFactors = FALSE
  )
  cat("   Discovery pairs:", nrow(discovery_pairs), "\n")

  saveRDS(response_subset, file.path(write_sceptre_fp, "response_matrix.rds"))
  saveRDS(grna_indicator_matrix, file.path(write_sceptre_fp, "grna_matrix.rds"))
  write.csv(cell_covariates_sceptre, file.path(write_sceptre_fp, "cell_covariates.csv"), row.names = FALSE)
  write.csv(grna_target_df_kept, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names = FALSE)
  saveRDS(formula_string, file.path(write_sceptre_fp, "formula_object.rds"))
  saveRDS(discovery_pairs, file.path(write_sceptre_fp, "discovery_pairs.rds"))
  cat("   sceptre written.\n")

  ## Section 5b: FR-Perturb format (COMMENTED OUT FOR NOW) -------------
  # write_frperturb_fp <- file.path(write_fp, "frperturb")
  # dir.create(write_frperturb_fp, showWarnings = FALSE, recursive = TRUE)
  #
  # library(reticulate)
  # library(SingleCellExperiment)
  # library(zellkonverter)
  #
  # # Setup conda environment (same as positive control)
  # env_name <- "r-anndata"
  # curr_envs <- conda_list()$name
  # if(!env_name %in% curr_envs) {
  #   conda_create("r-anndata", packages = c("python=3.12"))
  #   py_install(c("numpy", "scipy", "h5py", "anndata"),
  #              envname = "r-anndata", pip = TRUE)
  # }
  # use_condaenv("r-anndata", required = TRUE)
  # py_config()
  #
  # # Prepare data for FR-Perturb
  # cell_names <- rownames(cell_covariates_subset)
  # response_subset_frpert <- response_subset |> `colnames<-`(cell_names)
  #
  # # Compute response covariates
  # response_n_nonzero_subset <- Matrix::colSums(response_subset_frpert > 0)
  # response_n_umis_subset <- Matrix::colSums(response_subset_frpert)
  #
  # # Build covariate dataframe (include both *_full and *_subset, plus log-transformed)
  # cell_covs_frpert <- dplyr::select(
  #   cell_covariates_sceptre,
  #   # Keep *_full versions
  #   response_n_nonzero_full,
  #   response_n_umis_full,
  #   grna_n_nonzero_full,
  #   grna_n_umis_full,
  #   response_p_mito_full,
  #   prep_batch,
  #   # Also include *_subset versions
  #   grna_n_nonzero_subset,
  #   grna_n_umis_subset
  # ) |>
  #   mutate(
  #     response_n_nonzero_subset = response_n_nonzero_subset,
  #     response_n_umis_subset = response_n_umis_subset
  #   )
  #
  # # Create perturbation column (high MOI - use pseudo-targets with ":" delimiter)
  # # Verify no pseudo-targets contain ":"
  # stopifnot(!any(grepl(":", nt_grna_target_df$grna_target)))
  #
  # perturb_df <- cell_info |>
  #   group_by(cell_name) |>
  #   summarise(perturbation = paste0(grna_target, collapse = ":"))
  #
  # cell_covs_frpert <- left_join(
  #   cell_covs_frpert %>% mutate(cell_name = rownames(.)),
  #   perturb_df,
  #   by = "cell_name"
  # ) |>
  #   mutate(
  #     # Log-transform covariates (FR-Perturb doesn't take logs)
  #     log_grna_n_nonzero_subset = log(grna_n_nonzero_subset),
  #     log_grna_n_umis_subset = log(grna_n_umis_subset),
  #     log_response_n_nonzero = log(response_n_nonzero_subset),
  #     log_response_n_umis = log(response_n_umis_subset)
  #   )
  #
  # sce <- SingleCellExperiment(
  #   assays  = list(counts = response_subset_frpert),
  #   colData = cell_covs_frpert
  # )
  #
  # writeH5AD(
  #   sce,
  #   file = file.path(write_frperturb_fp, "response_matrix.h5ad"),
  #   X_name = "counts"
  # )
  # cat("   FR-perturb written.\n")
}

# Run the function
make_neg_control_gasperini(
  dataset_name,
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame
)

# ============================================================================
# UNIT TESTS (commented out for normal runs)
# ============================================================================

# # Test 1: Basic functionality - verify all output files created
# test_basic_functionality <- function() {
#   cat("\n=== Test 1: Basic functionality ===\n")
#
#   output_dir <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/neg-control/input_data/gasperini_neg_control"
#   )
#
#   # Check shared files
#   stopifnot(file.exists(file.path(output_dir, "cell_info.csv")))
#   stopifnot(file.exists(file.path(output_dir, "cell_covariates.csv")))
#
#   # Check SCEPTRE files
#   sceptre_dir <- file.path(output_dir, "sceptre")
#   stopifnot(file.exists(file.path(sceptre_dir, "response_matrix.rds")))
#   stopifnot(file.exists(file.path(sceptre_dir, "grna_matrix.rds")))
#   stopifnot(file.exists(file.path(sceptre_dir, "cell_covariates.csv")))
#   stopifnot(file.exists(file.path(sceptre_dir, "grna_target_data_frame.csv")))
#   stopifnot(file.exists(file.path(sceptre_dir, "formula_object.rds")))
#
#   # Check FR-Perturb files (COMMENTED OUT)
#   # frperturb_dir <- file.path(output_dir, "frperturb")
#   # stopifnot(file.exists(file.path(frperturb_dir, "response_matrix.h5ad")))
#
#   cat("✓ All expected SCEPTRE output files exist\n")
# }
#
# # Test 2: Data dimensions and structure
# test_data_dimensions <- function() {
#   cat("\n=== Test 2: Data dimensions and structure ===\n")
#
#   output_dir <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/neg-control/input_data/gasperini_neg_control"
#   )
#   sceptre_dir <- file.path(output_dir, "sceptre")
#
#   # Load data
#   response_matrix <- readRDS(file.path(sceptre_dir, "response_matrix.rds"))
#   grna_matrix <- readRDS(file.path(sceptre_dir, "grna_matrix.rds"))
#   grna_target_df <- read.csv(file.path(sceptre_dir, "grna_target_data_frame.csv"))
#   cell_info <- read.csv(file.path(output_dir, "cell_info.csv"))
#
#   # Check dimensions
#   expected_genes <- 1000
#   expected_guides <- 100
#
#   stopifnot(nrow(response_matrix) == expected_genes)
#   stopifnot(nrow(grna_matrix) == expected_guides)
#   stopifnot(nrow(grna_target_df) == expected_guides)
#
#   cat(sprintf("✓ Response matrix: %d genes × %d cells\n",
#               nrow(response_matrix), ncol(response_matrix)))
#   cat(sprintf("✓ gRNA matrix: %d guides × %d cells\n",
#               nrow(grna_matrix), ncol(grna_matrix)))
#   cat(sprintf("✓ Number of cells match: %d\n", ncol(response_matrix)))
#
#   # Check that rownames are set
#   stopifnot(!is.null(rownames(response_matrix)))
#   stopifnot(!is.null(rownames(grna_matrix)))
#   cat("✓ Rownames are set for matrices\n")
# }
#
# # Test 3: Pseudo-target naming
# test_pseudo_targets <- function() {
#   cat("\n=== Test 3: Pseudo-target naming ===\n")
#
#   output_dir <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/neg-control/input_data/gasperini_neg_control"
#   )
#
#   grna_target_df <- read.csv(file.path(output_dir, "sceptre/grna_target_data_frame.csv"))
#
#   # Check pseudo-target format
#   expected_pattern <- "^non-targeting[0-9]+$"
#   all_match <- all(grepl(expected_pattern, grna_target_df$grna_target))
#   stopifnot(all_match)
#   cat("✓ All pseudo-targets match pattern 'non-targeting[N]'\n")
#
#   # Check uniqueness
#   stopifnot(!any(duplicated(grna_target_df$grna_target)))
#   cat("✓ All pseudo-targets are unique\n")
#
#   # Check no colons (for FR-Perturb safety)
#   stopifnot(!any(grepl(":", grna_target_df$grna_target)))
#   cat("✓ No colons in pseudo-target names\n")
#
#   # Check sequential numbering
#   numbers <- as.integer(sub("non-targeting", "", grna_target_df$grna_target))
#   stopifnot(all(numbers == 1:100))
#   cat("✓ Pseudo-targets numbered sequentially 1-100\n")
# }
#
# # Test 4: Formula object
# test_formula_object <- function() {
#   cat("\n=== Test 4: Formula object ===\n")
#
#   output_dir <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/neg-control/input_data/gasperini_neg_control"
#   )
#
#   formula_obj <- readRDS(file.path(output_dir, "sceptre/formula_object.rds"))
#
#   # Check it's a formula
#   stopifnot(inherits(formula_obj, "formula"))
#   cat("✓ Formula object is of class 'formula'\n")
#
#   # Check formula structure
#   expected_formula <- ~ 1 + log(response_n_nonzero_full + 1) +
#     log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) +
#     log(grna_n_umis_full + 1) + prep_batch
#
#   stopifnot(identical(formula_obj, expected_formula))
#   cat("✓ Formula matches expected structure\n")
#
#   # Print formula
#   cat("Formula:", deparse(formula_obj), "\n")
# }
#
# # Test 5: Covariate columns
# test_covariates <- function() {
#   cat("\n=== Test 5: Covariate columns ===\n")
#
#   output_dir <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/neg-control/input_data/gasperini_neg_control"
#   )
#
#   # Check SCEPTRE covariates
#   cell_covs_sceptre <- read.csv(file.path(output_dir, "sceptre/cell_covariates.csv"))
#
#   expected_cols_sceptre <- c(
#     "response_n_nonzero_full", "response_n_umis_full",
#     "grna_n_nonzero_full", "grna_n_umis_full",
#     "response_p_mito_full", "prep_batch",
#     "grna_n_nonzero_subset", "grna_n_umis_subset"
#   )
#
#   stopifnot(all(expected_cols_sceptre %in% colnames(cell_covs_sceptre)))
#   cat("✓ SCEPTRE covariates have all expected columns\n")
#   cat("  Columns:", paste(colnames(cell_covs_sceptre), collapse=", "), "\n")
#
#   # Check FR-Perturb covariates (COMMENTED OUT)
#   # library(reticulate)
#   # use_condaenv("r-anndata", required = TRUE)
#   # anndata <- import("anndata")
#   # adata <- anndata$read_h5ad(file.path(output_dir, "frperturb/response_matrix.h5ad"))
#   #
#   # expected_cols_frperturb <- c(
#   #   "response_n_nonzero_full", "response_n_umis_full",
#   #   "grna_n_nonzero_full", "grna_n_umis_full",
#   #   "response_p_mito_full", "prep_batch",
#   #   "grna_n_nonzero_subset", "grna_n_umis_subset",
#   #   "response_n_nonzero_subset", "response_n_umis_subset",
#   #   "perturbation",
#   #   "log_grna_n_nonzero_subset", "log_grna_n_umis_subset",
#   #   "log_response_n_nonzero", "log_response_n_umis"
#   # )
#   #
#   # obs_cols <- names(adata$obs)
#   # missing_cols <- setdiff(expected_cols_frperturb, obs_cols)
#   # if (length(missing_cols) > 0) {
#   #   cat("Missing columns:", paste(missing_cols, collapse=", "), "\n")
#   #   stop("FR-Perturb covariates missing expected columns")
#   # }
#   # cat("✓ FR-Perturb covariates have all expected columns\n")
#   # cat("  Columns:", paste(obs_cols, collapse=", "), "\n")
# }
#
# # Test 6: Reproducibility
# test_reproducibility <- function() {
#   cat("\n=== Test 6: Reproducibility ===\n")
#
#   temp_dir <- tempdir()
#   temp_dataset1 <- file.path(temp_dir, "test_repro_1")
#   temp_dataset2 <- file.path(temp_dir, "test_repro_2")
#
#   # Run twice with same seed
#   response_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
#   grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
#   cell_covariates <- scep@covariate_data_frame
#   scep_assn_mat <- read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
#   grna_target_df_input <- scep@grna_target_data_frame
#
#   # First run
#   make_neg_control_gasperini(
#     temp_dataset1,
#     response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#     num_nt_guides = 50, num_genes = 500, random_seed = 123
#   )
#
#   # Second run
#   make_neg_control_gasperini(
#     temp_dataset2,
#     response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#     num_nt_guides = 50, num_genes = 500, random_seed = 123
#   )
#
#   # Compare outputs
#   base_path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
#                          "association/neg-control/input_data")
#
#   grna_df1 <- read.csv(file.path(base_path, temp_dataset1, "sceptre/grna_target_data_frame.csv"))
#   grna_df2 <- read.csv(file.path(base_path, temp_dataset2, "sceptre/grna_target_data_frame.csv"))
#
#   response_mat1 <- readRDS(file.path(base_path, temp_dataset1, "sceptre/response_matrix.rds"))
#   response_mat2 <- readRDS(file.path(base_path, temp_dataset2, "sceptre/response_matrix.rds"))
#
#   stopifnot(identical(grna_df1, grna_df2))
#   stopifnot(identical(rownames(response_mat1), rownames(response_mat2)))
#   cat("✓ Two runs with same seed produce identical outputs\n")
#
#   # Clean up
#   unlink(file.path(base_path, temp_dataset1), recursive = TRUE)
#   unlink(file.path(base_path, temp_dataset2), recursive = TRUE)
# }
#
# # Test 7: Different seed produces different selection
# test_different_seed <- function() {
#   cat("\n=== Test 7: Different seed produces different selection ===\n")
#
#   temp_dir <- tempdir()
#   temp_dataset1 <- file.path(temp_dir, "test_seed_1")
#   temp_dataset2 <- file.path(temp_dir, "test_seed_2")
#
#   response_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
#   grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
#   cell_covariates <- scep@covariate_data_frame
#   scep_assn_mat <- read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
#   grna_target_df_input <- scep@grna_target_data_frame
#
#   # Run with different seeds
#   make_neg_control_gasperini(
#     temp_dataset1,
#     response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#     num_nt_guides = 50, num_genes = 500, random_seed = 123
#   )
#
#   make_neg_control_gasperini(
#     temp_dataset2,
#     response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#     num_nt_guides = 50, num_genes = 500, random_seed = 456
#   )
#
#   base_path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
#                          "association/neg-control/input_data")
#
#   # Compare outputs - should be different
#   grna_df1 <- read.csv(file.path(base_path, temp_dataset1, "sceptre/grna_target_data_frame.csv"))
#   grna_df2 <- read.csv(file.path(base_path, temp_dataset2, "sceptre/grna_target_data_frame.csv"))
#
#   response_mat1 <- readRDS(file.path(base_path, temp_dataset1, "sceptre/response_matrix.rds"))
#   response_mat2 <- readRDS(file.path(base_path, temp_dataset2, "sceptre/response_matrix.rds"))
#
#   # Check that guides are different
#   guides_same <- identical(grna_df1$grna_id, grna_df2$grna_id)
#   stopifnot(!guides_same)
#   cat("✓ Different seeds produce different guide selection\n")
#
#   # Check that genes are different
#   genes_same <- identical(rownames(response_mat1), rownames(response_mat2))
#   stopifnot(!genes_same)
#   cat("✓ Different seeds produce different gene selection\n")
#
#   # Clean up
#   unlink(file.path(base_path, temp_dataset1), recursive = TRUE)
#   unlink(file.path(base_path, temp_dataset2), recursive = TRUE)
# }
#
# # Test 8: Edge case - request more guides than available
# test_edge_case_guides <- function() {
#   cat("\n=== Test 8: Edge case - more guides requested than available ===\n")
#
#   temp_dataset <- "test_edge_guides"
#
#   response_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
#   grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
#   cell_covariates <- scep@covariate_data_frame
#   scep_assn_mat <- read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
#   grna_target_df_input <- scep@grna_target_data_frame
#
#   # Request more guides than available (101 available, request 200)
#   output <- capture.output({
#     make_neg_control_gasperini(
#       temp_dataset,
#       response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#       num_nt_guides = 200, num_genes = 500, random_seed = 123
#     )
#   }, type = "message")
#
#   # Check that warning was issued
#   has_warning <- any(grepl("available", output))
#   stopifnot(has_warning)
#   cat("✓ Warning issued when requesting more guides than available\n")
#
#   # Check that all available guides were used
#   base_path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
#                          "association/neg-control/input_data")
#   grna_df <- read.csv(file.path(base_path, temp_dataset, "sceptre/grna_target_data_frame.csv"))
#   stopifnot(nrow(grna_df) == 101)  # All available guides used
#   cat(sprintf("✓ Used all %d available guides\n", nrow(grna_df)))
#
#   # Clean up
#   unlink(file.path(base_path, temp_dataset), recursive = TRUE)
# }
#
# # Test 9: Parameter variation
# test_parameter_variation <- function() {
#   cat("\n=== Test 9: Parameter variation ===\n")
#
#   temp_dataset <- "test_param_variation"
#
#   response_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
#   grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
#   cell_covariates <- scep@covariate_data_frame
#   scep_assn_mat <- read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
#   grna_target_df_input <- scep@grna_target_data_frame
#
#   # Run with custom parameters
#   custom_guides <- 30
#   custom_genes <- 200
#
#   make_neg_control_gasperini(
#     temp_dataset,
#     response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df_input,
#     num_nt_guides = custom_guides, num_genes = custom_genes, random_seed = 999
#   )
#
#   base_path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
#                          "association/neg-control/input_data")
#
#   # Check dimensions match custom parameters
#   response_mat <- readRDS(file.path(base_path, temp_dataset, "sceptre/response_matrix.rds"))
#   grna_mat <- readRDS(file.path(base_path, temp_dataset, "sceptre/grna_matrix.rds"))
#
#   stopifnot(nrow(response_mat) == custom_genes)
#   stopifnot(nrow(grna_mat) == custom_guides)
#
#   cat(sprintf("✓ Custom parameters applied: %d guides, %d genes\n",
#               custom_guides, custom_genes))
#
#   # Clean up
#   unlink(file.path(base_path, temp_dataset), recursive = TRUE)
# }
#
# # Run all tests
# run_all_tests <- function() {
#   cat("\n")
#   cat("================================================================================\n")
#   cat("RUNNING UNIT TESTS FOR make_neg_control_gasperini()\n")
#   cat("================================================================================\n")
#
#   tests <- list(
#     test_basic_functionality,
#     test_data_dimensions,
#     test_pseudo_targets,
#     test_formula_object,
#     test_covariates,
#     test_reproducibility,
#     test_different_seed,
#     test_edge_case_guides,
#     test_parameter_variation
#   )
#
#   n_passed <- 0
#   n_failed <- 0
#
#   for (i in seq_along(tests)) {
#     tryCatch({
#       tests[[i]]()
#       n_passed <- n_passed + 1
#     }, error = function(e) {
#       cat("✗ Test failed:", conditionMessage(e), "\n")
#       n_failed <- n_failed + 1
#     })
#   }
#
#   cat("\n")
#   cat("================================================================================\n")
#   cat(sprintf("TEST SUMMARY: %d passed, %d failed\n", n_passed, n_failed))
#   cat("================================================================================\n")
#
#   if (n_failed > 0) {
#     stop("Some tests failed")
#   }
# }
#
# # Uncomment to run tests:
# # run_all_tests()
