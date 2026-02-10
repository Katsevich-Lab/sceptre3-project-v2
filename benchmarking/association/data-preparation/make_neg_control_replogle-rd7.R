library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/replogle-rd7/sceptre-pipeline"
)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/replogle-rd7/sceptre-pipeline"
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

make_neg_control_replogle <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    num_genes = 5000,
    num_mixscale_nts = 20,
    random_seed = 243535
) {
  # Negative control data preparation for Replogle (LOW MOI):
  # 1. Use ALL NT guides (we don't have many NT cells to spare)
  # 2. Select genes and downsample
  # 3. Filter to cells expressing NO targeting guides and at least one NT guide
  # 4. Randomly assign single guide per cell (enforce low MOI)
  # 5. Filter cells/genes with all-0 expression
  # 6. Create pseudo-targets for each NT guide
  # 7. Write outputs for SCEPTRE and FR-Perturb

  # Section 0: Random sampling (MODIFIED FOR LOW MOI) -------------------
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  # 0a. Identify all NT guides
  all_nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == "non-targeting"]
  cat("Found", length(all_nt_guides), "NT guides total.\n")

  # Input validation
  if (length(all_nt_guides) == 0) {
    stop("No non-targeting guides found in grna_target_df")
  }

  # 0c. Identify all genes
  all_genes <- rownames(response_odm)
  cat("Found", length(all_genes), "genes total.\n")

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

  # Section 1: Identify NT-only cells (LOW MOI SPECIFIC) ----------------
  # Get cells expressing NO targeting guides AND at least one NT guide
  targeting_guides <- setdiff(rownames(scep_assn_mat), all_nt_guides)

  all_cell_idx <- intersect(
    which(Matrix::colSums(scep_assn_mat[targeting_guides, ]) == 0),
    which(Matrix::colSums(scep_assn_mat[all_nt_guides, ]) > 0)
  )
  cat("Found", length(all_cell_idx), "cells expressing no targeting guides and at least one NT guide.\n")

  if (length(all_cell_idx) == 0) {
    stop("No cells found expressing only NT guides")
  }

  # Section 2: Subset response matrix and filter cells -------------------
  # 2a. Initial subset
  response_subset <- odm_to_sparse_matrix(response_odm, genes, all_cell_idx)

  # 2b. Remove cells with no UMIs for selected genes
  cells_with_no_umis <- Matrix::colSums(response_subset) == 0
  cat("After initial subsetting,", sum(cells_with_no_umis),
      "cells with no UMIs are removed.\n")

  all_cell_idx <- all_cell_idx[!cells_with_no_umis]
  response_subset <- response_subset[, !cells_with_no_umis]

  stopifnot(Matrix::colSums(response_subset) > 0)
  stopifnot(ncol(response_subset) == length(all_cell_idx))

  # 2c. Remove genes with all-0 expression
  genes_with_no_umis <- Matrix::rowSums(response_subset) == 0
  response_subset <- response_subset[!genes_with_no_umis, ]
  cat(sum(genes_with_no_umis), "genes removed due to all-0 expression.\n")

  genes <- rownames(response_subset)
  num_genes <- length(genes)
  cat("Final dataset:", num_genes, "genes.\n")

  # Section 3: Enforce single guide per cell (LOW MOI) ------------------
  # Create binary assignment matrix and randomly pick one guide per cell
  grna_indicator_matrix <- scep_assn_mat[all_nt_guides, all_cell_idx] + 0L

  for(cell_idx in seq_len(ncol(grna_indicator_matrix))) {
    expressed_guides <- which(grna_indicator_matrix[, cell_idx] == 1)

    if(length(expressed_guides) == 0) {
      stop("Cell ", cell_idx, " has no expressed guides! Should not happen.")
    }

    if(length(expressed_guides) > 1) {
      guide_to_keep <- sample(expressed_guides, 1)
      guides_to_remove <- setdiff(expressed_guides, guide_to_keep)
      grna_indicator_matrix[guides_to_remove, cell_idx] <- 0
    }
  }

  # Validate: every cell has exactly one guide, every guide has at least one cell
  stopifnot(Matrix::colSums(grna_indicator_matrix) == 1)
  stopifnot(Matrix::rowSums(grna_indicator_matrix) >= 1)
  cat("Enforced single guide per cell (low MOI).\n")

  # Section 4: Create pseudo-targets -------------------------------------
  nt_grna_target_df <- data.frame(
    grna_id = all_nt_guides,
    grna_target = paste0("non-targeting", seq_along(all_nt_guides)),
    stringsAsFactors = FALSE
  )
  cat("Created", nrow(nt_grna_target_df), "pseudo-targets.\n")

  stopifnot(!any(duplicated(nt_grna_target_df$grna_target)))

  # Section 5: Create cell_info ------------------------------------------
  cell_info <- data.frame(
    cell_idx = all_cell_idx,
    expressed_guide = rownames(grna_indicator_matrix)[
      apply(grna_indicator_matrix, 2, function(col) which(col == 1))
    ],
    stringsAsFactors = FALSE
  ) |>
    left_join(nt_grna_target_df, by = c("expressed_guide" = "grna_id")) |>
    mutate(cell_name = rownames(cell_covariates)[cell_idx])

  # Section 6: Prepare covariates (NO BATCH for Replogle) ---------------
  cell_covariates_subset <- cell_covariates[all_cell_idx, ] |>
    transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full = response_n_umis,
      grna_n_nonzero_full = grna_n_nonzero,
      grna_n_umis_full = grna_n_umis
      # No prep_batch or response_p_mito for Replogle
    )
  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))

  # Section 7: Create output directory -----------------------------------
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/neg-control/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)

  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")

  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")

  ## Section 8a: SCEPTRE format -----------------------------------------
  write_sceptre_fp <- file.path(write_fp, "sceptre")
  dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)

  grna_target_df_kept <- nt_grna_target_df

  # NOTE: i am not doing subset covariates. Just full ones
  # Compute sceptre-specific covariates
  # grna_matrix_subset <- odm_to_sparse_matrix(
  #   grna_odm,
  #   all_nt_guides,
  #   all_cell_idx
  # )

  # grna_n_nonzero_subset <- Matrix::colSums(grna_matrix_subset != 0)
  # grna_n_umis_subset <- Matrix::colSums(grna_matrix_subset)
  cell_covariates_sceptre = cell_covariates_subset
  # cell_covariates_sceptre <- cbind(
  #   cell_covariates_subset,
  #   data.frame(
  #     grna_n_nonzero_subset = grna_n_nonzero_subset,
  #     grna_n_umis_subset = grna_n_umis_subset
  #   )
  # )

  # Save formula object (NO BATCH for Replogle)
  formula_object <- ~ 1 + log(response_n_nonzero_full + 1) +
    log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) +
    log(grna_n_umis_full + 1)

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
  saveRDS(formula_object, file.path(write_sceptre_fp, "formula_object.rds"))
  saveRDS(discovery_pairs, file.path(write_sceptre_fp, "discovery_pairs.rds"))
  cat("   sceptre written.\n")

  ## Section 8b: FR-Perturb format --------------------------------------
  write_frperturb_fp <- file.path(write_fp, "frperturb")
  dir.create(write_frperturb_fp, showWarnings = FALSE, recursive = TRUE)

  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)

  # Setup conda environment
  env_name <- "r-anndata"
  curr_envs <- conda_list()$name
  if(!env_name %in% curr_envs) {
    conda_create("r-anndata", packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"),
               envname = "r-anndata", pip = TRUE)
  }
  use_condaenv("r-anndata", required = TRUE)
  py_config()

  # Prepare data for FR-Perturb (low MOI - simple perturbation column)
  cell_names <- paste0("cell_idx_", all_cell_idx)
  response_subset_frpert <- response_subset |> `colnames<-`(cell_names)

  # Build covariate dataframe with log-transformed covariates
  cell_covs_frpert <- cell_covariates_subset |>
    mutate_all(log1p) %>%
    `names<-`(paste0(names(.), "_log1p")) |>
    mutate(perturbation = cell_info$grna_target)

  stopifnot(!any(is.na(cell_covs_frpert)))

  sce <- SingleCellExperiment(
    assays  = list(counts = response_subset_frpert),
    colData = cell_covs_frpert
  )

  writeH5AD(
    sce,
    file = file.path(write_frperturb_fp, "response_matrix.h5ad"),
    X_name = "counts"
  )
  cat("   FR-perturb written.\n")

  ## Section 8c: Mixscale format ----------------------------------------
  write_mixscale_fp <- file.path(write_fp, "mixscale")
  dir.create(write_mixscale_fp, showWarnings = FALSE, recursive = TRUE)

  # Prepare response matrix with cell names as column names
  response_mat_mixscale <- response_subset |>
    `colnames<-`(cell_names)

  # Create assignments vector: named vector mapping cell names to pseudo-targets
  assign_vec <- cell_info$grna_target |>
    setNames(cell_names)

  # Sample pseudo-targets for Mixscale control group
  if(num_mixscale_nts > nrow(nt_grna_target_df)) {
    warning("num_mixscale_nts (", num_mixscale_nts, ") exceeds available pseudo-targets (",
            nrow(nt_grna_target_df), "). Using all pseudo-targets.")
    mixscale_nt_targets <- nt_grna_target_df$grna_target
  } else {
    mixscale_nt_targets <- sample(nt_grna_target_df$grna_target, num_mixscale_nts)
  }
  cat("   Mixscale pseudo-targets to loop over:", length(mixscale_nt_targets), "\n")

  saveRDS(response_mat_mixscale, file.path(write_mixscale_fp, "response_matrix.rds"))
  saveRDS(assign_vec, file.path(write_mixscale_fp, "assignments.rds"))
  saveRDS(mixscale_nt_targets, file.path(write_mixscale_fp, "mixscale_nt_targets.rds"))
  cat("   Mixscale written.\n")
}

# response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
# grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
# cell_covariates = scep@covariate_data_frame
# scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
# grna_target_df = scep@grna_target_data_frame

# Run the function
num_genes = 100
make_neg_control_replogle(
  dataset_name = paste0("replogle-rd7_neg_control_", num_genes, "genes"),
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame,
  num_genes = num_genes,
  num_mixscale_nts = 6
)


num_genes = 5000
make_neg_control_replogle(
  dataset_name = paste0("replogle-rd7_neg_control_", num_genes, "genes"),
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame,
  num_genes = num_genes,
  num_mixscale_nts = 6
)
