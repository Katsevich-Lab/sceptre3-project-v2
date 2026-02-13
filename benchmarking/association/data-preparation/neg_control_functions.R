# neg_control_functions.R
# Negative control specific functions for data preparation
# Used by make_neg_control_replogle-rd7.R and tested by test-neg_control_functions.R

#' Select genes randomly from response ODM
#'
#' @param response_odm OnDisk Matrix of gene expression
#' @param num_genes Number of genes to select
#' @param random_seed Optional random seed
#' @return Character vector of selected gene names
select_genes_random <- function(response_odm, num_genes, random_seed = NULL) {
  if(!is.null(random_seed)) {
    set.seed(random_seed)
  }

  all_genes <- rownames(response_odm)
  cat("Found", length(all_genes), "genes total.\n")

  if (length(all_genes) < num_genes) {
    warning("Requested ", num_genes, " genes but only ",
            length(all_genes), " available. Using all available.")
    return(all_genes)
  } else {
    genes <- sample(all_genes, num_genes)
    cat("Selected", num_genes, "genes.\n")
    return(genes)
  }
}


#' Identify cells expressing only NT guides (no targeting guides)
#'
#' @param scep_assn_mat Sparse binary matrix of guide assignments (guides × cells)
#' @param grna_target_df Data frame with grna_id and grna_target columns
#' @param nt_target_name Target name for non-targeting guides (default "non-targeting")
#' @param targets_to_ignore Character vector of additional target names to ignore (e.g., "nt_off_target", "unknown")
#' @return Integer vector of cell indices
identify_nt_only_cells <- function(scep_assn_mat, grna_target_df,
                                   nt_target_name = "non-targeting",
                                   targets_to_ignore = character(0)) {
  # Identify NT guides
  nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == nt_target_name]
  cat("Found", length(nt_guides), "NT guides total.\n")

  if (length(nt_guides) == 0) {
    stop("No non-targeting guides found with target name '", nt_target_name, "'")
  }

  # Identify targeting guides (all guides NOT in nt_target_name or targets_to_ignore)
  non_targeting_targets <- c(nt_target_name, targets_to_ignore)
  targeting_guides <- grna_target_df$grna_id[!grna_target_df$grna_target %in% non_targeting_targets]
  stopifnot(targeting_guides %in% rownames(scep_assn_mat)) # these should all be present

  cat("Found", length(targeting_guides), "targeting guides (excluding",
      paste(non_targeting_targets, collapse=", "), ")\n")

  # Find cells with NO targeting guides AND at least one NT guide
  all_cell_idx <- intersect(
    which(Matrix::colSums(scep_assn_mat[targeting_guides, ,drop=FALSE]) == 0),
    which(Matrix::colSums(scep_assn_mat[nt_guides, ,drop = FALSE]) > 0)
  )

  cat("Found", length(all_cell_idx), "cells expressing no targeting guides and at least one NT guide.\n")

  if (length(all_cell_idx) == 0) {
    stop("No cells found expressing only NT guides")
  }

  return(all_cell_idx)
}


#' Filter all objects consistently based on expression
#'
#' @param response_matrix Sparse matrix (genes × cells)
#' @param grna_indicator Sparse binary matrix (guides × cells)
#' @param grna_target_df Data frame with grna_id and grna_target columns
#' @param cell_idx Integer vector of cell indices
#' @return List with filtered: response_matrix, grna_indicator, grna_target_df, cell_idx, genes
filter_all_objects <- function(response_matrix, grna_indicator, grna_target_df, cell_idx) {
  # Filter cells with all-0 UMI in response matrix
  cells_to_keep <- Matrix::colSums(response_matrix) > 0
  response_matrix <- response_matrix[, cells_to_keep, drop=FALSE]
  grna_indicator <- grna_indicator[, cells_to_keep, drop=FALSE]
  cell_idx <- cell_idx[cells_to_keep]

  # Filter genes with all-0 expression in response matrix
  genes_to_keep <- Matrix::rowSums(response_matrix) > 0
  response_matrix <- response_matrix[genes_to_keep, , drop=FALSE]
  genes <- rownames(response_matrix)

  # Filter guides with all-0 in gRNA indicator matrix
  guides_to_keep <- Matrix::rowSums(grna_indicator) > 0
  grna_indicator <- grna_indicator[guides_to_keep, , drop=FALSE]

  # Update grna_target_df to match kept guides
  grna_target_df <- grna_target_df[grna_target_df$grna_id %in% rownames(grna_indicator), ]

  # Validation
  stopifnot(ncol(response_matrix) == ncol(grna_indicator))
  stopifnot(ncol(response_matrix) == length(cell_idx))
  stopifnot(nrow(grna_indicator) == nrow(grna_target_df))
  stopifnot(all(Matrix::colSums(response_matrix) > 0))
  stopifnot(all(Matrix::rowSums(response_matrix) > 0))
  stopifnot(all(Matrix::colSums(grna_indicator) == 1))  # Low MOI already enforced
  stopifnot(all(Matrix::rowSums(grna_indicator) > 0))

  return(list(
    response_matrix = response_matrix,
    grna_indicator = grna_indicator,
    grna_target_df = grna_target_df,
    cell_idx = cell_idx,
    genes = genes
  ))
}


#' Create pseudo-target names for NT guides
#'
#' @param guide_ids Character vector of guide IDs
#' @return Data frame with grna_id and grna_target columns
create_pseudo_targets <- function(guide_ids) {
  nt_grna_target_df <- data.frame(
    grna_id = guide_ids,
    grna_target = paste0("non-targeting", seq_along(guide_ids)),
    stringsAsFactors = FALSE
  )

  cat("Created", nrow(nt_grna_target_df), "pseudo-targets.\n")

  stopifnot(!any(duplicated(nt_grna_target_df$grna_target)))

  return(nt_grna_target_df)
}


#' Prepare cell metadata and covariates
#'
#' @param grna_indicator_matrix Sparse binary matrix (guides × cells)
#' @param grna_target_df Data frame with grna_id and grna_target mapping
#' @param cell_covariates Data frame of all cell covariates
#' @param cell_idx Integer vector of cell indices to include
#' @return List with: cell_info (data frame), cell_covariates_subset (data frame)
prepare_cell_metadata <- function(grna_indicator_matrix, grna_target_df, cell_covariates, cell_idx) {
  # Create cell_info with guide assignments
  cell_info <- data.frame(
    cell_idx = cell_idx,
    expressed_guide = rownames(grna_indicator_matrix)[
      apply(grna_indicator_matrix, 2, function(col) which(col == 1))
    ],
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(grna_target_df, by = c("expressed_guide" = "grna_id")) |>
    dplyr::mutate(cell_name = rownames(cell_covariates)[cell_idx])

  # Subset and prepare covariates (NO BATCH for Replogle)
  cell_covariates_subset <- cell_covariates[cell_idx, ] |>
    dplyr::transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full = response_n_umis,
      grna_n_nonzero_full = grna_n_nonzero,
      grna_n_umis_full = grna_n_umis,
      batch = batch
    )

  return(list(
    cell_info = cell_info,
    cell_covariates_subset = cell_covariates_subset
  ))
}


#' Main negative control data preparation function
#'
#' @param dataset_name Name for the output dataset
#' @param response_odm OnDisk Matrix of gene expression
#' @param grna_odm OnDisk Matrix of gRNA expression (not used for neg control)
#' @param cell_covariates Data frame of cell-level covariates
#' @param scep_assn_mat Sparse binary matrix of guide assignments
#' @param grna_target_df Data frame mapping guide IDs to targets
#' @param num_genes Number of genes to randomly sample
#' @param include_batch Whether to include batch as a covariate
#' @param random_seed Random seed for reproducibility
make_neg_control_replogle <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    num_genes = 5000,
    include_batch = FALSE,
    random_seed = 243535
) {
  # Negative control data preparation for Replogle (LOW MOI):
  # 1. Select genes randomly
  # 2. Filter to cells expressing NO targeting guides and at least one NT guide
  # 3. Subset and filter response matrix
  # 4. Enforce low MOI (single guide per cell)
  # 5. Create pseudo-targets for each NT guide
  # 6. Prepare cell metadata
  # 7. Write outputs in SCEPTRE, Mixscale (all targets), and FR-Perturb formats

  # Section 0: Setup
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  # Section 1: Gene selection
  genes <- select_genes_random(response_odm, num_genes, random_seed)

  # Section 2: Cell filtering (NT-only cells)
  cat("\n=== Cell filtering ===\n")
  all_cell_idx <- identify_nt_only_cells(scep_assn_mat, grna_target_df)

  # Section 3: Response matrix subsetting (NO FILTERING YET)
  cat("\n=== Extracting response matrix ===\n")
  response_subset <- odm_to_sparse_matrix(response_odm, genes, all_cell_idx)
  cat("Extracted:", nrow(response_subset), "genes ×", ncol(response_subset), "cells (unfiltered)\n")

  # Section 4: gRNA indicator matrix and low MOI enforcement
  cat("\n=== Creating gRNA indicator and enforcing low MOI ===\n")
  nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == "non-targeting"]
  grna_indicator <- scep_assn_mat[nt_guides, all_cell_idx] + 0L
  cat("Created gRNA indicator:", nrow(grna_indicator), "guides ×", ncol(grna_indicator), "cells\n")

  grna_indicator <- enforce_single_guide_per_cell(grna_indicator, random_seed)
  cat("Enforced single guide per cell (low MOI)\n")

  # Section 5: Pseudo-targets
  cat("\n=== Creating pseudo-targets ===\n")
  nt_grna_target_df <- create_pseudo_targets(rownames(grna_indicator))

  # Section 6: Filter all objects consistently
  cat("\n=== Filtering all objects ===\n")
  cat("Before filtering:", nrow(response_subset), "genes ×", ncol(response_subset),
      "cells ×", nrow(grna_indicator), "guides\n")

  filtered <- filter_all_objects(
    response_matrix = response_subset,
    grna_indicator = grna_indicator,
    grna_target_df = nt_grna_target_df,
    cell_idx = all_cell_idx
  )

  response_subset <- filtered$response_matrix
  grna_indicator <- filtered$grna_indicator
  nt_grna_target_df <- filtered$grna_target_df
  all_cell_idx <- filtered$cell_idx
  genes <- filtered$genes

  cat("After filtering:", nrow(response_subset), "genes ×", ncol(response_subset),
      "cells ×", nrow(grna_indicator), "guides\n")

  # Section 7: Cell metadata (after filtering)
  cat("\n=== Preparing cell metadata ===\n")
  metadata <- prepare_cell_metadata(
    grna_indicator,
    nt_grna_target_df,
    cell_covariates,
    all_cell_idx
  )
  cell_info <- metadata$cell_info
  cell_covariates_subset <- metadata$cell_covariates_subset

  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))

  # Section 7: Create discovery pairs (Cartesian product)
  discovery_pairs <- expand.grid(
    grna_target = nt_grna_target_df$grna_target,
    response_id = genes,
    stringsAsFactors = FALSE
  )

  # Section 8: Write outputs
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/neg-control/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)

  # Write shared metadata files
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")

  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"),
            row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")

  # Formula object
  formula_string <- "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1)"
  if(include_batch) {
    formula_string <- paste(formula_string, "+ batch")
  } else {
    # it was still included when they were saved to disc,
    # but not included in any of the model specific datasets
    cell_covariates_subset$batch <- NULL
  }

  # needed for mixscale and frperturb
  cell_names <- paste0("cell_idx_", cell_info$cell_idx)

  # Write SCEPTRE format
  write_sceptre_output(
    response_matrix = response_subset,
    grna_matrix = grna_indicator,
    cell_covariates = cell_covariates_subset,
    grna_target_df = nt_grna_target_df,
    discovery_pairs = discovery_pairs,
    formula_object = formula_string,
    output_path = file.path(write_fp, "sceptre")
  )

  # Write Mixscale format
  write_mixscale_output(
    response_matrix = response_subset,
    cell_info = cell_info,
    output_path = file.path(write_fp, "mixscale"),
    all_targets = nt_grna_target_df$grna_target
  )


  # Prepare FR-Perturb covariates (Replogle: log1p all covariates)
  # regardless of whether batch is used, do not log it
  covariates_to_log1p <- setdiff(names(cell_covariates_subset), "batch")
  cell_covariates_frpert_with_perturbation <- prepare_frperturb_covariates(
    cell_covariates = cell_covariates_subset,
    grna_targets = cell_info$grna_target,
    covariates_to_log1p = covariates_to_log1p
  )

  write_frperturb_output(
    response_matrix = response_subset,
    cell_names = cell_names,
    cell_covariates_frpert = cell_covariates_frpert_with_perturbation,
    output_path = file.path(write_fp, "frperturb")
  )

  cat("\nNegative control dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")
}
