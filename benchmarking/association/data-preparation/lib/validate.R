# lib/validate.R
# Always-on invariant checks. validate_dataset() runs unconditionally inside
# build_dataset() right before writing, so a malformed dataset aborts loudly
# instead of being written. These are belt-and-suspenders on top of the filter
# step -- cheap, and they encode the contract every output must satisfy.


# ---- filter_all_objects ----------------------------------------------------
# Consistency filter: drop all-zero-response cells, all-zero-response genes, and
# all-zero guides, keeping every object aligned. Ported from legacy
# neg_control_functions.R; run UNCONDITIONALLY by build_dataset() (legacy only
# called it for neg-control -- that omission was bug #1/#2).
filter_all_objects <- function(response_matrix, grna_indicator, grna_target_df, cell_idx, low_moi = TRUE) {
  cells_to_keep <- Matrix::colSums(response_matrix) > 0
  response_matrix <- response_matrix[, cells_to_keep, drop = FALSE]
  grna_indicator  <- grna_indicator[, cells_to_keep, drop = FALSE]
  cell_idx        <- cell_idx[cells_to_keep]

  genes_to_keep   <- Matrix::rowSums(response_matrix) > 0
  response_matrix <- response_matrix[genes_to_keep, , drop = FALSE]
  genes           <- rownames(response_matrix)

  guides_to_keep  <- Matrix::rowSums(grna_indicator) > 0
  grna_indicator  <- grna_indicator[guides_to_keep, , drop = FALSE]

  grna_target_df  <- grna_target_df[grna_target_df$grna_id %in% rownames(grna_indicator), ]

  stopifnot(ncol(response_matrix) == ncol(grna_indicator))
  stopifnot(ncol(response_matrix) == length(cell_idx))
  stopifnot(nrow(grna_indicator) == nrow(grna_target_df))
  stopifnot(all(Matrix::colSums(response_matrix) > 0))
  stopifnot(all(Matrix::rowSums(response_matrix) > 0))
  stopifnot(all(Matrix::rowSums(grna_indicator) > 0))
  if (low_moi) {
    stopifnot(all(Matrix::colSums(grna_indicator) == 1))
  }

  list(
    response_matrix = response_matrix,
    grna_indicator  = grna_indicator,
    grna_target_df  = grna_target_df,
    cell_idx        = cell_idx,
    genes           = genes
  )
}


# ---- validate_dataset ------------------------------------------------------
#' Assert every invariant a finished dataset must satisfy, before writing.
validate_dataset <- function(response_matrix, grna_indicator, grna_target_df,
                             cell_info, cell_covariates_subset, cell_idx,
                             discovery_pairs, dataset, excluded_idx = integer(0)) {
  n_cells <- ncol(response_matrix)

  # --- alignment across all cell-indexed objects ---------------------------
  stopifnot(ncol(grna_indicator) == n_cells)
  stopifnot(nrow(cell_covariates_subset) == n_cells)
  stopifnot(nrow(cell_info) == n_cells)
  stopifnot(length(cell_idx) == n_cells)
  stopifnot(identical(cell_info$cell_idx, cell_idx))

  # --- guides <-> target df ------------------------------------------------
  stopifnot(nrow(grna_indicator) == nrow(grna_target_df))
  stopifnot(setequal(rownames(grna_indicator), grna_target_df$grna_id))

  # --- no degenerate rows/cols --------------------------------------------
  stopifnot(all(Matrix::rowSums(response_matrix) > 0))
  stopifnot(all(Matrix::colSums(response_matrix) > 0))
  stopifnot(all(Matrix::rowSums(grna_indicator) > 0))
  if (dataset$moi$low_moi) {
    stopifnot(all(Matrix::colSums(grna_indicator) == 1))
  }

  # --- cell indices are valid + no excluded cell leaked in -----------------
  stopifnot(!any(duplicated(cell_idx)))
  stopifnot(all(cell_idx >= 1L))
  stopifnot(length(intersect(cell_idx, excluded_idx)) == 0)

  # --- discovery pairs are a clean cartesian of SURVIVING targets x genes --
  if (!is.null(discovery_pairs)) {
    stopifnot(all(discovery_pairs$grna_target %in% grna_target_df$grna_target))
    stopifnot(all(discovery_pairs$response_id %in% rownames(response_matrix)))
    n_targets <- length(unique(discovery_pairs$grna_target))
    stopifnot(nrow(discovery_pairs) == n_targets * nrow(response_matrix))
  }

  invisible(TRUE)
}
