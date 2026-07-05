# lib/regime.R
# Two ORTHOGONAL axes of variation, kept deliberately separate:
#
#   * MOI regime  (moi_high / moi_low)  -- genuinely MOI-driven knobs:
#       cell-info structure (concatenated vs single guide/target per cell),
#       enforce-single-guide, whether mixscale runs (needs one target per
#       cell -> low-MOI only), and the per-cell perturbation string.
#
#   * Dataset config (dataset_gasperini / dataset_replogle) -- dataset-driven
#       knobs: assignment-method directory, whether an exclusion list is
#       required, batch-covariate name + whether it is modeled, cell-name
#       policy (real barcodes vs synthesized cell_idx_N), and a pointer to
#       the MOI regime that dataset happens to use.
#
# They coincide today (Gasperini<->high-MOI, Replogle<->low-MOI) but are not
# the same thing; splitting them keeps each knob attached to its real cause.
#
# Also holds the shared covariate-subsetting + enforce helpers and the two
# cell-info builders (ported from the legacy metadata functions).


# ===========================================================================
# Shared mechanics
# ===========================================================================

# ---- subset_cell_covariates -----------------------------------------------
# The 4 `_full` covariates are identical across builders; only the batch column
# NAME is dataset-specific (prep_batch vs batch), so it's a parameter. The batch
# column is always carried into the written cell_covariates.csv; whether it
# enters the formula is a separate dataset policy (model_includes_batch).
subset_cell_covariates <- function(cell_covariates, cell_idx, batch_col) {
  stopifnot(batch_col %in% names(cell_covariates))
  out <- cell_covariates[cell_idx, ] |>
    dplyr::transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full    = response_n_umis,
      grna_n_nonzero_full     = grna_n_nonzero,
      grna_n_umis_full        = grna_n_umis
    )
  out[[batch_col]] <- cell_covariates[[batch_col]][cell_idx]
  out
}


# ---- enforce_single_guide_per_cell ----------------------------------------
# Ported from legacy utils_data_prep.R. With "maximum" assignment every retained
# Replogle cell already has exactly one guide, so this is a validated no-op; kept
# as a guard.
enforce_single_guide_per_cell <- function(grna_indicator_matrix, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)

  result <- grna_indicator_matrix
  num_expressed_guides <- Matrix::colSums(result)
  if (any(num_expressed_guides == 0)) {
    stop("Some cells have no expressed guides! Should not happen here.")
  }

  idx_to_enforce <- which(num_expressed_guides > 1)
  for (cell_idx in idx_to_enforce) {
    expressed_guides <- which(result[, cell_idx] == 1)
    guide_to_keep <- sample(expressed_guides, 1)
    guides_to_remove <- setdiff(expressed_guides, guide_to_keep)
    result[guides_to_remove, cell_idx] <- 0
  }

  stopifnot(all(Matrix::colSums(result) == 1))
  return(result)
}


# ===========================================================================
# MOI-specific: cell-info builders (structure of guide/target per cell)
# ===========================================================================

# ---- build_cell_info_low_moi ----------------------------------------------
# One guide -> one target per cell. Adds the bug-5 fix (assert unique grna_id
# before the join). Uniform signature; concat_string is ignored.
build_cell_info_low_moi <- function(grna_indicator, grna_target_df,
                                    cell_covariates, cell_idx, concat_string = NULL) {
  stopifnot(all(Matrix::colSums(grna_indicator) == 1))
  stopifnot(!any(duplicated(grna_target_df$grna_id)))   # bug-5 guard

  data.frame(
    cell_idx = cell_idx,
    expressed_guide = rownames(grna_indicator)[
      apply(grna_indicator, 2, function(col) which(col == 1))
    ],
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(grna_target_df, by = c("expressed_guide" = "grna_id")) |>
    dplyr::mutate(cell_name = rownames(cell_covariates)[cell_idx])
}


# ---- build_cell_info_high_moi ---------------------------------------------
# Possibly several guides/targets per cell -> concatenated strings. Ported from
# the legacy prepare_cell_metadata_high_moi (dgCMatrix column-pointer walk).
build_cell_info_high_moi <- function(grna_indicator, grna_target_df,
                                     cell_covariates, cell_idx, concat_string) {
  if (!inherits(grna_indicator, "dgCMatrix")) {
    grna_indicator <- as(grna_indicator, "dgCMatrix")
  }

  guide_ids <- rownames(grna_indicator)
  if (is.null(guide_ids) || !all(guide_ids %in% grna_target_df$grna_id)) {
    stop("`grna_indicator` must have rownames giving guide IDs.")
  }

  target_by_guide <- setNames(as.character(grna_target_df$grna_target), grna_target_df$grna_id)
  target_ids <- unname(target_by_guide[guide_ids])
  if (anyNA(target_ids)) {
    missing_guides <- unique(guide_ids[is.na(target_ids)])
    stop("Some guides in `grna_indicator` are missing from `grna_target_df`: ",
         paste(head(missing_guides, 10), collapse = ", "))
  }

  p <- grna_indicator@p
  i <- grna_indicator@i + 1L
  if (any(diff(p) == 0L)) {
    stop("In build_cell_info_high_moi(), some cells have no assigned guides, which should not happen here.")
  }

  n_cells <- ncol(grna_indicator)
  expressed_guide_concat  <- character(n_cells)
  expressed_target_concat <- character(n_cells)
  for (col in seq_len(n_cells)) {
    rows <- i[seq.int(p[col] + 1L, p[col + 1L])]
    expressed_guide_concat[col]  <- paste(guide_ids[rows],  collapse = concat_string)
    expressed_target_concat[col] <- paste(target_ids[rows], collapse = concat_string)
  }

  data.frame(
    cell_idx = cell_idx,
    expressed_guide_concat,
    expressed_target_concat,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(cell_name = rownames(cell_covariates)[cell_idx])
}


# ===========================================================================
# MOI regimes
# ===========================================================================

moi_high <- list(
  name                 = "high_moi",
  low_moi              = FALSE,
  enforce_single_guide = FALSE,
  writes_mixscale      = FALSE,   # mixscale needs one target per cell
  # NT is NOT force-included: no high-MOI method needs an NT population (sceptre
  # uses complement; FR-Perturb centers on the grand average). A cell carrying a
  # sampled target guide may incidentally also express NT, but it enters via its
  # target guide -- NT guides aren't recorded and NT-only cells aren't pulled in.
  # (pos-Gasperini once force-included NT to boost its small N; re-enable via a
  #  per-selector force_nt override on that driver, not here.)
  force_nt_inclusion   = FALSE,
  build_cell_info      = build_cell_info_high_moi,
  perturbation_of      = function(cell_info) cell_info$expressed_target_concat
)

moi_low <- list(
  name                 = "low_moi",
  low_moi              = TRUE,
  enforce_single_guide = TRUE,    # validated no-op under "maximum"
  writes_mixscale      = TRUE,
  # Low-MOI: one guide per cell, so NT cells are DISTINCT cells that must be
  # pulled in during selection. sceptre (nt_cells), mixscale (nt.cell.class) and
  # FR-Perturb (--control-perturbation-name) all require an NT population.
  force_nt_inclusion   = TRUE,
  build_cell_info      = build_cell_info_low_moi,
  perturbation_of      = function(cell_info) cell_info$grna_target
)


# ===========================================================================
# Dataset configs
# ===========================================================================

.formula_base <-
  "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1)"

#' Model formula for a dataset (batch term appended only if that dataset models it).
dataset_formula <- function(dataset) {
  if (dataset$model_includes_batch) paste0(.formula_base, " + ", dataset$batch_col) else .formula_base
}

dataset_gasperini <- list(
  name                 = "gasperini",
  assignment_method    = "thresholding-5",
  excluded_required    = FALSE,
  moi                  = moi_high,
  batch_col            = "prep_batch",
  model_includes_batch = TRUE,
  frperturb_cell_names = function(cell_info) cell_info$cell_name   # real barcodes
)

dataset_replogle <- list(
  name                 = "replogle-rd7",
  assignment_method    = "maximum",
  excluded_required    = TRUE,
  moi                  = moi_low,
  batch_col            = "batch",
  model_includes_batch = FALSE,
  frperturb_cell_names = function(cell_info) paste0("cell_idx_", cell_info$cell_idx)  # no barcodes
)
