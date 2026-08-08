# grnafix offset + a standard sceptre mixture (glmpoisgrnafix_<mixture>).
#
# The "grnafix" offset is the circularity-corrected glmpois baseline already used
# by poisthresh_assign.R: per guide, the guide's own contribution is removed from
# the per-cell totals before fitting a Poisson GLM (log link):
#     other_umis_i = grna_n_umis_i    - y_i           (subtract own UMIs)
#     other_nnz_i  = grna_n_nonzero_i - 1{y_i > 0}     (subtract own presence)
#     log(mu0_i) = b0 + b1 log1p(other_umis_i) + b2 log1p(other_nnz_i)
#
# The mixture is whatever `assign_one_grna_pure_R` implements for `family` (here
# the standard two-component Poisson mixture, "pois"). Both pieces exist already;
# the only new thing is feeding the per-guide circularity design matrix in as the
# covariate_matrix, so that fit_baseline_glm_pure_R(g, X_guide) reproduces the
# grnafix offset and the unchanged EM runs on top of it.
#
# `assign_one_grna_pure_R` and `fit_baseline_glm_pure_R` are sourced from
# sceptre_assign_pure_R.R by the run boilerplate.

# Per-guide circularity-corrected design matrix (intercept + log1p grna covars).
grnafix_build_X <- function(other_umis, other_nnz) {
  X <- cbind(1, log1p(other_umis), log1p(other_nnz))
  colnames(X) <- c("(Intercept)", "log1p(other_grna_n_umis)", "log1p(other_grna_n_nonzero)")
  X
}

# Driver. `grna_n_nonzero` / `grna_n_umis` are the per-cell TOTALS over all guides
# (this guide included); each worker removes the guide's own contribution before
# fitting the offset, then runs the requested sceptre mixture on top of it.
grnafix_assign <- function(grna_matrix,
                           grna_n_nonzero,
                           grna_n_umis,
                           pi_guesses,
                           g_pert_guesses,
                           grna_ids               = rownames(grna_matrix),
                           family                 = "pois",
                           fix_curr_g_pert_bug    = FALSE,
                           n_nonzero_cells_cutoff = 10L,
                           backup_threshold       = 5L,
                           probability_threshold  = 0.8,
                           keep_fits              = FALSE,
                           cl                     = NULL) {
  force(grna_matrix); force(grna_n_nonzero); force(grna_n_umis)
  force(pi_guesses); force(g_pert_guesses); force(family); force(fix_curr_g_pert_bug)
  force(n_nonzero_cells_cutoff); force(backup_threshold); force(probability_threshold)
  force(keep_fits)

  worker <- function(grna_id) {
    g          <- as.numeric(grna_matrix[grna_id, ])
    other_umis <- pmax(grna_n_umis    - g,                  0)   # - y_i
    other_nnz  <- pmax(grna_n_nonzero - as.numeric(g > 0),  0)   # - 1{y_i > 0}
    X          <- grnafix_build_X(other_umis, other_nnz)
    assign_one_grna_pure_R(
      g                      = g,
      covariate_matrix       = X,
      pi_guesses             = pi_guesses,
      g_pert_guesses         = g_pert_guesses,
      n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
      backup_threshold       = backup_threshold,
      probability_threshold  = probability_threshold,
      keep_fits              = keep_fits,
      fix_curr_g_pert_bug    = fix_curr_g_pert_bug,
      family                 = family,
      offset_model_fit_fn    = fit_baseline_glm_pure_R
    )
  }

  per_guide <- if (is.null(cl)) lapply(grna_ids, worker) else parallel::parLapplyLB(cl, grna_ids, worker)
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  run_meta <- list(
    model                 = "glmpoisgrnafix-offset + sceptre mixture",
    family                = family,
    offset_covariates     = "log1p(grna_n_umis - y), log1p(grna_n_nonzero - 1{y>0})",
    fix_curr_g_pert_bug   = fix_curr_g_pert_bug,
    n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold
  )

  list(per_guide = per_guide, run_meta = run_meta)
}
