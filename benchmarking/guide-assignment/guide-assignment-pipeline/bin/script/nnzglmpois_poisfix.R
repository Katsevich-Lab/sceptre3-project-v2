# Single-predictor Poisson GLM offset: g ~ log1p(grna_n_nonzero) only (plus
# intercept). Standard vanilla Poisson MLE GLM (fit_baseline_glm_pure_R), just
# with the dataset formula overridden to this one covariate. Mixture: Poisson
# with the multiplicative (fixed) M-step g_mus_pert1 = g_mus_pert0 * exp(curr_g_pert).
#
# `grna_n_nonzero` is the per-cell count of detected gRNAs, added to cov_df in
# run_variant.R (Matrix::colSums(grna_matrix > 0)).

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = ~ log1p(grna_n_nonzero),
    cpus                = cpus,
    offset_model_fit_fn = fit_baseline_glm_pure_R,
    offset_spec         = list(
      name        = "fit_baseline_glm_pure_R",
      description = "Vanilla Poisson MLE GLM via stats::glm.fit",
      params      = list(formula = "~ log1p(grna_n_nonzero)")
    ),
    fix_curr_g_pert_bug = TRUE
  )
}
