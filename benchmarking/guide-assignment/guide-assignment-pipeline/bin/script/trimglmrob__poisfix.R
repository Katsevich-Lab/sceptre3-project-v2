# Trimmed robust Poisson GLM offset (top 0.1% by g dropped, then glmrob Mqle),
# Poisson mixture, multiplicative M-step (fixed).

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  TRIM_FRAC <- 0.001
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glmrob_trimmed_pure_R(g, X, trim_frac = TRIM_FRAC)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glmrob_trimmed_pure_R",
      description = "Robust Poisson GLM via glmrob (Mqle) on cells outside the top trim_frac of g",
      params      = list(trim_frac = TRIM_FRAC)
    ),
    worker_libraries    = c("Matrix", "robustbase"),
    fix_curr_g_pert_bug = TRUE
  )
}
