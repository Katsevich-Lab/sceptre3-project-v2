# Single-predictor Poisson GLM offset: g ~ log1p(grna_n_nonzero) only (plus
# intercept). Standard vanilla Poisson MLE GLM (fit_baseline_glm_pure_R), just
# with the dataset formula overridden to this one covariate. Mixture:
# sum-of-Poissons additive-signal (family = "poissum") with use_log_gamma = TRUE
# (fit with signal mean exp(gamma), classify with gamma; log-scale pathology).

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
    family              = "poissum",
    use_log_gamma       = TRUE
  )
}
