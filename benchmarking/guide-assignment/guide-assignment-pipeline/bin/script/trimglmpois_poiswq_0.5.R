# Trimmed Poisson GLM offset (top 0.1% by g dropped), multiplicative Pois-Pois
# mixture with a robust weighted-quantile gamma M-step. gamma_update_prob = 0.5
# uses the weighted MEDIAN of y/mu0 in place of the (outlier-sensitive)
# weighted mean.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  TRIM_FRAC <- 0.001
  GAMMA_UPDATE_PROB <- 0.5
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_trimmed_pure_R(g, X, trim_frac = TRIM_FRAC)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glm_trimmed_pure_R",
      description = "Poisson MLE GLM fit on cells outside the top trim_frac of g",
      params      = list(trim_frac = TRIM_FRAC)
    ),
    family              = "pois-wq",
    gamma_update_prob   = GAMMA_UPDATE_PROB
  )
}
