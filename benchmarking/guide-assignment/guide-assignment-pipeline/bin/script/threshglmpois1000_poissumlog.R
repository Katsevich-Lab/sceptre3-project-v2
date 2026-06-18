# Threshold Poisson GLM offset (glm.fit Poisson only on cells with g <= Y_MAX=1000),
# sum-of-Poissons additive-signal mixture (family = "poissum") with
# use_log_gamma = TRUE: the model is fit with signal mean exp(gamma), but cells
# are CLASSIFIED using gamma itself as the additive signal mean (mu1 = mu0 +
# gamma), mimicking the log-scale pathology.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  Y_MAX <- 1000
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_pois_threshold_pure_R(g, X, y_max = Y_MAX)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glm_pois_threshold_pure_R",
      description = "Poisson GLM (glm.fit) fit only on cells with g <= y_max",
      params      = list(y_max = Y_MAX)
    ),
    family              = "poissum",
    use_log_gamma       = TRUE
  )
}
