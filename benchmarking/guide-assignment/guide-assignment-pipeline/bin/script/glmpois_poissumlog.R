# Vanilla Poisson GLM offset, sum-of-Poissons additive-signal mixture
# (family = "poissum") with use_log_gamma = TRUE: the model is fit with signal
# mean exp(gamma), but cells are CLASSIFIED using gamma itself as the additive
# signal mean (mu1 = mu0 + gamma), mimicking the log-scale pathology.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = fit_baseline_glm_pure_R,
    offset_spec         = list(
      name        = "fit_baseline_glm_pure_R",
      description = "Vanilla Poisson MLE GLM via stats::glm.fit",
      params      = list()
    ),
    family              = "poissum",
    use_log_gamma       = TRUE
  )
}
