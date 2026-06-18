# Vanilla Poisson GLM offset, additive NB-signal Poisson mixture
# (family = "nbsum"): y_i ~ (1-pi) Pois(mu0_i) + pi [Pois(mu0_i) +
# NB(mu = exp(gamma), size = theta)]. Fit and classify both use exp(gamma) as
# the NB signal mean.

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
    family              = "nbsum"
  )
}
