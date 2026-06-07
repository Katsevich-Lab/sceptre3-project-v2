# Vanilla Poisson GLM offset, additive ZTP mixture (positive counts only):
# zero cells always get posterior = 0; positive cells fit
# (1-rho)*ZTP(mu0) + rho*ZTP(mu0+delta) with correct ZTP M-step.

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
    family              = "pois-additive-nonzero"
  )
}
