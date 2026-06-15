# Threshold Poisson GLM offset (glm.fit Poisson only on cells with g <= Y_MAX=5000),
# pois0-nb mixture: (1-pi) Pois(mu0) + pi NB(mu0*exp(gamma), theta), with cells
# g == 0 forced non-pert (positive-cells-only EM, faster). Poisson non-pert
# component, so no phi/estimate_phi_fn is needed. pi_max=0.5 guard inside.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  Y_MAX <- 5000
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
    family              = "pois0-nb"
  )
}
