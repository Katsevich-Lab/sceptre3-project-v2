# Capped Poisson GLM offset (stats::glm.fit Poisson on pmin(g, Y_MAX=100); all
# cells kept, large counts capped), pois-nb mixture: (1-pi) Pois(mu0) +
# pi NB(mu0*exp(gamma), theta). Poisson non-pert component, so no
# phi/estimate_phi_fn is needed. pi_max=0.5 minority guard is applied inside.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_pois_capped_pure_R(g, X)   # Y_MAX = 100 (default)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glm_pois_capped_pure_R",
      description = "Poisson GLM (glm.fit) fit on pmin(g, y_max)",
      params      = list(y_max = 100)
    ),
    family              = "pois-nb"
  )
}
