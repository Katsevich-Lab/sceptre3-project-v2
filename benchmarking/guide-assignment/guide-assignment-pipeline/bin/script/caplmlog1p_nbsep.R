# Capped log1p-OLS offset (stats::lm.fit of log1p(pmin(g, Y_MAX=100)); baseline
# mean recovered as expm1(prediction)), NB-separate mixture (phi0 != phi1).
# The OLS offset has no $theta, so the per-guide seed phi (used to initialize
# phi0 = phi1) comes from the means-based estimator; phi = 5 is the fallback.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_lm_log1p_capped_pure_R(g, X)   # Y_MAX = 100 (default)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_lm_log1p_capped_pure_R",
      description = "OLS (lm.fit) of log1p(pmin(g, y_max)); mu0 = expm1(fitted)",
      params      = list(y_max = 100)
    ),
    family              = "nb-separate",
    estimate_phi_fn     = estimate_phi_from_offset_means_nb,
    # Fallback if the means-based seed theta is non-finite/non-positive.
    phi                 = 5
  )
}
