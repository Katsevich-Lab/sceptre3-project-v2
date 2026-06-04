# Trimmed Poisson GLM offset (top 0.1% of cells by g dropped before fitting),
# combined with the NB-NB mixture EM that estimates class-specific
# overdispersions (phi0 != phi1). Initial seed phi (used to start both phi0
# and phi1) is estimated per-guide from the trimmed Poisson fit via
# sceptre's estimate_theta; fallback phi = 5 when that errors.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  TRIM_FRAC <- 0.1 / 100
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
    worker_libraries    = c("Matrix", "sceptre"),
    family              = "nb-separate",
    estimate_phi_fn     = estimate_phi_from_offset_fit_sceptre,
    # Fallback used when the initial estimate_phi_fn errors or returns a bad
    # value. Recorded per-guide as em_phi_source = "fallback".
    phi                 = 5
  )
}
