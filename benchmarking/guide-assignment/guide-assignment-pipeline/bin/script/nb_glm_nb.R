# NB-NB mixture EM with a MASS::glm.nb offset model. Theta is estimated
# jointly with the regression coefficients inside the offset fit, and then
# reused as the shared overdispersion (phi) for the NB-NB EM.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = fit_baseline_glm_nb,
    offset_spec         = list(
      name        = "fit_baseline_glm_nb",
      description = "NB GLM via MASS::glm.nb (joint MLE of beta and theta)",
      params      = list()
    ),
    worker_libraries    = c("Matrix", "MASS", "sceptre"),
    family              = "nb-shared",
    estimate_phi_fn     = estimate_phi_from_offset_fit_theta,
    # Fallback used when the glm.nb-derived theta is bad (e.g. fit didn't
    # converge cleanly). Recorded per-guide as em_phi_source = "fallback".
    phi                 = 5
  )
}
