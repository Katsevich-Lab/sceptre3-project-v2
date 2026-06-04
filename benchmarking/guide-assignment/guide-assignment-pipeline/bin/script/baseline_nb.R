# NB-NB mixture EM with vanilla Poisson GLM offset.
# phi (shared overdispersion) is estimated per-guide from the Poisson fit
# via sceptre's internal estimate_theta (Newton on the NB theta score).

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
    worker_libraries    = c("Matrix", "sceptre"),
    family              = "nb-shared",
    estimate_phi_fn     = estimate_phi_from_offset_fit_sceptre,
    # Fallback used when estimate_phi_fn errors / returns a bad value for a
    # given guide. Recorded per-guide as em_phi_source = "fallback".
    phi                 = 5
  )
}
