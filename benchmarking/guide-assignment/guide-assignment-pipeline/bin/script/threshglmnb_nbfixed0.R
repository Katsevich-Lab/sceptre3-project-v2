# Threshold NB GLM offset (fit MASS::glm.nb only to cells with g <= Y_MAX,
# default 100), NB-fixed0 mixture: the nonpert size theta0 is fixed at the
# offset fit's theta and only theta1 (perturbed size), gamma, and pi are
# estimated. Seed/fixed phi per guide is theta from the threshold NB GLM fit.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_nb_threshold_pure_R(g, X)   # Y_MAX = 100 (default)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glm_nb_threshold_pure_R",
      description = "NB GLM via MASS::glm.nb fit only on cells with g <= y_max",
      params      = list(y_max = 100)
    ),
    worker_libraries    = c("Matrix", "MASS", "sceptre"),
    family              = "nb-fixed0",
    estimate_phi_fn     = estimate_phi_from_offset_fit_theta,
    # Fallback when the glm.nb-derived theta is bad (e.g. fit didn't converge).
    phi                 = 5
  )
}
