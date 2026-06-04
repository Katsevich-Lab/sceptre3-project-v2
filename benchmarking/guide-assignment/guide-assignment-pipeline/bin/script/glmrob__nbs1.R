# Robust Poisson GLM offset (robustbase::glmrob, Mqle), NB-shared mixture
# with 1 phi update. Initial phi per guide from sceptre:::estimate_theta on
# the glmrob fit.

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = fit_baseline_glmrob_pure_R,
    offset_spec         = list(
      name        = "fit_baseline_glmrob_pure_R",
      description = "Robust Poisson GLM via robustbase::glmrob (Mqle)",
      params      = list()
    ),
    worker_libraries    = c("Matrix", "robustbase", "sceptre"),
    family              = "nb-shared",
    estimate_phi_fn     = estimate_phi_from_offset_fit_sceptre,
    n_phi_updates       = 1L,
    # Fallback when initial estimate_phi_fn errors / returns a bad value.
    phi                 = 5
  )
}
