# Robust Poisson GLM offset (robustbase::glmrob, MT), Poisson mixture,
# multiplicative M-step (fixed).

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
      description = "Robust Poisson GLM via robustbase::glmrob (MT)",
      params      = list()
    ),
    worker_libraries    = c("Matrix", "robustbase"),
    fix_curr_g_pert_bug = TRUE
  )
}
