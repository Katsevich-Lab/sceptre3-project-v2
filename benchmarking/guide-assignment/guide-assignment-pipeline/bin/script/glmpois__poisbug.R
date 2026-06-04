# Vanilla Poisson GLM offset, Poisson mixture, additive M-step (bugged):
# g_mus_pert1 = g_mus_pert0 + curr_g_pert.

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
    )
  )
}
