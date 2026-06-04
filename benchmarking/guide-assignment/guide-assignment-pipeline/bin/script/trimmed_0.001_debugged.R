# Trimmed Poisson GLM offset (top 0.1% of cells by g dropped before fitting)
# combined with the curr_g_pert M-step bug fix turned ON (multiplicative
# perturbation effect: g_mus_pert1 = g_mus_pert0 * exp(curr_g_pert)).

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
    fix_curr_g_pert_bug = TRUE
  )
}
