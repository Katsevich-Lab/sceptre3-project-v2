# Trimmed NB GLM offset (top 0.1% by g dropped, then MASS::glm.nb),
# Poisson mixture, additive M-step (bugged).

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  TRIM_FRAC <- 0.001
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_nb_trimmed_pure_R(g, X, trim_frac = TRIM_FRAC)
  }
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = formula,
    cpus                = cpus,
    offset_model_fit_fn = offset_model_fit_fn,
    offset_spec         = list(
      name        = "fit_baseline_glm_nb_trimmed_pure_R",
      description = "NB GLM via MASS::glm.nb fit on cells outside the top trim_frac of g",
      params      = list(trim_frac = TRIM_FRAC)
    ),
    worker_libraries    = c("Matrix", "MASS", "sceptre")
  )
}
