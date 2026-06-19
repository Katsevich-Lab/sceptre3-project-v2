# Intercept-only Poisson GLM offset (g ~ 1) with log1p(log1p(grna_n_nonzero))
# as a FIXED offset (coef pinned to 1):
#   log(mu0) = beta0 + log1p(log1p(grna_n_nonzero))
#   mu0      = exp(beta0) * (1 + log1p(grna_n_nonzero)).
# A shape-constrained, robust alternative to the free univariate fit
# exp(beta0) * (grna_n_nonzero + 1)^beta1. Mixture: sum-of-Poissons
# additive-signal (family = "poissum") with use_log_gamma = TRUE (fit with
# signal mean exp(gamma), classify with gamma; log-scale pathology).
#
# `grna_n_nonzero` is the per-cell count of detected gRNAs, added to cov_df in
# run_variant.R (Matrix::colSums(grna_matrix > 0)).

source(file.path(bin_dir, "script", "lib", "run_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_variant(
    grna_matrix         = grna_matrix,
    extra_covariates    = extra_covariates,
    formula             = ~ log1p(log1p(grna_n_nonzero)),
    cpus                = cpus,
    offset_model_fit_fn = fit_baseline_glm_nnz_offset_pure_R,
    offset_spec         = list(
      name        = "fit_baseline_glm_nnz_offset_pure_R",
      description = "Intercept-only Poisson GLM, log1p(log1p(grna_n_nonzero)) fixed offset",
      params      = list(formula = "~ log1p(log1p(grna_n_nonzero))")
    ),
    family              = "poissum",
    use_log_gamma       = TRUE
  )
}
