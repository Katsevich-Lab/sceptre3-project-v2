# Joint Poisson-Poisson mixture (no offset model), covariate = grna_n_nonzero.
#
#   Z_i  ~ Bernoulli(pi)                               # perturbation indicator
#   Y_i | Z_i ~ Poisson(mu_i)
#   log(mu_i) = b0 + b1 Z_i + b2 log1p(grna_n_nonzero_i - 1{y_i > 0})
#
# Four parameters: pi, b0, b1 (>= 0, pert effect), b2 (covariate coef). Fit
# jointly via EM over all cells; assign pert iff P(Z_i = 1 | data) >= 0.8.

source(file.path(bin_dir, "script", "lib", "run_poismix_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_poismix_variant(
    grna_matrix           = grna_matrix,
    cpus                  = cpus,
    n_em_rep              = 5L,
    n_fit_cutoff          = 10L,
    backup_threshold      = 5L,
    probability_threshold = 0.8,
    max_iter              = 200L,
    tol                   = 1e-8
  )
}
