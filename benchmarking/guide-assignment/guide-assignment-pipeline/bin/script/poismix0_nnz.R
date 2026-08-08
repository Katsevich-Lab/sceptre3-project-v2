# Joint Poisson-Poisson mixture, covariate = grna_n_nonzero, pois0 speedup:
# cells with y_i = 0 are forced to the non-pert component (like pois0nb).
#
#   Z_i  ~ Bernoulli(pi)   (forced 0 when y_i = 0)
#   Y_i | Z_i ~ Poisson(mu_i)
#   log(mu_i) = b0 + b1 Z_i + b2 log1p(grna_n_nonzero_i - 1{y_i > 0})
#
# Same four parameters and joint EM as script_poismix_nnz, but the E-step and
# the positive part of the M-step touch only positive cells; the zeros enter the
# (b0, b2) score in closed form via distinct-nnz bins. Assign pert iff
# P(Z_i = 1 | data) >= 0.8.

source(file.path(bin_dir, "script", "lib", "run_poismix_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_poismix_variant(
    grna_matrix           = grna_matrix,
    cpus                  = cpus,
    variant               = "pois0",
    n_em_rep              = 5L,
    n_fit_cutoff          = 10L,
    backup_threshold      = 5L,
    probability_threshold = 0.8,
    max_iter              = 200L,
    tol                   = 1e-8
  )
}
