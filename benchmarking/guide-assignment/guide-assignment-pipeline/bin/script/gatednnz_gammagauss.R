# Gated logistic / Gamma-Gaussian mixture (no offset model).
#
#   x_i  = log(grna_n_nonzero_i)              # gate covariate
#   z_i  = log(y_i), fit on y_i >= 2          # y in {0,1} -> non-pert
#   pi_i = ilogit(a0 + a1 x_i)
#   z_i ~ (1 - pi_i) Gamma(alpha, beta) + pi_i Normal(mu, sigma)   # Gaussian = high
#
# Assign pert iff posterior P(pert_i | data) >= 0.8.

source(file.path(bin_dir, "script", "lib", "run_gated_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_gated_variant(
    grna_matrix           = grna_matrix,
    cpus                  = cpus,
    n_em_rep              = 5L,
    n_fit_cutoff          = 10L,
    backup_threshold      = 5L,
    probability_threshold = 0.8,
    max_iter              = 200L,
    tol                   = 1e-6
  )
}
