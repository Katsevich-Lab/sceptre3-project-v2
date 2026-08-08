# glmpois offset (circularity-corrected grna_n_umis + grna_n_nonzero) with a
# fixed-shift Poisson mixture:  (1-pi) Pois(mu0) + pi Pois(mu0 + K),  K = 5.
# Only pi is estimated; K is hard-coded. Assign pert iff P(pert | y) >= 0.8.

source(file.path(bin_dir, "script", "lib", "run_poisthresh_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_poisthresh_variant(
    grna_matrix           = grna_matrix,
    cpus                  = cpus,
    K                     = 5,
    n_fit_cutoff          = 10L,
    backup_threshold      = 5L,
    probability_threshold = 0.8,
    max_iter              = 200L,
    tol                   = 1e-8
  )
}
