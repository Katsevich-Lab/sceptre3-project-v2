# grnafix offset (circularity-corrected glmpois: log1p(grna_n_umis - y) +
# log1p(grna_n_nonzero - 1{y>0})) with the standard Poisson mixture using the
# additive (bugged) M-step:  g_mus_pert1 = g_mus_pert0 + curr_g_pert.
# Both the offset (cf. glmpoisgrnafix_poisthresh) and the poisbug mixture (cf.
# glmpois_poisbug) have been used before -- this is the two combined.

source(file.path(bin_dir, "script", "lib", "run_grnafix_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_grnafix_variant(
    grna_matrix            = grna_matrix,
    cpus                   = cpus,
    family                 = "pois",
    fix_curr_g_pert_bug    = FALSE,   # additive (bugged) M-step
    n_nonzero_cells_cutoff = 10L,
    backup_threshold       = 5L,
    probability_threshold  = 0.8
  )
}
