# allgrnafix offset (circularity-corrected grna covars log1p(grna_n_umis - y) +
# log1p(grna_n_nonzero - 1{y>0}), PLUS the two uncorrected response covars
# log1p(response_n_umis_full) + log1p(response_n_nonzero_full)) with the standard
# Poisson mixture using the additive (bugged) M-step:
#   g_mus_pert1 = g_mus_pert0 + curr_g_pert.
# Same as glmpoisgrnafix_poisbug but the offset also keeps the response covariates
# (no circularity correction needed there -- they come from the gene matrix, not
# the grna matrix).

source(file.path(bin_dir, "script", "lib", "run_allgrnafix_variant.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  run_allgrnafix_variant(
    grna_matrix            = grna_matrix,
    extra_covariates       = extra_covariates,
    cpus                   = cpus,
    family                 = "pois",
    fix_curr_g_pert_bug    = FALSE,   # additive (bugged) M-step
    n_nonzero_cells_cutoff = 10L,
    backup_threshold       = 5L,
    probability_threshold  = 0.8
  )
}
