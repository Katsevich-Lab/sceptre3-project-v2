# Like baseline_nb (vanilla Poisson GLM offset + NB-NB mixture, initial phi
# from sceptre:::estimate_theta), but with one phi-update round:
#   EM -> update phi -> EM
# The saved (pi, g_pert) come from the second EM, matching the updated phi.
# Use n_phi_updates = K (in run_meta) for K updates; here K = 1.

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
    ),
    worker_libraries    = c("Matrix", "sceptre"),
    family              = "nb-shared",
    estimate_phi_fn     = estimate_phi_from_offset_fit_sceptre,
    n_phi_updates       = 1L,
    # Fallback used when initial estimate_phi_fn errors / returns a bad value
    # for a given guide. Recorded per-guide as em_phi_source = "fallback".
    phi                 = 5
  )
}
