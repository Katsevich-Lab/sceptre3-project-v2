# Trimmed Poisson GLM offset (top 0.1% of cells by g dropped before fitting),
# combined with the NB-NB mixture EM that estimates class-specific
# overdispersions (phi0 != phi1). Initial seed phi (used to start both phi0
# and phi1) is estimated per-guide from the trimmed Poisson fit via
# sceptre's estimate_theta; fallback phi = 5 when that errors.

IMPL_PATH <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  stopifnot(requireNamespace("sceptre", quietly = TRUE))  # for sceptre:::estimate_theta
  if (!file.exists(IMPL_PATH)) stop("Implementation not found at: ", IMPL_PATH)
  source(IMPL_PATH)

  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)

  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_trimmed_pure_R(g, X, trim_frac = 0.001)
  }
  attr(offset_model_fit_fn, "spec") <- list(
    name        = "fit_baseline_glm_trimmed_pure_R",
    description = "Poisson MLE GLM fit on cells outside the top trim_frac of g",
    params      = list(trim_frac = 0.001)
  )

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, IMPL_PATH)
    parallel::clusterEvalQ(cl, library(Matrix))
    parallel::clusterEvalQ(cl, library(sceptre))
  }

  results <- sceptre_assign_pure_R(
    grna_matrix          = grna_matrix,
    grna_ids             = rownames(grna_matrix),
    covariate_data_frame = cov_df,
    formula_object       = formula,
    cl                   = cl,
    offset_model_fit_fn  = offset_model_fit_fn,
    family               = "nb-separate",
    estimate_phi_fn      = estimate_phi_from_offset_fit_sceptre,
    # Fallback used when the initial estimate_phi_fn errors or returns a bad
    # value. Recorded per-guide as em_phi_source = "fallback".
    phi                  = 5
  )

  per_guide <- results$per_guide
  ns <- vapply(per_guide, function(r) length(r$assignments), integer(1))
  assignment_matrix <- Matrix::sparseMatrix(
    i = rep.int(seq_along(per_guide), ns),
    j = unlist(lapply(per_guide, `[[`, "assignments"), use.names = FALSE),
    x = TRUE,
    dims = c(nrow(grna_matrix), ncol(grna_matrix)),
    dimnames = list(rownames(grna_matrix), colnames(grna_matrix))
  )

  list(assignment_matrix = assignment_matrix, extras = results)
}
