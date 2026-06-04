# Trimmed Poisson GLM offset (top 0.1% of cells by g dropped before fitting)
# combined with the curr_g_pert M-step bug fix turned ON (multiplicative
# perturbation effect: g_mus_pert1 = g_mus_pert0 * exp(curr_g_pert)).

IMPL_PATH <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  if (!file.exists(IMPL_PATH)) stop("Implementation not found at: ", IMPL_PATH)
  source(IMPL_PATH)

  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)

  TRIM_FRAC <- 0.1 / 100
  offset_model_fit_fn <- function(g, X) {
    fit_baseline_glm_trimmed_pure_R(g, X, trim_frac = TRIM_FRAC)
  }
  attr(offset_model_fit_fn, "spec") <- list(
    name        = "fit_baseline_glm_trimmed_pure_R",
    description = "Poisson MLE GLM fit on cells outside the top trim_frac of g",
    params      = list(trim_frac = TRIM_FRAC)
  )

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, IMPL_PATH)
    parallel::clusterEvalQ(cl, library(Matrix))
  }

  results <- sceptre_assign_pure_R(
    grna_matrix          = grna_matrix,
    grna_ids             = rownames(grna_matrix),
    covariate_data_frame = cov_df,
    formula_object       = formula,
    cl                   = cl,
    offset_model_fit_fn  = offset_model_fit_fn,
    fix_curr_g_pert_bug  = TRUE
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
