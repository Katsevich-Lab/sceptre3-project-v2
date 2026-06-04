# Baseline: pure-R port of sceptre's default mixture-model gRNA assignment.
# Sources the shared implementation from bin/script/lib/.

# `bin_dir` is defined at the top level of bin/run_script.R, which source()s
# this file into the global environment.
IMPL_PATH <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  if (!file.exists(IMPL_PATH)) stop("Implementation not found at: ", IMPL_PATH)
  source(IMPL_PATH)

  # Mirror load_scep_objs()'s pure_R_cov_df augmentation
  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)

  # Offset model: vanilla Poisson MLE GLM (sceptre default).
  offset_model_fit_fn <- fit_baseline_glm_pure_R
  attr(offset_model_fit_fn, "spec") <- list(
    name        = "fit_baseline_glm_pure_R",
    description = "Vanilla Poisson MLE GLM via stats::glm.fit",
    params      = list()
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
    offset_model_fit_fn  = offset_model_fit_fn
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
