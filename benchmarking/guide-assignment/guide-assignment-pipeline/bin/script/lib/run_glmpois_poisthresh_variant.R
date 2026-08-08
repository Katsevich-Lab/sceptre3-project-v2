# Shared boilerplate for the glmpois_poisthresh<K> variant(s): vanilla glmpois
# offset + fixed-shift Poisson mixture. Unlike run_poisthresh_variant.R (which
# applies a per-guide circularity correction), the offset here is the standard
# glmpois baseline fit on the SHARED covariate matrix built once from the run's
# formula -- exactly the offset used by script_glmpois (glmpois_poisfix.R).
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s the
# variant into the global environment; the variant in turn source()s this file.

POISTHRESH_IMPL_PATH <- file.path(bin_dir, "script", "lib", "poisthresh_assign.R")
SCEPTRE_IMPL_PATH    <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")

run_glmpois_poisthresh_variant <- function(grna_matrix, extra_covariates, formula, cpus, K,
                                           n_fit_cutoff          = 10L,
                                           backup_threshold      = 5L,
                                           probability_threshold = 0.8,
                                           pi_init               = 0.05,
                                           max_iter              = 200L,
                                           tol                   = 1e-8,
                                           worker_libraries      = "Matrix",
                                           keep_fits             = FALSE,
                                           ...) {
  for (pkg in worker_libraries) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        paste0("requireNamespace('%s') failed.\n",
               "  .libPaths()  : %s\n",
               "  R.home()     : %s\n",
               "  R_LIBS_USER  : %s\n",
               "  installed?   : %s\n"),
        pkg,
        paste(.libPaths(), collapse = ":"),
        R.home(),
        Sys.getenv("R_LIBS_USER", "<unset>"),
        pkg %in% rownames(utils::installed.packages())
      ))
    }
  }
  if (!file.exists(SCEPTRE_IMPL_PATH))    stop("Implementation not found at: ", SCEPTRE_IMPL_PATH)
  if (!file.exists(POISTHRESH_IMPL_PATH)) stop("Implementation not found at: ", POISTHRESH_IMPL_PATH)
  source(SCEPTRE_IMPL_PATH)      # fit_baseline_glm_pure_R
  source(POISTHRESH_IMPL_PATH)   # poisthresh_em, glmpois_poisthresh_assign

  stopifnot(length(K) == 1L, is.finite(K), K > 0)

  # Build the shared design matrix exactly as sceptre_assign_pure_R does: augment
  # the covariate data frame with the grna totals, then model.matrix(formula).
  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)
  covariate_matrix <- stats::model.matrix(object = formula, data = cov_df)
  if (nrow(covariate_matrix) != nrow(cov_df)) {
    stop("Rows lost from covariate data frame after applying the formula (likely NAs).")
  }

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, SCEPTRE_IMPL_PATH)
    parallel::clusterCall(cl, source, POISTHRESH_IMPL_PATH)
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- glmpois_poisthresh_assign(
    grna_matrix           = grna_matrix,
    covariate_matrix      = covariate_matrix,
    K                     = K,
    grna_ids              = rownames(grna_matrix),
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    pi_init               = pi_init,
    max_iter              = max_iter,
    tol                   = tol,
    cl                    = cl,
    keep_fits             = keep_fits
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
