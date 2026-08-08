# Shared boilerplate for the poisthresh<K> variant(s): glmpois offset +
# fixed-shift Poisson mixture. Parallels run_poismix_variant.R. The per-cell
# totals grna_n_nonzero / grna_n_umis are computed here and passed to the
# driver, which removes each guide's own contribution (circularity correction).
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s the
# variant into the global environment; the variant in turn source()s this file.

POISTHRESH_IMPL_PATH <- file.path(bin_dir, "script", "lib", "poisthresh_assign.R")

run_poisthresh_variant <- function(grna_matrix, cpus, K,
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
  if (!file.exists(POISTHRESH_IMPL_PATH)) stop("Implementation not found at: ", POISTHRESH_IMPL_PATH)
  source(POISTHRESH_IMPL_PATH)

  stopifnot(length(K) == 1L, is.finite(K), K > 0)

  grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  grna_n_umis    <- Matrix::colSums(grna_matrix)

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, POISTHRESH_IMPL_PATH)
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- poisthresh_assign(
    grna_matrix           = grna_matrix,
    grna_n_nonzero        = grna_n_nonzero,
    grna_n_umis           = grna_n_umis,
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
