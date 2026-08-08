# Shared boilerplate for the gated logistic / Gamma-Gaussian variant(s).
#
# Parallels run_variant.R but for the gated model, which has NO offset model and
# no EM-family/offset machinery. The only covariate is grna_n_nonzero (computed
# here from the grna_matrix and shared across guides). A variant file declares
# the model hyperparameters and the worker libraries; everything else (cluster
# setup, assignment-matrix construction) lives here.
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s the
# variant into the global environment; the variant in turn source()s this file.

GATED_IMPL_PATH <- file.path(bin_dir, "script", "lib", "gated_assign.R")

run_gated_variant <- function(grna_matrix, cpus,
                              n_em_rep              = 5L,
                              pi_guess_range        = c(1e-5, 0.1),
                              mu_pert_guess_range   = log(c(10, 5000)),
                              n_fit_cutoff          = 10L,
                              backup_threshold      = 5L,
                              probability_threshold = 0.8,
                              max_iter              = 200L,
                              tol                   = 1e-6,
                              mu_margin             = 1e-3,
                              sigma_floor           = 0.05,
                              seed                  = 4L,
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
  if (!file.exists(GATED_IMPL_PATH)) stop("Implementation not found at: ", GATED_IMPL_PATH)
  source(GATED_IMPL_PATH)

  # Sole covariate: number of nonzero guides per cell (shared across guides).
  grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, GATED_IMPL_PATH)
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- gated_assign(
    grna_matrix           = grna_matrix,
    grna_n_nonzero        = grna_n_nonzero,
    grna_ids              = rownames(grna_matrix),
    n_em_rep              = n_em_rep,
    pi_guess_range        = pi_guess_range,
    mu_pert_guess_range   = mu_pert_guess_range,
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    max_iter              = max_iter,
    tol                   = tol,
    mu_margin             = mu_margin,
    sigma_floor           = sigma_floor,
    seed                  = seed,
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
