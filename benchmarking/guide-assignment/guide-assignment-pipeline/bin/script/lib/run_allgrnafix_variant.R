# Shared boilerplate for the glmpoisallgrnafix_<mixture> variant(s): allgrnafix
# offset (circularity-corrected grna covars + uncorrected response covars) + a
# standard sceptre mixture. Parallels run_grnafix_variant.R; the only additions
# are pulling the two response covariates out of extra_covariates and passing
# them to the driver.
#
# The per-cell grna totals grna_n_nonzero / grna_n_umis are computed here; the
# driver removes each guide's own contribution before fitting the offset. The
# response covariates are used uncorrected. The shared seeded EM starting guesses
# are generated exactly as sceptre_assign_pure_R() does, so the mixture behaves
# identically to the corresponding glmpois_<mixture> variant -- only the offset
# differs.
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s the
# variant into the global environment; the variant in turn source()s this file.

IMPL_PATH           <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")
EM_VARIANTS_PATH    <- file.path(bin_dir, "script", "lib", "em_variants.R")
ALLGRNAFIX_IMPL_PATH <- file.path(bin_dir, "script", "lib", "allgrnafix_assign.R")

run_allgrnafix_variant <- function(grna_matrix, extra_covariates, cpus,
                                   family                 = "pois",
                                   fix_curr_g_pert_bug    = FALSE,
                                   n_em_rep               = 5L,
                                   pi_guess_range         = c(1e-5, 0.1),
                                   g_pert_guess_range     = log(c(10, 5000)),
                                   n_nonzero_cells_cutoff = 10L,
                                   backup_threshold       = 5L,
                                   probability_threshold  = 0.8,
                                   seed                   = 4L,
                                   worker_libraries       = "Matrix",
                                   keep_fits              = FALSE,
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
  if (!file.exists(IMPL_PATH))            stop("Implementation not found at: ", IMPL_PATH)
  if (!file.exists(ALLGRNAFIX_IMPL_PATH)) stop("Implementation not found at: ", ALLGRNAFIX_IMPL_PATH)
  source(IMPL_PATH)                              # assign_one_grna_pure_R, fit_baseline_glm_pure_R
  if (file.exists(EM_VARIANTS_PATH)) source(EM_VARIANTS_PATH)
  source(ALLGRNAFIX_IMPL_PATH)                   # allgrnafix_assign

  # Response covariates (uncorrected), pulled from the cell covariate frame.
  needed <- c("response_n_umis_full", "response_n_nonzero_full")
  missing <- setdiff(needed, colnames(extra_covariates))
  if (length(missing) > 0L) {
    stop("run_allgrnafix_variant: extra_covariates is missing column(s): ",
         paste(missing, collapse = ", "))
  }
  response_umis <- as.numeric(extra_covariates[["response_n_umis_full"]])
  response_nnz  <- as.numeric(extra_covariates[["response_n_nonzero_full"]])
  if (length(response_umis) != ncol(grna_matrix)) {
    stop("run_allgrnafix_variant: extra_covariates has ", length(response_umis),
         " rows but grna_matrix has ", ncol(grna_matrix), " cells.")
  }

  grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  grna_n_umis    <- Matrix::colSums(grna_matrix)

  # Shared seeded EM starting guesses, identical to sceptre_assign_pure_R().
  # Save/restore the global RNG state so callers' streams aren't perturbed.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  pi_guesses     <- stats::runif(n_em_rep, pi_guess_range[1],     pi_guess_range[2])
  g_pert_guesses <- stats::runif(n_em_rep, g_pert_guess_range[1], g_pert_guess_range[2])
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, IMPL_PATH)
    if (file.exists(EM_VARIANTS_PATH)) parallel::clusterCall(cl, source, EM_VARIANTS_PATH)
    parallel::clusterCall(cl, source, ALLGRNAFIX_IMPL_PATH)
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- allgrnafix_assign(
    grna_matrix            = grna_matrix,
    grna_n_nonzero         = grna_n_nonzero,
    grna_n_umis            = grna_n_umis,
    response_umis          = response_umis,
    response_nnz           = response_nnz,
    pi_guesses             = pi_guesses,
    g_pert_guesses         = g_pert_guesses,
    grna_ids               = rownames(grna_matrix),
    family                 = family,
    fix_curr_g_pert_bug    = fix_curr_g_pert_bug,
    n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
    backup_threshold       = backup_threshold,
    probability_threshold  = probability_threshold,
    keep_fits              = keep_fits,
    cl                     = cl
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
