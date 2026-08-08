# Shared boilerplate for the glmpoisgrnafix_<mixture> variant(s): grnafix
# (circularity-corrected glmpois) offset + a standard sceptre mixture. Parallels
# run_poisthresh_variant.R, but the mixture/EM is reused verbatim from
# sceptre_assign_pure_R.R rather than the fixed-shift poisthresh EM.
#
# The per-cell totals grna_n_nonzero / grna_n_umis are computed here and passed to
# the driver, which removes each guide's own contribution (circularity correction)
# before fitting the offset. The shared seeded EM starting guesses are generated
# here exactly as sceptre_assign_pure_R() does, so the mixture behaves identically
# to the corresponding glmpois_<mixture> variant -- only the offset differs.
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s the
# variant into the global environment; the variant in turn source()s this file.

IMPL_PATH        <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")
EM_VARIANTS_PATH <- file.path(bin_dir, "script", "lib", "em_variants.R")
GRNAFIX_IMPL_PATH <- file.path(bin_dir, "script", "lib", "grnafix_assign.R")

run_grnafix_variant <- function(grna_matrix, cpus,
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
  if (!file.exists(IMPL_PATH))         stop("Implementation not found at: ", IMPL_PATH)
  if (!file.exists(GRNAFIX_IMPL_PATH)) stop("Implementation not found at: ", GRNAFIX_IMPL_PATH)
  source(IMPL_PATH)                              # assign_one_grna_pure_R, fit_baseline_glm_pure_R
  if (file.exists(EM_VARIANTS_PATH)) source(EM_VARIANTS_PATH)
  source(GRNAFIX_IMPL_PATH)                      # grnafix_assign

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
    parallel::clusterCall(cl, source, GRNAFIX_IMPL_PATH)
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- grnafix_assign(
    grna_matrix            = grna_matrix,
    grna_n_nonzero         = grna_n_nonzero,
    grna_n_umis            = grna_n_umis,
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
