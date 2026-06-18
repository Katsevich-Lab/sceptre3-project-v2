# Shared boilerplate for bin/script/{variant}.R files.
#
# A variant file declares: which offset-model fit function to use, an `offset_spec`
# describing it, any extra arguments to forward to sceptre_assign_pure_R(), and
# which packages workers need loaded. Everything else (cov_df augmentation,
# cluster setup, assignment-matrix construction) lives here.
#
# `bin_dir` is defined at the top level of bin/run_script.R, which source()s
# the variant into the global environment; the variant in turn source()s this
# file. We resolve IMPL_PATH here so workers can re-source the implementation.

IMPL_PATH        <- file.path(bin_dir, "script", "lib", "sceptre_assign_pure_R.R")
EM_VARIANTS_PATH <- file.path(bin_dir, "script", "lib", "em_variants.R")

run_variant <- function(grna_matrix, extra_covariates, formula, cpus,
                        offset_model_fit_fn,
                        offset_spec,
                        worker_libraries = "Matrix",
                        ...) {
  for (pkg in worker_libraries) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        paste0("requireNamespace('%s') failed.\n",
               "  .libPaths()  : %s\n",
               "  R.home()     : %s\n",
               "  R_LIBS_USER  : %s\n",
               "  R_LIBS_SITE  : %s\n",
               "  R_LIBS       : %s\n",
               "  installed?   : %s\n",
               "  find.package : %s\n"),
        pkg,
        paste(.libPaths(), collapse = ":"),
        R.home(),
        Sys.getenv("R_LIBS_USER", "<unset>"),
        Sys.getenv("R_LIBS_SITE", "<unset>"),
        Sys.getenv("R_LIBS",      "<unset>"),
        pkg %in% rownames(utils::installed.packages()),
        tryCatch(find.package(pkg), error = function(e) conditionMessage(e))
      ))
    }
  }
  if (!file.exists(IMPL_PATH)) stop("Implementation not found at: ", IMPL_PATH)
  source(IMPL_PATH)
  if (file.exists(EM_VARIANTS_PATH)) source(EM_VARIANTS_PATH)

  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)

  attr(offset_model_fit_fn, "spec") <- offset_spec

  cl <- NULL
  if (cpus > 1L) {
    cl <- parallel::makeCluster(cpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, IMPL_PATH)
    if (file.exists(EM_VARIANTS_PATH)) {
      parallel::clusterCall(cl, source, EM_VARIANTS_PATH)
    }
    for (pkg in worker_libraries) {
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    }
  }

  results <- sceptre_assign_pure_R(
    grna_matrix          = grna_matrix,
    grna_ids             = rownames(grna_matrix),
    covariate_data_frame = cov_df,
    formula_object       = formula,
    cl                   = cl,
    offset_model_fit_fn  = offset_model_fit_fn,
    ...
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
