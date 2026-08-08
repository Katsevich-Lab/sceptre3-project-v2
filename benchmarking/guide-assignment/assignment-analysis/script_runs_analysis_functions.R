# Analysis helpers for the script_* guide-assignment pipeline runs.
#
# These operate on "extras" objects: the contents of the
# assignment_extras_{method}_{dataset}.rds files the pipeline publishes
# (see guide-assignment-pipeline/bin/run_script.R). Each loaded extras object
# has the structure returned by sceptre_assign_pure_R():
#
#   extras$per_guide  -- named list, one element per gRNA (names = grna_ids);
#                        each element holds per-guide EM output, e.g. em_g_pert
#                        (the converged gamma), em_pi, em_log_lik, em_converged,
#                        em_phi_*, prob_quantiles, offset_model_summary, ...
#                        em_g_pert is NA for guides below the n_nonzero cutoff
#                        or where EM didn't converge.
#   extras$run_meta   -- run-level metadata constant across guides, incl. $family.
#
# NB on `em_g_pert` ("gamma"): its meaning depends on the family. It's a
# log-fold-change for nb-shared / nb-separate (and pois w/ fix_curr_g_pert_bug
# = TRUE), but an additive count shift for pois-additive / pois-additive-nonzero
# (and the default pois). Comparing gamma across two different families mixes
# scales; functions below warn when the two runs' families disagree.


# small null-coalescing helper
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---- field extraction -------------------------------------------------------

# Pull a scalar per-guide field into a named numeric vector (names = grna_ids).
# `field` is the name of an element of each extras$per_guide[[id]] (e.g.
# "em_g_pert", "em_pi", "em_log_lik").
extract_per_guide_scalar <- function(extras, field) {
  if (is.null(extras$per_guide)) {
    stop("`extras` has no $per_guide; is this an assignment_extras_*.rds object?")
  }
  vapply(extras$per_guide, function(r) {
    v <- r[[field]]
    if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
  }, numeric(1))
}

# Convenience wrapper for the gamma (perturbation effect) field.
extract_em_g_pert <- function(extras) extract_per_guide_scalar(extras, "em_g_pert")


# ---- gamma vs gamma scatterplot ---------------------------------------------

# Scatterplot of gamma (run x) vs gamma (run y), one point per gRNA shared by
# both runs. Returns a ggplot object.
#
#   extras_x, extras_y : loaded extras objects (run on the same dataset)
#   label_x, label_y   : axis labels (default to each run's family, if present)
#   log_scale          : if TRUE, log10 both axes (sensible for additive deltas,
#                        which are positive counts; gammas can be negative so
#                        the default is FALSE)
#   drop_na            : if TRUE (default), drop guides where either gamma is NA
plot_gamma_vs_gamma <- function(extras_x, extras_y,
                                label_x   = NULL,
                                label_y   = NULL,
                                log_scale = FALSE,
                                drop_na   = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_gamma_vs_gamma requires the ggplot2 package.")
  }

  fam_x <- extras_x$run_meta$family
  fam_y <- extras_y$run_meta$family
  if (!is.null(fam_x) && !is.null(fam_y) && !identical(fam_x, fam_y)) {
    warning("The two runs use different families ('", fam_x, "' vs '", fam_y,
            "'); em_g_pert is not on the same scale across families.")
  }
  if (is.null(label_x)) label_x <- paste0("gamma (", fam_x %||% "x", ")")
  if (is.null(label_y)) label_y <- paste0("gamma (", fam_y %||% "y", ")")

  gx <- extract_em_g_pert(extras_x)
  gy <- extract_em_g_pert(extras_y)

  common <- intersect(names(gx), names(gy))
  if (length(common) == 0L) {
    stop("No gRNA ids shared between the two extras objects -- ",
         "were they run on the same dataset?")
  }
  only_x <- setdiff(names(gx), names(gy))
  only_y <- setdiff(names(gy), names(gx))
  if (length(only_x) || length(only_y)) {
    warning(length(only_x), " guide(s) only in x and ", length(only_y),
            " only in y; plotting the ", length(common), " in common.")
  }

  df <- data.frame(
    grna_id = common,
    gamma_x = gx[common],
    gamma_y = gy[common],
    row.names = NULL
  )
  n_total <- nrow(df)
  if (drop_na) {
    df <- df[is.finite(df$gamma_x) & is.finite(df$gamma_y), , drop = FALSE]
  }
  n_dropped <- n_total - nrow(df)
  if (n_dropped > 0L) {
    message(n_dropped, " of ", n_total,
            " guides dropped (NA/non-finite gamma in at least one run).")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = gamma_x, y = gamma_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(
      x        = label_x,
      y        = label_y,
      title    = "Per-guide EM gamma: run x vs run y",
      subtitle = paste0(nrow(df), " guides")
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw()

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }
  p
}


# ---- offset-GLM coefficients: run x vs run y --------------------------------

# Pull each guide's offset-model GLM coefficients into a named list (names =
# grna_ids); each element is a named numeric coef vector, or NULL for guides
# below the n_nonzero cutoff (offset_model_summary is NULL there).
#
# `coef_names`: optional character vector to relabel coefficients positionally.
# The glm.fit (pois) fits store coefficients under the real covariate names
# ("(Intercept)", "loglib", ...), but the NB / glmrob fits go through a formula
# interface that renames the design columns to X1..Xp, so those store
# uninformative names. Pass coef_names = colnames(build_covariate_matrix(input_dir))
# to restore the real names (and to align a pois run against an NB run).
extract_offset_coefs <- function(extras, coef_names = NULL) {
  if (is.null(extras$per_guide)) {
    stop("`extras` has no $per_guide; is this an assignment_extras_*.rds object?")
  }
  lapply(extras$per_guide, function(r) {
    co <- r$offset_model_summary$coefficients
    if (is.null(co)) return(NULL)
    nm <- names(co)
    co <- as.numeric(co)
    if (!is.null(coef_names)) {
      if (length(co) != length(coef_names)) {
        stop("coef_names has length ", length(coef_names),
             " but a guide has ", length(co), " coefficients.")
      }
      nm <- coef_names
    }
    stats::setNames(co, nm)
  })
}

# Faceted scatterplot of per-guide offset-GLM coefficients: run x vs run y.
# One facet per coefficient that BOTH runs have (matched by name), one point
# per gRNA shared by both runs. Within a facet, x = beta_<coef> from run x and
# y = beta_<coef> from run y. Returns a ggplot object.
#
#   extras_x, extras_y : loaded extras objects (run on the same dataset)
#   label_x, label_y   : method names, used in the axis labels and title
#   coef_names         : optional positional relabeling of coefficients (see
#                        extract_offset_coefs); applied to both runs, so they
#                        must share the coefficient count
#   drop_na            : if TRUE (default), drop (guide, coef) points where
#                        either beta is NA/non-finite
plot_coef_vs_coef <- function(extras_x, extras_y,
                              label_x    = "x",
                              label_y    = "y",
                              coef_names = NULL,
                              drop_na    = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_coef_vs_coef requires the ggplot2 package.")
  }

  cx <- extract_offset_coefs(extras_x, coef_names)
  cy <- extract_offset_coefs(extras_y, coef_names)

  gx <- names(Filter(Negate(is.null), cx))
  gy <- names(Filter(Negate(is.null), cy))
  common_guides <- intersect(gx, gy)
  if (length(common_guides) == 0L) {
    stop("No gRNAs with offset coefficients in both runs (all below the ",
         "n_nonzero cutoff, or non-overlapping datasets?).")
  }

  # long df over (guide, coef) for the coefficients each guide shares
  rows <- lapply(common_guides, function(id) {
    bx <- cx[[id]]; by <- cy[[id]]
    coefs <- intersect(names(bx), names(by))
    if (length(coefs) == 0L) return(NULL)
    data.frame(
      grna_id = id,
      coef    = coefs,
      beta_x  = as.numeric(bx[coefs]),
      beta_y  = as.numeric(by[coefs]),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || nrow(df) == 0L) {
    stop("The two runs share no coefficient names. (NB / glmrob fits store ",
         "coefficients as X1..Xp; pass `coef_names` to relabel them consistently.)")
  }

  n_total <- nrow(df)
  if (drop_na) {
    df <- df[is.finite(df$beta_x) & is.finite(df$beta_y), , drop = FALSE]
  }
  n_dropped <- n_total - nrow(df)
  if (n_dropped > 0L) {
    message(n_dropped, " of ", n_total,
            " (guide, coef) points dropped (NA/non-finite beta in >= 1 run).")
  }

  df$coef <- factor(df$coef, levels = unique(df$coef))

  ggplot2::ggplot(df, ggplot2::aes(x = beta_x, y = beta_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::facet_wrap(~ coef, scales = "free") +
    ggplot2::labs(
      x        = paste0("beta (", label_x, ")"),
      y        = paste0("beta (", label_y, ")"),
      title    = paste0("Per-guide offset-GLM coefficients: ", label_x, " vs ", label_y),
      subtitle = paste0(length(unique(df$grna_id)), " guides; ",
                        nlevels(df$coef), " shared coefficient(s); ",
                        "facet = coefficient")
    ) +
    ggplot2::theme_bw()
}


# ---- reconstructing per-cell offsets and perturbed means --------------------
#
# Neither the per-cell offset (eta = X %*% coef) nor the perturbed-component
# mean (mu_pert) is stored in the extras objects (that needs keep_fits = TRUE).
# But both are reconstructable from what IS stored -- the per-guide offset
# coefficients (offset_model_summary$coefficients) and the converged gamma
# (em_g_pert) -- plus the dataset's covariate matrix X.
#
# Caveat: this reconstruction assumes coef is on the SAME scale/columns as X,
# which holds for glmpois / glmnb (glm.fit on the raw covariate_matrix). It does
# NOT hold for glmrob, which fits on a standardized X; those coefficients won't
# reconstruct eta this way.


# Build the (shared-across-guides) model matrix for a dataset, exactly as the
# pipeline does (run_variant.R augments cov_df with grna_n_nonzero / grna_n_umis,
# then sceptre_assign_pure_R calls model.matrix(formula, cov_df)).
build_covariate_matrix <- function(input_dir) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("build_covariate_matrix requires the Matrix package.")
  }
  extra_covariates <- utils::read.csv(file.path(input_dir, "cell_covariates.csv"))
  grna_matrix      <- readRDS(file.path(input_dir, "grna_matrix.rds"))
  formula          <- stats::as.formula(readRDS(file.path(input_dir, "formula_object.rds")))

  cov_df <- extra_covariates
  cov_df$grna_n_nonzero <- Matrix::colSums(grna_matrix > 0)
  cov_df$grna_n_umis    <- Matrix::colSums(grna_matrix)

  X <- stats::model.matrix(object = formula, data = cov_df)
  if (nrow(X) != nrow(cov_df)) {
    stop("Rows lost building the covariate matrix (NAs in covariates?).")
  }
  X
}


# Per-guide reconstructed offset on the LINK scale: eta^(g)_i = (X %*% coef^(g))_i.
# Returns a named list (grna_id -> numeric length n_cells). Guides with no offset
# fit (below the n_nonzero cutoff) are dropped. NA coefficients (rank-deficient
# fits) are zeroed, mirroring the pipeline's `coef[is.na(coef)] <- 0`.
reconstruct_offset_eta <- function(extras, covariate_matrix) {
  pg <- extras$per_guide
  out <- list()
  for (id in names(pg)) {
    coef <- pg[[id]]$offset_model_summary$coefficients
    if (is.null(coef)) next  # below-cutoff guide: no offset fit
    if (length(coef) != ncol(covariate_matrix)) {
      warning("Guide '", id, "': coef length (", length(coef),
              ") != ncol(X) (", ncol(covariate_matrix), "); skipping.")
      next
    }
    coef[is.na(coef)] <- 0
    out[[id]] <- as.numeric(covariate_matrix %*% coef)
  }
  out
}


# Name of the offset-model fit fn recorded in run_meta (its attr("spec")$name),
# or NA. Tells which link/scale the stored coefficients live on.
offset_fit_name <- function(extras) {
  fn <- extras$run_meta$offset_model_fit_fn
  if (is.null(fn)) return(NA_character_)
  spec <- attr(fn, "spec")
  if (is.null(spec) || is.null(spec$name)) NA_character_ else spec$name
}


# Per-guide offset on the COMMON log-mean scale eta = log(mu0), so offsets are
# comparable across different offset models. Reconstructs X %*% coef per guide
# (NA coefs zeroed), then:
#   * log-link GLM fits (pois / nb): X %*% coef IS log(mu0) -> return as-is.
#   * log1p-OLS fit (fit_baseline_lm_log1p_capped_pure_R): X %*% coef predicts
#     log1p(count), so mu0 = expm1(prediction) and eta = log(mu0) =
#     log(expm1(X %*% coef)), floored at mu_floor (matching the fit's pmax).
# Self-contained (does not call reconstruct_offset_eta) so it stays correct
# regardless of that function's scale convention. glmrob offsets are still NOT
# reconstructable (standardized X); see reconstruct_offset_eta's caveat.
reconstruct_offset_logmu0 <- function(extras, covariate_matrix, mu_floor = 1e-8) {
  nm       <- offset_fit_name(extras)
  is_log1p <- !is.na(nm) && grepl("log1p", nm)
  pg  <- extras$per_guide
  out <- list()
  for (id in names(pg)) {
    coef <- pg[[id]]$offset_model_summary$coefficients
    if (is.null(coef)) next
    if (length(coef) != ncol(covariate_matrix)) {
      warning("Guide '", id, "': coef length (", length(coef),
              ") != ncol(X) (", ncol(covariate_matrix), "); skipping.")
      next
    }
    coef[is.na(coef)] <- 0
    lp <- as.numeric(covariate_matrix %*% coef)
    out[[id]] <- if (is_log1p) log(pmax(expm1(lp), mu_floor)) else lp
  }
  out
}


# ---- offset vs offset scatter -----------------------------------------------

# Per-cell scatter of the offset eta = log(mu0) from run x vs run y, faceted by
# gRNA, for the guides in `grna_ids`. Both runs are put on the common log(mu0)
# scale (the log1p-OLS offset is converted via log(expm1(.))), so points on the
# dashed y = x line mean the two offset models agree. Useful for checking how a
# new offset (e.g. caplmlog1p) tracks a reference offset: a flat shelf at
# ~log(mu_floor) flags floored/degenerate cells, and a systematic sag below the
# line flags a downward level bias.
#
#   extras_x, extras_y : loaded extras objects, run on the SAME dataset
#   grna_ids           : guides to plot (one facet each)
#   input_dir          : dataset dir to rebuild X from (or pass covariate_matrix)
#   covariate_matrix   : precomputed X (skips input_dir if supplied)
#   title_             : custom plot title
#   label_x, label_y   : axis labels (default to each run's offset-fit name)
#   max_cells          : subsample this many cells (shared set across facets)
#   facet_scales       : facet_wrap scales ("fixed" default; "free" to zoom each)
#   seed               : RNG seed for the cell subsample
plot_offset_scatter <- function(extras_x, extras_y,
                                grna_ids,
                                input_dir        = NULL,
                                covariate_matrix = NULL,
                                title_           = NULL,
                                label_x          = NULL,
                                label_y          = NULL,
                                max_cells        = 3000L,
                                facet_scales     = "fixed",
                                seed             = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_offset_scatter requires the ggplot2 package.")
  }
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) stop("Provide either `covariate_matrix` or `input_dir`.")
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(label_x)) label_x <- paste0("offset = log(mu0)  [", offset_fit_name(extras_x) %||% "x", "]")
  if (is.null(label_y)) label_y <- paste0("offset = log(mu0)  [", offset_fit_name(extras_y) %||% "y", "]")

  off_x  <- reconstruct_offset_logmu0(extras_x, covariate_matrix)
  off_y  <- reconstruct_offset_logmu0(extras_y, covariate_matrix)
  common <- intersect(names(off_x), names(off_y))
  grna_ids <- intersect(grna_ids, common)
  if (length(grna_ids) == 0L) {
    stop("None of the requested grna_ids have offset fits in both runs.")
  }

  n_cells  <- nrow(covariate_matrix)
  set.seed(seed)
  cell_idx <- if (n_cells > max_cells) sort(sample.int(n_cells, max_cells)) else seq_len(n_cells)

  rows <- lapply(grna_ids, function(id) {
    data.frame(grna_id  = id,
               offset_x = off_x[[id]][cell_idx],
               offset_y = off_y[[id]][cell_idx],
               row.names = NULL)
  })
  df <- do.call(rbind, rows)
  df$grna_id <- factor(df$grna_id, levels = grna_ids)

  ggplot2::ggplot(df, ggplot2::aes(x = offset_x, y = offset_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.3, size = 0.6) +
    ggplot2::facet_wrap(~ grna_id, scales = facet_scales) +
    ggplot2::labs(
      x        = label_x,
      y        = label_y,
      title    = title_ %||% "Per-cell offset log(mu0): run x vs run y",
      subtitle = paste0(length(grna_ids), " guides; ", length(cell_idx), " cells/facet")
    ) +
    ggplot2::theme_bw()
}


# The perturbed-component mean as a function of the baseline mean mu0 = e^{eta},
# the converged effect gamma, and the family/bug toggle. This is exactly the
# `g_mus_pert1` line inside each EM kernel:
#   pois (bugged, fix=FALSE)            : mu0 + gamma          (additive count bump)
#   pois (fixed,  fix=TRUE)             : mu0 * exp(gamma)     (multiplicative)
#   pois-additive / -additive-nonzero  : mu0 + delta          (additive; gamma IS delta)
#   nb-shared / nb-separate            : mu0 * exp(gamma)     (multiplicative)
perturbed_mean <- function(mu0, gamma, family, fix_curr_g_pert_bug = FALSE) {
  if (family == "pois") {
    if (isTRUE(fix_curr_g_pert_bug)) mu0 * exp(gamma) else mu0 + gamma
  } else if (family %in% c("pois-additive", "pois-additive-nonzero")) {
    mu0 + gamma
  } else {  # nb-shared, nb-separate
    mu0 * exp(gamma)
  }
}


# ---- mu_pert vs mu_pert scatter ---------------------------------------------

# Per-cell scatter of the perturbed-component mean from run x vs run y, faceted
# by gRNA. Reconstructs mu_pert from the stored offset coefficients + gamma (see
# above). Points are coloured by the offset eta = log(mu0), which is the variable
# that drives the bugged-vs-fixed divergence: as eta -> -Inf, the bugged mean
# flattens to ~gamma while the fixed mean collapses to 0.
#
#   extras_x, extras_y : loaded extras objects, run on the SAME dataset
#   input_dir          : dataset dir to rebuild X from (or pass covariate_matrix)
#   covariate_matrix   : precomputed X (skips input_dir if supplied)
#   grna_ids           : which guides to plot (default: first n_guides shared)
#   n_guides           : number of guides to facet over when grna_ids is NULL
#   max_cells          : subsample this many cells for plotting (per facet, the
#                        same sampled cell set is reused across guides)
#   log_scale          : log10 both axes (default TRUE; means span decades)
#   seed               : RNG seed for the cell subsample
plot_mu_pert_scatter <- function(extras_x, extras_y,
                                 input_dir        = NULL,
                                 covariate_matrix = NULL,
                                 grna_ids         = NULL,
                                 n_guides         = 12L,
                                 max_cells        = 3000L,
                                 log_scale        = TRUE,
                                 label_x          = NULL,
                                 label_y          = NULL,
                                 seed             = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_mu_pert_scatter requires the ggplot2 package.")
  }
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide either `covariate_matrix` or `input_dir`.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }

  fam_x <- extras_x$run_meta$family; fix_x <- extras_x$run_meta$fix_curr_g_pert_bug
  fam_y <- extras_y$run_meta$family; fix_y <- extras_y$run_meta$fix_curr_g_pert_bug
  if (is.null(label_x)) {
    label_x <- paste0("mu_pert (", fam_x %||% "x",
                      if (isTRUE(fix_x)) ", fixed)" else if (identical(fam_x, "pois")) ", bugged)" else ")")
  }
  if (is.null(label_y)) {
    label_y <- paste0("mu_pert (", fam_y %||% "y",
                      if (isTRUE(fix_y)) ", fixed)" else if (identical(fam_y, "pois")) ", bugged)" else ")")
  }

  eta_x_all <- reconstruct_offset_eta(extras_x, covariate_matrix)
  eta_y_all <- reconstruct_offset_eta(extras_y, covariate_matrix)
  common    <- intersect(names(eta_x_all), names(eta_y_all))
  if (length(common) == 0L) {
    stop("No guides with offset fits shared between the two runs.")
  }
  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) stop("None of the requested grna_ids have offset fits in both runs.")
  }

  n_cells  <- nrow(covariate_matrix)
  set.seed(seed)
  cell_idx <- if (n_cells > max_cells) sort(sample.int(n_cells, max_cells)) else seq_len(n_cells)

  rows <- lapply(grna_ids, function(id) {
    eta_x   <- eta_x_all[[id]][cell_idx]
    gamma_x <- extras_x$per_guide[[id]]$em_g_pert
    gamma_y <- extras_y$per_guide[[id]]$em_g_pert
    mu0_x   <- exp(eta_x)
    mu0_y   <- exp(eta_y_all[[id]][cell_idx])
    data.frame(
      grna_id   = id,
      eta       = eta_x,                       # offset (same model => same for both here)
      mu_pert_x = perturbed_mean(mu0_x, gamma_x, fam_x, fix_x),
      mu_pert_y = perturbed_mean(mu0_y, gamma_y, fam_y, fix_y),
      row.names = NULL
    )
  })
  df <- do.call(rbind, rows)
  df <- df[is.finite(df$mu_pert_x) & is.finite(df$mu_pert_y), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mu_pert_x, y = mu_pert_y, colour = eta)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.4, size = 0.7) +
    ggplot2::scale_colour_viridis_c(name = "offset eta") +
    ggplot2::facet_wrap(~ grna_id, scales = "free") +
    ggplot2::labs(
      x     = label_x,
      y     = label_y,
      title = "Per-cell perturbed-component mean: run x vs run y",
      subtitle = paste0(length(grna_ids), " guides, ",
                        length(cell_idx), " cells each")
    ) +
    ggplot2::theme_bw()

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }
  p
}


# ---- per-guide 4-group (method x method) y vs offset ------------------------

# Short human label for a run, used in facet strips / legends.
run_label <- function(extras) {
  fam <- extras$run_meta$family %||% "?"
  if (identical(fam, "pois")) {
    if (isTRUE(extras$run_meta$fix_curr_g_pert_bug)) "pois(fixed)" else "pois(bugged)"
  } else {
    fam
  }
}

# Load the gRNA count matrix for a dataset (the y_i source). Kept separate from
# build_covariate_matrix so callers can reuse either independently.
load_grna_matrix <- function(input_dir) {
  readRDS(file.path(input_dir, "grna_matrix.rds"))
}


# Grid of scatterplots: each ROW is a guide, the 4 COLUMNS are the combinations
# of (run x: pert/unpert) x (run y: pert/unpert). Within each panel, observed
# guide count y_i (y-axis) vs offset o_i = eta_i (x-axis), over the cells that
# fall in that 2x2 cell for that guide. Each point's SIZE is log1p(y_i) and its
# COLOUR is the unobserved true perturbation status (true_perts[g, i] == 1) -- so
# in the disagreement columns you can read off who was right.
#
# A cell is "pert" for a run iff its index is in that run's stored assignments
# (per_guide[[g]]$assignments, i.e. posterior T1 >= probability_threshold). The
# offset o_i is reconstructed from run x's stored coefficients (same offset model
# as run y when comparing poisbug vs poisfix).
#
# `grna_matrix` is the gRNA count matrix (y_i source). `true_perts` is OPTIONAL:
# when supplied (a corresponding [gRNA x cell] matrix), points are coloured by
# the latent perturbation indicator (true_perts[g, i] == 1) so you can read off
# who was right in the disagreement columns. When omitted (real data, no ground
# truth), the colour is dropped and the plot still shows where the two methods
# agreed / disagreed.
#
#   extras_x, extras_y     : loaded extras objects, run on the SAME dataset
#   grna_matrix            : sparse gRNA count matrix (y_i source)
#   true_perts             : optional latent perturbation-status matrix (colour);
#                            NULL => no ground truth, no colour
#   input_dir              : dataset dir (loads X for o_i if covariate_matrix NULL)
#   covariate_matrix       : optional precomputed X (skips that reload)
#   grna_ids / n_guides    : which / how many guides (rows); default first n_guides
#   max_points_per_panel   : subsample cap per (guide, combo) panel
#   name_x, name_y         : run labels for the column strips (default run_label())
#   log_y                  : log1p y-axis (counts include zeros); default FALSE
#   seed                   : RNG seed for subsampling
plot_assignment_groups_y_vs_o <- function(extras_x, extras_y,
                                          grna_matrix,
                                          true_perts           = NULL,
                                          input_dir            = NULL,
                                          covariate_matrix     = NULL,
                                          grna_ids             = NULL,
                                          n_guides             = 6L,
                                          max_points_per_panel = 2000L,
                                          name_x               = NULL,
                                          name_y               = NULL,
                                          log_y                = TRUE,
                                          cap_y = NULL,
                                          seed                 = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_assignment_groups_y_vs_o requires the ggplot2 package.")
  }
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_x)) name_x <- run_label(extras_x)
  if (is.null(name_y)) name_y <- run_label(extras_y)
  has_truth <- !is.null(true_perts)

  eta_x_all <- reconstruct_offset_eta(extras_x, covariate_matrix)
  common    <- intersect(names(eta_x_all), names(extras_y$per_guide))
  # keep only guides that actually ran EM in run y too (offset fit present)
  common    <- common[vapply(common, function(id) {
    !is.null(extras_y$per_guide[[id]]$offset_model_summary$coefficients)
  }, logical(1))]
  if (length(common) == 0L) stop("No guides with offset fits shared between the two runs.")

  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) stop("None of the requested grna_ids have offset fits in both runs.")
  }

  # Row facet labels carry each guide's converged gamma for both runs.
  gx <- extract_em_g_pert(extras_x)
  gy <- extract_em_g_pert(extras_y)
  row_labels <- stats::setNames(
    sprintf("%s\n%s g=%.2f\n%s g=%.2f",
            grna_ids, name_x, gx[grna_ids], name_y, gy[grna_ids]),
    grna_ids
  )

  n_cells <- nrow(covariate_matrix)
  lev <- c(
    sprintf("%s: pert\n%s: pert",     name_x, name_y),
    sprintf("%s: pert\n%s: unpert",   name_x, name_y),
    sprintf("%s: unpert\n%s: pert",   name_x, name_y),
    sprintf("%s: unpert\n%s: unpert", name_x, name_y)
  )

  set.seed(seed)
  per_guide_df <- lapply(grna_ids, function(id) {
    o  <- eta_x_all[[id]]                            # offset (link scale)
    y  <- as.numeric(grna_matrix[id, ])              # observed guide count
    tp <- if (has_truth) as.numeric(true_perts[id, ]) == 1 else NA  # latent truth (NA if unknown)
    m1 <- logical(n_cells); m1[extras_x$per_guide[[id]]$assignments] <- TRUE
    m2 <- logical(n_cells); m2[extras_y$per_guide[[id]]$assignments] <- TRUE
    combo <- factor(
      ifelse(m1 & m2, lev[1],
      ifelse(m1 & !m2, lev[2],
      ifelse(!m1 & m2, lev[3], lev[4]))),
      levels = lev
    )
    gdf <- data.frame(grna_id = id, row_label = row_labels[[id]],
                      o = o, y = y, log1p_y = log1p(y),
                      true_pert = tp, combo = combo, row.names = NULL)
    # subsample within each (guide, combo) panel
    do.call(rbind, lapply(split(gdf, gdf$combo, drop = FALSE), function(d) {
      if (nrow(d) > max_points_per_panel) {
        d[sample.int(nrow(d), max_points_per_panel), , drop = FALSE]
      } else d
    }))
  })
  df <- do.call(rbind, per_guide_df)
  df$row_label <- factor(df$row_label, levels = unname(row_labels))

  
  if(!is.null(cap_y)) {
    df$y <- pmin(cap_y, df$y)
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = o, y = y,
                                        # size = log1p_y,
                                        shape = ifelse(y == 0, 1, 2) |> factor())) +
    ggplot2::geom_point(alpha = 0.5) +
    # ggplot2::scale_size_continuous(name = "log1p(y)", range = c(0.3, 3)) +
    ggplot2::facet_grid(row_label ~ combo, scales = "free_y") +
    ggplot2::labs(
      x     = "offset  o_i = eta_i  (link scale)",
      y     = "observed guide count  y_i",
      title = "y vs offset, split by (run x call) x (run y call)",
      subtitle = paste0(length(grna_ids), " guides; up to ",
                        max_points_per_panel, " cells per panel")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
    scale_shape_discrete(solid=F) +
    geom_hline(yintercept = 1)

  if (has_truth) {
    p <- p + ggplot2::aes(colour = true_pert) +
      ggplot2::scale_colour_manual(
        name   = "true perturbed",
        values = c(`FALSE` = "lightblue4", `TRUE` = "firebrick")
      )
  }

  if (log_y) {
    p <- p + ggplot2::scale_y_continuous(trans = "log1p")
  }
  p
}


# ---- mu_pert / mu_nonpert vs true cell mean (cleanser-sim guides) -----------

# For the cleanser simulation guides params_cleanser_mupow={mupow}_grna_*, scatter
# each method's reconstructed component means against the TRUE per-cell mean used
# to generate the data. Each ROW is a guide; the 4 COLUMNS are, per method, the
# non-perturbed mean mu_nonpert = e^{eta} and the perturbed mean mu_pert (via
# perturbed_mean()). Points coloured by true perturbation status
# (true_perts[g, i] == 1). The dashed y = x line is the target: a well-fit
# component mean lands on it for the cells it is meant to describe -- mu_nonpert
# on truly-unperturbed cells, mu_pert on truly-perturbed cells.
#
#   extras_x, extras_y : loaded extras objects, run on the SAME dataset
#   true_pert_means    : length-N (cells) vector of latent perturbed means;
#                        used as the x-axis in the mu_pert columns
#   true_nonpert_means : length-N (cells) vector of latent non-perturbed means;
#                        used as the x-axis in the mu_nonpert columns
#   true_perts         : matrix [gRNA x cell], latent perturbation indicator
#   mupow              : selects guides named params_cleanser_mupow={mupow}_grna_*
#   input_dir          : dataset dir to build X from (or pass covariate_matrix)
#   covariate_matrix   : optional precomputed X
#   grna_ids / n_guides: further restrict / cap guides (rows); default first n_guides
#   max_cells          : subsample this many cells per guide (shared across panels)
#   name_x, name_y     : run labels for the column strips (default run_label())
#   facet_scales       : passed to facet_grid (default "free_y")
#   log_scale          : log10 both axes (default TRUE)
#   seed               : RNG seed for the cell subsample
plot_mu_vs_true_mean <- function(extras_x, extras_y,
                                 true_pert_means, true_nonpert_means,
                                 true_perts, mupow,
                                 input_dir        = NULL,
                                 covariate_matrix = NULL,
                                 grna_ids         = NULL,
                                 n_guides         = 6L,
                                 max_cells        = 2000L,
                                 name_x           = NULL,
                                 name_y           = NULL,
                                 facet_scales     = "free_y",
                                 log_scale        = TRUE,
                                 seed             = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_mu_vs_true_mean requires the ggplot2 package.")
  }
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_x)) name_x <- run_label(extras_x)
  if (is.null(name_y)) name_y <- run_label(extras_y)

  fam_x <- extras_x$run_meta$family; fix_x <- extras_x$run_meta$fix_curr_g_pert_bug
  fam_y <- extras_y$run_meta$family; fix_y <- extras_y$run_meta$fix_curr_g_pert_bug

  eta_x_all <- reconstruct_offset_eta(extras_x, covariate_matrix)
  eta_y_all <- reconstruct_offset_eta(extras_y, covariate_matrix)
  common    <- intersect(names(eta_x_all), names(eta_y_all))
  if (length(common) == 0L) {
    stop("No guides have reconstructable offsets. The covariate_matrix (",
         ncol(covariate_matrix), " cols) doesn't match the offset-model ",
         "coefficients (see the coef-length warnings above). Omit ",
         "`covariate_matrix` and pass `input_dir` so X is rebuilt exactly as ",
         "the pipeline did, or pass `build_covariate_matrix(input_dir)`.")
  }

  # filter to the cleanser-sim guides for this mupow (literal prefix, no regex)
  prefix <- sprintf("params_cleanser_mupow=%s_grna_", mupow)
  common <- common[startsWith(common, prefix)]
  if (length(common) == 0L) {
    stop("No guides match prefix '", prefix, "' with offset fits in both runs.")
  }
  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) stop("None of the requested grna_ids match the mupow prefix.")
  }

  n_cells <- nrow(covariate_matrix)
  if (length(true_pert_means) != n_cells || length(true_nonpert_means) != n_cells) {
    stop("true_pert_means / true_nonpert_means must each be length ", n_cells,
         " (one per cell); got ", length(true_pert_means), " and ",
         length(true_nonpert_means), ".")
  }
  true_pert_means    <- as.numeric(true_pert_means)
  true_nonpert_means <- as.numeric(true_nonpert_means)

  panel_lev <- c(
    sprintf("%s: mu_nonpert", name_x),
    sprintf("%s: mu_pert",    name_x),
    sprintf("%s: mu_nonpert", name_y),
    sprintf("%s: mu_pert",    name_y)
  )

  set.seed(seed)
  per_guide_df <- lapply(grna_ids, function(id) {
    cell_idx <- if (n_cells > max_cells) sort(sample.int(n_cells, max_cells)) else seq_len(n_cells)
    tm_p   <- true_pert_means[cell_idx]      # x-axis for the mu_pert columns
    tm_np  <- true_nonpert_means[cell_idx]   # x-axis for the mu_nonpert columns
    tp     <- as.numeric(true_perts[id, ])[cell_idx] == 1
    mu0_x  <- exp(eta_x_all[[id]][cell_idx])
    mu0_y  <- exp(eta_y_all[[id]][cell_idx])
    mp_x   <- perturbed_mean(mu0_x, extras_x$per_guide[[id]]$em_g_pert, fam_x, fix_x)
    mp_y   <- perturbed_mean(mu0_y, extras_y$per_guide[[id]]$em_g_pert, fam_y, fix_y)
    meta   <- data.frame(grna_id = id, true_pert = tp, row.names = NULL)
    rbind(
      cbind(meta, true_mean = tm_np, mu = mu0_x, panel = panel_lev[1]),
      cbind(meta, true_mean = tm_p,  mu = mp_x,  panel = panel_lev[2]),
      cbind(meta, true_mean = tm_np, mu = mu0_y, panel = panel_lev[3]),
      cbind(meta, true_mean = tm_p,  mu = mp_y,  panel = panel_lev[4])
    )
  })
  df <- do.call(rbind, per_guide_df)
  df$panel   <- factor(df$panel, levels = panel_lev)
  df$grna_id <- factor(df$grna_id, levels = grna_ids)
  df <- df[is.finite(df$true_mean) & is.finite(df$mu), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = true_mean, y = mu, colour = true_pert)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.4, size = 0.7) +
    ggplot2::scale_colour_manual(
      name   = "true perturbed",
      values = c(`FALSE` = "lightblue4", `TRUE` = "firebrick")
    ) +
    ggplot2::facet_grid(grna_id ~ panel, scales = facet_scales) +
    ggplot2::labs(
      x        = "true per-cell mean",
      y        = "reconstructed component mean",
      title    = sprintf("Component means vs true mean (mupow = %s)", mupow),
      subtitle = paste0(length(grna_ids), " guides; up to ", max_cells, " cells each")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }
  p
}


# ---- Jaccard agreement vs guide tail mean -----------------------------------

# Per-guide Jaccard agreement between the assignment sets of two runs, plotted
# against each guide's "tail mean": mean(grna_matrix[g, grna_matrix[g, ] >= thresh]).
# One point per guide shared by both runs.
#
# Jaccard(A1, A2) = |A1 intersect A2| / |A1 union A2|, where A1 / A2 are the cell
# indices each run assigned to guide g (per_guide[[g]]$assignments). When both
# runs assign zero cells the union is empty and Jaccard is undefined -> set to
# `empty_union` (default NA, dropped).
#
#   extras_x, extras_y : loaded extras objects, run on the SAME dataset
#   grna_matrix        : gRNA count matrix (tail-mean source; rownames = grna_ids)
#   thresh             : tail threshold; tail mean averages counts >= thresh
#   grna_ids           : optional restriction (default all guides in both runs)
#   empty_union        : Jaccard value when both runs assign nothing (default NA)
#   log_x              : log10 x-axis (default TRUE; tail means are positive counts)
#   name_x, name_y     : run labels for the title (default run_label())
plot_jaccard_vs_tail_mean <- function(extras_x, extras_y, grna_matrix, thresh,
                                      grna_ids    = NULL,
                                      empty_union = NA_real_,
                                      log_x       = TRUE,
                                      name_x      = NULL,
                                      name_y      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_jaccard_vs_tail_mean requires the ggplot2 package.")
  }
  if (is.null(name_x)) name_x <- run_label(extras_x)
  if (is.null(name_y)) name_y <- run_label(extras_y)

  common <- intersect(names(extras_x$per_guide), names(extras_y$per_guide))
  common <- intersect(common, rownames(grna_matrix))
  if (length(common) == 0L) stop("No guides shared between the two runs and grna_matrix.")
  if (!is.null(grna_ids)) {
    common <- intersect(grna_ids, common)
    if (length(common) == 0L) stop("None of the requested grna_ids are shared.")
  }

  rows <- lapply(common, function(id) {
    a1  <- extras_x$per_guide[[id]]$assignments
    a2  <- extras_y$per_guide[[id]]$assignments
    uni <- length(union(a1, a2))
    jac <- if (uni == 0L) empty_union else length(intersect(a1, a2)) / uni
    g    <- as.numeric(grna_matrix[id, ])
    tail <- g[g >= thresh]
    tail_mean <- if (length(tail) == 0L) NA_real_ else mean(tail)
    data.frame(grna_id = id, tail_mean = tail_mean, jaccard = jac,
               n_union = uni, row.names = NULL)
  })
  df <- do.call(rbind, rows)

  n_total <- nrow(df)
  df <- df[is.finite(df$tail_mean) & is.finite(df$jaccard), , drop = FALSE]
  n_dropped <- n_total - nrow(df)
  if (n_dropped > 0L) {
    message(n_dropped, " of ", n_total,
            " guides dropped (no counts >= thresh, or empty assignment union).")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = tail_mean, y = jaccard)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x        = sprintf("tail mean  (mean of y_i >= %s)", format(thresh)),
      y        = "Jaccard agreement of assignments",
      title    = sprintf("Assignment agreement vs guide tail mean: %s vs %s",
                         name_x, name_y),
      subtitle = paste0(nrow(df), " guides")
    ) +
    ggplot2::theme_bw()

  if (log_x) p <- p + ggplot2::scale_x_log10()
  p
}


# ---- one-sided disagreement histograms --------------------------------------

# For each guide, count the cells where the two runs disagree, in each direction:
#   (1) run x says pert, run y says unpert   = |A1 \ A2|
#   (2) run y says pert, run x says unpert   = |A2 \ A1|
# where A1 / A2 are the stored assignment sets (per_guide[[g]]$assignments).
# Returns a 2-panel plot, one histogram per direction over all shared guides.
#
#   extras_x, extras_y : loaded extras objects, run on the SAME dataset
#   grna_ids           : optional restriction (default all guides in both runs)
#   bins               : histogram bins (default 30)
#   log_x              : log1p x-axis (counts include 0); default FALSE
#   facet_scales       : passed to facet_wrap (default "fixed" so the two panels
#                        share an x-axis -> asymmetry is read directly)
#   name_x, name_y     : run labels for the panel strips (default run_label())
plot_disagreement_histograms <- function(extras_x, extras_y,
                                         grna_ids     = NULL,
                                         bins         = 30L,
                                         log_x        = FALSE,
                                         facet_scales = "fixed",
                                         name_x       = NULL,
                                         name_y       = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_disagreement_histograms requires the ggplot2 package.")
  }
  if (is.null(name_x)) name_x <- run_label(extras_x)
  if (is.null(name_y)) name_y <- run_label(extras_y)

  common <- intersect(names(extras_x$per_guide), names(extras_y$per_guide))
  if (length(common) == 0L) stop("No guides shared between the two runs.")
  if (!is.null(grna_ids)) {
    common <- intersect(grna_ids, common)
    if (length(common) == 0L) stop("None of the requested grna_ids are shared.")
  }

  dir_lev <- c(
    sprintf("%s pert, %s unpert", name_x, name_y),
    sprintf("%s pert, %s unpert", name_y, name_x)
  )

  rows <- lapply(common, function(id) {
    a1 <- extras_x$per_guide[[id]]$assignments
    a2 <- extras_y$per_guide[[id]]$assignments
    data.frame(
      grna_id   = id,
      count     = c(length(setdiff(a1, a2)), length(setdiff(a2, a1))),
      direction = factor(dir_lev, levels = dir_lev),
      row.names = NULL
    )
  })
  df <- do.call(rbind, rows)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = count)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::facet_wrap(~ direction, scales = facet_scales) +
    ggplot2::labs(
      x        = "number of cells (one-sided disagreement)",
      y        = "number of guides",
      title    = sprintf("One-sided assignment disagreement: %s vs %s",
                         name_x, name_y),
      subtitle = paste0(length(common), " guides")
    ) +
    ggplot2::theme_bw()

  if (log_x) p <- p + ggplot2::scale_x_continuous(trans = "log1p")
  p
}


# ---- per-guide 4-group log1p(y) vs offset -----------------------------------

# Same 4-bin split as plot_assignment_groups_y_vs_o (each ROW a guide, the 4
# COLUMNS the combinations of run x: pert/unpert x run y: pert/unpert), but the
# y-axis is log1p(y_i) plotted directly and there is no ground-truth colouring.
# Within each panel: log1p(observed guide count) vs offset o_i = eta_i, over the
# cells that fall in that 2x2 bin for that guide. Row labels carry each guide's
# converged gamma for both runs.
#
#   extras_x, extras_y     : loaded extras objects, run on the SAME dataset
#   grna_matrix            : sparse gRNA count matrix (y_i source)
#   input_dir              : dataset dir (loads X for o_i if covariate_matrix NULL)
#   covariate_matrix       : optional precomputed X (skips that reload)
#   grna_ids / n_guides    : which / how many guides (rows); default first n_guides
#   max_points_per_panel   : subsample cap per (guide, combo) panel
#   name_x, name_y         : run labels for the column strips (default run_label())
#   facet_scales           : passed to facet_grid (default "free_y")
#   drop                   : facet_grid drop; FALSE (default) keeps all 4 columns
#                            even when a combo bin is empty
#   seed                   : RNG seed for subsampling
plot_groups_log1py_vs_o <- function(extras_x, extras_y,
                                    grna_matrix,
                                    input_dir            = NULL,
                                    covariate_matrix     = NULL,
                                    grna_ids             = NULL,
                                    n_guides             = 6L,
                                    max_points_per_panel = 2000L,
                                    name_x               = NULL,
                                    name_y               = NULL,
                                    facet_scales         = "free_y",
                                    drop                 = FALSE,
                                    title_ = "",
                                    seed                 = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_groups_log1py_vs_o requires the ggplot2 package.")
  }
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_x)) name_x <- run_label(extras_x)
  if (is.null(name_y)) name_y <- run_label(extras_y)

  eta_x_all <- reconstruct_offset_eta(extras_x, covariate_matrix)
  common    <- intersect(names(eta_x_all), names(extras_y$per_guide))
  common    <- common[vapply(common, function(id) {
    !is.null(extras_y$per_guide[[id]]$offset_model_summary$coefficients)
  }, logical(1))]
  if (length(common) == 0L) stop("No guides with offset fits shared between the two runs.")

  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) stop("None of the requested grna_ids have offset fits in both runs.")
  }

  # Row facet labels carry each guide's converged gamma for both runs.
  gx <- extract_em_g_pert(extras_x)
  gy <- extract_em_g_pert(extras_y)
  row_labels <- stats::setNames(
    sprintf("%s\n%s g=%.2f\n%s g=%.2f",
            grna_ids, name_x, gx[grna_ids], name_y, gy[grna_ids]),
    grna_ids
  )

  n_cells <- nrow(covariate_matrix)
  lev <- c(
    sprintf("%s: pert\n%s: pert",     name_x, name_y),
    sprintf("%s: pert\n%s: unpert",   name_x, name_y),
    sprintf("%s: unpert\n%s: pert",   name_x, name_y),
    sprintf("%s: unpert\n%s: unpert", name_x, name_y)
  )

  set.seed(seed)
  per_guide_df <- lapply(grna_ids, function(id) {
    o  <- eta_x_all[[id]]                            # offset (link scale)
    y  <- as.numeric(grna_matrix[id, ])              # observed guide count
    m1 <- logical(n_cells); m1[extras_x$per_guide[[id]]$assignments] <- TRUE
    m2 <- logical(n_cells); m2[extras_y$per_guide[[id]]$assignments] <- TRUE
    combo <- factor(
      ifelse(m1 & m2, lev[1],
      ifelse(m1 & !m2, lev[2],
      ifelse(!m1 & m2, lev[3], lev[4]))),
      levels = lev
    )
    gdf <- data.frame(grna_id = id, row_label = row_labels[[id]],
                      o = o, log1p_y = log1p(y), combo = combo, row.names = NULL)
    do.call(rbind, lapply(split(gdf, gdf$combo, drop = FALSE), function(d) {
      if (nrow(d) > max_points_per_panel) {
        d[sample.int(nrow(d), max_points_per_panel), , drop = FALSE]
      } else d
    }))
  })
  df <- do.call(rbind, per_guide_df)
  df$row_label <- factor(df$row_label, levels = unname(row_labels))

  ggplot2::ggplot(df, ggplot2::aes(x = o, y = log1p_y)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.7) +
    ggplot2::facet_grid(row_label ~ combo, scales = facet_scales, drop = drop) +
    ggplot2::labs(
      x        = "offset  o_i = eta_i  (link scale)",
      y        = "log1p(observed guide count)  log1p(y_i)",
      title    = title_,
      subtitle = paste0(length(grna_ids), " guides; up to ",
                        max_points_per_panel, " cells per panel")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}


# ---- per-cell zero-truncated upper-tail non-perturbed score t_i (Poisson) ----
#
# For a guide g, score each cell i by how far into the upper tail of the
# NON-PERTURBED Poisson law its observed count y_i sits, conditioned on being a
# detected (>= 1) count:
#
#   t_i := -log P_nonpert(Y_i >= max(y_i, 1) | Y_i >= 1)
#        = -log [ P(Pois(lambda_i) >= max(y_i, 1)) / (1 - dpois(0, lambda_i)) ],
#
# lambda_i = mu_{0i} = exp(o_i), the per-cell non-pert Poisson mean; o_i is the
# offset (GLM linear predictor) reconstructed from the stored offset coefficients
# via reconstruct_offset_logmu0 (correct for the log-link GLM and log1p-OLS
# offsets alike). Dividing by 1 - dpois(0, .) = 1 - e^{-lambda_i} renormalizes
# onto the zero-truncated law P(. | Y >= 1); max(y_i, 1) maps zero-count cells
# onto the y = 1 boundary, where t_i = 0 (as it is for any y_i in {0, 1}). t_i is
# monotone increasing in y_i for fixed lambda_i; larger t_i = deeper in the tail
# = stronger univariate evidence of perturbation. t_i >= 0 always.
#
# In R, P(Y >= k) = ppois(k - 1, lambda, lower.tail = FALSE).
#
# Returns a named list (grna_id -> numeric length n_cells), one t_i per cell, in
# covariate_matrix / grna_matrix column (cell) order -- mirroring
# reconstruct_offset_eta. Guides with no reconstructable offset (below the
# n_nonzero cutoff) are dropped. glmrob offsets (fit on standardized X) are not
# reconstructable and error out rather than return silently-wrong scores.
compute_nonpert_score_pois <- function(extras, grna_matrix,
                                       input_dir        = NULL,
                                       covariate_matrix = NULL,
                                       grna_ids         = NULL) {
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    covariate_matrix <- build_covariate_matrix(input_dir)
  }

  # glmrob fits on a standardized X, so X %*% coef does NOT reconstruct o_i.
  nm <- offset_fit_name(extras)
  if (!is.na(nm) && grepl("glmrob", nm)) {
    stop("Offset model '", nm, "' was fit on standardized covariates; o_i is not ",
         "reconstructable from the stored coefficients, so t_i is undefined.")
  }

  # o_i = log(mu_{0i}) on the common log-mean scale; lambda_i = exp(o_i).
  log_mu0 <- reconstruct_offset_logmu0(extras, covariate_matrix)

  ids <- names(log_mu0)
  if (!is.null(grna_ids)) {
    ids <- intersect(grna_ids, ids)
    if (length(ids) == 0L) stop("None of the requested grna_ids have a reconstructable offset.")
  }

  # log(1 - e^{-lambda}) = log(1 - dpois(0, lambda)), stable for all lambda > 0
  # (Machler's log1mexp split at log 2).
  log1m_p0 <- function(lambda) {
    ifelse(lambda <= log(2), log(-expm1(-lambda)), log1p(-exp(-lambda)))
  }

  out <- list()
  for (id in ids) {
    lambda <- exp(log_mu0[[id]])
    y      <- as.numeric(grna_matrix[id, ])
    if (length(y) != length(lambda)) {
      warning("Guide '", id, "': grna_matrix has ", length(y),
              " cells but offset has ", length(lambda), "; skipping.")
      next
    }
    k        <- pmax(y, 1)
    log_surv <- stats::ppois(k - 1, lambda, lower.tail = FALSE, log.p = TRUE)  # log P(Y >= k)
    out[[id]] <- -(log_surv - log1m_p0(lambda))                               # t_i
  }
  out
}


# ---- per-cell non-pert tail score t_i vs observed count y_i -----------------
#
# For a selection of guides, scatter the per-cell non-perturbed upper-tail score
#   t_i = -log P_nonpert(Y_i >= max(y_i, 1) | Y_i >= 1)
# (see compute_nonpert_score_pois) on the y-axis against the observed guide count
# y_i on the x-axis, one facet per guide. Each point is COLOURED by the cell's
# TRUE perturbation status (true_perts[g, i] == 1) and SHAPED by the run's
# ASSIGNED status (whether the cell is in per_guide[[g]]$assignments). t_i is the
# single-run analog of the offset scatter, so it works on real data too: omit
# `true_perts` and the colour aesthetic is dropped, leaving shape = assigned.
#
#   extras                 : a loaded extras object
#   grna_matrix            : sparse gRNA count matrix (y_i source; rownames = grna_ids)
#   true_perts             : optional latent perturbation matrix [gRNA x cell];
#                            NULL (default) => no ground truth, no colour aesthetic
#   input_dir              : dataset dir (builds X for o_i if covariate_matrix NULL)
#   covariate_matrix       : optional precomputed X (skips that reload)
#   grna_ids / n_guides    : which / how many guides (facets); default first n_guides
#   max_points_per_panel   : subsample cap per guide facet
#   name_                  : run label for the title (default run_label())
#   cap_y                  : optional cap on y_i (counts pmin()'d before plotting)
#   show_true              : restrict cells by TRUE status: "both" (default) keeps
#                            all; "np" only truly non-perturbed; "p" only truly
#                            perturbed. "np"/"p" require `true_perts`.
#   show_est               : restrict cells by ASSIGNED status (the run's own
#                            assignments), same values as `show_true`: "both"
#                            (default), "np" only assigned non-pert, "p" only
#                            assigned pert. Always available (no ground truth needed).
#   show_pos               : if TRUE, plot only cells with y_i > 0 (drop the
#                            zero-count cells); default FALSE. Applied on top of
#                            `show_true` / `show_est`, on the uncapped count.
#   log_x                  : log1p x-axis (counts include zeros); default TRUE
#   facet_scales           : passed to facet_wrap (default "free")
#   seed                   : RNG seed for subsampling
plot_t_score_vs_y <- function(extras,
                              grna_matrix,
                              true_perts           = NULL,
                              input_dir            = NULL,
                              covariate_matrix     = NULL,
                              grna_ids             = NULL,
                              n_guides             = 6L,
                              max_points_per_panel = 2000L,
                              name_                = NULL,
                              cap_y                = NULL,
                              show_true            = c("both", "np", "p"),
                              show_est             = c("both", "np", "p"),
                              show_pos             = FALSE,
                              log_x                = TRUE,
                              facet_scales         = "free",
                              seed                 = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_t_score_vs_y requires the ggplot2 package.")
  }
  show_true <- match.arg(show_true)
  show_est  <- match.arg(show_est)
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_)) name_ <- run_label(extras)
  has_truth <- !is.null(true_perts)
  if (show_true != "both" && !has_truth) {
    stop("show_true = '", show_true, "' filters by true perturbation status, ",
         "but `true_perts` was not supplied.")
  }

  # per-cell t_i = -log P_nonpert(Y >= max(y,1) | Y >= 1), one vector per guide.
  t_list <- compute_nonpert_score_pois(extras, grna_matrix,
                                       covariate_matrix = covariate_matrix)
  common <- intersect(names(t_list), rownames(grna_matrix))
  if (length(common) == 0L) {
    stop("No guides with reconstructable scores also present in grna_matrix.")
  }
  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) {
      stop("None of the requested grna_ids have reconstructable scores / are in grna_matrix.")
    }
  }

  # Facet labels carry each guide's converged gamma.
  g_pert <- extract_em_g_pert(extras)
  facet_labels <- stats::setNames(
    sprintf("%s\ng=%.2f", grna_ids, g_pert[grna_ids]),
    grna_ids
  )

  n_cells <- nrow(covariate_matrix)
  set.seed(seed)
  per_guide_df <- lapply(grna_ids, function(id) {
    t   <- t_list[[id]]                             # transformed score
    y   <- as.numeric(grna_matrix[id, ])           # observed guide count
    est <- logical(n_cells); est[extras$per_guide[[id]]$assignments] <- TRUE
    tp  <- if (has_truth) as.numeric(true_perts[id, ]) == 1 else NA
    gdf <- data.frame(grna_id = id, facet_label = facet_labels[[id]],
                      t = t, y = y, est_pert = est, true_pert = tp,
                      row.names = NULL)
    if (show_true == "np") {
      gdf <- gdf[!gdf$true_pert, , drop = FALSE]
    } else if (show_true == "p") {
      gdf <- gdf[gdf$true_pert, , drop = FALSE]
    }
    if (show_est == "np") {
      gdf <- gdf[!gdf$est_pert, , drop = FALSE]
    } else if (show_est == "p") {
      gdf <- gdf[gdf$est_pert, , drop = FALSE]
    }
    if (show_pos) {
      gdf <- gdf[gdf$y > 0, , drop = FALSE]
    }
    if (nrow(gdf) > max_points_per_panel) {
      gdf <- gdf[sample.int(nrow(gdf), max_points_per_panel), , drop = FALSE]
    }
    gdf
  })
  df <- do.call(rbind, per_guide_df)
  df$facet_label <- factor(df$facet_label, levels = unname(facet_labels))
  df$est_pert    <- factor(df$est_pert, levels = c(FALSE, TRUE))
  if (has_truth) df$true_pert <- factor(df$true_pert, levels = c(FALSE, TRUE))

  if (!is.null(cap_y)) df$y <- pmin(cap_y, df$y)

  # colour = true perturbation status (when available); shape = assigned status.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = y, y = t)) +
    ggplot2::facet_wrap(~ facet_label, scales = facet_scales) +
    ggplot2::labs(
      x     = if (is.null(cap_y)) "observed guide count  y_i"
              else sprintf("observed guide count  y_i  (capped at %s)", format(cap_y)),
      y     = "t_i  =  -log P_nonpert(Y_i >= y_i | Y_i >= 1)",
      title = sprintf("non-pert tail score t_i vs y_i (%s)", name_),
      subtitle = paste0(length(grna_ids), " guides; up to ",
                        max_points_per_panel, " cells per panel")
    ) +
    ggplot2::theme_bw()

  if (has_truth) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = true_pert, shape = est_pert),
                          alpha = 0.5) +
      ggplot2::scale_colour_manual(
        name   = "true perturbed",
        values = c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4"),
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name   = "assigned perturbed",
        values = c(`FALSE` = 1, `TRUE` = 2),
        drop   = FALSE
      )
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(shape = est_pert), alpha = 0.5) +
      ggplot2::scale_shape_manual(
        name   = "assigned perturbed",
        values = c(`FALSE` = 1, `TRUE` = 2),
        drop   = FALSE
      )
  }

  # x axis: observed counts, anchored at 0. log1p breaks (1-3-10 style) when log_x.
  if (log_x) {
    cand <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000,
              10000, 30000, 1e5, 3e5, 1e6)
    log1p_breaks <- function(limits) {
      b <- cand[cand <= max(limits, na.rm = TRUE)]
      if (length(b) < 2L) c(0, max(limits, na.rm = TRUE)) else b
    }
    p <- p + ggplot2::scale_x_continuous(
      trans  = "log1p",
      breaks = log1p_breaks,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    )
  }
  p <- p + ggplot2::expand_limits(x = 0, y = 0)
  p
}


# ---- per-guide (t_i vs y_i scatter) + (t_i histogram) panels ----------------
#
# Same per-cell non-pert tail score t_i and same cell-selection flags as
# plot_t_score_vs_y, but each guide gets a two-plot panel laid out side by side:
#   (1) scatter of t_i (y-axis) vs observed count y_i (x-axis), coloured by true
#       perturbation status and shaped by assigned status (as in plot_t_score_vs_y);
#   (2) histogram of t_i over the same selected cells.
# Guides are stacked vertically (one row per guide) and the figure carries a
# single shared legend at the bottom. Returns a cowplot object.
#
# `t_max` caps t_i (pmin) in BOTH the scatter and the histogram, since the tail
# score can run very large; capped cells pile up at t_max. The scatter is
# subsampled to `max_points_per_panel`, but the histogram uses ALL selected cells
# (so its shape is not thinned).
#
#   extras                 : a loaded extras object
#   grna_matrix            : sparse gRNA count matrix (y_i source; rownames = grna_ids)
#   true_perts             : optional latent perturbation matrix [gRNA x cell];
#                            NULL (default) => no ground truth, no colour aesthetic
#   input_dir              : dataset dir (builds X for o_i if covariate_matrix NULL)
#   covariate_matrix       : optional precomputed X (skips that reload)
#   grna_ids / n_guides    : which / how many guides (rows); default first n_guides
#   max_points_per_panel   : subsample cap for the SCATTER (histogram uses all)
#   name_                  : run label for the title (default run_label())
#   cap_y                  : optional cap on y_i in the scatter (pmin before plotting)
#   t_max                  : cap on t_i in BOTH plots (default 500)
#   bins                   : histogram bin count (default 50)
#   show_true / show_est / show_pos : cell-selection flags, identical to
#                            plot_t_score_vs_y
#   log_x                  : log1p x-axis for the scatter (counts); default TRUE
#   seed                   : RNG seed for the scatter subsample
plot_t_score_panels <- function(extras,
                                grna_matrix,
                                true_perts           = NULL,
                                input_dir            = NULL,
                                covariate_matrix     = NULL,
                                grna_ids             = NULL,
                                n_guides             = 6L,
                                max_points_per_panel = 2000L,
                                name_                = NULL,
                                cap_y                = NULL,
                                t_max                = 500,
                                bins                 = 50L,
                                show_true            = c("both", "np", "p"),
                                show_est             = c("both", "np", "p"),
                                show_pos             = FALSE,
                                log_x                = TRUE,
                                seed                 = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_t_score_panels requires the ggplot2 package.")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("plot_t_score_panels requires the cowplot package.")
  }
  show_true <- match.arg(show_true)
  show_est  <- match.arg(show_est)
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_)) name_ <- run_label(extras)
  has_truth <- !is.null(true_perts)
  if (show_true != "both" && !has_truth) {
    stop("show_true = '", show_true, "' filters by true perturbation status, ",
         "but `true_perts` was not supplied.")
  }

  t_list <- compute_nonpert_score_pois(extras, grna_matrix,
                                       covariate_matrix = covariate_matrix)
  common <- intersect(names(t_list), rownames(grna_matrix))
  if (length(common) == 0L) {
    stop("No guides with reconstructable scores also present in grna_matrix.")
  }
  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) {
      stop("None of the requested grna_ids have reconstructable scores / are in grna_matrix.")
    }
  }

  g_pert  <- extract_em_g_pert(extras)
  n_cells <- nrow(covariate_matrix)

  # log1p x breaks (1-3-10 style) for the count axis of the scatter.
  cand <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000,
            10000, 30000, 1e5, 3e5, 1e6)
  log1p_breaks <- function(limits) {
    b <- cand[cand <= max(limits, na.rm = TRUE)]
    if (length(b) < 2L) c(0, max(limits, na.rm = TRUE)) else b
  }

  # Selected cells for a guide (flags applied; t_i capped at t_max). No subsample.
  build_gdf <- function(id) {
    t   <- pmin(t_max, t_list[[id]])
    y   <- as.numeric(grna_matrix[id, ])
    est <- logical(n_cells); est[extras$per_guide[[id]]$assignments] <- TRUE
    tp  <- if (has_truth) as.numeric(true_perts[id, ]) == 1 else NA
    gdf <- data.frame(t = t, y = y, est_pert = est, true_pert = tp,
                      row.names = NULL)
    if (show_true == "np") {
      gdf <- gdf[!gdf$true_pert, , drop = FALSE]
    } else if (show_true == "p") {
      gdf <- gdf[gdf$true_pert, , drop = FALSE]
    }
    if (show_est == "np") {
      gdf <- gdf[!gdf$est_pert, , drop = FALSE]
    } else if (show_est == "p") {
      gdf <- gdf[gdf$est_pert, , drop = FALSE]
    }
    if (show_pos) {
      gdf <- gdf[gdf$y > 0, , drop = FALSE]
    }
    gdf$est_pert <- factor(gdf$est_pert, levels = c(FALSE, TRUE))
    if (has_truth) gdf$true_pert <- factor(gdf$true_pert, levels = c(FALSE, TRUE))
    gdf
  }

  make_scatter <- function(id, gdf) {
    d <- gdf
    if (nrow(d) > max_points_per_panel) {
      d <- d[sample.int(nrow(d), max_points_per_panel), , drop = FALSE]
    }
    if (!is.null(cap_y)) d$y <- pmin(cap_y, d$y)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = y, y = t)) +
      ggplot2::labs(
        x     = if (is.null(cap_y)) "observed guide count  y_i"
                else sprintf("y_i  (capped at %s)", format(cap_y)),
        y     = sprintf("t_i  (capped at %s)", format(t_max)),
        title = sprintf("%s   g=%.2f", id, g_pert[[id]])
      ) +
      ggplot2::theme_bw()
    if (has_truth) {
      p <- p +
        ggplot2::geom_point(ggplot2::aes(colour = true_pert, shape = est_pert),
                            alpha = 0.5) +
        ggplot2::scale_colour_manual(
          name = "true perturbed",
          values = c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4"), drop = FALSE) +
        ggplot2::scale_shape_manual(
          name = "assigned perturbed",
          values = c(`FALSE` = 1, `TRUE` = 2), drop = FALSE)
    } else {
      p <- p +
        ggplot2::geom_point(ggplot2::aes(shape = est_pert), alpha = 0.5) +
        ggplot2::scale_shape_manual(
          name = "assigned perturbed",
          values = c(`FALSE` = 1, `TRUE` = 2), drop = FALSE)
    }
    if (log_x) {
      p <- p + ggplot2::scale_x_continuous(
        trans = "log1p", breaks = log1p_breaks,
        expand = ggplot2::expansion(mult = c(0, 0.02)))
    }
    p + ggplot2::expand_limits(x = 0, y = 0)
  }

  make_hist <- function(id, gdf) {
    ggplot2::ggplot(gdf, ggplot2::aes(x = t)) +
      ggplot2::geom_histogram(bins = bins) +
      ggplot2::labs(
        x     = sprintf("t_i  (capped at %s)", format(t_max)),
        y     = "number of cells",
        title = sprintf("n = %d cells", nrow(gdf))
      ) +
      ggplot2::expand_limits(x = 0) +
      ggplot2::theme_bw()
  }

  set.seed(seed)
  rows <- lapply(grna_ids, function(id) {
    gdf <- build_gdf(id)
    cowplot::plot_grid(
      make_scatter(id, gdf) + ggplot2::theme(legend.position = "none"),
      make_hist(id, gdf),
      ncol = 2, rel_widths = c(2, 1)
    )
  })
  body <- cowplot::plot_grid(plotlist = rows, ncol = 1)

  # one shared legend pulled from a representative scatter
  legend <- cowplot::get_legend(
    make_scatter(grna_ids[[1]], build_gdf(grna_ids[[1]])) +
      ggplot2::theme(legend.position = "bottom")
  )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(sprintf("non-pert tail score t_i (%s)", name_),
                        fontface = "bold", x = 0.01, hjust = 0)

  cowplot::plot_grid(title, body, legend, ncol = 1,
                     rel_heights = c(0.04, 1, 0.06))
}


# ---- single-run y vs offset, coloured by estimated perturbation -------------

# For a selection of guides, scatter the observed guide count y_i (y-axis)
# against the offset o_i = eta_i (x-axis, link scale), one facet per guide.
# Each point is COLOURED by the run's ESTIMATED perturbation status (whether the
# cell is in per_guide[[g]]$assignments) and, when ground truth is available,
# SHAPED by the TRUE perturbation status (true_perts[g, i] == 1). This is the
# single-run analog of plot_assignment_groups_y_vs_o: it works on ONE extras
# object, so it is usable on real data where `true_perts` is unknown (omit it and
# the shape aesthetic is dropped, leaving estimated status as colour only).
#
#   extras                 : a loaded extras object
#   grna_matrix            : sparse gRNA count matrix (y_i source; rownames = grna_ids)
#   true_perts             : optional latent perturbation matrix [gRNA x cell];
#                            NULL (default) => no ground truth, no shape aesthetic
#   input_dir              : dataset dir (loads X for o_i if covariate_matrix NULL)
#   covariate_matrix       : optional precomputed X (skips that reload)
#   grna_ids / n_guides    : which / how many guides (facets); default first n_guides
#   max_points_per_panel   : subsample cap per guide facet
#   name_                  : run label for the title (default run_label())
#   log_y                  : log1p y-axis (counts include zeros); default TRUE
#   cap_y                  : optional cap on y_i (counts pmin()'d before plotting)
#   show_true              : restrict cells by TRUE perturbation status:
#                            "both" (default) keeps all cells; "np" keeps only
#                            truly non-perturbed (is_pert == 0); "p" keeps only
#                            truly perturbed (is_pert == 1). "np"/"p" require
#                            `true_perts` (errors otherwise).
#   show_est               : restrict cells by ESTIMATED perturbation status (the
#                            run's own assignments), same values as `show_true`:
#                            "both" (default) keeps all; "np" keeps only cells
#                            estimated non-perturbed; "p" only those estimated
#                            perturbed. Always available (no ground truth needed).
#   exclude_0              : if TRUE, drop all cells with y_i == 0 (applied on top
#                            of `show_true` / `show_est`); default FALSE
#   color_is_pert          : controls which status maps to which aesthetic. FALSE
#                            (default): colour = estimated status, shape = true
#                            status. TRUE: colour = true status, shape = estimated
#                            status. The TRUE-status aesthetic only appears when
#                            `true_perts` is supplied; TRUE here requires it.
#   facet_scales           : passed to facet_wrap (default "free")
#   seed                   : RNG seed for subsampling
plot_y_vs_o_one_run <- function(extras,
                                grna_matrix,
                                true_perts           = NULL,
                                input_dir            = NULL,
                                covariate_matrix     = NULL,
                                grna_ids             = NULL,
                                n_guides             = 6L,
                                max_points_per_panel = 2000L,
                                name_                = NULL,
                                log_y                = TRUE,
                                cap_y                = NULL,
                                show_true            = c("both", "np", "p"),
                                show_est             = c("both", "np", "p"),
                                exclude_0            = FALSE,
                                color_is_pert        = FALSE,
                                facet_scales         = "free",
                                seed                 = 1L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_y_vs_o_one_run requires the ggplot2 package.")
  }
  show_true <- match.arg(show_true)
  show_est  <- match.arg(show_est)
  if (is.null(covariate_matrix)) {
    if (is.null(input_dir)) {
      stop("Provide `covariate_matrix`, or `input_dir` to build it from.")
    }
    covariate_matrix <- build_covariate_matrix(input_dir)
  }
  if (is.null(name_)) name_ <- run_label(extras)
  has_truth <- !is.null(true_perts)
  if (show_true != "both" && !has_truth) {
    stop("show_true = '", show_true, "' filters by true perturbation status, ",
         "but `true_perts` was not supplied.")
  }

  eta_all <- reconstruct_offset_eta(extras, covariate_matrix)
  common  <- intersect(names(eta_all), rownames(grna_matrix))
  if (length(common) == 0L) {
    stop("No guides with offset fits also present in grna_matrix.")
  }
  if (is.null(grna_ids)) {
    grna_ids <- utils::head(common, n_guides)
  } else {
    grna_ids <- intersect(grna_ids, common)
    if (length(grna_ids) == 0L) {
      stop("None of the requested grna_ids have offset fits / are in grna_matrix.")
    }
  }

  # Facet labels carry each guide's converged gamma.
  g_pert <- extract_em_g_pert(extras)
  facet_labels <- stats::setNames(
    sprintf("%s\ng=%.2f", grna_ids, g_pert[grna_ids]),
    grna_ids
  )

  n_cells <- nrow(covariate_matrix)
  set.seed(seed)
  per_guide_df <- lapply(grna_ids, function(id) {
    o   <- eta_all[[id]]                            # offset (link scale)
    y   <- as.numeric(grna_matrix[id, ])            # observed guide count
    est <- logical(n_cells); est[extras$per_guide[[id]]$assignments] <- TRUE
    tp  <- if (has_truth) as.numeric(true_perts[id, ]) == 1 else NA
    gdf <- data.frame(grna_id = id, facet_label = facet_labels[[id]],
                      o = o, y = y, est_pert = est, true_pert = tp,
                      row.names = NULL)
    if (show_true == "np") {
      gdf <- gdf[!gdf$true_pert, , drop = FALSE]
    } else if (show_true == "p") {
      gdf <- gdf[gdf$true_pert, , drop = FALSE]
    }
    if (show_est == "np") {
      gdf <- gdf[!gdf$est_pert, , drop = FALSE]
    } else if (show_est == "p") {
      gdf <- gdf[gdf$est_pert, , drop = FALSE]
    }
    if (exclude_0) {
      gdf <- gdf[gdf$y != 0, , drop = FALSE]
    }
    if (nrow(gdf) > max_points_per_panel) {
      gdf <- gdf[sample.int(nrow(gdf), max_points_per_panel), , drop = FALSE]
    }
    gdf
  })
  df <- do.call(rbind, per_guide_df)
  df$facet_label <- factor(df$facet_label, levels = unname(facet_labels))
  df$est_pert    <- factor(df$est_pert, levels = c(FALSE, TRUE))
  if (has_truth) df$true_pert <- factor(df$true_pert, levels = c(FALSE, TRUE))

  if (!is.null(cap_y)) df$y <- pmin(cap_y, df$y)

  # color_is_pert decides which status drives colour vs shape. The TRUE-status
  # aesthetic only appears when ground truth is available; the ESTIMATED-status
  # aesthetic is always present.
  if (color_is_pert && !has_truth) {
    stop("color_is_pert = TRUE colours by true perturbation status, but ",
         "`true_perts` was not supplied.")
  }
  colour_var <- if (color_is_pert) "true_pert" else "est_pert"
  shape_var  <- if (color_is_pert) "est_pert"  else "true_pert"
  colour_lab <- if (color_is_pert) "true perturbed" else "estimated perturbed"
  shape_lab  <- if (color_is_pert) "estimated perturbed" else "true perturbed"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = o, y = y,
                                        colour = .data[[colour_var]])) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_colour_manual(
      name   = colour_lab,
      values = c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4"),  # ggplot2 default hues
      drop   = FALSE
    ) +
    ggplot2::facet_wrap(~ facet_label, scales = facet_scales) +
    ggplot2::labs(
      x     = "offset  o_i = eta_i  (link scale)",
      y     = "observed guide count  y_i",
      title = sprintf("y vs offset, coloured by %s (%s)", colour_lab, name_),
      subtitle = paste0(length(grna_ids), " guides; up to ",
                        max_points_per_panel, " cells per panel")
    ) +
    ggplot2::theme_bw()

  # The shape aesthetic needs the TRUE status; add it only when available.
  if (has_truth) {
    p <- p + ggplot2::aes(shape = .data[[shape_var]]) +
      ggplot2::scale_shape_manual(
        name   = shape_lab,
        values = c(`FALSE` = 1, `TRUE` = 2),
        drop   = FALSE
      )
  }

  # Y axis: always anchored at 0, with readable ticks. Breaks are supplied as
  # functions so they recompute per panel under facet_scales = "free"; we then
  # force 0 into every panel's range with expand_limits() rather than a global
  # `limits` (which would defeat the free scales).
  if (log_y) {
    # Under log1p the default breaks bunch up at the top; lay down 1-3-10 style
    # breaks spanning [0, panel max].
    cand <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000,
              10000, 30000, 1e5, 3e5, 1e6)
    log1p_breaks <- function(limits) {
      b <- cand[cand <= max(limits, na.rm = TRUE)]
      if (length(b) < 2L) c(0, max(limits, na.rm = TRUE)) else b
    }
    p <- p + ggplot2::scale_y_continuous(
      trans  = "log1p",
      breaks = log1p_breaks,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    )
  }
  p <- p + ggplot2::expand_limits(y = 0)
  p
}


# ============================================================================
# Gated logistic / Gamma-Gaussian variant (run_meta$model == "gated-gamma-gauss")
# ============================================================================
# These operate on extras from the gated variant, whose per_guide elements carry
# a $gated sublist (a0, a1, alpha, beta, mu, sigma, mean_pi, mean_gamma, n_fit)
# instead of an offset model. See guide-assignment-pipeline/bin/script/lib/
# gated_assign.R for the model.

# Parse the "guide type" (data-generating process) from gRNA ids of the form
# "<sim details>_<idx>", dropping the trailing _<idx> replicate index, so all
# iid copies of a setting share a type. E.g.
#   "sim_NPhighvar_Phighdisp_66" -> "sim_NPhighvar_Phighdisp".
guide_type_from_id <- function(ids) sub("_[0-9]+$", "", ids)

# Pull the gated-model per-guide parameters into a data frame (one row per
# guide, names = grna_ids). Columns: the 6 model params (a0, a1, alpha, beta,
# mu, sigma), plus mean_pi / mean_gamma / n_fit, the derived non-pert log-scale
# mean gamma_mean = alpha/beta, and em_converged. Guides below the fit cutoff or
# where EM failed have NA params.
extract_gated_params <- function(extras) {
  if (is.null(extras$per_guide)) {
    stop("`extras` has no $per_guide; is this an assignment_extras_*.rds object?")
  }
  get_field <- function(f) {
    vapply(extras$per_guide, function(r) {
      v <- r$gated[[f]]
      if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
    }, numeric(1))
  }
  df <- data.frame(
    grna_id    = names(extras$per_guide),
    a0         = get_field("a0"),
    a1         = get_field("a1"),
    alpha      = get_field("alpha"),
    beta       = get_field("beta"),
    mu         = get_field("mu"),
    sigma      = get_field("sigma"),
    mean_pi    = get_field("mean_pi"),
    mean_gamma = get_field("mean_gamma"),
    n_fit      = get_field("n_fit"),
    stringsAsFactors = FALSE
  )
  df$gamma_mean   <- df$alpha / df$beta
  df$em_converged <- extract_per_guide_scalar(extras, "em_converged")
  rownames(df) <- NULL
  df
}

# Histograms of each gated-model parameter, grouped by guide type. Guide type is
# parsed from the grna_id (guide_type_from_id), so all iid replicates of a
# simulation setting are pooled into one distribution.
#
#   extras   : gated extras object
#   params   : which params to plot (default the 6 model params). Any column of
#              extract_gated_params() is allowed, e.g. "gamma_mean", "mean_pi".
#   bins     : histogram bins (default 40)
#   overlay  : FALSE (default) -> facet_grid(guide_type ~ parameter): one
#              histogram per (type, param) cell, with the x-scale shared down
#              each parameter column (so types are directly comparable). TRUE ->
#              facet_wrap(~ parameter) with the types overlaid (fill =
#              guide_type, semi-transparent) -- compact but can clutter.
#   drop_na  : drop guides with NA params (below cutoff / EM failed); default TRUE
#
# Returns a ggplot object.
plot_gated_param_histograms <- function(extras,
                                        params  = c("a0", "a1", "alpha", "beta", "mu", "sigma"),
                                        bins    = 40,
                                        overlay = FALSE,
                                        drop_na = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_gated_param_histograms().")
  }
  df <- extract_gated_params(extras)
  unknown <- setdiff(params, names(df))
  if (length(unknown)) {
    stop("Unknown param(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(setdiff(names(df), "grna_id"), collapse = ", "))
  }
  df$guide_type <- guide_type_from_id(df$grna_id)

  long <- do.call(rbind, lapply(params, function(p) {
    data.frame(guide_type = df$guide_type, parameter = p, value = df[[p]],
               stringsAsFactors = FALSE)
  }))
  if (drop_na) long <- long[is.finite(long$value), , drop = FALSE]
  if (!nrow(long)) stop("No finite parameter values to plot (all NA?).")
  long$parameter <- factor(long$parameter, levels = params)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = value))
  if (overlay) {
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(fill = guide_type),
                              bins = bins, position = "identity", alpha = 0.5) +
      ggplot2::facet_wrap(~ parameter, scales = "free")
  } else {
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(fill = guide_type),
                              bins = bins, show.legend = FALSE) +
      ggplot2::facet_grid(guide_type ~ parameter, scales = "free")
  }
  p +
    ggplot2::labs(
      title = paste0("Gated-model parameters by guide type (",
                     extras$run_meta$model %||% "gated", ")"),
      x = "parameter value", y = "count", fill = "guide type"
    ) +
    ggplot2::theme_bw()
}
