# Pure-R re-implementation of sceptre's mixture-model gRNA assignment.
#
# Mirrors:
#   - R/mixture_model_functs.R::assign_grnas_to_cells_mixture (driver)
#   - R/mixture_model_functs.R::obtain_em_assignments       (per-guide logic)
#   - src/mixture_functs.cpp::run_reduced_em_algo_cpp        (EM inner loop)
#
# Everything is at the top level so each piece can be called and inspected
# independently. To parallelize with a PSOCK cluster, source this file on
# every worker first:
#
#   cl <- parallel::makeCluster(n)
#   parallel::clusterCall(cl, source, "experiments/sceptre_assign_pure_R.R")
#   res <- sceptre_assign_pure_R(..., cl = cl)
#   parallel::stopCluster(cl)


# ---- EM tolerance (mirrors compute_tolerance_cpp) ---------------------------
compute_em_tolerance_pure_R <- function(curr_log_lik, prev_log_lik) {
  if (curr_log_lik == -Inf || prev_log_lik == -Inf) {
    1.0
  } else {
    abs(curr_log_lik - prev_log_lik) /
      min(abs(curr_log_lik), abs(prev_log_lik))
  }
}


# ---- Offset-model fit contract ----------------------------------------------
# Each offset-model fit fn must return an object with:
#   $fitted.values         (length-N numeric; required, used by EM as g_mus_pert0)
#   $offset_model_summary  (named list; saved per-guide as `offset_model`).
# Variant-specific per-guide fields (n_trimmed, effective_thresh, ...) live in
# $offset_model_summary; variant-level run params (trim_frac, ...) belong on
# attr(offset_model_fit_fn, "spec") so they travel with run_meta.


# ---- Poisson GLM baseline fit ------------------------------------------------
# Returns the fitted mean for each cell under the "no perturbation" model.
# Matches the call in obtain_em_assignments.
fit_baseline_glm_pure_R <- function(g, covariate_matrix) {
  fit <- suppressWarnings(
    stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson())
  )
  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged
  )
  fit
}


# ---- Robust Poisson GLM baseline fit ----------------------------------------
# Drop-in replacement for fit_baseline_glm_pure_R that downweights cells with
# large Pearson residuals (perturbed cells, for heavy-tailed guides). Uses the
# transformation-based MT estimator from `robustbase::glmrob`.
#
# Motivation: the standard MLE fit is contaminated by perturbed cells, which
# inflates `g_mus_pert0` -- particularly for heavy-tailed guides -- and makes
# the mixture EM over-assign. A robust fit treats the perturbed cells as
# outliers from the "no perturbation" model and shrinks their influence.
#
# Returns an object with `fitted.values` and `coefficients` (just like
# stats::glm.fit), so it slots in unchanged downstream.
fit_baseline_glmrob_pure_R <- function(g, covariate_matrix) {
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("`fit_baseline_glmrob_pure_R` requires the `robustbase` package. ",
         "Install it with `install.packages(\"robustbase\")`.")
  }

  # Standardize non-intercept columns to keep glmrob's iteratively-reweighted
  # Fisher information well-conditioned. Iterative robust solvers invert
  # X' W X-like matrices that are fragile when columns have wildly different
  # scales (a common situation with sceptre's covariates: `(Intercept)` is 1,
  # `log(response_n_umis)` is 8-12, `grna_n_umis` can be 0-10000).
  # Standardization is invariant for fitted.values, so downstream code is
  # unaffected.
  X            <- covariate_matrix
  is_intercept <- apply(X, 2L, function(col) length(unique(col)) == 1L)
  centers      <- ifelse(is_intercept, 0, colMeans(X))
  scales       <- ifelse(is_intercept, 1, apply(X, 2L, stats::sd))
  scales[scales == 0] <- 1
  X            <- sweep(X, 2L, centers, "-")
  X            <- sweep(X, 2L, scales,  "/")

  # Use safe column names for the formula interface ("(Intercept)" etc. aren't
  # valid R identifiers).
  safe_names <- paste0("X", seq_len(ncol(X)))
  colnames(X) <- safe_names
  df <- as.data.frame(X)
  df$.y <- g
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )
  fit <- suppressWarnings(
    robustbase::glmrob(
      fml, data = df,
      family = stats::poisson(),
      method = "MT"
    )
  )
  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged,
    # Per-cell robustness weights (length n). w_r downweights large Pearson
    # residuals; w_x downweights high-leverage covariate rows.
    w_r          = fit$w.r,
    w_x          = fit$w.x
  )
  fit
}


# ---- Trimmed robust Poisson GLM baseline fit --------------------------------
# Robust analogue of fit_baseline_glm_trimmed_pure_R: drop the top `trim_frac`
# of cells by g, then run glmrob on the remainder, then evaluate fitted means
# on ALL cells. Standardization stats are computed from the kept rows so the
# MT solver sees a well-conditioned design (same rationale as the untrimmed
# glmrob fit); trimmed cells are then evaluated on that same scale.
fit_baseline_glmrob_trimmed_pure_R <- function(g, covariate_matrix, trim_frac = 0.05) {
  stopifnot(trim_frac >= 0, trim_frac < 1)
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("`fit_baseline_glmrob_trimmed_pure_R` requires the `robustbase` package.")
  }
  n      <- length(g)
  n_trim <- floor(trim_frac * n)
  if (n_trim == 0L) {
    return(fit_baseline_glmrob_pure_R(g, covariate_matrix))
  }

  keep <- rank(-g, ties.method = "first") > n_trim

  effective_trim_frac <- trim_frac
  if (all(g[keep] == 0)) {
    keep                <- rep(TRUE, n)
    effective_trim_frac <- 0
  }

  X            <- covariate_matrix
  is_intercept <- apply(X, 2L, function(col) length(unique(col)) == 1L)
  centers      <- ifelse(is_intercept, 0, colMeans(X[keep, , drop = FALSE]))
  scales       <- ifelse(is_intercept, 1, apply(X[keep, , drop = FALSE], 2L, stats::sd))
  scales[scales == 0] <- 1
  X            <- sweep(X, 2L, centers, "-")
  X            <- sweep(X, 2L, scales,  "/")

  safe_names <- paste0("X", seq_len(ncol(X)))
  colnames(X) <- safe_names
  df <- as.data.frame(X[keep, , drop = FALSE])
  df$.y <- g[keep]
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )
  fit <- suppressWarnings(
    robustbase::glmrob(
      fml, data = df,
      family = stats::poisson(),
      method = "MT"
    )
  )

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(X %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$trim_frac         <- effective_trim_frac
  fit$n_trimmed         <- sum(!keep)

  # Expand the kept-rows-only weights back to full length n, NA at trimmed
  # cells. This keeps w_r / w_x indexed by cell across guides regardless of
  # how many were trimmed.
  w_r_full <- rep(NA_real_, n); w_r_full[keep] <- fit$w.r
  w_x_full <- rep(NA_real_, n); w_x_full[keep] <- fit$w.x

  fit$offset_model_summary <- list(
    coefficients        = fit$coefficients,
    deviance            = fit$deviance,
    iter                = fit$iter,
    converged           = fit$converged,
    effective_trim_frac = effective_trim_frac,
    n_trimmed           = sum(!keep),
    # Per-cell robustness weights (length n). NA marks trimmed cells (not fit).
    # w_r downweights large Pearson residuals; w_x downweights high-leverage
    # covariate rows.
    w_r                 = w_r_full,
    w_x                 = w_x_full
  )
  fit
}


# ---- Trimmed Poisson GLM baseline fit ---------------------------------------
# Alternative robust baseline that's much more predictable than glmrob: drop
# the top `trim_frac` of cells by g, fit stats::glm.fit on the remainder, then
# evaluate fitted means on ALL cells (the trimmed cells receive a baseline
# prediction extrapolated from the bulk, not influenced by themselves).
#
# This avoids the M-estimator pathologies of glmrob (multiple solutions for
# Mqle, slow / sometimes brittle convergence for MT): the Poisson MLE on a
# subset of cells is a convex problem with a unique solution. The only knob is
# what fraction to drop. Its failure modes are easy to reason about: too high
# `trim_frac` -> fit becomes unstable from too few cells; too low -> still
# contaminated by perturbed cells.
#
# Drop-in replacement for fit_baseline_glm_pure_R. Returns an object with
# `fitted.values` (length N) and `coefficients` (length p), plus some
# trim-bookkeeping fields for inspection.
# ---- Threshold-based trimmed Poisson GLM baseline fit -----------------------
# Same idea as fit_baseline_glm_trimmed_pure_R, but interpretable: trim a cell
# if its guide count, RESCALED to a cell of average library size, exceeds
# `thresh`. Concretely:
#
#   trim cell i  iff  g_i * (mean(lib_size) / lib_size_i)  >=  thresh
#
# i.e. "if we observed cell i's guide UMIs in an average-sized cell, would
# that be >= thresh?". This naturally accounts for capture-efficiency
# differences: a big cell needs a proportionally larger absolute g_i to get
# flagged.
#
# `lib_size` is a length-N vector of per-cell library sizes (typically
# `response_n_umis` from the covariate data frame, but you can pass anything).
# `thresh` is the count threshold on the rescaled scale.
#
# Like the rank-based variant, falls back to fitting on all cells if the kept
# set is all zeros, and flags it via thresh = Inf in the returned fit.
fit_baseline_glm_threshold_pure_R <- function(g, covariate_matrix, lib_size,
                                              thresh = 10) {
  n <- length(g)
  stopifnot(length(lib_size) == n, all(lib_size > 0),
            length(thresh) == 1L, is.finite(thresh), thresh > 0)

  scaled_g <- g * (mean(lib_size) / lib_size)
  keep     <- scaled_g < thresh

  effective_thresh <- thresh
  if (all(g[keep] == 0)) {
    keep             <- rep(TRUE, n)
    effective_thresh <- Inf
  }

  fit <- suppressWarnings(stats::glm.fit(
    y      = g[keep],
    x      = covariate_matrix[keep, , drop = FALSE],
    family = stats::poisson()
  ))

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$thresh            <- effective_thresh
  fit$n_trimmed         <- sum(!keep)
  fit$offset_model_summary <- list(
    coefficients     = fit$coefficients,
    deviance         = fit$deviance,
    iter             = fit$iter,
    converged        = fit$converged,
    effective_thresh = effective_thresh,
    n_trimmed        = sum(!keep)
  )
  fit
}


fit_baseline_glm_trimmed_pure_R <- function(g, covariate_matrix, trim_frac = 0.05) {
  stopifnot(trim_frac >= 0, trim_frac < 1)
  n      <- length(g)
  n_trim <- floor(trim_frac * n)
  if (n_trim == 0L) {
    return(fit_baseline_glm_pure_R(g, covariate_matrix))
  }

  # Drop the top n_trim cells by g (deterministic; ties broken by index).
  keep <- rank(-g, ties.method = "first") > n_trim

  # Fall back to all cells if the kept set is uninformative (all zeros).
  # Poisson glm.fit on an all-zero response gives a degenerate intercept-only
  # fit at -Inf with NA slopes; not useful. Flag via effective_trim_frac = 0.
  effective_trim_frac <- trim_frac
  if (all(g[keep] == 0)) {
    keep <- rep(TRUE, n)
    effective_trim_frac <- 0
  }

  fit <- suppressWarnings(stats::glm.fit(
    y      = g[keep],
    x      = covariate_matrix[keep, , drop = FALSE],
    family = stats::poisson()
  ))

  # Evaluate fitted means on ALL cells using the trimmed-fit coefficients.
  coef <- fit$coefficients
  coef[is.na(coef)] <- 0                  # guard against rank-deficient fits
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$trim_frac         <- effective_trim_frac
  fit$n_trimmed         <- sum(!keep)
  fit$offset_model_summary <- list(
    coefficients        = fit$coefficients,
    deviance            = fit$deviance,
    iter                = fit$iter,
    converged           = fit$converged,
    effective_trim_frac = effective_trim_frac,
    n_trimmed           = sum(!keep)
  )
  fit
}


# ---- NB GLM baseline fit (MASS::glm.nb) -------------------------------------
# Iterative MLE of (beta, theta) under a NB GLM via MASS::glm.nb. Returns
# fitted means (used by EM as g_mus_pert0) plus the estimated theta on the
# fit object, so a paired `estimate_phi_fn` can read it directly without
# re-estimating.
fit_baseline_glm_nb <- function(g, covariate_matrix) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`fit_baseline_glm_nb` requires the `MASS` package.")
  }
  # MASS::glm.nb needs a formula + data interface, not glm.fit's x/y. The "-1"
  # only suppresses adding *another* intercept on top of `covariate_matrix`;
  # the original intercept column from the user's formula (named "(Intercept)"
  # by model.matrix) is renamed X1 here and stays as a predictor. So this is
  # equivalent to calling glm.fit(y = g, x = covariate_matrix, ...) as in the
  # other fit fns -- intercept iff the user's formula had one.
  safe_names <- paste0("X", seq_len(ncol(covariate_matrix)))
  df <- as.data.frame(covariate_matrix)
  colnames(df) <- safe_names
  df$.y <- g
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )

  # Try glm.nb directly (fast path). If it fails ("no valid set of coefficients
  # has been found: please supply starting values" on near-Poisson or
  # ill-conditioned designs), fall back to a Poisson glm.fit with theta
  # estimated separately via sceptre:::estimate_theta -- the same routine used
  # by estimate_phi_from_offset_fit_sceptre. Downstream
  # estimate_phi_from_offset_fit_theta then reads $theta and returns this
  # estimated value (em_phi_source = "estimated") rather than triggering the
  # scalar-phi fallback.
  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula = fml, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    if (!requireNamespace("sceptre", quietly = TRUE)) {
      stop("fit_baseline_glm_nb's Poisson fallback requires the `sceptre` package.")
    }
    fit <- suppressWarnings(
      stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson())
    )
    fit$theta <- tryCatch(
      sceptre:::estimate_theta(
        y     = g,
        mu    = fit$fitted.values,
        dfr   = fit$df.residual,
        limit = 50L,
        eps   = .Machine$double.eps^(1/4)
      )[[1]],
      error = function(e) NA_real_
    )
    fit$SE.theta  <- NA_real_
    fit$twologlik <- NA_real_
  }

  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged,
    theta        = fit$theta,
    SE_theta     = fit$SE.theta,
    twologlik    = fit$twologlik
  )
  fit
}


# ---- Trimmed NB GLM baseline fit (MASS::glm.nb on trimmed sample) -----------
# NB analogue of fit_baseline_glm_trimmed_pure_R: drop the top `trim_frac` of
# cells by g, run MASS::glm.nb on the remainder, then evaluate fitted means on
# ALL cells. Same Poisson + sceptre:::estimate_theta fallback as
# fit_baseline_glm_nb when glm.nb fails; the resulting `$theta` is preserved so
# estimate_phi_from_offset_fit_theta works downstream.
fit_baseline_glm_nb_trimmed_pure_R <- function(g, covariate_matrix, trim_frac = 0.05) {
  stopifnot(trim_frac >= 0, trim_frac < 1)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`fit_baseline_glm_nb_trimmed_pure_R` requires the `MASS` package.")
  }
  n      <- length(g)
  n_trim <- floor(trim_frac * n)
  if (n_trim == 0L) {
    return(fit_baseline_glm_nb(g, covariate_matrix))
  }

  keep <- rank(-g, ties.method = "first") > n_trim

  effective_trim_frac <- trim_frac
  if (all(g[keep] == 0)) {
    keep                <- rep(TRUE, n)
    effective_trim_frac <- 0
  }

  safe_names <- paste0("X", seq_len(ncol(covariate_matrix)))
  df <- as.data.frame(covariate_matrix[keep, , drop = FALSE])
  colnames(df) <- safe_names
  df$.y <- g[keep]
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )

  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula = fml, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    if (!requireNamespace("sceptre", quietly = TRUE)) {
      stop("fit_baseline_glm_nb_trimmed_pure_R's Poisson fallback requires the `sceptre` package.")
    }
    fit <- suppressWarnings(stats::glm.fit(
      y      = g[keep],
      x      = covariate_matrix[keep, , drop = FALSE],
      family = stats::poisson()
    ))
    fit$theta <- tryCatch(
      sceptre:::estimate_theta(
        y     = g[keep],
        mu    = fit$fitted.values,
        dfr   = fit$df.residual,
        limit = 50L,
        eps   = .Machine$double.eps^(1/4)
      )[[1]],
      error = function(e) NA_real_
    )
    fit$SE.theta  <- NA_real_
    fit$twologlik <- NA_real_
  }

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$trim_frac         <- effective_trim_frac
  fit$n_trimmed         <- sum(!keep)
  fit$offset_model_summary <- list(
    coefficients        = fit$coefficients,
    deviance            = fit$deviance,
    iter                = fit$iter,
    converged           = fit$converged,
    theta               = fit$theta,
    SE_theta            = fit$SE.theta,
    twologlik           = fit$twologlik,
    effective_trim_frac = effective_trim_frac,
    n_trimmed           = sum(!keep)
  )
  fit
}


# ---- Threshold-based NB GLM baseline fit ------------------------------------
# NB analogue of fit_baseline_glm_nb_trimmed_pure_R, but with a raw count
# threshold instead of a rank-based trim fraction: fit MASS::glm.nb only to
# cells with g <= y_max, then evaluate fitted means on ALL cells (high-count
# cells receive a baseline prediction extrapolated from the bulk, not
# influenced by themselves). This treats the extreme cells -- which are exactly
# the perturbed ones we don't want contaminating the "no perturbation" offset
# -- as out-of-sample.
#
# Same Poisson + sceptre:::estimate_theta fallback as fit_baseline_glm_nb when
# glm.nb fails; the resulting `$theta` is preserved so
# estimate_phi_from_offset_fit_theta works downstream. Falls back to fitting on
# all cells if the kept set is empty or all zeros, flagged via
# effective_y_max = Inf in the returned summary.
fit_baseline_glm_nb_threshold_pure_R <- function(g, covariate_matrix, y_max = 100) {
  stopifnot(length(y_max) == 1L, is.finite(y_max), y_max >= 0)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`fit_baseline_glm_nb_threshold_pure_R` requires the `MASS` package.")
  }
  n <- length(g)

  keep <- g <= y_max

  effective_y_max <- y_max
  if (!any(keep) || all(g[keep] == 0)) {
    keep            <- rep(TRUE, n)
    effective_y_max <- Inf
  }

  safe_names <- paste0("X", seq_len(ncol(covariate_matrix)))
  df <- as.data.frame(covariate_matrix[keep, , drop = FALSE])
  colnames(df) <- safe_names
  df$.y <- g[keep]
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )

  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula = fml, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    if (!requireNamespace("sceptre", quietly = TRUE)) {
      stop("fit_baseline_glm_nb_threshold_pure_R's Poisson fallback requires the `sceptre` package.")
    }
    fit <- suppressWarnings(stats::glm.fit(
      y      = g[keep],
      x      = covariate_matrix[keep, , drop = FALSE],
      family = stats::poisson()
    ))
    fit$theta <- tryCatch(
      sceptre:::estimate_theta(
        y     = g[keep],
        mu    = fit$fitted.values,
        dfr   = fit$df.residual,
        limit = 50L,
        eps   = .Machine$double.eps^(1/4)
      )[[1]],
      error = function(e) NA_real_
    )
    fit$SE.theta  <- NA_real_
    fit$twologlik <- NA_real_
  }

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$y_max             <- effective_y_max
  fit$n_trimmed         <- sum(!keep)
  fit$offset_model_summary <- list(
    coefficients    = fit$coefficients,
    deviance        = fit$deviance,
    iter            = fit$iter,
    converged       = fit$converged,
    theta           = fit$theta,
    SE_theta        = fit$SE.theta,
    twologlik       = fit$twologlik,
    effective_y_max = effective_y_max,
    n_trimmed       = sum(!keep)
  )
  fit
}


# ---- Capped NB GLM baseline fit ---------------------------------------------
# Like fit_baseline_glm_nb_threshold_pure_R, but instead of DROPPING cells with
# g > y_max it CAPS them: fit MASS::glm.nb on pmin(g, y_max). All n cells are
# kept -- capping just bounds the influence (leverage) of any single large
# count on the baseline fit, rather than excluding those cells entirely.
#
# Same Poisson + sceptre:::estimate_theta fallback as fit_baseline_glm_nb when
# glm.nb fails (run on the capped response); the resulting `$theta` is preserved
# so estimate_phi_from_offset_fit_theta works downstream. Note theta here is the
# dispersion of the CAPPED counts, which is generally larger (less overdispersed)
# than the raw-count theta.
fit_baseline_glm_nb_capped_pure_R <- function(g, covariate_matrix, y_max = 100) {
  stopifnot(length(y_max) == 1L, is.finite(y_max), y_max >= 0)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("`fit_baseline_glm_nb_capped_pure_R` requires the `MASS` package.")
  }
  n       <- length(g)
  n_capped <- sum(g > y_max)
  g_cap   <- pmin(g, y_max)

  safe_names <- paste0("X", seq_len(ncol(covariate_matrix)))
  df <- as.data.frame(covariate_matrix)
  colnames(df) <- safe_names
  df$.y <- g_cap
  fml <- stats::as.formula(
    paste0(".y ~ -1 + ", paste(safe_names, collapse = " + "))
  )

  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula = fml, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    if (!requireNamespace("sceptre", quietly = TRUE)) {
      stop("fit_baseline_glm_nb_capped_pure_R's Poisson fallback requires the `sceptre` package.")
    }
    fit <- suppressWarnings(stats::glm.fit(
      y      = g_cap,
      x      = covariate_matrix,
      family = stats::poisson()
    ))
    fit$theta <- tryCatch(
      sceptre:::estimate_theta(
        y     = g_cap,
        mu    = fit$fitted.values,
        dfr   = fit$df.residual,
        limit = 50L,
        eps   = .Machine$double.eps^(1/4)
      )[[1]],
      error = function(e) NA_real_
    )
    fit$SE.theta  <- NA_real_
    fit$twologlik <- NA_real_
  }

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$y_max             <- y_max
  fit$n_capped          <- n_capped
  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged,
    theta        = fit$theta,
    SE_theta     = fit$SE.theta,
    twologlik    = fit$twologlik,
    y_max        = y_max,
    n_capped     = n_capped
  )
  fit
}


# ---- Capped Poisson GLM baseline fit ----------------------------------------
# Poisson analogue of fit_baseline_glm_nb_capped_pure_R: cap the response
# (pmin(g, y_max)) and fit stats::glm.fit with a Poisson family on all n cells.
# Capping bounds the leverage of any single large count without dropping cells.
# No NB size is produced; pair with a Poisson family, or with nb-fixed0 using
# estimate_phi_from_offset_means_nb (theta0 from the means).
fit_baseline_glm_pois_capped_pure_R <- function(g, covariate_matrix, y_max = 100) {
  stopifnot(length(y_max) == 1L, is.finite(y_max), y_max >= 0)
  n_capped <- sum(g > y_max)
  g_cap    <- pmin(g, y_max)

  fit <- suppressWarnings(stats::glm.fit(
    y      = g_cap,
    x      = covariate_matrix,
    family = stats::poisson()
  ))

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0                  # guard against rank-deficient fits
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$y_max             <- y_max
  fit$n_capped          <- n_capped
  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged,
    y_max        = y_max,
    n_capped     = n_capped
  )
  fit
}


# ---- Threshold-based Poisson GLM baseline fit -------------------------------
# Poisson analogue of fit_baseline_glm_nb_threshold_pure_R: fit stats::glm.fit
# (Poisson) only on cells with g <= y_max, then evaluate fitted means on ALL
# cells. Drops the high-count tail (the perturbed cells) from the "no
# perturbation" offset fit, rather than capping them. No NB size is produced;
# pair with a Poisson family, or with nb-fixed0 via
# estimate_phi_from_offset_means_nb. Falls back to fitting on all cells if the
# kept set is empty or all zeros, flagged via effective_y_max = Inf.
fit_baseline_glm_pois_threshold_pure_R <- function(g, covariate_matrix, y_max = 100) {
  stopifnot(length(y_max) == 1L, is.finite(y_max), y_max >= 0)
  n <- length(g)

  keep <- g <= y_max

  effective_y_max <- y_max
  if (!any(keep) || all(g[keep] == 0)) {
    keep            <- rep(TRUE, n)
    effective_y_max <- Inf
  }

  fit <- suppressWarnings(stats::glm.fit(
    y      = g[keep],
    x      = covariate_matrix[keep, , drop = FALSE],
    family = stats::poisson()
  ))

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(covariate_matrix %*% coef)
  fit$fitted.values     <- exp(eta_all)
  fit$linear.predictors <- eta_all
  fit$y_max             <- effective_y_max
  fit$n_trimmed         <- sum(!keep)
  fit$offset_model_summary <- list(
    coefficients    = fit$coefficients,
    deviance        = fit$deviance,
    iter            = fit$iter,
    converged       = fit$converged,
    effective_y_max = effective_y_max,
    n_trimmed       = sum(!keep)
  )
  fit
}


# ---- Capped log1p-OLS baseline fit ------------------------------------------
# A non-GLM offset: cap the counts at y_max, take log1p, and run ordinary least
# squares (stats::lm.fit) of log1p(pmin(g, y_max)) on the covariate matrix. The
# OLS prediction estimates E[log1p(count)]; the baseline MEAN is recovered as
# mu0 = expm1(prediction) (the inverse of log1p), floored at a small positive
# value so log(mu0) is finite for the EM. fitted.values therefore holds mu0
# (the mean), consistent with the GLM offset fits.
#
# No NB size is produced (this isn't an NB fit), so $theta is absent: pair this
# offset with a Poisson family, or with nb-fixed0 using
# estimate_phi_from_offset_means_nb (which derives theta0 from the means, not a
# stored $theta).
fit_baseline_lm_log1p_capped_pure_R <- function(g, covariate_matrix, y_max = 100,
                                                mu_floor = 1e-8) {
  stopifnot(length(y_max) == 1L, is.finite(y_max), y_max >= 0)
  n        <- length(g)
  n_capped <- sum(g > y_max)
  z        <- log1p(pmin(g, y_max))

  fit <- stats::lm.fit(x = covariate_matrix, y = z)

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0                  # guard against rank-deficient fits
  eta_log1p <- as.numeric(covariate_matrix %*% coef)   # predicted log1p(count)
  mu0       <- pmax(expm1(eta_log1p), mu_floor)        # back to the count mean

  fit$fitted_log1p      <- eta_log1p      # OLS prediction on the log1p scale
  fit$fitted.values     <- mu0            # baseline mean mu0 (EM uses this)
  fit$linear.predictors <- log(mu0)       # eta = log(mu0), the actual offset
  fit$y_max             <- y_max
  fit$n_capped          <- n_capped

  # residual SD on the log1p scale, for diagnostics
  resid_df <- n - sum(!is.na(fit$coefficients))
  sigma    <- if (resid_df > 0) sqrt(sum((z - eta_log1p)^2) / resid_df) else NA_real_

  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    sigma_log1p  = sigma,
    y_max        = y_max,
    n_capped     = n_capped
  )
  fit
}


# ---- NB phi (= theta) estimators --------------------------------------------
# An `estimate_phi_fn` is called per guide as
#   estimate_phi_fn(g, offset_model_fit) -> scalar phi
# and the result is then handed to run_em_nb_nb_pure_R as the (shared) NB
# overdispersion. Use one of the helpers below, or write your own that
# matches that contract.

# Estimate NB theta by fixing mu from the Poisson offset fit and running
# sceptre's internal Newton estimator. sceptre:::estimate_theta returns
# `list(theta, method)`; we take the theta only.
estimate_phi_from_offset_fit_sceptre <- function(g, offset_model_fit,
                                                  limit = 50L,
                                                  eps   = .Machine$double.eps^(1/4)) {
  if (!requireNamespace("sceptre", quietly = TRUE)) {
    stop("estimate_phi_from_offset_fit_sceptre requires the sceptre package.")
  }
  sceptre:::estimate_theta(
    y     = g,
    mu    = offset_model_fit$fitted.values,
    dfr   = offset_model_fit$df.residual,
    limit = limit,
    eps   = eps
  )[[1]]
}

# Read theta directly off the offset fit. Use this when the offset model
# already estimated theta jointly (e.g. fit_baseline_glm_nb -> MASS::glm.nb).
estimate_phi_from_offset_fit_theta <- function(g, offset_model_fit) {
  offset_model_fit$theta
}

# Estimate the non-pert NB size from the per-cell baseline MEANS rather than the
# counts. The "data" are mu0_i = exp(offset_i) = offset_model_fit$fitted.values,
# treated as a sample from a single NB with a fixed mean m = exp(mean(offset))
# (the value you proposed). theta is then identified by method of moments:
#   Var = m + m^2 / theta   =>   theta = m^2 / (Var - m),
# where Var is the empirical variance of the mu0_i. Using the fitted means
# instead of g means perturbed cells -- whose large counts otherwise crush
# theta -- cannot contaminate the estimate.
#
# Why a mean is needed: NB is a two-parameter family and the size theta is not
# identified from a variance alone (Var mixes mean and dispersion), so a mean
# must be pinned. Here it is fixed at exp(mean(offset)) per your spec rather
# than estimated as the sample mean.
#
# Note: exp(offset) is a smooth function of the covariates and is usually far
# LESS variable than a count with that mean would be, so Var <= m is common; in
# that case theta is undefined/negative and we return theta_cap (effectively
# Poisson, i.e. a tight non-pert component). g is unused (kept for the
# estimate_phi_fn contract).
estimate_phi_from_offset_means_nb <- function(g, offset_model_fit,
                                              theta_cap = 1e6,
                                              min_theta = 1e-4) {
  mu0    <- offset_model_fit$fitted.values
  offset <- log(mu0)
  m      <- exp(mean(offset))      # fixed NB mean = exp(mean(offset))
  v      <- stats::var(mu0)        # empirical dispersion of the baseline means
  if (!is.finite(m) || !is.finite(v) || v <= m) {
    return(theta_cap)
  }
  theta <- m^2 / (v - m)
  min(max(theta, min_theta), theta_cap)
}


# ---- Phi update given (pi, gamma) from a converged EM -----------------------
# Maximize the marginal log-likelihood of the NB-NB mixture in phi alone,
# holding (pi, gamma) and the offset means fixed. Used by
# run_em_nb_nb_update_phi_pure_R to refine phi between EM passes.
# Numerically stable form: log1p / log + dnbinom(log = TRUE) + logaddexp.
est_phi_given_model_pure_R <- function(g, pi, gamma, g_mus_pert0,
                                       log_phi_lower = -10,
                                       log_phi_upper =  10) {
  mu0       <- g_mus_pert0
  mu1       <- exp(gamma) * g_mus_pert0
  log_1m_pi <- log1p(-pi)
  log_pi    <- log(pi)
  fn <- function(log_phi) {
    phi    <- exp(log_phi)
    log_p0 <- log_1m_pi + stats::dnbinom(g, mu = mu0, size = phi, log = TRUE)
    log_p1 <- log_pi    + stats::dnbinom(g, mu = mu1, size = phi, log = TRUE)
    mean(logaddexp(log_p0, log_p1))
  }
  opt <- stats::optimize(fn, lower = log_phi_lower, upper = log_phi_upper,
                         maximum = TRUE)
  exp(opt$maximum)
}


# ---- Pois-Pois mixture EM (mirrors run_reduced_em_algo_cpp) ------------------
# g_mus_pert0 is the per-cell baseline mean from the Poisson GLM.
# Note: matches the existing C++ exactly, including the M-step
# `curr_g_pert <- log(e1) - log(e2)` that is then *added* to g_mus_pert0
# in the next iteration's E-step. (Flagged as something worth revisiting.)
run_em_pois_pois_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                    log_g_factorial = lgamma(g + 1),
                                    max_iter = 50L, min_iter = 3L,
                                    ep_tol = 0.5e-4, fix_curr_g_pert_bug=FALSE) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0) == n,
            length(log_g_factorial) == n)

  log_g_mus_pert0 <- log(g_mus_pert0)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_

  for (i in seq_len(B)) {
    curr_pi      <- pi_guesses[i]
    curr_g_pert  <- g_pert_guesses[i]
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s         <- numeric(n)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      g_mus_pert1 <- if (fix_curr_g_pert_bug) {
        g_mus_pert0 * exp(curr_g_pert)
      } else {
        g_mus_pert0 + curr_g_pert
      }

      log_g_mus_pert1 <- log(g_mus_pert1)

      # mixture density (per-cell), and total log-likelihood
      p0 <- exp(log(1 - curr_pi) + g * log_g_mus_pert0 - g_mus_pert0 - log_g_factorial)
      p1 <- exp(log(curr_pi)     + g * log_g_mus_pert1 - g_mus_pert1 - log_g_factorial)
      s  <- p0 + p1
      s[s < 1e-100] <- 1e-100
      curr_log_lik <- sum(log(s))

      # E-step: posterior P(perturbed | g_i), computed via the "quotient" form
      # used in the C++ for numerical stability.
      quotient <- log(1 - curr_pi) - log(curr_pi) +
                  g * (log_g_mus_pert0 - log_g_mus_pert1) +
                  g_mus_pert1 - g_mus_pert0
      Ti1s <- 1 / (exp(quotient) + 1)

      # bail out if degenerate
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # update pi; flip if it crossed 0.5 (label-switching guard)
      curr_pi <- sum(Ti1s) / n
      if (curr_pi > 0.5) {
        Ti1s    <- 1 - Ti1s
        curr_pi <- 1 - curr_pi
      }

      # convergence check (must happen *before* the iteration counter bump
      # and *before* the M-step, to match the C++ control flow exactly)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for the perturbation effect
      e1 <- sum(Ti1s * g)
      e2 <- sum(Ti1s * g_mus_pert0)
      curr_g_pert <- log(e1) - log(e2)
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_g_pert
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert)
}


# ---- Weighted quantile + robust (weighted-quantile) gamma M-step ------------
# Inverse-CDF (type-1) weighted quantile: sort by x, normalize cumulative
# weights, return the smallest x whose cumulative weight reaches each prob.
weighted_quantile <- function(x, w, probs = 0.5) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]

  if (length(x) == 0L || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }

  o <- order(x)
  x <- x[o]
  w <- w[o]

  cw <- cumsum(w) / sum(w)

  vapply(probs, function(p) {
    x[which(cw >= p)[1L]]
  }, numeric(1))
}

# Robust replacement for the multiplicative Poisson M-step
#   exp(gamma) = sum(r * y) / sum(r * mu0).
# That closed form is the weighted MEAN of s_i = y_i / mu0_i with weights
# w_i = r_i * mu0_i; here we instead take a weighted QUANTILE of s (the
# weighted median by default), which resists a few high-count cells dragging
# gamma up. gamma is clamped to [gamma_min, gamma_max]; gamma_min = 0 keeps the
# perturbed mean at or above baseline. Note: a weighted quantile <= 0.5 can
# collapse to gamma_min when the perturbed mean is small (many perturbed cells
# still draw y = 0, so s = 0); raising `prob` steps over that zero-spike.
gamma_update_weighted_quantile <- function(y, mu0, r, prob = 0.5,
                                           gamma_min = 0, gamma_max = 30) {
  mu0 <- pmax(mu0, 1e-300)

  alpha_i <- y / mu0
  w <- r * mu0

  alpha <- weighted_quantile(alpha_i, w, probs = prob)

  if (!is.finite(alpha)) {
    return(gamma_min)
  }

  gamma <- log(alpha)

  min(max(gamma, gamma_min), gamma_max)
}


# ---- Pois-Pois mixture EM with a robust weighted-quantile gamma M-step ------
# Same scaffolding as run_em_pois_pois_pure_R (outer loop over B starts,
# sceptre's relative-tol convergence test, degeneracy bail-out), but:
#
#   1. E-step is always the multiplicative ("unbugged") model
#      mu_pert1 = mu0 * exp(gamma); the s = y/mu0, gamma = log(alpha)
#      formulation has no additive analogue, so there is no
#      `fix_curr_g_pert_bug` toggle.
#   2. M-step replaces exp(gamma) = sum(r*y)/sum(r*mu0) (a weighted mean of
#      y/mu0) with a weighted QUANTILE of y/mu0 -- see
#      gamma_update_weighted_quantile.
#
# The label-switch flip from the Poisson/NB EMs is intentionally dropped:
# gamma_min = 0 pins component 1 as the elevated (perturbed) component, so
# flipping would mislabel it.
#
# Caveat: a quantile M-step is an M-estimator, not the Q-maximizer, so the
# observed-data log-likelihood is not guaranteed to increase monotonically;
# the relative-tol convergence test can in principle oscillate (rare in
# practice). outer_log_lik is the observed Poisson-mixture loglik of the best
# converged start, used only to pick among starts. outer_g_pert returns the
# converged gamma (log-fold-change).
run_em_pois_pois_wquantile_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                              log_g_factorial = lgamma(g + 1),
                                              max_iter = 50L, min_iter = 3L,
                                              ep_tol = 0.5e-4,
                                              gamma_update_prob = 0.5,
                                              gamma_min = 0, gamma_max = 30) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0) == n,
            length(log_g_factorial) == n)

  log_g_mus_pert0 <- log(g_mus_pert0)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_

  for (i in seq_len(B)) {
    curr_pi      <- pi_guesses[i]
    curr_g_pert  <- g_pert_guesses[i]
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s         <- numeric(n)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      # multiplicative ("unbugged") perturbed mean
      g_mus_pert1     <- g_mus_pert0 * exp(curr_g_pert)
      log_g_mus_pert1 <- log(g_mus_pert1)

      p0 <- exp(log(1 - curr_pi) + g * log_g_mus_pert0 - g_mus_pert0 - log_g_factorial)
      p1 <- exp(log(curr_pi)     + g * log_g_mus_pert1 - g_mus_pert1 - log_g_factorial)
      s  <- p0 + p1
      s[s < 1e-100] <- 1e-100
      curr_log_lik <- sum(log(s))

      # E-step: posterior P(perturbed | g_i), quotient form for stability.
      quotient <- log(1 - curr_pi) - log(curr_pi) +
                  g * (log_g_mus_pert0 - log_g_mus_pert1) +
                  g_mus_pert1 - g_mus_pert0
      Ti1s <- 1 / (exp(quotient) + 1)

      # bail out if degenerate
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # update pi (no label-switch flip; see header)
      curr_pi <- sum(Ti1s) / n

      # convergence check (before the M-step, matching the Poisson EM)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # robust M-step: weighted quantile of y/mu0 in place of the weighted mean
      curr_g_pert <- gamma_update_weighted_quantile(
        y = g, mu0 = g_mus_pert0, r = Ti1s,
        prob = gamma_update_prob, gamma_min = gamma_min, gamma_max = gamma_max
      )
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_g_pert
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert)
}


# ---- Additive Poisson mixture EM (correct additive M-step) ------------------
# Mixture: g_i ~ (1 - pi) * Pois(g_mus_pert0_i) + pi * Pois(g_mus_pert0_i + delta),
# with delta >= 0 an additive count bump. Same outer-loop scaffolding as
# run_em_pois_pois_pure_R (B starting guesses, sceptre's relative-tol
# convergence test, label-switch guard, degeneracy bail-out), but:
#
#   1. E-step uses Pois(mu0 + delta) for the perturbed component (additive
#      mean), not Pois(mu0 * exp(gamma)).
#   2. M-step for delta is the actual MLE: the root of
#        sum(r * (g / (mu0 + delta) - 1))
#      via uniroot with adaptive bracket expansion. Boundary check at
#      delta = min_delta returns the boundary if the score is non-positive
#      there (kkt-style).
#
# Contrast with run_em_pois_pois_pure_R(fix_curr_g_pert_bug = FALSE), which
# pairs the additive mean with the multiplicative M-step `log(e1) - log(e2)`
# (the "bugged" variant). This kernel implements the corresponding correct
# additive M-step.
#
# g_pert_guesses is interpreted as log(delta) seeds so the same driver default
# g_pert_guess_range = log(c(10, 5000)) initializes delta in [10, 5000] count
# units, a sensible additive scale for sceptre guides.
run_em_pois_additive_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                        log_g_factorial = lgamma(g + 1),
                                        max_iter  = 50L, min_iter = 3L,
                                        ep_tol    = 0.5e-4,
                                        min_delta = 1e-10) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses)  == B,
            length(g_mus_pert0)     == n,
            length(log_g_factorial) == n)

  # Floor on baseline mean for log-stability in the Poisson density.
  mu0     <- pmax(g_mus_pert0, 1e-300)
  log_mu0 <- log(mu0)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_

  # 1D MLE for delta given responsibilities r. Closes over (g, mu0).
  # Adaptive bracket: hi starts at max(2*old_delta, max(g), ...) and doubles
  # until the score turns negative, or 100 doublings give up.
  mstep_delta <- function(r, old_delta) {
    if (sum(r) <= 0) return(max(old_delta, min_delta))
    score <- function(d) sum(r * (g / (mu0 + d) - 1))

    s_lo <- score(min_delta)
    if (!is.finite(s_lo) || s_lo <= 0) return(min_delta)

    hi   <- max(old_delta * 2, max(g), mean(g) + mean(mu0), 1)
    hi   <- max(hi, min_delta * 2)
    s_hi <- score(hi)
    k    <- 0L
    while (is.finite(s_hi) && s_hi > 0 && k < 100L) {
      hi <- hi * 2; s_hi <- score(hi); k <- k + 1L
    }
    if (!is.finite(s_hi) || s_hi > 0) return(hi)
    stats::uniroot(score, lower = min_delta, upper = hi, tol = 1e-8)$root
  }

  for (i in seq_len(B)) {
    curr_pi      <- pi_guesses[i]
    curr_delta   <- max(exp(g_pert_guesses[i]), min_delta)
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s         <- numeric(n)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      mu1     <- mu0 + curr_delta
      log_mu1 <- log(mu1)

      # mixture density per cell + total log-lik
      log_p0 <- log1p(-curr_pi) + g * log_mu0 - mu0 - log_g_factorial
      log_p1 <- log(curr_pi)    + g * log_mu1 - mu1 - log_g_factorial
      curr_log_lik <- sum(logaddexp(log_p0, log_p1))

      # E-step: posterior P(perturbed | g_i)
      Ti1s <- stats::plogis(log_p1 - log_p0)

      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # pi update + label-switch guard (mirrors run_em_pois_pois_pure_R)
      curr_pi <- sum(Ti1s) / n
      if (curr_pi > 0.5) {
        Ti1s    <- 1 - Ti1s
        curr_pi <- 1 - curr_pi
      }

      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for delta: 1D MLE via uniroot
      curr_delta <- max(mstep_delta(Ti1s, curr_delta), min_delta)
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_delta
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert)
}


# ---- Additive ZTP mixture EM (positive counts only) -------------------------
# Same additive model as run_em_pois_additive_pure_R, but only positive cells
# enter the mixture; zero cells always get posterior = 0.
#
#   g_i | g_i > 0  ~  (1 - rho) * ZTP(mu0_i) + rho * ZTP(mu0_i + delta)
#
# rho is the fraction of *positive* cells in the elevated component and can
# legitimately be large, so:
#   (a) No label-switch flip — delta >= 0 makes component ordering identifiable.
#   (b) rho starting guesses are spread evenly across (0.05, 0.95); the driver
#       default pi_guess_range = c(1e-5, 0.1) is too narrow for this model.
#
# M-step for delta: optimize() over log(delta) (unchanged from prototype —
# ZTP log-lik has no simple closed-form score root). outer_g_pert returns
# delta (additive count shift).
run_em_pois_additive_nonzero_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                                 log_g_factorial = NULL,
                                                 max_iter  = 50L, min_iter = 3L,
                                                 ep_tol    = 0.5e-4,
                                                 min_delta = 1e-10) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0)    == n)

  g   <- as.numeric(g)
  mu0 <- pmax(as.numeric(g_mus_pert0), 1e-300)

  pos   <- g > 0
  gp    <- g[pos]
  mu0p  <- mu0[pos]
  n_pos <- sum(pos)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_

  if (n_pos == 0L) {
    return(list(outer_Ti1s      = outer_Ti1s,
                outer_i         = 0L,
                outer_converged = TRUE,
                outer_log_lik   = 0,
                outer_pi        = NA_real_,
                outer_g_pert    = NA_real_))
  }

  # Stable log(1 - exp(-mu)) for mu > 0.
  log1mexp <- function(x) {
    out        <- numeric(length(x))
    small      <- x < log(2)
    out[small]  <- log(-expm1(-x[small]))
    out[!small] <- log1p(-exp(-x[!small]))
    out
  }

  log_ztpois <- function(y, mu) {
    stats::dpois(y, lambda = mu, log = TRUE) - log1mexp(mu)
  }

  delta_upper <- max(10, 10 * max(gp), 2 * max(mu0p + gp))

  mstep_delta <- function(r_pos, old_delta) {
    if (sum(r_pos) < 1e-8) return(old_delta)
    res <- tryCatch(
      stats::optimize(
        f        = function(eta) -sum(r_pos * log_ztpois(gp, mu0p + exp(eta))),
        interval = c(log(min_delta), log(delta_upper))
      ),
      error = function(e) NULL
    )
    if (is.null(res)) return(old_delta)
    max(exp(res$minimum), min_delta)
  }

  rho_guesses <- seq(0.05, 0.95, length.out = B)

  for (i in seq_len(B)) {
    curr_rho   <- rho_guesses[i]
    curr_delta <- max(exp(g_pert_guesses[i]), min_delta)
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s_pos     <- numeric(n_pos)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      log0 <- log1p(-curr_rho) + log_ztpois(gp, mu0p)
      log1 <- log(curr_rho)    + log_ztpois(gp, mu0p + curr_delta)
      curr_log_lik <- sum(logaddexp(log0, log1))
      Ti1s_pos <- stats::plogis(log1 - log0)

      if (all(Ti1s_pos <= 1e-100) || any(!is.finite(Ti1s_pos))) {
        curr_log_lik <- -Inf
        break
      }

      curr_rho <- min(max(mean(Ti1s_pos), 1e-10), 1 - 1e-10)

      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      curr_delta <- mstep_delta(Ti1s_pos, curr_delta)
    }

    if (converged && curr_log_lik > outer_log_lik) {
      Ti1s_full       <- numeric(n)
      Ti1s_full[pos]  <- Ti1s_pos
      outer_Ti1s      <- Ti1s_full
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_rho
      outer_g_pert    <- curr_delta
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert)
}


# ---- NB primitives (from experiments/pure_R_nb_em.R) ------------------------
# softplus(x) = log(1 + exp(x)), numerically stable.
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# logaddexp(x, y) = log(exp(x) + exp(y)), elementwise, numerically stable.
logaddexp <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log1p(exp(-abs(x - y)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

# Per-cell observed log-likelihood under the NB-NB mixture.
# Returns a length-n vector; sum() it to get the total log-lik.
# Uses the softplus parameterization t = offset - log(phi), exploiting the
# identity log dnbinom(y; phi, mu = exp(eta)) = -phi*softplus(eta - log(phi))
#                                              - y*softplus(-(eta - log(phi)))
# plus constants in (y, phi) that drop out of the EM (E-step quotient, and
# differences of log-liks across iterations).
obs_ll_vec <- function(gamma, pi, y, offset, phi) {
  t0 <- offset - log(phi)
  t1 <- gamma + offset - log(phi)
  a0 <- log1p(-pi) - phi * softplus(t0) - y * softplus(-t0)
  a1 <- log(pi)    - phi * softplus(t1) - y * softplus(-t1)
  logaddexp(a0, a1)
}

# Score of the Q-function in gamma. Setting to zero gives the M-step MLE.
# Derived from log dnbinom(y; phi, mu = exp(offset + gamma)) under weighting
# by prob_is_1.
gamma_score <- function(gamma_, y_plus_phi, offset_minus_log_phi, phi, prob_is_1) {
  sum(prob_is_1 * (y_plus_phi * plogis(-gamma_ - offset_minus_log_phi) - phi))
}


# ---- NB-NB mixture EM (proper joint MLE of gamma and pi, fixed phi) ---------
# Same control-flow scaffolding as run_em_pois_pois_pure_R (outer loop over B
# starting guesses, sceptre's convergence test, label-switch guard,
# degeneracy bail-out), but with two substantive differences:
#
#   1. E-step uses NB-NB densities (multiplicative model: pert mean =
#      mu_pert0 * exp(gamma)). There is no additive variant here -- gamma is
#      always a log-fold-change -- so there is no `fix_curr_g_pert_bug` toggle.
#   2. M-step is the *actual* MLE under the Q-function:
#        - pi <- mean(prob_is_1)        (closed form)
#        - gamma <- root of gamma_score (uniroot)
#      As phi -> Inf, the gamma root collapses to log(sum(prob*y)/sum(prob*mu0))
#      = sceptre's log(e1/e2), so this reduces to the Poisson EM in the
#      large-phi limit (with fix_curr_g_pert_bug = TRUE).
#
# Argument naming follows the Poisson version (g, g_mus_pert0, pi_guesses,
# g_pert_guesses) so it slots into assign_one_grna_pure_R. Internally we use
# offset = log(g_mus_pert0) and gamma = curr_g_pert.
run_em_nb_nb_pure_R <- function(g, g_mus_pert0, phi, pi_guesses, g_pert_guesses,
                                log_g_factorial = NULL,
                                max_iter = 50L, min_iter = 3L,
                                ep_tol = 0.5e-4) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0) == n,
            length(phi) == 1L, is.finite(phi), phi > 0)

  offset               <- log(g_mus_pert0)
  offset_minus_log_phi <- offset - log(phi)
  y_plus_phi           <- g + phi

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_

  for (i in seq_len(B)) {
    curr_pi      <- pi_guesses[i]
    curr_gamma   <- g_pert_guesses[i]
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s         <- numeric(n)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      # total observed log-likelihood at current (gamma, pi)
      ll_per_cell  <- obs_ll_vec(gamma = curr_gamma, pi = curr_pi,
                                 y = g, offset = offset, phi = phi)
      curr_log_lik <- sum(ll_per_cell)

      # E-step: posterior responsibility prob_is_1 in (0, 1)
      # Same softplus form as obs_ll_vec for numerical stability.
      t0   <- offset_minus_log_phi
      t1   <- curr_gamma + offset_minus_log_phi
      qlog <- -phi * softplus(t1) - g * softplus(-t1) +
               phi * softplus(t0) + g * softplus(-t0) +
               qlogis(curr_pi)
      Ti1s <- plogis(qlog)

      # bail out if degenerate
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # update pi; label-switch guard (sceptre style)
      curr_pi <- sum(Ti1s) / n
      if (curr_pi > 0.5) {
        Ti1s    <- 1 - Ti1s
        curr_pi <- 1 - curr_pi
      }

      # convergence check (sceptre's relative-tol form; same as Poisson)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for gamma: MLE via uniroot on the Q-function score.
      # Bracket centered on current gamma with extendInt for robustness.
      root_fit <- tryCatch(
        stats::uniroot(
          gamma_score,
          interval     = c(curr_gamma - 2, curr_gamma + 2),
          extendInt    = "downX",
          check.conv   = TRUE,
          tol          = 1e-8,
          y_plus_phi   = y_plus_phi,
          offset_minus_log_phi = offset_minus_log_phi,
          phi          = phi,
          prob_is_1    = Ti1s
        ),
        error = function(e) NULL
      )
      if (is.null(root_fit)) {
        curr_log_lik <- -Inf
        break
      }
      curr_gamma <- root_fit$root
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_gamma
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert)
}


# ---- NB-NB EM alternating with phi updates ----------------------------------
# Coordinate ascent over (pi, gamma) and phi:
#   for k = 1..(n_phi_updates + 1):
#     EM at current phi -> (pi, gamma)
#     if k <= n_phi_updates: phi <- est_phi_given_model_pure_R(...)
# Total: (n_phi_updates + 1) EM passes, n_phi_updates phi updates. The last
# EM uses the final phi, so the returned (pi, gamma) match it.
#
# If a given EM pass doesn't converge (or returns non-finite pi/gamma), the
# phi update for that round is skipped and curr_phi is retained. The full
# per-pass trajectory is returned for diagnostics.
run_em_nb_nb_update_phi_pure_R <- function(g, g_mus_pert0, phi, pi_guesses,
                                           g_pert_guesses, n_phi_updates,
                                           log_g_factorial = NULL,
                                           max_iter = 50L, min_iter = 3L,
                                           ep_tol = 0.5e-4) {
  stopifnot(length(n_phi_updates) == 1L, is.finite(n_phi_updates),
            n_phi_updates >= 0L)
  stopifnot(length(phi) == 1L, is.finite(phi), phi > 0)

  K <- n_phi_updates + 1L
  phi_traj       <- numeric(K)
  pi_traj        <- rep(NA_real_,    K)
  g_pert_traj    <- rep(NA_real_,    K)
  log_lik_traj   <- rep(NA_real_,    K)
  converged_traj <- rep(NA,          K)
  init_i_traj    <- rep(NA_integer_, K)

  curr_phi <- phi
  em_fit   <- NULL

  for (k in seq_len(K)) {
    phi_traj[k] <- curr_phi

    em_fit <- run_em_nb_nb_pure_R(
      g               = g,
      g_mus_pert0     = g_mus_pert0,
      phi             = curr_phi,
      pi_guesses      = pi_guesses,
      g_pert_guesses  = g_pert_guesses,
      log_g_factorial = log_g_factorial,
      max_iter        = max_iter,
      min_iter        = min_iter,
      ep_tol          = ep_tol
    )

    pi_traj[k]        <- em_fit$outer_pi
    g_pert_traj[k]    <- em_fit$outer_g_pert
    log_lik_traj[k]   <- em_fit$outer_log_lik
    converged_traj[k] <- em_fit$outer_converged
    init_i_traj[k]    <- em_fit$outer_i

    # Update phi unless this was the final EM call.
    if (k <= n_phi_updates) {
      if (isTRUE(em_fit$outer_converged) &&
          is.finite(em_fit$outer_pi) && is.finite(em_fit$outer_g_pert)) {
        new_phi <- tryCatch(
          est_phi_given_model_pure_R(
            g           = g,
            pi          = em_fit$outer_pi,
            gamma       = em_fit$outer_g_pert,
            g_mus_pert0 = g_mus_pert0
          ),
          error = function(e) NA_real_
        )
        if (!is.na(new_phi) && is.finite(new_phi) && new_phi > 0) {
          curr_phi <- new_phi
        }
      }
      # else: EM didn't converge or update failed -> keep curr_phi
    }
  }

  em_fit$phi_final  <- curr_phi
  em_fit$trajectory <- list(
    phi       = phi_traj,
    pi        = pi_traj,
    g_pert    = g_pert_traj,
    log_lik   = log_lik_traj,
    converged = converged_traj,
    init_i    = init_i_traj
  )
  em_fit
}


# ---- NB size (theta = phi) MLE via uniroot ----------------------------------
# 1D weighted MLE of the NB size parameter holding mu fixed:
#   max_phi  sum_i weights_i * log dnbinom(g_i; mu_i, phi)
# The score is monotone decreasing in phi (NB log-lik is concave in phi),
# so uniroot finds the unique root. Returns scalar phi; on uniroot failure
# returns init_phi unchanged.
nb_size_score <- function(log_phi, g, mu, weights) {
  phi <- exp(log_phi)
  sum(weights *
      (log(phi) - log(phi + mu) + (mu - g) / (mu + phi) +
       digamma(g + phi) - digamma(phi)))
}

update_nb_size_pure_R <- function(g, mu, weights, init_phi,
                                   log_phi_lower = -15, log_phi_upper = 15) {
  res <- tryCatch(
    stats::uniroot(nb_size_score,
                   interval  = c(log_phi_lower, log_phi_upper),
                   extendInt = "downX",
                   g = g, mu = mu, weights = weights),
    error = function(e) NULL
  )
  if (is.null(res)) init_phi else exp(res$root)
}


# ---- Joint (gamma, phi1) M-step via alternating uniroot ---------------------
# Block coordinate ascent: alternate the 1D MLEs for gamma (given phi1) and
# phi1 (given gamma) until either max_iter is hit or the change in
# (gamma, log phi1) falls below tol. Each individual M-step uses uniroot;
# on failure the prev value is retained for that step.
update_gamma_phi1_pure_R <- function(g, g_mus_pert0, Ti1s,
                                      gamma_init, phi1_init,
                                      max_iter = 100L, tol = 1e-6) {
  curr_gamma <- gamma_init
  curr_phi1  <- phi1_init
  k <- 0L
  for (k in seq_len(max_iter)) {
    # gamma update given phi1 (reuse the existing nb-shared score)
    y_plus_phi           <- g + curr_phi1
    offset_minus_log_phi <- log(g_mus_pert0) - log(curr_phi1)
    res_gamma <- tryCatch(
      stats::uniroot(gamma_score,
                     interval             = c(curr_gamma - 2, curr_gamma + 2),
                     extendInt            = "downX",
                     y_plus_phi           = y_plus_phi,
                     offset_minus_log_phi = offset_minus_log_phi,
                     phi                  = curr_phi1,
                     prob_is_1            = Ti1s),
      error = function(e) NULL
    )
    new_gamma <- if (is.null(res_gamma)) curr_gamma else res_gamma$root

    # phi1 update given gamma (1D NB size MLE)
    mu1      <- exp(new_gamma) * g_mus_pert0
    new_phi1 <- update_nb_size_pure_R(g, mu = mu1, weights = Ti1s,
                                       init_phi = curr_phi1)

    delta      <- abs(new_gamma - curr_gamma) + abs(log(new_phi1) - log(curr_phi1))
    curr_gamma <- new_gamma
    curr_phi1  <- new_phi1
    if (delta < tol) break
  }
  list(gamma = curr_gamma, phi1 = curr_phi1, n_iter = k)
}


# ---- NB-NB separate-phi mixture EM ------------------------------------------
# Same outer-loop scaffolding as run_em_nb_nb_pure_R (B starting guesses,
# sceptre's relative-tol convergence test, degeneracy bail-out), but the
# mixture has class-specific overdispersions:
#   nonpert cells:  g_i ~ NB(mean = g_mus_pert0_i,           size = phi0)
#   pert    cells:  g_i ~ NB(mean = exp(gamma)*g_mus_pert0_i, size = phi1)
#
# M-step is ECM in three conditional blocks: phi0, then (gamma, phi1)
# alternating, then pi (closed-form mean(Ti1s), clipped). No label-switch
# flip: with phi0 != phi1, flipping would also require swapping phi0 <-> phi1
# and re-anchoring gamma to the offset, which the conditional M-step can't
# do directly. Pi clipping keeps pi in (1e-8, 1 - 1e-8) for numerical
# stability instead.
#
# Each starting init seeds phi0 = phi1 = phi (the scalar input); the EM
# then separates them.
run_em_nb_separate_pure_R <- function(g, g_mus_pert0, phi,
                                       pi_guesses, g_pert_guesses,
                                       log_g_factorial = NULL,
                                       max_iter = 50L, min_iter = 3L,
                                       ep_tol = 0.5e-4,
                                       max_iter_gamma_phi1 = 100L,
                                       tol_gamma_phi1      = 1e-6) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0)    == n,
            length(phi) == 1L, is.finite(phi), phi > 0)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_
  outer_phi0      <- NA_real_
  outer_phi1      <- NA_real_

  for (i in seq_len(B)) {
    curr_pi      <- min(max(pi_guesses[i], 1e-8), 1 - 1e-8)
    curr_gamma   <- g_pert_guesses[i]
    curr_phi0    <- phi
    curr_phi1    <- phi
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s         <- numeric(n)
    converged    <- FALSE
    iteration    <- 1L

    repeat {
      # E-step + observed log-lik (single dnbinom pass per class)
      mu1          <- exp(curr_gamma) * g_mus_pert0
      log_dnb0     <- stats::dnbinom(g, mu = g_mus_pert0, size = curr_phi0, log = TRUE)
      log_dnb1     <- stats::dnbinom(g, mu = mu1,         size = curr_phi1, log = TRUE)
      log_p0       <- log1p(-curr_pi) + log_dnb0
      log_p1       <- log(curr_pi)    + log_dnb1
      curr_log_lik <- sum(logaddexp(log_p0, log_p1))
      Ti1s         <- stats::plogis(log_p1 - log_p0)

      # bail-out on degenerate posteriors
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # pi update (closed-form, clipped; no flip with separate phi)
      curr_pi <- min(max(mean(Ti1s), 1e-8), 1 - 1e-8)

      # convergence check (sceptre's relative-tol form)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for phi0 (weights = 1 - Ti1s; class mean = g_mus_pert0)
      curr_phi0 <- update_nb_size_pure_R(
        g        = g, mu = g_mus_pert0,
        weights  = 1 - Ti1s,
        init_phi = curr_phi0
      )

      # M-step for (gamma, phi1) (weights = Ti1s; alternating uniroot)
      gp <- update_gamma_phi1_pure_R(
        g           = g,
        g_mus_pert0 = g_mus_pert0,
        Ti1s        = Ti1s,
        gamma_init  = curr_gamma,
        phi1_init   = curr_phi1,
        max_iter    = max_iter_gamma_phi1,
        tol         = tol_gamma_phi1
      )
      curr_gamma <- gp$gamma
      curr_phi1  <- gp$phi1
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_gamma
      outer_phi0      <- curr_phi0
      outer_phi1      <- curr_phi1
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert,
       outer_phi0      = outer_phi0,
       outer_phi1      = outer_phi1)
}


# ---- NB-NB mixture EM with FIXED nonpert size (theta0), free theta1 ---------
# Sits between nb-shared (a single shared size) and nb-separate (both sizes
# free): the nonpert size theta0 is held FIXED at the value supplied by the
# offset model, and only the perturbed size theta1 is estimated.
#   nonpert cells: g_i ~ NB(mean = g_mus_pert0_i,            size = theta0)  [FIXED]
#   pert    cells: g_i ~ NB(mean = exp(gamma)*g_mus_pert0_i, size = theta1)  [free]
# Both the offset (o_i = log g_mus_pert0_i) and theta0 are treated as known.
# theta0 enters via the `phi` argument -- intended to come from an NB offset
# fit's $theta (e.g. fit_baseline_glm_nb_threshold_pure_R) through
# estimate_phi_from_offset_fit_theta. The EM fits (gamma, pi, theta1).
#
# Same outer-loop scaffolding as run_em_nb_separate_pure_R (B starts, sceptre's
# relative-tol convergence test, degeneracy bail-out), but with NO phi0 M-step.
# The (gamma, theta1) M-step maximizes the weighted perturbed-component NB
# log-likelihood -- valid because the nonpert term (theta0, mu0 both fixed) is
# constant in (gamma, theta1). It is solved jointly via box-constrained optim
# with gamma >= 0; the gamma >= 0 constraint pins component 1 as the elevated
# (perturbed) one, so (as in nb-separate) there is no label-switch flip. pi is
# the closed-form clipped mean(Ti1s).
run_em_nb_fixed0_pure_R <- function(g, g_mus_pert0, phi,
                                    pi_guesses, g_pert_guesses,
                                    log_g_factorial = NULL,
                                    max_iter = 50L, min_iter = 3L,
                                    ep_tol = 0.5e-4,
                                    max_gamma = 20, min_theta = 1e-4, max_theta = 1e6,
                                    mstep_maxit = 100L) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B,
            length(g_mus_pert0)    == n,
            length(phi) == 1L, is.finite(phi), phi > 0)

  theta0 <- phi
  mu0    <- g_mus_pert0

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_
  outer_phi0      <- NA_real_
  outer_phi1      <- NA_real_

  for (i in seq_len(B)) {
    curr_pi     <- min(max(pi_guesses[i], 1e-8), 1 - 1e-8)
    curr_gamma  <- max(0, min(g_pert_guesses[i], max_gamma))
    curr_theta1 <- min(max(theta0, min_theta), max_theta)
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s        <- numeric(n)
    converged   <- FALSE
    iteration   <- 1L

    repeat {
      # E-step + observed log-lik (nonpert size fixed at theta0)
      mu1      <- mu0 * exp(curr_gamma)
      log_dnb0 <- stats::dnbinom(g, mu = mu0, size = theta0,      log = TRUE)
      log_dnb1 <- stats::dnbinom(g, mu = mu1, size = curr_theta1, log = TRUE)
      log_p0   <- log1p(-curr_pi) + log_dnb0
      log_p1   <- log(curr_pi)    + log_dnb1
      curr_log_lik <- sum(logaddexp(log_p0, log_p1))
      Ti1s     <- stats::plogis(log_p1 - log_p0)

      # bail-out on degenerate posteriors
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # pi update (closed-form, clipped; no flip with gamma >= 0)
      curr_pi <- min(max(mean(Ti1s), 1e-8), 1 - 1e-8)

      # convergence check (sceptre's relative-tol form)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for (gamma, theta1): maximize the weighted perturbed-component
      # NB log-lik (the nonpert term is constant in these params). Joint
      # box-constrained optim over (gamma, log theta1), gamma >= 0.
      neg_q <- function(par) {
        mu_c <- mu0 * exp(par[1])
        -sum(Ti1s * stats::dnbinom(g, mu = mu_c, size = exp(par[2]), log = TRUE))
      }
      opt <- tryCatch(
        stats::optim(
          par     = c(curr_gamma, log(curr_theta1)),
          fn      = neg_q,
          method  = "L-BFGS-B",
          lower   = c(0,         log(min_theta)),
          upper   = c(max_gamma, log(max_theta)),
          control = list(maxit = mstep_maxit)
        ),
        error = function(e) NULL
      )
      if (is.null(opt)) {
        curr_log_lik <- -Inf
        break
      }
      curr_gamma  <- opt$par[1]
      curr_theta1 <- exp(opt$par[2])
    }

    if (converged && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_gamma
      outer_phi0      <- theta0       # fixed
      outer_phi1      <- curr_theta1  # estimated
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert,
       outer_phi0      = outer_phi0,
       outer_phi1      = outer_phi1)
}


# ---- Poisson-null / NB-pert mixture EM --------------------------------------
# Mixture with a POISSON non-pert component and an NB perturbed component:
#   nonpert cells: g_i ~ Pois(mu0_i)                        [no free parameter]
#   pert    cells: g_i ~ NB(mean = exp(gamma)*mu0_i, size = theta)
# where mu0_i = exp(o_i) is the fixed offset mean. The EM fits (gamma, pi,
# theta) -- like nb-fixed0, but the non-pert component is exactly Poisson
# rather than NB with a fixed theta0, so there is no phi/theta0 input at all.
#
# Same outer-loop scaffolding as run_em_nb_fixed0_pure_R: the (gamma, theta)
# M-step maximizes the weighted perturbed-component NB log-likelihood (valid
# because the Poisson non-pert term, mu0 fixed, is constant in gamma/theta),
# solved jointly via box-constrained optim with gamma >= 0 (which pins the
# perturbed component as the elevated one, so no label-switch flip). pi is the
# closed-form clipped mean(Ti1s). The non-pert Poisson log-density never changes
# across iterations, so it is precomputed once. outer_phi0 is reported as Inf
# (Poisson = NB with theta -> Inf); outer_phi1 is the estimated perturbed theta.
#
# Minority guard (pi_max): the Poisson non-pert component is rigid (variance =
# mean), so a perturbed NB with gamma -> 0 and tiny theta can collapse into a
# fat-tailed component that swallows the MAJORITY of cells at higher observed
# log-likelihood than the true minority-perturbation solution. To reject that
# degenerate optimum, only converged starts with pi <= pi_max (default 0.5) are
# eligible for selection -- the perturbation fraction for a guide is essentially
# never a majority. If no start qualifies, outer_converged stays FALSE and
# outer_Ti1s is all zeros (fail-safe: no assignments, rather than ~all cells).
run_em_pois_nb_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                  log_g_factorial = NULL,
                                  max_iter = 50L, min_iter = 3L,
                                  ep_tol = 0.5e-4,
                                  max_gamma = 20, min_theta = 1e-4, max_theta = 1e6,
                                  mstep_maxit = 100L, pi_max = 0.5) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B, length(g_mus_pert0) == n)

  mu0 <- g_mus_pert0
  # Poisson non-pert log-density: fixed (mu0 known), so compute once.
  log_f0 <- stats::dpois(g, lambda = mu0, log = TRUE)

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_
  outer_phi0      <- NA_real_
  outer_phi1      <- NA_real_

  for (i in seq_len(B)) {
    curr_pi    <- min(max(pi_guesses[i], 1e-8), 1 - 1e-8)
    curr_gamma <- max(0, min(g_pert_guesses[i], max_gamma))
    curr_theta <- min(max(1, min_theta), max_theta)
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s       <- numeric(n)
    converged  <- FALSE
    iteration  <- 1L

    repeat {
      # E-step + observed log-lik
      mu1    <- mu0 * exp(curr_gamma)
      log_f1 <- stats::dnbinom(g, mu = mu1, size = curr_theta, log = TRUE)
      log_p0 <- log1p(-curr_pi) + log_f0
      log_p1 <- log(curr_pi)    + log_f1
      curr_log_lik <- sum(logaddexp(log_p0, log_p1))
      Ti1s   <- stats::plogis(log_p1 - log_p0)

      # bail-out on degenerate posteriors
      if (all(Ti1s <= 1e-100) || any(!is.finite(Ti1s))) {
        curr_log_lik <- -Inf
        break
      }

      # pi update (closed-form, clipped; no flip with gamma >= 0)
      curr_pi <- min(max(mean(Ti1s), 1e-8), 1 - 1e-8)

      # convergence check (sceptre's relative-tol form)
      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for (gamma, theta): maximize the weighted perturbed-component
      # NB log-lik (the Poisson non-pert term is constant in these params).
      neg_q <- function(par) {
        mu_c <- mu0 * exp(par[1])
        -sum(Ti1s * stats::dnbinom(g, mu = mu_c, size = exp(par[2]), log = TRUE))
      }
      opt <- tryCatch(
        stats::optim(
          par     = c(curr_gamma, log(curr_theta)),
          fn      = neg_q,
          method  = "L-BFGS-B",
          lower   = c(0,         log(min_theta)),
          upper   = c(max_gamma, log(max_theta)),
          control = list(maxit = mstep_maxit)
        ),
        error = function(e) NULL
      )
      if (is.null(opt)) {
        curr_log_lik <- -Inf
        break
      }
      curr_gamma <- opt$par[1]
      curr_theta <- exp(opt$par[2])
    }

    # Minority guard: reject the degenerate "NB swallows the majority" optimum.
    if (converged && curr_pi <= pi_max && curr_log_lik > outer_log_lik) {
      outer_Ti1s      <- Ti1s
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_gamma
      outer_phi0      <- Inf          # Poisson non-pert = NB(theta -> Inf)
      outer_phi1      <- curr_theta   # estimated perturbed size
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert,
       outer_phi0      = outer_phi0,
       outer_phi1      = outer_phi1)
}


# ---- Poisson-null / NB-pert mixture EM, zeros forced non-pert ("pois0-nb") --
# Same model as run_em_pois_nb_pure_R -- (1-pi) Pois(mu0) + pi NB(mu0 e^gamma,
# theta) -- but cells with g == 0 are forced to the non-pert component: a 0-UMI
# cell cannot carry the guide. Because zeros never enter the perturbed
# component, the EM only ever touches the POSITIVE cells with the (expensive)
# dnbinom/optim calls; the zeros enter only through closed-form O(1) terms. For
# the sparse guide-count regime (often >99% zeros) this makes the per-iteration
# cost scale with n_pos rather than n.
#
# Bookkeeping that the hard-zero constraint makes exact and cheap:
#   * pi M-step: pi = sum(Ti1s_pos) / n  (zeros contribute weight 0). This is
#     the exact MLE of pi for the constrained model -- the n_zero "votes" for
#     non-pert enter the pi-score only through the n_zero * log(1 - pi) term,
#     whose derivative gives pi = sum(w)/n.
#   * observed log-lik: sum over positives of the 2-component mixture, plus the
#     zeros' (definitely-non-pert) contribution n_zero * log(1 - pi) + sum(log
#     dpois(0; mu0_zero)) = n_zero * log(1 - pi) - sum(mu0_zero). The second
#     piece is a constant (precomputed); only n_zero * log(1 - pi) updates, so
#     no per-iteration pass over the zeros is needed.
#
# Same (gamma, theta) joint box-constrained optim M-step (gamma >= 0, no flip)
# and the pi_max minority guard as run_em_pois_nb_pure_R. outer_phi0 = Inf
# (Poisson non-pert); outer_phi1 = estimated perturbed theta.
run_em_pois0_nb_pure_R <- function(g, g_mus_pert0, pi_guesses, g_pert_guesses,
                                   log_g_factorial = NULL,
                                   max_iter = 50L, min_iter = 3L,
                                   ep_tol = 0.5e-4,
                                   max_gamma = 20, min_theta = 1e-4, max_theta = 1e6,
                                   mstep_maxit = 100L, pi_max = 0.5) {
  n <- length(g)
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B, length(g_mus_pert0) == n)

  mu0      <- g_mus_pert0
  is_pos   <- g > 0
  n_zero   <- sum(!is_pos)
  y_pos    <- g[is_pos]
  mu0_pos  <- mu0[is_pos]
  # Fixed Poisson non-pert log-density on positives; zeros' non-pert log-density
  # dpois(0; mu0) = -mu0 enters the log-lik only as the constant below.
  log_f0_pos     <- stats::dpois(y_pos, lambda = mu0_pos, log = TRUE)
  loglik_zero_const <- -sum(mu0[!is_pos])

  outer_Ti1s      <- numeric(n)
  outer_i         <- 0L
  outer_log_lik   <- -Inf
  outer_converged <- FALSE
  outer_pi        <- NA_real_
  outer_g_pert    <- NA_real_
  outer_phi0      <- NA_real_
  outer_phi1      <- NA_real_

  for (i in seq_len(B)) {
    curr_pi    <- min(max(pi_guesses[i], 1e-8), 1 - 1e-8)
    curr_gamma <- max(0, min(g_pert_guesses[i], max_gamma))
    curr_theta <- min(max(1, min_theta), max_theta)
    prev_log_lik <- -Inf
    curr_log_lik <- -Inf
    Ti1s_pos   <- numeric(length(y_pos))
    converged  <- FALSE
    iteration  <- 1L

    repeat {
      # E-step over POSITIVE cells only (zeros are forced non-pert).
      mu1_pos    <- mu0_pos * exp(curr_gamma)
      log_f1_pos <- stats::dnbinom(y_pos, mu = mu1_pos, size = curr_theta, log = TRUE)
      log_p0_pos <- log1p(-curr_pi) + log_f0_pos
      log_p1_pos <- log(curr_pi)    + log_f1_pos
      Ti1s_pos   <- stats::plogis(log_p1_pos - log_p0_pos)

      # Observed log-lik: positives' mixture + zeros' (non-pert) contribution.
      curr_log_lik <- sum(logaddexp(log_p0_pos, log_p1_pos)) +
                      n_zero * log1p(-curr_pi) + loglik_zero_const

      # bail-out on degenerate posteriors
      if (length(Ti1s_pos) == 0L || all(Ti1s_pos <= 1e-100) ||
          any(!is.finite(Ti1s_pos))) {
        curr_log_lik <- -Inf
        break
      }

      # pi over ALL n cells (zeros contribute 0): exact MLE = sum(w_pos) / n
      curr_pi <- min(max(sum(Ti1s_pos) / n, 1e-8), 1 - 1e-8)

      tol <- compute_em_tolerance_pure_R(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_iter) {
        converged <- TRUE
        break
      }
      prev_log_lik <- curr_log_lik
      iteration    <- iteration + 1L
      if (iteration >= max_iter) break

      # M-step for (gamma, theta): POSITIVE cells only (zeros have weight 0).
      neg_q <- function(par) {
        mu_c <- mu0_pos * exp(par[1])
        -sum(Ti1s_pos * stats::dnbinom(y_pos, mu = mu_c, size = exp(par[2]), log = TRUE))
      }
      opt <- tryCatch(
        stats::optim(
          par     = c(curr_gamma, log(curr_theta)),
          fn      = neg_q,
          method  = "L-BFGS-B",
          lower   = c(0,         log(min_theta)),
          upper   = c(max_gamma, log(max_theta)),
          control = list(maxit = mstep_maxit)
        ),
        error = function(e) NULL
      )
      if (is.null(opt)) {
        curr_log_lik <- -Inf
        break
      }
      curr_gamma <- opt$par[1]
      curr_theta <- exp(opt$par[2])
    }

    # Minority guard: reject the degenerate "NB swallows the majority" optimum.
    if (converged && curr_pi <= pi_max && curr_log_lik > outer_log_lik) {
      Ti1s_full          <- numeric(n)
      Ti1s_full[is_pos]  <- Ti1s_pos     # zeros stay 0 (forced non-pert)
      outer_Ti1s      <- Ti1s_full
      outer_i         <- i
      outer_log_lik   <- curr_log_lik
      outer_converged <- TRUE
      outer_pi        <- curr_pi
      outer_g_pert    <- curr_gamma
      outer_phi0      <- Inf          # Poisson non-pert = NB(theta -> Inf)
      outer_phi1      <- curr_theta   # estimated perturbed size
    }
  }

  list(outer_Ti1s      = outer_Ti1s,
       outer_i         = outer_i,
       outer_converged = outer_converged,
       outer_log_lik   = outer_log_lik,
       outer_pi        = outer_pi,
       outer_g_pert    = outer_g_pert,
       outer_phi0      = outer_phi0,
       outer_phi1      = outer_phi1)
}


# ---- Per-guide assignment (mirrors obtain_em_assignments) --------------------
# Returns the assignments + a small per-guide summary by default. Set
# keep_fits = TRUE to also retain the full offset-model fit object and full
# EM result (each carries length-N vectors per guide -- expensive at scale).
#
# Per-guide output schema (same fields regardless of family / cutoff path):
#   assignments           : integer vector of cell indices
#   n_nonzero             : sum(g >= 1)
#   n_assigned            : length(assignments)
#   em_converged          : logical or NA (NA when EM was skipped below cutoff)
#   em_log_lik            : numeric or NA
#   em_init_i             : integer or NA  -- winning starting-guess index
#   em_pi                 : numeric or NA  -- converged mixing proportion
#   em_g_pert             : numeric or NA  -- converged perturbation effect
#                            (additive for pois w/ fix_curr_g_pert_bug=FALSE,
#                             log-fold-change otherwise; see run_meta)
#   em_phi_pert           : numeric or NA  -- NB only. For nb-shared: equals
#                            em_phi_nonpert (one shared phi). For nb-separate:
#                            the class-specific size for perturbed cells (= phi1).
#   em_phi_nonpert        : numeric or NA  -- NB only. For nb-separate: phi0.
#   em_phi_source         : character or NA  -- one of "estimated" (estimate_phi_fn
#                            succeeded), "fallback" (estimate_phi_fn errored, used
#                            scalar phi), or "fixed" (no estimate_phi_fn provided).
#                            Describes the *initial* phi only; refinement across
#                            phi updates is in em_trajectory. NA for pois / below cutoff.
#   em_trajectory         : list of length-(n_phi_updates+1) parallel vectors
#                            recording each EM pass:
#                              $phi       : phi used as input. For nb-separate
#                                           (length 1), this is the seed value
#                                           used to initialize phi0 = phi1.
#                              $pi        : output pi (NA if pass didn't converge)
#                              $g_pert    : output g_pert
#                              $log_lik   : output log-lik
#                              $converged : output converged flag
#                              $init_i    : output winning starting-guess index
#                            For pois / below-cutoff guides, all entries are NA.
#                            For nb-separate the trajectory is length 1 (n_phi_updates
#                            must be 0); the final phi0/phi1 live in em_phi_nonpert/
#                            em_phi_pert and are not in the trajectory.
#                            em_pi / em_g_pert / em_log_lik / em_converged / em_init_i
#                            and em_phi_pert / em_phi_nonpert are the LAST trajectory
#                            entry; the trajectory shows how they were arrived at.
#   prob_quantiles        : length-7 named numeric (0/1/10/50/90/99/100% of Ti1s)
#   n_above_prob_thresh   : sum(Ti1s >= probability_threshold)
#   offset_model_summary  : offset_fit$offset_model_summary, or NULL below cutoff
#   elapsed_sec           : wall time for this guide
#
# Below the n_nonzero cutoff, EM is skipped: assignments come from
# `which(g >= backup_threshold)` and all em_* / prob_* fields are NA.
assign_one_grna_pure_R <- function(g, covariate_matrix,
                                   pi_guesses, g_pert_guesses,
                                   n_nonzero_cells_cutoff = 10L,
                                   backup_threshold       = 5L,
                                   probability_threshold  = 0.8,
                                   keep_fits              = FALSE,
                                   fix_curr_g_pert_bug    = FALSE,
                                   use_log_gamma          = FALSE,
                                   family                 = c("pois", "pois-additive", "pois-additive-nonzero", "pois-wq", "nb-shared", "nb-separate", "nb-fixed0", "pois-nb", "pois0-nb", "poissum", "nbsum"),
                                   phi                    = NULL,
                                   estimate_phi_fn        = NULL,
                                   n_phi_updates          = 0L,
                                   gamma_update_prob      = 0.5,
                                   offset_model_fit_fn    = fit_baseline_glm_pure_R) {
  family <- match.arg(family)
  if (family %in% c("nb-shared", "nb-separate", "nb-fixed0")) {
    if (is.null(estimate_phi_fn) &&
        (is.null(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0)) {
      stop("`family = \"", family, "\"` requires either `estimate_phi_fn` or ",
           "a single positive finite `phi` (NB size parameter).")
    }
  }
  stopifnot(length(n_phi_updates) == 1L, is.finite(n_phi_updates),
            n_phi_updates >= 0L)
  n_phi_updates <- as.integer(n_phi_updates)
  if (family %in% c("nb-separate", "nb-fixed0", "pois-nb", "pois0-nb", "nbsum") && n_phi_updates > 0L) {
    stop("`n_phi_updates` must be 0 for `family = \"", family, "\"` ",
         "(theta1 is updated inside each Q-step; theta0 is fixed/estimated jointly).")
  }

  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)

  # defaults for fields that may not get set (below-cutoff path)
  offset_model_fit     <- NULL
  offset_model_summary <- NULL
  em_fit               <- NULL
  em_converged         <- NA
  em_log_lik           <- NA_real_
  em_init_i            <- NA_integer_
  em_pi                <- NA_real_
  em_g_pert            <- NA_real_
  em_phi_pert          <- NA_real_
  em_phi_nonpert       <- NA_real_
  em_phi_source        <- NA_character_
  em_trajectory        <- list(
    phi       = rep(NA_real_,    n_phi_updates + 1L),
    pi        = rep(NA_real_,    n_phi_updates + 1L),
    g_pert    = rep(NA_real_,    n_phi_updates + 1L),
    log_lik   = rep(NA_real_,    n_phi_updates + 1L),
    converged = rep(NA,          n_phi_updates + 1L),
    init_i    = rep(NA_integer_, n_phi_updates + 1L)
  )
  prob_quantile_probs  <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles       <- stats::setNames(
    rep(NA_real_, length(prob_quantile_probs)),
    paste0(prob_quantile_probs * 100, "%")
  )
  n_above_prob_thresh  <- NA_integer_

  if (n_nonzero >= n_nonzero_cells_cutoff) {
    offset_model_fit     <- offset_model_fit_fn(g, covariate_matrix)
    offset_model_summary <- offset_model_fit$offset_model_summary
    g_mus_pert0          <- offset_model_fit$fitted.values
    log_g_factorial      <- lgamma(g + 1)

    # For NB: phi is either estimated per-guide from the offset fit, or the
    # caller-supplied scalar. For nb-shared this is THE shared phi; for
    # nb-separate it's the seed value used to initialize both phi0 and phi1.
    # If estimation errors / returns a non-finite value, fall back to the
    # caller's `phi` scalar and record that via em_phi_source.
    phi_used <- NA_real_
    if (family %in% c("nb-shared", "nb-separate", "nb-fixed0")) {
      if (!is.null(estimate_phi_fn)) {
        phi_estimate <- tryCatch(
          estimate_phi_fn(g, offset_model_fit),
          error = function(e) NULL
        )
        if (!is.null(phi_estimate) && length(phi_estimate) == 1L &&
            is.finite(phi_estimate) && phi_estimate > 0) {
          phi_used      <- phi_estimate
          em_phi_source <- "estimated"
        } else {
          phi_used      <- phi
          em_phi_source <- "fallback"
        }
      } else {
        phi_used      <- phi
        em_phi_source <- "fixed"
      }
      if (is.null(phi_used) || !is.finite(phi_used) || phi_used <= 0) {
        stop("NB EM needs a positive finite phi but got ", phi_used,
             " (estimate_phi_fn failed/returned bad value and no usable ",
             "scalar phi was supplied).")
      }
    }

    em_fit <- if (family == "pois") {
      run_em_pois_pois_pure_R(
        g                   = g,
        g_mus_pert0         = g_mus_pert0,
        pi_guesses          = pi_guesses,
        g_pert_guesses      = g_pert_guesses,
        log_g_factorial     = log_g_factorial,
        fix_curr_g_pert_bug = fix_curr_g_pert_bug
      )
    } else if (family == "pois-additive") {
      run_em_pois_additive_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "pois-additive-nonzero") {
      run_em_pois_additive_nonzero_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "pois-wq") {
      run_em_pois_pois_wquantile_pure_R(
        g                 = g,
        g_mus_pert0       = g_mus_pert0,
        pi_guesses        = pi_guesses,
        g_pert_guesses    = g_pert_guesses,
        log_g_factorial   = log_g_factorial,
        gamma_update_prob = gamma_update_prob
      )
    } else if (family == "nb-shared") {
      # Always route nb-shared through the update_phi wrapper. With
      # n_phi_updates = 0 it runs one EM pass (same as run_em_nb_nb_pure_R)
      # but returns the length-1 trajectory; with > 0 it alternates EM and
      # phi updates.
      run_em_nb_nb_update_phi_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        phi             = phi_used,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        n_phi_updates   = n_phi_updates,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "nb-separate") {
      run_em_nb_separate_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        phi             = phi_used,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "nb-fixed0") {
      run_em_nb_fixed0_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        phi             = phi_used,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "pois-nb") {
      run_em_pois_nb_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "pois0-nb") {
      run_em_pois0_nb_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
      )
    } else if (family == "poissum") {
      # Sum-of-Poissons additive-signal mixture (fit_pois_additive_signal_em in
      # em_variants.R). Single EM fit (no multi-start), offset = log(mu0).
      # use_log_gamma routes the log-scale classification. Adapt its richer
      # return to the outer_* contract the dispatch below consumes; assignments
      # are driven by posterior_classification, which honours use_log_gamma.
      fit <- fit_pois_additive_signal_em(
        y             = g,
        offset        = log(g_mus_pert0),
        use_log_gamma = use_log_gamma
      )
      list(
        outer_Ti1s      = fit$posterior_classification,
        outer_i         = 1L,
        outer_converged = fit$converged,
        outer_log_lik   = fit$loglik,
        outer_pi        = fit$pi,
        outer_g_pert    = fit$gamma
      )
    } else {  # "nbsum"
      # Additive NB-signal Poisson mixture (fit_pois_additive_nb_signal_em in
      # em_variants.R): nonpert Pois(mu0), pert Pois(mu0) + NB(exp(gamma), theta).
      # Single EM fit, offset = log(mu0); use_log_gamma routes the log-scale
      # classification. theta is surfaced as the perturbed NB size (outer_phi1),
      # with outer_phi0 = Inf for the Poisson baseline, mirroring pois-nb.
      fit <- fit_pois_additive_nb_signal_em(
        y             = g,
        offset        = log(g_mus_pert0),
        use_log_gamma = use_log_gamma
      )
      list(
        outer_Ti1s      = fit$posterior_classification,
        outer_i         = 1L,
        outer_converged = fit$converged,
        outer_log_lik   = fit$loglik,
        outer_pi        = fit$pi,
        outer_g_pert    = fit$gamma,
        outer_phi1      = fit$theta,
        outer_phi0      = Inf
      )
    }

    # Always trust the EM result (even if it didn't converge); the user
    # diagnoses flaky fits via em_converged + em_log_lik.
    assignments         <- which(em_fit$outer_Ti1s >= probability_threshold)
    prob_quantiles      <- stats::quantile(em_fit$outer_Ti1s, probs = prob_quantile_probs)
    n_above_prob_thresh <- length(assignments)

    em_converged   <- em_fit$outer_converged
    em_log_lik     <- em_fit$outer_log_lik
    em_init_i      <- em_fit$outer_i
    em_pi          <- em_fit$outer_pi
    em_g_pert      <- em_fit$outer_g_pert
    if (family == "nb-shared") {
      # phi_final reflects all phi updates (= phi_used if n_phi_updates == 0).
      em_phi_pert    <- em_fit$phi_final
      em_phi_nonpert <- em_fit$phi_final
      em_trajectory  <- em_fit$trajectory
      # em_phi_source already set above when phi_used was determined
    } else if (family %in% c("nb-separate", "nb-fixed0", "pois-nb", "pois0-nb", "nbsum")) {
      # nb-separate:      phi0/phi1 both estimated inside each Q-step.
      # nb-fixed0:        phi0 (= phi_used) held fixed; only phi1 estimated.
      # pois-nb/pois0-nb: phi0 = Inf (Poisson non-pert); only phi1 estimated.
      # nbsum:            phi0 = Inf (Poisson baseline); phi1 = theta (NB signal).
      em_phi_pert    <- em_fit$outer_phi1
      em_phi_nonpert <- em_fit$outer_phi0
      # Length-1 trajectory: one EM pass, no outer phi-update loop.
      # $phi records the seed/fixed value (initializes phi0/phi1 for
      # nb-separate; is the fixed phi0 for nb-fixed0).
      em_trajectory$phi[1L]       <- phi_used
      em_trajectory$pi[1L]        <- em_fit$outer_pi
      em_trajectory$g_pert[1L]    <- em_fit$outer_g_pert
      em_trajectory$log_lik[1L]   <- em_fit$outer_log_lik
      em_trajectory$converged[1L] <- em_fit$outer_converged
      em_trajectory$init_i[1L]    <- em_fit$outer_i
    }
  } else {
    assignments <- which(g >= backup_threshold)
  }

  elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  out <- list(
    assignments          = assignments,
    n_nonzero            = n_nonzero,
    n_assigned           = length(assignments),
    em_converged         = em_converged,
    em_log_lik           = em_log_lik,
    em_init_i            = em_init_i,
    em_pi                = em_pi,
    em_g_pert            = em_g_pert,
    em_phi_pert          = em_phi_pert,
    em_phi_nonpert       = em_phi_nonpert,
    em_phi_source        = em_phi_source,
    em_trajectory        = em_trajectory,
    prob_quantiles       = prob_quantiles,
    n_above_prob_thresh  = n_above_prob_thresh,
    offset_model_summary = offset_model_summary,
    elapsed_sec          = elapsed_sec
  )
  if (keep_fits) {
    out$offset_model_fit <- offset_model_fit
    out$em_fit           <- em_fit
  }
  out
}


# ---- Driver ------------------------------------------------------------------
# Does the one-time setup (model matrix, starting guesses) and then loops
# over guides, optionally in parallel via parLapplyLB.
sceptre_assign_pure_R <- function(grna_matrix,
                                  covariate_data_frame,
                                  formula_object,
                                  grna_ids               = rownames(grna_matrix),
                                  n_em_rep               = 5L,
                                  pi_guess_range         = c(1e-5, 0.1),
                                  g_pert_guess_range     = log(c(10, 5000)),
                                  n_nonzero_cells_cutoff = 10L,
                                  backup_threshold       = 5L,
                                  probability_threshold  = 0.8,
                                  seed                   = 4L,
                                  cl                     = NULL,
                                  keep_fits              = FALSE,
                                  fix_curr_g_pert_bug    = FALSE,
                                  use_log_gamma          = FALSE,
                                  family                 = c("pois", "pois-additive", "pois-additive-nonzero", "pois-wq", "nb-shared", "nb-separate", "nb-fixed0", "pois-nb", "pois0-nb", "poissum", "nbsum"),
                                  phi                    = NULL,
                                  estimate_phi_fn        = NULL,
                                  n_phi_updates          = 0L,
                                  gamma_update_prob      = 0.5,
                                  offset_model_fit_fn    = fit_baseline_glm_pure_R) {
  family <- match.arg(family)
  if (family %in% c("nb-shared", "nb-separate", "nb-fixed0")) {
    if (is.null(estimate_phi_fn) &&
        (is.null(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0)) {
      stop("`family = \"", family, "\"` requires either `estimate_phi_fn` or ",
           "a single positive finite `phi` (NB size parameter).")
    }
  }
  stopifnot(length(n_phi_updates) == 1L, is.finite(n_phi_updates),
            n_phi_updates >= 0L)
  n_phi_updates <- as.integer(n_phi_updates)
  if (family %in% c("nb-separate", "nb-fixed0", "pois-nb", "pois0-nb", "nbsum") && n_phi_updates > 0L) {
    stop("`n_phi_updates` must be 0 for `family = \"", family, "\"` ",
         "(theta1 is updated inside each Q-step; theta0 is fixed/estimated jointly).")
  }

  # 0. force all promises that will be captured by the worker closure, so
  # they're shipped to PSOCK workers as values rather than as unevaluated
  # expressions that try to look up symbols (e.g. `scep_sims`) on the worker.
  force(grna_matrix)
  force(n_nonzero_cells_cutoff); force(backup_threshold); force(probability_threshold)
  force(keep_fits); force(fix_curr_g_pert_bug); force(use_log_gamma)
  force(phi); force(estimate_phi_fn)
  force(n_phi_updates); force(gamma_update_prob); force(offset_model_fit_fn)

  # 1. design matrix
  covariate_matrix <- stats::model.matrix(object = formula_object,
                                          data   = covariate_data_frame)
  if (nrow(covariate_matrix) != nrow(covariate_data_frame)) {
    stop("Rows lost from `covariate_data_frame` after applying `formula_object` ",
         "(likely due to NAs).")
  }

  # 2. starting guesses; shared across all guides, identical to sceptre's
  # get_random_starting_guesses(set.seed(4)) behavior. Save/restore the global
  # RNG state so callers' streams aren't perturbed.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  pi_guesses     <- stats::runif(n_em_rep, pi_guess_range[1],     pi_guess_range[2])
  g_pert_guesses <- stats::runif(n_em_rep, g_pert_guess_range[1], g_pert_guess_range[2])
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }

  # 3. per-guide worker (closure captures everything it needs)
  worker <- function(grna_id) {
    g <- as.numeric(grna_matrix[grna_id, ])
    assign_one_grna_pure_R(
      g                      = g,
      covariate_matrix       = covariate_matrix,
      pi_guesses             = pi_guesses,
      g_pert_guesses         = g_pert_guesses,
      n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
      backup_threshold       = backup_threshold,
      probability_threshold  = probability_threshold,
      keep_fits              = keep_fits,
      fix_curr_g_pert_bug    = fix_curr_g_pert_bug,
      use_log_gamma          = use_log_gamma,
      family                 = family,
      phi                    = phi,
      estimate_phi_fn        = estimate_phi_fn,
      n_phi_updates          = n_phi_updates,
      gamma_update_prob      = gamma_update_prob,
      offset_model_fit_fn    = offset_model_fit_fn
    )
  }

  # 4. dispatch
  per_guide <- if (is.null(cl)) {
    lapply(grna_ids, worker)
  } else {
    parallel::parLapplyLB(cl, grna_ids, worker)
  }
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  # 5. run-level metadata: parameters that are constant across guides for this
  # run. The fit fn itself is saved here -- its closure captures variant
  # hyperparameters (e.g. trim_frac), and an optional attr(., "spec") can carry
  # a human-readable label.
  run_meta <- list(
    family                 = family,
    offset_model_fit_fn    = offset_model_fit_fn,
    estimate_phi_fn        = estimate_phi_fn,
    phi                    = phi,
    n_phi_updates          = n_phi_updates,
    gamma_update_prob      = gamma_update_prob,
    fix_curr_g_pert_bug    = fix_curr_g_pert_bug,
    use_log_gamma          = use_log_gamma,
    n_em_rep               = n_em_rep,
    pi_guess_range         = pi_guess_range,
    g_pert_guess_range     = g_pert_guess_range,
    n_nonzero_cells_cutoff = n_nonzero_cells_cutoff,
    backup_threshold       = backup_threshold,
    probability_threshold  = probability_threshold,
    seed                   = seed
  )

  list(per_guide = per_guide, run_meta = run_meta)
}
