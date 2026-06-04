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
# Cantoni-Ronchetti Mqle estimator from `robustbase::glmrob`.
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
  # Fisher information well-conditioned. Mqle inverts X' W X directly, which
  # is fragile when columns have wildly different scales (a common situation
  # with sceptre's covariates: `(Intercept)` is 1, `log(response_n_umis)` is
  # 8-12, `grna_n_umis` can be 0-10000). Standardization is invariant for
  # fitted.values, so downstream code is unaffected.
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
      method = "Mqle"
    )
  )
  fit$offset_model_summary <- list(
    coefficients = fit$coefficients,
    deviance     = fit$deviance,
    iter         = fit$iter,
    converged    = fit$converged
  )
  fit
}


# ---- Trimmed robust Poisson GLM baseline fit --------------------------------
# Robust analogue of fit_baseline_glm_trimmed_pure_R: drop the top `trim_frac`
# of cells by g, then run glmrob on the remainder, then evaluate fitted means
# on ALL cells. Standardization stats are computed from the kept rows so the
# Mqle solver sees a well-conditioned design (same rationale as the untrimmed
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
      method = "Mqle"
    )
  )

  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  eta_all <- as.numeric(X %*% coef)
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


# ---- Trimmed Poisson GLM baseline fit ---------------------------------------
# Alternative robust baseline that's much more predictable than glmrob: drop
# the top `trim_frac` of cells by g, fit stats::glm.fit on the remainder, then
# evaluate fitted means on ALL cells (the trimmed cells receive a baseline
# prediction extrapolated from the bulk, not influenced by themselves).
#
# This avoids the multiple-solutions pathology of Mqle: the Poisson MLE on a
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
                                   family                 = c("pois", "nb-shared", "nb-separate"),
                                   phi                    = NULL,
                                   estimate_phi_fn        = NULL,
                                   n_phi_updates          = 0L,
                                   offset_model_fit_fn    = fit_baseline_glm_pure_R) {
  family <- match.arg(family)
  if (family %in% c("nb-shared", "nb-separate")) {
    if (is.null(estimate_phi_fn) &&
        (is.null(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0)) {
      stop("`family = \"", family, "\"` requires either `estimate_phi_fn` or ",
           "a single positive finite `phi` (NB size parameter).")
    }
  }
  stopifnot(length(n_phi_updates) == 1L, is.finite(n_phi_updates),
            n_phi_updates >= 0L)
  n_phi_updates <- as.integer(n_phi_updates)
  if (family == "nb-separate" && n_phi_updates > 0L) {
    stop("`n_phi_updates` must be 0 for `family = \"nb-separate\"` ",
         "(phi0/phi1 are updated inside each Q-step).")
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
    if (family %in% c("nb-shared", "nb-separate")) {
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
    } else {  # "nb-separate"
      run_em_nb_separate_pure_R(
        g               = g,
        g_mus_pert0     = g_mus_pert0,
        phi             = phi_used,
        pi_guesses      = pi_guesses,
        g_pert_guesses  = g_pert_guesses,
        log_g_factorial = log_g_factorial
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
    } else if (family == "nb-separate") {
      # phi0/phi1 are estimated jointly with (pi, gamma) inside each Q-step.
      em_phi_pert    <- em_fit$outer_phi1
      em_phi_nonpert <- em_fit$outer_phi0
      # Length-1 trajectory: one EM pass, no outer phi-update loop.
      # $phi records the seed value used to initialize both phi0 and phi1.
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
                                  family                 = c("pois", "nb-shared", "nb-separate"),
                                  phi                    = NULL,
                                  estimate_phi_fn        = NULL,
                                  n_phi_updates          = 0L,
                                  offset_model_fit_fn    = fit_baseline_glm_pure_R) {
  family <- match.arg(family)
  if (family %in% c("nb-shared", "nb-separate")) {
    if (is.null(estimate_phi_fn) &&
        (is.null(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0)) {
      stop("`family = \"", family, "\"` requires either `estimate_phi_fn` or ",
           "a single positive finite `phi` (NB size parameter).")
    }
  }
  stopifnot(length(n_phi_updates) == 1L, is.finite(n_phi_updates),
            n_phi_updates >= 0L)
  n_phi_updates <- as.integer(n_phi_updates)
  if (family == "nb-separate" && n_phi_updates > 0L) {
    stop("`n_phi_updates` must be 0 for `family = \"nb-separate\"` ",
         "(phi0/phi1 are updated inside each Q-step).")
  }

  # 0. force all promises that will be captured by the worker closure, so
  # they're shipped to PSOCK workers as values rather than as unevaluated
  # expressions that try to look up symbols (e.g. `scep_sims`) on the worker.
  force(grna_matrix)
  force(n_nonzero_cells_cutoff); force(backup_threshold); force(probability_threshold)
  force(keep_fits); force(fix_curr_g_pert_bug); force(phi); force(estimate_phi_fn)
  force(n_phi_updates); force(offset_model_fit_fn)

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
      family                 = family,
      phi                    = phi,
      estimate_phi_fn        = estimate_phi_fn,
      n_phi_updates          = n_phi_updates,
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
    fix_curr_g_pert_bug    = fix_curr_g_pert_bug,
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
