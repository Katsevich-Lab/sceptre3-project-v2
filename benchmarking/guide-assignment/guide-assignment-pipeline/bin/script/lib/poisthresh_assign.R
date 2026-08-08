# Fixed-shift Poisson mixture with a glmpois offset (poisthresh<K>).
#
# Stagewise per guide:
#   1. glmpois offset model (Poisson GLM, log link) on circularity-corrected
#      covariates -- this guide's own contribution is removed from the totals:
#        other_umis_i = grna_n_umis_i    - y_i           (subtract own UMIs)
#        other_nnz_i  = grna_n_nonzero_i - 1{y_i > 0}     (subtract own presence)
#        log(mu0_i) = b0 + b1 log1p(other_umis_i) + b2 log1p(other_nnz_i)
#      (log1p matches the project's canonical log(x + 1) covariate transform.)
#   2. Two-component Poisson mixture with a HARD-CODED additive shift K (not
#      estimated):
#        y_i ~ (1 - pi) Pois(mu0_i) + pi Pois(mu0_i + K)
#      Only pi is fit, by EM. Component 1 (the +K bump) is always the high one,
#      so there is no label-switch ambiguity. The per-cell mixture log-lik is
#      concave in pi, so a single EM run reaches the global MLE.
#
# Assign cell i perturbed iff posterior P(pert_i | y_i) >= probability_threshold.

# logaddexp(a, b) = log(exp(a) + exp(b)), stable.
poisthresh_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  m + log1p(exp(-abs(a - b)))
}

# glmpois offset fit. X is the (n x p) design (intercept + transformed covars).
# Falls back to an intercept-only marginal mean if the GLM fails to converge.
poisthresh_fit_offset <- function(g, X) {
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(y = g, x = X, family = stats::poisson())),
    error = function(e) NULL
  )
  ok <- !is.null(fit) && isTRUE(fit$converged) && all(is.finite(fit$coefficients))
  if (ok) {
    return(list(mu0 = pmax(fit$fitted.values, 1e-300),
                coef = fit$coefficients, deviance = fit$deviance, converged = TRUE))
  }
  m <- max(mean(g), 1e-8)                       # intercept-only fallback
  list(mu0 = rep(m, length(g)),
       coef = stats::setNames(c(log(m), rep(NA_real_, ncol(X) - 1L)), colnames(X)),
       deviance = NA_real_, converged = FALSE)
}

# EM for pi only (mu0 and K fixed -> component densities are fixed, computed
# once). Returns the MLE pi, the per-cell responsibilities, and the log-lik.
poisthresh_em <- function(g, mu0, K, pi_init = 0.05, max_iter = 200L, tol = 1e-8) {
  lf0 <- stats::dpois(g, lambda = mu0,     log = TRUE)
  lf1 <- stats::dpois(g, lambda = mu0 + K, log = TRUE)
  pi_hat <- min(max(pi_init, 1e-10), 1 - 1e-10)
  ll_old <- -Inf; converged <- FALSE; iter <- 0L; ll <- -Inf
  for (iter in seq_len(max_iter)) {
    lp0 <- log1p(-pi_hat) + lf0
    lp1 <- log(pi_hat)    + lf1
    ll  <- sum(poisthresh_logaddexp(lp0, lp1))
    if (is.finite(ll) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) {
      converged <- TRUE; ll_old <- ll; break
    }
    ll_old <- ll
    r <- stats::plogis(lp1 - lp0)               # P(pert | y_i)
    pi_hat <- min(max(mean(r), 1e-10), 1 - 1e-10)
  }
  # final E-step at the converged pi so r matches the reported pi
  r <- stats::plogis((log(pi_hat) + lf1) - (log1p(-pi_hat) + lf0))
  r[!is.finite(r)] <- 0
  list(pi = pi_hat, r = r, loglik = ll_old, converged = converged, iter = iter)
}

assign_one_grna_poisthresh <- function(g, other_umis, other_nnz, K,
                                       n_fit_cutoff          = 10L,
                                       backup_threshold      = 5L,
                                       probability_threshold = 0.8,
                                       pi_init               = 0.05,
                                       max_iter              = 200L,
                                       tol                   = 1e-8,
                                       keep_fits             = FALSE) {
  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)
  n_fit     <- length(g)

  em_fit       <- NULL
  em_converged <- NA
  em_log_lik   <- NA_real_
  poisthresh   <- list(pi = NA_real_, K = K, b0 = NA_real_, b_umis = NA_real_,
                       b_nnz = NA_real_, mean_mu0 = NA_real_, mean_r = NA_real_,
                       offset_converged = NA, n_fit = n_fit)
  offset_model_summary <- NULL
  prob_quantile_probs  <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles <- stats::setNames(rep(NA_real_, length(prob_quantile_probs)),
                                    paste0(prob_quantile_probs * 100, "%"))
  n_above_prob_thresh <- NA_integer_

  if (n_nonzero >= n_fit_cutoff) {
    X <- cbind(1, log1p(other_umis), log1p(other_nnz))
    colnames(X) <- c("(Intercept)", "log1p(other_grna_n_umis)", "log1p(other_grna_n_nonzero)")
    off <- poisthresh_fit_offset(g, X)
    em_fit <- tryCatch(poisthresh_em(g, off$mu0, K, pi_init, max_iter, tol),
                       error = function(e) NULL)
    if (!is.null(em_fit) && is.finite(em_fit$loglik)) {
      r_post      <- em_fit$r
      assignments <- which(r_post >= probability_threshold)
      prob_quantiles      <- stats::quantile(r_post, probs = prob_quantile_probs)
      n_above_prob_thresh <- length(assignments)
      poisthresh <- list(
        pi = em_fit$pi, K = K,
        b0 = unname(off$coef[1L]), b_umis = unname(off$coef[2L]), b_nnz = unname(off$coef[3L]),
        mean_mu0 = mean(off$mu0), mean_r = mean(r_post),
        offset_converged = off$converged, n_fit = n_fit
      )
      offset_model_summary <- list(coef = off$coef, deviance = off$deviance,
                                   converged = off$converged)
      em_converged <- em_fit$converged
      em_log_lik   <- em_fit$loglik
    } else {
      assignments  <- which(g >= backup_threshold)
      em_converged <- FALSE
    }
  } else {
    assignments <- which(g >= backup_threshold)
  }

  elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  out <- list(
    assignments          = assignments,
    n_nonzero            = n_nonzero,
    n_fit                = n_fit,
    n_assigned           = length(assignments),
    em_converged         = em_converged,
    em_log_lik           = em_log_lik,
    em_init_i            = NA_integer_,        # single EM run (pi-likelihood is concave)
    poisthresh           = poisthresh,
    prob_quantiles       = prob_quantiles,
    n_above_prob_thresh  = n_above_prob_thresh,
    offset_model_summary = offset_model_summary,
    elapsed_sec          = elapsed_sec
  )
  if (keep_fits) out$em_fit <- em_fit
  out
}

# Driver. `grna_n_nonzero` / `grna_n_umis` are the per-cell TOTALS over all
# guides (this guide included); each worker removes the guide's own contribution
# before fitting the offset (the circularity correction).
poisthresh_assign <- function(grna_matrix,
                              grna_n_nonzero,
                              grna_n_umis,
                              K,
                              grna_ids              = rownames(grna_matrix),
                              n_fit_cutoff          = 10L,
                              backup_threshold      = 5L,
                              probability_threshold = 0.8,
                              pi_init               = 0.05,
                              max_iter              = 200L,
                              tol                   = 1e-8,
                              cl                    = NULL,
                              keep_fits             = FALSE) {
  force(grna_matrix); force(grna_n_nonzero); force(grna_n_umis); force(K)
  force(n_fit_cutoff); force(backup_threshold); force(probability_threshold)
  force(pi_init); force(max_iter); force(tol); force(keep_fits)

  worker <- function(grna_id) {
    g          <- as.numeric(grna_matrix[grna_id, ])
    other_umis <- pmax(grna_n_umis    - g,                  0)   # - y_i
    other_nnz  <- pmax(grna_n_nonzero - as.numeric(g > 0),  0)   # - 1{y_i > 0}
    assign_one_grna_poisthresh(
      g = g, other_umis = other_umis, other_nnz = other_nnz, K = K,
      n_fit_cutoff = n_fit_cutoff, backup_threshold = backup_threshold,
      probability_threshold = probability_threshold, pi_init = pi_init,
      max_iter = max_iter, tol = tol, keep_fits = keep_fits
    )
  }

  per_guide <- if (is.null(cl)) lapply(grna_ids, worker) else parallel::parLapplyLB(cl, grna_ids, worker)
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  run_meta <- list(
    model                 = "pois-thresh-glmpois-offset",
    offset_covariates     = "log1p(grna_n_umis - y), log1p(grna_n_nonzero - 1{y>0})",
    mixture               = "(1-pi) Pois(mu0) + pi Pois(mu0 + K)",
    K                     = K,
    n_free_params         = 1L,                # pi only (K fixed)
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    pi_init               = pi_init,
    max_iter              = max_iter,
    tol                   = tol
  )

  list(per_guide = per_guide, run_meta = run_meta)
}


# ---- Standard glmpois offset + fixed-shift mixture (glmpois_poisthresh<K>) ----
# Same fixed-K Poisson mixture as above, but the offset is the VANILLA glmpois
# baseline: a Poisson GLM on the shared covariate matrix built once from the
# run's formula (response + grna covariates, log1p), with NO circularity
# correction. This reuses `poisthresh_em`; only the offset model differs.
# `fit_baseline_glm_pure_R` is sourced from sceptre_assign_pure_R.R.

assign_one_grna_glmpois_poisthresh <- function(g, covariate_matrix, K,
                                               n_fit_cutoff          = 10L,
                                               backup_threshold      = 5L,
                                               probability_threshold = 0.8,
                                               pi_init               = 0.05,
                                               max_iter              = 200L,
                                               tol                   = 1e-8,
                                               keep_fits             = FALSE) {
  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)
  n_fit     <- length(g)

  em_fit       <- NULL
  em_converged <- NA
  em_log_lik   <- NA_real_
  poisthresh   <- list(pi = NA_real_, K = K, coef = NA_real_, mean_mu0 = NA_real_,
                       mean_r = NA_real_, offset_converged = NA, n_fit = n_fit)
  offset_model_summary <- NULL
  prob_quantile_probs  <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles <- stats::setNames(rep(NA_real_, length(prob_quantile_probs)),
                                    paste0(prob_quantile_probs * 100, "%"))
  n_above_prob_thresh <- NA_integer_

  if (n_nonzero >= n_fit_cutoff) {
    off          <- fit_baseline_glm_pure_R(g, covariate_matrix)
    mu0          <- pmax(off$fitted.values, 1e-300)
    off_conv     <- isTRUE(off$converged)
    em_fit <- tryCatch(poisthresh_em(g, mu0, K, pi_init, max_iter, tol),
                       error = function(e) NULL)
    if (!is.null(em_fit) && is.finite(em_fit$loglik)) {
      r_post      <- em_fit$r
      assignments <- which(r_post >= probability_threshold)
      prob_quantiles      <- stats::quantile(r_post, probs = prob_quantile_probs)
      n_above_prob_thresh <- length(assignments)
      poisthresh <- list(
        pi = em_fit$pi, K = K, coef = off$coefficients,
        mean_mu0 = mean(mu0), mean_r = mean(r_post),
        offset_converged = off_conv, n_fit = n_fit
      )
      offset_model_summary <- off$offset_model_summary
      em_converged <- em_fit$converged
      em_log_lik   <- em_fit$loglik
    } else {
      assignments  <- which(g >= backup_threshold)
      em_converged <- FALSE
    }
  } else {
    assignments <- which(g >= backup_threshold)
  }

  elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  out <- list(
    assignments          = assignments,
    n_nonzero            = n_nonzero,
    n_fit                = n_fit,
    n_assigned           = length(assignments),
    em_converged         = em_converged,
    em_log_lik           = em_log_lik,
    em_init_i            = NA_integer_,        # single EM run (pi-likelihood is concave)
    poisthresh           = poisthresh,
    prob_quantiles       = prob_quantiles,
    n_above_prob_thresh  = n_above_prob_thresh,
    offset_model_summary = offset_model_summary,
    elapsed_sec          = elapsed_sec
  )
  if (keep_fits) out$em_fit <- em_fit
  out
}

# Driver. `covariate_matrix` is the shared design built once from the formula
# (vanilla glmpois covariates); the SAME matrix is used for every guide, so
# there is no circularity correction.
glmpois_poisthresh_assign <- function(grna_matrix,
                                      covariate_matrix,
                                      K,
                                      grna_ids              = rownames(grna_matrix),
                                      n_fit_cutoff          = 10L,
                                      backup_threshold      = 5L,
                                      probability_threshold = 0.8,
                                      pi_init               = 0.05,
                                      max_iter              = 200L,
                                      tol                   = 1e-8,
                                      cl                    = NULL,
                                      keep_fits             = FALSE) {
  force(grna_matrix); force(covariate_matrix); force(K)
  force(n_fit_cutoff); force(backup_threshold); force(probability_threshold)
  force(pi_init); force(max_iter); force(tol); force(keep_fits)

  worker <- function(grna_id) {
    g <- as.numeric(grna_matrix[grna_id, ])
    assign_one_grna_glmpois_poisthresh(
      g = g, covariate_matrix = covariate_matrix, K = K,
      n_fit_cutoff = n_fit_cutoff, backup_threshold = backup_threshold,
      probability_threshold = probability_threshold, pi_init = pi_init,
      max_iter = max_iter, tol = tol, keep_fits = keep_fits
    )
  }

  per_guide <- if (is.null(cl)) lapply(grna_ids, worker) else parallel::parLapplyLB(cl, grna_ids, worker)
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  run_meta <- list(
    model                 = "pois-thresh-glmpois-standard-offset",
    offset_covariates     = "vanilla glmpois covariate_matrix (formula, log1p), no circularity correction",
    mixture               = "(1-pi) Pois(mu0) + pi Pois(mu0 + K)",
    K                     = K,
    n_free_params         = 1L,                # pi only (K fixed)
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    pi_init               = pi_init,
    max_iter              = max_iter,
    tol                   = tol
  )

  list(per_guide = per_guide, run_meta = run_meta)
}
