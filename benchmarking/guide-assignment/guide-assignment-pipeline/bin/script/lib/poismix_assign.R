# Joint Poisson-Poisson mixture guide-assignment model (no offset model).
#
# Like the gated variant, there is NO separate offset model: the covariate
# coefficient is fit *jointly* with the mixture inside the EM (rather than
# stagewise, where one first fits an offset model and then runs the mixture
# with the offsets fixed).
#
# Per guide, over ALL cells:
#   Z_i  ~ Bernoulli(pi)                                   # perturbation indicator
#   Y_i | Z_i ~ Poisson(mu_i)                              # UMI count
#   log(mu_i) = b0 + b1 Z_i + b2 x_i
#   x_i = log1p(grna_n_nonzero_i - 1{y_i > 0})             # decircularized covariate
#
# Four parameters: pi, b0 (background log-mean), b1 (perturbation effect),
# b2 (covariate coefficient, shared across components). The label convention
# b1 >= 0 forces Z = 1 to be the higher-mean (perturbed) component.
#
# A cell is assigned perturbed iff its posterior responsibility
#   r_i = P(Z_i = 1 | y_i, x_i) = pi f1(y_i) / [(1 - pi) f0(y_i) + pi f1(y_i)]
# is >= probability_threshold (default 0.8).
#
# Initialization mirrors sceptre's multi-start scheme (random (pi, g_pert)
# seed grid generated once by the driver under an isolated seed, shared across
# guides; pick the best converged start).

# ---- helpers -----------------------------------------------------------------

ETA_CAP <- 30  # clamp linear predictors to avoid exp() overflow

# Weighted Poisson Fisher scoring for the M-step regression. The complete-data
# weighted log-likelihood is
#   sum_i (1 - r_i)[y_i eta0_i - exp(eta0_i)] + r_i[y_i eta1_i - exp(eta1_i)]
# with eta0 = b0 + b2 x, eta1 = b0 + b1 + b2 x. This is the row-doubled weighted
# Poisson GLM the prototype fits with glm(), solved here directly on the
# 3-parameter design [1, z, x] without materializing the 2n expanded rows.
# A tiny ridge keeps the 3x3 information matrix invertible when a component
# carries ~no weight (b1 unidentified) or x is near-constant.
# Weighted Poisson expected complete-data log-lik (drops the y-only constant),
# the objective the M-step maximizes. Used for backtracking line search.
poismix_q <- function(beta, y, x, w0, w1) {
  eta0 <- pmin(pmax(beta[1L] + beta[3L] * x, -ETA_CAP), ETA_CAP)
  eta1 <- pmin(pmax(beta[1L] + beta[2L] + beta[3L] * x, -ETA_CAP), ETA_CAP)
  sum(w0 * (y * eta0 - exp(eta0))) + sum(w1 * (y * eta1 - exp(eta1)))
}

poismix_irls_poisson <- function(y, x, r, beta_start = c(0, 0, 0),
                                 n_iter = 50L, tol = 1e-8, ridge = 1e-6) {
  beta <- beta_start
  w0 <- 1 - r; w1 <- r
  q_cur <- poismix_q(beta, y, x, w0, w1)
  for (i in seq_len(n_iter)) {
    b0 <- beta[1L]; b1 <- beta[2L]; b2 <- beta[3L]
    eta0 <- pmin(pmax(b0 + b2 * x, -ETA_CAP), ETA_CAP)
    eta1 <- pmin(pmax(b0 + b1 + b2 * x, -ETA_CAP), ETA_CAP)
    mu0 <- exp(eta0); mu1 <- exp(eta1)
    Wm0 <- w0 * mu0; Wm1 <- w1 * mu1
    res0 <- w0 * (y - mu0); res1 <- w1 * (y - mu1)

    # gradient (score) w.r.t. (b0, b1, b2)
    g <- c(sum(res0) + sum(res1),
           sum(res1),
           sum(x * res0) + sum(x * res1))
    # Fisher information (X' W X), 3x3 symmetric
    H11 <- sum(Wm0) + sum(Wm1)
    H12 <- sum(Wm1)
    H13 <- sum(x * Wm0) + sum(x * Wm1)
    H22 <- sum(Wm1)
    H23 <- sum(x * Wm1)
    H33 <- sum(x * x * Wm0) + sum(x * x * Wm1)
    H <- matrix(c(H11, H12, H13,
                  H12, H22, H23,
                  H13, H23, H33), 3L, 3L)
    diag(H) <- diag(H) + ridge
    delta <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(delta) || any(!is.finite(delta))) break
    # backtracking line search: guarantee a monotone increase in Q
    step <- 1; accepted <- FALSE
    repeat {
      cand <- beta + step * delta
      q_new <- poismix_q(cand, y, x, w0, w1)
      if (is.finite(q_new) && q_new >= q_cur - 1e-9) { accepted <- TRUE; break }
      step <- step / 2
      if (step < 1e-10) break
    }
    if (!accepted) break
    beta <- beta + step * delta; q_cur <- q_new
    if (max(abs(step * delta)) < tol) break
  }
  beta
}

# One M-step. `r` = current responsibilities (P(pert)). Enforces b1 >= 0 by
# relabeling so Z = 1 is the high-mean component.
poismix_m_step <- function(y, x, r, beta_start = c(0, 0, 0),
                          eps = 1e-10, ridge = 1e-6) {
  pi_hat <- min(max(mean(r), eps), 1 - eps)
  beta <- poismix_irls_poisson(y, x, r, beta_start = beta_start, ridge = ridge)
  b0 <- beta[1L]; b1 <- beta[2L]; b2 <- beta[3L]
  if (is.finite(b1) && b1 < 0) {            # label flip: Z = 1 must be high mean
    b0 <- b0 + b1; b1 <- -b1; pi_hat <- 1 - pi_hat
  }
  list(pi = pi_hat, b0 = b0, b1 = b1, b2 = b2)
}

# One E-step. Returns responsibilities + total log-likelihood (log-sum-exp).
poismix_e_step <- function(y, x, par, eps = 1e-10) {
  pi_hat <- min(max(par$pi, eps), 1 - eps)
  eta0 <- pmin(pmax(par$b0 + par$b2 * x, -ETA_CAP), ETA_CAP)
  eta1 <- pmin(pmax(par$b0 + par$b1 + par$b2 * x, -ETA_CAP), ETA_CAP)
  lf0 <- stats::dpois(y, lambda = exp(eta0), log = TRUE)
  lf1 <- stats::dpois(y, lambda = exp(eta1), log = TRUE)
  ln0 <- log1p(-pi_hat) + lf0
  ln1 <- log(pi_hat) + lf1
  m <- pmax(ln0, ln1)
  ll_cell <- m + log(exp(ln0 - m) + exp(ln1 - m))
  r <- exp(ln1 - ll_cell)
  r[!is.finite(r)] <- 0
  list(r = r, ll = sum(ll_cell))
}

# Run EM from a single initialization (responsibilities `r_init`). Mirrors
# gated_em_one_start: M-step then E-step, report the matched (par, r, ll).
poismix_em_one_start <- function(y, x, r_init, max_iter = 200L, tol = 1e-8,
                                eps = 1e-10, ridge = 1e-6) {
  r <- r_init
  ll_old <- -Inf
  converged <- FALSE
  par <- NULL
  beta_warm <- c(0, 0, 0)
  iter <- 0L
  ll <- -Inf
  for (iter in seq_len(max_iter)) {
    par <- poismix_m_step(y, x, r, beta_start = beta_warm, eps = eps, ridge = ridge)
    beta_warm <- c(par$b0, par$b1, par$b2)
    es <- poismix_e_step(y, x, par, eps = eps)
    r <- es$r
    ll <- es$ll
    if (is.finite(ll) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) {
      converged <- TRUE; ll_old <- ll; break
    }
    ll_old <- ll
  }
  list(par = par, r = r, loglik = ll_old, converged = converged, iter = iter)
}

# Build initial params from a sceptre-style random seed:
#   pi0      -- initial (constant) prior perturbation probability
#   g_pert0  -- initial perturbed-component mean on the log(count) scale
# The background intercept b0 is method-of-moments on cells below exp(g_pert0)
# (the presumed non-pert background); b2 starts at 0; b1 = g_pert0 - b0 (>= 0).
poismix_init_par_from_seed <- function(y, x, pi0, g_pert0, eps = 1e-10) {
  pi0 <- min(max(pi0, eps), 1 - eps)
  # background from cells below the seeded perturbed mean (exclude pert outliers)
  lo <- y[y < exp(g_pert0)]
  if (length(lo) < 5L) lo <- y
  b0 <- log(max(mean(lo), 0.01))
  b1 <- max(g_pert0 - b0, 0.1)
  list(pi = pi0, b0 = b0, b1 = b1, b2 = 0)
}

# Multi-start EM, sceptre-style. `pi_guesses` / `g_pert_guesses` are generated
# once by the driver under an isolated seed and shared across guides. Selection
# mirrors run_em_pois_pois_pure_R's `converged && loglik > best`.
run_em_poismix <- function(y, x, pi_guesses, g_pert_guesses,
                          max_iter = 200L, tol = 1e-8, eps = 1e-10, ridge = 1e-6) {
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B)
  best <- NULL; best_i <- NA_integer_
  better <- function(cand, cur) {
    if (is.null(cur)) return(TRUE)
    if (cand$converged != cur$converged) return(cand$converged)
    is.finite(cand$loglik) &&
      (!is.finite(cur$loglik) || cand$loglik > cur$loglik)
  }
  for (i in seq_len(B)) {
    par0 <- poismix_init_par_from_seed(y, x, pi_guesses[i], g_pert_guesses[i], eps)
    r0 <- poismix_e_step(y, x, par0, eps = eps)$r
    fit <- poismix_em_one_start(y, x, r0, max_iter, tol, eps, ridge)
    if (better(fit, best)) { best <- fit; best_i <- i }
  }
  list(
    pi = best$par$pi, b0 = best$par$b0, b1 = best$par$b1, b2 = best$par$b2,
    r = best$r, loglik = best$loglik,
    converged = best$converged, iter = best$iter, start_i = best_i
  )
}

# ---- per-guide -----------------------------------------------------------------

assign_one_grna_poismix <- function(g, grna_n_nonzero,
                                    pi_guesses, g_pert_guesses,
                                    n_fit_cutoff          = 10L,
                                    backup_threshold      = 5L,
                                    probability_threshold = 0.8,
                                    max_iter              = 200L,
                                    tol                   = 1e-8,
                                    eps                   = 1e-10,
                                    ridge                 = 1e-6,
                                    keep_fits             = FALSE) {
  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)

  # Guide-specific decircularized covariate (subtract this guide's own presence).
  other_nnz <- pmax(grna_n_nonzero - as.numeric(g > 0), 0)
  x <- log1p(other_nnz)
  y <- as.numeric(g)
  ok <- is.finite(y) & is.finite(x)
  n_fit <- sum(ok)

  # defaults for the below-cutoff / failure paths
  em_fit       <- NULL
  em_converged <- NA
  em_log_lik   <- NA_real_
  em_init_i    <- NA_integer_
  poismix      <- list(pi = NA_real_, b0 = NA_real_, b1 = NA_real_, b2 = NA_real_,
                       mean_mu_np = NA_real_, mean_mu_pert = NA_real_,
                       mean_r = NA_real_, n_fit = n_fit)
  prob_quantile_probs <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles <- stats::setNames(rep(NA_real_, length(prob_quantile_probs)),
                                    paste0(prob_quantile_probs * 100, "%"))
  n_above_prob_thresh <- NA_integer_

  # Need enough informative (nonzero) cells to identify the perturbed component.
  if (n_nonzero >= n_fit_cutoff && n_fit >= n_fit_cutoff) {
    yf <- y[ok]; xf <- x[ok]
    fit_idx <- which(ok)
    em_fit <- tryCatch(
      run_em_poismix(yf, xf, pi_guesses, g_pert_guesses, max_iter, tol, eps, ridge),
      error = function(e) NULL
    )
    if (!is.null(em_fit) && is.finite(em_fit$loglik)) {
      r_post <- em_fit$r
      assignments <- fit_idx[r_post >= probability_threshold]
      prob_quantiles <- stats::quantile(r_post, probs = prob_quantile_probs)
      n_above_prob_thresh <- length(assignments)
      eta0 <- pmin(pmax(em_fit$b0 + em_fit$b2 * xf, -ETA_CAP), ETA_CAP)
      eta1 <- pmin(pmax(em_fit$b0 + em_fit$b1 + em_fit$b2 * xf, -ETA_CAP), ETA_CAP)
      poismix <- list(
        pi = em_fit$pi, b0 = em_fit$b0, b1 = em_fit$b1, b2 = em_fit$b2,
        mean_mu_np = mean(exp(eta0)), mean_mu_pert = mean(exp(eta1)),
        mean_r = mean(r_post), n_fit = n_fit
      )
      em_converged <- em_fit$converged
      em_log_lik   <- em_fit$loglik
      em_init_i    <- em_fit$start_i
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
    em_init_i            = em_init_i,
    poismix              = poismix,
    prob_quantiles       = prob_quantiles,
    n_above_prob_thresh  = n_above_prob_thresh,
    offset_model_summary = NULL,            # no offset model in this variant
    elapsed_sec          = elapsed_sec
  )
  if (keep_fits) out$em_fit <- em_fit
  out
}

# ---- driver --------------------------------------------------------------------
# `grna_n_nonzero` is the shared per-cell covariate (number of nonzero guides in
# each cell); each guide subtracts its own presence (decircularization) before
# taking log1p.
poismix_assign <- function(grna_matrix,
                          grna_n_nonzero,
                          grna_ids              = rownames(grna_matrix),
                          n_em_rep              = 5L,
                          pi_guess_range        = c(1e-5, 0.1),
                          g_pert_guess_range    = log(c(10, 5000)),
                          n_fit_cutoff          = 10L,
                          backup_threshold      = 5L,
                          probability_threshold = 0.8,
                          max_iter              = 200L,
                          tol                   = 1e-8,
                          eps                   = 1e-10,
                          ridge                 = 1e-6,
                          seed                  = 4L,
                          cl                    = NULL,
                          keep_fits             = FALSE) {
  force(grna_matrix); force(grna_n_nonzero)
  force(n_fit_cutoff); force(backup_threshold); force(probability_threshold)
  force(max_iter); force(tol); force(eps); force(ridge); force(keep_fits)

  # Starting guesses generated ONCE under an isolated seed and shared across all
  # guides, exactly like sceptre_assign_pure_R's driver. The global RNG state is
  # saved/restored so callers' streams aren't perturbed.
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
  force(pi_guesses); force(g_pert_guesses)

  worker <- function(grna_id) {
    g <- as.numeric(grna_matrix[grna_id, ])
    assign_one_grna_poismix(
      g                     = g,
      grna_n_nonzero        = grna_n_nonzero,
      pi_guesses            = pi_guesses,
      g_pert_guesses        = g_pert_guesses,
      n_fit_cutoff          = n_fit_cutoff,
      backup_threshold      = backup_threshold,
      probability_threshold = probability_threshold,
      max_iter              = max_iter,
      tol                   = tol,
      eps                   = eps,
      ridge                 = ridge,
      keep_fits             = keep_fits
    )
  }

  per_guide <- if (is.null(cl)) {
    lapply(grna_ids, worker)
  } else {
    parallel::parLapplyLB(cl, grna_ids, worker)
  }
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  run_meta <- list(
    model                 = "pois-pois-joint-nnz",
    covariate             = "log1p(grna_n_nonzero - 1{y>0})",
    response              = "y (all cells)",
    n_params              = 4L,
    n_em_rep              = n_em_rep,
    pi_guess_range        = pi_guess_range,
    g_pert_guess_range    = g_pert_guess_range,
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    max_iter              = max_iter,
    tol                   = tol,
    ridge                 = ridge,
    seed                  = seed
  )

  list(per_guide = per_guide, run_meta = run_meta)
}

# ==============================================================================
# pois0 speedup variant: cells with y_i = 0 are forced to the non-pert
# component (responsibility 0), exactly as run_em_pois0_nb_pure_R does. This is
# an excellent approximation because the perturbed component has a high mean, so
# P(pert | y = 0) = pi dpois(0; mu1) / [...] is essentially 0.
#
# The constraint lets the E-step and the *positive* part of the M-step touch
# only the positive cells. Unlike pois0nb (fixed offset mu0), our non-pert mean
# mu0_i = exp(b0 + b2 x_i) is estimated jointly, so the zeros still contribute
# to the (b0, b2) score. But every zero has x = log1p(grna_n_nonzero) with
# integer grna_n_nonzero, so the zeros collapse to a few distinct x values; we
# precompute (value, count) bins once and add their closed-form contribution to
# the 3x3 score/information each Fisher step. Cost per EM iteration drops from
# O(n) to O(n_pos + n_distinct_nnz).
#
# Zeros' observed-data contribution (definitely non-pert):
#   sum_zero log[(1 - pi) dpois(0; mu0_i)] = n_zero log(1 - pi) - sum_zero mu0_i
# with sum_zero mu0_i = sum_k cnt_k exp(b0 + b2 v_k).
# ==============================================================================

# Build the per-guide fit object: positive cells (full mixture) + zero bins.
poismix0_make_fit <- function(g, grna_n_nonzero) {
  is_pos <- g > 0
  y_pos  <- as.numeric(g[is_pos])
  # decircularized covariate; positives subtract their own presence (-1).
  x_pos  <- log1p(pmax(grna_n_nonzero[is_pos] - 1, 0))
  # zeros: y = 0 so other_nnz = grna_n_nonzero; bin by distinct nnz value.
  nnz_zero <- grna_n_nonzero[!is_pos]
  if (length(nnz_zero)) {
    tab <- table(nnz_zero)
    zv  <- log1p(as.numeric(names(tab)))
    zc  <- as.numeric(tab)
  } else {
    zv <- numeric(0); zc <- numeric(0)
  }
  list(y_pos = y_pos, x_pos = x_pos, zv = zv, zc = zc,
       n_total = length(g), n_pos = length(y_pos), n_zero = sum(!is_pos),
       pos_idx = which(is_pos))
}

# Objective for the pois0 M-step: positives' weighted Poisson log-lik plus the
# zeros' (non-pert) contribution -sum_zero mu0 (binned). Drops the y-only const.
poismix0_q <- function(beta, fit, w0, w1) {
  y <- fit$y_pos; x <- fit$x_pos
  eta0 <- pmin(pmax(beta[1L] + beta[3L] * x, -ETA_CAP), ETA_CAP)
  eta1 <- pmin(pmax(beta[1L] + beta[2L] + beta[3L] * x, -ETA_CAP), ETA_CAP)
  q <- sum(w0 * (y * eta0 - exp(eta0))) + sum(w1 * (y * eta1 - exp(eta1)))
  if (length(fit$zv)) {
    q <- q - sum(fit$zc * exp(pmin(pmax(beta[1L] + beta[3L] * fit$zv, -ETA_CAP), ETA_CAP)))
  }
  q
}

# Weighted Poisson Fisher scoring over positives (both components) plus the
# binned zeros (non-pert component only). b1 kept >= 0 by projection, with
# backtracking line search for monotone increase in the objective.
poismix0_irls <- function(fit, r, beta_start = c(0, 0, 0),
                          n_iter = 50L, tol = 1e-8, ridge = 1e-6) {
  beta <- beta_start
  y <- fit$y_pos; x <- fit$x_pos; zv <- fit$zv; zc <- fit$zc
  w0 <- 1 - r; w1 <- r
  q_cur <- poismix0_q(beta, fit, w0, w1)
  for (i in seq_len(n_iter)) {
    b0 <- beta[1L]; b1 <- beta[2L]; b2 <- beta[3L]
    # positives, both components
    eta0 <- pmin(pmax(b0 + b2 * x, -ETA_CAP), ETA_CAP)
    eta1 <- pmin(pmax(b0 + b1 + b2 * x, -ETA_CAP), ETA_CAP)
    mu0 <- exp(eta0); mu1 <- exp(eta1)
    Wm0 <- w0 * mu0; Wm1 <- w1 * mu1
    res0 <- w0 * (y - mu0); res1 <- w1 * (y - mu1)
    g1 <- sum(res0) + sum(res1); g2 <- sum(res1); g3 <- sum(x * res0) + sum(x * res1)
    H11 <- sum(Wm0) + sum(Wm1); H12 <- sum(Wm1); H13 <- sum(x * Wm0) + sum(x * Wm1)
    H22 <- sum(Wm1); H23 <- sum(x * Wm1); H33 <- sum(x * x * Wm0) + sum(x * x * Wm1)
    # binned zeros, non-pert component only (y = 0, z = 0, weight = count)
    if (length(zv)) {
      mu0z <- exp(pmin(pmax(b0 + b2 * zv, -ETA_CAP), ETA_CAP))
      Zmu <- sum(zc * mu0z); Zvmu <- sum(zc * zv * mu0z); Zv2mu <- sum(zc * zv * zv * mu0z)
      g1 <- g1 - Zmu; g3 <- g3 - Zvmu
      H11 <- H11 + Zmu; H13 <- H13 + Zvmu; H33 <- H33 + Zv2mu
    }
    g <- c(g1, g2, g3)
    H <- matrix(c(H11, H12, H13, H12, H22, H23, H13, H23, H33), 3L, 3L)
    diag(H) <- diag(H) + ridge
    delta <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(delta) || any(!is.finite(delta))) break
    step <- 1; accepted <- FALSE
    repeat {
      cand <- beta + step * delta
      if (cand[2L] < 0) cand[2L] <- 0        # project b1 >= 0 (zeros anchor comp0 low)
      q_new <- poismix0_q(cand, fit, w0, w1)
      if (is.finite(q_new) && q_new >= q_cur - 1e-9) { accepted <- TRUE; break }
      step <- step / 2
      if (step < 1e-10) break
    }
    if (!accepted) break
    beta <- cand; q_cur <- q_new
    if (max(abs(step * delta)) < tol) break
  }
  beta
}

poismix0_m_step <- function(fit, r, beta_start = c(0, 0, 0), eps = 1e-10, ridge = 1e-6) {
  pi_hat <- min(max(sum(r) / fit$n_total, eps), 1 - eps)   # zeros weight 0
  beta <- poismix0_irls(fit, r, beta_start = beta_start, ridge = ridge)
  list(pi = pi_hat, b0 = beta[1L], b1 = beta[2L], b2 = beta[3L])
}

# E-step over positives; zeros forced non-pert. loglik adds the zeros' constant.
poismix0_e_step <- function(fit, par, eps = 1e-10) {
  pi_hat <- min(max(par$pi, eps), 1 - eps)
  y <- fit$y_pos; x <- fit$x_pos
  eta0 <- pmin(pmax(par$b0 + par$b2 * x, -ETA_CAP), ETA_CAP)
  eta1 <- pmin(pmax(par$b0 + par$b1 + par$b2 * x, -ETA_CAP), ETA_CAP)
  lf0 <- stats::dpois(y, lambda = exp(eta0), log = TRUE)
  lf1 <- stats::dpois(y, lambda = exp(eta1), log = TRUE)
  ln0 <- log1p(-pi_hat) + lf0
  ln1 <- log(pi_hat) + lf1
  m <- pmax(ln0, ln1)
  ll_pos <- m + log(exp(ln0 - m) + exp(ln1 - m))
  r <- exp(ln1 - ll_pos)
  r[!is.finite(r)] <- 0
  # zeros: n_zero log(1 - pi) - sum_zero mu0
  zsum <- if (length(fit$zv)) {
    sum(fit$zc * exp(pmin(pmax(par$b0 + par$b2 * fit$zv, -ETA_CAP), ETA_CAP)))
  } else 0
  ll <- sum(ll_pos) + fit$n_zero * log1p(-pi_hat) - zsum
  list(r = r, ll = ll)
}

poismix0_em_one_start <- function(fit, r_init, max_iter = 200L, tol = 1e-8,
                                  eps = 1e-10, ridge = 1e-6) {
  r <- r_init; ll_old <- -Inf; converged <- FALSE; par <- NULL
  beta_warm <- c(0, 0, 0); iter <- 0L; ll <- -Inf
  for (iter in seq_len(max_iter)) {
    par <- poismix0_m_step(fit, r, beta_start = beta_warm, eps = eps, ridge = ridge)
    beta_warm <- c(par$b0, par$b1, par$b2)
    es <- poismix0_e_step(fit, par, eps = eps)
    r <- es$r; ll <- es$ll
    if (is.finite(ll) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) {
      converged <- TRUE; ll_old <- ll; break
    }
    ll_old <- ll
  }
  list(par = par, r = r, loglik = ll_old, converged = converged, iter = iter)
}

poismix0_init_par_from_seed <- function(fit, pi0, g_pert0, eps = 1e-10) {
  pi0 <- min(max(pi0, eps), 1 - eps)
  # background = mean over presumed-non-pert cells: positives below the seeded
  # perturbed mean plus ALL zeros (which pull the mean toward the true low level
  # and prevent the pert outliers from inflating b0).
  lo <- fit$y_pos[fit$y_pos < exp(g_pert0)]
  bg_mean <- sum(lo) / (length(lo) + fit$n_zero)
  b0 <- log(max(bg_mean, 1e-3))
  b1 <- max(g_pert0 - b0, 0.1)
  list(pi = pi0, b0 = b0, b1 = b1, b2 = 0)
}

run_em_poismix0 <- function(fit, pi_guesses, g_pert_guesses,
                            max_iter = 200L, tol = 1e-8, eps = 1e-10, ridge = 1e-6) {
  B <- length(pi_guesses)
  stopifnot(length(g_pert_guesses) == B)
  best <- NULL; best_i <- NA_integer_
  better <- function(cand, cur) {
    if (is.null(cur)) return(TRUE)
    if (cand$converged != cur$converged) return(cand$converged)
    is.finite(cand$loglik) && (!is.finite(cur$loglik) || cand$loglik > cur$loglik)
  }
  for (i in seq_len(B)) {
    par0 <- poismix0_init_par_from_seed(fit, pi_guesses[i], g_pert_guesses[i], eps)
    r0 <- poismix0_e_step(fit, par0, eps = eps)$r
    fit_i <- poismix0_em_one_start(fit, r0, max_iter, tol, eps, ridge)
    if (better(fit_i, best)) { best <- fit_i; best_i <- i }
  }
  list(pi = best$par$pi, b0 = best$par$b0, b1 = best$par$b1, b2 = best$par$b2,
       r = best$r, loglik = best$loglik,
       converged = best$converged, iter = best$iter, start_i = best_i)
}

assign_one_grna_poismix0 <- function(g, grna_n_nonzero,
                                     pi_guesses, g_pert_guesses,
                                     n_fit_cutoff          = 10L,
                                     backup_threshold      = 5L,
                                     probability_threshold = 0.8,
                                     max_iter              = 200L,
                                     tol                   = 1e-8,
                                     eps                   = 1e-10,
                                     ridge                 = 1e-6,
                                     keep_fits             = FALSE) {
  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)

  em_fit       <- NULL
  em_converged <- NA
  em_log_lik   <- NA_real_
  em_init_i    <- NA_integer_
  poismix      <- list(pi = NA_real_, b0 = NA_real_, b1 = NA_real_, b2 = NA_real_,
                       mean_mu_np = NA_real_, mean_mu_pert = NA_real_,
                       mean_r = NA_real_, n_fit = length(g))
  prob_quantile_probs <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles <- stats::setNames(rep(NA_real_, length(prob_quantile_probs)),
                                    paste0(prob_quantile_probs * 100, "%"))
  n_above_prob_thresh <- NA_integer_

  if (n_nonzero >= n_fit_cutoff) {
    fit <- poismix0_make_fit(g, grna_n_nonzero)
    em_fit <- tryCatch(
      run_em_poismix0(fit, pi_guesses, g_pert_guesses, max_iter, tol, eps, ridge),
      error = function(e) NULL
    )
    if (!is.null(em_fit) && is.finite(em_fit$loglik)) {
      r_post <- em_fit$r                                  # over positives; zeros = 0
      assignments <- fit$pos_idx[r_post >= probability_threshold]
      prob_quantiles <- stats::quantile(r_post, probs = prob_quantile_probs)
      n_above_prob_thresh <- length(assignments)
      eta0 <- pmin(pmax(em_fit$b0 + em_fit$b2 * fit$x_pos, -ETA_CAP), ETA_CAP)
      eta1 <- pmin(pmax(em_fit$b0 + em_fit$b1 + em_fit$b2 * fit$x_pos, -ETA_CAP), ETA_CAP)
      poismix <- list(
        pi = em_fit$pi, b0 = em_fit$b0, b1 = em_fit$b1, b2 = em_fit$b2,
        mean_mu_np = mean(exp(eta0)), mean_mu_pert = mean(exp(eta1)),
        mean_r = mean(r_post), n_fit = length(g)
      )
      em_converged <- em_fit$converged
      em_log_lik   <- em_fit$loglik
      em_init_i    <- em_fit$start_i
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
    n_fit                = length(g),
    n_assigned           = length(assignments),
    em_converged         = em_converged,
    em_log_lik           = em_log_lik,
    em_init_i            = em_init_i,
    poismix              = poismix,
    prob_quantiles       = prob_quantiles,
    n_above_prob_thresh  = n_above_prob_thresh,
    offset_model_summary = NULL,
    elapsed_sec          = elapsed_sec
  )
  if (keep_fits) out$em_fit <- em_fit
  out
}

poismix0_assign <- function(grna_matrix,
                            grna_n_nonzero,
                            grna_ids              = rownames(grna_matrix),
                            n_em_rep              = 5L,
                            pi_guess_range        = c(1e-5, 0.1),
                            g_pert_guess_range    = log(c(10, 5000)),
                            n_fit_cutoff          = 10L,
                            backup_threshold      = 5L,
                            probability_threshold = 0.8,
                            max_iter              = 200L,
                            tol                   = 1e-8,
                            eps                   = 1e-10,
                            ridge                 = 1e-6,
                            seed                  = 4L,
                            cl                    = NULL,
                            keep_fits             = FALSE) {
  force(grna_matrix); force(grna_n_nonzero)
  force(n_fit_cutoff); force(backup_threshold); force(probability_threshold)
  force(max_iter); force(tol); force(eps); force(ridge); force(keep_fits)

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  set.seed(seed)
  pi_guesses     <- stats::runif(n_em_rep, pi_guess_range[1],     pi_guess_range[2])
  g_pert_guesses <- stats::runif(n_em_rep, g_pert_guess_range[1], g_pert_guess_range[2])
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
  force(pi_guesses); force(g_pert_guesses)

  worker <- function(grna_id) {
    g <- as.numeric(grna_matrix[grna_id, ])
    assign_one_grna_poismix0(
      g = g, grna_n_nonzero = grna_n_nonzero,
      pi_guesses = pi_guesses, g_pert_guesses = g_pert_guesses,
      n_fit_cutoff = n_fit_cutoff, backup_threshold = backup_threshold,
      probability_threshold = probability_threshold, max_iter = max_iter,
      tol = tol, eps = eps, ridge = ridge, keep_fits = keep_fits
    )
  }

  per_guide <- if (is.null(cl)) lapply(grna_ids, worker) else parallel::parLapplyLB(cl, grna_ids, worker)
  per_guide <- stats::setNames(per_guide, as.character(grna_ids))

  run_meta <- list(
    model                 = "pois-pois-joint-nnz-pois0",
    covariate             = "log1p(grna_n_nonzero - 1{y>0})",
    response              = "y (positives fit mixture; y==0 forced non-pert)",
    n_params              = 4L,
    n_em_rep              = n_em_rep,
    pi_guess_range        = pi_guess_range,
    g_pert_guess_range    = g_pert_guess_range,
    n_fit_cutoff          = n_fit_cutoff,
    backup_threshold      = backup_threshold,
    probability_threshold = probability_threshold,
    max_iter              = max_iter,
    tol                   = tol,
    ridge                 = ridge,
    seed                  = seed
  )

  list(per_guide = per_guide, run_meta = run_meta)
}
