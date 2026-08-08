# Gated logistic / Gamma-Gaussian mixture guide-assignment model.
#
# This is structurally different from sceptre_assign_pure_R(): there is NO
# offset model. The single covariate is grna_n_nonzero, entering through a
# logistic *gate* on the (cell-specific) prior probability of perturbation.
#
# Per guide, restricted to cells with y_i >= 2 (the "fit cells"):
#   x_i  = log(grna_n_nonzero_i)                       # decircularized covariate
#                                                       # (the -1 from this guide
#                                                       #  is ignored for now)
#   z_i  = log(y_i)
#   pi_i = ilogit(a0 + a1 * x_i)                        # P(cell i perturbed)
#   z_i ~ (1 - pi_i) * Gamma(alpha, beta)  +  pi_i * Normal(mu, sigma)
#          \________ non-pert (low) _______/    \__ pert (high) __/
#
# Six parameters: a0, a1 (gate), alpha, beta (Gamma rate-param), mu, sigma.
#
# Cells with y_i in {0, 1} are not fit and are assigned non-pert. A cell is
# assigned perturbed iff its posterior responsibility
#   gamma_i = P(pert_i | z_i, x_i)
#           = pi_i N(z_i) / [ (1 - pi_i) Gamma(z_i) + pi_i N(z_i) ]
# is >= probability_threshold (default 0.8). This posterior plays the exact role
# the old mixture posterior did; pi_i is just cell-specific now.
#
# Identifiability: the Gaussian is forced to be the HIGH component (pert) via a
# constrained M-step -- if a step yields mu < Gamma mean (alpha/beta), mu is
# clamped up to (alpha/beta + mu_margin).

# ---- helpers -----------------------------------------------------------------

# Weighted Gamma shape MLE: solve log(alpha) - digamma(alpha) = s for alpha,
# where s = log(weighted mean of z) - weighted mean of log(z) >= 0. Closed-form
# seed (Minka / Choi-Wette) followed by a few Newton steps.
gated_gamma_shape_mle <- function(s) {
  s <- max(s, 1e-8)
  a <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
  for (i in seq_len(10L)) {
    f  <- log(a) - digamma(a) - s
    fp <- 1 / a - trigamma(a)
    a_new <- a - f / fp
    if (!is.finite(a_new) || a_new <= 0) a_new <- a / 2
    if (abs(a_new - a) < 1e-8) { a <- a_new; break }
    a <- a_new
  }
  a
}

# Weighted (fractional-response) logistic regression of soft labels `y` (the
# responsibilities) on design [1, x]. Hand-rolled 2-parameter IRLS so there are
# no glm warnings about non-integer successes and no per-call overhead.
gated_irls_logit <- function(x, y, n_iter = 25L, tol = 1e-8) {
  a0 <- 0; a1 <- 0
  for (i in seq_len(n_iter)) {
    eta <- a0 + a1 * x
    p   <- stats::plogis(eta)
    w   <- pmax(p * (1 - p), 1e-8)
    u   <- eta + (y - p) / w                 # IRLS working response
    s0  <- sum(w);        sx  <- sum(w * x)
    sxx <- sum(w * x * x); b0 <- sum(w * u);  b1 <- sum(w * x * u)
    det <- s0 * sxx - sx * sx
    if (!is.finite(det) || abs(det) < 1e-12) break
    a0_new <- (sxx * b0 - sx * b1) / det
    a1_new <- (s0  * b1 - sx * b0) / det
    if (!is.finite(a0_new) || !is.finite(a1_new)) break
    d <- abs(a0_new - a0) + abs(a1_new - a1)
    a0 <- a0_new; a1 <- a1_new
    if (d < tol) break
  }
  c(a0, a1)
}

# One M-step. `gamma` = current responsibilities (pert). Returns the 6 params.
gated_m_step <- function(z, x, gamma, mu_margin = 1e-3, sigma_floor = 0.05) {
  w1 <- gamma; w0 <- 1 - gamma
  W1 <- sum(w1); W0 <- sum(w0)

  # Gaussian (pert / high component)
  if (W1 > 1e-8) {
    mu    <- sum(w1 * z) / W1
    sig2  <- sum(w1 * (z - mu)^2) / W1
  } else {
    mu    <- max(z)
    sig2  <- stats::var(z)
  }
  sigma <- sqrt(max(sig2, sigma_floor^2))

  # Gamma (non-pert / low component)
  if (W0 > 1e-8) {
    zbar    <- sum(w0 * z) / W0
    meanlog <- sum(w0 * log(z)) / W0
    s       <- log(zbar) - meanlog
    alpha   <- gated_gamma_shape_mle(s)
    beta    <- alpha / zbar
  } else {
    alpha <- 1; beta <- 1 / mean(z)
  }

  # hard constraint: Gaussian mean must exceed Gamma mean (pert = high).
  gmean <- alpha / beta
  if (!is.finite(mu) || mu < gmean + mu_margin) mu <- gmean + mu_margin

  ab <- gated_irls_logit(x, gamma)
  list(a0 = ab[1L], a1 = ab[2L], alpha = alpha, beta = beta, mu = mu, sigma = sigma)
}

# One E-step. Returns responsibilities + total log-likelihood (log-sum-exp).
gated_e_step <- function(z, x, par) {
  eta      <- par$a0 + par$a1 * x
  log_pi   <- stats::plogis(eta, log.p = TRUE)                    # log pi
  log_1mpi <- stats::plogis(eta, lower.tail = FALSE, log.p = TRUE) # log(1 - pi)
  la <- log_1mpi + stats::dgamma(z, shape = par$alpha, rate = par$beta, log = TRUE)
  lb <- log_pi   + stats::dnorm(z, mean = par$mu, sd = par$sigma, log = TRUE)
  m  <- pmax(la, lb)
  ll_cell <- m + log(exp(la - m) + exp(lb - m))
  gamma   <- exp(lb - ll_cell)
  gamma[!is.finite(gamma)] <- 0
  list(gamma = gamma, ll = sum(ll_cell))
}

# Run EM from a single initialization (responsibilities `gamma_init`).
gated_em_one_start <- function(z, x, gamma_init, max_iter = 200L,
                               tol = 1e-6, mu_margin = 1e-3, sigma_floor = 0.05) {
  gamma     <- gamma_init
  ll_old    <- -Inf
  converged <- FALSE
  par       <- NULL
  iter      <- 0L
  for (iter in seq_len(max_iter)) {
    par <- gated_m_step(z, x, gamma, mu_margin, sigma_floor)
    es  <- gated_e_step(z, x, par)
    gamma <- es$gamma
    ll    <- es$ll
    if (is.finite(ll) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) {
      converged <- TRUE; ll_old <- ll; break
    }
    ll_old <- ll
  }
  list(par = par, gamma = gamma, loglik = ll_old,
       converged = converged, iter = iter)
}

# Build an initial parameter set from a sceptre-style random seed:
#   pi0      -- initial (constant) prior perturbation probability -> gate intercept
#   mu_pert0 -- initial perturbed-component mean on the log(y) scale
# The Gamma (non-pert) is seeded by method-of-moments on the cells below the
# seeded perturbed mean (the presumed background); the gate slope starts at 0.
gated_init_par_from_seed <- function(z, pi0, mu_pert0, mu_margin = 1e-3,
                                     sigma_floor = 0.05) {
  a0 <- stats::qlogis(pi0)
  a1 <- 0
  lo <- z[z < mu_pert0]
  if (length(lo) < 5L) lo <- z
  m <- mean(lo); v <- max(stats::var(lo), 1e-4)
  alpha <- max(m * m / v, 0.05)
  beta  <- max(m / v, 0.05)
  mu    <- mu_pert0
  sigma <- max(stats::sd(z), sigma_floor)
  gmean <- alpha / beta
  if (!is.finite(mu) || mu < gmean + mu_margin) mu <- gmean + mu_margin
  list(a0 = a0, a1 = a1, alpha = alpha, beta = beta, mu = mu, sigma = sigma)
}

# Multi-start EM, sceptre-style. Each start is a random (pi, mu_pert) seed
# (`pi_guesses` / `mu_pert_guesses`, generated once by the driver under an
# isolated seed and shared across guides). Each start runs EM to convergence;
# selection mirrors run_em_pois_pois_pure_R's `converged && loglik > best`: a
# converged start always beats a non-converged one, ties break on log-lik.
run_em_gated_gamma_gauss <- function(z, x, pi_guesses, mu_pert_guesses,
                                     max_iter = 200L, tol = 1e-6,
                                     mu_margin = 1e-3, sigma_floor = 0.05) {
  B <- length(pi_guesses)
  stopifnot(length(mu_pert_guesses) == B)
  best <- NULL; best_i <- NA_integer_
  better <- function(cand, cur) {
    if (is.null(cur)) return(TRUE)
    if (cand$converged != cur$converged) return(cand$converged)   # converged wins
    is.finite(cand$loglik) &&
      (!is.finite(cur$loglik) || cand$loglik > cur$loglik)
  }
  for (i in seq_len(B)) {
    par0   <- gated_init_par_from_seed(z, pi_guesses[i], mu_pert_guesses[i],
                                       mu_margin, sigma_floor)
    gamma0 <- gated_e_step(z, x, par0)$gamma
    fit    <- gated_em_one_start(z, x, gamma0, max_iter, tol, mu_margin, sigma_floor)
    if (better(fit, best)) { best <- fit; best_i <- i }
  }
  list(
    a0 = best$par$a0, a1 = best$par$a1,
    alpha = best$par$alpha, beta = best$par$beta,
    mu = best$par$mu, sigma = best$par$sigma,
    gamma = best$gamma, loglik = best$loglik,
    converged = best$converged, iter = best$iter, start_i = best_i
  )
}

# ---- per-guide -----------------------------------------------------------------

assign_one_grna_gated <- function(g, grna_n_nonzero,
                                  pi_guesses, mu_pert_guesses,
                                  n_fit_cutoff          = 10L,
                                  backup_threshold      = 5L,
                                  probability_threshold = 0.8,
                                  max_iter              = 200L,
                                  tol                   = 1e-6,
                                  mu_margin             = 1e-3,
                                  sigma_floor           = 0.05,
                                  keep_fits             = FALSE) {
  t_start   <- Sys.time()
  n_nonzero <- sum(g >= 1)
  fit_idx   <- which(g >= 2L)               # continuous model fit on y >= 2
  n_fit     <- length(fit_idx)

  # defaults for the below-cutoff / failure paths
  em_fit       <- NULL
  em_converged <- NA
  em_log_lik   <- NA_real_
  em_init_i    <- NA_integer_
  gated        <- list(a0 = NA_real_, a1 = NA_real_, alpha = NA_real_,
                       beta = NA_real_, mu = NA_real_, sigma = NA_real_,
                       mean_pi = NA_real_, mean_gamma = NA_real_, n_fit = n_fit)
  prob_quantile_probs <- c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
  prob_quantiles <- stats::setNames(rep(NA_real_, length(prob_quantile_probs)),
                                    paste0(prob_quantile_probs * 100, "%"))
  n_above_prob_thresh <- NA_integer_

  if (n_fit >= n_fit_cutoff) {
    z <- log(g[fit_idx])
    x <- log(grna_n_nonzero[fit_idx])
    em_fit <- tryCatch(
      run_em_gated_gamma_gauss(z, x, pi_guesses, mu_pert_guesses,
                               max_iter, tol, mu_margin, sigma_floor),
      error = function(e) NULL
    )
    if (!is.null(em_fit) && is.finite(em_fit$loglik)) {
      gamma_post  <- em_fit$gamma
      assignments <- fit_idx[gamma_post >= probability_threshold]
      prob_quantiles <- stats::quantile(gamma_post, probs = prob_quantile_probs)
      n_above_prob_thresh <- length(assignments)
      eta_fit <- em_fit$a0 + em_fit$a1 * x
      gated <- list(
        a0 = em_fit$a0, a1 = em_fit$a1,
        alpha = em_fit$alpha, beta = em_fit$beta,
        mu = em_fit$mu, sigma = em_fit$sigma,
        mean_pi = mean(stats::plogis(eta_fit)),
        mean_gamma = mean(gamma_post),
        n_fit = n_fit
      )
      em_converged <- em_fit$converged
      em_log_lik   <- em_fit$loglik
      em_init_i    <- em_fit$start_i
    } else {
      # EM failed -> fall back to a hard threshold
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
    gated                = gated,
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
# each cell); the same vector is used for every guide (the guide-specific -1
# decircularization is ignored for now).
gated_assign <- function(grna_matrix,
                         grna_n_nonzero,
                         grna_ids              = rownames(grna_matrix),
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
                         cl                    = NULL,
                         keep_fits             = FALSE) {
  force(grna_matrix); force(grna_n_nonzero)
  force(n_fit_cutoff); force(backup_threshold); force(probability_threshold)
  force(max_iter); force(tol); force(mu_margin)
  force(sigma_floor); force(keep_fits)

  # Starting guesses: generated ONCE here under an isolated seed and shared
  # across all guides, exactly like sceptre_assign_pure_R's driver. The global
  # RNG state is saved/restored so callers' streams aren't perturbed.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  pi_guesses      <- stats::runif(n_em_rep, pi_guess_range[1],      pi_guess_range[2])
  mu_pert_guesses <- stats::runif(n_em_rep, mu_pert_guess_range[1], mu_pert_guess_range[2])
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
  force(pi_guesses); force(mu_pert_guesses)

  worker <- function(grna_id) {
    g <- as.numeric(grna_matrix[grna_id, ])
    assign_one_grna_gated(
      g                     = g,
      grna_n_nonzero        = grna_n_nonzero,
      pi_guesses            = pi_guesses,
      mu_pert_guesses       = mu_pert_guesses,
      n_fit_cutoff          = n_fit_cutoff,
      backup_threshold      = backup_threshold,
      probability_threshold = probability_threshold,
      max_iter              = max_iter,
      tol                   = tol,
      mu_margin             = mu_margin,
      sigma_floor           = sigma_floor,
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
    model                 = "gated-gamma-gauss",
    gate_covariate        = "log(grna_n_nonzero)",
    response              = "log(y), y >= 2",
    n_params              = 6L,
    y_min_fit             = 2L,
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
    seed                  = seed
  )

  list(per_guide = per_guide, run_meta = run_meta)
}
