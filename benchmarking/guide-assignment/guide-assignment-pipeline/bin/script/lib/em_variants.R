# Experimental EM-algorithm variants for gRNA mixture-model assignment.
#
# Companion to sceptre_assign_pure_R.R, which holds the offset-model fits, the
# per-guide driver, and the "production" EM kernels (run_em_*_pure_R). New /
# experimental EM kernels live here so that file stays manageable. Both files
# are sourced by run_variant.R in the main process AND on every PSOCK worker, so
# any kernel defined here is visible wherever the production kernels are.
#
# Contract note: the kernels in sceptre_assign_pure_R.R return the per-guide EM
# trajectory under fixed names (outer_Ti1s, outer_pi, outer_g_pert,
# outer_log_lik, outer_converged, outer_i, trajectory) that the driver's family
# dispatch (assign_one_grna_pure_R) consumes. Prototypes here may use their own
# richer return shape while being developed; adapt them to the outer_* contract
# before wiring into the family dispatch.


# ---- Additive-signal Poisson mixture EM (single fit) -------------------------
#
# A CONSISTENT additive-signal Poisson mixture. Each cell is
#
#   y_i ~ (1 - pi) Pois(mu0_i)  +  pi Pois(mu0_i + lambda),   mu0_i = exp(offset_i),
#
# i.e. the perturbed component adds an independent Poisson signal with mean
# lambda = exp(gamma) on top of the baseline mean mu0_i (additive on the count
# scale). Both the E-step and the M-step use this same additive mean, so unlike
# sceptre's pois kernel this model is internally consistent: pi and gamma are a
# proper joint EM fit. lambda's M-step is a 1-D score-root solve (uniroot) --
# the score sum_i w_i [y_i / (mu0_i + lambda) - 1] is strictly decreasing in
# lambda, so the bracket-and-root is well defined; there is no closed form
# because mu0_i varies across cells.
#
# use_log_gamma = TRUE: when CLASSIFYING cells, use gamma itself -- not
# exp(gamma) -- as the additive signal mean (mu1 = mu0 + gamma). This MIMICS THE
# LOG-SCALE PATHOLOGY in sceptre's bugged pois kernel, where a log-scale effect
# is plugged into the perturbed mean additively, collapsing it to just above
# baseline. It is NOT a faithful replication of sceptre's assignments: sceptre's
# gamma is the multiplicative-MLE update log(sum w*y / sum w*mu0), whereas this
# gamma is log(lambda) for the correctly-fit additive bump -- different
# quantities. For exact sceptre-bug behaviour use
# run_em_pois_pois_pure_R(fix_curr_g_pert_bug = FALSE) instead. The fit-based
# posterior (using exp(gamma)) is always returned alongside, so the two
# classifications can be compared.
fit_pois_additive_signal_em <- function(
  y,
  offset,
  pi_init = 0.001,
  gamma_init = NULL,
  max_iter = 200,
  tol = 1e-8,
  min_pi = 1e-10,
  min_signal_mean = 1e-10,
  max_signal_mean = 1e8,
  use_log_gamma = FALSE,
  classify_threshold = 0.5,
  verbose = FALSE
) {
  stopifnot(length(y) == length(offset))
  stopifnot(all(y >= 0, na.rm = TRUE))
  stopifnot(all(is.finite(y)))
  stopifnot(all(is.finite(offset)))

  y <- as.numeric(y)
  offset <- as.numeric(offset)

  n <- length(y)
  mu0 <- exp(offset)

  logsumexp2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }

  # Initialize signal mean lambda = exp(gamma).
  if (is.null(gamma_init)) {
    # Crude excess-count initialization.
    lambda_init <- mean(pmax(y - mu0, 0))
    lambda_init <- max(lambda_init, min_signal_mean)
    gamma <- log(lambda_init)
  } else {
    gamma <- gamma_init
  }

  signal_mean <- function(gamma) {
    exp(gamma)
  }

  pi <- min(max(pi_init, min_pi), 1 - min_pi)
  lambda <- min(max(signal_mean(gamma), min_signal_mean), max_signal_mean)
  gamma <- log(lambda)

  compute_posterior <- function(pi, signal_mean_for_pert) {
    mu1 <- mu0 + signal_mean_for_pert

    log_f0 <- dpois(y, lambda = mu0, log = TRUE)
    log_f1 <- dpois(y, lambda = mu1, log = TRUE)

    log_num0 <- log1p(-pi) + log_f0
    log_num1 <- log(pi) + log_f1
    log_den <- logsumexp2(log_num0, log_num1)

    exp(log_num1 - log_den)
  }

  compute_loglik <- function(pi, signal_mean_for_pert) {
    mu1 <- mu0 + signal_mean_for_pert

    log_f0 <- dpois(y, lambda = mu0, log = TRUE)
    log_f1 <- dpois(y, lambda = mu1, log = TRUE)

    sum(logsumexp2(
      log1p(-pi) + log_f0,
      log(pi) + log_f1
    ))
  }

  update_lambda <- function(w) {
    # Maximize:
    # sum_i w_i [ y_i log(mu0_i + lambda) - (mu0_i + lambda) ]
    #
    # Score:
    # sum_i w_i [ y_i / (mu0_i + lambda) - 1 ].

    score <- function(lambda) {
      sum(w * (y / (mu0 + lambda) - 1))
    }

    lower <- min_signal_mean

    if (score(lower) <= 0) {
      return(lower)
    }

    upper <- max(1, 2 * lower)
    while (score(upper) > 0 && upper < max_signal_mean) {
      upper <- upper * 2
    }

    upper <- min(upper, max_signal_mean)

    if (score(upper) > 0) {
      return(max_signal_mean)
    }

    uniroot(score, lower = lower, upper = upper)$root
  }

  history <- data.frame(
    iter = integer(),
    loglik = numeric(),
    pi = numeric(),
    gamma = numeric(),
    signal_mean = numeric(),
    n_eff = numeric()
  )

  old_loglik <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # E-step under the actual fitted model: signal mean = exp(gamma)
    lambda <- exp(gamma)
    w <- compute_posterior(pi, signal_mean_for_pert = lambda)

    # M-step for pi
    pi <- mean(w)
    pi <- min(max(pi, min_pi), 1 - min_pi)

    # M-step for gamma through lambda = exp(gamma)
    lambda <- update_lambda(w)
    gamma <- log(lambda)

    loglik <- compute_loglik(pi, signal_mean_for_pert = lambda)

    history <- rbind(
      history,
      data.frame(
        iter = iter,
        loglik = loglik,
        pi = pi,
        gamma = gamma,
        signal_mean = lambda,
        n_eff = sum(w)
      )
    )

    if (verbose) {
      message(sprintf(
        "iter %03d | loglik %.6f | pi %.4g | gamma %.4f | signal_mean %.4f | n_eff %.2f",
        iter, loglik, pi, gamma, lambda, sum(w)
      ))
    }

    if (is.finite(old_loglik)) {
      if (abs(loglik - old_loglik) < tol * (1 + abs(old_loglik))) {
        converged <- TRUE
        break
      }
    }

    old_loglik <- loglik
  }

  # Standard fitted posterior: uses exp(gamma) as signal mean.
  signal_mean_fit <- exp(gamma)
  posterior_fit <- compute_posterior(
    pi,
    signal_mean_for_pert = signal_mean_fit
  )

  # Optional bug-style / log-gamma classification:
  # use gamma itself as the additive signal mean.
  if (use_log_gamma) {
    signal_mean_classify <- gamma

    if (signal_mean_classify <= 0) {
      warning(
        "use_log_gamma = TRUE but fitted gamma <= 0; ",
        "clipping classification signal mean to min_signal_mean."
      )
      signal_mean_classify <- min_signal_mean
    }
  } else {
    signal_mean_classify <- signal_mean_fit
  }

  posterior_classification <- compute_posterior(
    pi,
    signal_mean_for_pert = signal_mean_classify
  )

  called <- posterior_classification > classify_threshold

  list(
    pi = pi,
    gamma = gamma,
    signal_mean_fit = signal_mean_fit,
    signal_mean_classification = signal_mean_classify,
    mu0 = mu0,
    mu1_fit = mu0 + signal_mean_fit,
    mu1_classification = mu0 + signal_mean_classify,
    posterior_fit = posterior_fit,
    posterior_classification = posterior_classification,
    called = called,
    loglik = tail(history$loglik, 1),
    history = history,
    converged = converged,
    n_iter = nrow(history),
    use_log_gamma = use_log_gamma
  )
}


# ---- Additive NB-signal Poisson mixture EM (single fit) ----------------------
#
# Like fit_pois_additive_signal_em, but the perturbed signal is NEGATIVE
# BINOMIAL rather than Poisson. Each cell is
#
#   y_i ~ (1 - pi) Pois(mu0_i)  +  pi [ Pois(mu0_i) + NB(mu = exp(gamma), size = theta) ],
#
# i.e. a perturbed cell is the independent SUM of the baseline Poisson and an NB
# signal with mean exp(gamma) and dispersion (size) theta. The perturbed density
# is the Poisson-NB convolution P(Y = y) = sum_{s=0}^{y} NB(s) Pois(y - s),
# evaluated exactly via log_pois_plus_nb (O(y_i) terms per cell). theta -> Inf
# recovers the Poisson-signal model (poissum).
#
# EM: E-step posterior as usual; pi M-step = mean(w); (gamma, theta) M-step
# maximizes the weighted convolution log-likelihood sum_i w_i log f1_i jointly
# via L-BFGS-B on (log signal_mean, log theta). The baseline Pois(mu0_i) has no
# free parameters (mu0 is fixed by the offset).
#
# use_log_gamma = TRUE: classify using gamma itself -- not exp(gamma) -- as the
# NB signal mean, the same log-scale pathology mimic as
# fit_pois_additive_signal_em.
#
# NOTE: the exact convolution inside a numerical (gamma, theta) optimization on
# every EM iteration is considerably slower than poissum; expect it to be costly
# for guides with large counts.
fit_pois_additive_nb_signal_em <- function(
  y,
  offset,
  pi_init = 0.001,
  gamma_init = NULL,
  theta_init = 1,
  max_iter = 200,
  tol = 1e-8,
  min_pi = 1e-10,
  min_signal_mean = 1e-10,
  max_signal_mean = 1e8,
  min_theta = 1e-4,
  max_theta = 1e6,
  use_log_gamma = FALSE,
  classify_threshold = 0.5,
  verbose = FALSE
) {
  stopifnot(length(y) == length(offset))
  stopifnot(all(y >= 0, na.rm = TRUE))
  stopifnot(all(is.finite(y)))
  stopifnot(all(is.finite(offset)))
  stopifnot(all(y == floor(y)))

  y <- as.integer(y)
  offset <- as.numeric(offset)

  n <- length(y)
  mu0 <- exp(offset)

  logsumexp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }

  logsumexp2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }

  log_pois_plus_nb_one <- function(y_i, mu0_i, signal_mean, theta) {
    s <- 0:y_i

    log_terms <- stats::dpois(
      x = y_i - s,
      lambda = mu0_i,
      log = TRUE
    ) +
      stats::dnbinom(
        x = s,
        size = theta,
        mu = signal_mean,
        log = TRUE
      )

    logsumexp(log_terms)
  }

  log_pois_plus_nb <- function(signal_mean, theta) {
    vapply(
      seq_along(y),
      function(i) {
        log_pois_plus_nb_one(
          y_i = y[i],
          mu0_i = mu0[i],
          signal_mean = signal_mean,
          theta = theta
        )
      },
      numeric(1)
    )
  }

  if (is.null(gamma_init)) {
    signal_mean_init <- mean(pmax(y - mu0, 0))
    signal_mean_init <- max(signal_mean_init, min_signal_mean)
    gamma <- log(signal_mean_init)
  } else {
    gamma <- gamma_init
  }

  signal_mean <- min(max(exp(gamma), min_signal_mean), max_signal_mean)
  gamma <- log(signal_mean)

  theta <- min(max(theta_init, min_theta), max_theta)
  pi <- min(max(pi_init, min_pi), 1 - min_pi)

  compute_posterior <- function(pi, signal_mean_for_pert, theta) {
    log_f0 <- stats::dpois(
      x = y,
      lambda = mu0,
      log = TRUE
    )

    log_f1 <- log_pois_plus_nb(
      signal_mean = signal_mean_for_pert,
      theta = theta
    )

    log_num0 <- log1p(-pi) + log_f0
    log_num1 <- log(pi) + log_f1
    log_den <- logsumexp2(log_num0, log_num1)

    exp(log_num1 - log_den)
  }

  compute_loglik <- function(pi, signal_mean_for_pert, theta) {
    log_f0 <- stats::dpois(
      x = y,
      lambda = mu0,
      log = TRUE
    )

    log_f1 <- log_pois_plus_nb(
      signal_mean = signal_mean_for_pert,
      theta = theta
    )

    sum(logsumexp2(
      log1p(-pi) + log_f0,
      log(pi) + log_f1
    ))
  }

  history <- data.frame(
    iter = integer(),
    loglik = numeric(),
    pi = numeric(),
    gamma = numeric(),
    signal_mean = numeric(),
    theta = numeric(),
    n_eff = numeric()
  )

  old_loglik <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    signal_mean <- exp(gamma)

    # E-step
    w <- compute_posterior(
      pi = pi,
      signal_mean_for_pert = signal_mean,
      theta = theta
    )

    # M-step for pi
    pi <- mean(w)
    pi <- min(max(pi, min_pi), 1 - min_pi)

    # M-step for gamma and theta.
    # Maximize weighted log convolution likelihood:
    # sum_i w_i log P(Pois(mu0_i) + NB(exp(gamma), theta) = y_i)
    objective <- function(par) {
      signal_mean_curr <- exp(par[1])
      theta_curr <- exp(par[2])

      log_f1 <- log_pois_plus_nb(
        signal_mean = signal_mean_curr,
        theta = theta_curr
      )

      -sum(w * log_f1)
    }

    opt <- optim(
      par = c(log(signal_mean), log(theta)),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(log(min_signal_mean), log(min_theta)),
      upper = c(log(max_signal_mean), log(max_theta)),
      control = list(maxit = 100)
    )

    gamma <- opt$par[1]
    signal_mean <- exp(gamma)
    theta <- exp(opt$par[2])

    loglik <- compute_loglik(
      pi = pi,
      signal_mean_for_pert = signal_mean,
      theta = theta
    )

    history <- rbind(
      history,
      data.frame(
        iter = iter,
        loglik = loglik,
        pi = pi,
        gamma = gamma,
        signal_mean = signal_mean,
        theta = theta,
        n_eff = sum(w)
      )
    )

    if (verbose) {
      message(sprintf(
        paste0(
          "iter %03d | loglik %.6f | pi %.4g | ",
          "gamma %.4f | signal_mean %.4f | theta %.4f | n_eff %.2f"
        ),
        iter, loglik, pi, gamma, signal_mean, theta, sum(w)
      ))
    }

    if (is.finite(old_loglik)) {
      if (abs(loglik - old_loglik) < tol * (1 + abs(old_loglik))) {
        converged <- TRUE
        break
      }
    }

    old_loglik <- loglik
  }

  # Standard fitted posterior: uses exp(gamma) as NB signal mean.
  signal_mean_fit <- exp(gamma)

  posterior_fit <- compute_posterior(
    pi = pi,
    signal_mean_for_pert = signal_mean_fit,
    theta = theta
  )

  # Optional bug-style final classification:
  # use gamma itself as the NB signal mean.
  if (use_log_gamma) {
    signal_mean_classification <- gamma

    if (signal_mean_classification <= 0) {
      warning(
        "use_log_gamma = TRUE but fitted gamma <= 0; ",
        "clipping classification signal mean to min_signal_mean."
      )
      signal_mean_classification <- min_signal_mean
    }
  } else {
    signal_mean_classification <- signal_mean_fit
  }

  posterior_classification <- compute_posterior(
    pi = pi,
    signal_mean_for_pert = signal_mean_classification,
    theta = theta
  )

  called <- posterior_classification > classify_threshold

  list(
    pi = pi,
    gamma = gamma,
    theta = theta,
    signal_mean_fit = signal_mean_fit,
    signal_mean_classification = signal_mean_classification,
    mu0 = mu0,
    posterior_fit = posterior_fit,
    posterior_classification = posterior_classification,
    called = called,
    loglik = tail(history$loglik, 1),
    history = history,
    converged = converged,
    n_iter = nrow(history),
    use_log_gamma = use_log_gamma
  )
}
