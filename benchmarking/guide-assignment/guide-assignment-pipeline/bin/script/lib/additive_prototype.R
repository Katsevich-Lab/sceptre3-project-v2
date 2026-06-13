fit_additive_pois_mix <- function(y,
                                  mu0,
                                  init_pi = NULL,
                                  init_delta = NULL,
                                  max_iter = 200L,
                                  tol = 1e-8,
                                  min_delta = 1e-10,
                                  pi_bounds = c(1e-6, 1 - 1e-6),
                                  verbose = FALSE) {
  # y: guide UMI count vector for one guide across cells
  # mu0: background/null mean for each cell, on count scale
  #      e.g. mu0 = exp(offset), or your existing g_mus_pert0
  
  y <- as.numeric(y)
  mu0 <- as.numeric(mu0)
  
  if (length(y) != length(mu0)) {
    stop("y and mu0 must have the same length.")
  }
  if (any(!is.finite(y)) || any(!is.finite(mu0))) {
    stop("y and mu0 must be finite.")
  }
  if (any(y < 0) || any(mu0 < 0)) {
    stop("y and mu0 must be nonnegative.")
  }
  
  # Avoid exact zero background means for log/Poisson numerical stability.
  mu0 <- pmax(mu0, 1e-300)
  
  n <- length(y)
  
  clamp <- function(x, lo, hi) {
    pmin(pmax(x, lo), hi)
  }
  
  logspace2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  # Crude initialization.
  if (is.null(init_pi)) {
    # Cells that look surprising under the background model.
    p_tail <- ppois(y - 1, lambda = mu0, lower.tail = FALSE)
    init_pi <- mean(y > 0 & p_tail < 0.01)
    
    if (!is.finite(init_pi) || init_pi <= 0) {
      init_pi <- min(max(mean(y > 0) / 2, 0.01), 0.25)
    }
  }
  
  pi <- clamp(init_pi, pi_bounds[1], pi_bounds[2])
  
  if (is.null(init_delta)) {
    idx <- y > 0
    
    if (any(idx)) {
      init_delta <- stats::quantile(y[idx], probs = 0.9,
                                    names = FALSE, type = 1) -
        stats::median(mu0[idx])
    } else {
      init_delta <- min_delta
    }
    
    init_delta <- max(as.numeric(init_delta), min_delta, 0.1)
  }
  
  delta <- max(as.numeric(init_delta), min_delta)
  
  mstep_delta <- function(r, y, mu0, old_delta) {
    if (sum(r) <= 0) {
      return(max(old_delta, min_delta))
    }
    
    score <- function(d) {
      sum(r * (y / (mu0 + d) - 1))
    }
    
    # Boundary optimum near zero.
    s_lo <- score(min_delta)
    if (!is.finite(s_lo) || s_lo <= 0) {
      return(min_delta)
    }
    
    # Find an upper bracket where the score is negative.
    hi <- max(old_delta * 2, max(y), mean(y) + mean(mu0), 1)
    hi <- max(hi, min_delta * 2)
    
    s_hi <- score(hi)
    k <- 0L
    
    while (is.finite(s_hi) && s_hi > 0 && k < 100L) {
      hi <- hi * 2
      s_hi <- score(hi)
      k <- k + 1L
    }
    
    if (!is.finite(s_hi) || s_hi > 0) {
      # Very rare fallback: Q is still increasing over the searched range.
      return(hi)
    }
    
    stats::uniroot(score, lower = min_delta, upper = hi, tol = tol)$root
  }
  
  trace <- matrix(NA_real_, nrow = max_iter, ncol = 4)
  colnames(trace) <- c("iter", "loglik", "pi", "delta")
  
  loglik_old <- -Inf
  converged <- FALSE
  r <- rep(pi, n)
  
  for (iter in seq_len(max_iter)) {
    # E-step.
    log0 <- log1p(-pi) + stats::dpois(y, lambda = mu0, log = TRUE)
    log1 <- log(pi) + stats::dpois(y, lambda = mu0 + delta, log = TRUE)
    
    logden <- logspace2(log0, log1)
    r <- exp(log1 - logden)
    
    loglik <- sum(logden)
    
    trace[iter, ] <- c(iter, loglik, pi, delta)
    
    if (verbose) {
      message(sprintf(
        "iter %03d  loglik %.6f  pi %.5g  delta %.5g  eta %.5g",
        iter, loglik, pi, delta, log(delta)
      ))
    }
    
    # M-step.
    pi_new <- clamp(mean(r), pi_bounds[1], pi_bounds[2])
    delta_new <- mstep_delta(r, y, mu0, delta)
    
    if (iter > 1L) {
      ll_change <- abs(loglik - loglik_old) / (abs(loglik_old) + 1)
      par_change <- max(
        abs(pi_new - pi),
        abs(log(delta_new) - log(delta))
      )
      
      if (ll_change < tol && par_change < sqrt(tol)) {
        pi <- pi_new
        delta <- delta_new
        converged <- TRUE
        break
      }
    }
    
    pi <- pi_new
    delta <- delta_new
    loglik_old <- loglik
  }
  
  # Final E-step using final parameters.
  log0 <- log1p(-pi) + stats::dpois(y, lambda = mu0, log = TRUE)
  log1 <- log(pi) + stats::dpois(y, lambda = mu0 + delta, log = TRUE)
  logden <- logspace2(log0, log1)
  r <- exp(log1 - logden)
  loglik <- sum(logden)
  
  trace <- as.data.frame(trace[seq_len(iter), , drop = FALSE])
  
  list(
    pi = pi,
    delta = delta,
    eta = log(delta),
    mu0 = mu0,
    mu1 = mu0 + delta,
    posterior = r,
    loglik = loglik,
    iter = iter,
    converged = converged,
    trace = trace
  )
}



fit_additive_pois_mix_nonzero <- function(y,
                                          mu0,
                                          init_rho = NULL,
                                          init_delta = NULL,
                                          max_iter = 200L,
                                          tol = 1e-8,
                                          min_delta = 1e-8,
                                          rho_bounds = c(1e-6, 1 - 1e-6),
                                          delta_upper = NULL,
                                          zero_truncated = TRUE,
                                          threshold = 0.5,
                                          verbose = FALSE) {
  # Positive-only additive mixture:
  #
  # For y_i > 0:
  #   y_i ~ (1 - rho) Pois_+(mu0_i) + rho Pois_+(mu0_i + delta)
  #
  # For y_i == 0:
  #   posterior_i = 0
  #
  # Pois_+ means zero-truncated Poisson.
  #
  # If zero_truncated = FALSE, this literally subsets to y > 0
  # and uses ordinary Poisson densities. I do not recommend that
  # as the default, but it is useful for comparison.
  
  y <- as.numeric(y)
  mu0 <- as.numeric(mu0)
  
  if (length(y) != length(mu0)) {
    stop("y and mu0 must have the same length.")
  }
  if (any(!is.finite(y)) || any(!is.finite(mu0))) {
    stop("y and mu0 must be finite.")
  }
  if (any(y < 0) || any(mu0 < 0)) {
    stop("y and mu0 must be nonnegative.")
  }
  
  mu0 <- pmax(mu0, 1e-300)
  
  n <- length(y)
  pos <- y > 0
  n_pos <- sum(pos)
  
  posterior <- numeric(n)
  assigned <- rep(FALSE, n)
  
  if (n_pos == 0L) {
    return(list(
      rho_pos = 0,
      delta = NA_real_,
      eta = NA_real_,
      posterior = posterior,
      assigned = assigned,
      mu0 = mu0,
      mu1 = rep(NA_real_, n),
      n = n,
      n_pos = n_pos,
      loglik_pos = NA_real_,
      iter = 0L,
      converged = TRUE,
      trace = data.frame()
    ))
  }
  
  yp <- y[pos]
  mu0p <- mu0[pos]
  
  clamp <- function(x, lo, hi) {
    pmin(pmax(x, lo), hi)
  }
  
  logspace2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  # Stable log(1 - exp(-x)) for x > 0.
  log1mexp <- function(x) {
    out <- numeric(length(x))
    small <- x < log(2)
    
    out[small] <- log(-expm1(-x[small]))
    out[!small] <- log1p(-exp(-x[!small]))
    
    out
  }
  
  log_component <- function(y, mu) {
    out <- stats::dpois(y, lambda = mu, log = TRUE)
    
    if (zero_truncated) {
      out <- out - log1mexp(mu)
    }
    
    out
  }
  
  if (is.null(init_rho)) {
    # Crude initialization: among positive counts, mark cells whose
    # count is surprising under the zero-truncated null.
    q0 <- -expm1(-mu0p)
    p_tail <- stats::ppois(yp - 1, lambda = mu0p, lower.tail = FALSE) /
      pmax(q0, 1e-300)
    
    p_tail <- clamp(p_tail, 0, 1)
    
    init_rho <- mean(p_tail < 0.01)
    
    if (!is.finite(init_rho) || init_rho <= 0) {
      init_rho <- min(max(mean(yp >= 2) / 2, 0.01), 0.25)
    }
  }
  
  rho <- clamp(init_rho, rho_bounds[1], rho_bounds[2])
  
  if (is.null(init_delta)) {
    init_delta <- stats::quantile(yp, probs = 0.9,
                                  names = FALSE, type = 1) -
      stats::median(mu0p)
    
    if (!is.finite(init_delta)) {
      init_delta <- 1
    }
    
    init_delta <- max(init_delta, 0.1)
  }
  
  delta <- max(as.numeric(init_delta), min_delta)
  
  if (is.null(delta_upper)) {
    delta_upper <- max(
      10,
      10 * max(yp),
      10 * delta,
      2 * max(mu0p + yp)
    )
  }
  
  delta_upper <- max(delta_upper, delta, min_delta * 10)
  
  mstep_delta <- function(r) {
    if (sum(r) < 1e-8) {
      return(delta)
    }
    
    neg_q_eta <- function(eta) {
      d <- exp(eta)
      -sum(r * log_component(yp, mu0p + d))
    }
    
    opt <- stats::optimize(
      f = neg_q_eta,
      interval = c(log(min_delta), log(delta_upper))
    )
    
    max(exp(opt$minimum), min_delta)
  }
  
  trace <- matrix(NA_real_, nrow = max_iter, ncol = 4)
  colnames(trace) <- c("iter", "loglik_pos", "rho_pos", "delta")
  
  loglik_old <- -Inf
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    # E-step on positive-count cells only.
    log0 <- log1p(-rho) + log_component(yp, mu0p)
    log1 <- log(rho) + log_component(yp, mu0p + delta)
    
    logden <- logspace2(log0, log1)
    r <- exp(log1 - logden)
    
    loglik <- sum(logden)
    
    trace[iter, ] <- c(iter, loglik, rho, delta)
    
    if (verbose) {
      message(sprintf(
        "iter %03d  loglik_pos %.6f  rho_pos %.5g  delta %.5g  eta %.5g",
        iter, loglik, rho, delta, log(delta)
      ))
    }
    
    # M-step.
    rho_new <- clamp(mean(r), rho_bounds[1], rho_bounds[2])
    delta_new <- mstep_delta(r)
    
    if (iter > 1L) {
      ll_change <- abs(loglik - loglik_old) / (abs(loglik_old) + 1)
      par_change <- max(
        abs(rho_new - rho),
        abs(log(delta_new) - log(delta))
      )
      
      if (ll_change < tol && par_change < sqrt(tol)) {
        rho <- rho_new
        delta <- delta_new
        converged <- TRUE
        break
      }
    }
    
    rho <- rho_new
    delta <- delta_new
    loglik_old <- loglik
  }
  
  # Final E-step with final parameters.
  log0 <- log1p(-rho) + log_component(yp, mu0p)
  log1 <- log(rho) + log_component(yp, mu0p + delta)
  
  logden <- logspace2(log0, log1)
  r <- exp(log1 - logden)
  
  posterior[pos] <- r
  posterior[!pos] <- 0
  
  assigned <- posterior > threshold
  assigned[!pos] <- FALSE
  
  trace <- as.data.frame(trace[seq_len(iter), , drop = FALSE])
  
  list(
    rho_pos = rho,
    delta = delta,
    eta = log(delta),
    posterior = posterior,
    assigned = assigned,
    mu0 = mu0,
    mu1 = mu0 + delta,
    n = n,
    n_pos = n_pos,
    positive_fraction = n_pos / n,
    call_fraction = mean(assigned),
    loglik_pos = sum(logden),
    iter = iter,
    converged = converged,
    zero_truncated = zero_truncated,
    threshold = threshold,
    trace = trace
  )
}