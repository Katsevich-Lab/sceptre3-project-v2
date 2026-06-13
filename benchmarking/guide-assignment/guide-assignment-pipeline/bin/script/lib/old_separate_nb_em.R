

e_step <- function(y, curr_params, exp_offset) {
  q_j <- dnbinom(x = y, mu = exp(curr_params$gamma) * exp_offset, size = curr_params$phi1, log = TRUE) -
    dnbinom(x = y, mu = exp_offset, size = curr_params$phi0, log = TRUE) +
    qlogis(curr_params$pi)
  plogis(q_j)
}

update_phi0 <- function(y, curr_params, exp_offset, prob_is_1) {
  
  deriv_wrt_phi0 <- function(phi0) {
    mu0 = exp_offset
    dQ_dphi0 = sum((1-prob_is_1) *
                     (log(phi0) - log(phi0 + mu0) + (mu0 - y) / (mu0 + phi0) +
                        digamma(y + phi0) - digamma(phi0)))
    dQ_dphi0
  }
  fn <- function(eta0) {
    phi0 = exp(eta0)
    sum((1-prob_is_1) *
          dnbinom(x = y, mu = exp_offset,
                  size = phi0, log = TRUE))
  }
  gr = function(eta0) {
    phi0 = exp(eta0)
    deriv_wrt_phi0(phi0) * phi0
  }
  
  opt <- optim(
    par = log(curr_params$phi0),
    fn = fn, gr =gr, 
    method = "L-BFGS-B",
    lower=-15, upper=15,
    control = list(fnscale = -1)
  )
  list(phi0 = exp(opt$par))
}
# update_phi0(y, curr_params, exp_offset,prob_is_1 = e_step(y, curr_params, exp_offset))


update_gamma_and_phi1 <- function(y, curr_params, exp_offset, prob_is_1, tol=1e-6, max_iter=100) {
  # idea: we optimize over phi1, and for each value of phi1 we get gamma
  # gamma is convex so it should be really fast
  
  get_gamma <- function(phi1, gamma_init) {
    fn <- function(gamma) {
      mu1 = exp(gamma)
      sum(prob_is_1 *
            dnbinom(x = y, mu = mu1 * exp_offset,
                    size = phi1, log = TRUE))
    }
    gr = function(gamma) {
      mu1 = exp(gamma)
      sum(prob_is_1 * phi1 * (
        (y - mu1 * exp_offset) / (phi1 + mu1 * exp_offset)
      ))
    }
    opt <- optim(
      par = gamma_init,
      fn = fn, gr =gr, 
      method = "L-BFGS-B",
      lower=-20, upper=20,
      control = list(fnscale = -1)
    )
    list(gamma = opt$par)
  }
  
  get_phi1 <- function(gamma_, phi1_init) {
    mu1_with_offset =exp(gamma_) * exp_offset
    deriv_wrt_phi1 <- function(phi1) {
      dQ_dphi1 = sum(prob_is_1 *
                       (log(phi1) - log(phi1 + mu1_with_offset) + (mu1_with_offset - y) / (mu1_with_offset + phi1) +
                          digamma(y + phi1) - digamma(phi1))
      )
      dQ_dphi1
    }
    fn <- function(eta1) {
      phi1 = exp(eta1)
      sum(prob_is_1 *
            dnbinom(x = y, mu = mu1_with_offset,
                    size = phi1, log = TRUE))
    }
    gr = function(eta1) {
      phi1 = exp(eta1)
      deriv_wrt_phi1(phi1) * phi1
    }
    opt <- optim(
      par = log(phi1_init),
      fn = fn, gr =gr, 
      method = "L-BFGS-B",
      lower=-15, upper=15,
      control = list(fnscale = -1)
    )
    list(
      phi1 = exp(opt$par)
    )
  }
  
  fn_joint <- function(gamma_phi1) {
    new_gamma = get_gamma(gamma_phi1[2], gamma_phi1[1])
    new_phi1 = get_phi1(new_gamma$gamma, gamma_phi1[2])
    c(new_gamma$gamma, new_phi1$phi1)
  }
  
  curr_gamma_phi1 <- c(curr_params$gamma, curr_params$phi1)
  for(i in 1:max_iter) {
    new_gamma_phi1 <- fn_joint(curr_gamma_phi1)
    if(abs(new_gamma_phi1[1] - curr_gamma_phi1[1]) +
       abs(log(new_gamma_phi1[2]) - log(curr_gamma_phi1[2])) < tol) {
      break
    }
    curr_gamma_phi1 <- new_gamma_phi1
  }
  list(
    gamma = new_gamma_phi1[1],
    phi1 = new_gamma_phi1[2],
    n_iter_update_gamma_and_phi1 = i,
    conv = i < max_iter
  )
}
q_step <- function(y, curr_params, exp_offset, prob_is_1, tol_gamma_phi1 = 1e-6, max_iter_gamma_phi1 = 100) {
  phi0 <- update_phi0(y, curr_params, exp_offset,prob_is_1 = prob_is_1)$phi0
  new_gamma_phi1 <- update_gamma_and_phi1(y=y, curr_params, exp_offset, prob_is_1, tol = tol_gamma_phi1, max_iter = max_iter_gamma_phi1)
  new_pi <- min(max(mean(prob_is_1), 1e-8), 1 - 1e-8)  # clipping to a box for safety
  list(
    gamma = new_gamma_phi1$gamma,
    phi1 = new_gamma_phi1$phi1,
    phi0 = phi0,
    pi = new_pi
  )
}


obs_loglik_nb <- function(y, curr_params, exp_offset) {
  
  logf0 <- dnbinom(
    y, mu = exp_offset, size = curr_params$phi0, log = TRUE
  )
  
  logf1 <- dnbinom(
    x = y,
    mu = exp(curr_params$gamma) * exp_offset,
    size = curr_params$phi1,
    log = TRUE
  )
  
  a <- log1p(-curr_params$pi) + logf0
  b <- log(curr_params$pi) + logf1
  m <- pmax(a, b)
  
  sum(m + log(exp(a - m) + exp(b - m)))
}
# has_converged_nb <- function(y, exp_offset, curr_params, new_params, tol = 1e-8) {
#   ll_curr <- obs_loglik_nb(y, curr_params, exp_offset)
#   ll_new  <- obs_loglik_nb(y, new_params, exp_offset)
# 
#   abs(ll_new - ll_curr) < tol * (1 + abs(ll_curr))
# }
fit_nb_mixture <- function(y, init_params, exp_offset, tol=1e-8, max_iter=200, tol_gamma_phi1=1e-6, max_iter_gamma_phi1=100, details=0) {
  
  curr_params = init_params
  curr_params$pi <-  min(max(curr_params$pi, 1e-8), 1 - 1e-8) # this gets done to updates too, so it'll always be the case
  for (i in seq_len(max_iter)) {
    prob_is_1 <- e_step(y, curr_params, exp_offset)
    new_params <- new_params <- q_step(
      y = y, curr_params = curr_params,
      exp_offset = exp_offset, prob_is_1 = prob_is_1,
      tol_gamma_phi1 = tol_gamma_phi1, max_iter_gamma_phi1 = max_iter_gamma_phi1
    )
    ll_curr <- obs_loglik_nb(y, curr_params, exp_offset)
    ll_new  <- obs_loglik_nb(y, new_params, exp_offset)
    
    has_converged <- abs(ll_new - ll_curr) < tol * (1 + abs(ll_curr))
    
    if (details >= 1) {
      cat(sprintf(
        "iter=%d  loglik=%.8f  pi=%.6f  gamma=%.6f  phi0=%.6f  phi1=%.6f\n",
        i, ll_new,
        new_params$pi,
        new_params$gamma,
        new_params$phi0,
        new_params$phi1
      ))
    }
    
    if (has_converged) {
      curr_params <- new_params
      break
    }
    curr_params <- new_params
  }
  curr_params
}

