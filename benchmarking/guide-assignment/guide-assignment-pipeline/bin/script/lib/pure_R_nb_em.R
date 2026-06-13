
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# logaddexp <- function(x, y) {
#   m <- pmax(x, y)
#   out <- m + log(exp(x - m) + exp(y - m))
#   out[is.infinite(m) & m < 0] <- -Inf
#   out
# }
logaddexp <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log1p(exp(-abs(x - y)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}
obs_ll_vec <- function(gamma, pi, y, offset, phi) {
  t0 <- offset - log(phi)
  t1 <- gamma + offset - log(phi)
  
  a0 <- log1p(-pi) - phi * softplus(t0) - y * softplus(-t0)
  a1 <- log(pi)    - phi * softplus(t1) - y * softplus(-t1)
  
  logaddexp(a0, a1)
}

e_step <- function(curr_params, y, offset, phi) {
  t0 <- offset - log(phi)
  t1 <- curr_params$gamma + offset - log(phi)
  
  q <- -phi * softplus(t1) - y * softplus(-t1) + phi * softplus(t0) + y * softplus(-t0) + qlogis(curr_params$pi)
  prob_is_1 <- plogis(q)
  prob_is_1
}


gamma_gradient <- function(gamma_, y_plus_phi, offset_minus_log_phi, phi, prob_is_1) {
  sum(prob_is_1 * (y_plus_phi * plogis(-gamma_ - offset_minus_log_phi) - phi))
}

q_step <- function(curr_params, prob_is_1, y, offset, phi) {
  # 1. optimize over gamma
  # 2. update pi
  
  # pre-compute args just for `gamma_gradient()`
  y_plus_phi = y + phi
  offset_minus_log_phi = offset - log(phi)
  
  update_gamma <- uniroot(
    gamma_gradient,
    interval = c(curr_params$gamma - 2, curr_params$gamma + 2),
    extendInt = "downX",
    check.conv = TRUE,
    tol = 1e-8,
    y_plus_phi=y_plus_phi, offset_minus_log_phi=offset_minus_log_phi,
    phi=phi, prob_is_1=prob_is_1
  )
  
  
  gamma_new <- update_gamma$root
  pi_new <- mean(prob_is_1)
  list(gamma = gamma_new, pi = pi_new)
}

# optimization ideaS:
# - probably a lot of stuff involving y and phi that i could precompute
fit_nb_mix_given_overdisp <- function(y, init_params, offset, phi, tol=1e-6, max_iter = 200, verbose=FALSE) {
  
  stopifnot(setequal(names(init_params), c("gamma", "pi")))
  curr_params = init_params
  curr_avg_obs_loglik <- mean(obs_ll_vec(gamma = curr_params$gamma, pi = curr_params$pi, y = y, offset = offset, phi = phi))
  
  for(i in seq_len(max_iter)) {
    
    prob_is_1 <- e_step(curr_params = curr_params, y = y, offset = offset, phi = phi)
    new_params <- q_step(curr_params, prob_is_1 = prob_is_1, y = y, offset = offset, phi = phi)
    
    # check change in avg observed log likelihood
    # doing this instead of a relative change, since const(y, phi) disappear here
    new_avg_obs_loglik <- mean(obs_ll_vec(gamma = new_params$gamma, pi = new_params$pi, y = y, offset = offset, phi = phi))
    
    if(verbose) {
      cat(sprintf(
        "iter=%d  avg loglik=%.8f  pi=%.6f  gamma=%.6f \n",
        i, new_avg_obs_loglik,
        new_params$pi,
        new_params$gamma
      ))
    }
    
    if(abs(new_avg_obs_loglik - curr_avg_obs_loglik) < tol) {
      break
    }
    curr_params = new_params
    curr_avg_obs_loglik = new_avg_obs_loglik
  }
  list(
    gamma = curr_params$gamma, pi = curr_params$pi,
    n_iter = i, converged = i < max_iter,
    prob_is_1 = prob_is_1
  )
} 
