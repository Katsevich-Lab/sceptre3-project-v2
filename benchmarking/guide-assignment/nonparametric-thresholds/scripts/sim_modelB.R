#!/usr/bin/env Rscript
# =============================================================================
# MODEL B  --  owned mechanistic (CellBender-style) generative model.
#
# A transparent, lab-owned DGP with exactly the load-bearing realistic features,
# calibrated to the survey/barnyard and DELIBERATELY misspecified vs every
# method (no method models all of: decoupled depths, a shared rank-1 ambient
# pool, MOI co-occurrence, doublets).  It reaches regimes the barnyard can't
# (controlled MOI/signal, known FDR) and is NOT fishash's package, so a
# fishash-family win here is not home-field.
#
# Per (g,c):
#   infection pool   p_g     ~ Dirichlet(alpha_inf)             (guide popularity)
#   ambient pool     chi_g   ~ Dirichlet(alpha_amb)             (calibrated Gini)
#   guide capture    w_g     ~ lognormal(0, w_sigma)
#   cell signal depth d_c    ~ lognormal(0, d_sigma)            (capture)
#   droplet amb depth kappa_c~ lognormal(log kappa_bar, k_sigma) -- DECOUPLED from d_c
#   integrations     n_c     ~ Poisson(MOI); guides ~ multinomial(p_g)  => co-occurrence
#   signal counts            ~ NB(mu = mu_pert * d_c * w_g, size = theta_sig)
#   ambient counts           ~ Poisson(chi_g * kappa_c)   [NB(1/phi_noise) if phi_noise>0]
#   doublets: a fraction delta of cells are summed with a fresh partner cell;
#             the OBSERVED integration truth is the union, doublet is labeled.
#
# The two-setting axis mirrors barnyard purity: delta>0 = "doublets in" (~0.90),
# delta=0 = "clean" (~0.99).  NB: with full truth, a doublet's guides are
# correctly labeled signal (union), so -- unlike the barnyard's species-only
# truth -- doublets do NOT corrupt the ambient here; ambient stays Poisson,
# matching real *post-doublet* ambient.  That semantic difference is the point.
# Outputs: results/sim_framework/datasets/B__*  +  manifest_modelB.csv
# =============================================================================
suppressPackageStartupMessages({ library(Matrix) })
source(file.path(getwd(), "scripts", "sim_lib.R"))

# generate `C` singlet cells; returns signal (sparse), ambient (sparse), Z (lgl)
# Perturbation: if pert_rate is given, each guide g perturbs each cell independently
# at its own rate (Bernoulli) -> guide g hits ~pert_rate_g * C cells, matching the
# real per-guide signal fraction (co-occurrence emerges from overlap).  Otherwise
# the MOI/multinomial model (guides per cell ~ Poisson(moi)).
.gen_cells <- function(C, G, p_g, chi_g, w_g, mu_pert, theta_sig,
                       d_sigma, kappa_bar, k_sigma, moi, phi_noise, pert_rate = NULL) {
  # Lognormal capture/depth with the -sigma^2/2 mean correction so each factor's
  # expectation is exactly 1 (or kappa_bar).  Without the correction E inflates
  # by exp(sigma^2/2) -- negligible at sigma=0.05, ~8% at 0.4.
  d_c     <- exp(rnorm(C, -0.5 * d_sigma^2, d_sigma))                  # E[d_c]    = 1
  kappa_c <- kappa_bar * exp(rnorm(C, -0.5 * k_sigma^2, k_sigma))      # E[kappa_c]= kappa_bar
  # vector-or-scalar guards: silent length-mismatch would corrupt indexing
  stopifnot(length(mu_pert)   %in% c(1L, G),
            length(theta_sig) %in% c(1L, G),
            is.null(pert_rate) || length(pert_rate) %in% c(1L, G))
  if (!is.null(pert_rate)) {
    pr <- if (length(pert_rate) == 1) rep(pert_rate, G) else pert_rate
    n_pert <- rbinom(G, C, pmin(pmax(pr, 0), 1))                      # perturbed cells per guide
    zi <- rep.int(seq_len(G), n_pert)
    zj <- unlist(lapply(seq_len(G), function(g) if (n_pert[g] > 0) sample.int(C, n_pert[g]) else integer(0)), use.names = FALSE)
  } else {
    n_c <- rpois(C, moi)
    zi <- integer(0); zj <- integer(0)
    for (cc in which(n_c > 0)) {
      k <- min(n_c[cc], G); g_hit <- sample.int(G, k, replace = FALSE, prob = p_g)
      zi <- c(zi, g_hit); zj <- c(zj, rep.int(cc, k))
    }
  }
  mu_g <- if (length(mu_pert) == 1)  rep(mu_pert, G)  else mu_pert    # per-guide signal level (vector ok)
  th_g <- if (length(theta_sig) == 1) rep(theta_sig, G) else theta_sig # per-guide NB dispersion (vector ok)
  mu_sig <- mu_g[zi] * d_c[zj] * w_g[zi]
  sig <- rnbinom(length(zi), mu = pmax(mu_sig, 1e-6), size = pmax(th_g[zi], 0.02))
  S <- sparseMatrix(i = zi, j = zj, x = as.numeric(sig), dims = c(G, C))
  Z <- sparseMatrix(i = zi, j = zj, x = TRUE, dims = c(G, C))
  lam <- outer(chi_g, kappa_c)                                        # rank-1 ambient field
  amb_x <- if (phi_noise > 0) rnbinom(length(lam), mu = as.numeric(lam), size = 1 / phi_noise)
           else rpois(length(lam), as.numeric(lam))
  B <- drop0(Matrix(amb_x, nrow = G, ncol = C, sparse = TRUE))
  list(S = S, B = B, Z = Z)
}

build_modelB <- function(G = 100, C = 3500, moi = 1, mu_pert = 10, theta_sig = 1,
                         alpha_inf = 3, alpha_amb = 3, w_sigma = 0.3, d_sigma = 0.4,
                         kappa_bar = 2.5, k_sigma = 0.5, doublet_rate = 0, phi_noise = 0,
                         pert_rate = NULL, seed = 1) {
  set.seed(seed)
  p_g   <- rdirichlet1(rep(alpha_inf, G))
  chi_g <- rdirichlet1(rep(alpha_amb, G))
  w_g   <- exp(rnorm(G, -0.5 * w_sigma^2, w_sigma))     # E[w_g] = 1 (lognormal mean correction)
  base <- .gen_cells(C, G, p_g, chi_g, w_g, mu_pert, theta_sig, d_sigma, kappa_bar, k_sigma, moi, phi_noise, pert_rate)
  S <- base$S; B <- base$B; Z <- base$Z; is_dbl <- logical(C)
  D <- round(doublet_rate * C)
  if (D > 0) {
    didx <- sort(sample(C, D)); is_dbl[didx] <- TRUE
    part <- .gen_cells(D, G, p_g, chi_g, w_g, mu_pert, theta_sig, d_sigma, kappa_bar, k_sigma, moi, phi_noise, pert_rate)
    add_at <- function(M, Madd, cols) { M <- as(M, "TsparseMatrix"); Ma <- as(Madd, "TsparseMatrix")
      sparseMatrix(i = c(M@i + 1L, Ma@i + 1L), j = c(M@j + 1L, cols[Ma@j + 1L]),
                   x = c(M@x, Ma@x), dims = dim(M), use.last.ij = FALSE) }
    S <- add_at(S, part$S, didx); B <- add_at(B, part$B, didx)
    Zd <- as(Z, "TsparseMatrix"); Pd <- as(part$Z, "TsparseMatrix")
    Z <- sparseMatrix(i = c(Zd@i + 1L, Pd@i + 1L), j = c(Zd@j + 1L, didx[Pd@j + 1L]),
                      x = TRUE, dims = c(G, C))
  }
  Y <- drop0(S + B)
  rn <- paste0("g", seq_len(G)); cn <- paste0("c", seq_len(C))
  dimnames(Y) <- dimnames(Z) <- list(rn, cn)
  list(counts = as(Y, "CsparseMatrix"), Z = as(Z, "lgCMatrix"),
       is_doublet = is_dbl, ambient_nnz = sum(B > 0), signal_nnz = sum(S > 0))
}

# Sourcing this file always defines build_modelB() / .gen_cells() only.
# Model B regimes are built by sim_regimes.R (one regime per real dataset).
# The MODELB_FUNCS_ONLY flag is kept for backward compatibility with any caller
# that still sets it before source()ing.
