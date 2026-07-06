#!/usr/bin/env Rscript
# ============================================================================
# Canonical noise-conditioned contingency test for guide assignment, modular so the
# two improvements can be ablated cleanly against fishash:
#
#   cell_margin = "observed" : draws = raw library N_{:,c}            (== fishash)
#               = "ambient"  : draws = noise_col[c] - noise[g,c] + y  (decontaminated DEPTH)
#   tail        = "hyper"    : hypergeometric (equidispersed, == fishash)
#               = "betabinom": beta-binomial conditional tail (overdispersion rho from noise)
#
# Noise field = fishash::impute_masked_counts (rank-1 Poisson completion of the assigned-masked
# matrix), iterated with fishash's mask schedule. FDR by BH (fdr="BH") or Guo-Sarkar block (fdr="GS").
# With cell_margin="observed", tail="hyper", fdr="GS" this reproduces fishash.
# ============================================================================
suppressPackageStartupMessages({library(fishash); library(Matrix); library(extraDistr)})

contingency_assign <- function(counts, q = 0.05, refit = 10, min_count = 2,
                               cell_margin = c("ambient","observed"),
                               tail = c("hyper","nb"), fdr = c("GS","BH"),
                               rho_fixed = NULL, init_assigned = NULL) {
  # init_assigned: optional logical/0-1 sparse matrix (same dims as counts) seeding the assigned
  # set A for the FIRST background estimate. Default NULL = "no integrations" (background_1 = raw
  # counts, so the initial depth d_c is the raw library size). Passing e.g. (counts > 2) seeds the
  # first depth with the CLEANSER-style <=2 ambient estimate. Only affects iteration 1's rate field;
  # the schedule then takes over. Used to probe initialization sensitivity of the fixed point.
  # rho_fixed: if supplied, use this GLOBAL negative-binomial dispersion for tail="nb"
  # instead of the per-entry Pearson-residual estimate (which is signal-contaminated).
  # The clean way to set it is from the aggregate ambient count-2 excess over Poisson
  # (soup-dominated): obs/exp_Poisson(count=2) ~ 1 + rho, so rho ~ excess - 1.
  cell_margin <- match.arg(cell_margin); tail <- match.arg(tail); fdr <- match.arg(fdr)
  counts  <- as(counts, "CsparseMatrix")
  G <- nrow(counts); C <- ncol(counts)
  obs_col <- Matrix::colSums(counts)
  n_rows <- sum(Matrix::rowSums(counts) > 0); n_cols <- sum(obs_col > 0)
  n_entries <- as.numeric(n_rows) * as.numeric(n_cols)
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  if (!is.null(init_assigned)) init_assigned <- as(init_assigned, "CsparseMatrix")
  mask <- init_assigned; prev <- NULL; assigned <- NULL; rho <- 0
  for (it in seq_len(refit + 1)) {
    background <- if (it == 1 && is.null(init_assigned)) counts
                  else fishash::impute_masked_counts(counts, mask)
    nr <- Matrix::rowSums(background); nc <- Matrix::colSums(background); Tn <- sum(background)
    bgc <- background[cbind(i, j)]                      # noise estimate at each nonzero (g,c)
    K <- nr[i] - bgc + y                                # row (guide) margin: decontaminated (fishash)
    if (cell_margin == "ambient") { draws <- nc[j] - bgc + y; pop <- Tn - bgc + y }
    else                          { draws <- obs_col[j];      pop <- Tn - nc[j] + obs_col[j] }
    if (tail == "hyper") {
      logp <- stats::phyper(y - 1, K, pop - K, draws, lower.tail = FALSE, log.p = TRUE)
    } else {                                            # overdispersed conditional tail (NB; the
      # numerically-stable Poisson limit of the beta-binomial, since pop >> draws and p is small).
      # Conditional mean lambda = draws*K/pop; dispersion phi from unassigned Pearson residuals.
      lambda <- draws * (K / pop)
      un <- if (is.null(assigned)) rep(TRUE, length(y)) else !assigned[cbind(i, j)]
      # estimate phi from CONFIDENTLY-NOISE entries only: unassigned, estimable mean, and count not
      # in the upper tail (excludes signal that survived masking, esp. at high MOI). median-based.
      if (!is.null(rho_fixed)) {
        rho <- rho_fixed                                # clean global dispersion (count-2 calibrated)
      } else {
        ycap <- if (any(un)) as.numeric(quantile(y[un], 0.98, na.rm = TRUE)) else Inf
        sel <- un & (lambda >= 0.5) & (y <= ycap)      # drop the upper-tail (signal leakage)
        if (sum(sel) >= 20) {
          d <- ((y[sel] - lambda[sel])^2 - lambda[sel]) / lambda[sel]^2
          d <- d[d <= quantile(d, 0.95, na.rm = TRUE)] # trim residual outliers, keep OD signal
          rho <- max(mean(d, na.rm = TRUE), 0)
        } else rho <- 0
      }
      logp <- if (rho > 0) stats::pnbinom(y - 1, mu = lambda, size = 1/rho, lower.tail = FALSE, log.p = TRUE)
              else          stats::phyper(y - 1, K, pop - K, draws, lower.tail = FALSE, log.p = TRUE)
    }
    if (fdr == "BH") {
      A <- stats::p.adjust(exp(logp), "BH") < q
    } else {                                            # Guo-Sarkar block (per-cell), like fishash
      ml <- sparseMatrix(i = i, j = j, x = logp, dims = dim(counts))
      colmin <- sparseMatrixStats::colMins(ml)
      B <- sum(stats::p.adjust(pmin(exp(colmin) * n_rows, 1), "BH") <= q)
      cutoff <- log(q) - log(n_entries) + log(max(B, 1))
      A <- logp <= cutoff
    }
    A <- A & !is.na(A) & (y >= min_count)
    assigned <- sparseMatrix(i = i[A], j = j[A], x = TRUE, dims = dim(counts), dimnames = dimnames(counts))
    mask <- if (it > 3 && !is.null(mask)) (mask | assigned) else assigned
    if (it > 1 && sum(abs(prev - assigned)) == 0) break
    prev <- assigned
  }
  list(assigned = assigned, rho = rho)
}
