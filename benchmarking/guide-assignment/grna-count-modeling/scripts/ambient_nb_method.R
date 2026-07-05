#!/usr/bin/env Rscript
# ============================================================================
# Candidate flagship method: a DECONTAMINATED, OVERDISPERSION-ADAPTIVE ambient test.
#
# Idea: assign (g,c) when its count exceeds the rank-1 ambient null lambda_{g,c}=gamma_g*kappa_c,
# but (1) estimate that null from the NOISE (mask assigned entries and re-estimate the margins =
# fishash's Simpson/mean fix, done transparently), and (2) test with a NEGATIVE-BINOMIAL tail whose
# single overdispersion phi is estimated from the ambient residuals by method-of-moments (CLEANSER's
# variance fix, but cross-guide-shared and closed-form). phi -> 0 recovers the plain Poisson ambient
# test, so the method auto-adapts: Poisson on equidispersed (CROP-seq) noise, NB on overdispersed
# (direct-capture) noise -- no manual --cs/--dc switch. One knob (FDR q); phi auto-estimated.
#
# This is the only design that addresses BOTH failure modes (Simpson mean + overdispersion variance)
# as a fast closed-form test. ambient_nb_assign(N, q, n_iter, overdispersion).
# ============================================================================
suppressPackageStartupMessages(library(Matrix))

ambient_nb_assign <- function(N, q = 0.05, n_iter = 3, overdispersion = TRUE, phi_method = "robust") {
  N <- as(N, "TsparseMatrix")
  i <- N@i + 1L; j <- N@j + 1L; y <- as.numeric(N@x)
  G <- nrow(N); C <- ncol(N); nz <- length(y)
  if (nz == 0) return(list(assigned = sparseMatrix(i=integer(0), j=integer(0), dims=dim(N)), phi=0))
  assigned <- rep(FALSE, nz); phi <- 0
  for (it in seq_len(n_iter)) {
    un <- !assigned
    # rank-1 ambient mean from the UNASSIGNED (noise) margins  -> Simpson/mean decontamination
    rg <- as.numeric(tapply(y[un], factor(i[un], levels = seq_len(G)), sum)); rg[is.na(rg)] <- 0
    cc <- as.numeric(tapply(y[un], factor(j[un], levels = seq_len(C)), sum)); cc[is.na(cc)] <- 0
    tot <- sum(y[un]); if (tot <= 0) break
    lambda <- rg[i] * cc[j] / tot                      # lambda at each nonzero entry
    # overdispersion phi: ROBUST quasi-likelihood dispersion = TRIMMED mean of per-entry
    # NB dispersion estimates d = ((N-lambda)^2 - lambda)/lambda^2 over unassigned entries with
    # estimable mean (lambda >= 0.5), trimming the top 5% (signal contamination). phi->0 => Poisson.
    if (overdispersion && phi_method == "robust") {       # trimmed per-entry dispersion
      sel <- un & (lambda >= 0.5)
      if (sum(sel) >= 20) {
        d <- ((y[sel] - lambda[sel])^2 - lambda[sel]) / lambda[sel]^2
        cap <- as.numeric(quantile(d, 0.95, na.rm = TRUE))
        phi <- max(mean(d[d <= cap], na.rm = TRUE), 0)
      } else phi <- 0
    } else if (overdispersion && phi_method == "global") { # global moment estimator
      lam2_all <- (sum(rg^2) * sum(cc^2)) / tot^2
      sumLam2_un <- lam2_all - sum(lambda[assigned]^2)
      sumLam_un  <- tot - sum(lambda[assigned])
      ss <- sum(y[un]^2) - 2 * sum(y[un] * lambda[un]) + sumLam2_un
      phi <- if (sumLam2_un > 0) max((ss - sumLam_un) / sumLam2_un, 0) else 0
    }
    # one-sided upper-tail test  P(count >= y)
    p <- if (phi > 0) stats::pnbinom(y - 1, mu = lambda, size = 1/phi, lower.tail = FALSE)
         else          stats::ppois(y - 1, lambda, lower.tail = FALSE)
    padj <- stats::p.adjust(p, "BH")
    new_assigned <- !is.na(padj) & padj < q
    converged <- it > 1 && all(new_assigned == assigned)
    assigned <- new_assigned
    if (converged) break
  }
  list(assigned = sparseMatrix(i = i[assigned], j = j[assigned], x = TRUE,
                               dims = dim(N), dimnames = dimnames(N)), phi = phi)
}
