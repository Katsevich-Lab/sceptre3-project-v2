# =============================================================================
# Simple nonparametric per-gRNA thresholding schemes for guide assignment.
#
# Both methods operate on the raw integer count vector of a single gRNA across
# cells (no covariate adjustment, no mixture model -- deliberately simple). They
# return an integer threshold `t`; a cell is assigned (perturbed) iff
#   count >= t.
# When no usable threshold is found (gRNA is unimodal / has no perturbed mode)
# the method abstains: t = Inf, so no cell is assigned.
#
#   otsu_threshold_log1p()       Otsu's between-class-variance split on log(1+y).
#   smoothed_valley_threshold()  smoothed discrete PMF on log(y+1); cut at the
#                                empirical valley (min PMF) between the two
#                                dominant modes:  t = argmin_{m1 < k < m2} p_hat(k).
#
# Otsu finds the best variance-separating split; the smoothed-valley method
# finds the empirical dip. They coincide for cleanly separated bimodal data and
# diverge under skew / overlap.
# =============================================================================

# ---- Otsu threshold on log(1+y) --------------------------------------------
# Maximize the between-class variance of z = log(1 + count) over all integer
# split levels. Returns the smallest count assigned to the "perturbed" (upper)
# class, so assignment is `count >= t`.
otsu_threshold_log1p <- function(counts) {
  counts <- as.numeric(counts)
  n <- length(counts)
  if (n == 0L) return(list(t = Inf, ok = FALSE, reason = "empty"))

  tab      <- table(counts)               # frequency of each distinct count
  levels_k <- as.numeric(names(tab))      # distinct counts, ascending
  w        <- as.numeric(tab)             # frequency per level
  K        <- length(levels_k)
  if (K < 2L) return(list(t = Inf, ok = FALSE, reason = "single_level"))

  z <- log1p(levels_k)                    # intensity on the log scale
  p <- w / sum(w)                         # normalized weights

  # Cumulative class-0 weight and weighted z-sum for each split index i
  # (class0 = levels 1..i, class1 = levels i+1..K).
  w0     <- cumsum(p)
  mu_cum <- cumsum(p * z)
  muT    <- mu_cum[K]

  i   <- seq_len(K - 1L)                   # candidate splits (exclude last)
  w0i <- w0[i]
  w1i <- 1 - w0i
  mu0 <- mu_cum[i] / w0i
  mu1 <- (muT - mu_cum[i]) / w1i
  sigma_b <- w0i * w1i * (mu0 - mu1)^2     # between-class variance

  best <- which.max(sigma_b)
  t    <- levels_k[best + 1L]              # first count in the upper class
  list(t = t, ok = TRUE, sigma_b = sigma_b[best],
       n_levels = K, reason = "otsu")
}

# ---- Smoothed-valley threshold ---------------------------------------------
# Build a smoothed discrete PMF p_hat(k) over the integer support on the
# z = log(1 + k) scale (Gaussian KDE of log1p(counts), interpolated back to the
# integer grid). Find the two dominant modes (highest local maxima) and cut at
# the integer of minimum smoothed PMF strictly between them.
#
#   bw            kernel bandwidth in z (log1p) units; NULL => nrd0 of z.
#   min_dip_ratio require valley height < min_dip_ratio * (smaller mode height)
#                 for the distribution to count as genuinely bimodal; else
#                 abstain. 1 = accept any dip; lower = stricter.
smoothed_valley_threshold <- function(counts, bw = NULL, min_dip_ratio = 0.95) {
  counts <- as.numeric(counts)
  n      <- length(counts)
  kmax   <- if (n) max(counts) else 0
  if (kmax < 1) return(list(t = Inf, ok = FALSE, reason = "all_zero"))

  z <- log1p(counts)
  if (is.null(bw)) bw <- stats::bw.nrd0(z)
  if (!is.finite(bw) || bw <= 0) bw <- 0.5

  # Smooth density of z via FFT KDE, then evaluate on the integer grid.
  dens <- stats::density(z, bw = bw, n = 512,
                         from = 0, to = log1p(kmax) + 3 * bw)
  k     <- 0:kmax
  zk    <- log1p(k)
  phat  <- stats::approx(dens$x, dens$y, xout = zk, rule = 2)$y
  phat[phat < 0] <- 0
  if (sum(phat) > 0) phat <- phat / sum(phat)

  # Local maxima over the integer grid (strict interior peaks + endpoint check).
  is_peak <- logical(length(k))
  if (length(k) >= 2L) {
    is_peak[1]               <- phat[1] > phat[2]
    is_peak[length(k)]       <- phat[length(k)] > phat[length(k) - 1L]
    if (length(k) >= 3L) {
      mid <- 2:(length(k) - 1L)
      is_peak[mid] <- phat[mid] >= phat[mid - 1L] & phat[mid] > phat[mid + 1L]
    }
  }
  peak_idx <- which(is_peak)
  if (length(peak_idx) < 2L)
    return(list(t = Inf, ok = FALSE, reason = "unimodal",
                n_modes = length(peak_idx), phat = phat, k = k))

  # Two dominant modes = the two tallest peaks, ordered in k.
  top2 <- peak_idx[order(phat[peak_idx], decreasing = TRUE)[1:2]]
  m1   <- min(top2)
  m2   <- max(top2)

  # Valley = integer of minimum smoothed PMF strictly between the modes.
  between <- (m1 + 1L):(m2 - 1L)          # indices into k (1-based, k starts 0)
  if (length(between) < 1L)
    return(list(t = Inf, ok = FALSE, reason = "adjacent_modes",
                phat = phat, k = k))
  v_rel    <- between[which.min(phat[between])]
  valley_k <- k[v_rel]
  valley_h <- phat[v_rel]

  # Bimodality / dip check: the trough must be a real dip below both modes.
  mode_h <- min(phat[m1], phat[m2])
  if (valley_h >= min_dip_ratio * mode_h)
    return(list(t = Inf, ok = FALSE, reason = "shallow_dip",
                valley_k = valley_k, mode1 = k[m1], mode2 = k[m2],
                phat = phat, k = k))

  list(t = valley_k, ok = TRUE, reason = "valley",
       mode1 = k[m1], mode2 = k[m2], valley_h = valley_h, mode_h = mode_h,
       bw = bw, phat = phat, k = k)
}

# ---- Ambient-proportion test (principled per-(cell,guide) assignment) -------
# The whole method in one sentence: assign guide g to cell c when its UMI count
# significantly exceeds what ambient contamination alone would give, where the
# expected ambient count is (cell's total guide UMIs) x (guide's share of all
# guide UMIs); keep the calls at a controlled false-discovery rate q.
#
#   pi_g       = sum_c y[g,c] / sum_{g,c} y[g,c]     (guide g's ambient share)
#   lambda_gc  = n_c * pi_g                          (expected ambient count)
#   p_gc       = P(Pois(lambda_gc) >= y_gc)          (upper-tail vs ambient null)
#   assign if BH-adjusted p_gc < q.
#
# Unlike a single per-guide count threshold, the cutoff adapts per cell (via the
# library size n_c) and per guide (via pi_g); the one knob q is the ambient FDR
# and is interpretable + label-free. This is the contingency-table / hypergeo-
# metric family (geomux, Teyssier 2024; fishash, Kamm et al. 2026). The exact
# hypergeometric is the default (robust); model="poisson" is the simpler
# intuition/approximation -- identical on real many-guide data, but it can
# under-call when there are very few guides (pi_g self-inflated by its own
# signal). n_iter>1 re-estimates pi_g from the uncalled background (off by
# default; it doesn't help on real data).
ambient_test_assign <- function(grna_matrix, q = 0.05,
                                model = c("hypergeometric", "poisson"), n_iter = 1) {
  model <- match.arg(model)
  tg  <- methods::as(grna_matrix, "TsparseMatrix")
  i   <- tg@i + 1L; j <- tg@j + 1L; y <- tg@x        # nonzero (guide, cell, count)
  n_c <- Matrix::colSums(grna_matrix)                # per-cell library size (full)
  N   <- sum(y)                                      # full total (population)
  G   <- nrow(grna_matrix)
  called <- rep(FALSE, length(y)); pval <- rep(1, length(y))
  for (it in seq_len(n_iter)) {
    bg <- !called                                    # background = uncalled pairs
    K_bg <- as.numeric(tapply(y[bg], factor(i[bg], levels = seq_len(G)), sum)); K_bg[is.na(K_bg)] <- 0
    pi_g <- K_bg / sum(y[bg])                        # ambient composition
    if (model == "poisson") {
      pval <- stats::ppois(y - 1, n_c[j] * pi_g[i], lower.tail = FALSE)
    } else {
      K_eff <- round(pi_g * N)                       # background profile rescaled to full N
      pval  <- stats::phyper(y - 1, K_eff[i], N - K_eff[i], n_c[j], lower.tail = FALSE)
    }
    called <- stats::p.adjust(pval, "BH") < q
    called[is.na(called)] <- FALSE
  }
  A <- Matrix::sparseMatrix(i = i[called], j = j[called], x = TRUE,
                            dims = dim(grna_matrix), dimnames = dimnames(grna_matrix))
  list(assignment_matrix = A, pval = pval, called = called, i = i, j = j, y = y)
}

# ---- Apply a per-gRNA threshold function across a count matrix --------------
# grna_matrix: gRNAs (rows) x cells (cols), sparse or dense.
# threshold_fn: function(counts_vector) -> list with $t (integer threshold).
# Returns a sparse logical assignment matrix (count >= t), plus a per-gRNA
# data frame of thresholds / diagnostics.
assign_by_threshold <- function(grna_matrix, threshold_fn,
                                min_nonzero = 0L, verbose = FALSE) {
  n_grna <- nrow(grna_matrix)
  n_cell <- ncol(grna_matrix)
  # Row-oriented access is fast on a dgRMatrix.
  gm <- methods::as(grna_matrix, "RsparseMatrix")

  i_idx <- integer(0)
  j_idx <- integer(0)
  diag_rows <- vector("list", n_grna)

  for (g in seq_len(n_grna)) {
    counts <- gm[g, ]
    counts <- as.numeric(counts)
    n_nz   <- sum(counts > 0)
    if (n_nz < min_nonzero) {
      res <- list(t = Inf, ok = FALSE, reason = "too_few_nonzero")
    } else {
      res <- threshold_fn(counts)
    }
    assigned <- if (is.finite(res$t)) which(counts >= res$t) else integer(0)
    if (length(assigned)) {
      i_idx <- c(i_idx, rep.int(g, length(assigned)))
      j_idx <- c(j_idx, assigned)
    }
    diag_rows[[g]] <- data.frame(
      grna_id    = rownames(grna_matrix)[g] %||% as.character(g),
      t          = res$t,
      ok         = isTRUE(res$ok),
      reason     = res$reason %||% NA_character_,
      n_nonzero  = n_nz,
      n_assigned = length(assigned),
      stringsAsFactors = FALSE
    )
    if (verbose && g %% 200 == 0) cat("  ", g, "/", n_grna, "gRNAs\n")
  }

  assignment_matrix <- Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = TRUE,
    dims = c(n_grna, n_cell),
    dimnames = list(rownames(grna_matrix), colnames(grna_matrix))
  )
  list(assignment_matrix = assignment_matrix,
       diagnostics = do.call(rbind, diag_rows))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
