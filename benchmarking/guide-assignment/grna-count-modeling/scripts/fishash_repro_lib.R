#!/usr/bin/env Rscript
# ============================================================================
# Fishash+ (Poisson-tail) run through the fishash authors' OWN simulations.
#
# "Fishash+ Poisson" = the depth-decontaminated ambient test, but scored with a
# POISSON upper tail against the denoised rank-1 ambient rate  lambda_{g,c} = a_g * d_c
# (canonical_model.qmd, unconditional form), INSTEAD of fishash's hypergeometric
# contingency table. The rank-1 ambient factors a_g (guide share) and d_c (cell
# ambient DEPTH) come from the SAME masked rank-1 Poisson completion fishash uses
# (fishash::impute_masked_counts), iterated with fishash's exact mask schedule and
# GS FDR. Because the Poisson rate uses the fitted ambient depth d_c (not the raw
# library N_{:,c} that fishash plugs in as hypergeometric "draws"), the depth fix
# is intrinsic -- this is the Poisson-form analogue of depth_fix.
#
# Everything else (mask schedule, min_count, exclude_empty, GS FDR at q) is copied
# byte-for-byte from fishash::fishash so the ONLY differences vs fishash are
#   (tail: hypergeometric -> Poisson)  and  (draws: raw library -> denoised depth).
#
# Scoring = the authors' own per-entry pooled confusion metric (their
# bin/process_assignments.R): Precision, Recall on the "full" and "nonzero" subsets.
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(SummarizedExperiment)
  library(sparseMatrixStats); library(fishash)
})

# ---- single stable location for the (regenerable) sim dataset cache ----------
# Datasets are byte-exact from the authors' seeds, so this is a CACHE, not a precious
# artifact. Override with $FISH_DATA; default is a gitignored dir under results/.
# (All scripts setwd() to the project root before sourcing, so the relative path resolves.)
sim_data_dir <- function(scenario = NULL) {
  base <- Sys.getenv("FISH_DATA", file.path("results", "fishash_repro", "datasets"))
  if (is.null(scenario)) base else file.path(base, scenario)
}

# ---- rank-1 Poisson factors of the masked matrix -----------------------------
# Reproduces the alternating fit inside fishash::impute_masked_counts EXACTLY, but
# returns the factors (a = guide_freqs summing to 1, d = cell_sizes = ambient depth)
# so we can form the full rank-1 ambient rate lambda_{g,c} = a_g * d_c at EVERY entry
# (impute_masked_counts only fills the *masked* entries; unmasked ones keep observed).
rank1_factors <- function(counts, mask = NULL, eps = 1e-4, max_iter = 10) {
  counts <- as(counts, "CsparseMatrix")
  G <- nrow(counts); C <- ncol(counts)
  if (is.null(mask)) {                       # no masking -> closed-form independence fit
    rs <- Matrix::rowSums(counts); cs <- Matrix::colSums(counts); tot <- sum(counts)
    a <- rs / tot; a[!is.finite(a)] <- 0
    return(list(a = as.numeric(a), d = as.numeric(cs)))   # a_g*d_c = rs*cs/tot
  }
  mask <- as(mask, "CsparseMatrix")
  counts0 <- counts - counts * mask          # zero out masked (assigned) entries
  rs0 <- Matrix::rowSums(counts0); cs0 <- Matrix::colSums(counts0)
  skip_row <- rs0 == 0; skip_col <- cs0 == 0
  cell_sizes <- rep(1, C); guide_freqs <- rep(0, G)
  prev_freqs <- NULL; prev_sizes <- NULL
  for (it in seq_len(max_iter)) {
    gf <- sum(cell_sizes) - Matrix::rowSums(mask %*% Matrix::Diagonal(x = cell_sizes))
    gf <- rs0 / gf; gf[skip_row] <- 0; gf <- gf / sum(gf)
    cs <- 1 - Matrix::colSums(Matrix::Diagonal(x = gf) %*% mask)
    cs <- cs0 / cs; cs[skip_col] <- 0
    guide_freqs <- gf; cell_sizes <- cs
    if (it > 1) {
      dg <- max(abs(log(gf) - log(prev_freqs)), na.rm = TRUE)
      ds <- max(abs(log(cs) - log(prev_sizes)), na.rm = TRUE)
      if (is.finite(dg) && is.finite(ds) && dg < eps && ds < eps) break
    }
    prev_freqs <- gf; prev_sizes <- cs
  }
  list(a = as.numeric(guide_freqs), d = as.numeric(cell_sizes))
}

# ---- Guo-Sarkar block FDR cutoff (copied from fishash_internal) ---------------
.gs_cutoff <- function(i, j, log_pval, dims, n_rows, n_entries, q) {
  mat_logpval <- sparseMatrix(i = i, j = j, x = log_pval, dims = dims, index1 = TRUE)
  colmin <- sparseMatrixStats::colMins(mat_logpval)
  B <- sum(stats::p.adjust(pmin(exp(colmin) * n_rows, 1), method = "BH") <= q)
  log(q) - log(n_entries) + log(B)
}

# ---- Fishash+ Poisson --------------------------------------------------------
# tail = "poisson": lambda = a_g*d_c, p = P(Pois(lambda) >= y)         [ambient depth d_c]
# cell_margin: "ambient" (ours; d_c = denoised depth) is the ONLY sensible Poisson
#   form -- kept as an arg only to allow the raw-depth ablation lambda = a_g*N_{:,c}.
fishash_plus_poisson <- function(counts, q = 0.05, refit = 10, min_count = 2,
                                 cell_margin = c("ambient", "observed"),
                                 exclude_empty = TRUE, return_score = FALSE) {
  cell_margin <- match.arg(cell_margin)
  counts <- as(counts, "CsparseMatrix")
  col_sums <- Matrix::colSums(counts); row_sums <- Matrix::rowSums(counts)
  n_rows <- if (exclude_empty) sum(row_sums > 0) else nrow(counts)
  n_cols <- if (exclude_empty) sum(col_sums > 0) else ncol(counts)
  n_entries <- as.numeric(n_rows) * as.numeric(n_cols)
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  mask <- NULL; prev <- NULL; assigned <- NULL; log_pval <- NULL
  for (it in seq_len(refit + 1)) {
    fac <- if (it == 1) rank1_factors(counts, NULL) else rank1_factors(counts, mask)
    d_use <- if (cell_margin == "ambient") fac$d[j] else col_sums[j]   # depth: denoised vs raw
    # floor guards d_c=0 (fully-masked cells): their residual counts get forced significant.
    # Benign in practice -- such entries are the cell's own dominant (true-positive) guides.
    lambda <- pmax(fac$a[i] * d_use, 1e-9)
    log_pval <- stats::ppois(y - 1, lambda, lower.tail = FALSE, log.p = TRUE)
    cutoff <- .gs_cutoff(i, j, log_pval, dim(counts), n_rows, n_entries, q)
    A <- (log_pval <= cutoff) & (y >= min_count)
    assigned <- sparseMatrix(i = i[A], j = j[A], x = TRUE, dims = dim(counts),
                             dimnames = dimnames(counts))
    mask <- if (it > 3 && !is.null(mask)) (mask | assigned) else assigned
    if (it > 1 && sum(abs(prev - assigned)) == 0) break
    prev <- assigned
  }
  if (!return_score) return(assigned)
  # continuous score (converged Poisson-tail -log p) for PR curves; count<min_count can't be positive
  score <- -log_pval; score[y < min_count] <- -Inf
  list(assigned = assigned, i = i, j = j, score = score)
}

# average precision (full-PR-curve AUPRC) of a per-entry score vs boolean truth
average_precision <- function(score, truth) {
  o <- order(score, decreasing = TRUE)
  tp <- cumsum(truth[o]); P <- sum(truth)
  if (P == 0) return(NA_real_)
  precision <- tp / seq_along(tp); recall <- tp / P
  sum(diff(c(0, recall)) * precision)
}

# ---- the authors' exact scoring (bin/process_assignments.R) ------------------
# per-entry pooled confusion, on subset = "full" (all entries) and "nonzero" (counts>0)
.get_confusion <- function(true, est) {
  tp <- sum(est * true); fn <- sum(true) - tp; fp <- sum(est) - tp
  stopifnot(length(est) == length(true))
  tot <- length(est); tn <- tot - tp - fn - fp
  data.frame(TN = tn, FN = fn, FP = fp, TP = tp,
             Precision = tp / (tp + fp), Recall = tp / (tp + fn))
}

# full data.frame form (both "full" and "nonzero" subsets) -- used by the method-panel runner.
# Named score_subsets (NOT score_assignment) to avoid clashing with sim_lib.R::score_assignment,
# which the head-to-head runners also source.
score_subsets <- function(sim, assn, method, sim_label) {
  gt <- assay(sim, "ground_truth"); nz <- assay(sim, "counts") > 0
  cbind(
    data.frame(method = method, sim_label = sim_label, subset = c("full", "nonzero")),
    rbind(.get_confusion(gt, assn), .get_confusion(gt[nz], assn[nz]))
  )
}

# compact vector form c(F1, Precision, Recall) on the "full" subset -- used by the head-to-head
# runners (single self-contained scorer, replacing the previously duplicated inline `score1`).
score_entries <- function(sim, A) {
  gt <- assay(sim, "ground_truth")
  tp <- sum(A * gt); fp <- sum(A) - tp; fn <- sum(gt) - tp
  P <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  R <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  c(F1 = 2 * P * R / (P + R), Precision = P, Recall = R)
}

f1_score <- function(precision, recall) 2 * precision * recall / (precision + recall)
