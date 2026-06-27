#!/usr/bin/env Rscript
# =============================================================================
# sim_lib.R  --  shared foundation for the simulation-framework investigation
#
# Three things every downstream script needs:
#   (1) METHOD REGISTRY  run_methods(): run the assignment-method panel on a
#       guides x cells count matrix, return a named list of logical assignment
#       matrices.  (R methods only; geomux runs out-of-process -- see geomux_*.)
#   (2) SCORING          score_assignment() / score_panel(): per-guide
#       precision / recall / Jaccard vs a ground-truth matrix, plus a pooled
#       realized-FDR / recall (the FDR-control view).
#   (3) REALISM GATE     realism_stats(): method-agnostic summaries that decide
#       whether a simulated matrix "looks like" real data.  Two kinds:
#         - ambient-SHAPE   (needs a known ambient mask): rank-1-conditional
#           Fano (index of dispersion) + max ambient count.  Real direct-capture
#           ambient is Poisson (Fano ~ 1); this is the check the lab's
#           comprehensive sim fails.
#         - STRUCTURAL      (one fixed estimator on any matrix): library/ambient
#           depth ratio, ambient-pool Gini, MOI, signal fraction, mode-gap.
#
# Paths follow the repo convention: HERE resolves to the folder root.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(sparseMatrixStats)
})

.here <- function() {
  # 1. explicit override; 2. marker file relative to cwd (the usual case: run
  # from the folder root); 3. derive from the running script's --file path.
  if (nzchar(Sys.getenv("SIM_HERE"))) return(normalizePath(Sys.getenv("SIM_HERE")))
  if (file.exists(file.path(getwd(), "scripts", "contingency_method.R"))) return(normalizePath(getwd()))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  h <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
  if (basename(h) == "scripts") dirname(h) else h
}
HERE <- .here()
GA   <- normalizePath(file.path(HERE, ".."))
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
source(file.path(HERE, "scripts", "contingency_method.R"))

# -----------------------------------------------------------------------------
# (1) METHOD REGISTRY
# -----------------------------------------------------------------------------
# Each entry maps a method name -> function(counts) returning a logical
# guides x cells assignment matrix (sparse).  geomux is handled separately
# because it needs its own python venv; run_methods() can splice geomux calls
# in via the `geomux_calls` argument (a guides x cells logical matrix or NULL).
# bulletproof coercion to a logical CsparseMatrix from any sparse type
# (numeric dgC, logical lgC, or pattern ngC/ngR as sceptre returns).
.as_lgl <- function(A, dn = NULL) {
  Tt <- as(as(A, "CsparseMatrix"), "TsparseMatrix")
  keep <- if (.hasSlot(Tt, "x")) Tt@x != 0 else rep(TRUE, length(Tt@i))
  Matrix::sparseMatrix(i = Tt@i[keep] + 1L, j = Tt@j[keep] + 1L, x = TRUE,
                       dims = dim(Tt), dimnames = if (!is.null(dn)) dn else dimnames(Tt))
}

# ---- sceptre mixture (the current/"legacy" SCEPTRE assignment) --------------
# Builds a minimal sceptre object (dummy 1-gene response) and runs the mixture
# assignment with log grna-covariates -- exactly run_sceptre.R's call.
sceptre_assign <- function(counts, moi = "high") {
  gm <- as(counts, "CsparseMatrix"); G <- nrow(gm); C <- ncol(gm)
  if (is.null(rownames(gm))) rownames(gm) <- paste0("g", seq_len(G))
  if (is.null(colnames(gm))) colnames(gm) <- paste0("c", seq_len(C))
  set.seed(1)
  resp <- Matrix::sparseMatrix(i = rep(1L, 5), j = sample(C, 5), x = 1,
                               dims = c(1, C), dimnames = list("gene1", colnames(gm)))
  tdf <- data.frame(grna_id = rownames(gm), grna_target = paste0("t", seq_len(G)),
                    stringsAsFactors = FALSE)
  tdf$grna_target[1] <- "non-targeting"
  suppressMessages(utils::capture.output({
    obj <- sceptre::import_data(response_matrix = resp, grna_matrix = gm,
                                grna_target_data_frame = tdf, moi = moi)
    obj <- sceptre::set_analysis_parameters(obj, formula = ~1)
    obj <- sceptre::assign_grnas(obj, method = "mixture",
             formula_object = ~ log(grna_n_umis + 1) + log(grna_n_nonzero + 1))
  }))
  A <- sceptre::get_grna_assignments(obj)
  # sceptre returns guide rownames but drops cell colnames, preserving cell ORDER
  # (we skip run_qc, so no cells are dropped) -- exactly as the pipeline's
  # run_sceptre.R relies on.  Re-attach cell names by position, align guides by name.
  stopifnot(ncol(A) == C)
  if (is.null(colnames(A))) colnames(A) <- colnames(gm)
  A[rownames(gm), colnames(gm)]
}

# Fixed global UMI-count threshold: assign (g,c) iff count >= t.  The simplest
# possible baseline (Louis: "y>>20 perturbed, y<<20 not; only y~10 is hard").
.fixed_thresh <- function(counts, tcut) {
  ct <- as(counts, "TsparseMatrix"); keep <- ct@x >= tcut
  Matrix::sparseMatrix(i = ct@i[keep] + 1L, j = ct@j[keep] + 1L, x = TRUE,
                       dims = dim(counts), dimnames = dimnames(counts))
}

run_methods <- function(counts, methods = c("thresh10","ambient","otsu","valley","fishash","depth_fix","sceptre"),
                        q = 0.05, external = list(), verbose = TRUE) {
  counts <- as(counts, "CsparseMatrix")
  dn <- dimnames(counts)
  out <- list()
  tt <- function(nm, expr) { t0 <- proc.time()[3]; v <- force(expr)
    if (verbose) cat(sprintf("    [%-9s] %5.1fs\n", nm, proc.time()[3]-t0)); v }
  for (m in methods) {
    out[[m]] <- switch(m,
      ambient   = tt("ambient",  .as_lgl(ambient_test_assign(counts, q = q,
                       model = "hypergeometric", n_iter = 1)$assignment_matrix, dn)),
      otsu      = tt("otsu",     .as_lgl(assign_by_threshold(counts, otsu_threshold_log1p)$assignment_matrix, dn)),
      valley    = tt("valley",   .as_lgl(assign_by_threshold(counts, smoothed_valley_threshold)$assignment_matrix, dn)),
      fishash   = tt("fishash",  .as_lgl(SummarizedExperiment::assay(
                       fishash::fishash(counts, padj_cutoff = q), "assigned"), dn)),
      depth_fix = tt("depth_fix",.as_lgl(contingency_assign(counts, q = q, cell_margin = "ambient",
                       tail = "hyper", fdr = "GS")$assigned, dn)),
      sceptre   = tt("sceptre",  .as_lgl(sceptre_assign(counts), dn)),
      {  # default: fixed-threshold methods named "thresh<N>" (e.g. thresh10)
        if (grepl("^thresh[0-9]+$", m)) tt(m, .as_lgl(.fixed_thresh(counts, as.integer(sub("thresh","",m))), dn))
        else stop("unknown method: ", m)
      })
  }
  # splice in out-of-process methods (geomux, crispat) passed as logical matrices
  for (nm in names(external)) if (!is.null(external[[nm]])) out[[nm]] <- .as_lgl(external[[nm]], dn)
  out
}

# -----------------------------------------------------------------------------
# (2) SCORING
# -----------------------------------------------------------------------------
# Per-guide precision/recall/Jaccard vs a logical truth matrix, plus a pooled
# realized FDR and recall.  Guides with no truth and no calls are dropped from
# the per-guide means (NA); pooled FDR uses all calls.
score_assignment <- function(A, truth) {
  A <- as(A, "lgCMatrix"); truth <- as(truth, "lgCMatrix")
  G <- nrow(A)
  Ar <- as(A, "RsparseMatrix"); Tr <- as(truth, "RsparseMatrix")
  per <- do.call(rbind, lapply(seq_len(G), function(g) {
    p <- which(as.logical(Ar[g, ])); t <- which(as.logical(Tr[g, ]))
    inter <- length(intersect(p, t)); uni <- length(union(p, t))
    data.frame(guide = g,
               n_pred = length(p), n_true = length(t),
               jac  = if (uni == 0) NA_real_ else inter / uni,
               prec = if (length(p)) inter / length(p) else NA_real_,
               rec  = if (length(t)) inter / length(t) else NA_real_)
  }))
  npred <- sum(A); ntrue <- sum(truth); tp <- sum(A & truth)
  list(per_guide = per,
       summary = data.frame(
         jaccard   = mean(per$jac,  na.rm = TRUE),
         precision = mean(per$prec, na.rm = TRUE),
         recall    = mean(per$rec,  na.rm = TRUE),
         fdr_pooled    = if (npred) (npred - tp) / npred else 0,
         recall_pooled = if (ntrue) tp / ntrue else NA_real_,
         n_pred = npred, n_true = ntrue))
}

# Score a named list of assignment matrices against one truth matrix.
score_panel <- function(assignments, truth) {
  do.call(rbind, lapply(names(assignments), function(m) {
    s <- score_assignment(assignments[[m]], truth)$summary; s$method <- m; s
  }))
}

# -----------------------------------------------------------------------------
# (3) REALISM GATE
# -----------------------------------------------------------------------------
.gini <- function(x) { x <- sort(x); n <- length(x); if (n == 0 || sum(x) == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n }

# Rank-1-conditional index of dispersion (Fano) of the AMBIENT counts.
#   Given a logical ambient_mask (TRUE where the (g,c) count is provably/known
#   ambient), fit a rank-1 Poisson mean  lambda_gc = a_g * b_c / S  by moments
#   over the masked entries, then Fano = mean_{lambda>=floor} (y-lambda)^2 /
#   lambda.  Conditioning on the rank-1 mean removes mean heterogeneity (gamma_g
#   spans, kappa_c spans) so Fano isolates TRUE overdispersion: Poisson -> ~1.
# Returns Fano, the max ambient count, and the fraction of ambient entries used.
ambient_dispersion <- function(counts, ambient_mask, lambda_floor = 0.3) {
  counts <- as(counts, "CsparseMatrix"); ambient_mask <- as(ambient_mask, "lgCMatrix")
  # ambient counts as a (possibly zero-laden) matrix over masked entries
  Amb <- counts * ambient_mask                       # ambient counts where masked, else 0
  a_g <- Matrix::rowSums(Amb); b_c <- Matrix::colSums(Amb); S <- sum(Amb)
  nz <- Amb@x[Amb@x > 0]                              # nonzero ambient counts
  vmr_nz <- if (length(nz) > 1 && mean(nz) > 0) var(nz) / mean(nz) else NA_real_  # robust secondary
  if (S == 0) return(list(fano = NA_real_, max_ct = 0, n_used = 0, vmr_nz = vmr_nz))
  mm <- as(ambient_mask, "TsparseMatrix"); i <- mm@i + 1L; j <- mm@j + 1L  # all masked cells (incl zeros)
  y  <- counts[cbind(i, j)]
  lam <- a_g[i] * b_c[j] / S
  keep <- is.finite(lam) & lam >= lambda_floor
  fano <- if (sum(keep) < 50) NA_real_ else mean((y[keep] - lam[keep])^2 / lam[keep])
  list(fano = fano, max_ct = max(y, 0), n_used = sum(keep), vmr_nz = vmr_nz)
}

# Structural summaries via ONE fixed estimator (the ambient test) so real and
# simulated matrices are judged through the same lens.  Mirrors
# ambient_survey_compute.R::metrics() but trimmed to the gate quantities.
structural_stats <- function(counts, maxc = 25000, seed = 1) {
  gm <- as(counts, "CsparseMatrix"); storage.mode(gm@x) <- "double"
  if (ncol(gm) > maxc) { set.seed(seed); gm <- gm[, sort(sample(ncol(gm), maxc))] }
  gm <- gm[, Matrix::colSums(gm) > 0, drop = FALSE]
  lib <- Matrix::colSums(gm)
  A <- suppressMessages(ambient_test_assign(gm, q = 0.05, model = "hypergeometric",
                                            n_iter = 1)$assignment_matrix)
  A <- as(A, "CsparseMatrix") * 1
  sig_cell  <- Matrix::colSums(gm * A)
  amb_depth <- pmax(lib - sig_cell, 0)
  moi <- Matrix::colSums(A)
  tot_g <- Matrix::rowSums(gm); sig_g <- Matrix::rowSums(gm * A); amb_g <- tot_g - sig_g
  well <- amb_g >= 10
  dr <- lib / pmax(amb_depth, 0.5)
  assigned <- sig_cell > 0
  data.frame(
    n_guides = nrow(gm), n_cells = ncol(gm),
    moi_med = median(moi), assigned_frac = mean(assigned),
    signal_frac_med = median((sig_cell / lib)[assigned]),
    depth_ratio_med = median(dr), depth_ratio_p90 = as.numeric(quantile(dr, .90)),
    depth_cor = suppressWarnings(cor(log1p(lib), log1p(amb_depth))),
    comp_gini = if (sum(well) > 1) .gini(amb_g[well]) else NA_real_)
}

# Per-guide bimodality / mode-gap (the realism_separation axis): on how many
# guides is the perturbed mode visibly separated, and by how much (log units)?
separation_stats <- function(counts, min_nonzero = 30, max_guides = 400) {
  gm <- as(counts, "RsparseMatrix")
  G <- nrow(gm); idx <- seq_len(G)
  if (G > max_guides) { set.seed(1); idx <- sort(sample(G, max_guides)) }
  gaps <- c(); ok <- logical(0)
  for (g in idx) {
    cv <- as.numeric(gm[g, ]); if (sum(cv > 0) < min_nonzero) next
    v <- smoothed_valley_threshold(cv)
    ok <- c(ok, isTRUE(v$ok))
    if (isTRUE(v$ok)) gaps <- c(gaps, log1p(v$mode2) - log1p(v$mode1))
  }
  data.frame(frac_bimodal = if (length(ok)) mean(ok) else NA_real_,
             mode_gap_med = if (length(gaps)) median(gaps) else NA_real_)
}

# One call: full realism row for a matrix.  ambient_mask optional (only then are
# the ambient-shape stats computed).
realism_stats <- function(counts, ambient_mask = NULL, label = NA_character_) {
  st <- structural_stats(counts)
  sp <- separation_stats(counts)
  ad <- if (!is.null(ambient_mask)) ambient_dispersion(counts, ambient_mask)
        else list(fano = NA_real_, max_ct = NA_real_, vmr_nz = NA_real_)
  data.frame(label = label, st, sp,
             ambient_fano = ad$fano, ambient_vmr_nz = ad$vmr_nz, ambient_max = ad$max_ct)
}

# -----------------------------------------------------------------------------
# helpers shared by the simulators
# -----------------------------------------------------------------------------
# Build a logical truth matrix from integration matrix Z and observed counts:
# the *recoverable* truth is "integrated AND observed (count>0)".
truth_observed <- function(Z, counts) as(as(Z, "lgCMatrix") & (as(counts, "CsparseMatrix") > 0), "lgCMatrix")

# --- dataset on-disk layout (shared by Model A/B builders and the orchestrator) -
SIMFW  <- function() file.path(HERE, "results", "sim_framework")
DSDIR  <- function() file.path(SIMFW(), "datasets")
save_dataset <- function(id, counts, Z, meta) {
  d <- file.path(DSDIR(), id); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  saveRDS(as(counts, "CsparseMatrix"), file.path(d, "counts.rds"))
  saveRDS(as(Z, "lgCMatrix"),          file.path(d, "Z.rds"))
  saveRDS(meta,                        file.path(d, "meta.rds"))
  invisible(d)
}
load_dataset <- function(id) {
  d <- file.path(DSDIR(), id)
  list(counts = readRDS(file.path(d, "counts.rds")), Z = readRDS(file.path(d, "Z.rds")),
       meta = readRDS(file.path(d, "meta.rds")), id = id, dir = d)
}
# simple Dirichlet draw (rgamma-based) so we avoid an extra package dependency
rdirichlet1 <- function(alpha) { g <- rgamma(length(alpha), shape = alpha, rate = 1); g / sum(g) }

cat("sim_lib.R loaded: run_methods(), score_assignment/score_panel(), realism_stats().\n")
