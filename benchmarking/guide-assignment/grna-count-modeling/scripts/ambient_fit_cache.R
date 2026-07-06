#!/usr/bin/env Rscript
# ============================================================================
# Run fishash+ (contingency_assign, cell_margin="ambient") on ALL guides of ALL real
# datasets and CACHE the fitted ambient-model parameters so downstream per-guide views
# (e.g. the marginal depth-mixed Poisson panel) never re-run the fit.
#
# Per dataset we save results/ambient_ceiling/fit_cache/<ds>_fit.rds :
#   list(dataset, guides, cells, G, C, Rg, Cc, Tn, rho, guide_mean)
# where Rg[g], Cc[c], Tn are fishash's rank-1 Poisson completion of the signal-masked matrix
# (impute_masked_counts on the depth_fix assignment): lambda_{gc} = (Rg[g]/Tn) * Cc[c].
# The (large) assignment matrix is saved separately as <ds>_assigned.rds.
#
# Also rebuilds results/ambient_ceiling/per_guide_ambient_ceiling.csv (adds guide_mean) and
# dataset_summary.csv, identical in spirit to ambient_ceiling_survey.R.
# fit_cache/ is gitignored (large; regenerable by re-running this script).
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(fishash); library(extraDistr); library(sparseMatrixStats)
})
# run from the grna-count-modeling/ folder root
source("scripts/contingency_method.R"); source("scripts/datasets.R")
OUT   <- "results/ambient_ceiling"; CACHE <- file.path(OUT, "fit_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

datasets <- dataset_paths()

per_guide <- list(); per_ds <- list()
for (nm in names(datasets)) {
  p <- datasets[[nm]]; t0 <- proc.time()["elapsed"]
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), nm)); flush.console()
  if (!file.exists(p)) { cat("  MISSING, skipping\n"); next }
  counts <- load_grna_matrix(p)
  G <- nrow(counts); C <- ncol(counts)
  gname <- rownames(counts); if (is.null(gname)) gname <- as.character(seq_len(G))
  cname <- colnames(counts)
  cat(sprintf("  %d guides x %d cells (nnz=%d); fishash+ ...\n", G, C, length(counts@x))); flush.console()

  res <- tryCatch(fishashplus_assign(counts),
                  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) { rm(counts); gc(); next }

  # rank-1 ambient field from the signal-masked completion (the fitted parameters)
  bg <- fishash::impute_masked_counts(counts, res$assigned)
  Rg <- Matrix::rowSums(bg); Cc <- Matrix::colSums(bg); Tn <- sum(bg)
  guide_mean <- Matrix::rowMeans(counts)
  names(Rg) <- gname; names(guide_mean) <- gname; names(Cc) <- cname

  saveRDS(list(dataset = nm, guides = gname, cells = cname, G = G, C = C,
               Rg = Rg, Cc = Cc, Tn = Tn, rho = res$rho, guide_mean = guide_mean),
          file.path(CACHE, paste0(nm, "_fit.rds")))
  saveRDS(res$assigned, file.path(CACHE, paste0(nm, "_assigned.rds")))

  # per-guide ambient ceiling (max unassigned nonzero count), same definition as the survey
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  asg <- as.logical(res$assigned[cbind(i, j)]); asg[is.na(asg)] <- FALSE; amb <- !asg
  ceil  <- tapply(y[amb], factor(i[amb], levels = seq_len(G)), max)
  nz    <- tapply(y,      factor(i,      levels = seq_len(G)), length)
  n_amb <- tapply(amb,    factor(i,      levels = seq_len(G)), sum)
  minA  <- tapply(ifelse(asg, y, NA), factor(i, levels = seq_len(G)),
                  function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
  df <- data.frame(dataset = nm, guide = gname,
                   n_cells_nonzero = as.integer(ifelse(is.na(nz), 0L, nz)),
                   n_assigned = as.integer(ifelse(is.na(nz),0L,nz) - ifelse(is.na(n_amb),0L,n_amb)),
                   ambient_ceiling = as.numeric(ceil), min_assigned = as.numeric(minA),
                   guide_mean = as.numeric(guide_mean), stringsAsFactors = FALSE)
  df <- df[df$n_cells_nonzero > 0, , drop = FALSE]
  per_guide[[nm]] <- df
  per_ds[[nm]] <- data.frame(dataset = nm, guides = nrow(df), cells = C,
    med_ambient_ceiling = median(df$ambient_ceiling, na.rm = TRUE),
    p90_ambient_ceiling = as.numeric(quantile(df$ambient_ceiling, 0.90, na.rm = TRUE)),
    max_ambient_ceiling = max(df$ambient_ceiling, na.rm = TRUE), rho = res$rho)
  cat(sprintf("  cached in %.1fs; med ceiling=%.0f max=%.0f\n", proc.time()["elapsed"] - t0,
              per_ds[[nm]]$med_ambient_ceiling, per_ds[[nm]]$max_ambient_ceiling)); flush.console()
  rm(counts, bg, res, ct); gc()
}
write.csv(do.call(rbind, per_guide), file.path(OUT, "per_guide_ambient_ceiling.csv"), row.names = FALSE)
write.csv(do.call(rbind, per_ds),    file.path(OUT, "dataset_summary.csv"), row.names = FALSE)
cat("\nAll fits cached to", CACHE, "\n")
