#!/usr/bin/env Rscript
# ============================================================================
# Apply fishash+ (contingency_assign, cell_margin="ambient") to every REAL gRNA
# count dataset and tabulate, per guide, the LARGEST count that can plausibly be
# ambient = max{ count[g,c] : count[g,c] > 0 and NOT assigned as real }.
#
# fishash+ = fishash's noise-conditioned contingency test with the cell margin
# supplied from the denoised ambient depth (the "depth_fix" recommendation).
# Output: per-guide CSV + pooled/faceted histogram of the per-guide ambient ceiling.
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2)
})

# run from the grna-count-modeling/ folder root. NOTE: the per-guide ceiling COMPUTE here is superseded
# by scripts/ambient_fit_cache.R (which also caches the fit params); this script is kept for its pooled/
# faceted ambient-ceiling histograms.
source("scripts/contingency_method.R"); source("scripts/datasets.R")
OUT <- "results/ambient_ceiling"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
datasets <- dataset_paths()

per_guide <- list()
per_ds    <- list()

for (nm in names(datasets)) {
  p <- datasets[[nm]]
  t0 <- proc.time()["elapsed"]
  cat(sprintf("[%s] loading %s\n", format(Sys.time(), "%H:%M:%S"), nm)); flush.console()
  if (!file.exists(p)) { cat("  MISSING, skipping\n"); next }
  counts <- readRDS(p)
  counts <- as(counts, "CsparseMatrix")
  storage.mode(counts@x) <- "double"
  G <- nrow(counts); C <- ncol(counts)
  cat(sprintf("  dims %d guides x %d cells, nnz=%d ; running fishash+ ...\n",
              G, C, length(counts@x))); flush.console()

  res <- tryCatch(
    contingency_assign(counts, q = 0.05, refit = 10, min_count = 2,
                       cell_margin = "ambient", tail = "hyper", fdr = "GS"),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next

  # nonzero triplets (guide i, cell j, count y), matching contingency_assign
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  asg <- as.logical(res$assigned[cbind(i, j)]); asg[is.na(asg)] <- FALSE

  # per-guide summaries over its NONZERO entries
  gname <- rownames(counts); if (is.null(gname)) gname <- as.character(seq_len(G))
  amb   <- !asg                                            # ambient (unassigned) nonzeros
  ceil  <- tapply(y[amb], factor(i[amb], levels = seq_len(G)), max)   # largest ambient count / guide
  nz    <- tapply(y,      factor(i,      levels = seq_len(G)), length)# nonzero cells / guide
  n_amb <- tapply(amb,    factor(i,      levels = seq_len(G)), sum)   # ambient nonzeros / guide
  minA  <- tapply(ifelse(asg, y, NA), factor(i, levels = seq_len(G)),
                  function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)) # smallest assigned

  df <- data.frame(
    dataset          = nm,
    guide            = gname,
    n_cells_nonzero  = as.integer(ifelse(is.na(nz), 0L, nz)),
    n_assigned       = as.integer(ifelse(is.na(nz), 0L, nz) - ifelse(is.na(n_amb), 0L, n_amb)),
    ambient_ceiling  = as.numeric(ceil),          # NA if guide has no ambient nonzero (all assigned / all zero)
    min_assigned     = as.numeric(minA),
    stringsAsFactors = FALSE)
  df <- df[df$n_cells_nonzero > 0, , drop = FALSE]  # keep only observed guides
  per_guide[[nm]] <- df

  per_ds[[nm]] <- data.frame(
    dataset = nm, guides = nrow(df), cells = C,
    med_ambient_ceiling = median(df$ambient_ceiling, na.rm = TRUE),
    p90_ambient_ceiling = as.numeric(quantile(df$ambient_ceiling, 0.90, na.rm = TRUE)),
    max_ambient_ceiling = max(df$ambient_ceiling, na.rm = TRUE),
    rho = res$rho)
  cat(sprintf("  done in %.1fs ; median ceiling=%.0f  p90=%.0f  max=%.0f\n",
              proc.time()["elapsed"] - t0,
              per_ds[[nm]]$med_ambient_ceiling, per_ds[[nm]]$p90_ambient_ceiling,
              per_ds[[nm]]$max_ambient_ceiling)); flush.console()
  # incremental save so a late crash doesn't lose earlier work
  saveRDS(list(per_guide = per_guide, per_ds = per_ds), file.path(OUT, "_progress.rds"))
}

all_guides <- do.call(rbind, per_guide)
ds_summary <- do.call(rbind, per_ds)
write.csv(all_guides, file.path(OUT, "per_guide_ambient_ceiling.csv"), row.names = FALSE)
write.csv(ds_summary, file.path(OUT, "dataset_summary.csv"), row.names = FALSE)

# ---- histograms of the per-guide ambient ceiling -------------------------
plotdf <- all_guides[is.finite(all_guides$ambient_ceiling) & all_guides$ambient_ceiling >= 1, ]
plotdf$log10_ceiling <- log10(plotdf$ambient_ceiling)
brks  <- c(1, 2, 3, 5, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
lbls  <- as.character(brks)

# pooled
p_pooled <- ggplot(plotdf, aes(x = ambient_ceiling)) +
  geom_histogram(bins = 45, fill = "#3b6ea5", colour = "white", linewidth = 0.15) +
  scale_x_log10(breaks = brks, labels = lbls) +
  labs(title = "Largest count fishash+ still calls ambient, per guide",
       subtitle = sprintf("Pooled over %d guides across %d real datasets",
                          nrow(plotdf), length(unique(plotdf$dataset))),
       x = "Per-guide ambient ceiling  (max unassigned count, log scale)",
       y = "Number of guides") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "ambient_ceiling_hist_pooled.png"), p_pooled,
       width = 8, height = 5, dpi = 150)

# faceted per dataset
p_facet <- ggplot(plotdf, aes(x = ambient_ceiling)) +
  geom_histogram(bins = 30, fill = "#3b6ea5", colour = "white", linewidth = 0.1) +
  scale_x_log10(breaks = c(1,10,100,1000,10000), labels = c("1","10","100","1k","10k")) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(title = "Per-guide ambient ceiling by dataset (fishash+)",
       x = "Largest count called ambient (log scale)", y = "Guides") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "ambient_ceiling_hist_facet.png"), p_facet,
       width = 13, height = 9, dpi = 150)

cat("\n==== DATASET SUMMARY ====\n")
print(ds_summary, row.names = FALSE)
cat("\nWrote:\n  ", file.path(OUT, "per_guide_ambient_ceiling.csv"),
    "\n  ", file.path(OUT, "dataset_summary.csv"),
    "\n  ", file.path(OUT, "ambient_ceiling_hist_pooled.png"),
    "\n  ", file.path(OUT, "ambient_ceiling_hist_facet.png"), "\n")
