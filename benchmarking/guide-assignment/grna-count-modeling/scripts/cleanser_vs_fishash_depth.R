#!/usr/bin/env Rscript
# ============================================================================
# Per-cell ambient depth: fishash+ (rank-1 completion) vs CLEANSER (sum of <=2 UMIs).
# Scatter per dataset.
#   fishash+ depth  d_c^F = Cc[c]                 (cached rank-1 cell factor; colSums of noise)
#   CLEANSER depth  d_c^C = colSums(counts * (counts <= 2))   (sub-threshold, signal-free UMIs)
# fishash+ depth read from cache (no refit); CLEANSER depth is a colSums over the raw matrix.
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
set.seed(1)
OUT   <- "results/ambient_ceiling"; CACHE <- file.path(OUT, "fit_cache")

source("scripts/datasets.R")
paths <- dataset_paths()

NPLOT <- 12000L                      # per-dataset cap for the scatter (corr uses ALL cells)
pts <- list(); summ <- list()
for (nm in names(paths)) {
  fit <- readRDS(file.path(CACHE, paste0(nm, "_fit.rds")))
  counts <- load_grna_matrix(paths[[nm]])
  Cc <- fit$Cc; if (!is.null(colnames(counts)) && !is.null(names(Cc)))
    Cc <- Cc[match(colnames(counts), names(Cc))]
  cl <- kappa_leq2(counts)                                # CLEANSER per-cell depth (sum of <=2)
  fh <- as.numeric(Cc)                                    # fishash+ per-cell depth (cached)
  keep <- is.finite(cl) & is.finite(fh) & cl > 0 & fh > 0
  cr_p <- suppressWarnings(cor(log10(fh[keep]), log10(cl[keep])))
  cr_s <- suppressWarnings(cor(fh[keep], cl[keep], method = "spearman"))
  summ[[nm]] <- data.frame(dataset = nm, cells = ncol(counts), n_both_pos = sum(keep),
    cl_zero = sum(cl == 0), pearson_log = cr_p, spearman = cr_s,
    med_ratio_cl_over_fh = median((cl[keep]) / (fh[keep])))
  idx <- which(keep); if (length(idx) > NPLOT) idx <- sample(idx, NPLOT)
  pts[[nm]] <- data.frame(dataset = nm, fish = fh[idx], clean = cl[idx])
  cat(sprintf("%-20s r_log=%.3f  spearman=%.3f  med(CL/FH)=%.2f  (cl_zero=%d)\n",
              nm, cr_p, cr_s, summ[[nm]]$med_ratio_cl_over_fh, summ[[nm]]$cl_zero)); flush.console()
  rm(counts); gc()
}
P <- do.call(rbind, pts); S <- do.call(rbind, summ)
write.csv(S, file.path(OUT, "cleanser_vs_fishash_depth_summary.csv"), row.names=FALSE)

lev <- S$dataset[order(S$pearson_log)]                    # worst agreement first
P$dataset <- factor(P$dataset, levels = lev)
lab <- setNames(sprintf("%s  (r=%.2f)", S$dataset, S$pearson_log), S$dataset)

p <- ggplot(P, aes(fish, clean)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.15, size = 0.3, colour = "#3b6ea5") +
  facet_wrap(~ dataset, scales = "free", ncol = 4, labeller = labeller(dataset = lab)) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Per-cell ambient depth: fishash+ vs CLEANSER (sum of <=2 UMIs)",
       subtitle = "dashed = y=x.  Below the line = CLEANSER underestimates depth (its <=2 cap discards ambient UMIs >2).  facet r = Pearson on log10.",
       x = "fishash+ ambient depth  d_c  (rank-1 cell factor, log)",
       y = "CLEANSER ambient depth (sum of counts <=2, log)") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "cleanser_vs_fishash_depth_scatter.png"), p, width = 14, height = 10, dpi = 150)

cat("\n=== summary ===\n"); print(S, row.names = FALSE)
cat("\nwrote:\n  ", file.path(OUT, "cleanser_vs_fishash_depth_scatter.png"),
    "\n  ", file.path(OUT, "cleanser_vs_fishash_depth_summary.csv"), "\n")
