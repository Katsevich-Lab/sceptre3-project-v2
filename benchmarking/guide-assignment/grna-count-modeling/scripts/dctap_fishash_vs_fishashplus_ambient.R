#!/usr/bin/env Rscript
# ============================================================================
# DC-TAP K562 high-MOI: per-guide fitted-ambient overlap plots, fishash vs fishash+.
#
# For each positive-control guide, overlay on the observed gRNA-count histogram the
# FITTED AMBIENT count distribution implied by each method's rank-1 ambient rate field:
#   fishash  = contingency_assign(cell_margin="observed")
#   fishash+ = contingency_assign(cell_margin="ambient")   (depth_fix)
#
# CRITICAL CONSTRUCTION: the fitted "number of ambient cells at count k" for guide g is
#   pred_ambient(k) = sum over { cells c : NOT assigned[g,c] } of dpois(k, lambda_gc)
# i.e. summed over that method's AMBIENT (unassigned) cells ONLY. This integrates to
# (#ambient cells for g), not the total cell count. Summing over all cells would treat
# every signal cell as ambient and over-predict the low counts. We verify the integral.
#
# lambda_gc from the returned background matrix bg of impute_masked_counts:
#   Rg = rowSums(bg); Cc = colSums(bg); Tn = sum(bg);  lambda_gc = (Rg[g]/Tn) * Cc[c]
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2); library(tidyr); library(scales)
})

HERE <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
source(file.path(HERE, "contingency_method.R"))

OUT_DIR <- normalizePath(file.path(HERE, "..", "results", "ambient_ceiling"), mustWork = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATA <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey/dctap_k562_highmoi_GSE303901/sceptre/grna_matrix_aligned.rds")
KMAX <- 60L

# ---- positive-control guides (12) ----------------------------------------
pc <- list(
  GATA1 = c("GATA1-TSS-1","GATA1-TSS-2","e-GATA1-1","e-GATA1-2"),
  HDAC6 = c("HDAC6-TSS-1","HDAC6-TSS-2","e-HDAC6-1","e-HDAC6-2"),
  MYC   = c("MYC-TSS-1","MYC-TSS-2","MYC-e2-1","MYC-e2-2")
)
pc_df <- do.call(rbind, lapply(names(pc), function(g)
  data.frame(gene = g, guide = pc[[g]], stringsAsFactors = FALSE)))

# ---- load matrix ----------------------------------------------------------
counts <- as(readRDS(DATA), "CsparseMatrix")
storage.mode(counts@x) <- "double"
rn <- rownames(counts)
stopifnot(all(pc_df$guide %in% rn))
if (!all(pc_df$guide %in% rn)) {
  cat("Available rownames:\n"); print(rn)
  stop("Some positive-control guides not found.")
}
cat(sprintf("Loaded matrix: %d guides x %d cells\n", nrow(counts), ncol(counts)))

# ---- run both methods -----------------------------------------------------
# Both return $assigned; we separately recompute the FINAL background field using the
# same assigned-mask that the method converged to, so lambda_gc matches the method's fit.
run_method <- function(cell_margin) {
  cat(sprintf("Running contingency_assign(cell_margin='%s') ...\n", cell_margin))
  t0 <- Sys.time()
  res <- contingency_assign(counts, q = 0.05, refit = 10, min_count = 2,
                            cell_margin = cell_margin, tail = "hyper", fdr = "GS")
  cat(sprintf("  done in %.1fs; assigned entries = %d\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(res$assigned)))
  assigned <- as(res$assigned, "CsparseMatrix")
  # Final background from the converged assigned mask (rank-1 ambient completion).
  bg <- fishash::impute_masked_counts(counts, assigned)
  Rg <- Matrix::rowSums(bg); Cc <- Matrix::colSums(bg); Tn <- sum(bg)
  list(assigned = assigned, Rg = Rg, Cc = Cc, Tn = Tn)
}
fh  <- run_method("observed")   # fishash
fhp <- run_method("ambient")    # fishash+

# ---- per-guide fitted-ambient curves --------------------------------------
# For guide g: over unassigned cells, lambda_gc = (Rg[g]/Tn)*Cc[c].
# pred_ambient(k) = sum_{c: !assigned[g,c]} dpois(k, lambda_gc), k = 0..KMAX with 60+ top bin.
kbins <- 0:KMAX
counts_full <- as.matrix(counts[pc_df$guide, , drop = FALSE])  # 12 x C dense (small)

fitted_ambient_curve <- function(m, guide) {
  gi <- match(guide, rn)
  assigned_g <- as.logical(m$assigned[gi, ])          # length C
  amb <- !assigned_g                                   # ambient (unassigned) cells
  lam <- (m$Rg[gi] / m$Tn) * m$Cc                      # length C, per-cell ambient rate
  lam_amb <- lam[amb]
  # pred over k=0..KMAX-1 exactly; final bin KMAX = 60+ gets upper-tail mass (>=60).
  pred <- numeric(length(kbins))
  for (ki in seq_along(kbins)) {
    k <- kbins[ki]
    if (k < KMAX) pred[ki] <- sum(dpois(k, lam_amb))
    else          pred[ki] <- sum(stats::ppois(KMAX - 1L, lam_amb, lower.tail = FALSE)) # P(X>=60)
  }
  list(pred = pred, n_ambient = sum(amb), total_pred = sum(pred),
       ambient_counts = counts_full[match(guide, pc_df$guide), amb])
}

# ---- verification: integrates-to-#ambient-cells ---------------------------
chk_g <- pc_df$guide[1]
chk <- fitted_ambient_curve(fh, chk_g)
cat(sprintf("\n[CHECK] guide '%s' (fishash): sum_k pred_ambient = %.2f  vs  #ambient cells = %d  (ratio %.4f)\n",
            chk_g, chk$total_pred, chk$n_ambient, chk$total_pred / chk$n_ambient))
chk2 <- fitted_ambient_curve(fhp, chk_g)
cat(sprintf("[CHECK] guide '%s' (fishash+): sum_k pred_ambient = %.2f  vs  #ambient cells = %d  (ratio %.4f)\n\n",
            chk_g, chk2$total_pred, chk2$n_ambient, chk2$total_pred / chk2$n_ambient))
# Note: the 60+ tail bin makes the pred integral ~= #ambient (Poisson mass fully summed).

# ---- observed histogram over ALL cells ------------------------------------
obs_hist <- function(guide) {
  x <- counts_full[match(guide, pc_df$guide), ]
  xb <- pmin(as.integer(round(x)), KMAX)               # cap at 60+
  tabulate(xb + 1L, nbins = KMAX + 1L)                 # index 1 == count 0
}

# ---- ambient ceiling: largest count k with >=1 UNASSIGNED cell ------------
ambient_ceiling <- function(m, guide) {
  gi <- match(guide, rn)
  assigned_g <- as.logical(m$assigned[gi, ])
  x <- counts_full[match(guide, pc_df$guide), ]
  amb_counts <- x[!assigned_g]
  amb_counts <- amb_counts[amb_counts > 0]
  if (length(amb_counts) == 0) return(0)
  max(amb_counts)
}

# ---- assemble plotting frame + overlap CSV --------------------------------
plot_rows <- list()
csv_rows  <- list()
for (r in seq_len(nrow(pc_df))) {
  g <- pc_df$guide[r]; gene <- pc_df$gene[r]
  h  <- obs_hist(g)
  cf <- fitted_ambient_curve(fh,  g)
  cp <- fitted_ambient_curve(fhp, g)

  facet_lab <- sprintf("%s  (%s)", g, gene)
  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    guide = g, gene = gene, facet = facet_lab,
    count = kbins, observed = h, fishash = cf$pred, `fishash+` = cp$pred,
    check.names = FALSE, stringsAsFactors = FALSE)

  # histogram-intersection overlap of the two fitted-ambient curves (0..1)
  inter <- sum(pmin(cf$pred, cp$pred))
  overlap <- inter / mean(c(cf$total_pred, cp$total_pred))
  csv_rows[[length(csv_rows) + 1]] <- data.frame(
    guide = g, gene = gene,
    n_ambient_fishash = cf$n_ambient, n_ambient_fishashplus = cp$n_ambient,
    total_pred_fishash = cf$total_pred, total_pred_fishashplus = cp$total_pred,
    hist_intersection = overlap,
    ambient_ceiling_fishash = ambient_ceiling(fh, g),
    ambient_ceiling_fishashplus = ambient_ceiling(fhp, g),
    stringsAsFactors = FALSE)
}
pdat <- do.call(rbind, plot_rows)
overlap_df <- do.call(rbind, csv_rows)

# stable facet ordering: GATA1, HDAC6, MYC then guide order given
facet_levels <- unlist(lapply(seq_len(nrow(pc_df)), function(r)
  sprintf("%s  (%s)", pc_df$guide[r], pc_df$gene[r])))
pdat$facet <- factor(pdat$facet, levels = facet_levels)

# ---- write CSV ------------------------------------------------------------
csv_path <- file.path(OUT_DIR, "dctap_fishash_vs_fishashplus_ambient_overlap.csv")
write.csv(overlap_df, csv_path, row.names = FALSE)
cat("Wrote overlap CSV:\n  ", csv_path, "\n", sep = "")
print(overlap_df, row.names = FALSE)

# ---- figure ---------------------------------------------------------------
long <- tidyr::pivot_longer(pdat, cols = c("fishash", "fishash+"),
                            names_to = "method", values_to = "pred")
long$method <- factor(long$method, levels = c("fishash", "fishash+"))

log1p_trans <- scales::trans_new("log1p", transform = log1p, inverse = expm1)

p <- ggplot(pdat, aes(x = count)) +
  geom_col(aes(y = observed), fill = "grey78", colour = NA, width = 1) +
  geom_line(data = long, aes(x = count, y = pred, colour = method), linewidth = 0.7) +
  facet_wrap(~ facet, ncol = 4) +
  scale_colour_manual(values = c("fishash" = "#0072B2", "fishash+" = "#c0392b"),
                      name = NULL) +
  scale_y_continuous(trans = log1p_trans, breaks = c(0, 1, 10, 100, 1000)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60), labels = c("0","20","40","60+")) +
  labs(
    title = "DC-TAP high-MOI: fitted ambient count distribution, fishash vs fishash+",
    subtitle = "Grey = observed gRNA-count histogram (all cells). Lines = fitted ambient cells at each count, summed over each method's UNASSIGNED cells only.",
    x = "gRNA count", y = "number of cells (log1p)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 9))

png_path <- file.path(OUT_DIR, "dctap_fishash_vs_fishashplus_ambient_fit.png")
ggsave(png_path, p, width = 13, height = 9, dpi = 140)
cat("Wrote figure:\n  ", png_path, "\n", sep = "")
