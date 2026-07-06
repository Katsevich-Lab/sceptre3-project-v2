#!/usr/bin/env Rscript
# Show the gRNA count distribution of the replogle-rd7 guide with the LARGEST count called
# ambient by fishash+ (per-guide ambient ceiling = 7: guide 8832_TFAM_P1P2_ENSG00000108064),
# overlaid with the fitted MARGINAL DEPTH-MIXED POISSON ambient model.
#
# Ambient model (same as scripts/barnyard_marginal_fit.R): each cell c contributes an ambient
# count ~ Poisson(lambda_{gc}), lambda_{gc} = a_g d_c = (R_g / T) C_c, where R_g,C_c,T come from
# fishash's rank-1 Poisson completion of the SIGNAL-MASKED matrix (impute_masked_counts on the
# depth_fix assignment). Marginal predicted #cells with count k = sum_c dpois(k, lambda_{gc}).
# Baseline overlay: a single homogeneous Poisson at the guide's mean (ignores depth spread).
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2); library(tidyr) })
source("scripts/contingency_method.R")
source("scripts/datasets.R")   # load_grna_matrix()

OUT <- "results/ambient_ceiling"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GUIDE <- "8832_TFAM_P1P2_ENSG00000108064"

counts <- load_grna_matrix("replogle-rd7")   # as(CsparseMatrix) + storage.mode double
cat(sprintf("replogle-rd7: %d guides x %d cells\n", nrow(counts), ncol(counts)))

## ---- run fishash+ (depth_fix) and extract the denoised rank-1 ambient rate field ----
res <- contingency_assign(counts, q = 0.05, refit = 10, min_count = 2,
                          cell_margin = "ambient", tail = "hyper", fdr = "GS")
bg  <- fishash::impute_masked_counts(counts, res$assigned)    # rank-1 completion, signal masked
Rg  <- Matrix::rowSums(bg); Cc <- Matrix::colSums(bg); Tn <- sum(bg)

g   <- match(GUIDE, rownames(counts)); stopifnot(!is.na(g))
y   <- as.numeric(counts[g, ])                                # observed count in every cell
lam <- (Rg[g] / Tn) * Cc                                      # lambda_{gc} = a_g d_c, all cells
asg <- as.logical(res$assigned[g, ]); asg[is.na(asg)] <- FALSE
amb_ceiling <- max(y[y > 0 & !asg])                           # largest count left unassigned
cat(sprintf("guide %s: %d cells nonzero, %d assigned; ambient ceiling = %d; min assigned = %d\n",
            GUIDE, sum(y > 0), sum(asg), amb_ceiling, min(y[asg])))

## ---- marginal: observed vs depth-mixed Poisson vs single Poisson --------------------
KMAX <- 25L                                                   # individual bins 0..24, then 25+
ks   <- 0:KMAX
mu   <- mean(y); n <- length(y)
obs  <- vapply(ks, function(k) if (k < KMAX) sum(y == k) else sum(y >= k), numeric(1))
dm   <- vapply(ks, function(k) if (k < KMAX) sum(dpois(k, lam)) else sum(1 - ppois(k-1, lam)), numeric(1))
sp   <- vapply(ks, function(k) if (k < KMAX) n*dpois(k, mu)    else n*(1 - ppois(k-1, mu)),   numeric(1))
klab <- ifelse(ks < KMAX, as.character(ks), paste0(KMAX, "+"))

D <- data.frame(k = ks, klab = factor(klab, levels = klab), observed = obs, depthmix = dm, single = sp)
Dl <- pivot_longer(D, c(depthmix, single), names_to = "model", values_to = "pred")
Dl$model <- factor(Dl$model, c("depthmix","single"),
  labels = c("fitted marginal depth-mixed Poisson (a_g d_c)", "single Poisson (guide mean)"))

p <- ggplot(D, aes(klab, observed)) +
  geom_col(fill = "grey82", width = 0.72) +
  geom_vline(xintercept = amb_ceiling + 1, linetype = "dashed", colour = "grey40") +
  annotate("text", x = amb_ceiling + 1, y = max(obs), hjust = -0.08, vjust = 1,
           label = sprintf("ambient ceiling = %d", amb_ceiling), size = 4, colour = "grey30") +
  geom_line(data = Dl, aes(klab, pred, colour = model, group = model), linewidth = 0.8) +
  geom_point(data = Dl, aes(klab, pred, colour = model), size = 2.2) +
  scale_y_continuous(trans = scales::trans_new("log1p", log1p, expm1),
                     breaks = c(0, 1, 10, 100, 1000, 10000, 1e5), labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = c("fitted marginal depth-mixed Poisson (a_g d_c)" = "#D55E00",
                                 "single Poisson (guide mean)" = "#0072B2")) +
  labs(title = sprintf("replogle-rd7  %s", GUIDE),
       subtitle = sprintf("guide with the largest count fishash+ still calls ambient (ceiling %d).  Depth-mixed Poisson = ambient model; observed above ceiling = signal.", amb_ceiling),
       x = "gRNA UMI count", y = "number of cells (log1p; axis floor at 0)", colour = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "top",
       axis.text.x = element_text(size = 9))
ggsave(file.path(OUT, "replogle_top_ambient_guide_fit.png"), p, width = 11, height = 6, dpi = 140)

## ---- goodness-of-fit over the AMBIENT region (k <= ceiling) -------------------------
amb_bins <- ks <= amb_ceiling
chi_dm <- sum((obs[amb_bins] - dm[amb_bins])^2 / pmax(dm[amb_bins], 1e-6))
chi_sp <- sum((obs[amb_bins] - sp[amb_bins])^2 / pmax(sp[amb_bins], 1e-6))
cat(sprintf("\nmarginal var/mean = %.2f ; ambient-region (k<=%d) chi^2:  depth-mix = %.1f   single = %.1f\n",
            var(y)/mean(y), amb_ceiling, chi_dm, chi_sp))
cat("\nobserved vs models (counts of cells):\n")
print(data.frame(k = klab, observed = obs, depthmix = round(dm,1), single = round(sp,1)), row.names = FALSE)
cat("\nwrote", file.path(OUT, "replogle_top_ambient_guide_fit.png"), "\n")
