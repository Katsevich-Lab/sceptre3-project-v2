## Clean RAMAC ambient-mode fit for the collaborator writeup: observed counts vs the depth-mixed
## Poisson with the denoised ambient rate a_g * d_c (NO crude-depth comparison). Shows the model
## explains counts 0-1 but leaves a real excess at count 2 (obs 53 vs predicted ~9). Reads cached
## per-cell RAMAC counts + denoised rank-1 (a, d); no reprocessing of the raw matrix.
suppressPackageStartupMessages({library(ggplot2)})
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

r  <- readRDS("results/replogle_ambient_poisson/ramac_depthfix_res.rds")
ad <- readRDS("results/replogle_ambient_poisson/replogle_rank1_ad.rds")
v  <- r$v; lam <- ad$a[r$RAMAC] * ad$d              # per-cell ambient rate a_RAMAC * d_c
amb <- which(v <= 7)                                # RAMAC's ambient cells (below the empty gap)

K <- 5L
ks <- 0:K
obs  <- vapply(ks, function(k) if (k < K) sum(v[amb] == k) else sum(v[amb] >= k), numeric(1))
pred <- vapply(ks, function(k) if (k < K) sum(dpois(k, lam[amb])) else sum(1 - ppois(k-1, lam[amb])), numeric(1))
df <- data.frame(klab = factor(ifelse(ks < K, as.character(ks), paste0(K, "+")),
                               levels = c(as.character(0:(K-1)), paste0(K, "+"))),
                 observed = obs, model = pred)

p <- ggplot(df, aes(klab, observed)) +
  geom_col(fill = "grey80", width = 0.66) +
  geom_line(aes(y = model, group = 1), colour = "#D55E00", linewidth = 0.9) +
  geom_point(aes(y = model), colour = "#D55E00", size = 3) +
  annotate("segment", x = 3, xend = 3, y = 12, yend = 45, colour = "grey30",
           arrow = grid::arrow(length = grid::unit(0.16, "cm"), ends = "both")) +
  annotate("text", x = 3.35, y = 24, hjust = 0, size = 3.5, colour = "grey20",
           label = "count 2:\n53 observed\n~9 predicted") +
  scale_y_log10() +
  labs(title = "RAMAC's ambient mode leaves a real excess at count 2",
       subtitle = "Replogle guide RAMAC, cells below its empty gap. Bars = observed; orange = depth-mixed Poisson (a_g d_c).\nThe model fits counts 0-1; a separate signal mode (70-6072 UMIs) lies above an empty gap.",
       x = "gRNA UMI count", y = "number of cells (log)") +
  theme_bw(base_size = 11.5)
ggsave(file.path(OUT, "ramac_ambient_fit.png"), p, width = 8.5, height = 4.8, dpi = 130)
cat(sprintf("count-2: obs=%d model=%.1f (%.1fx excess)\n", obs[3], pred[3], obs[3]/pred[3]))
cat("wrote ramac_ambient_fit.png\n")
