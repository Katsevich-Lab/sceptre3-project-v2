#!/usr/bin/env Rscript
# Compare depth_fix vs fishash at d_sigma=0.05 (primary, snapshotted in
# scores_B_dsig005.csv) vs the current scores_all.csv (which now holds the
# d_sigma=0.3 sweep results for Model B).  Outputs:
#   results/sim_framework/dsigma_sweep.csv -- per-regime delta table
#   console summary -- overall depth_fix-vs-fishash gap at both settings
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
OUT <- file.path(getwd(), "results", "sim_framework")
pri <- read.csv(file.path(OUT, "scores_B_dsig005.csv"))             # d_sigma=0.05 Model B
all <- read.csv(file.path(OUT, "scores_all.csv"))                   # current = sweep run
swe <- all[all$model == "B", ]                                       # Model B at d_sigma=0.30

bind <- bind_rows(
  pri %>% mutate(d_sigma = 0.05),
  swe %>% mutate(d_sigma = 0.30)
) %>% filter(method %in% c("depth_fix", "fishash"))

per_regime <- bind %>%
  group_by(regime, d_sigma, method) %>%
  summarise(jaccard = mean(jaccard, na.rm = TRUE),
            FDR     = mean(fdr_pooled, na.rm = TRUE),
            recall  = mean(recall, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = c(jaccard, FDR, recall)) %>%
  mutate(gap_jac = jaccard_depth_fix - jaccard_fishash) %>%
  arrange(regime, d_sigma)
write.csv(per_regime, file.path(OUT, "dsigma_sweep.csv"), row.names = FALSE)

overall <- bind %>%
  group_by(d_sigma, method) %>%
  summarise(jaccard = round(mean(jaccard, na.rm = TRUE), 3),
            FDR     = round(mean(fdr_pooled, na.rm = TRUE), 3),
            recall  = round(mean(recall, na.rm = TRUE), 3), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = c(jaccard, FDR, recall)) %>%
  mutate(gap_jaccard = jaccard_depth_fix - jaccard_fishash)
cat("=== depth_fix vs fishash on Model B regimes, d_sigma=0.05 vs 0.30 ===\n")
print(as.data.frame(overall), row.names = FALSE)

cat("\n=== per-regime delta (jaccard_gap = depth_fix - fishash) ===\n")
delta <- per_regime %>% select(regime, d_sigma, gap_jac) %>%
  pivot_wider(names_from = d_sigma, values_from = gap_jac, names_prefix = "d_sigma_") %>%
  mutate(shift = d_sigma_0.3 - d_sigma_0.05)
print(as.data.frame(delta), row.names = FALSE)

# headline plot: depth_fix-vs-fishash Jaccard gap by regime at the two d_sigmas
pd <- per_regime %>% select(regime, d_sigma, gap_jac)
pd$d_sigma <- factor(paste0("d_sigma=", pd$d_sigma))
p <- ggplot(pd, aes(reorder(regime, gap_jac), gap_jac, fill = d_sigma)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c("d_sigma=0.05" = "#185FA5", "d_sigma=0.3" = "#D98E04")) +
  coord_flip() +
  labs(title = "depth_fix - fishash Jaccard gap by regime: d_sigma=0.05 vs d_sigma=0.3",
       subtitle = "if depth_fix's edge is real (not a model-fairness artifact), the orange bars should not collapse below blue",
       y = "depth_fix - fishash (Jaccard)", x = NULL, fill = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(OUT, "fig_dsigma_sweep.png"), p, width = 10, height = 5.5, dpi = 140)
cat("\nwrote dsigma_sweep.csv + fig_dsigma_sweep.png\n")
