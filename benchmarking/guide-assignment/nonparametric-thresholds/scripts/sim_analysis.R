#!/usr/bin/env Rscript
# Aggregate scores_all.csv into the document's tables + figures.
# Model A = barnyard semi-synthetic (ground truth; chemistry x purity x mu x MOI).
# Model B = data-derived regimes (one per real non-barnyard dataset).
# Counts are computed from the data so figure titles stay in sync after rebuilds.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
HERE <- getwd(); OUT <- file.path(HERE, "results", "sim_framework")
s   <- read.csv(file.path(OUT, "scores_all.csv"))
reg <- read.csv(file.path(OUT, "regimes.csv"))           # regime -> separation, mu, moi, ambient

ord <- c("thresh3","thresh10","thresh20","otsu","valley","ambient","geomux","fishash","depth_fix","sceptre","crispat")
ord <- ord[ord %in% s$method]; s$method <- factor(s$method, levels = ord)
fam <- c(thresh3="fixed threshold", thresh10="fixed threshold", thresh20="fixed threshold",
         otsu="per-guide threshold", valley="per-guide threshold",
         ambient="ambient/contingency", geomux="ambient/contingency", fishash="ambient/contingency",
         depth_fix="ambient/contingency", sceptre="mixture", crispat="mixture")
s$family <- fam[as.character(s$method)]
pal <- c("fixed threshold"="#9A8C98","per-guide threshold"="#C9ADA7",
         "ambient/contingency"="#1D9E75","mixture"="#185FA5")

n_datasets <- length(unique(s$id))
n_regimes  <- length(unique(s$regime[s$model == "B"]))

# ---- overall ranking --------------------------------------------------------
overall <- s %>% group_by(method) %>%
  summarise(jaccard=round(mean(jaccard,na.rm=T),3), precision=round(mean(precision,na.rm=T),3),
            recall=round(mean(recall,na.rm=T),3), FDR=round(mean(fdr_pooled,na.rm=T),3), .groups="drop") %>%
  arrange(-jaccard)
write.csv(overall, file.path(OUT, "method_overall.csv"), row.names = FALSE)
cat(sprintf("=== overall method ranking (mean over all %d datasets) ===\n", n_datasets))
print(as.data.frame(overall), row.names = FALSE)

# =============================================================================
# MODEL B: the data-derived regimes, ordered easy -> hard by mode separation
# =============================================================================
B <- s %>% filter(model == "B") %>% left_join(reg, by = "regime")
reg_order <- reg %>% arrange(desc(separation)) %>% pull(regime)
B$regime <- factor(B$regime, levels = rev(reg_order))                 # hard -> easy on x
Bagg <- B %>% group_by(regime, separation, family, method) %>%
  summarise(jaccard = mean(jaccard, na.rm = TRUE), FDR = mean(fdr_pooled, na.rm = TRUE), .groups = "drop")

sel <- c("thresh3","thresh10","ambient","geomux","fishash","depth_fix","sceptre","crispat")
mcol <- c(thresh3="#C77",thresh10="#955",ambient="#3a3",geomux="#7c4",fishash="#1D9E75",
          depth_fix="#0a6",sceptre="#185FA5",crispat="#7aa")
sweep <- Bagg %>% filter(method %in% sel) %>%
  pivot_longer(c(jaccard, FDR), names_to = "metric", values_to = "value")
sweep$metric <- factor(sweep$metric, levels = c("jaccard", "FDR"))
f_sweep <- ggplot(sweep, aes(regime, value, colour = method, group = method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
  geom_hline(data = data.frame(metric = factor("FDR", levels = c("jaccard", "FDR"))),
             aes(yintercept = 0.05), linetype = "dashed", colour = "grey40") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") + scale_colour_manual(values = mcol) +
  labs(title = sprintf("Method performance across %d data-derived regimes (one per real dataset)", n_regimes),
       subtitle = "regimes ordered hard -> easy by mode separation",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1),
                                  legend.position = "top", legend.title = element_blank())
ggsave(file.path(OUT, "fig_regime_sweep.png"), f_sweep, width = 11, height = 7, dpi = 140)

# heatmap: method x regime, fill = Jaccard
hm <- B %>% group_by(regime, method) %>% summarise(jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop")
f_hm <- ggplot(hm, aes(regime, method, fill = jaccard)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.3) +
  scale_fill_viridis_c(option = "D", limits = c(0, 1)) +
  labs(title = sprintf("Jaccard by method x regime (%d real-data-derived regimes)", n_regimes),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(OUT, "fig_regime_heatmap.png"), f_hm, width = 12, height = 5, dpi = 140)

# regime-level table for the doc.  The left_join above suffixes shared columns
# from `reg` with .y; we use the regimes.csv versions explicitly for clarity.
regtab <- B %>% group_by(regime, separation, mu_pert, moi, target_pct3, sim_pct3) %>%
  summarise(best = ord[which.max(tapply(jaccard, method, mean, na.rm = TRUE))],
            depth_fix = round(mean(jaccard[method == "depth_fix"], na.rm = TRUE), 2),
            thresh3   = round(mean(jaccard[method == "thresh3"],   na.rm = TRUE), 2),
            .groups = "drop") %>% arrange(desc(separation))
write.csv(regtab, file.path(OUT, "regime_method_table.csv"), row.names = FALSE)

# =============================================================================
# MODEL A: barnyard semi-synthetic (FDR control, purity, high-MOI)
# =============================================================================
A <- s %>% filter(model == "A")
A$chemistry <- factor(A$chemistry, levels = c("CROP-seq", "direct-capture"))
A$moi_level <- factor(A$moi_level, levels = c("low", "high", "vhigh"))

# FDR control by chemistry
af <- A %>% group_by(chemistry, family, method) %>%
  summarise(FDR = mean(fdr_pooled, na.rm = TRUE), .groups = "drop")
f_fdr <- ggplot(af, aes(method, FDR, fill = family)) + geom_col() +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  facet_wrap(~chemistry) + scale_fill_manual(values = pal) +
  labs(title = "Model A (real barnyard ambient): realized FDR (dashed = nominal 0.05)",
       y = "mean FDR", x = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   legend.position = "top", legend.title = element_blank())
ggsave(file.path(OUT, "fig_fdr.png"), f_fdr, width = 10, height = 4.6, dpi = 140)

# purity sensitivity (direct-capture)
ap <- A %>% filter(chemistry == "direct-capture") %>% group_by(method, purity) %>%
  summarise(FDR = mean(fdr_pooled, na.rm = TRUE), .groups = "drop") %>%
  mutate(purity = factor(paste0("purity ", purity)))
f_pur <- ggplot(ap, aes(method, FDR, fill = purity)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("purity 0.9" = "#D98E04", "purity 0.99" = "#185FA5")) +
  labs(title = "Model A purity sensitivity (direct-capture): FDR at 0.90 vs 0.99",
       subtitle = "the doublet false-positive floor at 0.90 collapses at 0.99 for calibrated methods",
       y = "mean FDR", x = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   legend.position = "top", legend.title = element_blank())
ggsave(file.path(OUT, "fig_purity.png"), f_pur, width = 9, height = 4.5, dpi = 140)

# depth_fix vs fishash by MOI (the high-MOI masking story; Model A only)
moicmp <- s %>% filter(method %in% c("depth_fix", "fishash"), model == "A") %>%
  group_by(moi_level = factor(moi_level, levels = c("low", "high", "vhigh")), method) %>%
  summarise(jaccard = round(mean(jaccard, na.rm = TRUE), 3),
            recall  = round(mean(recall,  na.rm = TRUE), 3), .groups = "drop")
write.csv(moicmp, file.path(OUT, "depthfix_vs_fishash_moi.csv"), row.names = FALSE)

cat(sprintf("\nwrote method_overall.csv, regime_method_table.csv, depthfix_vs_fishash_moi.csv (%d regimes / %d datasets)\n",
            n_regimes, n_datasets))
cat("figures: fig_regime_sweep.png, fig_regime_heatmap.png, fig_fdr.png, fig_purity.png\n")
