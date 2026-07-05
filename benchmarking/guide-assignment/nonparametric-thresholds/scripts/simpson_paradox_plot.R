#!/usr/bin/env Rscript
# Plot the Simpson's-paradox reproduction (reads results/simpson_paradox_repro.csv,
# produced by simpson_paradox_repro.R). One panel for realized FDR, one for recall;
# uncorrected (plain Fisher, refit=0) vs corrected (Simpson correction, refit=10).
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
RES  <- file.path(HERE, "results")

df <- read.csv(file.path(RES, "simpson_paradox_repro.csv"))
df$regime <- factor(df$regime, levels = c("low", "medium", "high"),
  labels = c("low\n(noise ⊥ signal)", "medium", "high\n(noise ∝ signal)"))

long <- df %>%
  pivot_longer(c(unc_fdr, cor_fdr, unc_rec, cor_rec),
               names_to = "k", values_to = "value") %>%
  mutate(method = ifelse(grepl("^unc", k), "uncorrected (refit=0)", "corrected (refit=10)"),
         metric = ifelse(grepl("fdr$", k), "realized FDR", "recall"),
         method = factor(method, levels = c("uncorrected (refit=0)", "corrected (refit=10)")))

hline <- data.frame(metric = "realized FDR", yint = 0.05)

p <- ggplot(long, aes(regime, value, colour = method)) +
  geom_hline(data = hline, aes(yintercept = yint), linetype = "dashed", colour = "grey40") +
  stat_summary(fun = mean, geom = "line", aes(group = method), linewidth = 0.8,
               position = position_dodge(0.4)) +
  stat_summary(fun = mean, geom = "point", size = 2.6, position = position_dodge(0.4)) +
  geom_point(alpha = 0.25, size = 1, position = position_dodge(0.4)) +
  facet_wrap(~metric, scales = "free_y") +
  scale_colour_manual(values = c("uncorrected (refit=0)" = "#d7301f",
                                 "corrected (refit=10)"   = "#238b45")) +
  labs(title = "Reproducing fishash's Simpson's-paradox phenomenon (paper Fig. 7)",
       subtitle = paste("authors' simulator (20 guides, 4000 cells, SNR 4, MOI 0.3);",
                        "real fishash package; dashed = nominal 5% FDR"),
       x = "signal–noise composition correlation", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 8.5))

ggsave(file.path(RES, "simpson_paradox_repro.png"), p, width = 8, height = 4, dpi = 150)
cat("wrote", file.path(RES, "simpson_paradox_repro.png"), "\n")
