#!/usr/bin/env Rscript
# Place the 5 recent (2024-2026) survey datasets on the same mode-separation axis
# as the real Gasperini/Replogle data and the simulations. Truth-free statistic:
# per-guide gap = log1p(perturbed mode) - log1p(background mode) from the smoothed
# valley detector. Answers: do the (Replogle-calibrated) sims represent real data?

suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
SURV <- path.expand("~/data/external/perturbseq-survey")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

sep_of <- function(counts) { v <- smoothed_valley_threshold(counts)
  if (isTRUE(v$ok)) c(gap = log1p(v$mode2) - log1p(v$mode1), m2 = v$mode2, bim = 1)
  else c(gap = NA_real_, m2 = NA_real_, bim = 0) }

surv <- c(`DC-TAP K562 MOI=1 (low)` = "dctap_k562_lowmoi_GSE303901",
          `DC-TAP K562 MOI=14 (high)` = "dctap_k562_highmoi_GSE303901",
          `EndoC T2D (low)` = "endoc_t2d_perturbseq_GSE273677",
          `Erythroid multiome (high)` = "perturb_multiome_erythroid_GSE274113",
          `CD4 genome-scale (low)` = "tcell_cd4_perturbseq_GSE314342")
set.seed(3); rows <- list()
for (sn in names(surv)) {
  gm <- readRDS(file.path(SURV, surv[[sn]], "grna_matrix.rds"))
  sel <- if (nrow(gm) > 400) sort(sample(nrow(gm), 400)) else seq_len(nrow(gm))
  sub <- as(gm[sel, , drop = FALSE], "RsparseMatrix")
  s <- t(vapply(seq_along(sel), function(r) sep_of(as.numeric(sub[r, ])), numeric(3)))
  rows[[sn]] <- data.frame(group = sn, kind = "real (2024-26)", gap = s[, "gap"], bim = s[, "bim"])
  cat("done", sn, "\n")
}
surv_df <- bind_rows(rows)

# combine with the earlier real+sim separation results if present
prev_fp <- file.path(HERE, "results", "realism_separation.rds")
combined <- surv_df
if (file.exists(prev_fp)) {
  prev <- readRDS(prev_fp) |>
    mutate(kind = ifelse(kind == "real", "real (Gasp/Repl)", "simulation")) |>
    select(group, kind, gap, bim = bimodal) |> mutate(bim = as.numeric(bim))
  combined <- bind_rows(prev, surv_df)
}
summ <- combined %>% group_by(group, kind) %>%
  summarize(pct_bimodal = round(100*mean(bim, na.rm=TRUE)), median_gap = round(median(gap, na.rm=TRUE),2), .groups="drop") %>%
  arrange(median_gap)
cat("\n===== mode separation across real (incl. 2024-26) + simulated =====\n")
print(as.data.frame(summ), row.names = FALSE)
write.csv(summ, file.path(HERE, "results", "survey_separation_summary.csv"), row.names = FALSE)

pd <- combined %>% filter(!is.na(gap)) %>% mutate(group = factor(group, levels = summ$group))
p <- ggplot(pd, aes(group, gap, fill = kind)) +
  geom_boxplot(outlier.size = .3, width = .7) + coord_flip() +
  scale_fill_manual(values = c(`real (2024-26)`="#d7301f", `real (Gasp/Repl)`="#fdae61", simulation="grey80")) +
  labs(title = "Mode separation across real datasets (incl. 5 recent) and simulations",
       subtitle = "gap = log1p(signal mode) - log1p(background mode); higher = more separable",
       x = NULL, y = "per-guide mode gap (log1p scale)", fill = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(HERE, "results", "Survey_SeparationAxis.png"), p, width = 10, height = 6, dpi = 120)
cat("\nWrote results/Survey_SeparationAxis.png\n")
