#!/usr/bin/env Rscript
# Figure backing the guide-margin contamination claim: per guide, how much its
# GLOBAL count is inflated by signal over its AMBIENT level (total_g / ambient_g),
# and the spread of that across guides (q90/q10) = guide_contam_spread.
#   Panel A: guide_contam_spread across all 17 datasets (from survey CSV).
#   Panel B: per-guide total vs ambient counts for a few datasets (the spread, shown).
suppressPackageStartupMessages({
  library(Matrix); library(sparseMatrixStats); library(ggplot2); library(dplyr); library(patchwork)
})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
GA   <- normalizePath(file.path(HERE, ".."))
SURV <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
DATA <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
REPRO<- file.path(HERE, "external", "repro_work")
D    <- file.path(HERE, "results", "ambient_intuition")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

# ---------- Panel A: the metric across all datasets (from the survey table) ----
t <- read.csv(file.path(D, "survey_depth_composition.csv"))
t$lab <- substr(gsub("_", " ", sub("_figshare.*$","", sub("_GSE[0-9]+$","", t$dataset))), 1, 24)
tA <- t[order(t$guide_contam_spread), ]; tA$lab <- factor(tA$lab, levels = tA$lab)
pA <- ggplot(tA, aes(guide_contam_spread, lab)) +
  geom_col(aes(fill = log10(moi)), width = 0.75) +
  geom_vline(xintercept = median(t$guide_contam_spread), linetype = "dashed", colour = "grey45") +
  scale_fill_viridis_c(option = "C", end = 0.9, name = "log10 MOI") +
  labs(title = "A. Does all-UMI share equal ambient share?  Spread across 17 datasets",
       subtitle = "q90/q10 of (all-UMI share / ambient share) across guides; 1 = identical for every guide; dashed = median (1.9x)",
       x = "spread of  share_g / gamma_g   across guides", y = NULL) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8), axis.text.y = element_text(size = 7.5),
    legend.position = c(0.85, 0.25), legend.key.size = unit(0.35, "cm"))

# ---------- Panel B: per-guide total vs ambient counts, a few datasets ----------
sel <- list(
  `ipsc (8.2x)`        = function() readRDS(file.path(SURV, "ipsc_crispri_hipsci_figshare27989294/grna_matrix.rds")),
  `replogle (4.0x)`    = function() readRDS(file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds")),
  `schraivogel (1.3x)` = function() { library(fishash); data(crispat_schraivogel); crispat_schraivogel },
  `barnyard cropseq (1.4x)` = function() as(readMM(file.path(REPRO, "mix0hr_Cropseq_grna_counts.mtx")), "CsparseMatrix"))

per_guide <- function(gm) {
  gm <- as(gm, "CsparseMatrix"); storage.mode(gm@x) <- "double"
  if (ncol(gm) > 25000) { set.seed(1); gm <- gm[, sort(sample(ncol(gm), 25000))] }
  gm <- gm[, colSums(gm) > 0, drop = FALSE]
  A <- suppressMessages(ambient_test_assign(gm, q = 0.05, model = "hypergeometric", n_iter = 1)$assignment_matrix) * 1
  tot <- rowSums(gm); sig <- rowSums(gm * A); amb <- tot - sig
  data.frame(total = tot, ambient = amb)[tot > 0, ]
}
rows <- list()
for (nm in names(sel)) {
  pg <- tryCatch(per_guide(sel[[nm]]()), error = function(e) { cat(nm, "ERR", conditionMessage(e), "\n"); NULL })
  if (!is.null(pg)) { pg$dataset <- nm; rows[[nm]] <- pg; cat(nm, ":", nrow(pg), "guides\n") }
}
pg <- do.call(rbind, rows)
# share = representation among ALL UMIs ; gamma = representation among AMBIENT UMIs
pg <- pg %>% group_by(dataset) %>%
  mutate(share = total / sum(total), gamma = ambient / sum(ambient)) %>%
  ungroup() %>% filter(ambient >= 10)            # well-observed guides only
pg$dataset <- factor(pg$dataset, levels = names(sel))
ann <- pg %>% group_by(dataset) %>% summarise(
  spread = quantile(share / gamma, .9) / quantile(share / gamma, .1),
  within2x = mean(abs(log2(share / gamma)) <= 1), .groups = "drop")

pB <- ggplot(pg, aes(gamma, share)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(alpha = 0.35, size = 0.8, colour = "#7F77DD") +
  geom_text(data = ann, aes(x = 10^(-0.5), y = Inf,
            label = sprintf("spread %.1fx\n%.0f%% within 2x", spread, 100 * within2x)),
            hjust = 1, vjust = 1.3, size = 2.7, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "B. Per guide: representation among ALL UMIs vs among AMBIENT UMIs",
       subtitle = "on the dashed line = a guide's all-UMI share equals its ambient share (naive composition is a good proxy)",
       x = "ambient representation   gamma_g = amb_g / (all ambient UMIs)",
       y = "all-UMI representation   share_g = N_g / N") +
  theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8), strip.text = element_text(size = 8.5))

ggsave(file.path(D, "guide_margin.png"), pA / pB + plot_layout(heights = c(1.1, 1)),
       width = 12, height = 9, dpi = 150)
cat("wrote", file.path(D, "guide_margin.png"), "\n")
