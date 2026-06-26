#!/usr/bin/env Rscript
# Which simulation setting has mode-separation like REAL Gasperini / Replogle?
#
# For each guide we measure the separation between the background and perturbed
# count modes on the log(1+count) scale, using a statistic computable IDENTICALLY
# on real and simulated data with NO ground truth: the smoothed-valley mode
# detector returns the two dominant modes (m1 < m2); we record
#   - bimodal?  (is a clean second mode detectable at all)
#   - gap = log1p(m2) - log1p(m1)         (separation magnitude)
#   - depth = 1 - valley_height/mode_height (how clean the trough is)
#   - m2     (perturbed-mode count level, to flag absurd scales)
# For the sims we ALSO use the ground truth to report the "true overlap" (the
# fraction of truly-perturbed cells that fall below the detected cut = the
# irreducible misses) as a validation that the truth-free statistic tracks real
# difficulty.

suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
HERE  <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
DATA  <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(HERE, "..", "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
extract_part <- function(x, k) vapply(strsplit(x, "_", fixed = TRUE),
  function(z) if (length(z) >= k) z[[k]] else NA_character_, character(1))

# per-guide separation row from the smoothed-valley detector (no truth used)
sep_row <- function(counts) {
  v <- smoothed_valley_threshold(counts)
  bimodal <- isTRUE(v$ok)
  data.frame(
    bimodal    = bimodal,
    gap        = if (bimodal) log1p(v$mode2) - log1p(v$mode1) else NA_real_,
    depth      = if (bimodal) 1 - v$valley_h / v$mode_h else NA_real_,
    m2         = if (bimodal) v$mode2 else NA_real_,
    t          = v$t,
    n_nonzero  = sum(counts > 0)
  )
}

# config: real datasets (sampled, no truth) + sims (all guides, truth available)
cfg <- list(
  list(label = "Gasperini (real)", id = "gasperini",            grp = "real", sample = 400),
  list(label = "Replogle (real)",  id = "replogle-rd7_medium",  grp = "real", sample = 400),
  list(label = "old_sims",         id = "sims_sum_repeat_old",  grp = "sim",  by_p = TRUE),
  list(label = "2np_3p",           id = "sims_sum_2np_3p",      grp = "sim",  by_p = TRUE),
  list(label = "1np_3p_disp",      id = "sims_sum_1np_3p_disp", grp = "sim",  by_p = TRUE),
  list(label = "gasp_calib (NEW)", id = "sims_gasperini_calibrated", grp = "sim", by_p = TRUE)
)

set.seed(11)
rows <- list()
for (cc in cfg) {
  gm <- readRDS(file.path(DATA, cc$id, "sceptre", "grna_matrix.rds"))
  gn <- rownames(gm)
  sel <- if (!is.null(cc$sample) && cc$sample < length(gn)) sort(sample(seq_along(gn), cc$sample)) else seq_along(gn)
  sub <- as(gm[sel, , drop = FALSE], "RsparseMatrix")
  truth <- NULL
  tfp <- file.path(DATA, cc$id, "true_pert_matrix.rds")
  if (file.exists(tfp)) truth <- readRDS(tfp)

  for (r in seq_along(sel)) {
    g <- gn[sel[r]]
    counts <- as.numeric(sub[r, ])
    sr <- sep_row(counts)
    # group label: sims split by perturbation level (P part of the guide name)
    grp_label <- if (isTRUE(cc$by_p)) paste0(cc$label, ":", extract_part(g, 3)) else cc$label
    overlap <- NA_real_
    if (!is.null(truth)) {
      tr <- as.numeric(truth[g, ]) == 1 & counts > 0
      if (sum(tr) > 0 && is.finite(sr$t)) overlap <- mean(counts[tr] < sr$t)
    }
    rows[[length(rows) + 1]] <- cbind(
      data.frame(dataset = cc$label, group = grp_label, kind = cc$grp, overlap = overlap), sr)
  }
  cat("processed", cc$label, "(", length(sel), "guides )\n")
}
df <- bind_rows(rows)
saveRDS(df, file.path(HERE, "results", "realism_separation.rds"))

# ---- summary table ----------------------------------------------------------
summ <- df %>% group_by(group, kind) %>% summarize(
  n_guides       = n(),
  pct_bimodal    = round(100 * mean(bimodal), 1),
  median_gap     = round(median(gap, na.rm = TRUE), 2),
  median_depth   = round(median(depth, na.rm = TRUE), 2),
  median_m2_count= round(median(m2, na.rm = TRUE)),
  median_overlap = round(median(overlap, na.rm = TRUE), 3),
  .groups = "drop") %>% arrange(median_gap)
cat("\n===== per-guide separation by dataset =====\n")
print(as.data.frame(summ), row.names = FALSE)
write.csv(summ, file.path(HERE, "results", "realism_separation_summary.csv"), row.names = FALSE)

# ---- plot: per-guide mode-gap distribution, real highlighted ----------------
pd <- df %>% filter(bimodal) %>%
  mutate(group = factor(group, levels = summ$group))   # ordered by median gap
real_med <- df %>% filter(kind == "real", bimodal) %>% group_by(dataset) %>%
  summarize(m = median(gap, na.rm = TRUE), .groups = "drop")
p <- ggplot(pd, aes(group, gap, fill = kind)) +
  geom_hline(data = real_med, aes(yintercept = m, color = dataset), linetype = 2, linewidth = 0.6) +
  geom_boxplot(outlier.size = 0.4, width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c(real = "#d1495b", sim = "grey75")) +
  scale_color_manual(values = c(`Gasperini (real)` = "#d1495b", `Replogle (real)` = "#1f6f8b")) +
  labs(title = "Mode separation: real vs. simulated (per guide)",
       subtitle = "gap = log1p(perturbed mode) - log1p(background mode); dashed = real medians",
       x = NULL, y = "mode gap on log(1+count) scale", fill = NULL, color = "real reference") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(HERE, "results", "realism_separation.png"), p, width = 11, height = 5.5, dpi = 120)
cat("\nWrote results/realism_separation.png\n")
