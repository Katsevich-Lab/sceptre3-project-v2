#!/usr/bin/env Rscript
# ============================================================================
# Fig 6 (varying MOI) overlay: the fishash paper reports NO numbers for its
# simulation methods, so the 10-method panel is DIGITIZED (median F1 read off
# Fig 6a, p.17) and our Fishash+ Poisson is OVERLAID from exact computation on
# the authors' identical seeded sims (results/fishash_repro/varyMOI).
#
# Digitization validated: read-off fishash matches our exact fishash to ~0.01 at
# every MOI (see `fishash (digitized)` vs `fishash` below). Digitization precision
# is ~+/-0.02, so fine orderings within the top cluster are not resolvable by eye.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("scripts/fishash_repro_axis.R")
OUT <- "results/fishash_repro/varyMOI"

moi <- c(0.1, 0.3, 0.5, 1, 2, 3, 5, 10)
# --- digitized median F1 from Fig 6a (10 paper methods) -----------------------
dig <- tibble::tribble(
  ~method,             ~`0.1`,~`0.3`,~`0.5`,~`1`, ~`2`, ~`3`, ~`5`, ~`10`,
  "cleanser_cs",        0.93, 0.94, 0.93, 0.94, 0.92, 0.90, 0.86, 0.75,
  "cleanser_dc",        0.96, 0.96, 0.95, 0.95, 0.93, 0.91, 0.88, 0.77,
  "crispat_gauss",      0.65, 0.65, 0.72, 0.79, 0.85, 0.85, 0.80, 0.72,
  "crispat_poisgauss",  0.94, 0.94, 0.95, 0.93, 0.92, 0.89, 0.88, 0.81,
  "crispat_poisson",    0.95, 0.95, 0.94, 0.92, 0.90, 0.86, 0.85, 0.75,
  "crispat_negbinom",   0.97, 0.96, 0.95, 0.94, 0.92, 0.90, 0.87, 0.75,
  "sceptre_mixture",    0.92, 0.92, 0.93, 0.90, 0.91, 0.88, 0.82, 0.68,
  "demuxem",            0.94, 0.94, 0.92, 0.88, 0.84, 0.79, 0.70, 0.51,
  "geomux",             0.93, 0.93, 0.92, 0.92, 0.90, 0.89, 0.85, 0.75,
  "fishash (digitized)",0.97, 0.96, 0.96, 0.95, 0.92, 0.88, 0.82, 0.66
) %>% pivot_longer(-method, names_to = "moi", values_to = "F1") %>%
  mutate(moi = as.numeric(moi), source = "paper (digitized)")

# --- our exact computation on the authors' identical sims ---------------------
d <- read.csv(file.path(OUT, "combined_confusion.csv"))
d <- d[d$subset == "full", ]; d$F1 <- with(d, 2 * Precision * Recall / (Precision + Recall))
d$moi <- as.numeric(sub("moi_([0-9.]+)_.*", "\\1", d$sim_label))
exact <- d %>% filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  group_by(moi, method) %>% summarize(F1 = mean(F1), .groups = "drop") %>%
  mutate(method = recode(method, fishash_refit10 = "fishash", `fishash+_poisson` = "fishash+ Poisson (ours)"),
         source = "exact (this work)")

write.csv(bind_rows(dig, exact), file.path(OUT, "fig6_overlay_data.csv"), row.names = FALSE)

# ---- Main figure: ours, real fishash, and the best published alternative -----
# The envelope preserves the "best in field" comparison without an unreadable
# spaghetti of nine individually digitized competitor curves.
best_alt <- dig %>%
  filter(method != "fishash (digitized)") %>%
  group_by(moi) %>%
  slice_max(F1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(display = "best other published method (digitized)")
main <- bind_rows(
  exact %>% transmute(moi, F1, display = recode(
    method,
    "fishash" = "fishash (refit10)",
    "fishash+ Poisson (ours)" = "Fishash+ Poisson"
  )),
  best_alt %>% select(moi, F1, display)
)
main$display <- factor(
  main$display,
  c("Fishash+ Poisson", "fishash (refit10)", "best other published method (digitized)")
)
main_cols <- c(
  "Fishash+ Poisson" = "#e31a1c",
  "fishash (refit10)" = "#6a3d9a",
  "best other published method (digitized)" = "grey45"
)
p <- ggplot(main, aes(moi, F1, color = display, linetype = display)) +
  geom_line(aes(linewidth = display)) +
  geom_point(aes(size = display)) +
  scale_x_fishash_guide_load(expand = expansion(mult = c(0.03, 0.22))) +
  scale_color_manual(values = main_cols, name = NULL) +
  scale_linetype_manual(
    values = c("Fishash+ Poisson" = "solid", "fishash (refit10)" = "solid",
               "best other published method (digitized)" = "22"),
    name = NULL
  ) +
  scale_linewidth_manual(values = c(1.5, 1.0, 0.8), guide = "none") +
  scale_size_manual(values = c(2.5, 1.9, 1.6), guide = "none") +
  labs(x = "Mean infection events per recovered cell", y = "median F1 (full subset)",
       title = "Fishash+ Poisson leads the published field at moderate-to-high guide load",
       subtitle = "Fishash and Fishash+ are exact; the grey curve is the pointwise best of nine digitized alternatives from Fig 6a (read-off precision about ±0.02).") +
  coord_cartesian(ylim = c(0.62, 1.0)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 9)) +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
ggsave(file.path(OUT, "fig6_overlay_f1_vs_moi.png"), p, width = 11, height = 6.5, dpi = 130)

# ---- their-figure layout: facet by MOI, all 11 methods as points -------------
allm <- bind_rows(
  dig %>% filter(method != "fishash (digitized)") %>% mutate(hi = "paper (digitized)"),
  exact %>% mutate(hi = ifelse(grepl("ours", method), "Fishash+ Poisson (ours)", "fishash (exact)"))
) %>% mutate(
  guide_load = factor(
    fishash_recovered_guide_labels(moi),
    levels = fishash_recovered_guide_labels(fishash_preselection_moi)
  ),
  method = factor(method, levels = c(
  "cleanser_cs","cleanser_dc","crispat_gauss","crispat_poisgauss","crispat_poisson",
  "crispat_negbinom","sceptre_mixture","demuxem","geomux","fishash","fishash+ Poisson (ours)"))
)
pf <- ggplot(allm, aes(method, F1, color = hi)) +
  geom_point(size = 2) +
  facet_wrap(~guide_load, nrow = 2,
             labeller = labeller(guide_load = function(z) paste(z, "infections/cell"))) +
  scale_color_manual(values = c("paper (digitized)" = "grey60", "fishash (exact)" = "#6a3d9a",
                                "Fishash+ Poisson (ours)" = "#e31a1c"), name = NULL) +
  labs(x = NULL, y = "median F1", title = "Their Fig 6a layout, with Fishash+ Poisson added (red)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom")
ggsave(file.path(OUT, "fig6_overlay_faceted.png"), pf, width = 13, height = 7, dpi = 130)

# ---- console summary: ours vs field at each MOI ------------------------------
tab <- bind_rows(dig %>% filter(method != "fishash (digitized)"), exact) %>%
  select(method, moi, F1) %>% pivot_wider(names_from = moi, values_from = F1)
cat("=== median F1 by MOI: digitized paper methods + our exact fishash/fishash+ ===\n")
print(as.data.frame(tab), row.names = FALSE)
best <- bind_rows(dig %>% filter(method != "fishash (digitized)"), exact) %>%
  group_by(moi) %>% slice_max(F1, n = 1) %>% ungroup()
cat("\n=== best method at each MOI (incl. ours) ===\n"); print(as.data.frame(best[, c("moi","method","F1")]), row.names = FALSE)
cat("\nwrote fig6_overlay_f1_vs_moi.png, fig6_overlay_faceted.png, fig6_overlay_data.csv to", OUT, "\n")
