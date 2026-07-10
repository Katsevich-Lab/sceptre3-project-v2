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

# ---- Figure: median F1 vs MOI, our method overlaid on the digitized field ----
others <- dig %>% filter(method != "fishash (digitized)")
key    <- dig %>% filter(method %in% c("cleanser_dc", "crispat_poisgauss"))   # paper's high-MOI winners
p <- ggplot() +
  geom_line(data = others, aes(moi, F1, group = method), color = "grey75", linewidth = 0.5) +
  geom_line(data = key, aes(moi, F1, color = method), linewidth = 0.7, linetype = "22") +
  geom_point(data = key, aes(moi, F1, color = method), size = 1.5) +
  geom_line(data = exact %>% filter(method == "fishash"),
            aes(moi, F1), color = "#6a3d9a", linewidth = 1.0) +
  geom_line(data = exact %>% filter(grepl("ours", method)),
            aes(moi, F1), color = "#e31a1c", linewidth = 1.6) +
  geom_point(data = exact %>% filter(grepl("ours", method)),
             aes(moi, F1), color = "#e31a1c", size = 2.4) +
  scale_x_log10(breaks = moi, expand = expansion(mult = c(0.03, 0.22))) +
  scale_color_manual(values = c(cleanser_dc = "#1f78b4", crispat_poisgauss = "#33a02c"), name = "paper's high-MOI winners") +
  annotate("text", x = 11, y = 0.769, label = "fishash+ Poisson\n(ours)", color = "#e31a1c", hjust = 0, size = 3.8, fontface = 2, lineheight = 0.9) +
  annotate("text", x = 11, y = 0.66, label = "fishash", color = "#6a3d9a", hjust = 0, size = 3.8) +
  annotate("text", x = 0.1, y = 0.50, label = "grey = other 8 paper methods (digitized)", color = "grey55", hjust = 0, size = 3.3) +
  labs(x = "MOI", y = "median F1 (full subset)",
       title = "Fishash+ Poisson overlaid on the fishash paper's varying-MOI panel (Fig 6a)",
       subtitle = "10 paper methods digitized from Fig 6a; fishash & Fishash+ Poisson computed exactly on the authors' identical seeded sims") +
  coord_cartesian(ylim = c(0.45, 1.0)) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))
ggsave(file.path(OUT, "fig6_overlay_f1_vs_moi.png"), p, width = 11, height = 6.5, dpi = 130)

# ---- their-figure layout: facet by MOI, all 11 methods as points -------------
allm <- bind_rows(
  dig %>% filter(method != "fishash (digitized)") %>% mutate(hi = "paper (digitized)"),
  exact %>% mutate(hi = ifelse(grepl("ours", method), "Fishash+ Poisson (ours)", "fishash (exact)"))
) %>% mutate(method = factor(method, levels = c(
  "cleanser_cs","cleanser_dc","crispat_gauss","crispat_poisgauss","crispat_poisson",
  "crispat_negbinom","sceptre_mixture","demuxem","geomux","fishash","fishash+ Poisson (ours)")))
pf <- ggplot(allm, aes(method, F1, color = hi)) +
  geom_point(size = 2) +
  facet_wrap(~moi, nrow = 2, labeller = label_both) +
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
