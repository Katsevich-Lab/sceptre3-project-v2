#!/usr/bin/env Rscript
# ============================================================================
# Fig 5 (varying number of guides, both regimes) overlay. 10 paper methods
# DIGITIZED from Fig 5a (p.15); Fishash+ Poisson & fishash computed EXACTLY on
# the authors' identical seeded sims (varyNguides = high-gRNA, varyNguidesLow =
# low-gRNA). Both Fig-5 regimes are MOI 0.3 (low MOI), so this is the no-self-
# masking regime where the depth fix costs a little precision rather than winning
# -- included for completeness / to show no regression at scale.
# The authors' 80,000-guide tier is their runtime/memory scenario; we overlay on
# {20,200,2000,20000} where we computed our method.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
OUT <- "results/fishash_repro/varyNguides"

ng <- c(20, 200, 2000, 20000)
# --- digitized median F1 from Fig 5a ------------------------------------------
mk <- function(v) setNames(v, ng)
dig <- bind_rows(
  # HIGH-gRNA regime (top row of Fig 5a)
  tibble::tribble(
    ~method,             ~`20`,~`200`,~`2000`,~`20000`,
    "cleanser_cs",        0.84, 0.93, 0.95, 0.95,
    "cleanser_dc",        0.68, 0.96, 0.95, 0.95,
    "crispat_gauss",      0.72, 0.73, 0.90, 0.90,
    "crispat_poisgauss",  0.68, 0.94, 0.93, 0.93,
    "crispat_poisson",    0.91, 0.95, 0.95, 0.88,
    "crispat_negbinom",   0.90, 0.96, 0.96, 0.90,
    "sceptre_mixture",    0.82, 0.93, 0.93, 0.84,
    "demuxem",            0.87, 0.92, 0.93, 0.95,
    "geomux",             0.92, 0.90, 0.96, 0.96,
    "fishash (digitized)",0.94, 0.96, 0.97, 0.97
  ) %>% mutate(regime = "high_gRNA"),
  # LOW-gRNA regime (bottom row of Fig 5a) -- harder (SNR 1), reads are coarser
  tibble::tribble(
    ~method,             ~`20`,~`200`,~`2000`,~`20000`,
    "cleanser_cs",        0.68, 0.76, 0.75, 0.70,
    "cleanser_dc",        0.64, 0.81, 0.76, 0.66,
    "crispat_gauss",      0.57, 0.74, 0.79, 0.82,
    "crispat_poisgauss",  0.62, 0.79, 0.80, 0.70,
    "crispat_poisson",    0.80, 0.83, 0.84, 0.42,
    "crispat_negbinom",   0.74, 0.84, 0.68, 0.30,
    "sceptre_mixture",    0.57, 0.82, 0.83, 0.81,
    "demuxem",            0.77, 0.80, 0.83, 0.83,
    "geomux",             0.75, 0.83, 0.80, 0.82,
    "fishash (digitized)",0.68, 0.84, 0.87, 0.88
  ) %>% mutate(regime = "low_gRNA")
) %>% pivot_longer(-c(method, regime), names_to = "ng", values_to = "F1") %>%
  mutate(ng = as.numeric(ng))

# --- our exact computation ----------------------------------------------------
load_exact <- function(dir, regime) {
  d <- read.csv(file.path("results/fishash_repro", dir, "combined_confusion.csv"))
  d <- d[d$subset == "full", ]; d$F1 <- with(d, 2 * Precision * Recall / (Precision + Recall))
  d$ng <- as.integer(sub("nguides_([0-9]+)_.*", "\\1", d$sim_label))
  d %>% filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
    group_by(ng, method) %>% summarize(F1 = mean(F1), .groups = "drop") %>%
    mutate(regime = regime,
           method = recode(method, fishash_refit10 = "fishash", `fishash+_poisson` = "fishash+ Poisson (ours)"))
}
exact <- bind_rows(load_exact("varyNguides", "high_gRNA"), load_exact("varyNguidesLow", "low_gRNA"))

# validation: digitized fishash vs exact fishash
val <- dig %>% filter(method == "fishash (digitized)") %>%
  inner_join(exact %>% filter(method == "fishash") %>% select(ng, regime, F1exact = F1), by = c("ng", "regime"))
cat("=== digitization check: |digitized fishash - exact fishash| ===\n")
print(as.data.frame(val %>% mutate(absdiff = round(abs(F1 - F1exact), 3)) %>% select(regime, ng, F1, F1exact, absdiff)), row.names = FALSE)

# ---- figure: F1 vs nguides, facet by regime, ours overlaid -------------------
others <- dig %>% filter(method != "fishash (digitized)")
rlab <- c(high_gRNA = "high-gRNA regime (100 UMI, SNR 4)", low_gRNA = "low-gRNA regime (20 UMI, SNR 1)")
p <- ggplot() +
  geom_line(data = others, aes(ng, F1, group = method), color = "grey78", linewidth = 0.5) +
  geom_line(data = exact %>% filter(method == "fishash"), aes(ng, F1), color = "#6a3d9a", linewidth = 1.0) +
  geom_line(data = exact %>% filter(grepl("ours", method)), aes(ng, F1), color = "#e31a1c", linewidth = 1.6) +
  geom_point(data = exact %>% filter(grepl("ours", method)), aes(ng, F1), color = "#e31a1c", size = 2.2) +
  facet_wrap(~regime, nrow = 1, labeller = labeller(regime = rlab)) +
  scale_x_log10(breaks = ng) +
  labs(x = "number of guides", y = "median F1 (full subset)",
       title = "Fishash+ Poisson overlaid on the fishash paper's varying-guides panel (Fig 5a)",
       subtitle = "Both regimes are MOI 0.3 (low MOI): red (ours) tracks just below fishash at the top of the field — the small low-MOI precision cost, no self-masking to recover. grey = 8 other digitized paper methods.") +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  theme_bw(base_size = 12) + theme(plot.subtitle = element_text(size = 8.5))
ggsave(file.path(OUT, "fig5_overlay_f1_vs_nguides.png"), p, width = 12, height = 5.5, dpi = 130)
write.csv(bind_rows(dig %>% rename(), exact), file.path(OUT, "fig5_overlay_data.csv"), row.names = FALSE)
cat("\nwrote fig5_overlay_f1_vs_nguides.png + fig5_overlay_data.csv to", OUT, "\n")
