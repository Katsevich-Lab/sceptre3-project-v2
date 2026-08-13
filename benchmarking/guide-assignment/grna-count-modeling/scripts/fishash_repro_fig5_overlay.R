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

# ---- figure: decision-relevant comparison only -------------------------------
rlab <- c(high_gRNA = "high-gRNA regime (100 UMI, SNR 4)", low_gRNA = "low-gRNA regime (20 UMI, SNR 1)")
p <- ggplot(exact, aes(ng, F1, color = method, group = method)) +
  geom_line(aes(linewidth = method)) +
  geom_point(aes(size = method)) +
  facet_wrap(~regime, nrow = 1, labeller = labeller(regime = rlab)) +
  scale_x_log10(breaks = ng) +
  scale_color_manual(values = c("fishash" = "#6a3d9a", "fishash+ Poisson (ours)" = "#e31a1c"),
                     labels = c("fishash (refit10)", "Fishash+ Poisson"), name = NULL) +
  scale_linewidth_manual(values = c("fishash" = 1.0, "fishash+ Poisson (ours)" = 1.5),
                         guide = "none") +
  scale_size_manual(values = c("fishash" = 1.8, "fishash+ Poisson (ours)" = 2.4),
                    guide = "none") +
  labs(x = "number of guides", y = "median F1 (full subset)",
       title = "Fishash+ Poisson stays close to fishash as the guide library grows",
       subtitle = "Both regimes average 1.04 infection events per recovered cell (pre-selection MOI 0.3), leaving little co-occurrence to recover.") +
  coord_cartesian(ylim = c(0.68, 1.0)) +
  theme_bw(base_size = 12) +
  theme(plot.subtitle = element_text(size = 8.5), legend.position = "bottom")
ggsave(file.path(OUT, "fig5_overlay_f1_vs_nguides.png"), p, width = 12, height = 5.5, dpi = 130)
write.csv(bind_rows(dig %>% rename(), exact), file.path(OUT, "fig5_overlay_data.csv"), row.names = FALSE)
cat("\nwrote fig5_overlay_f1_vs_nguides.png + fig5_overlay_data.csv to", OUT, "\n")
