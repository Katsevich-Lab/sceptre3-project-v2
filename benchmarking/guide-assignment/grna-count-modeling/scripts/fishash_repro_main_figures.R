#!/usr/bin/env Rscript
# Consolidated main-text figures for fishash_repro.qmd.
#
# Figure 1 puts every decision-relevant accuracy comparator in one panel and
# pairs it with the sceptre runtime comparison. Figure 2 compares F1 levels at
# two robustness boundaries on common axes. A detailed decomposition shows
# precision, recall, and F1 under NB ambient.

# Run from benchmarking/guide-assignment/grna-count-modeling.
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})
source("scripts/fishash_repro_axis.R")

OUT <- "results/fishash_repro"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

f1_score <- function(precision, recall) {
  ifelse(precision + recall == 0, NA_real_, 2 * precision * recall / (precision + recall))
}

read_method_panel <- function(scenario) {
  d <- read.csv(file.path(OUT, scenario, "combined_confusion.csv"))
  d %>%
    filter(subset == "full") %>%
    mutate(
      scenario = scenario,
      moi = as.numeric(sub("moi_([0-9.]+)_.*", "\\1", sim_label)),
      iter = as.integer(sub(".*_iter_([0-9]+)$", "\\1", sim_label)),
      F1 = f1_score(Precision, Recall)
    )
}

read_guide_panel <- function(scenario) {
  d <- read.csv(file.path(OUT, scenario, "combined_confusion.csv"))
  d %>%
    filter(subset == "full") %>%
    mutate(
      scenario = scenario,
      ng = as.integer(sub("nguides_([0-9]+)_.*", "\\1", sim_label)),
      iter = as.integer(sub(".*_iter_([0-9]+)$", "\\1", sim_label)),
      F1 = f1_score(Precision, Recall)
    )
}

report_theme <- function(base_size = 19) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.5, "lines"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.9), color = "grey30"),
      strip.background = element_rect(fill = "grey95", color = "grey75"),
      strip.text = element_text(face = "bold")
    )
}

method_levels <- c(
  "fishash+",
  "fishash (refit10)",
  "sceptre mixture",
  "best remaining published method (digitized)"
)
method_cols <- c(
  "fishash+" = "#D62728",
  "fishash (refit10)" = "#6A3D9A",
  "sceptre mixture" = "#1F78B4",
  "best remaining published method (digitized)" = "grey45"
)
method_types <- c(
  "fishash+" = "solid",
  "fishash (refit10)" = "solid",
  "sceptre mixture" = "solid",
  "best remaining published method (digitized)" = "22"
)
method_labels <- c(
  "fishash+" = "fishash+",
  "fishash (refit10)" = "fishash",
  "sceptre mixture" = "sceptre",
  "best remaining published method (digitized)" = "best other (digitized)"
)

# ---- Figure 1A: one accuracy comparison on the 10 common seeded replicates --
vary_moi <- read_method_panel("varyMOI")

accuracy_exact <- vary_moi %>%
  filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  mutate(display = recode(
    method,
    fishash_refit10 = "fishash (refit10)",
    `fishash+_poisson` = "fishash+"
  )) %>%
  select(moi, iter, display, F1)

sceptre_exact <- read.csv(file.path(OUT, "sceptre_h2h", "h2h_accuracy_reps.csv")) %>%
  filter(iter <= 10, method == "sceptre_mixture") %>%
  transmute(moi, iter, display = "sceptre mixture", F1)

accuracy_summary <- bind_rows(accuracy_exact, sceptre_exact) %>%
  group_by(moi, display) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  )

digitized <- read.csv(file.path(OUT, "varyMOI", "fig6_overlay_data.csv")) %>%
  filter(
    source == "paper (digitized)",
    !method %in% c("fishash (digitized)", "sceptre_mixture")
  )

best_remaining <- digitized %>%
  group_by(moi) %>%
  slice_max(F1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    moi,
    display = "best remaining published method (digitized)",
    F1,
    q25 = NA_real_,
    q75 = NA_real_,
    n = NA_integer_,
    competitor_identity = method
  )

accuracy_plot_data <- bind_rows(
  accuracy_summary %>% mutate(competitor_identity = NA_character_),
  best_remaining
) %>%
  mutate(display = factor(display, method_levels))

write.csv(
  accuracy_plot_data %>% mutate(display = as.character(display)),
  file.path(OUT, "main_headline_accuracy.csv"),
  row.names = FALSE
)

p_accuracy <- ggplot(
  accuracy_plot_data,
  aes(moi, F1, color = display, group = display)
) +
  geom_ribbon(
    data = accuracy_plot_data %>% filter(!is.na(q25)),
    aes(ymin = q25, ymax = q75, fill = display),
    color = NA, alpha = 0.08, show.legend = FALSE
  ) +
  geom_line(
    data = accuracy_plot_data %>% filter(display == "fishash+"),
    linewidth = 1.35
  ) +
  geom_line(
    data = accuracy_plot_data %>% filter(!display %in% c(
      "fishash+", "best remaining published method (digitized)"
    )),
    linewidth = 0.95
  ) +
  geom_line(
    data = accuracy_plot_data %>% filter(display == "best remaining published method (digitized)"),
    linewidth = 0.8, linetype = "22"
  ) +
  geom_point(
    data = accuracy_plot_data %>% filter(display == "fishash+"),
    size = 2.5
  ) +
  geom_point(
    data = accuracy_plot_data %>% filter(display != "fishash+"),
    size = 1.8
  ) +
  scale_x_fishash_guide_load(secondary = FALSE, expand = expansion(mult = c(0.03, 0.05))) +
  scale_color_manual(
    values = method_cols, limits = method_levels, labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  scale_fill_manual(values = method_cols, limits = method_levels, drop = FALSE) +
  coord_cartesian(ylim = c(0.62, 1.0)) +
  labs(
    x = "MOI",
    y = "Median F1 (full subset)",
    title = "Accuracy across MOI"
  ) +
  report_theme() +
  guides(color = guide_legend(override.aes = list(linetype = method_types)))

# ---- Figure 1B: runtime at increasing library size --------------------------
runtime_summary <- read.csv(file.path(OUT, "sceptre_h2h", "sceptre_h2h.csv")) %>%
  filter(scenario == "varyNguides") %>%
  mutate(display = recode(
    method,
    `fishash+ Poisson` = "fishash+",
    sceptre_mixture = "sceptre mixture"
  )) %>%
  group_by(x, display) %>%
  summarize(secs = median(secs), .groups = "drop") %>%
  mutate(display = factor(display, method_levels))

speedups <- runtime_summary %>%
  mutate(display = as.character(display)) %>%
  pivot_wider(names_from = display, values_from = secs) %>%
  mutate(
    speedup = `sceptre mixture` / `fishash+`,
    label = paste0(round(speedup), "×"),
    hjust = case_when(x == min(x) ~ 0, x == max(x) ~ 1, TRUE ~ 0.5)
  )

write.csv(
  runtime_summary %>% mutate(display = as.character(display)),
  file.path(OUT, "main_headline_runtime.csv"),
  row.names = FALSE
)

p_runtime <- ggplot(
  runtime_summary,
  aes(x, secs, color = display, group = display)
) +
  geom_line(
    data = runtime_summary %>% filter(display == "fishash+"),
    linewidth = 1.2
  ) +
  geom_line(
    data = runtime_summary %>% filter(display == "sceptre mixture"),
    linewidth = 0.95
  ) +
  geom_point(size = 2.2) +
  geom_text(
    data = speedups,
    aes(x = x, y = `sceptre mixture`, label = label, hjust = hjust),
    inherit.aes = FALSE, vjust = -0.7, size = 5, fontface = "bold",
    color = "grey25"
  ) +
  scale_x_log10(breaks = c(20, 200, 2000, 20000), labels = scales::label_comma()) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.2))) +
  scale_color_manual(
    values = method_cols, limits = method_levels, labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  labs(
    x = "Number of guides (20k cells)",
    y = "Seconds per dataset (1 CPU)",
    title = "Runtime scaling"
  ) +
  report_theme() +
  guides(color = "none")

headline <- (p_accuracy | p_runtime) +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(
    title = "fishash+ combines field-leading accuracy with practical runtime",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 23),
      plot.tag = element_text(face = "bold", size = 20)
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(OUT, "main_headline.png"), headline, width = 13, height = 8.1, dpi = 160)

# ---- Figure 2A: matched methods under ambient-variance misspecification ------
overdispersed_summary <- read_method_panel("varyMOI_od") %>%
  filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  mutate(display = recode(
    method,
    fishash_refit10 = "fishash (refit10)",
    `fishash+_poisson` = "fishash+"
  )) %>%
  group_by(display, moi) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(display = factor(display, method_levels[1:2]))

write.csv(
  overdispersed_summary %>% mutate(display = as.character(display)),
  file.path(OUT, "main_boundary_overdispersion.csv"),
  row.names = FALSE
)

p_variance <- ggplot(
  overdispersed_summary,
  aes(moi, F1, color = display, group = display)
) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = display),
    alpha = 0.08, color = NA, show.legend = FALSE
  ) +
  geom_line(
    data = overdispersed_summary %>% filter(display == "fishash+"),
    linewidth = 1.35
  ) +
  geom_line(
    data = overdispersed_summary %>% filter(display != "fishash+"),
    linewidth = 0.95
  ) +
  geom_point(
    data = overdispersed_summary %>% filter(display == "fishash+"),
    size = 2.5
  ) +
  geom_point(
    data = overdispersed_summary %>% filter(display != "fishash+"),
    size = 1.8
  ) +
  scale_x_fishash_guide_load(
    secondary = FALSE,
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_color_manual(
    values = method_cols, limits = method_levels[1:2], labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  scale_fill_manual(values = method_cols, limits = method_levels[1:2], drop = FALSE) +
  coord_cartesian(ylim = c(0.15, 1.0)) +
  labs(
    x = "MOI",
    y = "Median F1 (full subset)",
    title = "NB-ambient setting",
    subtitle = "Matched fishash methods."
  ) +
  report_theme(base_size = 17) +
  theme(
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 30, hjust = 1),
    axis.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    legend.key.width = grid::unit(1.8, "lines"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(
      size = 16, color = "grey30", lineheight = 1.05,
      margin = margin(b = 8)
    )
  ) +
  guides(color = guide_legend(override.aes = list(linetype = method_types[1:2])))

# ---- Figure 2B: harder ambient-mean setting ---------------------------------
derived_exact <- bind_rows(
  read_method_panel("varyMOILow") %>%
    filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
    transmute(
      moi,
      display = recode(
        method,
        fishash_refit10 = "fishash (refit10)",
        `fishash+_poisson` = "fishash+"
      ),
      F1
    ),
  read.csv(file.path(OUT, "endo_contrast", "sceptre_scores.csv")) %>%
    filter(scn == "varyMOILow", method == "sceptre_mixture") %>%
    transmute(moi, display = "sceptre mixture", F1)
) %>%
  group_by(moi, display) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(competitor_identity = NA_character_)

harder_plot_data <- derived_exact %>%
  mutate(display = factor(as.character(display), method_levels[1:3]))

write.csv(
  harder_plot_data %>% mutate(display = as.character(display)),
  file.path(OUT, "main_boundary_exogenous.csv"),
  row.names = FALSE
)

p_harder <- ggplot(
  harder_plot_data,
  aes(moi, F1, color = display, group = display)
) +
  geom_ribbon(
    data = harder_plot_data %>% filter(!is.na(q25)),
    aes(ymin = q25, ymax = q75, fill = display),
    color = NA, alpha = 0.08, show.legend = FALSE
  ) +
  geom_line(
    data = harder_plot_data %>% filter(display == "fishash+"),
    linewidth = 1.35
  ) +
  geom_line(
    data = harder_plot_data %>% filter(display != "fishash+"),
    linewidth = 0.95
  ) +
  geom_point(
    data = harder_plot_data %>% filter(display == "fishash+"),
    size = 2.5
  ) +
  geom_point(
    data = harder_plot_data %>% filter(display != "fishash+"),
    size = 1.8
  ) +
  scale_x_fishash_guide_load(
    secondary = FALSE,
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_color_manual(
    values = method_cols, limits = method_levels[1:3], labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  scale_fill_manual(values = method_cols, limits = method_levels[1:3], drop = FALSE) +
  coord_cartesian(ylim = c(0.15, 1.0)) +
  labs(
    x = "MOI",
    y = "Median F1 (full subset)",
    title = "Increased exogenous setting",
    subtitle = "75% exogenous ambient\n20 guide UMIs/cell; SNR 1."
  ) +
  report_theme(base_size = 17) +
  theme(
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 30, hjust = 1),
    axis.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    legend.key.width = grid::unit(1.8, "lines"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(
      size = 16, color = "grey30", lineheight = 1.05,
      margin = margin(b = 8)
    )
  ) +
  guides(color = guide_legend(override.aes = list(linetype = method_types[1:3])))

boundary <- (p_variance | p_harder) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = "Two distinct robustness conditions",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 23),
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(file.path(OUT, "main_boundary_conditions.png"), boundary, width = 11.5, height = 7.5, dpi = 160)

# ---- Supporting analysis: four-way guide-count comparison ------------------
guide_exact <- bind_rows(
  read_guide_panel("varyNguides"),
  read_guide_panel("varyNguidesLow")
) %>%
  filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  filter(
    scenario == "varyNguidesLow" |
      (scenario == "varyNguides" & ng < 20000 & iter <= 5) |
      (scenario == "varyNguides" & ng == 20000 & iter == 1)
  ) %>%
  mutate(
    regime = recode(
      scenario,
      varyNguides = "high-gRNA (100 UMI, SNR 4)",
      varyNguidesLow = "low-gRNA (20 UMI, SNR 1)"
    ),
    display = recode(
      method,
      fishash_refit10 = "fishash (refit10)",
      `fishash+_poisson` = "fishash+"
    )
  ) %>%
  group_by(regime, ng, display) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "exact", competitor_identity = NA_character_)

guide_sceptre_high <- read.csv(file.path(OUT, "sceptre_h2h", "sceptre_h2h.csv")) %>%
  filter(scenario == "varyNguides", method == "sceptre_mixture") %>%
  group_by(x) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  transmute(
    regime = "high-gRNA (100 UMI, SNR 4)",
    ng = x,
    display = "sceptre mixture",
    q25, q75, n, F1,
    source = "exact",
    competitor_identity = NA_character_
  )

guide_digitized <- read.csv(file.path(OUT, "varyNguides", "fig5_overlay_data.csv")) %>%
  filter(method %in% c(
    "cleanser_cs", "cleanser_dc", "crispat_gauss", "crispat_poisgauss",
    "crispat_poisson", "crispat_negbinom", "sceptre_mixture", "demuxem",
    "geomux", "fishash (digitized)"
  ))

guide_sceptre_low <- guide_digitized %>%
  filter(regime == "low_gRNA", method == "sceptre_mixture") %>%
  transmute(
    regime = "low-gRNA (20 UMI, SNR 1)",
    ng,
    display = "sceptre mixture",
    q25 = NA_real_, q75 = NA_real_, n = NA_integer_, F1,
    source = "paper (digitized)",
    competitor_identity = NA_character_
  )

guide_best <- guide_digitized %>%
  filter(!method %in% c("fishash (digitized)", "sceptre_mixture")) %>%
  group_by(regime, ng) %>%
  slice_max(F1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    regime = recode(
      regime,
      high_gRNA = "high-gRNA (100 UMI, SNR 4)",
      low_gRNA = "low-gRNA (20 UMI, SNR 1)"
    ),
    ng,
    display = "best remaining published method (digitized)",
    q25 = NA_real_, q75 = NA_real_, n = NA_integer_, F1,
    source = "paper (digitized)",
    competitor_identity = method
  )

guide_fourway <- bind_rows(guide_exact, guide_sceptre_high, guide_sceptre_low, guide_best) %>%
  mutate(
    regime = factor(
      regime,
      c("high-gRNA (100 UMI, SNR 4)", "low-gRNA (20 UMI, SNR 1)")
    ),
    display = factor(display, method_levels)
  )

write.csv(
  guide_fourway %>% mutate(regime = as.character(regime), display = as.character(display)),
  file.path(OUT, "appendix_guide_count.csv"),
  row.names = FALSE
)

p_guide <- ggplot(guide_fourway, aes(ng, F1, color = display, group = display)) +
  geom_ribbon(
    data = guide_fourway %>% filter(!is.na(q25)),
    aes(ymin = q25, ymax = q75, fill = display),
    alpha = 0.08, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_line(
    data = guide_fourway %>% filter(display == "best remaining published method (digitized)"),
    linewidth = 0.8, linetype = "22"
  ) +
  geom_point(size = 1.9) +
  facet_wrap(~regime, nrow = 1) +
  scale_x_log10(breaks = c(20, 200, 2000, 20000), labels = scales::label_comma()) +
  scale_color_manual(
    values = method_cols, limits = method_levels, labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  scale_fill_manual(values = method_cols, limits = method_levels, drop = FALSE) +
  coord_cartesian(ylim = c(0.52, 1.0)) +
  labs(
    x = "Number of guides (MOI = 1.04)",
    y = "Median F1 (full subset)",
    title = "Four-way comparison as the guide library grows",
    subtitle = paste(
      "High-gRNA exact curves use five common seeds through 2,000 guides and one at 20,000.",
      "Low-gRNA sceptre and both grey curves are digitized from the paper.",
      sep = "\n"
    )
  ) +
  report_theme() +
  guides(color = guide_legend(override.aes = list(linetype = method_types)))

ggsave(file.path(OUT, "appendix_guide_count.png"), p_guide, width = 12.5, height = 7.4, dpi = 160)

# ---- Supporting analysis: method levels in the paired mean regimes ---------
mean_fish <- bind_rows(
  vary_moi,
  read_method_panel("varyMOILow")
) %>%
  filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  mutate(
    regime = recode(
      scenario,
      varyMOI = "25% exogenous (high-gRNA)",
      varyMOILow = "75% exogenous (low-gRNA)"
    ),
    display = recode(
      method,
      fishash_refit10 = "fishash (refit10)",
      `fishash+_poisson` = "fishash+"
    )
  ) %>%
  group_by(regime, moi, display) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "exact", competitor_identity = NA_character_)

mean_sceptre <- read.csv(file.path(OUT, "endo_contrast", "sceptre_scores.csv")) %>%
  group_by(scn, moi) %>%
  summarize(
    q25 = quantile(F1, 0.25, na.rm = TRUE),
    q75 = quantile(F1, 0.75, na.rm = TRUE),
    n = sum(!is.na(F1)),
    F1 = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  transmute(
    regime = recode(
      scn,
      varyMOI = "25% exogenous (high-gRNA)",
      varyMOILow = "75% exogenous (low-gRNA)"
    ),
    moi,
    display = "sceptre mixture",
    q25, q75, n, F1,
    source = "exact",
    competitor_identity = NA_character_
  )

mean_best <- best_remaining %>%
  transmute(
    regime = "25% exogenous (high-gRNA)",
    moi,
    display,
    q25, q75, n, F1,
    source = "paper (digitized)",
    competitor_identity
  )

mean_fourway <- bind_rows(mean_fish, mean_sceptre, mean_best) %>%
  mutate(
    regime = factor(
      regime,
      c("25% exogenous (high-gRNA)", "75% exogenous (low-gRNA)")
    ),
    display = factor(display, method_levels)
  )

write.csv(
  mean_fourway %>% mutate(regime = as.character(regime), display = as.character(display)),
  file.path(OUT, "appendix_mean_regimes.csv"),
  row.names = FALSE
)

p_mean_levels <- ggplot(mean_fourway, aes(moi, F1, color = display, group = display)) +
  geom_ribbon(
    data = mean_fourway %>% filter(!is.na(q25)),
    aes(ymin = q25, ymax = q75, fill = display),
    alpha = 0.08, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_line(
    data = mean_fourway %>% filter(display == "best remaining published method (digitized)"),
    linewidth = 0.8, linetype = "22"
  ) +
  geom_point(size = 1.9) +
  facet_wrap(~regime, nrow = 1) +
  scale_x_fishash_guide_load(secondary = FALSE) +
  scale_color_manual(
    values = method_cols, limits = method_levels, labels = method_labels,
    drop = FALSE, name = NULL
  ) +
  scale_fill_manual(values = method_cols, limits = method_levels, drop = FALSE) +
  coord_cartesian(ylim = c(0.15, 1.0)) +
  labs(
    x = "MOI",
    y = "Median F1 (full subset)",
    title = "Four-way comparison across ambient-mean regimes",
    subtitle = paste(
      "Exact curves use 10 simulations; grey is available only in the published 25% regime.",
      "The 75% varying-MOI sweep was constructed for this analysis and has no published competitor values.",
      sep = "\n"
    )
  ) +
  report_theme() +
  guides(color = guide_legend(override.aes = list(linetype = method_types)))

ggsave(file.path(OUT, "appendix_mean_regimes.png"), p_mean_levels, width = 12.5, height = 7.8, dpi = 160)

# ---- Supporting analysis: full-curve ranking -------------------------------
auprc <- read.csv(
  file.path(OUT, "varyMOI", "auprc_by_moi.csv"),
  check.names = FALSE
) %>%
  pivot_longer(-moi, names_to = "raw_method", values_to = "AUPRC") %>%
  filter(raw_method != "raw_count") %>%
  mutate(display = ifelse(raw_method == "fishash", "fishash (refit10)", "fishash+")) %>%
  group_by(moi, display) %>%
  summarize(
    q25 = quantile(AUPRC, 0.25, na.rm = TRUE),
    q75 = quantile(AUPRC, 0.75, na.rm = TRUE),
    n = sum(!is.na(AUPRC)),
    AUPRC = median(AUPRC, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(auprc, file.path(OUT, "appendix_auprc.csv"), row.names = FALSE)

p_auprc <- ggplot(auprc, aes(moi, AUPRC, color = display, fill = display, group = display)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.09, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 2.5) +
  scale_x_fishash_guide_load(secondary = FALSE) +
  scale_color_manual(
    values = method_cols[c("fishash+", "fishash (refit10)")],
    labels = method_labels[c("fishash+", "fishash (refit10)")],
    name = NULL
  ) +
  scale_fill_manual(values = method_cols[c("fishash+", "fishash (refit10)")]) +
  coord_cartesian(ylim = c(0.9, 1.0)) +
  labs(
    x = "MOI",
    y = "Median full AUPRC",
    title = "fishash+ improves full-curve ranking as MOI increases",
    subtitle = "Medians and IQRs across the 10 author-seeded simulations."
  ) +
  report_theme() +
  guides(fill = "none")

ggsave(
  file.path(OUT, "varyMOI", "fig6c_auprc_vs_moi.png"),
  p_auprc, width = 11, height = 7.2, dpi = 160
)

# ---- Main text: detailed overdispersion decomposition -----------------------
overdispersion <- read_method_panel("varyMOI_od") %>%
  filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  mutate(display = recode(
    method,
    fishash_refit10 = "fishash (refit10)",
    `fishash+_poisson` = "fishash+"
  )) %>%
  select(moi, iter, display, Precision, Recall, F1) %>%
  pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "value") %>%
  group_by(moi, display, stat) %>%
  summarize(
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    n = sum(!is.na(value)),
    value = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))

write.csv(overdispersion, file.path(OUT, "appendix_overdispersion.csv"), row.names = FALSE)

p_overdispersion <- ggplot(
  overdispersion,
  aes(moi, value, color = display, fill = display, group = display)
) +
  geom_hline(
    data = data.frame(stat = factor("Precision", c("Precision", "Recall", "F1"))),
    aes(yintercept = 0.95), inherit.aes = FALSE,
    linetype = "22", color = "grey35", linewidth = 0.7
  ) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.09, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.3) +
  facet_wrap(~stat, nrow = 1) +
  scale_x_fishash_guide_load(
    secondary = FALSE,
    breaks = fishash_preselection_moi[c(1, 3, 4, 5, 7, 8)]
  ) +
  scale_color_manual(
    values = method_cols[c("fishash+", "fishash (refit10)")],
    labels = method_labels[c("fishash+", "fishash (refit10)")],
    name = NULL
  ) +
  scale_fill_manual(values = method_cols[c("fishash+", "fishash (refit10)")]) +
  coord_cartesian(ylim = c(0.45, 1.0)) +
  labs(
    x = "MOI",
    y = "Median (full subset)",
    title = "NB ambient exposes a low-MOI precision failure"
  ) +
  report_theme() +
  guides(fill = "none")

ggsave(
  file.path(OUT, "varyMOI_od", "figC2_overdispersed_prf.png"),
  p_overdispersion, width = 13, height = 7.4, dpi = 160
)

# ---- Supporting analysis: mechanism ablation -------------------------------
ablation_levels <- c(
  "fishash+", "fitted-depth contingency", "fishash (refit10)",
  "Poisson tail, raw cell total", "fishash, no refit"
)
ablation_cols <- c(
  "fishash+" = "#D62728",
  "fitted-depth contingency" = "#E68613",
  "fishash (refit10)" = "#6A3D9A",
  "Poisson tail, raw cell total" = "grey40",
  "fishash, no refit" = "#2CA02C"
)

ablation <- vary_moi %>%
  filter(method %in% c(
    "fishash_refit10", "fishash+_poisson", "fishash+_contingency",
    "fishash_poisson_rawd", "fishash_refit0"
  )) %>%
  mutate(display = recode(
    method,
    fishash_refit10 = "fishash (refit10)",
    `fishash+_poisson` = "fishash+",
    `fishash+_contingency` = "fitted-depth contingency",
    fishash_poisson_rawd = "Poisson tail, raw cell total",
    fishash_refit0 = "fishash, no refit"
  )) %>%
  select(moi, iter, display, Precision, Recall, F1) %>%
  pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "value") %>%
  group_by(moi, display, stat) %>%
  summarize(
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    n = sum(!is.na(value)),
    value = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    display = factor(display, ablation_levels),
    stat = factor(stat, c("Precision", "Recall", "F1"))
  )

write.csv(
  ablation %>% mutate(display = as.character(display), stat = as.character(stat)),
  file.path(OUT, "appendix_ablation.csv"), row.names = FALSE
)

p_ablation <- ggplot(ablation, aes(moi, value, color = display, group = display)) +
  geom_line(aes(linewidth = display)) +
  geom_point(size = 2.1) +
  facet_wrap(~stat, nrow = 1) +
  scale_x_fishash_guide_load(
    secondary = FALSE,
    breaks = fishash_preselection_moi[c(1, 3, 4, 5, 7, 8)]
  ) +
  scale_color_manual(values = ablation_cols, limits = ablation_levels, name = NULL) +
  scale_linewidth_manual(
    values = c(1.35, 1.0, 1.0, 0.85, 0.85),
    limits = ablation_levels, guide = "none"
  ) +
  coord_cartesian(ylim = c(0.45, 1.0)) +
  labs(
    x = "MOI",
    y = "Median (full subset)",
    title = "The fitted ambient depth, not the tail family, drives the gain",
    subtitle = "Medians across 10 author-seeded simulations."
  ) +
  report_theme() +
  guides(color = guide_legend(nrow = 2))

ggsave(
  file.path(OUT, "varyMOI", "varyMOI_ablation_prf_vs_x.png"),
  p_ablation, width = 13, height = 8.2, dpi = 160
)

message(
  "wrote report and supporting analysis figures plus source CSVs to ", OUT
)
