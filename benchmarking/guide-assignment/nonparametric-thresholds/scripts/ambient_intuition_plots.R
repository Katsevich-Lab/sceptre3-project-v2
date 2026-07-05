#!/usr/bin/env Rscript
# Intuition plots for ambient DEPTH (per cell) and COMPOSITION (per guide) on
# crispat_schraivogel. Reads results/ambient_intuition/*.csv (from the compute
# script). Two figures: depth.png and composition.png.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
have_pw <- requireNamespace("patchwork", quietly = TRUE)
if (have_pw) library(patchwork)

HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
D <- file.path(HERE, "results", "ambient_intuition")

cells  <- read.csv(file.path(D, "per_cell.csv"))
cells  <- cells[cells$libsize > 0, ]              # drop empty droplets for log scatters
guides <- read.csv(file.path(D, "per_guide.csv"))
prof   <- read.csv(file.path(D, "cell_profiles.csv"))
gs     <- read.csv(file.path(D, "guide_scatter.csv"))
hl     <- readRDS(file.path(D, "highlights.rds"))

sig_cell <- cells[cells$cell == hl$cell_signal, ]
amb_cell <- cells[cells$cell == hl$cell_ambient, ]
g_inf    <- guides[guides$guide == hl$guide_inflated, ]
short <- function(g) sub("CROPseq_dCas9_DS_", "", sub("\\..*$|-P1P2$", "", g))
th <- theme_bw(base_size = 11) + theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5), legend.position = "none")

# ============================ FIGURE 1: DEPTH ===============================
rho1 <- cor(cells$libsize, cells$depth_fish)
p1 <- ggplot(cells, aes(libsize, depth_fish + 0.5)) +
  geom_point(aes(colour = top_frac), alpha = 0.10, size = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(data = sig_cell, aes(libsize, depth_fish + 0.5), colour = "#d7301f", size = 3) +
  geom_point(data = amb_cell, aes(libsize, depth_fish + 0.5), colour = "#0F6E56", size = 3) +
  annotate("text", x = sig_cell$libsize, y = sig_cell$depth_fish + 1.4,
           label = "signal-dominated cell", colour = "#d7301f", size = 3, hjust = 1.05) +
  annotate("text", x = amb_cell$libsize, y = amb_cell$depth_fish + 2.5,
           label = "ambient-only cell", colour = "#0F6E56", size = 3, hjust = -0.05) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_viridis_c(option = "C", end = 0.9) +
  labs(title = "Library size is NOT ambient depth",
       subtitle = sprintf("each point = a cell;  dashed = y=x;  r = %.2f;  median libsize/depth = %d x",
                           rho1, round(median(cells$libsize / pmax(cells$depth_fish,1)))),
       x = "naive gRNA library size (UMIs)", y = "fishash ambient depth  kappa_c") + th

rho2 <- cor(cells$depth_clean, cells$depth_fish)
p2 <- ggplot(cells, aes(depth_clean + 0.5, depth_fish + 0.5)) +
  geom_point(alpha = 0.08, size = 0.4, colour = "#185FA5") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Two ambient-depth estimators agree",
       subtitle = sprintf("CLEANSER sub-threshold sum vs fishash rank-1 depth;  r = %.2f", rho2),
       x = "CLEANSER depth (sum of gRNA UMIs <= 2)", y = "fishash ambient depth  kappa_c") + th

# panel 3: one signal cell's per-guide counts vs the ambient floor gamma_g*kappa_c
pf <- prof %>% arrange(desc(count_signalcell)) %>% mutate(rank = row_number())
p3 <- ggplot(pf, aes(rank)) +
  geom_segment(aes(xend = rank, y = amb_signalcell + 0.3, yend = count_signalcell + 0.3),
               colour = "grey75") +
  geom_point(aes(y = amb_signalcell + 0.3), colour = "#0F6E56", size = 1.3) +
  geom_point(aes(y = count_signalcell + 0.3,
                 colour = assigned_signalcell, size = assigned_signalcell)) +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = "#d7301f")) +
  scale_size_manual(values = c(`FALSE` = 0.9, `TRUE` = 2.6)) +
  scale_y_log10() + expand_limits(x = 95) +
  annotate("text", x = 30, y = pf$count_signalcell[1] * 0.55,
           label = sprintf("assigned guide:\n%s\n%d UMIs", short(pf$guide[1]), pf$count_signalcell[1]),
           colour = "#d7301f", size = 3, hjust = 0) +
  annotate("text", x = 30, y = 1.5, label = "green = ambient floor",
           colour = "#0F6E56", size = 3, hjust = 0) +
  labs(title = "What the test sees in ONE cell",
       subtitle = sprintf("library %d UMIs but ambient depth %.1f; guides ranked by count",
                          sig_cell$libsize, sig_cell$depth_fish),
       x = "guide (ranked)", y = "UMI count (+0.3, log)") + th

# ========================= FIGURE 2: COMPOSITION ============================
gg <- guides %>% arrange(desc(comp_fish)) %>% mutate(rank = row_number())
top_amb <- gg$guide[1]
pc1 <- ggplot(gg, aes(rank)) +
  geom_col(aes(y = comp_fish), fill = "#7F77DD", width = 0.85) +
  geom_point(aes(y = share_naive), colour = "#993C1D", size = 1) +
  annotate("text", x = 60, y = max(gg$comp_fish)*0.85,
           label = "bars = ambient composition (fishash)\ndots = naive global share",
           size = 3, hjust = 0, colour = "grey30") +
  annotate("text", x = 3, y = gg$comp_fish[1] + 0.004, label = short(top_amb), size = 2.8, hjust = 0) +
  labs(title = "Composition of the ambient gRNA pool",
       subtitle = "86 guides ranked by ambient share; naive share (dots) tracks it closely here",
       x = "guide (ranked by ambient share)", y = "share of background UMIs") + th

pc2 <- ggplot(guides, aes(comp_fish, share_naive)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(alpha = 0.6, size = 1.6, colour = "#185FA5") +
  geom_point(data = g_inf, aes(comp_fish, share_naive), colour = "#d7301f", size = 3) +
  annotate("text", x = g_inf$comp_fish, y = g_inf$share_naive*1.12,
           label = sprintf("%s\nshare %.1f%% vs ambient %.1f%%", short(hl$guide_inflated),
                           100*g_inf$share_naive, 100*g_inf$comp_fish),
           colour = "#d7301f", size = 2.8, hjust = 0.5) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Per guide: naive share vs ambient rate",
       subtitle = "above the line = global share overstates the guide's true ambient rate",
       x = "ambient composition  gamma_g (fishash)", y = "naive global share  N_g./N") + th

# panel: the highlighted guide across all cells, count vs ambient depth
gamma_g <- g_inf$comp_fish
gs2 <- gs %>% mutate(grp = ifelse(assigned, "assigned (signal)", "unassigned (noise)"))
line_df <- data.frame(d = exp(seq(log(0.5), log(max(gs$depth_fish)), length = 100)))
line_df$x <- line_df$d + 0.5; line_df$y <- gamma_g * line_df$d + 0.5
pc3 <- ggplot(gs2, aes(depth_fish + 0.5, count + 0.5)) +
  geom_jitter(data = filter(gs2, !assigned), width = 0.06, height = 0.08,
              alpha = 0.10, size = 0.4, colour = "grey45") +
  geom_jitter(data = filter(gs2, assigned), width = 0.06, height = 0.0,
              alpha = 0.7, size = 1.0, colour = "#d7301f") +
  geom_line(data = line_df, aes(x, y), colour = "#0F6E56", linewidth = 0.9) +
  annotate("text", x = 2, y = 150, label = "signal cells:\nfar above ambient line",
           colour = "#d7301f", size = 3, hjust = 0) +
  annotate("text", x = 20, y = 1.5, label = "ambient line  gamma_g * kappa_c",
           colour = "#0F6E56", size = 3, hjust = 0) +
  scale_x_log10() + scale_y_log10() +
  labs(title = sprintf("One guide (%s) across all cells", short(hl$guide_inflated)),
       subtitle = "noise hugs the ambient line; signal cells break far above it",
       x = "ambient depth  kappa_c", y = "guide UMI count") + th

if (have_pw) {
  ggsave(file.path(D, "depth.png"), p1 + p2 + p3 + plot_layout(nrow = 1),
         width = 14, height = 4.3, dpi = 150)
  ggsave(file.path(D, "composition.png"), pc1 + pc2 + pc3 + plot_layout(nrow = 1),
         width = 14, height = 4.3, dpi = 150)
} else {
  for (nm in c("p1","p2","p3","pc1","pc2","pc3"))
    ggsave(file.path(D, paste0(nm, ".png")), get(nm), width = 5, height = 4.3, dpi = 150)
}
cat("wrote depth.png and composition.png to", D, "\n")
