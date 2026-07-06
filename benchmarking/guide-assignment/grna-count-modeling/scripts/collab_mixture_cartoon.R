## ============================================================================
## Collaborator writeup — CONCEPTUAL CARTOON (schematic, NOT real data).
##
## The central picture: a single guide's observed gRNA-count distribution is a
## MIXTURE OF THREE SOURCES, drawn here by hand on a log-count axis:
##   1. Ambient contamination  (blue)  — steep decay peaking at count ~1-2; the
##      tall bulk, soup landing in every cell.
##   2. Weak integrations      (gold)  — a broad, lower hump at low-to-moderate
##      counts (~2-25) that OVERLAPS the ambient tail; this is why the low mode
##      is "partly real" — real integrations hide INSIDE it.
##   3. Strong integrations    (green) — a clean hump at high counts (~100-1000),
##      well separated from the low mode by an empty gap.
##
## The reader should see: low mode = ambient + weak (overlapping), an empty gap,
## then the high mode = strong; plus a thin outline of the observed TOTAL (sum of
## the three) = what actually gets measured.
##
## Everything below is constructed BY HAND — no dataset is read. Run from
## grna-count-modeling/.  Packages: ggplot2 (+ scales).
## ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(scales)})
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- Okabe-Ito palette (matches the rest of the writeup) --------------------
COL_AMB    <- "#0072B2"   # blue  — ambient contamination
COL_WEAK   <- "#E69F00"   # gold  — weak integrations
COL_STRONG <- "#009E73"   # green — strong integrations
COL_TOTAL  <- "grey25"    # observed total outline

## ---- x-grid on the LOG-count axis -------------------------------------------
## Work in log10(count); count runs ~1..2000 so both modes + the gap are visible.
lx <- seq(log10(1), log10(2500), length.out = 1600)   # log10 count
x  <- 10^lx

## ---- three component densities, defined by hand in log-count space ----------
## Gaussians-in-log give clean, tunable humps; heights are schematic (relative
## frequency), so overall normalization does not matter.
lognorm_bump <- function(lx, center_count, log_sd, height) {
  height * exp(-0.5 * ((lx - log10(center_count)) / log_sd)^2)
}

## 1. Ambient: tall, narrow, peaked at the very lowest counts (~1.3), with a
##    long-ish right tail into the low-moderate counts (the tail weak sits in).
amb <- lognorm_bump(lx, center_count = 1.3, log_sd = 0.24, height = 1.00)

## 2. Weak integrations: broad, much lower hump centred around count ~7, spanning
##    ~2-25 so it clearly OVERLAPS the ambient right tail (partly hidden inside it).
weak <- lognorm_bump(lx, center_count = 7, log_sd = 0.42, height = 0.28)

## 3. Strong integrations: a clean, moderate hump at high counts (~300), well to
##    the right of everything else, separated by an empty gap.
strong <- lognorm_bump(lx, center_count = 300, log_sd = 0.28, height = 0.55)

total <- amb + weak + strong

## long-form data frame for the filled component ribbons
comp_levels <- c("Ambient contamination", "Weak integrations", "Strong integrations")
df <- rbind(
  data.frame(lx = lx, x = x, y = amb,    comp = comp_levels[1]),
  data.frame(lx = lx, x = x, y = weak,   comp = comp_levels[2]),
  data.frame(lx = lx, x = x, y = strong, comp = comp_levels[3])
)
df$comp <- factor(df$comp, levels = comp_levels)
df_total <- data.frame(lx = lx, x = x, y = total)

## ---- locate the empty gap (between weak's right shoulder and strong's left) --
## purely for the bracket/label geometry; pick where total density is lowest
## between the two modes.
in_gap <- x > 25 & x < 150
gap_x  <- x[in_gap]
gap_lo <- min(gap_x); gap_hi <- max(gap_x)
gap_mid <- 10^mean(log10(c(gap_lo, gap_hi)))

pal <- c(
  "Ambient contamination" = COL_AMB,
  "Weak integrations"     = COL_WEAK,
  "Strong integrations"   = COL_STRONG
)

ytop <- max(total) * 1.06
## bracket geometry (drawn just under the axis top, above the curves)
low_l <- 1; low_r <- 27; hi_l <- 130; hi_r <- 2200
br_y  <- ytop * 0.985                 # bracket rail height
tick  <- ytop * 0.026                 # bracket end-tick drop
lab_y <- ytop * 0.905                 # bracket text height

brackets <- data.frame(
  x    = c(low_l, hi_l),
  xend = c(low_r, hi_r),
  y    = br_y
)
bracket_ticks <- data.frame(
  x    = c(low_l, low_r, hi_l, hi_r),
  yend = br_y - tick,
  y    = br_y
)
bracket_labs <- data.frame(
  x   = c(10^mean(log10(c(low_l, low_r))), 10^mean(log10(c(hi_l, hi_r)))),
  y   = lab_y,
  lab = c("low mode", "high mode")
)

p <- ggplot() +
  ## filled component densities (semi-transparent, stacked by draw order)
  geom_area(data = df, aes(x = x, y = y, fill = comp),
            position = "identity", alpha = 0.55, colour = NA) +
  ## crisp component outlines so overlapping humps stay legible
  geom_line(data = df, aes(x = x, y = y, colour = comp),
            linewidth = 0.7, show.legend = FALSE) +
  ## thin outline of the OBSERVED total (what actually gets measured)
  geom_line(data = df_total, aes(x = x, y = y),
            colour = COL_TOTAL, linewidth = 0.7, linetype = "22") +
  ## ---- direct, colour-matched labels on each component (no legend) ----------
  annotate("text", x = 1.02, y = 0.83, hjust = 0,
           label = "Ambient\ncontamination", colour = COL_AMB,
           fontface = "bold", size = 4.3, lineheight = 0.9) +
  annotate("text", x = 6.7, y = 0.635, hjust = 0.5,
           label = "Weak integrations", colour = "#B37400",
           fontface = "bold", size = 4.3) +
  annotate("text", x = 6.7, y = 0.565, hjust = 0.5,
           label = "(real, but hidden\ninside the low mode)", colour = "#B37400",
           fontface = "italic", size = 3.6, lineheight = 0.9) +
  annotate("segment", x = 6.7, xend = 6.7, y = 0.335, yend = 0.50,
           colour = COL_WEAK, linewidth = 0.5,
           arrow = arrow(length = unit(0.13, "cm"), ends = "first")) +
  annotate("text", x = 300, y = 0.60, hjust = 0.5,
           label = "Strong integrations", colour = "#00785A",
           fontface = "bold", size = 4.3) +
  ## label the observed-total envelope on its low-mode right shoulder, with a
  ## short leader to the dashed curve (kept clear of the bracket labels above)
  annotate("segment", x = 17, xend = 13.5, y = 0.46, yend = 0.335,
           colour = COL_TOTAL, linewidth = 0.4) +
  annotate("text", x = 17.5, y = 0.475,
           label = "observed total\n(what gets measured)", colour = COL_TOTAL,
           fontface = "italic", size = 3.7, hjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000),
                labels = label_number(accuracy = 1),
                expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "What we think a guide's counts are: a mixture of three sources",
    x = "gRNA UMI count in the cell  (log scale)",
    y = "relative frequency  (number of cells)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank()
  )

ggsave(file.path(OUT, "mixture_cartoon.png"), p, width = 9, height = 5, dpi = 130)
cat("wrote", file.path(OUT, "mixture_cartoon.png"), "\n")
