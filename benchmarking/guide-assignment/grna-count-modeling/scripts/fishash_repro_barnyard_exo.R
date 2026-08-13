#!/usr/bin/env Rscript
# ============================================================================
# What is a REALISTIC exogenous-noise fraction? -- estimated from the barnyard.
#
# In a species-mixing (barnyard) sample, a guide of the WRONG species in a cell is
# contamination. Decompose it into the guidebender model's two noise types (its exo
# model, paper Eq 6, is barcode-swap / index-hop; endo Eq 5 is diffuse ambient soup):
#   - endogenous soup  -> diffuse, Poisson, wrong-species counts ~1-2 in every cell
#   - exogenous swap    -> a few cells with ONE wrong-species guide at high count
# We gate to GEX-pure cells (frac_homo > 0.99), which EXCLUDES doublets (a doublet has
# real wrong-species mRNA -> frac_homo drops well below 0.99). So this is the swap-only,
# post-doublet exo fraction -- the apples-to-apples comparison to guidebender's exo,
# which is also doublet-free (guidebender has no doublet term at all).
#
# Finding (use the NON-co-cultured mix0hr control): exo ~= 0% (CROP-seq) to ~10-15%
# (direct-capture), concentrated in ~a handful of high-count cells -- FAR below the paper's
# 25% (high-gRNA) / 75% (low-gRNA) settings. So realistic data is endo-dominated: the regime
# where fishash+ (depth_fix) wins. CAVEAT: the 72h co-cultured mix72hr samples show much more
# (~55% direct-capture) -- that is co-culture cross-contamination / residual aggregates, NOT a
# clean index-hop read, so mix0hr is the reliable estimate.
# Data: external/repro_work/mix{0hr,72hr}_{DirectCapture,Cropseq}_*. Run from folder root.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})
REPRO <- "external/repro_work"; OUT <- "results/fishash_repro"
samples <- c("mix0hr_DirectCapture", "mix72hr_DirectCapture", "mix0hr_Cropseq", "mix72hr_Cropseq")

rows <- list()
for (s in samples) {
  mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (!file.exists(mtx)) next
  counts <- as(readMM(mtx), "CsparseMatrix")
  g    <- read.csv(file.path(REPRO, paste0(s, "_guides.csv")))
  meta <- read.csv(file.path(REPRO, paste0(s, "_meta.csv")))
  fh <- meta$homo_sum_gex / (meta$homo_sum_gex + meta$mus_sum_gex)
  hc <- which(fh > 0.99)                                     # GEX-pure human cells (doublet-free)
  M  <- counts[g$guide_type == "mus_guide", hc, drop = FALSE] # wrong-species (noise) counts
  tot <- sum(M); x <- M@x
  ex <- function(t) 100 * sum(x[x >= t]) / tot               # % of noise UMIs in count>=t events
  rows[[s]] <- data.frame(
    sample = s, pure_human_cells = length(hc), noise_umis = tot, median_nonzero = median(x),
    n_events_ge10 = sum(x >= 10),
    exo_pct_ge10 = round(ex(10), 1), exo_pct_ge5 = round(ex(5), 1), exo_pct_ge3 = round(ex(3), 1))
}
df <- do.call(rbind, rows); rownames(df) <- NULL
write.csv(df, file.path(OUT, "barnyard_exo_fraction.csv"), row.names = FALSE)

# Plot the count>=10 estimate as the cleanest barcode-swap signal. Vertical
# ranges show sensitivity to using count>=3, >=5, or >=10 as the cutoff.
plot_df <- df %>%
  mutate(
    chemistry = ifelse(grepl("DirectCapture", sample), "direct capture", "CROP-seq"),
    phase = ifelse(grepl("mix0hr", sample), "0 h control", "72 h co-culture"),
    label = factor(
      paste(phase, chemistry, sep = "\n"),
      levels = c("0 h control\nCROP-seq", "0 h control\ndirect capture",
                 "72 h co-culture\nCROP-seq", "72 h co-culture\ndirect capture")
    ),
    lower = pmin(exo_pct_ge3, exo_pct_ge5, exo_pct_ge10),
    upper = pmax(exo_pct_ge3, exo_pct_ge5, exo_pct_ge10)
  )
p <- ggplot(plot_df, aes(label, exo_pct_ge10, color = chemistry)) +
  geom_hline(yintercept = c(25, 75), linetype = "22", color = "grey65") +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = paste0(exo_pct_ge10, "%")), nudge_y = 4.5,
            show.legend = FALSE, fontface = "bold", size = 3.5) +
  annotate("text", x = 4.42, y = 25, label = "paper: 25% exogenous",
           hjust = 1, vjust = -0.5, color = "grey40", size = 3.2) +
  annotate("text", x = 4.42, y = 75, label = "paper: 75% exogenous",
           hjust = 1, vjust = -0.5, color = "grey40", size = 3.2) +
  scale_color_manual(values = c("CROP-seq" = "#1f78b4", "direct capture" = "#e31a1c")) +
  coord_cartesian(ylim = c(0, 82), clip = "off") +
  labs(
    x = NULL,
    y = "estimated exogenous share of wrong-species guide UMIs",
    color = NULL,
    title = "Barnyard controls are less exogenous than the paper's stress regimes",
    subtitle = "Points use count >=10 events; vertical ranges span count cutoffs 3, 5, and 10. The 72 h direct-capture sample is confounded by co-culture cross-contamination."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.5),
    plot.margin = margin(5.5, 32, 5.5, 5.5)
  )
ggsave(file.path(OUT, "barnyard_exo_fraction.png"), p, width = 10, height = 6, dpi = 130)

cat("=== realistic exogenous fraction from the barnyard (GEX-pure cells; doublets excluded) ===\n")
cat("  exo_pct_ge10 = clearest index-hop signal; ge3 = swap + Poisson soup-tail (upper bound)\n\n")
print(df, row.names = FALSE)
cat("\nRealistic exo (non-co-cultured mix0hr) ~ 0% (CROP-seq) to ~10-15% (direct-capture) vs paper's 25% / 75%.\n")
cat("mix72hr (72h co-culture) is confounded by cross-contamination -- not a clean index-hop read.\n")
cat("wrote", file.path(OUT, "barnyard_exo_fraction.csv"), "and barnyard_exo_fraction.png\n")
