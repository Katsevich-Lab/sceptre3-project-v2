#!/usr/bin/env Rscript
# ============================================================================
# method_comparison_2.qmd — Figure 1 (per-dataset ambient-noise landscape) + the shared
# per-dataset table the other mc2_ figure scripts read for capture-modality colouring.
#
# Run from the grna-count-modeling/ folder root:  Rscript scripts/mc2_dataset_landscape.R
#
# INPUTS  (produced upstream by scripts/writeup_compute.R -> results/ambient_ceiling/writeup/):
#   depth_summary.csv  ceilings.csv  assign_counts.csv  jaccard_ff.csv  depths_sampled.csv
# OUTPUTS (results/ambient_ceiling/):
#   dataset_landscape.csv         per-dataset metrics + authoritative capture modality (shared input)
#   dataset_landscape.png         Figure 1 (scatter: median ambient depth vs median ambient ceiling)
#   dataset_landscape_cross.png   same, with 5-95% variability crosses (not embedded in the writeup)
# Everything is fishash+ (= depth_fix); no reassignment happens here — it reads cached summaries.
# ============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(ggrepel)})
W <- "results/ambient_ceiling/writeup"

depth  <- read.csv(file.path(W, "depth_summary.csv"))
ceil   <- read.csv(file.path(W, "ceilings.csv"))
acount <- read.csv(file.path(W, "assign_counts.csv"))
jacc   <- read.csv(file.path(W, "jaccard_ff.csv"))

# per-guide ambient ceiling -> per-dataset summaries (median / p90 / max)
ceil_ds <- ceil %>% group_by(dataset) %>% summarise(
  med_ceil_fplus = median(ceil_fishashplus, na.rm=TRUE),
  p90_ceil_fplus = quantile(ceil_fishashplus, 0.90, na.rm=TRUE),
  max_ceil_fplus = max(ceil_fishashplus, na.rm=TRUE),
  med_ceil_fish  = median(ceil_fishash, na.rm=TRUE), .groups="drop")

df <- depth %>%
  left_join(ceil_ds, by="dataset") %>%
  left_join(acount %>% select(dataset, n_fplus=fishashplus, n_fish=fishash,
                              moi_fplus=moi_fishashplus, moi_fish=moi_fishash), by="dataset") %>%
  left_join(jacc %>% select(dataset, jaccard_ff=jaccard), by="dataset") %>%
  mutate(decontam_ratio = round(med_lib / med_depth_fishashplus, 2),
         assign_per_cell = round(n_fplus / cells, 3),
         med_depth_fishashplus = round(med_depth_fishashplus, 2),
         med_depth_cleanser = round(med_depth_cleanser, 2))

# authoritative gRNA capture chemistry (GEO records + methods + CLEANSER paper + build_inputs.py)
chem <- c(
  gasperini="CROP-seq", gastric_organoid="CROP-seq", a549="CROP-seq",
  barnyard_lrb100_0hr="CROP-seq", barnyard_lrb100_72hr="CROP-seq",
  `replogle-rd7`="direct capture", mccutcheon="direct capture", cd8_tcell="direct capture",
  dctap_k562_highmoi="direct capture", dctap_k562_lowmoi="direct capture",
  endoc_t2d="direct capture", invivo_cortex="direct capture", ipsc="direct capture",
  barnyard_mch2_0hr="direct capture", barnyard_mch2_72hr="direct capture",
  erythroid_multiome="other", cd4_tcell="other")
conf <- c(ipsc="med-high", erythroid_multiome="med", cd4_tcell="med")   # rest: high
df$modality  <- chem[df$dataset]
df$chem_conf <- ifelse(df$dataset %in% names(conf), conf[df$dataset], "high")
write.csv(df, "results/ambient_ceiling/dataset_landscape.csv", row.names=FALSE)

# 5-95% bands for the cross variant
dep <- read.csv(file.path(W, "depths_sampled.csv"))
xb <- dep %>% filter(fishashplus > 0) %>% group_by(dataset) %>%
  summarise(x5=quantile(fishashplus,.05), x95=quantile(fishashplus,.95), .groups="drop")
yb <- ceil %>% filter(!is.na(ceil_fishashplus)) %>% group_by(dataset) %>%
  summarise(y5=quantile(ceil_fishashplus,.05), y95=quantile(ceil_fishashplus,.95), .groups="drop")
df <- df %>% left_join(xb, by="dataset") %>% left_join(yb, by="dataset")
df$x <- df$med_depth_fishashplus; df$y <- df$med_ceil_fplus
pal <- c("CROP-seq"="#eb6834","direct capture"="#2a78d6","other"="#1d9e75")
df$modality <- factor(df$modality, levels=names(pal))

mk <- function(cross) {
  p <- ggplot(df, aes(x, y, color=modality))
  if (cross) p <- p +
    geom_errorbar(aes(xmin=x5, xmax=x95, y=y), orientation="y", width=0, linewidth=0.5, alpha=0.5) +
    geom_errorbar(aes(ymin=y5, ymax=y95, x=x), width=0, linewidth=0.5, alpha=0.5)
  xlab <- if (cross) "median ambient gRNA depth across cells  (fishash+)  |  arms = 5-95% across cells"
          else       "median ambient gRNA depth across cells  (fishash+)"
  ylab <- if (cross) "median largest ambient UMI count over guides  (fishash+)  |  arms = 5-95% across guides"
          else       "median over guides of largest ambient UMI count  (fishash+)"
  ttl  <- if (cross) "Per-dataset ambient-noise landscape (fishash+), with 5-95% variability crosses"
          else       "Per-dataset ambient-noise landscape (fishash+)"
  p + geom_point(size=2.7) +
    geom_text_repel(aes(label=dataset), size=2.9, max.overlaps=Inf,
                    min.segment.length=0, box.padding=0.5, show.legend=FALSE) +
    scale_x_log10() + scale_y_log10() + scale_color_manual(values=pal, drop=FALSE) +
    labs(x=xlab, y=ylab, color="capture modality", title=ttl) +
    theme_bw(base_size=11) +
    theme(legend.position=c(0.99,0.02), legend.justification=c(1,0),
          legend.background=element_rect(fill=alpha("white",0.7), color=NA),
          plot.title=element_text(size=11))
}
ggsave("results/ambient_ceiling/dataset_landscape.png",       mk(FALSE), width=8,   height=6,   dpi=150)
ggsave("results/ambient_ceiling/dataset_landscape_cross.png", mk(TRUE),  width=8.5, height=6.3, dpi=150)
cat("wrote dataset_landscape.csv + dataset_landscape.png + dataset_landscape_cross.png\n")
