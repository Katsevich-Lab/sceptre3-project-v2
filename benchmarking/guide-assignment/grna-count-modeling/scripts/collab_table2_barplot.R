#!/usr/bin/env Rscript
# Bar-graph of the barnyard "Table 2" (Liu et al. 2025) per-cell species-assignment accuracy,
# with our fishash+ row added. Reads results/collaborator_writeup/barnyard_table2_augmented.csv
# (fishash & fishash+ = our exact recomputation; other methods transcribed from the paper).
# Renders results/collaborator_writeup/barnyard_table2_barplot.png.
suppressPackageStartupMessages({library(ggplot2)})

OUT <- "results/collaborator_writeup"
aug <- read.csv(file.path(OUT, "barnyard_table2_augmented.csv"), comment.char="#",
                stringsAsFactors=FALSE)

# nice sample facet labels, ordered
samp_levels <- c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")
samp_labels <- c("CROP-seq, mix 0 hr","CROP-seq, mix 72 hr",
                 "direct-capture, mix 0 hr","direct-capture, mix 72 hr")
aug$sample <- factor(aug$sample, levels=samp_levels, labels=samp_labels)

# consistent method order across facets: by mean accuracy across the 4 samples (ascending
# so highest ends up at the TOP after coord_flip)
mean_acc <- tapply(aug$accuracy, aug$method, mean)
method_order <- names(sort(mean_acc))
aug$method <- factor(aug$method, levels=method_order)

# highlight colours: fishash+ accent, fishash secondary, rest grey
aug$hl <- ifelse(aug$method=="fishash+", "fishash+",
          ifelse(aug$method=="fishash",  "fishash", "other"))
aug$hl <- factor(aug$hl, levels=c("fishash+","fishash","other"))
pal <- c("fishash+"="#D55E00", "fishash"="#0072B2", "other"="grey70")

p <- ggplot(aug, aes(x=method, y=accuracy, fill=hl)) +
  geom_col(width=0.72) +
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=0.28,
                linewidth=0.35, colour="grey25") +
  facet_wrap(~sample, ncol=2) +
  scale_fill_manual(values=pal, name=NULL,
                    breaks=c("fishash+","fishash"),
                    labels=c("fishash+ (our method)","fishash (our recomputation)")) +
  coord_flip(ylim=c(0.40, 1.0)) +
  labs(
    x=NULL, y="Per-cell species-assignment accuracy"
  ) +
  theme_bw(base_size=16) +
  theme(
    legend.position="top",
    legend.margin=margin(b=-4),
    panel.grid.major.y=element_blank(),
    strip.background=element_rect(fill="grey92", colour=NA),
    strip.text=element_text(face="bold", size=14),
    axis.text.y=element_text(size=14)
  )

png_path <- file.path(OUT, "barnyard_table2_barplot.png")
ggsave(png_path, p, width=10, height=7.2, dpi=150)
cat("Wrote", png_path, "\n")
