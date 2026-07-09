#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr)})
have_pw <- requireNamespace("patchwork", quietly=TRUE)

land <- read.csv("results/ambient_ceiling/dataset_landscape.csv")   # has modality
mod  <- setNames(land$modality, land$dataset)

qs <- function(v){ v<-v[!is.na(v)]; q<-quantile(v,c(.05,.25,.5,.75,.95))
  data.frame(ymin=q[1],lower=q[2],middle=q[3],upper=q[4],ymax=q[5]) }

# per-cell fishash+ ambient depth (positive)
dep <- read.csv("results/ambient_ceiling/writeup/depths_sampled.csv") %>% filter(fishashplus>0)
sx  <- dep %>% group_by(dataset) %>% reframe(qs(fishashplus))
# per-guide fishash+ ambient ceiling
cei <- read.csv("results/ambient_ceiling/writeup/ceilings.csv")
sy  <- cei %>% group_by(dataset) %>% reframe(qs(ceil_fishashplus))

# shared row order: decreasing by median ambient depth (deepest at top)
ord <- sy %>% left_join(sx, by="dataset", suffix=c("_c","_d")) %>%
  arrange(middle_d) %>% pull(dataset)
for (d in list("sx","sy")) NULL
sx$dataset <- factor(sx$dataset, levels=ord); sy$dataset <- factor(sy$dataset, levels=ord)
sx$modality <- factor(mod[as.character(sx$dataset)], levels=c("CROP-seq","direct capture","other"))
sy$modality <- factor(mod[as.character(sy$dataset)], levels=c("CROP-seq","direct capture","other"))
pal <- c("CROP-seq"="#eb6834","direct capture"="#2a78d6","other"="#1d9e75")

strip <- function(s, xlab, showy, brks=waiver(), logscale=TRUE) {
  sc <- if (logscale) scale_y_log10(breaks=brks) else scale_y_continuous(breaks=brks)
  ggplot(s, aes(x=dataset, fill=modality, color=modality)) +
    geom_boxplot(aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax),
                 stat="identity", width=0.6, alpha=0.35, linewidth=0.4) +
    coord_flip() + sc +
    scale_fill_manual(values=pal, drop=FALSE) + scale_color_manual(values=pal, drop=FALSE) +
    labs(x=NULL, y=xlab) + theme_bw(base_size=11) +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          axis.text.y=if(showy) element_text() else element_blank(),
          legend.position="none")
}
pL <- strip(sx, "ambient depth per cell",   TRUE)
pR <- strip(sy, "ambient ceiling per guide", FALSE, brks=c(1,2,3,4,5,10,30))

if (have_pw) {
  library(patchwork)
  g <- (pL | pR) + plot_layout(widths=c(1,1), guides="collect") +
    plot_annotation(title="Per-dataset ambient variability (fishash+)",
      theme=theme(plot.title=element_text(size=12, hjust=0.5))) &
    theme(legend.position="bottom")
  ggsave("results/ambient_ceiling/dataset_variability_strips.png", g, width=10, height=6, dpi=150)
  cat("wrote dataset_variability_strips.png (patchwork)\n")
} else {
  ggsave("results/ambient_ceiling/dataset_variability_ceiling.png", pL, width=6, height=6, dpi=150)
  ggsave("results/ambient_ceiling/dataset_variability_depth.png",   pR, width=6, height=6, dpi=150)
  cat("wrote two PNGs (no patchwork)\n")
}
