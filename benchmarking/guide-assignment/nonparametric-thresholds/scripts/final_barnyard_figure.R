#!/usr/bin/env Rscript
# Consolidated barnyard comparison: the 10 methods from Fishash Table 2
# (transcribed) + Otsu, smoothed valley, and the ambient test (computed here,
# q=0.01). Metric = Table-2 accuracy (>=1 correct-species guide & 0 wrong).
suppressPackageStartupMessages({library(ggplot2); library(tidyr); library(dplyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
cols <- c("Cropseq mix0hr", "Cropseq mix72hr", "DirectCapture mix0hr", "DirectCapture mix72hr")
pub <- tribble(~method, ~`Cropseq mix0hr`, ~`Cropseq mix72hr`, ~`DirectCapture mix0hr`, ~`DirectCapture mix72hr`,
  "cleanser cs",0.9414,0.9410,0.8612,0.9221,  "cleanser dc",0.9178,0.8459,0.8614,0.9308,
  "crispat gauss",0.9367,0.9346,0.8200,0.8945, "crispat negbinom",0.4556,0.4815,0.8865,0.9381,
  "crispat poisgauss",0.9501,0.9502,0.8739,0.9331, "crispat poisson",0.9086,0.9169,0.8763,0.8956,
  "demuxem",0.9364,0.9364,0.8975,0.9408, "fishash",0.9563,0.9533,0.9100,0.9491,
  "geomux",0.9578,0.9555,0.7471,0.5252, "sceptre mixture",0.9491,0.9472,0.8953,0.9408)
mine <- tribble(~method, ~`Cropseq mix0hr`, ~`Cropseq mix72hr`, ~`DirectCapture mix0hr`, ~`DirectCapture mix72hr`,
  "Otsu",0.8368,0.9064,0.9432,0.9453, "smoothed valley",0.8668,0.9250,0.9446,0.9417,
  "AMBIENT TEST (q=.01)",0.9394,0.9655,0.9600,0.9562)
df <- bind_rows(pub %>% mutate(grp="published"), mine %>% mutate(grp=ifelse(grepl("AMBIENT",method),"ambient","nonparam"))) %>%
  pivot_longer(all_of(cols), names_to="sample", values_to="acc")
ord <- df %>% group_by(method) %>% summarize(m=mean(acc)) %>% arrange(m) %>% pull(method)
df$method <- factor(df$method, levels=ord)
p <- ggplot(df, aes(acc, method, color=grp)) +
  geom_point(size=2) + facet_wrap(~sample, nrow=1) +
  scale_color_manual(values=c(ambient="#d7301f", nonparam="#2c7fb8", published="grey55"), guide="none") +
  labs(title="Barnyard species-assignment accuracy: ambient test (red) vs published methods (grey) + Otsu/valley (blue)",
       subtitle="metric = fraction of cells with >=1 correct-species guide and 0 wrong-species guides (Fishash Table 2)",
       x="accuracy", y=NULL) +
  theme_bw(base_size=9) + theme(panel.grid.minor=element_blank())
ggsave(file.path(HERE,"results","Final_Barnyard_AllMethods.png"), p, width=12, height=4.2, dpi=120)
# mean-across-samples ranking
rank <- df %>% group_by(method,grp) %>% summarize(mean_acc=round(mean(acc),4),.groups="drop") %>% arrange(-mean_acc)
cat("Mean accuracy across the 4 barnyard samples (ranked):\n"); print(as.data.frame(rank), row.names=FALSE)
write.csv(rank, file.path(HERE,"results","Final_Barnyard_ranking.csv"), row.names=FALSE)
cat("\nWrote results/Final_Barnyard_AllMethods.png\n")
