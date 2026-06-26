#!/usr/bin/env Rscript
# Score the geomux PR sweep, combine with ours/fishash/otsu/valley, and plot
# precision-recall curves by difficulty regime.
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
EXP  <- file.path(HERE, "results", "comprehensive_bench")

gm <- readRDS(file.path(EXP,"gm.rds")); truth <- readRDS(file.path(EXP,"truth.rds"))
guides <- read.csv(file.path(EXP,"guides.csv")); Ncell <- ncol(gm)
mu <- as.integer(sub("mu","", vapply(strsplit(guides$group,"_"), `[`, character(1), 2)))
regime <- ifelse(mu <= 15, "hard (mu<=15)", ifelse(mu >= 150, "easy (mu>=150)", "moderate (mu=50)"))
tt <- as(truth,"TsparseMatrix"); gg <- as(gm,"TsparseMatrix")
truth_key <- intersect(tt@i*Ncell + tt@j, gg@i*Ncell + gg@j)
truth_g <- truth_key %/% Ncell + 1L
pr_one <- function(key, g) do.call(rbind, lapply(unique(regime), function(rg){
  gid <- which(regime==rg); A <- key[g %in% gid]; T <- truth_key[truth_g %in% gid]
  tp <- length(intersect(A,T))
  data.frame(regime=rg, precision=if(length(A))tp/length(A) else NA, recall=if(length(T))tp/length(T) else NA)}))

rows <- list()
for (f in list.files(file.path(EXP,"geomux_pr"), full.names=TRUE)) {
  d <- read.csv(f); gi <- match(d$guide, rownames(gm)); ci <- match(d$cell_barcode, colnames(gm))
  ok <- !is.na(gi)&!is.na(ci); pr <- pr_one((gi[ok]-1L)*Ncell+(ci[ok]-1L), gi[ok])
  pr$method <- "geomux"; pr$cutoff <- as.numeric(sub("q","",sub(".csv","",basename(f)))); rows[[length(rows)+1]] <- pr
}
geo <- do.call(rbind, rows)
cols <- c("regime","precision","recall","method","cutoff")
df <- rbind(read.csv(file.path(HERE,"results","pr_curves.csv"))[,cols], geo[,cols])
df$regime <- factor(df$regime, levels=c("hard (mu<=15)","moderate (mu=50)","easy (mu>=150)"))
write.csv(df, file.path(HERE,"results","pr_curves_all.csv"), row.names=FALSE)

curves <- df %>% filter(method %in% c("ours (ambient)","fishash","geomux")) %>% arrange(method, regime, cutoff)
pts    <- df %>% filter(method %in% c("otsu","valley"))
p <- ggplot(curves, aes(recall, precision, color=method)) +
  geom_path() + geom_point(size=1.3) +
  geom_point(data=pts, aes(recall, precision, color=method), shape=8, size=2.6) +
  facet_wrap(~regime, nrow=1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  scale_color_manual(values=c(`ours (ambient)`="#238b45", fishash="#d7301f", geomux="#7a0177", otsu="#2c7fb8", valley="#fb9a29")) +
  labs(title="Precision-recall by difficulty regime (comprehensive sim, full ground truth)",
       subtitle="curves = FDR-knob sweep (ours BH q; fishash padj; geomux p-value); stars = Otsu/valley operating points",
       x="recall", y="precision", color=NULL) + theme_bw(base_size=10) + theme(legend.position="top")
ggsave(file.path(HERE,"results","PR_curves.png"), p, width=11, height=4.2, dpi=120)

# matched-precision readout: recall at precision>=0.9 for each method/regime
cat("\n=== recall at highest cutoff with precision>=0.90 (per regime) ===\n")
df %>% filter(method %in% c("ours (ambient)","fishash","geomux"), precision>=0.90) %>%
  group_by(method, regime) %>% summarize(max_recall_at_prec90=round(max(recall),3), .groups="drop") %>%
  arrange(regime, -max_recall_at_prec90) %>% as.data.frame() %>% print(row.names=FALSE)
cat("\nWrote results/PR_curves.png + pr_curves_all.csv\n")
