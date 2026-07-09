#!/usr/bin/env Rscript
# Per-guide assignment concordance (Jaccard of assigned cells) for two comparisons, per dataset:
#   fishash vs fishash+   and   CLEANSER vs fishash+ .  CLEANSER uses chemistry-correct calls
#   (--cs where available, else --dc). Non-barnyard, clean-chemistry datasets only.
suppressMessages({library(Matrix); library(ggplot2); library(dplyr)})
source("scripts/datasets.R")
DS <- c("dctap_k562_highmoi","dctap_k562_lowmoi","endoc_t2d","invivo_cortex",
        "mccutcheon","a549","cd8_tcell","gastric_organoid")     # depth desc
ADIR <- "results/ambient_ceiling/writeup/assign"

load_cl <- function(gm, ds) {                                   # prefer --cs (chemistry-correct)
  cf <- c(file.path("results/ambient_ceiling/cleanser_calls_cs", paste0(ds,".csv")),
          file.path("results/ambient_ceiling/cleanser_calls",    paste0(ds,".csv")))
  cf <- cf[file.exists(cf)][1]; if (is.na(cf)) return(NULL)
  cl <- read.csv(cf); gg <- match(cl$guide,rownames(gm)); cc <- match(cl$cell,colnames(gm)); ok <- !is.na(gg)&!is.na(cc)
  as(sparseMatrix(i=gg[ok], j=cc[ok], x=TRUE, dims=dim(gm), dimnames=dimnames(gm)), "lgCMatrix")
}
jac_pg <- function(A, B, minu=3) { inter <- rowSums(A & B); uni <- rowSums(A | B)
  j <- inter/uni; j[uni < minu] <- NA; j }

rows <- list(); summ <- list()
for (ds in DS) {
  gm <- load_grna_matrix(ds)
  Af <- as(readRDS(file.path(ADIR,paste0(ds,"_fishash.rds"))),    "lgCMatrix")
  Ap <- as(readRDS(file.path(ADIR,paste0(ds,"_fishashplus.rds"))),"lgCMatrix")
  Ac <- load_cl(gm, ds); if (is.null(Ac)) next
  jfp <- jac_pg(Af, Ap); jcp <- jac_pg(Ac, Ap)
  rows[[length(rows)+1]] <- data.frame(dataset=ds, comparison="fishash vs fishash+", jaccard=jfp)
  rows[[length(rows)+1]] <- data.frame(dataset=ds, comparison="CLEANSER vs fishash+", jaccard=jcp)
  summ[[length(summ)+1]] <- data.frame(dataset=ds,
    med_fp=median(jfp,na.rm=TRUE), med_cp=median(jcp,na.rm=TRUE),
    n_guides=sum(!is.na(jfp)))
}
df <- bind_rows(rows) %>% filter(!is.na(jaccard))
df$dataset    <- factor(df$dataset, levels=DS)
df$comparison <- factor(df$comparison, levels=c("fishash vs fishash+","CLEANSER vs fishash+"))
pal <- c("fishash vs fishash+"="#2a78d6","CLEANSER vs fishash+"="#e34948")

dg <- position_dodge(0.8)
p <- ggplot(df, aes(dataset, jaccard, fill=comparison, color=comparison)) +
  geom_violin(position=dg, width=0.8, alpha=0.15, linewidth=0.2, scale="width") +
  geom_boxplot(position=dg, width=0.2, alpha=0.28, outlier.size=0.3, linewidth=0.22, show.legend=FALSE) +
  stat_summary(fun=median, fun.min=median, fun.max=median, geom="crossbar",   # bold median marker
               position=dg, width=0.34, linewidth=0.5, color="grey12", fill=NA, show.legend=FALSE) +
  stat_summary(fun=median, geom="text", aes(label=sprintf("%.2f", after_stat(y))),
               position=dg, size=2.7, color="grey12", vjust=1.7, show.legend=FALSE) +
  scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
  scale_y_continuous(limits=c(0,1.02), breaks=seq(0,1,0.25)) +
  labs(x=NULL, y="per-guide assignment concordance (Jaccard of assigned cells)",
       fill=NULL, color=NULL,
       title="Per-guide method concordance across datasets  (bold bar = median)") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", axis.text.x=element_text(angle=25, hjust=1),
        plot.title=element_text(size=12))
ggsave("results/ambient_ceiling/method_concordance.png", p, width=10.5, height=5.6, dpi=150)

st <- bind_rows(summ)
write.csv(st, "results/ambient_ceiling/method_concordance_summary.csv", row.names=FALSE)
cat("\n== median per-guide Jaccard ==\n"); print(st, row.names=FALSE)
cat("\nwrote results/ambient_ceiling/method_concordance.png\n")
