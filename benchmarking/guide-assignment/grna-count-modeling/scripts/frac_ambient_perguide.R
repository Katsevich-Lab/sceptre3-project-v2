#!/usr/bin/env Rscript
# ============================================================================
# Per-guide, per gRNA-count bin: the FRACTION of cells at that count each method leaves AMBIENT
# (unassigned). Purely empirical -- no model, no depth. Grey bars show the count distribution
# (scaled) for context; the three lines are the ambient fraction for fishash / fishash+ / CLEANSER.
# At low counts everything is ambient (fraction ~1); as counts rise each method starts assigning
# (fraction falls). A line that stays high = the method keeps calling that count ambient.
# Usage: Rscript scripts/frac_ambient_perguide.R <highmoi|lowmoi>
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2); library(dplyr); library(tidyr) })
source("scripts/contingency_method.R"); source("scripts/datasets.R")
OUT <- "results/ambient_ceiling"
which <- commandArgs(trailingOnly=TRUE)[1]; if(is.na(which)) which<-"highmoi"
cfg <- list(
  highmoi=list(dir=dctap_sceptre_dir("highmoi"), ds="dctap_k562_highmoi", lab="DC-TAP high-MOI", tag="dctap_highmoi"),
  lowmoi =list(dir=dctap_sceptre_dir("lowmoi"),  ds="dctap_k562_lowmoi",  lab="DC-TAP low-MOI",  tag="dctap_lowmoi"))[[which]]

gm  <- load_grna_matrix(file.path(cfg$dir,"grna_matrix_aligned.rds"))
tdf <- read.csv(file.path(cfg$dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
sym_of <- setNames(tdf$grna_target_symbol, tdf$grna_id)
pos <- tdf$grna_id[tdf$grna_target_symbol!="non-targeting"]

A_f <- as(fishash_assign(gm)$assigned,"lgCMatrix")
A_p <- as(fishashplus_assign(gm)$assigned,"lgCMatrix")
A_c <- as(load_cleanser_assignment(gm, cfg$ds),"lgCMatrix")
methods <- list(fishash=A_f, `fishash+`=A_p, CLEANSER=A_c)

KMAX <- 50L; ks <- 0:KMAX
bar<-list(); lin<-list()
for (g in pos) {
  y <- as.numeric(gm[g,])
  n <- vapply(ks, function(k) if(k<KMAX) sum(y==k) else sum(y>=k), numeric(1))
  bar[[g]] <- data.frame(guide=g, gene=sym_of[g], k=ks, n=n)
  for (mn in names(methods)) {
    amb <- !as.logical(methods[[mn]][g,]); amb[is.na(amb)]<-TRUE
    na <- vapply(ks, function(k) if(k<KMAX) sum(y==k & amb) else sum(y>=k & amb), numeric(1))
    lin[[paste(g,mn)]] <- data.frame(guide=g, gene=sym_of[g], method=mn, k=ks, n_ambient=na)
  }
}
BAR<-bind_rows(bar); LIN<-bind_rows(lin)
lev<-sprintf("%s (%s)", pos, sym_of[pos])
BAR$facet<-factor(sprintf("%s (%s)",BAR$guide,BAR$gene),levels=lev)
LIN$facet<-factor(sprintf("%s (%s)",LIN$guide,LIN$gene),levels=lev)
LIN$method<-factor(LIN$method,levels=c("fishash","fishash+","CLEANSER"))

p <- ggplot() +
  geom_col(data=BAR, aes(k, n), fill="grey82", width=0.9) +
  geom_line(data=LIN, aes(k, n_ambient, colour=method), linewidth=0.8) +
  facet_wrap(~facet, ncol=4) +
  scale_colour_manual(values=c(fishash="#7f8c8d",`fishash+`="#c0392b",CLEANSER="#2980b9"), name=NULL) +
  scale_y_continuous(trans=scales::trans_new("log1p",log1p,expm1), breaks=c(0,1,10,100,1000),
                     labels=scales::label_comma()) +
  coord_cartesian(xlim=c(0,KMAX)) +
  labs(title=sprintf("Number of cells at each gRNA count called ambient (%s)", cfg$lab),
       subtitle="Grey bars = observed cells per count. Lines = number of those cells each method leaves ambient (unassigned). Where a line tracks the bars, the method calls that count ambient; where it drops below, it assigns as signal.",
       x="gRNA UMI count", y="number of cells (log1p; floor 0)") +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave(file.path(OUT,sprintf("num_ambient_%s.png",cfg$tag)), p, width=13, height=8, dpi=140)
cat("wrote", file.path(OUT,sprintf("num_ambient_%s.png",cfg$tag)), "\n")
