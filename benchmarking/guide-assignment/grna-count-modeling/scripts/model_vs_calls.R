#!/usr/bin/env Rscript
# ============================================================================
# Two distinct objects per method, side by side, per guide: the rank-1 Poisson ambient MODEL
# (dashed = expected # ambient cells at each count, over the cells the method calls ambient) vs the
# actual ambient CALLS (solid = # cells actually left ambient). Self-consistent iff they coincide.
#
#   fishash   : rank-1 field from impute_masked_counts(counts, A_fishash)  [it computes this for the
#               guide-side Simpson correction, but CLASSIFIES the cell side against the raw library]
#   fishash+  : rank-1 field from impute_masked_counts(counts, A_fishash+) [and classifies against it]
#   CLEANSER  : rank-1 field on the <=2 (sub-threshold) matrix -- the model implied by its L_i depth
#               [it classifies via an NB mixture, decoupled from L_i]
# Usage:  Rscript scripts/model_vs_calls.R <highmoi|lowmoi>
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2); library(tidyr); library(dplyr) })
source("scripts/contingency_method.R"); source("scripts/datasets.R")
OUT <- "results/ambient_ceiling"
which <- commandArgs(trailingOnly=TRUE)[1]; if(is.na(which)) which<-"highmoi"
cfg <- list(
  highmoi=list(dir=dctap_sceptre_dir("highmoi"), ds="dctap_k562_highmoi", lab="DC-TAP high-MOI", tag="dctap_highmoi"),
  lowmoi =list(dir=dctap_sceptre_dir("lowmoi"),  ds="dctap_k562_lowmoi",  lab="DC-TAP low-MOI",  tag="dctap_lowmoi"))[[which]]

gm  <- load_grna_matrix(file.path(cfg$dir,"grna_matrix_aligned.rds"))
tdf <- read.csv(file.path(cfg$dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
sym_of <- setNames(tdf$grna_target_symbol, tdf$grna_id)

# assignments
A_f <- as(fishash_assign(gm)$assigned,"lgCMatrix")
A_p <- as(fishashplus_assign(gm)$assigned,"lgCMatrix")
A_c <- as(load_cleanser_assignment(gm, cfg$ds),"lgCMatrix")

# Expected-ambient model per method, built from the depth each method ACTUALLY USES as the cell
# exposure: the guide's ambient total R_g spread across cells in proportion to that depth.
#   lambda_{g,c} = R_g * d_c / sum_c d_c    (so sum_c lambda = R_g, the guide's ambient total)
#   fishash   : d_c = raw library size            (its actual cell exposure)
#   fishash+  : d_c = denoised rank-1 depth
#   CLEANSER  : d_c = L_i = sum of <=2 counts
# R_g = guide ambient total from the corresponding signal-free matrix (rank-1 completion / <=2 matrix).
mk <- function(A, Rg, d) list(A=A, Rg=Rg, d=as.numeric(d), sumd=sum(d))
bgf <- fishash::impute_masked_counts(gm, A_f); bgp <- fishash::impute_masked_counts(gm, A_p)
sub <- gm; sub@x[sub@x>2] <- 0                                # CLEANSER's <=2 sub-threshold matrix
rawlib <- Matrix::colSums(gm)
FLD <- list(fishash   = mk(A_f, Matrix::rowSums(bgf), rawlib),
            `fishash+`= mk(A_p, Matrix::rowSums(bgp), Matrix::colSums(bgp)),
            CLEANSER  = mk(A_c, Matrix::rowSums(sub), Matrix::colSums(sub)))

sel <- intersect(c("GATA1-TSS-2","e-HDAC6-1","HDAC6-TSS-1","MYC-e2-1"), rownames(gm))
KMAX <- 60L; ks <- 0:KMAX
rows_bar<-list(); rows_line<-list(); gof<-list()
for (g in sel) {
  y <- as.numeric(gm[g,]); gidx<-match(g,rownames(gm))
  obs <- vapply(ks, function(k) if(k<KMAX) sum(y==k) else sum(y>=k), numeric(1))
  rows_bar[[g]] <- data.frame(guide=g, gene=sym_of[g], k=ks, observed=obs)
  for (mn in names(FLD)) {
    f<-FLD[[mn]]; lam <- (f$Rg[gidx]/f$sumd)*f$d
    amb <- !as.logical(f$A[g,]); amb[is.na(amb)]<-TRUE
    model <- vapply(ks, function(k) if(k<KMAX) sum(dpois(k,lam[amb])) else sum(1-ppois(k-1,lam[amb])), numeric(1))
    cal   <- vapply(ks, function(k) if(k<KMAX) sum(y[amb]==k) else sum(y[amb]>=k), numeric(1))
    rows_line[[paste(g,mn)]] <- rbind(
      data.frame(guide=g,gene=sym_of[g],method=mn,k=ks,object="MODEL (rank-1 Poisson)",n=model),
      data.frame(guide=g,gene=sym_of[g],method=mn,k=ks,object="CALLS (actual assignment)",n=cal))
    gof[[paste(g,mn)]] <- data.frame(dataset=cfg$tag, guide=g, gene=sym_of[g], method=mn,
      n_ambient=sum(cal), L1_frac=sum(abs(model-cal))/sum(cal))
  }
}
BAR<-bind_rows(rows_bar); LIN<-bind_rows(rows_line); GOF<-bind_rows(gof)
lev<-sprintf("%s (%s)", sel, sym_of[sel])
mlev<-c("fishash","fishash+","CLEANSER")
BAR$facet<-factor(sprintf("%s (%s)",BAR$guide,BAR$gene),levels=lev)
LIN$facet<-factor(sprintf("%s (%s)",LIN$guide,LIN$gene),levels=lev); LIN$method<-factor(LIN$method,levels=mlev)
BARg<-BAR[rep(seq_len(nrow(BAR)),3),]; BARg$method<-factor(rep(mlev,each=nrow(BAR)),levels=mlev)

p <- ggplot(BARg, aes(k, observed)) +
  geom_col(fill="grey82", width=0.9) +
  geom_line(data=LIN, aes(k, n, colour=method, linetype=object), linewidth=0.7) +
  facet_grid(method ~ facet) +
  scale_colour_manual(values=c(fishash="#7f8c8d",`fishash+`="#c0392b",CLEANSER="#2980b9"), guide="none") +
  scale_linetype_manual(values=c("MODEL (rank-1 Poisson)"=2,"CALLS (actual assignment)"=1), name=NULL) +
  scale_y_continuous(trans=scales::trans_new("log1p",log1p,expm1), breaks=c(0,1,10,100,1000),
                     labels=scales::label_comma()) +
  coord_cartesian(xlim=c(0,KMAX)) +
  labs(title=sprintf("Expected ambient (from depth USED) vs observed, per guide (%s)", cfg$lab),
       subtitle="Grey bars = observed cells/count. Dashed = expected # ambient cells at each count from each method's depth USED (fishash=raw library, fishash+=rank-1, CLEANSER=≤2), over its called-ambient cells; solid = # cells it actually leaves ambient. Observed above dashed = signal the method should call.",
       x="gRNA UMI count", y="number of cells (log1p; floor 0)") +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave(file.path(OUT,sprintf("model_vs_calls_%s.png",cfg$tag)), p, width=13, height=8.5, dpi=140)

## overlay version: one panel per guide, the three methods' EXPECTED AMBIENT curves on the histogram
ov <- LIN %>% filter(object=="MODEL (rank-1 Poisson)")
pE <- ggplot(BAR, aes(k, observed)) +
  geom_col(fill="grey82", width=0.9) +
  geom_line(data=ov, aes(k, n, colour=method), linewidth=0.8) +
  facet_wrap(~facet, ncol=2, scales="free_y") +
  scale_colour_manual(values=c(fishash="#7f8c8d",`fishash+`="#c0392b",CLEANSER="#2980b9"),
                      name="expected ambient (from depth used):") +
  scale_y_continuous(trans=scales::trans_new("log1p",log1p,expm1), breaks=c(0,1,10,100,1000),
                     labels=scales::label_comma()) +
  coord_cartesian(xlim=c(0,KMAX)) +
  labs(title=sprintf("Expected ambient count per bar, by depth used (%s)", cfg$lab),
       subtitle="Grey bars = observed cells/count. Lines = expected # ambient cells at each count implied by each method's depth used\n(fishash=raw library, fishash+=rank-1, CLEANSER=≤2). Observed ABOVE a line = signal that method should call; a line that stays high (fishash) expects ambient at high counts and self-masks real signal.",
       x="gRNA UMI count", y="number of cells (log1p; floor 0)") +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave(file.path(OUT,sprintf("expected_ambient_%s.png",cfg$tag)), pE, width=12, height=8, dpi=140)

write.csv(GOF, file.path(OUT,sprintf("model_vs_calls_%s_gof.csv",cfg$tag)), row.names=FALSE)
cat(sprintf("=== %s : MODEL vs CALLS L1 discrepancy (0 = self-consistent) ===\n", cfg$lab))
print(GOF %>% group_by(method) %>% summarize(median_L1=median(L1_frac), .groups="drop") %>% as.data.frame(), row.names=FALSE, digits=3)
cat("\nwrote", file.path(OUT,sprintf("model_vs_calls_%s.png",cfg$tag)),"\n")
