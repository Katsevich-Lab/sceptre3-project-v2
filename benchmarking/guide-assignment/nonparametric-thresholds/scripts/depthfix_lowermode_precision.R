## ============================================================================
## The downstream cost of the detection floor, on WELL-SEPARATED (clean-gap) guides.
## Runs the real depth_fix once, then produces three linked outputs:
##  (1) lowermode_fraction.csv  -- for every clean-gap guide, what fraction of the
##      cells depth_fix calls fall in the below-gap LOWER MODE (counts 2..gap_lo)
##      rather than the signal mode (>=gap_hi). ~19% pooled; up to 92% for some.
##  (2) example_wellsep_assignments.png -- count histograms of 5 example guides,
##      coloured by how many cells in each bin depth_fix called.
##  (3) example_effectsizes.{csv,png} -- on-target knockdown (target expr / NT
##      baseline) estimated two ways: gap-based (signal mode only) vs ALL depth_fix
##      calls. Including the lower-mode calls dilutes the effect in proportion to
##      the lower-mode fraction (catastrophically for TOP2A: 0.05 -> 0.82).
## Run from nonparametric-thresholds/.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats); library(ondisc); library(sceptre)})
source("scripts/contingency_method.R")
OUT <- "results/global_ambient_poisson"
Dm <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")

mc   <- as(readRDS(Dm), "CsparseMatrix")
resp <- initialize_odm_from_backing_file(file.path(Dp,"response.odm"))
so   <- readRDS(file.path(Dp,"sceptre_object.rds")); lib <- exp(so@covariate_matrix[,"log(response_n_umis)"])
tdf  <- so@grna_target_data_frame
nt   <- tdf$grna_id[tdf$grna_target=="non-targeting"]; ntmask <- rep(FALSE, ncol(mc))
for(g in nt) ntmask <- ntmask | (as.numeric(mc[g,]) >= 30)

cat("running depth_fix (min_count=2, ambient/hyper/GS) once...\n")
A <- as(contingency_assign(mc, q=0.05, refit=10, min_count=2, cell_margin="ambient", tail="hyper", fdr="GS")$assigned, "CsparseMatrix")

## ---- (1) lower-mode fraction of assignments, per clean-gap guide ------------
pg <- read.csv(file.path(OUT,"perguide_replogle_rd7.csv"), stringsAsFactors=FALSE)
pg <- pg[pg$gap_lo>=2, ]; ridx <- match(pg$guide, rownames(mc))
lm <- data.frame()
for(i in seq_len(nrow(pg))){
  r<-ridx[i]; v<-as.numeric(mc[r,]); asg<-as.logical(A[r,]); n_tot<-sum(asg); if(!n_tot) next
  lm <- rbind(lm, data.frame(guide=pg$guide[i], gap_lo=pg$gap_lo[i], gap_hi=pg$gap_hi[i],
    n_assigned=n_tot, n_lower=sum(asg&v>=2&v<=pg$gap_lo[i]), n_count2=sum(asg&v==2),
    n_signal=sum(asg&v>=pg$gap_hi[i]),
    frac_lower=sum(asg&v>=2&v<=pg$gap_lo[i])/n_tot, frac_count2=sum(asg&v==2)/n_tot))
}
write.csv(lm, file.path(OUT,"lowermode_fraction.csv"), row.names=FALSE)
cat(sprintf("[lower-mode] %d clean-gap guides | pooled frac in lower mode = %.1f%% (count-2 %.1f%%)\n",
    nrow(lm), 100*sum(lm$n_lower)/sum(lm$n_assigned), 100*sum(lm$n_count2)/sum(lm$n_assigned)))
cat(sprintf("[lower-mode] per-guide frac_lower: median %.1f%% [IQR %.1f-%.1f%%] max %.1f%% ; >25%%: %d/%d guides\n",
    100*median(lm$frac_lower), 100*quantile(lm$frac_lower,.25), 100*quantile(lm$frac_lower,.75),
    100*max(lm$frac_lower), sum(lm$frac_lower>0.25), nrow(lm)))

## ---- 5 example guides spanning the range -----------------------------------
picks <- c("1290_CCDC6_P1P2_ENSG00000108091","7424_RPAP3_P1P2_ENSG00000005175",
           "968_C12orf60_P1P2_ENSG00000182993","9676_USP8_P1P2_ENSG00000138592",
           "9175_TOP2A_P1P2_ENSG00000131747")
gene <- function(g) sub("^[0-9]+_([^_]+)_.*","\\1",g)

## ---- (2) coloured count histograms -----------------------------------------
be <- function(maxc, thr=15, nlog=26){ b<-c(1:thr, round(exp(seq(log(thr+1),log(maxc+1),length.out=nlog)))); b<-sort(unique(b[b>=1&b<=maxc])); c(b,maxc+1) }
png(file.path(OUT,"example_wellsep_assignments.png"), width=1500, height=950, res=140)
par(mfrow=c(2,3), mar=c(4.2,4.3,3.6,1), mgp=c(2.4,0.7,0), cex.axis=0.95, cex.lab=1.05)
for(g in picks){
  r<-match(g,rownames(mc)); v<-as.numeric(mc[r,]); asg<-as.logical(A[r,]); info<-lm[lm$guide==g,]
  mx<-max(v); E<-be(mx); mb<-length(E)-1; lo<-E[1:mb]; hi<-E[2:(mb+1)]-1; FL<-0.5
  tot<-sapply(1:mb,function(i) sum(v>=lo[i]&v<=hi[i])); asn<-sapply(1:mb,function(i) sum(asg&v>=lo[i]&v<=hi[i]))
  plot(NA,log="xy",xlim=c(0.9,mx+1),ylim=c(FL,max(tot)*1.4),xlab="gRNA UMI count",ylab="number of cells",
       main=sprintf("%s   (gap %d-%d)\n%.0f%% of %d calls in lower mode",gene(g),info$gap_lo,info$gap_hi,100*info$frac_lower,info$n_assigned))
  rect(info$gap_lo+0.5,FL,info$gap_hi-0.5,max(tot)*1.4,col=rgb(1,0.6,0,0.09),border=NA)
  rect(lo,FL,hi+1,pmax(tot,FL),col="grey84",border="grey60",lwd=0.5)
  rect(lo,FL,hi+1,pmax(asn,FL),col=rgb(0.78,0.1,0.1,0.85),border=NA)
  abline(v=c(info$gap_lo+0.5,info$gap_hi-0.5),col="darkorange3",lty=3,lwd=0.8)
}
plot.new()
legend("topleft",bty="n",cex=1.15,fill=c("grey84",rgb(0.78,0.1,0.1,0.85),rgb(1,0.6,0,0.3)),border=c("grey60",NA,NA),
       legend=c("all cells with this count","called 'contains gRNA' (depth_fix)","the empty gap"))
text(0,0.5,adj=0,cex=1.0,labels=paste0("Well-separated Replogle guides. A clean guide (CCDC6)\n",
 "colours almost only the signal mode; as the lower-mode\n",
 "fraction rises (RPAP3->C12orf60->USP8->TOP2A), depth_fix\n",
 "increasingly calls below-gap cells -- mostly background\n","(~85-90%% no knockdown). TOP2A: 92%% of calls are lower mode."))
dev.off()

## ---- (3) effect sizes: gap-based vs all depth_fix calls ---------------------
set.seed(1); boot <- function(x,base) quantile(replicate(1000, mean(sample(x,length(x),TRUE))/base), c(.025,.975))
es <- data.frame()
for(g in picks){
  r<-match(g,rownames(mc)); tgt<-sub(".*_(ENSG[0-9]+)$","\\1",g)
  v<-as.numeric(mc[r,]); asg<-as.logical(A[r,]); cp<-as.numeric(resp[tgt,])/lib*1e4
  base<-mean(cp[ntmask & v==0]); gh<-lm$gap_hi[lm$guide==g]
  gc<-cp[v>=gh]; dc<-cp[asg]; eg<-mean(gc)/base; ed<-mean(dc)/base; cig<-boot(gc,base); cid<-boot(dc,base)
  es<-rbind(es,data.frame(gene=gene(g),frac_lower=lm$frac_lower[lm$guide==g],base_cp10k=round(base,2),
    n_gap=sum(v>=gh),n_depthfix=sum(asg),eff_gap=round(eg,3),eff_gap_lo=round(cig[1],3),eff_gap_hi=round(cig[2],3),
    eff_df=round(ed,3),eff_df_lo=round(cid[1],3),eff_df_hi=round(cid[2],3),effect_lost=round((ed-eg)/(1-eg),2)))
}
es <- es[order(es$frac_lower),]; write.csv(es, file.path(OUT,"example_effectsizes.csv"), row.names=FALSE)
cat("\n=== on-target effect: gap-based vs all depth_fix calls (expr / NT baseline; lower=stronger KD) ===\n")
print(es[,c("gene","frac_lower","base_cp10k","n_gap","n_depthfix","eff_gap","eff_df","effect_lost")], row.names=FALSE)

png(file.path(OUT,"example_effectsizes.png"), width=1150, height=680, res=140)
par(mar=c(5,4.8,3.6,1),mgp=c(2.9,0.8,0),cex.axis=0.95,cex.lab=1.1)
nD<-nrow(es); x<-seq_len(nD); w<-0.19
plot(NA,xlim=c(.5,nD+.5),ylim=c(0,1.15),xaxt="n",xlab="",ylab="target expression / NT baseline",
     main="On-target effect size: gap-based threshold vs all depth_fix calls\n(1.0 = no knockdown; lower bar = stronger knockdown)")
abline(h=1,col="grey45",lty=2)
rect(x-w-.02,0,x-.02,es$eff_gap,col="steelblue4",border=NA); rect(x+.02,0,x+w+.02,es$eff_df,col="firebrick",border=NA)
segments(x-w/2-.02,es$eff_gap_lo,x-w/2-.02,es$eff_gap_hi,lwd=1.4); segments(x+w/2+.02,es$eff_df_lo,x+w/2+.02,es$eff_df_hi,lwd=1.4)
text(x-w/2-.02,es$eff_gap,sprintf("%.2f",es$eff_gap),pos=3,cex=0.72,col="steelblue4")
text(x+w/2+.02,es$eff_df,sprintf("%.2f",es$eff_df),pos=3,cex=0.72,col="firebrick")
axis(1,at=x,labels=sprintf("%s\n(%.0f%% lower)",es$gene,100*es$frac_lower),padj=0.5,cex.axis=0.82)
legend("topleft",bty="n",cex=0.9,fill=c("steelblue4","firebrick"),
       legend=c("gap-based (signal mode only)","all depth_fix calls (adds lower mode)"))
dev.off()
cat("\nwrote lowermode_fraction.csv, example_wellsep_assignments.png, example_effectsizes.{csv,png}\n")
