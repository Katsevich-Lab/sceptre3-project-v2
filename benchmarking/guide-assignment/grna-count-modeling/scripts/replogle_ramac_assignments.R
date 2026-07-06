## ============================================================================
## What depth_fix actually ASSIGNS on RAMAC, and a precision cross-reference.
## Self-contained (reads the raw gRNA matrix; needs fishash for contingency_assign).
##  (1) plug-in Poisson carrier mask (p<1e-6)      -> ramac_labeled_cells.png
##  (2) real depth_fix assignment (hyper + GS FDR)  -> ramac_depthfix_assignments.png
##  (3) cross-reference: do RAMAC's low-count calls land in cells that clearly
##      belong to ANOTHER guide (doublet/spillover)? -> ramac_crossref.txt
## Run from the grna-count-modeling/ folder.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats)})
source("scripts/contingency_method.R")
source("scripts/datasets.R")   # load_grna_matrix()
OUT  <- "results/replogle_ambient_poisson"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
mc <- load_grna_matrix("replogle-rd7"); G <- nrow(mc); C <- ncol(mc)
RAMAC <- which(grepl("RAMAC", rownames(mc)))[1]; v <- as.numeric(mc[RAMAC,])

## --- rank-1 denoised fit (for the plug-in mask + a_g) ---
gv <- mc@i+1L; cv <- rep.int(seq_len(C), diff(mc@p)); xv <- mc@x
rowN <- as.numeric(rowSums(mc)); colN <- as.numeric(colSums(mc))
fs <- function(val,idx,n){ o<-numeric(n); s<-rowsum(val,idx); o[as.integer(rownames(s))]<-s[,1]; o }
d <- colN; d[d==0]<-1e-6; a <- rowN/sum(rowN); mask <- logical(length(xv))
for(o in 1:6){ for(it in 1:15){
  if(any(mask)){mrN<-fs(xv[mask],gv[mask],G);mrD<-fs(d[cv[mask]],gv[mask],G)}else{mrN<-numeric(G);mrD<-numeric(G)}
  a<-(rowN-mrN)/(sum(d)-mrD); a[!is.finite(a)|a<0]<-0; a<-a/sum(a)
  if(any(mask)){mcN<-fs(xv[mask],cv[mask],C);mcA<-fs(a[gv[mask]],cv[mask],C)}else{mcN<-numeric(C);mcA<-numeric(C)}
  d<-(colN-mcN)/(1-mcA); d[!is.finite(d)|d<0]<-0 }
  mask <- ppois(xv-1, a[gv]*d[cv], lower.tail=FALSE) < 1e-6 }
mu <- a[RAMAC]*d; plug <- (v>0) & (ppois(v-1, mu, lower.tail=FALSE) < 1e-6)

## --- real depth_fix assignment ---
cat("running depth_fix on full matrix...\n")
A <- contingency_assign(mc, q=0.05, refit=10, min_count=2, cell_margin="ambient", tail="hyper", fdr="GS")$assigned
asg <- as.logical(A[RAMAC,]); n_other_asg <- Matrix::colSums(A) - asg
mc0 <- mc; mc0[RAMAC,] <- 0; max_other <- colMaxs(as(mc0,"CsparseMatrix"))

## --- figures ---
be <- function(maxc, thr=12, nlog=32){ b<-c(0:thr, round(exp(seq(log(thr+1),log(maxc+1),length.out=nlog)))); b<-sort(unique(b[b>=0&b<=maxc])); c(b,maxc+1) }
histpanel <- function(sel, file, main, legtxt){
  mxk<-max(v); E<-be(mxk); mb<-length(E)-1; lo<-E[1:mb]; hi<-E[2:(mb+1)]-1; wid<-hi-lo+1
  tot<-sapply(1:mb,function(i) sum(v>=lo[i]&v<=hi[i])); sub<-sapply(1:mb,function(i) sum(sel&v>=lo[i]&v<=hi[i]))
  dt<-tot/C/wid; ds<-sub/C/wid; dl<-lo+0.5; dr<-hi+1.5; FLOOR<-3e-8
  png(file.path(OUT,file), width=1050, height=720, res=140)
  par(mar=c(4.4,4.8,3.2,1), mgp=c(2.7,0.8,0), cex.axis=1.05, cex.lab=1.15)
  plot(NA, log="xy", xlim=c(0.8,mxk+1), ylim=c(FLOOR,2), xlab="gRNA UMI count + 1", ylab="density  P(count)/count", main=main)
  rect(dl,FLOOR,dr,pmax(dt,FLOOR),col="grey85",border="grey60",lwd=0.6)
  rect(dl,FLOOR,dr,pmax(ds,FLOOR),col=rgb(0.75,0.1,0.1,0.8),border=NA)
  rect(8,FLOOR,70,2,col=rgb(1,0.6,0,0.08),border=NA); text(24,0.25,"empty gap",col="darkorange3",cex=0.9)
  legend("topright", bty="n", cex=0.95, fill=c("grey85","firebrick"), border=c("grey60",NA), legend=legtxt)
  dev.off()
}
histpanel(plug, "ramac_labeled_cells.png", "RAMAC: cells labeled 'contains RAMAC' (plug-in Poisson mask, p<1e-6)",
          c("ambient (not labeled)","labeled: contains RAMAC"))
histpanel(asg,  "ramac_depthfix_assignments.png", "RAMAC: cells assigned by depth_fix\n(hypergeometric + Guo-Sarkar FDR, min_count=2)",
          c("not assigned (soup)","assigned: contains RAMAC (depth_fix)"))

## --- cross-reference ---
sink(file.path(OUT,"ramac_crossref.txt"), split=TRUE)
cat("=== RAMAC assignments: depth_fix vs plug-in, and doublet cross-reference ===\n\n")
cat(sprintf("depth_fix total assignments (all guides): %d over %d cells = %.2f guides/cell\n",
            sum(A), C, sum(A)/C))
cat(sprintf("RAMAC assigned by depth_fix: %d (count %d..%d)\n", sum(asg), min(v[asg]), max(v[asg])))
cat(sprintf("RAMAC labeled by plug-in mask: %d (count %d..%d)\n", sum(plug), min(v[plug]), max(v[plug])))
cat(sprintf("agreement on count>=3: both=%d ; depth_fix-only (all count-2)=%d\n\n", sum(asg&plug), sum(asg&!plug)))
cat("Cross-reference: does the cell ALSO have another guide with raw count>=30?\n")
grp <- function(sel,lab){ if(sum(sel)==0){return()}; mo<-max_other[sel]
  cat(sprintf("  %-42s n=%3d  another-guide>=30: %3d (%2.0f%%)  median-other=%.0f\n",
              lab, sum(sel), sum(mo>=30), 100*mean(mo>=30), median(mo))) }
grp(v==2 & asg,  "count-2, ASSIGNED by depth_fix")
grp(v==2 & !asg, "count-2, NOT assigned (soup)")
grp(v>=3 & v<=7, "count 3-7 (gap shoulder, all assigned)")
grp(v>=70,       "count >=70 (signal mode)")
low <- asg & v<=7
cat(sprintf("\nRAMAC depth_fix calls with count<=7: %d ; of these, %d are in cells with another guide>=30 (probable spillover),\n  %d are clean (no other integration -> plausibly real weak RAMAC).\n",
            sum(low), sum(max_other[low]>=30), sum(max_other[low]<30)))
sink()
cat("\nwrote ramac_labeled_cells.png, ramac_depthfix_assignments.png, ramac_crossref.txt\n")
