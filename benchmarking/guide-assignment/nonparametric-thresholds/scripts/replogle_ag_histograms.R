## Rank Replogle guides by the fitted rank-1 ambient share a_g, and show the
## per-guide count histograms (across cells) for the extremes, with the
## depth-mixed ambient Poisson (mean_c Pois(k; a_g d_c)) overlaid.
suppressMessages({library(Matrix); library(sparseMatrixStats)})
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
OUT  <- "results/replogle_ambient_poisson"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
mc <- as(readRDS(DATA), "CsparseMatrix"); G <- nrow(mc); C <- ncol(mc)

## --- rank-1 denoised fit (same as replogle_ambient_analysis.R) ---
adf <- file.path(OUT, "replogle_rank1_ad.rds")
if (file.exists(adf)) { ad <- readRDS(adf); a <- ad$a; d <- ad$d } else {
  gv <- mc@i+1L; cv <- rep.int(seq_len(C), diff(mc@p)); xv <- mc@x
  rowN <- as.numeric(rowSums(mc)); colN <- as.numeric(colSums(mc))
  fastsum <- function(val,idx,n){ o<-numeric(n); s<-rowsum(val,idx); o[as.integer(rownames(s))]<-s[,1]; o }
  d <- colN; d[d==0]<-1e-6; a <- rowN/sum(rowN); mask <- logical(length(xv))
  for(outer in 1:6){
    for(inner in 1:15){
      if(any(mask)){mrN<-fastsum(xv[mask],gv[mask],G);mrD<-fastsum(d[cv[mask]],gv[mask],G)}else{mrN<-numeric(G);mrD<-numeric(G)}
      a<-(rowN-mrN)/(sum(d)-mrD); a[!is.finite(a)|a<0]<-0; a<-a/sum(a)
      if(any(mask)){mcN<-fastsum(xv[mask],cv[mask],C);mcA<-fastsum(a[gv[mask]],cv[mask],C)}else{mcN<-numeric(C);mcA<-numeric(C)}
      d<-(colN-mcN)/(1-mcA); d[!is.finite(d)|d<0]<-0
    }
    mask <- ppois(xv-1, a[gv]*d[cv], lower.tail=FALSE) < 1e-6
  }
  saveRDS(list(a=a,d=d), adf)
}

## --- rank guides by a_g ---
nz <- diff(mc@p); nz <- tabulate(mc@i+1L, nbins=G)          # #cells with count>0 per guide
info <- data.frame(idx=1:G, guide=rownames(mc), a_g=a, n_nonzero=nz,
                   max_ct=sparseMatrixStats::rowMaxs(mc), total=as.numeric(rowSums(mc)))
o <- order(info$a_g)
cat("=== a_g summary ===\n"); print(summary(a))
cat("\n=== LOWEST a_g (5) ===\n"); print(head(info[o,c("guide","a_g","n_nonzero","max_ct","total")],5), row.names=FALSE, digits=3)
cat("\n=== HIGHEST a_g (5) ===\n"); print(head(info[rev(o),c("guide","a_g","n_nonzero","max_ct","total")],5), row.names=FALSE, digits=3)

## pick 4 guides spanning a_g: highest, ~66th pct, ~33rd pct, lowest (lowest among
## guides with >=200 nonzero cells so the histogram is not a dead spike).
hi  <- rev(o)[1]
q66 <- o[round(0.66*G)]; q33 <- o[round(0.33*G)]
elig <- info$idx[info$n_nonzero >= 200]
lo  <- elig[which.min(info$a_g[elig])]
picks <- c(hi, q66, q33, lo)
tags  <- c("HIGHEST a_g", "~66th pct", "~33rd pct", "LOWEST a_g (>=200 nz cells)")

short <- function(nm) sub("_ENSG.*","", sub("^[0-9]+_","", nm))
be <- function(maxc, thr=12, nlog=30){ b<-c(0:thr, round(exp(seq(log(thr+1),log(maxc+1),length.out=nlog)))); b<-sort(unique(b[b>=0&b<=maxc])); c(b,maxc+1) }

png(file.path(OUT,"replogle_ag_histograms.png"), width=1500, height=1000, res=140)
par(mfrow=c(2,2), mar=c(4.2,4.4,3.0,1), mgp=c(2.4,0.7,0), cex.axis=0.95)
FLOOR <- 3e-8
for(i in seq_along(picks)){
  g <- picks[i]; v <- as.numeric(mc[g,]); mxk <- max(v)
  tabk <- tabulate(v+1, nbins=mxk+1); E<-be(mxk); mb<-length(E)-1; lo_<-E[1:mb]; hi_<-E[2:(mb+1)]-1; wid<-hi_-lo_+1
  dens <- sapply(1:mb,function(j) sum(tabk[(lo_[j]+1):(hi_[j]+1)]))/C/wid; dl<-lo_+0.5; dr<-hi_+1.5
  ks <- unique(round(c(0:20, exp(seq(log(21),log(mxk+1),length=120))-1))); ks<-ks[ks>=0&ks<=mxk]
  amb <- sapply(ks, function(k) mean(dpois(k, a[g]*d)))       # depth-mixed ambient Poisson (all cells)
  plot(NA, log="xy", xlim=c(0.8,mxk+1), ylim=c(FLOOR,2), xlab="gRNA UMI count + 1", ylab="density  P(count)/count",
       main=sprintf("%s: %s\na_g=%.2e  (mean ambient rate a_g*d = %.3f)", tags[i], short(info$guide[g]), a[g], a[g]*mean(d)))
  rect(dl, FLOOR, dr, pmax(dens,FLOOR), col="grey83", border="grey55", lwd=0.6)
  lines(ks+1, amb, col="royalblue", lwd=2.2)
  if(i==1) legend("topright", bty="n", cex=0.9, legend=c("observed counts","depth-mixed ambient Poisson (a_g d_c)"),
                  fill=c("grey83",NA), border=c("grey55",NA), col=c(NA,"royalblue"), lty=c(NA,1), lwd=c(NA,2.2))
}
dev.off(); cat("\nwrote replogle_ag_histograms.png\n")
