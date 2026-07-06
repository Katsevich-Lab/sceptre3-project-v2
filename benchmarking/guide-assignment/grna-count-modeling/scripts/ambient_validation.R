## ============================================================================
## Cross-dataset validation of the ambient generative-model assumptions.
## For each dataset: fit the fishash-style rank-1 DENOISED ambient (a_g, d_c) by
## carrier-masked IPF, then measure whether the soup is Poisson via the aggregate
## dispersion var/mean vs fitted rate mu=a_g d_c (zeros handled analytically).
## Also reports library/denoised-depth ratio and signal separability.
## This is the anti-overfitting backbone: are Replogle's findings general?
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats)})

## ---- rank-1 Poisson denoised ambient fit (mask carriers, IPF over noise) -----
source("scripts/ambient_lib.R")   # fit_rank1_denoised()

## ---- aggregate soup dispersion var/mean vs rate mu (zeros analytic) ----------
soup_dispersion <- function(fit, nb=40){
  a<-fit$a; d<-fit$d; gv<-fit$gv; cv<-fit$cv; xv<-fit$xv; mask<-fit$mask
  mu_nz <- a[gv]*d[cv]
  la<-log(a[a>0]); ld<-log(d[d>0])
  ha<-hist(la,breaks=80,plot=FALSE); hd<-hist(ld,breaks=200,plot=FALSE)
  LM<-outer(ha$mids,hd$mids,"+"); W<-outer(ha$counts,hd$counts,"*")
  mubrk<-seq(min(LM)-1e-9,max(LM)+1e-9,length.out=nb+1)
  acc<-function(idx,val){ s<-tapply(val,factor(idx,levels=1:nb),sum); s[is.na(s)]<-0; as.numeric(s) }
  nslot<-acc(findInterval(as.vector(LM),mubrk),as.vector(W))
  if(any(mask)){ lmm<-log(mu_nz[mask]); nmask<-acc(findInterval(lmm,mubrk),rep(1,sum(mask))) } else nmask<-numeric(nb)
  keep<-!mask & mu_nz>0; lmz<-log(mu_nz[keep]); Nz<-xv[keep]; bi<-findInterval(lmz,mubrk)
  sN<-acc(bi,Nz); sN2<-acc(bi,Nz^2)
  nsl<-pmax(nslot-nmask,0); mean_b<-sN/nsl; disp<-(sN2/nsl-mean_b^2)/mean_b
  muctr<-exp((mubrk[-1]+mubrk[-(nb+1)])/2)
  data.frame(mu=muctr, nslot=nsl, mean=mean_b, disp=disp)
}

## ---- one-dataset summary -----------------------------------------------------
summarize_ambient <- function(counts, name){
  counts <- as(counts,"CsparseMatrix"); if(nrow(counts)>ncol(counts)) counts<-t(counts)
  fit <- fit_rank1_denoised(counts)
  disp <- soup_dispersion(fit)
  # dispersion at "moderate" rate bins (mu 0.01-0.3) where power exists
  modr <- disp$mu>=0.01 & disp$mu<=0.3 & disp$nslot>1e4 & is.finite(disp$disp)
  hir  <- disp$mu>=0.05 & disp$nslot>1e3 & is.finite(disp$disp)
  # depth ratio: library size vs denoised depth
  L <- fit$colN; dd <- fit$d
  drat <- median(L[L>0]/pmax(dd[L>0],0.5))
  # signal separability: fraction of masked (carrier) entries; gap presence via
  # per-guide: median over guides of (min signal count / max ambient count)
  list(summary=data.frame(dataset=name, n_guides=nrow(counts), n_cells=ncol(counts),
         frac_carrier=mean(fit$mask),
         disp_mod=round(if(any(modr)) weighted.mean(disp$disp[modr],disp$nslot[modr]) else NA,3),
         disp_hi =round(if(any(hir)) max(disp$disp[hir],na.rm=TRUE) else NA,3),
         depth_ratio_lib_over_denoised=round(drat,1),
         med_denoised_depth=round(median(dd[dd>0]),1)),
       disp=cbind(dataset=name, disp))
}
