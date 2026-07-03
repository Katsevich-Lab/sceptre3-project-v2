## ============================================================================
## Replogle ambient-Poisson analysis — the depth-fit, RAMAC case study, and the
## decisive soup-Poisson test. Self-contained (reads the raw Replogle gRNA matrix;
## no scratchpad dependencies). Companion to scripts/replogle_ambient_poisson.R
## (which does the Poisson+NB mixture fit + marginal histograms). Together these
## reproduce every figure in results/replogle_ambient_poisson/ and back the
## writeup replogle_ambient_poisson.qmd.
##
## Question (from Louis): on Replogle, is the ambient (soup) background too heavy
## to be Poisson, even after marginalizing over ambient depth?
## Answer (this script): the soup is *very nearly* Poisson (aggregate var/mean
## ~1.00 at low rates -> ~1.03-1.11 at high rates). The apparent heaviness of the
## raw margin is integration/doublet contamination; the residual few-percent tail
## excess is small and partly weak-integration. The crude depth proxy (library -
## max) FAKES a perfect Poisson fit via a contamination-inflated tail; only the
## fishash-style rank-1 DENOISED depth gives the honest answer.
## Run from the nonparametric-thresholds/ folder.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats)})

DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
OUT  <- "results/replogle_ambient_poisson"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

mc <- as(readRDS(DATA), "CsparseMatrix")            # guides x cells
G <- nrow(mc); C <- ncol(mc)
RAMAC <- which(grepl("RAMAC", rownames(mc)))[1]     # the clean-gap case-study guide
cat(sprintf("Replogle: %d guides x %d cells, %d nonzeros. RAMAC = row %d (%s)\n",
            G, C, length(mc@x), RAMAC, rownames(mc)[RAMAC]))

## ---- (1) depth proxies: crude (library - carrier) --------------------------
L  <- colSums(mc); mx <- colMaxs(mc); d_crude <- L - mx      # soup depth ~ library minus the one carrier

## ---- (2) fishash-style rank-1 DENOISED depth (carrier-masked IPF) -----------
## noise entries N_{g,c} ~ Poisson(a_g d_c); mask carriers (Poisson upper tail
## p<1e-6); refit rank-1 a_g,d_c by IPF over unmasked entries; iterate. Gauge sum_g a_g=1.
gv <- mc@i + 1L; cv <- rep.int(seq_len(C), diff(mc@p)); xv <- mc@x
rowN <- as.numeric(rowSums(mc)); colN <- as.numeric(colSums(mc))
fastsum <- function(val, idx, n){ o<-numeric(n); s<-rowsum(val, idx); o[as.integer(rownames(s))]<-s[,1]; o }
d <- colN; d[d==0] <- 1e-6; a <- rowN/sum(rowN); mask <- logical(length(xv))
for (outer in 1:6){
  for (inner in 1:15){
    if(any(mask)){ mrN<-fastsum(xv[mask],gv[mask],G); mrD<-fastsum(d[cv[mask]],gv[mask],G) } else { mrN<-numeric(G); mrD<-numeric(G) }
    a <- (rowN-mrN)/(sum(d)-mrD); a[!is.finite(a)|a<0]<-0; a<-a/sum(a)
    if(any(mask)){ mcN<-fastsum(xv[mask],cv[mask],C); mcA<-fastsum(a[gv[mask]],cv[mask],C) } else { mcN<-numeric(C); mcA<-numeric(C) }
    d <- (colN-mcN)/(1-mcA); d[!is.finite(d)|d<0]<-0
  }
  mu <- a[gv]*d[cv]; mask <- ppois(xv-1, mu, lower.tail=FALSE) < 1e-6
  cat(sprintf("  IPF outer %d: masked (carrier) entries = %d (%.3f%%)\n", outer, sum(mask), 100*mean(mask)))
}
d_denoised <- d
cat("\n--- depth comparison (per cell) ---\n")
cat(sprintf("crude  d=L-max : mean=%.1f median=%.1f E[d^2]/E[d]=%.1f\n", mean(d_crude), median(d_crude), mean(d_crude^2)/mean(d_crude)))
cat(sprintf("denoised rank-1: mean=%.1f median=%.1f E[d^2]/E[d]=%.1f  cor(crude,denoised)=%.3f\n",
            mean(d_denoised), median(d_denoised), mean(d_denoised^2)/mean(d_denoised), cor(d_crude, d_denoised)))

## ---- (3) RAMAC: the clean bimodal gap --------------------------------------
v <- as.numeric(mc[RAMAC,]); mxk <- max(v)
cat(sprintf("\nRAMAC counts 1..90 that occur:\n")); print(table(v[v>=1 & v<=90]))
be <- function(maxc, thr=12, nlog=32){ b<-c(0:thr, round(exp(seq(log(thr+1),log(maxc+1),length.out=nlog)))); b<-sort(unique(b[b>=0&b<=maxc])); c(b,maxc+1) }
tabk <- tabulate(v+1, nbins=mxk+1); E<-be(mxk); mb<-length(E)-1; lo<-E[1:mb]; hi<-E[2:(mb+1)]-1; wid<-hi-lo+1
dens <- sapply(1:mb,function(i) sum(tabk[(lo[i]+1):(hi[i]+1)]))/C/wid; dl<-lo+0.5; dr<-hi+1.5
png(file.path(OUT,"gap_hist_RAMAC.png"), width=1000, height=720, res=140)
par(mar=c(4.4,4.8,3.2,1), mgp=c(2.7,0.8,0), cex.axis=1.05, cex.lab=1.15); FLOOR<-3e-8
plot(NA, log="xy", xlim=c(0.8,mxk+1), ylim=c(FLOOR,2), xlab="gRNA UMI count + 1", ylab="density   P(count) / count",
     main="Replogle guide RAMAC — marginal UMI counts\n(clean empty gap between ambient and signal modes)")
rect(dl, FLOOR, dr, pmax(dens,FLOOR), col="grey82", border="grey55", lwd=0.7)
rect(8, FLOOR, 70, 2, col=rgb(1,0.6,0,0.10), border=NA); text(24, 0.25, "empty gap\n(no cells\ncount 8-69)", col="darkorange3", cex=0.9)
dev.off()

## ---- (4) RAMAC lower mode: single vs depth-mixed Poisson (crude depth) ------
amb <- v <= 7; va <- v[amb]; ks <- 0:7; obs <- sapply(ks, function(k) sum(va==k))
lam1 <- mean(va); da_c <- d_crude[amb]; ag_c <- sum(va)/sum(da_c); lam2 <- ag_c*da_c
pm1 <- dpois(ks, lam1); pm2 <- sapply(ks, function(k) mean(dpois(k, lam2)))
png(file.path(OUT,"ramac_poisson_lowmode.png"), width=1000, height=720, res=140)
par(mar=c(4.4,4.8,3.2,1), mgp=c(2.7,0.8,0), cex.axis=1.05, cex.lab=1.15); ymn<-1e-7
plot(ks, pmax(obs/length(va),ymn), log="y", type="n", xlim=c(-0.3,7.3), ylim=c(ymn,1),
     xlab="gRNA UMI count", ylab="P(count)", main="RAMAC ambient (lower mode): observed vs Poisson fits")
for(i in seq_along(ks)) segments(ks[i], ymn, ks[i], max(obs[i]/length(va),ymn), col="grey70", lwd=8, lend=1)
points(ks, pmax(pm1,ymn), pch=1, col="royalblue", cex=1.6, lwd=2); lines(ks, pmax(pm1,ymn), col="royalblue", lty=3)
points(ks, pmax(pm2,ymn), pch=16, col="firebrick", cex=1.4); lines(ks, pmax(pm2,ymn), col="firebrick")
legend("topright", bty="n", cex=1.0, legend=c("observed","single Poisson","depth-mixed Poisson (crude depth)"),
       col=c("grey60","royalblue","firebrick"), pch=c(15,1,16), lty=c(NA,3,1), lwd=c(NA,2,2))
dev.off()

## ---- (5) RAMAC linear-x histogram, crude depth-mixed Poisson ---------------
exp_c <- sapply(ks, function(k) mean(dpois(k, lam2)))*length(va)
png(file.path(OUT,"ramac_ambient_linear.png"), width=1000, height=720, res=140)
par(mar=c(4.4,4.8,3.4,1), mgp=c(2.8,0.8,0), cex.axis=1.1, cex.lab=1.2); ymin<-0.5
plot(NA, xlim=c(-0.5,7.5), ylim=c(ymin,2e6), log="y", xaxt="n", xlab="gRNA UMI count", ylab="number of cells",
     main="RAMAC ambient counts vs depth-mixed Poisson fit (crude depth)")
axis(1, at=ks); rect(ks-0.34, ymin, ks+0.34, pmax(obs,ymin), col="grey82", border="grey45", lwd=0.8)
lines(ks, pmax(exp_c,ymin), col="firebrick", lwd=1.8); points(ks, pmax(exp_c,ymin), pch=16, col="firebrick", cex=1.5)
text(ks[obs>0], obs[obs>0], obs[obs>0], pos=3, cex=0.85, col="grey25", offset=0.35)
legend("topright", bty="n", cex=1.0, legend=c("observed cells","depth-mixed Poisson (a_g*d_c)"),
       col=c("grey70","firebrick"), pch=c(15,16), lty=c(NA,1), lwd=c(NA,1.8))
dev.off()

## ---- (6) RAMAC crude vs DENOISED depth (the reversal) ----------------------
fitpred <- function(dvec){ dd<-dvec[amb]; ag<-sum(va)/sum(dd); list(exp=sapply(ks,function(k) mean(dpois(k,ag*dd))*length(va)), maxrate=ag*max(dd)) }
fc <- fitpred(d_crude); fd <- fitpred(d_denoised)
cat(sprintf("\nRAMAC ambient fit: crude max-rate=%.3f, denoised max-rate=%.3f\n", fc$maxrate, fd$maxrate))
cat("  k  observed  crude-depth  denoised-depth\n")
for(i in seq_along(ks)) cat(sprintf("  %d  %8d  %10.1f  %13.1f\n", ks[i], obs[i], fc$exp[i], fd$exp[i]))
png(file.path(OUT,"ramac_denoised_fit.png"), width=1050, height=730, res=140)
par(mar=c(4.4,4.8,3.4,1), mgp=c(2.8,0.8,0), cex.axis=1.1, cex.lab=1.2); ymin<-0.5
plot(NA, xlim=c(-0.5,7.5), ylim=c(ymin,2e6), log="y", xaxt="n", xlab="gRNA UMI count", ylab="number of cells",
     main="RAMAC ambient: depth-mixed Poisson with crude vs denoised depth")
axis(1, at=ks); rect(ks-0.34, ymin, ks+0.34, pmax(obs,ymin), col="grey82", border="grey45", lwd=0.8)
lines(ks, pmax(fc$exp,ymin), col="royalblue", lwd=1.8); points(ks, pmax(fc$exp,ymin), pch=1, col="royalblue", cex=1.5, lwd=2)
lines(ks, pmax(fd$exp,ymin), col="firebrick", lwd=1.8); points(ks, pmax(fd$exp,ymin), pch=16, col="firebrick", cex=1.4)
text(ks[obs>0], obs[obs>0], obs[obs>0], pos=3, cex=0.8, col="grey25", offset=0.3)
legend("topright", bty="n", cex=1.0, legend=c("observed cells","crude depth (L-max)","denoised rank-1 depth"),
       col=c("grey70","royalblue","firebrick"), pch=c(15,1,16), lty=c(NA,1,1), lwd=c(NA,1.8,1.8))
dev.off()

## ---- (7) DECISIVE TEST: non-circular + aggregate soup dispersion -----------
sink(file.path(OUT,"soup_poisson_test.txt"), split=TRUE)
cat("=== Replogle soup-Poisson decisive test (denoised rank-1 fit) ===\n\n")
mu_nz <- a[gv]*d_denoised[cv]; mask <- ppois(xv-1, mu_nz, lower.tail=FALSE) < 1e-6
cat(sprintf("carriers masked: %d nonzero entries (%.3f%%)\n\n", sum(mask), 100*mean(mask)))
## Part 1: RAMAC counts in cells whose SOLE carrier is a DIFFERENT guide (non-circular)
mgc <- cv[mask]; mgg <- gv[mask]; tabc <- tabulate(mgc, nbins=C); singlet <- tabc==1
carrierG <- integer(C); iss <- singlet[mgc]; carrierG[mgc[iss]] <- mgg[iss]
clean <- singlet & carrierG != RAMAC & carrierG > 0
vcln <- v[clean]; dcln <- d_denoised[clean]; Nc <- length(vcln)
aR <- sum(vcln)/sum(dcln); lam <- aR*dcln
cat("Part 1 -- RAMAC counts in", Nc, "cells whose sole carrier is a DIFFERENT guide (non-circular)\n")
cat(sprintf("  a_RAMAC=%.3e  max ambient rate=%.3f\n", aR, aR*max(dcln)))
cat(sprintf("  %3s %10s %14s %8s\n","k","observed","exp Poisson","o/e"))
for (k in 0:7){ o<-sum(vcln==k); e<-mean(dpois(k,lam))*Nc; cat(sprintf("  %3d %10d %14.1f %8.2f\n", k, o, e, o/max(e,1e-9))) }
## Part 2: aggregate soup dispersion var/mean vs rate mu (zeros analytic)
la<-log(a[a>0]); ld<-log(d_denoised[d_denoised>0]); ha<-hist(la,breaks=100,plot=FALSE); hd<-hist(ld,breaks=250,plot=FALSE)
LM<-outer(ha$mids,hd$mids,"+"); W<-outer(ha$counts,hd$counts,"*"); nb<-45
mubrk<-seq(min(LM)-1e-9,max(LM)+1e-9,length.out=nb+1)
acc<-function(idx,val){ s<-tapply(val,factor(idx,levels=1:nb),sum); s[is.na(s)]<-0; as.numeric(s) }
nslot<-acc(findInterval(as.vector(LM),mubrk),as.vector(W))
nmask<-acc(findInterval(log(mu_nz[mask]),mubrk),rep(1,sum(mask)))
keep<-!mask & mu_nz>0; bi<-findInterval(log(mu_nz[keep]),mubrk); Nz<-xv[keep]
sN<-acc(bi,Nz); sN2<-acc(bi,Nz^2); nsl<-pmax(nslot-nmask,0); mean_b<-sN/nsl; disp<-(sN2/nsl-mean_b^2)/mean_b
muctr<-exp((mubrk[-1]+mubrk[-(nb+1)])/2)
cat("\nPart 2 -- aggregate soup dispersion var/mean by fitted rate mu (Poisson => 1)\n")
cat(sprintf("  %10s %12s %10s %10s\n","mu","n_slots","mean","var/mean"))
for (b in which(nsl>1e4 & is.finite(disp))) if (b%%3==0 || muctr[b]>0.05)
  cat(sprintf("  %10.4f %12.3e %10.5f %10.3f\n", muctr[b], nsl[b], mean_b[b], disp[b]))
sink()
png(file.path(OUT,"soup_dispersion_vs_rate.png"), width=1000, height=720, res=140)
par(mar=c(4.6,4.8,3.2,1), mgp=c(2.8,0.8,0), cex.axis=1.05, cex.lab=1.15); ok<-nsl>1e4 & is.finite(disp) & mean_b>0
plot(muctr[ok], disp[ok], log="x", type="b", pch=16, col="firebrick", ylim=c(0, max(2, max(disp[ok],na.rm=TRUE))),
     xlab="fitted ambient rate  mu = a_g * d_c", ylab="soup dispersion  var / mean",
     main="Replogle aggregate soup: dispersion vs rate\n(Poisson => 1 at all rates)")
abline(h=1, col="grey50", lty=2); dev.off()

saveRDS(list(d_crude=summary(d_crude), d_denoised=summary(d_denoised),
             cor=cor(d_crude,d_denoised), disp=data.frame(mu=muctr,nsl=nsl,mean=mean_b,disp=disp)),
        file.path(OUT,"replogle_ambient_analysis.rds"))
cat("\nwrote figures + soup_poisson_test.txt to", OUT, "\n")
