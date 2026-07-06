## ============================================================================
## Is the Replogle ambient (soup) distribution too heavy for a Poisson?
## ----------------------------------------------------------------------------
## Louis's concern: for Replogle, the ambient guide-count distributions look
## "too heavy to come from a Poisson, potentially even if you account for the
## extra heaviness you get by marginalizing over the ambient depths."
##
## Test, per guide, using the canonical model (canonical_model.qmd):
##     N_gc = Z_gc + X_gc * W_gc,   Z ~ Poisson(a_g d_c),   W ~ NB(mu_g s_c, theta)
## Fit the Poisson(ambient) + NB(signal) mixture by soft EM and ask whether the
## AMBIENT (X=0) component is genuinely overdispersed relative to Poisson, or
## whether the apparent heaviness is signal (real integrations, incl. low-count
## doublets where the guide is a secondary integration) leaking into the low
## counts.
##
## Decisive test = within-depth-stratum dispersion (var/mean):
##   * raw   (all cells)  -> what Louis sees; huge (~1500-3700x)
##   * soup  (P(signal)<0.5, decontaminated) -> is the SOUP itself Poisson?
##   * depth-mixing prediction 1 + Var(a_g d_c)/E(a_g d_c) for the soup
##
## Depth proxy: d_c = library_size_c - (largest single guide count in cell c).
##   In this LOW-MOI screen the cell's library size is dominated by its one
##   carrier guide, so library size is NOT the soup depth; subtracting the
##   carrier recovers the soup/ambient depth (d_c median 32 vs L_c median 1174).
##
## Data: replogle-rd7 sceptre grna_matrix.rds  (2666 guides x 616,184 cells)
## Output: results/replogle_ambient_poisson/{png,csv}
## ============================================================================
## Run from the grna-count-modeling/ folder (repo convention).
suppressMessages({library(Matrix); library(sparseMatrixStats)})

DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
OUT  <- "results/replogle_ambient_poisson"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

## ---- load + depth proxies --------------------------------------------------
m  <- as(readRDS(DATA), "CsparseMatrix")            # guides x cells (dgCMatrix)
L  <- colSums(m); mx <- colMaxs(m); d <- L - mx     # soup-depth proxy
rs <- rowSums(m)

## ---- pick guides across the abundance spectrum -----------------------------
o     <- order(rs)
picks <- o[round(quantile(seq_along(o), c(.15,.35,.55,.75,.90,.98)))]

## ---- weighted NB dispersion (theta) via MoM on Pearson residuals -----------
theta_mom <- function(x, w, mu){
  W <- sum(w); vv <- sum(w*(x-mu)^2)/W; ex <- vv - sum(w*mu)/W
  if (ex <= 0) return(Inf); (sum(w*mu^2)/W)/ex
}

## ---- soft EM for the per-guide Poisson(depth) + NB(free mean) mixture -------
fit_guide <- function(v, d, iters=25){
  intg <- v >= 15; pi <- max(mean(intg), 1e-6)
  ag <- sum(v[!intg])/sum(d[!intg]); mu <- max(mean(v[intg]), 1); th <- 1
  for (it in 1:iters){
    lam <- ag*d + 1e-12
    l0  <- dpois(v, lam); l1 <- dnbinom(v, mu=mu, size=th)
    r   <- (pi*l1)/((1-pi)*l0 + pi*l1 + 1e-300)
    pi  <- mean(r); ag <- sum((1-r)*v)/sum((1-r)*d)
    mu  <- sum(r*v)/sum(r); th <- theta_mom(v, r, mu)
  }
  list(pi=pi, ag=ag, mu=mu, th=th, r=r, lam=ag*d)
}

## ---- pooled within-depth-decile dispersion over a subset of cells ----------
within_disp <- function(v, d, keep){
  va <- v[keep]; da <- d[keep]
  br <- unique(quantile(da, seq(0,1,0.1))); br[1] <- -Inf; br[length(br)] <- Inf
  b  <- cut(da, br, labels=FALSE)
  M  <- do.call(rbind, tapply(seq_along(va), b,
        function(ix){x<-va[ix]; c(n=length(x), mu=mean(x), v=var(x))}))
  ok <- is.finite(M[,"v"]/M[,"mu"]) & M[,"mu"] > 0
  sum(M[ok,"n"]*(M[ok,"v"]/M[ok,"mu"]))/sum(M[ok,"n"])
}

## ---- run -------------------------------------------------------------------
res <- list(); rows <- list()
for (g in picks){
  nm <- rownames(m)[g]; v <- as.numeric(m[g,])
  f  <- fit_guide(v, d)
  keep <- f$r < 0.5
  lam  <- f$lam[keep]
  rows[[nm]] <- data.frame(
    guide=nm, pi=f$pi, ambient_rate_per_depth=f$ag, signal_mean=f$mu,
    signal_theta=f$th,
    disp_raw_all_cells      = within_disp(v, d, rep(TRUE, length(v))),
    disp_soup_only          = within_disp(v, d, keep),
    disp_depthmix_prediction= 1 + var(lam)/mean(lam))
  res[[nm]] <- list(fit=f, v=v)
}
summ <- do.call(rbind, rows)
write.csv(summ, file.path(OUT, "replogle_ambient_poisson_summary.csv"), row.names=FALSE)
print(summ, digits=3)

## ---- figure ----------------------------------------------------------------
nms   <- names(res)
short <- sub("_ENSG.*","", sub("^[0-9]+_","", nms)); short[grepl("non-target",nms)] <- "non-targeting"
png(file.path(OUT, "replogle_ambient_poisson.png"), width=1500, height=1000, res=140)
layout(matrix(c(1,2,3, 4,5,6, 7,7,7), nrow=3, byrow=TRUE), heights=c(1,1,1.15))
par(mar=c(3.4,3.6,2.4,0.8), mgp=c(2.1,0.7,0), cex.axis=0.9)
for (i in seq_along(nms)){
  f <- res[[nms[i]]]$fit; v <- res[[nms[i]]]$v; mxk <- max(v)
  ov <- sort(unique(v)); obs <- as.numeric(table(v))/length(v)
  ks <- unique(round(c(0:20, exp(seq(log(21), log(mxk+1), length=120))-1)))
  ks <- ks[ks>=0 & ks<=mxk]
  ambP <- sapply(ks, function(k) mean((1-f$pi)*dpois(k, f$lam)))
  sigN <- sapply(ks, function(k) f$pi*dnbinom(k, mu=f$mu, size=f$th))
  plot(ov+1, obs, log="xy", pch=16, cex=0.5, col="grey30",
       xlim=c(1,mxk+1), ylim=c(1e-7,1), xlab="gRNA count + 1", ylab="P(count)",
       main=sprintf("%s  (pi=%.1e, theta=%.2f)", short[i], f$pi, f$th))
  lines(ks+1, ambP+sigN, col="firebrick", lwd=2)
  lines(ks+1, ambP, col="royalblue", lwd=1.6, lty=2)
  lines(ks+1, sigN, col="darkorange", lwd=1.6, lty=3)
  abline(v=16, col="grey80", lty=3)
  if (i==1) legend("topright", bty="n", cex=0.8,
     legend=c("observed","full mixture","ambient Poisson","signal NB"),
     col=c("grey30","firebrick","royalblue","darkorange"),
     pch=c(16,NA,NA,NA), lty=c(NA,1,2,3), lwd=c(NA,2,1.6,1.6))
}
par(mar=c(4.2,4.2,2.2,1))
x <- seq_along(nms)
plot(NA, xlim=c(0.5,length(nms)+0.5), ylim=c(0.6, max(summ$disp_raw_all_cells)*1.5),
     log="y", xaxt="n", xlab="",
     ylab="within-depth-stratum dispersion (var/mean)",
     main="Where the heaviness lives: soup is Poisson; integrations are not")
axis(1, at=x, labels=short, las=2, cex.axis=0.85); abline(h=1, col="grey60", lty=2)
points(x, summ$disp_raw_all_cells,       pch=17, col="firebrick", cex=1.6)
points(x, summ$disp_soup_only,           pch=16, col="royalblue", cex=1.6)
points(x, summ$disp_depthmix_prediction, pch=1,  col="royalblue", cex=1.6)
text(x, summ$disp_raw_all_cells, sprintf("%.0f", summ$disp_raw_all_cells), pos=3, cex=0.7, col="firebrick")
legend("right", bty="n", cex=0.9,
  legend=c("ALL cells (raw ambient margin)","soup only (P(signal)<0.5)",
           "depth-mixing prediction","Poisson (var/mean=1)"),
  col=c("firebrick","royalblue","royalblue","grey60"),
  pch=c(17,16,1,NA), lty=c(NA,NA,NA,2))
dev.off()
cat("\nwrote", file.path(OUT, "replogle_ambient_poisson.png"), "\n")

## ---- marginal UMI-count histogram (log-log) with null + alternative densities
## Done rigorously: log-spaced INTEGER bins; bar height = mass / bin-width (in
## counts) = density w.r.t. count, directly comparable to the model PMF curves
## (also density w.r.t. count). Unit-width bins at low counts, log-spaced above.
sel <- intersect(c("2576_EIF5A_P1P2_ENSG00000132507",
                   "9610_UQCRC1_P1_ENSG00000010256",
                   "4215_IST1_P1P2_ENSG00000182149"), nms)
bin_edges <- function(maxc, thr=12, nlog=30){
  b <- c(0:thr, round(exp(seq(log(thr+1), log(maxc+1), length.out=nlog))))
  b <- sort(unique(b[b>=0 & b<=maxc])); c(b, maxc+1)   # bin i = counts [E[i], E[i+1]-1]
}
png(file.path(OUT, "replogle_marginal_hist.png"), width=1680, height=580, res=150)
par(mfrow=c(1,length(sel)), mar=c(4.0,4.4,3.2,1), mgp=c(2.5,0.7,0), cex.axis=0.95, cex.lab=1.05)
FLOOR <- 3e-8
for (nm in sel){
  f <- res[[nm]]$fit; v <- res[[nm]]$v; nc <- length(v); mxk <- max(v)
  sh <- short[match(nm, nms)]
  tabk <- tabulate(v+1, nbins=mxk+1)
  E <- bin_edges(mxk); mb <- length(E)-1; lo <- E[1:mb]; hi <- E[2:(mb+1)]-1; wid <- hi-lo+1
  mass <- sapply(1:mb, function(i) sum(tabk[(lo[i]+1):(hi[i]+1)]))/nc
  dens <- mass/wid; dl <- lo+0.5; dr <- hi+1.5             # display extent in (count+1)
  ks <- unique(round(c(0:25, exp(seq(log(26), log(mxk+1), length=200))-1))); ks <- ks[ks>=0 & ks<=mxk]
  null <- sapply(ks, function(k) mean((1-f$pi)*dpois(k, f$lam)))       # H0: ambient Poisson
  alt  <- sapply(ks, function(k) f$pi*dnbinom(k, mu=f$mu, size=f$th))  # H1: signal NB
  plot(NA, log="xy", xlim=c(0.8, mxk+1), ylim=c(FLOOR, 2),
       xlab="UMI count + 1", ylab="density   P(count) / count",
       main=sprintf("%s\nP(integrated)=%.1e,  NB theta=%.2f", sh, f$pi, f$th))
  rect(dl, FLOOR, dr, pmax(dens, FLOOR), col="grey85", border="grey55", lwd=0.6)
  lines(ks+1, null, col="royalblue", lwd=2.4); lines(ks+1, alt, col="darkorange", lwd=2.4)
  lines(ks+1, null+alt, col="firebrick", lwd=2.0, lty=2)
  if (nm==sel[1]) legend("topright", bty="n", cex=0.9,
     legend=c("observed (density hist)","null: ambient Poisson","alt: signal NB","full mixture"),
     fill=c("grey85",NA,NA,NA), border=c("grey55",NA,NA,NA),
     col=c(NA,"royalblue","darkorange","firebrick"), lty=c(NA,1,1,2), lwd=c(NA,2.4,2.4,2))
}
dev.off()
cat("wrote", file.path(OUT, "replogle_marginal_hist.png"), "\n")
