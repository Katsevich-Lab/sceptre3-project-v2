## Detailed per-dataset figures for the comprehensive writeup.
suppressMessages({library(Matrix)})
ROOT <- "results/method_decision"; DEN <- file.path(ROOT,"de_native")
methc <- c(thresh3="#7f8c8d", ambient="#e67e22", fishash="#7570b3", depth_fix="#c0392b", crispat="#16a085", cleanser="#8e44ad", sceptre="#2c3e50")

read_p <- function(ds,m,f){ p<-file.path(DEN,ds,m,f); if(!file.exists(p)) return(numeric(0)); r<-read.csv(p); v<-r$p_value; v[!is.na(v)] }

## ---- (1) Calibration QQ plots, one panel per genome-wide dataset -------------
gw <- c("replogle","gasperini","a549","endoc")
gw_lab <- c("replogle (low MOI)","gasperini (high MOI)","a549 (high MOI)","endoc (high MOI)")
png(file.path(ROOT,"detail_calibration_qq.png"), width=1200, height=1100, res=140)
par(mfrow=c(2,2), mar=c(4.2,4.2,2.6,1), mgp=c(2.4,0.7,0))
for(i in seq_along(gw)){
  ds <- gw[i]; ms <- names(methc)
  plot(NA, xlim=c(0,3.2), ylim=c(0,3.6), xlab="-log10 expected (uniform)", ylab="-log10 observed",
       main=paste0("Neg-control calibration QQ: ", gw_lab[i]))
  abline(0,1,col="grey60",lty=2)
  for(m in ms){
    p <- read_p(ds,m,"calibration.csv"); if(!length(p)) next
    p <- sort(p); n<-length(p); ex <- -log10(ppoints(n)); ob <- -log10(p)
    # subsample for plotting speed
    idx <- unique(round(seq(1,n,length.out=min(n,600))))
    lines(rev(ex)[idx], rev(ob)[idx], col=methc[m], lwd=1.8)
  }
  legend("topleft", legend=names(methc), col=methc, lwd=2, bty="n", cex=0.8)
}
dev.off(); cat("wrote detail_calibration_qq.png\n")

## ---- (2) Per-dataset ambient dispersion facets ------------------------------
disps <- readRDS(file.path(ROOT,"ambient_dispersion_curves.rds"))
nm <- names(disps)
png(file.path(ROOT,"detail_ambient_dispersion_facets.png"), width=1400, height=1000, res=135)
par(mfrow=c(3,4), mar=c(3.6,3.6,2.2,0.8), mgp=c(2.1,0.6,0), cex.axis=0.85)
for(k in nm){
  d <- disps[[k]]; d <- d[d$nslot>1e3 & is.finite(d$disp),]; d<-d[order(d$mu),]
  plot(d$mu, d$disp, log="x", type="b", pch=16, cex=0.5, col="firebrick", lwd=1.5,
       xlim=c(1e-4,1), ylim=c(0.8, max(2, max(d$disp,na.rm=TRUE))),
       xlab="rate mu=a_g d_c", ylab="var/mean",
       main=sub("_"," ",substr(k,1,20)))
  abline(h=1, col="grey50", lty=2)
}
dev.off(); cat("wrote detail_ambient_dispersion_facets.png\n")

## ---- (3) Per-dataset DE calibration + power (all 7 DE datasets) --------------
agg <- read.csv(file.path(ROOT,"de_aggregate.csv"))
dss <- c("replogle","gasperini","a549","endoc","dctap_lowmoi","dctap_highmoi","cd8_tcell")
png(file.path(ROOT,"detail_de_per_dataset.png"), width=1500, height=800, res=135)
par(mfrow=c(2,4), mar=c(6,4,2.6,0.8), mgp=c(2.4,0.7,0))
ms <- names(methc)
for(ds in dss){
  s <- agg[agg$dataset==ds,]; y <- s$t1e05[match(ms,s$method)]
  if(all(is.na(y))){ plot.new(); title(paste0(ds,"\n(no calibration:\ndegenerate)")); next }
  b<-barplot(y, col=methc[ms], names.arg=ms, las=2, cex.names=0.7, ylim=c(0,max(0.25,max(y,na.rm=TRUE))),
     main=paste0(ds," (",s$moi[1]," MOI)"), ylab="neg-ctrl t1e @0.05")
  abline(h=0.05,lty=2,col="grey30"); text(b,y,sprintf("%.3f",y),pos=3,cex=0.6,xpd=NA)
}
# 8th panel: assignment counts (log) for gasperini as example of over/under-calling
plot.new()
dev.off(); cat("wrote detail_de_per_dataset.png\n")

## ---- (4) Assignment counts per method per dataset (over/under-calling) -------
png(file.path(ROOT,"detail_assignment_counts.png"), width=1200, height=600, res=135)
par(mar=c(7,4.5,2.6,1), mgp=c(2.6,0.7,0))
dss2 <- c("replogle","gasperini","a549","endoc","dctap_lowmoi","dctap_highmoi")
M <- sapply(dss2, function(ds){ s<-agg[agg$dataset==ds,]; s$n_calls[match(names(methc),s$method)] })
rownames(M) <- names(methc)
barplot(M, beside=TRUE, col=methc, las=2, log="y", cex.names=0.8,
        ylab="# assigned (g,c) [log]", main="Assignment counts per method per dataset")
legend("topright", legend=names(methc), fill=methc, bty="n", cex=0.75, ncol=2)
dev.off(); cat("wrote detail_assignment_counts.png\n")
