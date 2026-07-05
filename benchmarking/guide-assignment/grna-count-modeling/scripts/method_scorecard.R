## Summary scorecard across the decisive axes. crispat + CLEANSER included.
## DE axes computed on the COMMON dataset subset (a549, endoc, replogle) so all
## methods are compared on identical data (crispat/cleanser were not run on
## gasperini for cost). Upstream sim from method_overall.csv (+ CLEANSER's prior
## re-scored value, flagged).
OUT <- "results/method_decision"
sim <- read.csv("results/sim_framework/method_overall.csv")
agg <- read.csv(file.path(OUT,"de_aggregate.csv"))
COMMON <- c("a549","endoc","replogle")
de_common <- function(m, col){ v <- agg[agg$method==m & agg$dataset %in% COMMON, col]; mean(v, na.rm=TRUE) }

meths <- c("depth_fix","fishash","crispat","cleanser","thresh3","ambient","sceptre")
lab   <- c("depth_fix\n(ours)","fishash\n(SOTA)","crispat\n(P-Gauss)","CLEANSER\n(MCMC)","thresh3","ambient","sceptre\n(current)")
simj  <- sim$jaccard[match(meths, sim$method)]
simj[meths=="cleanser"] <- 0.92         # prior re-scoring (0.89/0.95); flagged in caption
simfdr<- sim$FDR[match(meths, sim$method)]
det1e <- sapply(meths, de_common, col="t1e05")
depow <- sapply(meths, de_common, col="pow_frac05")

png(file.path(OUT,"method_scorecard.png"), width=1500, height=850, res=140)
par(mfrow=c(2,2), mar=c(5.8,4.5,3,1), mgp=c(2.6,0.7,0))
cols <- c("#c0392b","#7570b3","#16a085","#8e44ad","#7f8c8d","#e67e22","#2c3e50")
bp <- function(v, main, ylab, ref=NA, reflab=""){
  b <- barplot(v, col=cols, names.arg=lab, las=2, cex.names=0.7, main=main, ylab=ylab,
               ylim=c(0, max(v,ref,na.rm=TRUE)*1.2))
  if(!is.na(ref)){ abline(h=ref, lty=2, lwd=2, col="grey30"); text(par("usr")[1], ref, reflab, pos=3, cex=0.7, col="grey30") }
  text(b, v, sprintf("%.3f",v), pos=3, cex=0.72, xpd=NA)
}
bp(simj,  "Upstream accuracy (simulations)", "per-guide Jaccard vs truth")
bp(simfdr,"Upstream FDR control (sims)", "realized FDR", ref=0.05, reflab="nominal 0.05")
bp(det1e, "Downstream calibration (a549+endoc+replogle)", "neg-control type-I error @0.05", ref=0.05, reflab="nominal 0.05")
bp(depow, "Downstream power (a549+endoc+replogle)", "frac on-target p<0.05")
dev.off()
cat("wrote method_scorecard.png\n")
print(data.frame(method=meths, sim_jaccard=round(simj,3), sim_FDR=round(simfdr,3),
                 de_t1e=round(det1e,3), de_power=round(depow,3)), row.names=FALSE)
