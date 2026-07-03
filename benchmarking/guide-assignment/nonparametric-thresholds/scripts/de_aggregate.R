## Aggregate the DE grid: read every calibration.csv + power.csv, compute rich
## calibration and power metrics, and rank methods across datasets (anti-overfitting:
## a method must be good CONSISTENTLY, not on one dataset).
suppressMessages(library(Matrix))
ROOT <- "results/method_decision/de_native"
OUT  <- "results/method_decision"

read_p <- function(f){ if(!file.exists(f)) return(numeric(0)); r<-read.csv(f); p<-r$p_value; p[!is.na(p)] }
read_fc <- function(f){ if(!file.exists(f)) return(numeric(0)); r<-read.csv(f)
  if(!"fold_change" %in% names(r)) return(numeric(0)); fc<-r$fold_change; fc[is.finite(fc) & fc>0] }

rows <- list()
for(ds in list.dirs(ROOT, recursive=FALSE)){
  for(md in list.dirs(ds, recursive=FALSE)){
    cp <- read_p(file.path(md,"calibration.csv"))
    pp <- read_p(file.path(md,"power.csv"))
    fc <- read_fc(file.path(md,"power.csv"))
    sm <- if(file.exists(file.path(md,"summary.csv"))) read.csv(file.path(md,"summary.csv")) else NULL
    if(is.null(sm)) next
    # calibration metrics
    calm <- if(length(cp)>10){
      exp_med <- 0.5
      data.frame(cal_n=length(cp), t1e05=mean(cp<0.05), t1e01=mean(cp<0.01),
        t1e001=mean(cp<0.001), ks=as.numeric(ks.test(cp,"punif")$statistic),
        # QQ central slope on -log10 (uniform => slope 1); >1 = inflated
        mlog10_mean=mean(-log10(cp+1e-300)),   # uniform => ~0.434
        med_dev=abs(median(cp)-0.5))
    } else data.frame(cal_n=length(cp), t1e05=NA,t1e01=NA,t1e001=NA,ks=NA,mlog10_mean=NA,med_dev=NA)
    # power metrics
    powm <- if(length(pp)>0){
      data.frame(pow_n=length(pp), pow_med_p=median(pp), pow_frac05=mean(pp<0.05),
        pow_mlog10=mean(-log10(pp+1e-300)),
        # median knockdown magnitude among on-target pairs (|log2 fold change|)
        med_abs_lfc=if(length(fc)) median(abs(log2(fc))) else NA)
    } else data.frame(pow_n=0,pow_med_p=NA,pow_frac05=NA,pow_mlog10=NA,med_abs_lfc=NA)
    rows[[paste(basename(ds),basename(md))]] <- cbind(
      data.frame(dataset=basename(ds), method=basename(md), moi=sm$moi[1], n_calls=sm$n_calls[1]),
      calm, powm)
  }
}
tab <- do.call(rbind, rows); rownames(tab)<-NULL
## trustworthy for DE calibration = genome-wide gene expression (excludes dctap's
## 93-gene targeted panel and cd8's degenerate ultra-shallow guide capture).
TRUST <- c("a549","endoc","replogle","gasperini")
tab$trustworthy <- tab$dataset %in% TRUST
write.csv(tab, file.path(OUT,"de_aggregate.csv"), row.names=FALSE)
cat("=== per (dataset, method) ===\n")
print(tab[order(tab$dataset,tab$method), c("dataset","method","moi","n_calls","cal_n","t1e05","ks","mlog10_mean","pow_n","pow_med_p","pow_frac05","med_abs_lfc")], row.names=FALSE, digits=3)

## rank methods within each dataset (1=best), then average ranks across datasets.
## calibration: closeness of t1e05 to 0.05 (|t1e05-0.05|) and KS (lower=better).
## power: pow_med_p (lower=better), med_abs_lfc (higher=better).
meths <- unique(tab$method); dss <- unique(tab$dataset)
rk <- function(vals, lower_better=TRUE){ r<-rank(if(lower_better) vals else -vals, na.last="keep"); r }
agg <- data.frame(method=meths)
cal_rank <- pow_rank <- matrix(NA, length(meths), length(dss), dimnames=list(meths,dss))
for(d in dss){
  sub <- tab[tab$dataset==d,]
  # calibration score: |t1e05 - 0.05| + ks (both lower better), combined rank
  cscore <- abs(sub$t1e05-0.05) + sub$ks
  pscore <- sub$pow_med_p                      # lower = more power
  cal_rank[sub$method, d] <- rk(cscore, TRUE)
  pow_rank[sub$method, d] <- rk(pscore, TRUE)
}
agg$mean_cal_rank <- round(rowMeans(cal_rank[meths,,drop=FALSE], na.rm=TRUE),2)
agg$mean_pow_rank <- round(rowMeans(pow_rank[meths,,drop=FALSE], na.rm=TRUE),2)
agg$mean_t1e05 <- round(tapply(tab$t1e05, tab$method, mean, na.rm=TRUE)[meths],3)
agg$mean_ks    <- round(tapply(tab$ks, tab$method, mean, na.rm=TRUE)[meths],3)
agg$mean_pow_med_p <- round(tapply(tab$pow_med_p, tab$method, mean, na.rm=TRUE)[meths],3)
agg$mean_lfc   <- round(tapply(tab$med_abs_lfc, tab$method, mean, na.rm=TRUE)[meths],3)
agg <- agg[order(agg$mean_cal_rank+agg$mean_pow_rank),]
write.csv(agg, file.path(OUT,"de_method_ranking.csv"), row.names=FALSE)
cat("\n=== METHOD RANKING across ALL datasets (lower rank = better) ===\n")
print(agg, row.names=FALSE)

## TRUSTWORTHY-ONLY summary (genome-wide GE datasets)
tt <- tab[tab$trustworthy,]
cat("\n=== TRUSTWORTHY datasets only (genome-wide GE):", paste(unique(tt$dataset),collapse=", "), "===\n")
ts <- aggregate(cbind(t1e05,t1e01,ks,pow_med_p,pow_frac05,med_abs_lfc)~method, data=tt, FUN=function(x) round(mean(x,na.rm=TRUE),3))
print(ts[order(abs(ts$t1e05-0.05)),], row.names=FALSE)
write.csv(ts, file.path(OUT,"de_trustworthy_summary.csv"), row.names=FALSE)

## per-dataset t1e05 table (trustworthy)
cat("\n=== per-dataset t1e05 (trustworthy) ===\n")
w <- reshape(tt[,c("dataset","method","t1e05")], idvar="method", timevar="dataset", direction="wide")
print(w, row.names=FALSE)

## calibration figure: t1e05 by method, per trustworthy dataset
png(file.path(OUT,"de_calibration_compare.png"), width=1100, height=700, res=140)
par(mar=c(7,4.5,3,1), mgp=c(2.6,0.7,0))
meths2 <- intersect(c("thresh3","ambient","fishash","depth_fix","crispat","cleanser","sceptre"), unique(tt$method))
dss2 <- unique(tt$dataset); cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")[seq_along(dss2)]
plot(NA, xlim=c(0.5,length(meths2)+0.5), ylim=c(0, max(0.15, max(tt$t1e05,na.rm=TRUE))),
     xaxt="n", xlab="", ylab="negative-control type-I error @0.05",
     main="Downstream DE calibration by method (genome-wide-GE datasets)")
axis(1, at=seq_along(meths2), labels=meths2, las=2)
abline(h=0.05, col="grey40", lty=2, lwd=2)
for(k in seq_along(dss2)){
  s <- tt[tt$dataset==dss2[k],]; y <- s$t1e05[match(meths2, s$method)]
  points(seq_along(meths2), y, col=cols[k], pch=19, cex=1.6, type="b", lwd=1.5)
}
legend("topright", legend=dss2, col=cols, pch=19, lwd=1.5, bty="n", cex=0.9)
legend("topleft", legend="nominal 0.05", lty=2, col="grey40", bty="n", cex=0.9)
dev.off()
cat("wrote de_calibration_compare.png\n")
