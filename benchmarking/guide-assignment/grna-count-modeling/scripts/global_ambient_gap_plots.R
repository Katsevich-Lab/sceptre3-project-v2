## ============================================================================
## Aggregate + plot the GLOBAL clean-gap ambient-Poisson analysis.
## Reads results/global_ambient_poisson/perguide_*.csv (from global_ambient_gap.R),
## joins the aggregate soup var/mean (ambient_validation.csv), and produces the
## cross-dataset figures + summary table backing global_ambient_poisson.qmd.
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix)})
OUT <- "results/global_ambient_poisson"

files <- list.files(OUT, pattern="^perguide_.*\\.csv$", full.names=TRUE)
pg <- do.call(rbind, lapply(files, read.csv, stringsAsFactors=FALSE))
pg <- pg[is.finite(pg$c2_oe_den) & is.finite(pg$max_amb_rate), ]

## ---- per-dataset summary ----------------------------------------------------
agg <- do.call(rbind, lapply(split(pg, pg$dataset), function(d) data.frame(
  dataset=d$dataset[1], n_gap=nrow(d),
  med_c2_single=round(median(d$c2_oe_single,na.rm=TRUE),2),
  med_c2_den   =round(median(d$c2_oe_den,na.rm=TRUE),2),
  med_tail_den =round(median(d$tail_oe_den,na.rm=TRUE),2),
  med_max_rate =round(median(d$max_amb_rate),3),
  med_a_g      =signif(median(d$a_g),3))))
## join aggregate soup var/mean (the decisive low-rate-weighted test), if present
av_path <- "results/method_decision/ambient_validation.csv"
if(file.exists(av_path)){
  av <- read.csv(av_path, stringsAsFactors=FALSE)
  key <- c(replogle_rd7="replogle_rd7", dctap_k562_lowmoi="dctap_k562_lowmoi_GSE303901",
           dctap_k562_highmoi="dctap_k562_highmoi_GSE303901", endoc_t2d="endoc_t2d_perturbseq_GSE273677",
           a549_crispri="a549_crispri_perturbseq_GSE236304", invivo_cortex="invivo_cortex_perturbseq_GSE249416",
           erythroid_multiome="perturb_multiome_erythroid_GSE274113", gastric_organoid="gastric_organoid_cropseq_GSE280506",
           cd8_tcell_perturbcite="cd8_tcell_perturbcite_GSE279498", ipsc_crispri="ipsc_crispri_hipsci_figshare27989294",
           tcell_cd4="tcell_cd4_perturbseq_GSE314342")
  agg$soup_vmr_mod <- av$disp_mod[match(key[agg$dataset], av$dataset)]
} else agg$soup_vmr_mod <- NA
agg <- agg[order(agg$med_c2_den), ]
write.csv(agg, file.path(OUT,"summary_datasets.csv"), row.names=FALSE)
cat("=== per-dataset summary (sorted by denoised count-2 excess) ===\n"); print(agg, row.names=FALSE)

## dataset display order + colors (by median denoised excess)
ds <- agg$dataset; nD <- length(ds)
pal <- setNames(hcl(h=seq(15,375,length.out=nD+1)[1:nD], c=90, l=55), ds)
RAMAC_row <- pg[grepl("RAMAC", pg$guide), ]

## ---- FIG 1: per-guide denoised count-2 excess by dataset (box + points) ------
png(file.path(OUT,"global_c2_excess_by_dataset.png"), width=1150, height=760, res=140)
par(mar=c(9.5,4.8,3.4,1), mgp=c(3.2,0.8,0), cex.axis=0.95, cex.lab=1.12)
bx <- lapply(ds, function(x) pmax(pg$c2_oe_den[pg$dataset==x], 0.3))
plot(NA, xlim=c(0.5,nD+0.5), ylim=c(0.4, max(unlist(bx))*1.3), log="y", xaxt="n",
     xlab="", ylab="count-2 excess  observed / expected  (denoised depth)",
     main="Clean-gap ambient count-2 excess across guides & datasets\n(RAMAC generalized; denoised-depth-mixed Poisson null)")
abline(h=1, col="grey45", lty=2)
for(i in seq_len(nD)){
  y <- bx[[i]]; x0 <- i
  set.seed(i); jj <- x0 + runif(length(y),-0.22,0.22)
  points(jj, y, pch=16, col=adjustcolor(pal[ds[i]],0.35), cex=0.5)
  segments(x0-0.30, median(y), x0+0.30, median(y), col="black", lwd=3)
  q <- quantile(y,c(.25,.75)); segments(x0,q[1],x0,q[2],col="black",lwd=1.2)
}
## mark RAMAC
if(nrow(RAMAC_row)){ ix <- match("replogle_rd7", ds)
  points(ix, RAMAC_row$c2_oe_den[1], pch=23, bg="gold", col="black", cex=1.7, lwd=1.5)
  text(ix, RAMAC_row$c2_oe_den[1], "RAMAC", pos=4, offset=0.6, cex=0.85, font=2) }
axis(1, at=seq_len(nD), labels=sprintf("%s\n(n=%d)", ds, agg$n_gap), las=2, cex.axis=0.72)
legend("topleft", bty="n", cex=0.8, pch=c(NA,23), lty=c(2,NA), col=c("grey45","black"),
       pt.bg=c(NA,"gold"), legend=c("Poisson (=1): soup explains the 2s","RAMAC (the original case study)"))
dev.off()

## ---- FIG 2: the unifying scatter — excess collapses to 1 as ambient rate rises
png(file.path(OUT,"global_c2_excess_vs_rate.png"), width=1080, height=760, res=140)
par(mar=c(4.8,4.9,3.4,1), mgp=c(2.9,0.8,0), cex.axis=1.0, cex.lab=1.14)
plot(pg$max_amb_rate, pmax(pg$c2_oe_den,0.3), log="xy", pch=16, cex=0.55,
     col=adjustcolor(pal[pg$dataset],0.5),
     xlab="max fitted ambient rate over the guide's cells   a_g * max d_c",
     ylab="count-2 excess  (denoised depth)",
     main="The excess is weak-integration leakage, not soup\n(where the soup can actually make a 2, the Poisson fits)")
abline(h=1, col="grey45", lty=2)
## loess-ish running median in log-rate bins
lr <- log10(pg$max_amb_rate); br <- seq(min(lr),max(lr),length.out=16)
bi <- findInterval(lr, br); mm <- tapply(pg$c2_oe_den, bi, median)
xm <- tapply(pg$max_amb_rate, bi, median)
lines(xm, mm, col="black", lwd=2.5)
if(nrow(RAMAC_row)) points(RAMAC_row$max_amb_rate[1], RAMAC_row$c2_oe_den[1],
                           pch=23, bg="gold", col="black", cex=1.7, lwd=1.5)
legend("topright", bty="n", cex=0.72, pch=16, col=adjustcolor(pal[ds],0.7),
       legend=sprintf("%s", ds), ncol=1)
legend("bottomleft", bty="n", cex=0.82, lwd=c(2.5,1), lty=c(1,2), col=c("black","grey45"),
       legend=c("running median","Poisson (=1)"))
dev.off()

## ---- FIG 3: denoising universally shrinks the excess (single -> denoised) ----
png(file.path(OUT,"global_single_vs_denoised.png"), width=1000, height=740, res=140)
par(mar=c(9.5,4.8,3.4,1), mgp=c(3.2,0.8,0), cex.axis=0.95, cex.lab=1.12)
o <- order(agg$med_c2_single)
plot(NA, xlim=c(0.5,nD+0.5), ylim=c(0.7, max(agg$med_c2_single)*1.25), log="y", xaxt="n",
     xlab="", ylab="median count-2 excess",
     main="Denoising the depth universally shrinks the excess\n(crude/library depth always overstates the soup)")
abline(h=1, col="grey45", lty=2)
for(k in seq_len(nD)){ i<-o[k]
  arrows(k, agg$med_c2_single[i], k, agg$med_c2_den[i], length=0.09, lwd=2.2, col=pal[agg$dataset[i]])
  points(k, agg$med_c2_single[i], pch=1, cex=1.3, col="grey30")
  points(k, agg$med_c2_den[i], pch=16, cex=1.3, col=pal[agg$dataset[i]])
}
axis(1, at=seq_len(nD), labels=agg$dataset[o], las=2, cex.axis=0.72)
legend("topleft", bty="n", cex=0.9, pch=c(1,16), col=c("grey30","black"),
       legend=c("single Poisson (library-implied)","denoised-depth-mixed Poisson"))
dev.off()

cat("\nwrote 3 figures + summary_datasets.csv to", OUT, "\n")
