## ============================================================================
## ORTHOGONAL readout: does the guide's TARGET GENE get knocked down in
## strong- vs weak-integration vs background cells? This is the phenotype channel
## that the count-based analyses (rate argument, cross-reference) cannot see.
## For each clean-gap targeting guide in Replogle (CRISPRi), split cells by the
## guide's UMI count and measure the target gene's normalized expression:
##   strong (>=gap_hi) : must show knockdown if these are real integrations
##   weak   (2..gap_lo): the count-2 "excess" cells under test
##   count-1, and count-0 (per-gene baseline)
## If strong knocks down but weak does NOT, the weak cells are background/soup,
## not functional integrations. Aggregated across guides for power (each guide's
## target is lowly expressed, so one guide is underpowered -- e.g. RAMAC).
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix); library(ondisc); library(sceptre)})
source("scripts/datasets.R")   # replogle_rd7_pipeline_dir()
D <- replogle_rd7_pipeline_dir()
OUT <- "results/global_ambient_poisson"

grna <- initialize_odm_from_backing_file(file.path(D,"grna.odm"))
resp <- initialize_odm_from_backing_file(file.path(D,"response.odm"))
so   <- readRDS(file.path(D,"sceptre_object.rds"))
lib  <- exp(so@covariate_matrix[,"log(response_n_umis)"])
tdf  <- so@grna_target_data_frame
gene_ids <- rownames(resp)

## NT-cell mask: cells that are a strong (>=30) carrier of a non-targeting guide.
## This is the canonical Perturb-seq control (unperturbed baseline). We use it as
## the primary per-gene baseline; the complement (guide-absent) agrees to ~2%.
nt_guides <- tdf$grna_id[tdf$grna_target=="non-targeting"]
ntmask <- rep(FALSE, length(lib))
for(gg in nt_guides) ntmask <- ntmask | (as.numeric(grna[gg,]) >= 30)
cat("NT-strong control cells:", sum(ntmask), "\n")

## clean-gap guides (gap_lo>=2) with a target gene present in the expression matrix
pg <- read.csv(file.path(OUT,"perguide_replogle_rd7.csv"), stringsAsFactors=FALSE)
pg <- pg[pg$gap_lo>=2, ]
pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target!="non-targeting" & pg$target %in% gene_ids, ]
cat("clean-gap targeting guides with target gene in GEX:", nrow(pg), "\n")

rows <- vector("list", nrow(pg))
for(i in seq_len(nrow(pg))){
  g <- as.numeric(grna[pg$guide[i], ])
  e <- as.numeric(resp[pg$target[i], ])
  cp <- e/lib*1e4
  base <- mean(cp[ntmask & g==0])                      # per-gene baseline: NT control cells
  if(!is.finite(base) || base<=0) next
  grp <- list(strong=g>=pg$gap_hi[i], weak=g>=2 & g<=pg$gap_lo[i], one=g==1)
  mn <- sapply(grp, function(m) if(sum(m)) mean(cp[m]) else NA)
  rows[[i]] <- data.frame(guide=pg$guide[i], target=pg$target[i],
    base_cp10k=round(base,3), n_weak=sum(grp$weak), n_strong=sum(grp$strong),
    kd_strong=round(mn["strong"]/base,3), kd_weak=round(mn["weak"]/base,3),
    kd_one=round(mn["one"]/base,3), row.names=NULL)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT,"replogle_knockdown_by_integration.csv"), row.names=FALSE)

## power-positive guides: target expressed enough AND strong-integration knockdown demonstrably works
pw <- res[res$base_cp10k>=0.5 & is.finite(res$kd_strong) & res$kd_strong<=0.7 &
          is.finite(res$kd_weak) & res$n_weak>=5, ]
cat(sprintf("\npower-positive guides (base>=0.5 CP10K, strong knockdown<=0.7, n_weak>=5): %d\n", nrow(pw)))
summ <- function(x) sprintf("median=%.2f  mean=%.2f  [IQR %.2f-%.2f]", median(x,na.rm=T), mean(x,na.rm=T),
                            quantile(x,.25,na.rm=T), quantile(x,.75,na.rm=T))
cat("  kd_strong (real integrations): ", summ(pw$kd_strong), "\n")
cat("  kd_weak   (count-2..gap cells):", summ(pw$kd_weak),   "\n")
cat("  kd_one    (count-1 cells):     ", summ(pw$kd_one),    "\n")
cat(sprintf("\n  => strong cells knock the target DOWN to %.0f%% of baseline;\n     weak/count-2 cells sit at %.0f%% (100%% = no knockdown = background).\n",
    100*median(pw$kd_strong,na.rm=T), 100*median(pw$kd_weak,na.rm=T)))

## pooled fold-change with power weighting (weight by strong knockdown depth & n_weak)
cat(sprintf("\n  pooled weak knockdown (n_weak-weighted mean): %.2f of baseline\n",
    weighted.mean(pw$kd_weak, pw$n_weak, na.rm=TRUE)))

## ---- FIGURE: knockdown by integration class, per guide + aggregate ----------
png(file.path(OUT,"replogle_knockdown_by_integration.png"), width=1080, height=720, res=140)
par(mar=c(4.6,4.8,3.4,1), mgp=c(2.9,0.8,0), cex.axis=1.0, cex.lab=1.12)
xs <- c(strong=1, weak=2, one=3)
plot(NA, xlim=c(0.5,3.5), ylim=c(0, max(2, quantile(c(pw$kd_strong,pw$kd_weak,pw$kd_one),.95,na.rm=T))),
     xaxt="n", xlab="", ylab="target expression  /  baseline (guide-absent cells)",
     main=sprintf("Replogle CRISPRi: target knockdown by integration class\n(%d power-positive clean-gap guides; 1.0 = no knockdown = background)", nrow(pw)))
abline(h=1, col="grey45", lty=2)
for(k in 1:nrow(pw)){
  y <- c(pw$kd_strong[k], pw$kd_weak[k], pw$kd_one[k])
  lines(xs, y, col=adjustcolor("grey70",0.35), lwd=0.6)
  points(xs, y, pch=16, col=adjustcolor(c("firebrick","darkorange2","grey50"),0.4), cex=0.6)
}
med <- c(median(pw$kd_strong,na.rm=T), median(pw$kd_weak,na.rm=T), median(pw$kd_one,na.rm=T))
segments(xs-0.18, med, xs+0.18, med, lwd=4, col=c("firebrick","darkorange2","grey30"))
text(xs, med, sprintf("%.2f", med), pos=3, offset=0.8, font=2, cex=1.0)
axis(1, at=xs, labels=c("strong\n(real integ., >=gap_hi)","weak\n(count 2..gap_lo)","count-1"), padj=0.5, cex.axis=0.85)
## RAMAC highlighted
ir <- which(grepl("RAMAC", pw$guide))
if(length(ir)) { lines(xs, c(pw$kd_strong[ir],pw$kd_weak[ir],pw$kd_one[ir]), col="blue", lwd=1.6)
  points(xs, c(pw$kd_strong[ir],pw$kd_weak[ir],pw$kd_one[ir]), pch=23, bg="gold", cex=1.3) }
legend("topleft", bty="n", cex=0.8, legend=c("thick bar = median across guides","gold = RAMAC","dashed = no knockdown (background)"))
dev.off()
cat("\nwrote replogle_knockdown_by_integration.{csv,png}\n")
