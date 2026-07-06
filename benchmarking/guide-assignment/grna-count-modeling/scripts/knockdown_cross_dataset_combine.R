## ============================================================================
## Combine the cross-dataset knockdown-by-integration-class results (§7 extension).
## Computes replogle on the SAME generic thresholds (strong>=30 / weak 2..7) as the
## survey datasets for apples-to-apples, joins the survey CSVs, and makes the
## combined table + figure. Also reports, per dataset, the implied functional
## fraction of the weak band:  f = (1 - kd_weak) / (1 - kd_strong).
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix); library(ondisc); library(sceptre)})
OUT <- "results/global_ambient_poisson"
STRONG<-30; WEAKLO<-2; WEAKHI<-7; MIN_STRONG<-5; MIN_WEAK<-5

## ---- replogle, generic thresholds, via the odm pipeline ---------------------
Dr <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")
grna <- initialize_odm_from_backing_file(file.path(Dr,"grna.odm"))
resp <- initialize_odm_from_backing_file(file.path(Dr,"response.odm"))
so   <- readRDS(file.path(Dr,"sceptre_object.rds"))
lib  <- exp(so@covariate_matrix[,"log(response_n_umis)"])
tdf  <- so@grna_target_data_frame; gene_ids <- rownames(resp)
## NT-cell control mask (canonical baseline)
nt_guides <- tdf$grna_id[tdf$grna_target=="non-targeting"]
ntmask <- rep(FALSE, length(lib)); for(gg in nt_guides) ntmask <- ntmask | (as.numeric(grna[gg,])>=30)
pg <- read.csv(file.path(OUT,"perguide_replogle_rd7.csv"), stringsAsFactors=FALSE)  # clean-gap guides only
## use the clean-gap guide list (already targeting candidates) but generic count classes
pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target!="non-targeting" & pg$target %in% gene_ids, ]
rows <- list()
for(i in seq_len(nrow(pg))){
  gc <- as.numeric(grna[pg$guide[i],]); e <- as.numeric(resp[pg$target[i],]); cp <- e/lib*1e4
  base<-mean(cp[ntmask & gc==0]); if(!is.finite(base)||base<=0) next
  strong<-gc>=STRONG; weak<-gc>=WEAKLO&gc<=WEAKHI; one<-gc==1
  if(sum(strong)<MIN_STRONG || sum(weak)<MIN_WEAK) next
  rows[[length(rows)+1]] <- data.frame(dataset="replogle_rd7", guide=pg$guide[i], target=pg$target[i],
    base_cp10k=round(base,3), n_weak=sum(weak), n_strong=sum(strong),
    kd_strong=round(mean(cp[strong])/base,3), kd_weak=round(mean(cp[weak])/base,3),
    kd_one=round(mean(cp[one])/base,3))
}
rep_res <- do.call(rbind, rows); write.csv(rep_res, file.path(OUT,"knockdown_replogle_rd7.csv"), row.names=FALSE)

## ---- gather all datasets -----------------------------------------------------
files <- c(file.path(OUT,"knockdown_replogle_rd7.csv"),
           list.files(OUT, pattern="^knockdown_(endoc|a549|cd8|tcell)\\.csv$", full.names=TRUE))
rd <- function(f){ x <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
                   if(is.null(x) || !nrow(x)) NULL else x }
all <- do.call(rbind, Filter(Negate(is.null), lapply(files, rd)))
pp <- function(d) d[d$base_cp10k>=0.5 & is.finite(d$kd_strong) & d$kd_strong<=0.7 &
                    is.finite(d$kd_weak) & d$n_weak>=5, ]
summ <- do.call(rbind, lapply(split(all, all$dataset), function(d){
  pw <- pp(d)
  data.frame(dataset=d$dataset[1], n_cand=nrow(d), n_power=nrow(pw),
    kd_strong=if(nrow(pw)) round(median(pw$kd_strong),2) else NA,
    kd_weak  =if(nrow(pw)) round(median(pw$kd_weak),2) else NA,
    kd_one   =if(nrow(pw)) round(median(pw$kd_one),2) else NA,
    weak_functional_frac=if(nrow(pw)) round((1-median(pw$kd_weak))/(1-median(pw$kd_strong)),2) else NA)
}))
summ <- summ[order(summ$kd_weak), ]
write.csv(summ, file.path(OUT,"knockdown_cross_dataset_summary.csv"), row.names=FALSE)
cat("=== cross-dataset knockdown by integration class (generic strong>=30 / weak 2-7) ===\n")
print(summ, row.names=FALSE)

## a549 diagnostic: why no power-positive?
a5 <- all[all$dataset=="a549", ]
if(nrow(a5)) cat(sprintf("\n[a549 diag] %d candidates; best (lowest) kd_strong=%.2f, median kd_strong=%.2f, median base=%.2f CP10K\n",
   nrow(a5), min(a5$kd_strong,na.rm=T), median(a5$kd_strong,na.rm=T), median(a5$base_cp10k,na.rm=T)))

## ---- FIGURE: strong vs weak median knockdown per dataset --------------------
ds <- summ$dataset[!is.na(summ$kd_weak)]
png(file.path(OUT,"knockdown_cross_dataset.png"), width=1080, height=680, res=140)
par(mar=c(6.5,4.8,3.4,1), mgp=c(3.0,0.8,0), cex.axis=0.95, cex.lab=1.1)
nD <- length(ds); x <- seq_len(nD)
plot(NA, xlim=c(0.5,nD+0.5), ylim=c(0,1.15), xaxt="n", xlab="",
     ylab="target expression / baseline (median)",
     main="Knockdown by integration class across datasets\n(strong integrations vs the count-2..7 'weak' band; 1.0 = background)")
abline(h=1, col="grey45", lty=2)
sm <- summ[match(ds,summ$dataset),]
segments(x, sm$kd_strong, x, sm$kd_weak, col="grey70", lwd=1.5)
points(x, sm$kd_strong, pch=16, col="firebrick", cex=1.8)
points(x, sm$kd_weak,   pch=16, col="darkorange2", cex=1.8)
points(x, sm$kd_one,    pch=1,  col="grey40", cex=1.3)
text(x, sm$kd_weak, sprintf("%.0f%% real", 100*sm$weak_functional_frac), pos=3, offset=0.6, cex=0.72, col="darkorange3")
text(x, sm$kd_strong, sprintf("%.2f", sm$kd_strong), pos=1, offset=0.5, cex=0.72, col="firebrick")
axis(1, at=x, labels=sprintf("%s\n(n=%d)", ds, sm$n_power), las=1, cex.axis=0.72)
legend("bottomright", bty="n", cex=0.82, pch=c(16,16,1), col=c("firebrick","darkorange2","grey40"),
       legend=c("strong (real integrations)","weak (count 2-7)","count-1"))
dev.off()
cat("\nwrote knockdown_cross_dataset_summary.csv + knockdown_cross_dataset.png\n")
