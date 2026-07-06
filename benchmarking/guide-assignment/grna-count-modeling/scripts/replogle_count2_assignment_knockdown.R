## ============================================================================
## POSITIVE-CONTROL test of depth_fix's count-2 recall.
## §7 showed the count-2 band is a MIXTURE (~15% real weak integrations, ~85%
## background) -- so its aggregate knockdown is diluted. But depth_fix does not
## take all count-2 cells: it ASSIGNS only those whose depth-adjusted evidence is
## strong (2 UMIs in a *shallow* cell -> low ambient rate -> significant). If that
## assignment is meaningful, the count-2 cells depth_fix ASSIGNS should be enriched
## for real integrations and show GREATER on-target knockdown than the count-2
## cells it REJECTS. This uses the on-target knockdown as a positive control for
## the method's detection-floor recall.
##
## Groups per clean-gap targeting guide:
##   count-2 ASSIGNED   (depth_fix calls it) | count-2 REJECTED (depth_fix does not)
##   strong (>=gap_hi, reference)            | NT-control baseline
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats); library(ondisc); library(sceptre)})
source("scripts/contingency_method.R")
source("scripts/datasets.R")   # load_replogle_rd7_de()
OUT <- "results/global_ambient_poisson"

rd <- load_replogle_rd7_de()                            # guides x cells (aligned to response.odm)
mc <- rd$mc; resp <- rd$resp; so <- rd$so
lib <- exp(so@covariate_matrix[,"log(response_n_umis)"]); tdf <- so@grna_target_data_frame
gene_ids <- rownames(resp)

## NT-control mask (canonical baseline)
nt_guides <- tdf$grna_id[tdf$grna_target=="non-targeting"]
ntmask <- rep(FALSE, ncol(mc)); for(gg in nt_guides) ntmask <- ntmask | (as.numeric(mc[gg,])>=30)

## --- run the real depth_fix assignment on the full matrix ---
cat("running depth_fix (contingency_assign, ambient/hyper/GS, min_count=2) on full matrix...\n")
t0 <- proc.time()[3]
A <- contingency_assign(mc, q=0.05, refit=10, min_count=2, cell_margin="ambient", tail="hyper", fdr="GS")$assigned
A <- as(A, "CsparseMatrix")
cat(sprintf("  done in %.0fs; total assignments = %d (%.2f guides/cell)\n", proc.time()[3]-t0, sum(A), sum(A)/ncol(mc)))

## clean-gap targeting guides with measured target
pg <- read.csv(file.path(OUT,"perguide_replogle_rd7.csv"), stringsAsFactors=FALSE)
pg <- pg[pg$gap_lo>=2, ]; pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target!="non-targeting" & pg$target %in% gene_ids, ]
ridx <- match(pg$guide, rownames(mc))

## accumulate per-cell NT-normalised target expression for each group (pooled power)
poolA <- poolU <- poolS <- poolC1 <- numeric(0)          # assigned2 / rejected2 / strong / count1
## ungated pools (all guides with a strong band) for the robustness disclosure
ugA <- ugU <- numeric(0); n_ungated <- 0L
rows <- list()
for(i in seq_len(nrow(pg))){
  r <- ridx[i]; v <- as.numeric(mc[r,]); e <- as.numeric(resp[pg$target[i],]); cp <- e/lib*1e4
  base <- mean(cp[ntmask & v==0]); if(!is.finite(base)||base<=0) next
  asg <- as.logical(A[r,])
  a2 <- which(v==2 & asg); u2 <- which(v==2 & !asg)
  st <- which(v>=pg$gap_hi[i]); c1 <- which(v==1)
  if(length(st)){ ugA <- c(ugA, cp[a2]/base); ugU <- c(ugU, cp[u2]/base); n_ungated <- n_ungated+1L }
  # power gate: target expressed, strong-integration knockdown demonstrably works
  if(base<0.5 || !length(st) || mean(cp[st])/base > 0.7) next
  poolA <- c(poolA, cp[a2]/base); poolU <- c(poolU, cp[u2]/base)
  poolS <- c(poolS, cp[st]/base); poolC1 <- c(poolC1, cp[c1]/base)
  rows[[length(rows)+1]] <- data.frame(guide=pg$guide[i], target=pg$target[i], base_cp10k=round(base,3),
    n_assigned2=length(a2), n_rejected2=length(u2),
    kd_assigned2=if(length(a2)) round(mean(cp[a2])/base,3) else NA,
    kd_rejected2=if(length(u2)) round(mean(cp[u2])/base,3) else NA,
    kd_strong=round(mean(cp[st])/base,3))
}
res <- do.call(rbind, rows); write.csv(res, file.path(OUT,"replogle_count2_assignment_knockdown.csv"), row.names=FALSE)

kd <- function(x) if(length(x)) mean(x) else NA
cat(sprintf("\n=== count-2 knockdown, split by depth_fix assignment (%d power-positive guides) ===\n", nrow(res)))
cat(sprintf("  pooled cells: assigned2=%d  rejected2=%d  strong=%d  count1=%d\n",
    length(poolA), length(poolU), length(poolS), length(poolC1)))
cat(sprintf("  POOLED knockdown (mean target expr / NT baseline):\n"))
cat(sprintf("    count-2 ASSIGNED by depth_fix : %.3f\n", kd(poolA)))
cat(sprintf("    count-2 REJECTED (soup)       : %.3f\n", kd(poolU)))
cat(sprintf("    strong (real integrations)    : %.3f\n", kd(poolS)))
cat(sprintf("    count-1 (background)          : %.3f\n", kd(poolC1)))
fA <- (1-kd(poolA))/(1-kd(poolS)); fU <- (1-kd(poolU))/(1-kd(poolS))
cat(sprintf("  implied functional fraction: assigned2 = %.0f%%   rejected2 = %.0f%%\n", 100*fA, 100*fU))
cat(sprintf("  [robustness] UNGATED (all %d guides w/ strong band): assigned2=%.3f rejected2=%.3f (diff %+.3f)\n",
    n_ungated, kd(ugA), kd(ugU), kd(ugA)-kd(ugU)))

## bootstrap CI on the assigned-vs-rejected difference (pooled cells)
set.seed(1)
db <- replicate(2000, mean(sample(poolA, length(poolA), TRUE)) - mean(sample(poolU, length(poolU), TRUE)))
cat(sprintf("  assigned2 - rejected2 = %.3f  (95%% CI [%.3f, %.3f])\n",
    kd(poolA)-kd(poolU), quantile(db,.025), quantile(db,.975)))

## ---- FIGURE ----
png(file.path(OUT,"replogle_count2_assignment_knockdown.png"), width=1000, height=680, res=140)
par(mar=c(5.2,4.8,3.4,1), mgp=c(2.9,0.8,0), cex.axis=0.95, cex.lab=1.1)
vals <- c(kd(poolS), kd(poolA), kd(poolU), kd(poolC1))
cols <- c("firebrick","seagreen4","darkorange2","grey55")
b <- barplot(vals, col=cols, border=NA, ylim=c(0,1.15), ylab="target expression / NT baseline (pooled)",
  names.arg=c(sprintf("strong\n(n=%d)",length(poolS)), sprintf("count-2\nASSIGNED\n(n=%d)",length(poolA)),
              sprintf("count-2\nREJECTED\n(n=%d)",length(poolU)), sprintf("count-1\n(n=%d)",length(poolC1))),
  main="depth_fix's count-2 assignments are NOT enriched for knockdown\n(Replogle CRISPRi; assigned ~= rejected ~= background, far from strong)")
abline(h=1, col="grey45", lty=2)
text(b, vals, sprintf("%.2f", vals), pos=3, offset=0.4, font=2, cex=0.9)
dev.off()
cat("\nwrote replogle_count2_assignment_knockdown.{csv,png}\n")
