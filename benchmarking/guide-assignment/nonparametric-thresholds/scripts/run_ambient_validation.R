## Run the ambient-assumption validation across ALL datasets (anti-overfitting).
suppressMessages({library(Matrix); library(sparseMatrixStats)})
source(file.path(getwd(),"scripts","ambient_validation.R"))
OUT <- "results/method_decision"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

SURV <- path.expand("~/data/external/perturbseq-survey")
REPL <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
grna_path <- function(dd){
  for(p in c(file.path(dd,"sceptre","grna_matrix_aligned.rds"), file.path(dd,"grna_matrix.rds")))
    if(file.exists(p)) return(p); NA
}
datasets <- list()
for(dd in list.dirs(SURV, recursive=FALSE)){ p<-grna_path(dd); if(!is.na(p)) datasets[[basename(dd)]] <- p }
datasets[["replogle_rd7"]] <- REPL

summ <- list(); disps <- list()
for(nm in names(datasets)){
  cat("=== ambient validation:", nm, "===\n")
  g <- tryCatch(as(readRDS(datasets[[nm]]),"CsparseMatrix"), error=function(e) NULL)
  if(is.null(g)){ cat("  load fail\n"); next }
  if(nrow(g)>ncol(g)) g <- t(g)
  # cap cells for the very large ones (rank-1 fit is cheap but keep memory sane)
  if(ncol(g) > 300000){ set.seed(1); g <- g[, sort(sample(ncol(g), 300000))] }
  r <- tryCatch(summarize_ambient(g, nm), error=function(e){cat("  ERR:",conditionMessage(e),"\n"); NULL})
  if(is.null(r)) next
  summ[[nm]] <- r$summary; disps[[nm]] <- r$disp
  print(r$summary, row.names=FALSE)
  saveRDS(disps, file.path(OUT,"ambient_dispersion_curves.rds"))
  st <- do.call(rbind, summ); write.csv(st, file.path(OUT,"ambient_validation.csv"), row.names=FALSE)
}
st <- do.call(rbind, summ); rownames(st)<-NULL
write.csv(st, file.path(OUT,"ambient_validation.csv"), row.names=FALSE)
cat("\n==== AMBIENT VALIDATION SUMMARY ====\n"); print(st, row.names=FALSE)

## figure: dispersion vs rate, all datasets overlaid
allc <- do.call(rbind, disps)
png(file.path(OUT,"ambient_dispersion_all.png"), width=1100, height=760, res=140)
par(mar=c(4.5,4.6,2.5,1), mgp=c(2.6,0.7,0))
dsn <- unique(allc$dataset); cols <- rainbow(length(dsn))
plot(NA, log="x", xlim=c(1e-4, 1), ylim=c(0.8, 2),
     xlab="fitted ambient rate  mu = a_g d_c", ylab="soup dispersion var/mean",
     main="Ambient soup dispersion across datasets (Poisson => 1)")
abline(h=1, col="grey50", lty=2)
for(k in seq_along(dsn)){
  s <- allc[allc$dataset==dsn[k] & allc$nslot>1e3 & is.finite(allc$disp),]
  s <- s[order(s$mu),]; lines(s$mu, s$disp, col=cols[k], lwd=1.8, type="b", pch=16, cex=0.5)
}
legend("topleft", legend=sub("_"," ",substr(dsn,1,18)), col=cols, lwd=2, cex=0.6, bty="n")
dev.off()
cat("wrote ambient_dispersion_all.png\n")
