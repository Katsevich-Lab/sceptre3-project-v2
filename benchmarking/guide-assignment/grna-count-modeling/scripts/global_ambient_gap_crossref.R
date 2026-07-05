## ============================================================================
## CROSS-REFERENCE the count-2 cells: what ELSE is in them?
## The global writeup claims the below-gap count-2 excess is "weak integration,
## not soup". The rate argument only proves NOT-soup. This script does what the
## RAMAC doc did (doublet cross-reference) for EVERY clean-gap guide, at scale:
## for each count-2 ambient cell, look at the OTHER guides in that cell and
## classify it. Also does the same for count-1 cells and for the guide's own
## signal-mode cells (baselines). Answers directly: are these cells dominated by
## a different strong guide (doublet/spillover), by this guide (own weak
## integration), or by nothing (unresolved background)?
##
## Usage:  Rscript scripts/global_ambient_gap_crossref.R <name> <matrix_rds>
## Run from nonparametric-thresholds/.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats)})

args <- commandArgs(trailingOnly=TRUE)
NAME <- args[1]; PATH <- path.expand(args[2])
OUT  <- "results/global_ambient_poisson"
STRONG <- 30                                   # "clear integration of another guide" (RAMAC doc threshold)

pg_file <- file.path(OUT, sprintf("perguide_%s.csv", NAME))
if(!file.exists(pg_file)){ cat("[",NAME,"] no perguide file; skip\n"); quit(status=0) }
pg <- read.csv(pg_file, stringsAsFactors=FALSE)
pg <- pg[pg$gap_lo >= 2, ]                      # need a count-2 ambient bar
if(!nrow(pg)){ cat("[",NAME,"] no gap_lo>=2 guides; skip\n"); quit(status=0) }

mc <- as(readRDS(PATH), "CsparseMatrix")
if(nrow(mc) > ncol(mc)) mc <- as(t(mc), "CsparseMatrix")
C <- ncol(mc)
colmax  <- colMaxs(mc)                                   # top guide count per cell
nstrong <- as.numeric(Matrix::colSums(mc >= STRONG))     # # guides >= STRONG per cell

gv <- mc@i + 1L; cv <- rep.int(seq_len(C), diff(mc@p)); xv <- mc@x
ridx <- match(pg$guide, rownames(mc)); names(ridx) <- pg$guide
gap_hi <- setNames(pg$gap_hi, pg$guide)

## classify a set of cells given the focal guide's ACTUAL per-cell count (fcnt).
## Excludes the focal guide itself from the 'other strong' tally -- otherwise a
## signal cell (focal >= STRONG) trivially counts itself and p_other_strong -> 1.
## returns: other-carrier>=STRONG / other moderately stronger / focal-is-max
## (colmax <= focal count) ; plus median of the strongest OTHER guide count.
classify <- function(cells, fcnt){
  if(!length(cells)) return(NULL)
  cm <- colmax[cells]; ns <- nstrong[cells]
  other_strong <- (ns - (fcnt >= STRONG)) >= 1   # subtract the focal guide if it is itself strong
  focal_is_max <- cm <= fcnt                      # nobody beats the focal's own count
  other_mod    <- !other_strong & !focal_is_max  # another guide stronger than focal but < STRONG
  other_max <- ifelse(cm > fcnt, cm, NA)          # strongest other guide, when one exists
  c(n=length(cells),
    p_other_strong=mean(other_strong),
    p_other_mod   =mean(other_mod),
    p_focal_max   =mean(focal_is_max),
    med_other_max =median(other_max, na.rm=TRUE))
}

acc2 <- acc1 <- accS <- list()
for(i in seq_len(nrow(pg))){
  r <- ridx[i]; if(is.na(r)) next
  pos <- which(gv==r); vg <- xv[pos]; cg <- cv[pos]
  c2 <- cg[vg==2]; c1 <- cg[vg==1]; sigsel <- vg >= gap_hi[i]; cS <- cg[sigsel]
  acc2[[i]] <- classify(c2, 2); acc1[[i]] <- classify(c1, 1); accS[[i]] <- classify(cS, vg[sigsel])
}
poolw <- function(L){                            # cell-weighted pool across guides
  L <- L[!vapply(L, is.null, logical(1))]
  M <- do.call(rbind, L); w <- M[,"n"]
  c(n_cells=sum(w),
    p_other_strong=round(weighted.mean(M[,"p_other_strong"], w),3),
    p_other_mod   =round(weighted.mean(M[,"p_other_mod"], w),3),
    p_focal_max   =round(weighted.mean(M[,"p_focal_max"], w),3),
    med_other_max =round(weighted.mean(M[,"med_other_max"], w, na.rm=TRUE),0))
}
row2 <- poolw(acc2); row1 <- poolw(acc1); rowS <- poolw(accS)
out <- data.frame(dataset=NAME, cell_set=c("count-2 (tested)","count-1","signal-mode (real integ.)"),
  rbind(row2, row1, rowS), row.names=NULL)
cat(sprintf("\n===== %s : what else is in the cell? (STRONG=%d) =====\n", NAME, STRONG))
print(out, row.names=FALSE)
write.csv(out, file.path(OUT, sprintf("crossref_%s.csv", NAME)), row.names=FALSE)
