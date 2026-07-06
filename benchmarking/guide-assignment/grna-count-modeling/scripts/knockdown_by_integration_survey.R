## ============================================================================
## Extend the §7 phenotype channel (knockdown by integration class) to OTHER
## datasets. For each dataset with an aligned expression matrix + a resolvable
## guide->target-gene map, split cells by the guide's UMI count using GENERIC
## thresholds (strong >= STRONG / weak in WEAKLO..WEAKHI / count-1) -- a clean gap
## is NOT required (it is too rare once you also demand targeting+measured), and
## the phenotype question ("do low-count cells carry a functional guide?") does
## not need one. Measure the target gene's normalized expression vs the
## guide-absent baseline. Companion to scripts/replogle_knockdown_by_integration.R
## (the gap-based version on replogle's odm pipeline).
##
## Usage: Rscript scripts/knockdown_by_integration_survey.R <name> <sceptre_dir>
##   where <sceptre_dir> holds grna_matrix_aligned.rds, response_matrix.rds,
##   grna_target_data_frame.csv.
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix)})
args <- commandArgs(trailingOnly=TRUE); NAME <- args[1]; DIR <- path.expand(args[2])
OUT <- "results/global_ambient_poisson"
STRONG <- 30; WEAKLO <- 2; WEAKHI <- 7; MIN_STRONG <- 5; MIN_WEAK <- 5

g <- as(readRDS(file.path(DIR,"grna_matrix_aligned.rds")), "CsparseMatrix")
r <- as(readRDS(file.path(DIR,"response_matrix.rds")), "CsparseMatrix")
td <- read.csv(file.path(DIR,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
if(ncol(g)!=ncol(r)) stop("cell count mismatch")
lib <- as.numeric(Matrix::colSums(r))
tmap <- setNames(td$grna_target, td$grna_id)

## nonzero triplets of guide matrix
gv <- g@i+1L; cv <- rep.int(seq_len(ncol(g)), diff(g@p)); xv <- g@x
pos_by_g <- split(seq_along(gv), gv)

## NT-cell control mask: strong (>=30) carrier of a non-targeting guide. Canonical
## Perturb-seq baseline. Fall back to the complement (guide-absent) only when the
## dataset has no NTs (e.g. genome-wide tcell_cd4).
nt_ids <- td$grna_id[td$grna_target=="non-targeting"]
nt_rows <- match(nt_ids, rownames(g)); nt_rows <- nt_rows[!is.na(nt_rows)]
if(length(nt_rows)){
  ntmask <- Matrix::colSums(g[nt_rows, , drop=FALSE] >= 30) >= 1
  BASELINE <- "NT cells"
} else { ntmask <- rep(TRUE, ncol(g)); BASELINE <- "complement (no NTs in this dataset)" }
cat(sprintf("[%s] baseline = %s (%d control cells)\n", NAME, BASELINE, sum(ntmask)))

## pass 1: targeting guides with a resolvable, measured target AND enough
## strong + weak cells to be powered (generic count classes; no gap needed)
cand <- list()
for(gi in seq_len(nrow(g))){
  gid <- rownames(g)[gi]; tgt <- tmap[gid]
  if(is.na(tgt) || tgt=="non-targeting" || !(tgt %in% rownames(r))) next
  pos <- pos_by_g[[as.character(gi)]]; if(is.null(pos)) next
  vg <- xv[pos]
  if(sum(vg>=STRONG) < MIN_STRONG || sum(vg>=WEAKLO & vg<=WEAKHI) < MIN_WEAK) next
  cand[[length(cand)+1]] <- list(gi=gi, gid=gid, tgt=tgt, pos=pos)
}
cat(sprintf("[%s] targeting guides w/ measured target & >=%d strong & >=%d weak cells: %d\n",
            NAME, MIN_STRONG, MIN_WEAK, length(cand)))
if(!length(cand)){ cat("[",NAME,"] nothing to do\n"); quit(status=0) }

## pass 2: pull the needed target-gene rows (sparse), compute knockdown per guide
utg <- unique(vapply(cand, function(z) z$tgt, character(1)))
Es <- r[utg, , drop=FALSE]; rownames(Es) <- utg
E <- as.matrix(Es)                                       # genes(needed) x cells
rows <- vector("list", length(cand))
for(k in seq_along(cand)){ z <- cand[[k]]
  gc <- numeric(ncol(g)); gc[cv[z$pos]] <- xv[z$pos]     # this guide's counts (dense)
  cp <- E[z$tgt, ]/lib*1e4
  base <- mean(cp[ntmask & gc==0]); if(!is.finite(base)||base<=0) next
  strong <- gc>=STRONG; weak <- gc>=WEAKLO & gc<=WEAKHI; one <- gc==1
  mn <- c(strong=if(sum(strong)) mean(cp[strong]) else NA,
          weak  =if(sum(weak))   mean(cp[weak])   else NA,
          one   =if(sum(one))    mean(cp[one])    else NA)
  rows[[k]] <- data.frame(dataset=NAME, guide=z$gid, target=z$tgt,
    base_cp10k=round(base,3), n_weak=sum(weak), n_strong=sum(strong),
    kd_strong=round(mn["strong"]/base,3), kd_weak=round(mn["weak"]/base,3),
    kd_one=round(mn["one"]/base,3), row.names=NULL)
}
res <- do.call(rbind, rows[!vapply(rows,is.null,logical(1))])
write.csv(res, file.path(OUT, sprintf("knockdown_%s.csv", NAME)), row.names=FALSE)

pw <- res[res$base_cp10k>=0.5 & is.finite(res$kd_strong) & res$kd_strong<=0.7 &
          is.finite(res$kd_weak) & res$n_weak>=5, ]
cat(sprintf("[%s] power-positive guides: %d\n", NAME, nrow(pw)))
if(nrow(pw)){
  cat(sprintf("[%s] median kd  strong=%.2f  weak=%.2f  count1=%.2f  (n_weak-wtd weak=%.2f)\n",
    NAME, median(pw$kd_strong,na.rm=T), median(pw$kd_weak,na.rm=T), median(pw$kd_one,na.rm=T),
    weighted.mean(pw$kd_weak, pw$n_weak, na.rm=TRUE)))
}
