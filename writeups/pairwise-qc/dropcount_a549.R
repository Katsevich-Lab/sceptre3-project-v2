#!/usr/bin/env Rscript
# Drop fractions on a549 CRISPRi (high-MOI, complement control). Full trans grid:
# 53 targets x all genes. Compute observed n_nonzero_trt (for current filter) and the
# margin-based min attainable knockdown p (for nominal testability).
suppressMessages({ library(Matrix); library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
dir <- "/Users/ekatsevi/data/external/perturbseq-survey/a549_crispri_perturbseq_GSE236304/sceptre"
g <- readRDS(file.path(dir,"grna_matrix_aligned.rds")); r <- readRDS(file.path(dir,"response_matrix.rds"))
tdf <- read.csv(file.path(dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
N <- ncol(r); assigned <- g >= 5
tgts <- unique(tdf$grna_target[tdf$grna_target!="non-targeting"])
nonzero_r <- Matrix::rowSums(r > 0)                 # total non-zero cells per gene (all cells)
genes <- rownames(r)

## per (target, gene): observed non-zero trt count + margins
rows <- vector("list", length(tgts))
for (i in seq_along(tgts)) {
  gids <- tdf$grna_id[tdf$grna_target==tgts[i]]; gids <- gids[gids %in% rownames(g)]
  Tg <- which(Matrix::colSums(assigned[gids,,drop=FALSE]) > 0)   # treatment cells (complement = rest)
  n_trt <- length(Tg)
  nz_trt <- Matrix::rowSums(r[, Tg, drop=FALSE] > 0)             # observed non-zero trt per gene
  n_cntrl <- N - n_trt
  nz_cntrl <- nonzero_r - nz_trt                                 # non-zero complement per gene
  p_bar <- nz_cntrl / n_cntrl
  rows[[i]] <- data.table(target=tgts[i], n_trt=n_trt,
                          n_nonzero_trt=nz_trt, n_nonzero_cntrl=nz_cntrl,
                          p_bar=p_bar,
                          Ea = n_trt*p_bar,
                          minp_kd = (1-p_bar)^n_trt,
                          minp_2s = pmin(1, 2*pmin((1-p_bar)^n_trt, p_bar^n_trt)))
}
d <- rbindlist(rows)
cat(sprintf("a549: %s (target,gene) pairs | %d targets x %d genes | N=%d cells\n",
            format(nrow(d),big.mark=","), length(tgts), length(genes), N))
cat(sprintf("median n_trt=%.0f  median p_bar(gene)=%.3f\n", median(d$n_trt), median(d$p_bar)))

thr <- function(m) sprintf("%.1f%%", 100*mean(m))
cat("\n--- fraction of pairs THROWN OUT ---\n")
cat(sprintf(" current sceptre (obs<7 either side)          : %s\n", thr(d$n_nonzero_trt<7 | d$n_nonzero_cntrl<7)))
cat(sprintf(" nominal testability 1-sided knockdown  a=.10 : %s\n", thr(d$minp_kd > 0.10)))
cat(sprintf(" nominal testability 1-sided knockdown  a=.05 : %s\n", thr(d$minp_kd > 0.05)))
cat(sprintf(" nominal testability 2-sided            a=.10 : %s\n", thr(d$minp_2s > 0.10)))
cat(sprintf(" [ref] expected-count E[a] < 7                : %s\n", thr(d$Ea < 7)))

cur_drop <- d$n_nonzero_trt<7 | d$n_nonzero_cntrl<7; nom_keep <- d$minp_kd <= 0.10
cat(sprintf("\n current drops %s; of those nominal-1sided KEEPS %s; nominal drops but current keeps %s\n",
            thr(cur_drop), thr(nom_keep[cur_drop]), thr(!nom_keep & !cur_drop)))
saveRDS(d, file.path(scr,"a549_dropcount.rds")); cat("saved.\n")
