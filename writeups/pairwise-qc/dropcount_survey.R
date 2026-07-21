#!/usr/bin/env Rscript
# General drop-fraction computation for a perturbseq-survey dataset.
# Auto-detects MOI (=> control group), computes current-sceptre vs nominal-testability drops.
suppressMessages({ library(Matrix); library(data.table) })
ds  <- commandArgs(trailingOnly=TRUE)[1]
dir <- file.path("/Users/ekatsevi/data/external/perturbseq-survey", ds, "sceptre")
g <- readRDS(file.path(dir,"grna_matrix_aligned.rds")); r <- readRDS(file.path(dir,"response_matrix.rds"))
tdf <- read.csv(file.path(dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
N <- ncol(r); assigned <- g >= 5
nt_names <- c("non-targeting","non_targeting","NT","safe_harbor")
tdf$is_nt <- tolower(tdf$grna_target) %in% tolower(nt_names)
tgts <- unique(tdf$grna_target[!tdf$is_nt & tdf$grna_target!="" & !is.na(tdf$grna_target)])
moi <- median(Matrix::colSums(assigned))
lowmoi <- moi <= 1.5
set.seed(1); if (length(tgts) > 80) tgts <- sample(tgts, 80)   # cap targets for speed
nonzero_all <- Matrix::rowSums(r > 0)

## control-cell pool
nt_ids <- tdf$grna_id[tdf$is_nt]; nt_ids <- nt_ids[nt_ids %in% rownames(g)]
if (lowmoi && length(nt_ids) > 0) {
  NTcells <- which(Matrix::colSums(assigned[nt_ids,,drop=FALSE]) > 0)
  nz_ctrl_pool <- Matrix::rowSums(r[, NTcells, drop=FALSE] > 0); n_ctrl_pool <- length(NTcells)
  ctrl <- "NT cells"
} else { nz_ctrl_pool <- nonzero_all; n_ctrl_pool <- N; ctrl <- "complement" }

rows <- vector("list", length(tgts))
for (i in seq_along(tgts)) {
  gids <- tdf$grna_id[tdf$grna_target==tgts[i]]; gids <- gids[gids %in% rownames(g)]
  if (!length(gids)) next
  Tg <- which(Matrix::colSums(assigned[gids,,drop=FALSE]) > 0); n_trt <- length(Tg)
  if (n_trt < 1) next
  nz_trt <- Matrix::rowSums(r[, Tg, drop=FALSE] > 0)
  if (ctrl=="complement") { n_c <- N - n_trt; nz_c <- nonzero_all - nz_trt }
  else                    { n_c <- n_ctrl_pool; nz_c <- nz_ctrl_pool }
  p_bar <- pmin(0.999, nz_c / n_c)
  rows[[i]] <- data.table(n_trt=n_trt, n_nonzero_trt=nz_trt, n_nonzero_cntrl=nz_c,
                          minp_kd=(1-p_bar)^n_trt, minp_2s=pmin(1,2*pmin((1-p_bar)^n_trt, p_bar^n_trt)),
                          Ea=n_trt*p_bar)
}
d <- rbindlist(rows)
thr <- function(m) 100*mean(m)
cat(sprintf("%-34s MOI=%s ctrl=%-10s targets=%d genes=%d N=%d pairs=%s medNtrt=%.0f\n",
            ds, ifelse(lowmoi,"low","high"), ctrl, length(tgts), nrow(r), N,
            format(nrow(d),big.mark=","), median(d$n_trt)))
cat(sprintf("   THROWN OUT  current=%.1f%%  testable1sided(.10)=%.1f%%  (.05)=%.1f%%  testable2sided(.10)=%.1f%%  E[a]<7=%.1f%%\n",
            thr(d$n_nonzero_trt<7 | d$n_nonzero_cntrl<7), thr(d$minp_kd>.10), thr(d$minp_kd>.05),
            thr(d$minp_2s>.10), thr(d$Ea<7)))
cur <- d$n_nonzero_trt<7 | d$n_nonzero_cntrl<7
cat(sprintf("   of current's drops, nominal-1sided rescues %.1f%% ; nominal drops-but-current-keeps %.1f%%\n",
            thr((d$minp_kd<=.10)[cur]), thr(d$minp_kd>.10 & !cur)))
