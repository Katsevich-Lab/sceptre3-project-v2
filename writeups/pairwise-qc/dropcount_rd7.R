#!/usr/bin/env Rscript
# How many pairs would the NOMINAL TESTABILITY filter throw out, vs current sceptre?
# Replogle RD7 (genome-wide CRISPRi, low-MOI). Uses rho_g = n_trt/n_cntrl (data-driven,
# denominator-free): E[a] = rho_g * n_nonzero_cntrl; min attainable knockdown p = (1-p_bar)^n_trt.
suppressMessages({ library(arrow); library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
dir <- "/Users/ekatsevi/data/projects/sceptre3/nf_pipelines/rd7_trans/sceptre_outputs/trans_results"
rho <- as.data.table(readRDS(file.path(scr,"rd7_analyzed.rds"))$rho)[, .(grna_target, rho)]
setkey(rho, grna_target)

## data-driven control-cell count = max non-zero control count (ubiquitous gene ~ all NT cells)
n_cntrl <- open_dataset(dir) |> dplyr::summarise(m = max(n_nonzero_cntrl)) |> dplyr::collect() |> as.numeric()
cat(sprintf("n_cntrl (data-driven NT cells) = %d\n", n_cntrl))

## sample ~24 files spread across the 437 for a representative slice of the genome-wide grid
files <- sort(list.files(dir, "result_[0-9]+\\.parquet$", full.names=TRUE))
files <- files[seq(1, length(files), length.out=24)]
d <- as.data.table(open_dataset(files) |>
       dplyr::select(grna_target, n_nonzero_trt, n_nonzero_cntrl) |> dplyr::collect())
d <- d[rho, on="grna_target", nomatch=0]
d[, n_trt   := pmax(1, round(rho * n_cntrl))]
d[, p_bar   := pmin(0.999, n_nonzero_cntrl / n_cntrl)]
d[, Ea      := rho * n_nonzero_cntrl]                    # expected non-zero trt cells (denominator-free)
d[, minp_kd := (1 - p_bar)^n_trt]                        # min attainable one-sided (knockdown) p
d[, a_max   := pmin(n_trt, n_nonzero_trt + n_nonzero_cntrl)]
# min attainable two-sided p ~ 2*min(knockdown extreme, activation extreme)
d[, minp_2s := pmin(1, 2*pmin(minp_kd, p_bar^n_trt))]
cat(sprintf("sampled pairs: %s  (median n_trt=%.0f, median p_bar=%.3f)\n",
            format(nrow(d), big.mark=","), median(d$n_trt), median(d$p_bar)))

thr <- function(mask) sprintf("%.1f%%", 100*mean(mask))
cat("\n--- fraction of pairs THROWN OUT ---\n")
cat(sprintf(" current sceptre (obs n_nonzero_trt<7 or n_nonzero_cntrl<7): %s\n",
            thr(d$n_nonzero_trt < 7 | d$n_nonzero_cntrl < 7)))
cat(sprintf(" nominal testability, 1-sided knockdown  (minp_kd > 0.10)  : %s\n", thr(d$minp_kd > 0.10)))
cat(sprintf(" nominal testability, 1-sided knockdown  (minp_kd > 0.05)  : %s\n", thr(d$minp_kd > 0.05)))
cat(sprintf(" nominal testability, 2-sided            (minp_2s > 0.10)  : %s\n", thr(d$minp_2s > 0.10)))
cat(sprintf(" [ref] expected-count E[a] < 7 (my earlier fixed rule)     : %s\n", thr(d$Ea < 7)))

cat("\n--- overlap: of pairs current drops, how many does nominal-1sided KEEP? ---\n")
cur_drop <- d$n_nonzero_trt < 7 | d$n_nonzero_cntrl < 7
nom_keep <- d$minp_kd <= 0.10
cat(sprintf(" current drops %s of pairs; of those, nominal KEEPS %s (rescued)\n",
            thr(cur_drop), thr(nom_keep[cur_drop])))
cat(sprintf(" pairs nominal drops but current keeps: %s\n", thr(!nom_keep & !cur_drop)))
saveRDS(d[, .(grna_target, n_nonzero_trt, n_nonzero_cntrl, n_trt, p_bar, Ea, minp_kd, minp_2s)],
        file.path(scr,"rd7_dropcount.rds"))
cat("saved.\n")
