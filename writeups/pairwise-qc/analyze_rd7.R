#!/usr/bin/env Rscript
# Robust Replogle RD7 analysis. Estimate each target's rho_g = n_trt/n_cntrl
# directly from trans (null) pairs (ratio of sums of n_nonzero_trt to n_nonzero_cntrl),
# so the margin filter m = rho_g * n_nonzero_cntrl needs NO external cell counts and
# is automatically consistent with the analysis's effective cell set.
suppressMessages({ library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
ex  <- readRDS(file.path(scr, "rd7_extract.rds"))
cis <- as.data.table(ex$cis); trans <- as.data.table(ex$trans)

## estimate rho_g from trans pairs (pooled ratio of sums; robust to sparse real trans effects)
rho <- trans[n_nonzero_cntrl > 0, .(rho = sum(n_nonzero_trt) / sum(as.numeric(n_nonzero_cntrl)),
                                    n = .N), by = grna_target]
cat("rho_g (= n_trt/n_cntrl) quantiles:\n"); print(quantile(rho$rho, c(0,.25,.5,.75,1)))

setkey(rho, grna_target)
addm <- function(d){
  d <- merge(d, rho[, .(grna_target, rho)], by = "grna_target", all.x = TRUE)
  d[, m_expected := rho * n_nonzero_cntrl]
  d[, pass_old := n_nonzero_trt >= 7 & n_nonzero_cntrl >= 7]
  d[, pass_new := m_expected    >= 7 & n_nonzero_cntrl >= 7]
  d
}
cis <- addm(cis); trans <- addm(trans)
cis <- cis[is.finite(p_value) & !is.na(rho)]
saveRDS(list(cis=cis, trans=trans, rho=rho), file.path(scr, "rd7_analyzed.rds"))

cat("\n==== CIS positive controls (real knockdowns), n =", nrow(cis), "====\n")
cat(sprintf("median log2FC = %.2f ; %% with log2FC<-1 : %.0f%%\n",
            median(cis$log_2_fold_change), 100*mean(cis$log_2_fold_change < -1)))

## cross-tab pass_old vs pass_new among cis
tab <- cis[, .N, by = .(pass_old, pass_new)][order(-pass_old, -pass_new)]
cat("\ncross-tab of cis pairs (pass_old x pass_new):\n"); print(tab)

## the key discordant cell: current DROPS but margin KEEPS
disc_cell <- cis[pass_old == FALSE & pass_new == TRUE]
cat(sprintf("\nCIS pairs DROPPED by current filter but KEPT by margin: %d\n", nrow(disc_cell)))
cat(sprintf("  their median log2FC = %.2f ; median n_nonzero_trt = %.0f ; median m = %.0f\n",
            median(disc_cell$log_2_fold_change), median(disc_cell$n_nonzero_trt), median(disc_cell$m_expected)))
cat(sprintf("  fraction individually significant (p<0.05): %.2f\n", mean(disc_cell$p_value < 0.05)))

## opposite discordant cell: current KEEPS but margin DROPS
disc_cell2 <- cis[pass_old == TRUE & pass_new == FALSE]
cat(sprintf("\nCIS pairs KEPT by current but DROPPED by margin: %d\n", nrow(disc_cell2)))
if(nrow(disc_cell2)>0) cat(sprintf("  their median log2FC = %.2f ; median n_nonzero_trt = %.0f ; median m = %.1f ; frac p<0.05 = %.2f\n",
            median(disc_cell2$log_2_fold_change), median(disc_cell2$n_nonzero_trt),
            median(disc_cell2$m_expected), mean(disc_cell2$p_value<0.05)))

## discovery counts (BH<0.1) among cis, each filter
disc <- function(sub){ p.adjust(sub$p_value,"BH") -> a; sum(a<0.1) }
cat(sprintf("\nCIS discoveries (BH<0.1): current=%d (of %d kept) | margin=%d (of %d kept)\n",
            disc(cis[pass_old==TRUE]), sum(cis$pass_old),
            disc(cis[pass_new==TRUE]), sum(cis$pass_new)))

## Stratify by n_trt (rho as proxy): does the filter flaw bite in the LOW-cell regime?
cis[, ntrt_est := rho * (n_nonzero_cntrl / p_hat_cntrl)]  # not needed; use rho tertiles
cis[, rho_grp := cut(rho, quantile(rho, c(0,1/3,2/3,1)), include.lowest=TRUE,
                     labels=c("few trt cells","medium","many trt cells"))]
strat <- cis[, .(n=.N,
                 cur_keep = mean(pass_old), mar_keep = mean(pass_new),
                 cur_disc = disc(.SD[pass_old==TRUE]), mar_disc = disc(.SD[pass_new==TRUE])),
             by = rho_grp][order(rho_grp)]
cat("\nStratified by treatment-cell abundance (rho tertile):\n"); print(strat)

cat("\n==== TRANS entanglement (null-ish) ====\n")
tr <- trans[is.finite(p_value) & !is.na(rho)]
tr[, z := sign(log_2_fold_change) * qnorm(pmin(pmax(1 - p_value/2, 1e-8), 1-1e-8))]
cat(sprintf("cor(n_nonzero_trt, z)=%+.3f | cor(m, z)=%+.3f\n",
            cor(tr$n_nonzero_trt, tr$z, use="complete.obs"), cor(tr$m_expected, tr$z, use="complete.obs")))
cat(sprintf("mean z | n_nonzero_trt<7 : %+.3f  vs  >=7 : %+.3f\n",
            tr[n_nonzero_trt<7, mean(z)], tr[n_nonzero_trt>=7, mean(z)]))
cat("DONE\n")
