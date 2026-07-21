#!/usr/bin/env Rscript
# Evaluate the "testability" filter: given the fixed margins, take the MOST FAVORABLE
# 2x2 table and its minimum attainable Fisher/hypergeometric p-value; drop the pair
# if even that best case can't reach significance. (This is Tarone 1990 / significant-
# pattern-mining "testability".) Compare to the expected-count filter and the current one.
suppressMessages(library(data.table))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
d <- as.data.table(readRDS(file.path(scr, "sim_qc_res_both.rds")))

## recover margins
d[, n_cntrl := round(n_nonzero_cntrl / p_hat_cntrl)]
d[, N   := n_trt + n_cntrl]
d[, m1  := n_nonzero_trt + n_nonzero_cntrl]      # column margin: total non-zero (a+c)
d[, m0  := N - m1]                               # total zero

## minimum attainable ONE-SIDED (knockdown) p-value = most extreme table:
## put as few non-zero cells in treatment as possible: a_min = max(0, n_trt - m0)
d[, a_min := pmax(0L, n_trt - m0)]
d[, p_min := phyper(a_min, m1, m0, n_trt)]       # P(A <= a_min) at the extreme
## expected non-zero treatment cells (Cochran / my earlier proposal)
d[, E_a := n_trt * m1 / N]

cat("=== p_min  vs  exp(-E[a]) approximation ===\n")
d2 <- d[is.finite(p_min) & p_min>0]
cat(sprintf("median |log p_min - (-E_a)| : %.3f  (approx p_min ~ exp(-E[a]))\n",
            median(abs(log(d2$p_min) + d2$E_a))))
cat(sprintf("current thresh 7 on E[a]  <=>  p_min <= %.2e\n", exp(-7)))

## three filters (control side always easier; treatment side binds)
delta_nom  <- 0.05
delta_impl <- exp(-7)                             # matches current 'expected count >= 7'
delta_bonf <- 0.05 / nrow(d)                      # multiplicity-aware (Tarone-style)
d[, pass_current  := n_nonzero_trt >= 7 & n_nonzero_cntrl >= 7]
d[, pass_expected := E_a >= 7 & (n_cntrl*m1/N) >= 7]
d[, pass_test_nom  := p_min <= delta_nom  & n_nonzero_cntrl > 0]
d[, pass_test_impl := p_min <= delta_impl & n_nonzero_cntrl > 0]
d[, pass_test_bonf := p_min <= delta_bonf & n_nonzero_cntrl > 0]

## null calibration (ancillary filter must keep uniform)
nl <- d[kappa==0 & is.finite(p_value)]
cal <- function(m) sprintf("n=%4d  P(p<.05)=%.3f", sum(m), mean(nl$p_value[m] < .05))
cat("\n=== null calibration under each filter ===\n")
cat(" current      :", cal(nl$pass_current), "\n")
cat(" expected>=7  :", cal(nl$pass_expected), "\n")
cat(" testable@.05 :", cal(nl$pass_test_nom), "\n")
cat(" testable@e^-7:", cal(nl$pass_test_impl), "\n")
cat(" testable@Bonf:", cal(nl$pass_test_bonf), "\n")

## retention of TRUE effects by knockdown strength
te <- d[panel=="grid" & kappa>0]
cat("\n=== retention of true effects (mean pass rate) by kappa ===\n")
ag <- te[, .(current=mean(pass_current), expected=mean(pass_expected),
             test_impl=mean(pass_test_impl), test_bonf=mean(pass_test_bonf)), by=kappa][order(kappa)]
print(ag)

## discoveries (BH<0.1) within each filtered set
disc <- function(mask){ s <- d[mask & is.finite(p_value)]
  c(true=sum(p.adjust(s$p_value,"BH")<.1 & s$kappa>0), false=sum(p.adjust(s$p_value,"BH")<.1 & s$kappa==0)) }
cat("\n=== discoveries (true/false, BH<0.1) ===\n")
for(nm in c("pass_current","pass_expected","pass_test_impl","pass_test_bonf"))
  cat(sprintf(" %-16s tested=%4d  true=%3d false=%3d\n", nm, sum(d[[nm]]), disc(d[[nm]])[1], disc(d[[nm]])[2]))

## agreement between expected-count and testability@e^-7
cat(sprintf("\nexpected>=7 vs testable@e^-7 disagree on %d / %d pairs\n",
            sum(d$pass_expected != d$pass_test_impl), nrow(d)))
