#!/usr/bin/env Rscript
## Analyze the DE-pushthrough operating-point grid (scripts/de_opgrid.R outputs).
## Central question: how much does the assignment OPERATING POINT (correction, q) move
## downstream sceptre calibration (neg-control t1e) and power, vs how much the assignment
## itself changes? Figures + summary -> results/error_control/de/.
## NOTE: hard-QC "singleton" ops are EXCLUDED (byte-identical to keep-all — unresolved anomaly,
## likely sceptre's low-MOI DE already drops multiplet cells; flagged, not concluded).
suppressPackageStartupMessages({library(ggplot2)})
DE <- "results/error_control/de"
d <- read.csv(file.path(DE,"de_opgrid_master.csv"), stringsAsFactors=FALSE)
d <- d[d$qc!="singleton",]                                   # exclude the anomalous hard-QC ops
d$dataset <- factor(d$dataset, levels=c("replogle","endoc","gasperini","a549","dctap_lowmoi","dctap_highmoi","cd8_tcell"))
d$op <- factor(d$op, levels=unique(d$op[order(d$corr, d$q)]))
has_cal <- !is.na(d$t1e_05)

## -------- summary: within-dataset spread of assignment vs DE --------
agg <- do.call(rbind, lapply(split(d, d$dataset), function(x){
  data.frame(dataset=x$dataset[1], moi=x$moi[1], n_ops=nrow(x),
    calls_fold = max(x$n_calls,na.rm=T)/min(x$n_calls,na.rm=T),
    t1e_min=min(x$t1e_05,na.rm=T), t1e_max=max(x$t1e_05,na.rm=T),
    t1e_span=max(x$t1e_05,na.rm=T)-min(x$t1e_05,na.rm=T),
    pow_min=min(x$pow_frac_sig,na.rm=T), pow_max=max(x$pow_frac_sig,na.rm=T),
    pow_span=max(x$pow_frac_sig,na.rm=T)-min(x$pow_frac_sig,na.rm=T)) }))
agg[] <- lapply(agg, function(c) if(is.numeric(c)) round(c,3) else c)
write.csv(agg, file.path(DE,"de_opgrid_summary.csv"), row.names=FALSE)
cat("=== per-dataset: assignment fold-change vs DE span (singleton excluded) ===\n"); print(agg, row.names=FALSE)

## -------- Fig 1: calibration t1e per operating point, faceted by dataset --------
dd <- d[has_cal,]
p1 <- ggplot(dd, aes(op, t1e_05, group=1)) +
  geom_hline(yintercept=0.05, colour="grey55", linetype=2) +
  geom_line(colour="grey70") + geom_point(aes(colour=corr, shape=corr), size=2.4) +
  facet_wrap(~dataset, scales="free_y", ncol=3) +
  scale_x_discrete(labels=function(x) sub("_none","",sub("_q","·q",x))) +
  labs(title="Downstream DE calibration is flat across the assignment operating point",
       subtitle="Negative-control Type-I error (t1e at 0.05) per operating point. Dashed = nominal 0.05. Within each dataset the level barely moves; between datasets it differs (confounding/NT-count/panel).",
       x=NULL, y="neg-control t1e (target 0.05)", colour="correction", shape="correction") +
  theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=55,hjust=1,size=6))
ggsave(file.path(DE,"de_opgrid_calibration_by_op.png"), p1, width=11, height=6.5, dpi=140)

## -------- Fig 2: assignment count varies; DE metrics don't (single faceted plot) --------
d2 <- rbind(
  transform(d[has_cal,], metric="neg-control t1e (target 0.05)", value=t1e_05),
  transform(d[!is.na(d$pow_frac_sig),], metric="power (frac sig)", value=pow_frac_sig))
p2 <- ggplot(d2, aes(n_calls, value, colour=dataset)) +
  geom_hline(data=data.frame(metric="neg-control t1e (target 0.05)", y=0.05), aes(yintercept=y),
             linetype=2, colour="grey55") +
  geom_line(aes(group=dataset), alpha=.4) + geom_point(size=1.9) +
  facet_wrap(~metric, scales="free_y") + scale_x_log10() +
  labs(title="The assignment changes across operating points; the downstream DE outcome does not",
       subtitle="Each line = one dataset swept across corrections + q. Calls move (x); calibration and power (y) stay put. Between-dataset differences dwarf within-dataset operating-point differences.",
       x="# assignment calls (log10)", y=NULL, colour="dataset") + theme_bw(base_size=11)
ggsave(file.path(DE,"de_opgrid_calls_vs_de.png"), p2, width=12, height=5, dpi=140)

## -------- Fig 3: BH q-sweep (the cost-asymmetry test) --------
qs <- d[d$corr=="BH" & d$qc=="none",]
qs2 <- rbind(
  transform(qs[!is.na(qs$t1e_05),], metric="neg-control t1e", value=t1e_05),
  transform(qs[!is.na(qs$pow_frac_sig),], metric="power (frac sig)", value=pow_frac_sig))
p3 <- ggplot(qs2, aes(q, value, colour=dataset)) +
  geom_hline(data=data.frame(metric="neg-control t1e",y=0.05), aes(yintercept=y), linetype=2, colour="grey55") +
  geom_line() + geom_point(size=1.8) + facet_wrap(~metric, scales="free_y") +
  scale_x_log10(breaks=c(.01,.05,.1,.2,.4)) +
  labs(title="BH q-sweep: loosening the FDR threshold barely changes downstream DE",
       subtitle="The cost-asymmetry lever would show on weak effects; here on-target power is saturated, so both calibration and detected power are ~flat 0.01→0.40.",
       x="BH q (log)", y=NULL) + theme_bw(base_size=11)
ggsave(file.path(DE,"de_opgrid_qsweep.png"), p3, width=11, height=4.5, dpi=140)

cat("\nwrote de_opgrid_summary.csv + 3 figures to", DE, "\n")
