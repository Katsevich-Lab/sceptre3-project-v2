suppressMessages({ library(ggplot2); library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
fig <- file.path(scr, "figs")
d <- as.data.table(readRDS(file.path(scr, "a549_res_C80.rds")))
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.background = element_rect(fill = "grey92"))
col_old <- "#D55E00"; col_new <- "#0072B2"
nl <- d[type == "null" & is.finite(p_value)]
nl[, z := qnorm(pmin(pmax(p_value, 1e-6), 1 - 1e-6))]   # left-sided signed z

## entanglement: mean signed z vs filter-stat bins
nl[, trt_bin := cut(n_nonzero_trt, c(-1,2,4,6,9,14,25,1e6),
      labels=c("0-2","3-4","5-6","7-9","10-14","15-25","26+"))]
nl[, m_bin := cut(m_expected, c(-1,2,4,6,9,14,25,1e6),
      labels=c("0-2","3-4","5-6","7-9","10-14","15-25","26+"))]
e1 <- nl[, .(z=mean(z), se=sd(z)/sqrt(.N), n=.N), by=trt_bin][!is.na(trt_bin)][order(trt_bin)]
e2 <- nl[, .(z=mean(z), se=sd(z)/sqrt(.N), n=.N), by=m_bin][!is.na(m_bin)][order(m_bin)]
e1[, stat := "current filter: observed n_nonzero_trt"]; setnames(e1,"trt_bin","bin")
e2[, stat := "margin filter: expected m"]; setnames(e2,"m_bin","bin")
ent <- rbind(e1,e2)
c1 <- cor(nl$n_nonzero_trt, nl$z); c2 <- cor(nl$m_expected, nl$z)
pA <- ggplot(ent, aes(bin, z, colour=stat, group=stat)) +
  geom_hline(yintercept=0, colour="grey60") +
  geom_vline(xintercept=3.5, linetype=2, colour="grey40") +
  geom_line(linewidth=1) + geom_pointrange(aes(ymin=z-2*se, ymax=z+2*se)) +
  annotate("text", x=1.4, y=Inf, vjust=1.4, label="dropped by\ncurrent filter (<7)", size=3, colour="grey30") +
  scale_colour_manual(values=c(col_old,col_new)) +
  labs(x="filter statistic value (binned)", y="mean signed test statistic z\n(null pairs; negative = repression)",
       colour=NULL,
       title="Real data (a549 CRISPRi): the current filter is entangled with the test statistic",
       subtitle=sprintf("Null pairs with FEW observed nonzero treatment cells carry systematically negative z. cor(n_nonzero_trt,z)=%+.2f vs cor(m,z)=%+.2f.\nThe current threshold (dashed) preferentially removes the repression tail; the margin filter shows no such gradient.", c1, c2)) +
  th + theme(axis.text.x=element_text(size=9))
ggsave(file.path(fig,"figR1_a549_entanglement.png"), pA, width=9, height=5, dpi=140)

## calibration QQ (null), under both filters
qq <- function(pv,lab){pv<-sort(pv[is.finite(pv)]); data.table(theo=ppoints(length(pv)),obs=pv,set=lab)}
q <- rbindlist(list(qq(nl$p_value,"all null"),
                    qq(nl[pass_old==TRUE]$p_value,"kept by current filter"),
                    qq(nl[pass_new==TRUE]$p_value,"kept by margin filter")))
pB <- ggplot(q, aes(-log10(theo),-log10(obs),colour=set)) +
  geom_abline(slope=1,colour="grey60") + geom_line(linewidth=0.9) +
  scale_colour_manual(values=c("all null"="grey30","kept by current filter"=col_old,"kept by margin filter"=col_new)) +
  labs(x=expression(-log[10]~"expected p"), y=expression(-log[10]~"observed p"), colour=NULL,
       title="Real data: margin filter preserves null calibration",
       subtitle="a549 CRISPRi null (trans) pairs; both filters leave the retained nulls uniform.") + th
ggsave(file.path(fig,"figR2_a549_calibration.png"), pB, width=6, height=5, dpi=140)

cat(sprintf("a549 real-data: cor(n_nonzero_trt,z)=%+.3f cor(m,z)=%+.3f\n", c1,c2))
cat(sprintf("null mean z among n_nonzero_trt<7: %+.3f ; >=7: %+.3f\n",
            nl[n_nonzero_trt<7, mean(z)], nl[n_nonzero_trt>=7, mean(z)]))
cat("figures: figR1_a549_entanglement.png, figR2_a549_calibration.png\n")
