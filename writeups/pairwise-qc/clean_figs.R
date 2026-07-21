suppressMessages({ library(ggplot2); library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
fig <- file.path(scr, "figs")
dir.create(fig, showWarnings = FALSE)
th  <- theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = "grey85"),
        strip.text = element_text(size = 10.5),
        legend.position = "top", legend.title = element_blank(),
        plot.margin = margin(6,10,4,6))
col_old <- "#C4471A"; col_new <- "#0B5CAD"

## FIG1 clean — mechanism
both <- as.data.table(readRDS(file.path(scr, "sim_qc_res_both.rds")))
both[, z := sign(log_2_fold_change) * qnorm(pmin(pmax(1 - p_value/2, 1e-6), 1-1e-9))]
nl <- both[kappa == 0 & is.finite(z)]
c1 <- cor(nl$n_nonzero_trt, nl$z); c2 <- cor(nl$m_expected, nl$z)
m1 <- melt(nl[, .(z, a = n_nonzero_trt, b = m_expected)], id.vars="z", variable.name="stat", value.name="value")
m1[, stat := factor(stat, levels=c("a","b"), labels=c(
  sprintf("current:  observed non-zero trt count   (r = %+.2f)", c1),
  sprintf("proposed:  expected count  n_trt · n_nonzero_tot / N   (r = %+.2f)", c2)))]
p1 <- ggplot(m1, aes(value, z)) +
  geom_point(alpha=0.16, size=0.7, colour="grey25") +
  geom_smooth(method="lm", se=FALSE, colour=col_old, linewidth=1) +
  geom_vline(xintercept=7, linetype=2, colour=col_new) +
  geom_hline(yintercept=0, colour="grey65", linewidth=0.3) +
  facet_wrap(~stat, scales="free_x") +
  labs(x="filter statistic value", y="signed test statistic z  (null pairs)") + th
ggsave(file.path(fig,"fig1_mechanism.png"), p1, width=10, height=4.3, dpi=150)

## FIG2 clean — retention
te <- both[panel=="grid" & kappa>0]
ret <- melt(te[, .(current=mean(pass_old), margin=mean(pass_new)), by=.(kappa,p,n_trt)],
            id.vars=c("kappa","p","n_trt"), variable.name="filter", value.name="retention")
ret[, filter := factor(filter, levels=c("current","margin"),
      labels=c("current  (observed count ≥ 7)","proposed  (expected count ≥ 7)"))]
ret[, p_lab := paste0("baseline expr. ", scales::percent(p,1))]
ret[, nt_lab := paste0("n_trt = ", n_trt)]
p2 <- ggplot(ret, aes(kappa, retention, colour=filter)) +
  geom_line(linewidth=1) + geom_point(size=1.9) +
  facet_grid(p_lab ~ nt_lab) +
  scale_colour_manual(values=c(col_old,col_new)) +
  scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
  labs(x="true knockdown strength  κ", y="true-effect pairs retained by QC") + th
ggsave(file.path(fig,"fig2_retention.png"), p2, width=8.6, height=5.6, dpi=150)

## FIGR1 clean — a549 entanglement
d <- as.data.table(readRDS(file.path(scr,"a549_res_C80.rds")))
nla <- d[type=="null" & is.finite(p_value)]; nla[, z := qnorm(pmin(pmax(p_value,1e-6),1-1e-6))]
brks <- c(-1,2,4,6,9,14,25,1e6); labs7 <- c("0-2","3-4","5-6","7-9","10-14","15-25","26+")
nla[, trt_bin := cut(n_nonzero_trt, brks, labels=labs7)]
nla[, m_bin  := cut(m_expected,   brks, labels=labs7)]
e1 <- nla[!is.na(trt_bin), .(z=mean(z), se=sd(z)/sqrt(.N)), by=trt_bin][order(trt_bin)]
e2 <- nla[!is.na(m_bin),  .(z=mean(z), se=sd(z)/sqrt(.N)), by=m_bin][order(m_bin)]
setnames(e1,"trt_bin","bin"); setnames(e2,"m_bin","bin")
e1[, stat:="current:  observed non-zero trt count"]; e2[, stat:="proposed:  expected count"]
ent <- rbind(e1,e2)
ca1 <- cor(nla$n_nonzero_trt, nla$z); ca2 <- cor(nla$m_expected, nla$z)
pA <- ggplot(ent, aes(bin, z, colour=stat, group=stat)) +
  geom_hline(yintercept=0, colour="grey65") +
  geom_vline(xintercept=3.5, linetype=2, colour="grey55") +
  geom_line(linewidth=1) + geom_pointrange(aes(ymin=z-2*se, ymax=z+2*se), fatten=2.4) +
  annotate("text", x=1.5, y=-0.75, label="dropped by\ncurrent filter", size=3.1, colour="grey35", lineheight=.9) +
  scale_colour_manual(values=c(col_old,col_new)) +
  labs(x="filter statistic value (binned)", y="mean signed z  (null pairs)") + th
ggsave(file.path(fig,"figR1_a549_entanglement.png"), pA, width=8.4, height=4.4, dpi=150)
cat("clean figs written. sim r:", round(c1,2), round(c2,2), "| a549 r:", round(ca1,3), round(ca2,3), "\n")
