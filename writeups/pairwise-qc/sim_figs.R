#!/usr/bin/env Rscript
suppressMessages({ library(ggplot2); library(data.table) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
fig <- file.path(scr, "figs"); dir.create(fig, showWarnings = FALSE)
th  <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"), legend.position = "top")
col_old <- "#D55E00"; col_new <- "#0072B2"

both <- as.data.table(readRDS(file.path(scr, "sim_qc_res_both.rds")))
left <- as.data.table(readRDS(file.path(scr, "sim_qc_res_left.rds")))
# signed z from two-sided output
both[, z := sign(log_2_fold_change) * qnorm(pmin(pmax(1 - p_value/2, 1e-6), 1-1e-9))]

## ---------- Fig 1: MECHANISM — filter vs signed z under the null ----------
nl <- both[kappa == 0 & is.finite(z)]
c1 <- cor(nl$n_nonzero_trt, nl$z); c2 <- cor(nl$m_expected, nl$z)
m1 <- melt(nl[, .(z, `current filter:\nobserved n_nonzero_trt` = n_nonzero_trt,
                  `margin filter:\nexpected m = n_trt·p_ctrl` = m_expected)],
           id.vars = "z", variable.name = "stat", value.name = "value")
labs <- c(sprintf("current filter: observed n_nonzero_trt   (cor = %+.2f)", c1),
          sprintf("margin filter: expected m = n_trt·p_ctrl   (cor = %+.2f)", c2))
m1[, stat := factor(stat, labels = labs)]
p1 <- ggplot(m1, aes(value, z)) +
  geom_point(alpha = 0.18, size = 0.7, colour = "grey20") +
  geom_smooth(method = "lm", se = FALSE, colour = col_old, linewidth = 0.9) +
  geom_vline(xintercept = 7, linetype = 2, colour = col_new) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  facet_wrap(~stat, scales = "free_x") +
  labs(x = "filter statistic value", y = "signed test statistic z (null pairs)",
       title = "Under the null, the current filter is entangled with the test statistic; the margin filter is not",
       subtitle = "Each point is a null (κ=0) pair. Dashed line = threshold of 7. Removing pairs left of it removes the negative-z (repression) tail for the current filter only.") + th
ggsave(file.path(fig, "fig1_mechanism.png"), p1, width = 11, height = 4.6, dpi = 140)

## ---------- Fig 2: RETENTION of true effects vs knockdown strength ----------
te <- both[panel == "grid" & kappa > 0]
ret <- te[, .(old = mean(pass_old), new = mean(pass_new)), by = .(kappa, p, n_trt)]
retm <- melt(ret, id.vars = c("kappa","p","n_trt"), variable.name = "filter", value.name = "retention")
retm[, filter := factor(filter, levels = c("old","new"),
        labels = c("current (n_nonzero_trt≥7)","margin (m≥7)"))]
retm[, p_lab := paste0("baseline expr. ", scales::percent(p, 1))]
retm[, nt_lab := paste0("n_trt = ", n_trt)]
p2 <- ggplot(retm, aes(kappa, retention, colour = filter)) +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  facet_grid(p_lab ~ nt_lab) +
  scale_colour_manual(values = c(col_old, col_new)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "true knockdown strength κ", y = "fraction of true-effect pairs retained by QC",
       colour = NULL,
       title = "The current filter discards true effects in proportion to how strong they are",
       subtitle = "Margin filter retention is flat in κ (it correctly drops only the genuinely underpowered p=15%, n_trt=40 cell where expected m=6<7).") + th
ggsave(file.path(fig, "fig2_retention.png"), p2, width = 9, height = 6.2, dpi = 140)

## ---------- Fig 3: WHY — statistic value vs kappa ----------
sm <- te[, .(`observed n_nonzero_trt` = mean(n_nonzero_trt),
             `expected m (margin)` = mean(m_expected)), by = .(kappa, p, n_trt)]
sm0 <- both[panel=="grid" & kappa==0, .(`observed n_nonzero_trt`=mean(n_nonzero_trt),
             `expected m (margin)`=mean(m_expected)), by=.(kappa,p,n_trt)]
sm <- rbind(sm0, sm)
smm <- melt(sm, id.vars = c("kappa","p","n_trt"), variable.name = "stat", value.name = "val")
smm[, p_lab := paste0("baseline expr. ", scales::percent(p, 1))]
smm[, nt_lab := paste0("n_trt = ", n_trt)]
p3 <- ggplot(smm, aes(kappa, val, colour = stat)) +
  geom_hline(yintercept = 7, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 1.7) +
  facet_grid(p_lab ~ nt_lab, scales = "free_y") +
  scale_colour_manual(values = c("observed n_nonzero_trt" = col_old, "expected m (margin)" = col_new)) +
  labs(x = "true knockdown strength κ", y = "mean filter statistic", colour = NULL,
       title = "Mechanism: the knockdown depletes the observed treatment-cell count, but not the margin",
       subtitle = "Dashed line = threshold 7. observed n_nonzero_trt falls below 7 as κ grows; expected m stays put.") + th
ggsave(file.path(fig, "fig3_why.png"), p3, width = 9, height = 6.2, dpi = 140)

## ---------- Fig 4: CALIBRATION (null QQ) ----------
qq_dt <- function(pv, lab) { pv <- sort(pv[is.finite(pv)]);
  data.table(theo = ppoints(length(pv)), obs = pv, set = lab) }
nlb <- both[kappa == 0]
qq <- rbindlist(list(
  qq_dt(nlb$p_value, "all null pairs"),
  qq_dt(nlb[pass_old == TRUE]$p_value, "kept by current filter"),
  qq_dt(nlb[pass_new == TRUE]$p_value, "kept by margin filter")))
p4 <- ggplot(qq, aes(-log10(theo), -log10(obs), colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("all null pairs"="grey30",
      "kept by current filter"=col_old, "kept by margin filter"=col_new)) +
  labs(x = expression(-log[10]~"expected p (uniform)"),
       y = expression(-log[10]~"observed p"), colour = NULL,
       title = "Both filters preserve null calibration",
       subtitle = "The margin filter keeps the retained-null p-values uniform (it is exactly ancillary by construction).") + th
ggsave(file.path(fig, "fig4_calibration.png"), p4, width = 6.4, height = 5.2, dpi = 140)

## ---------- Fig 5: BOTTOM LINE — discoveries at BH<0.1 ----------
disc <- function(dt, side) {
  f <- function(sub) { sub <- sub[is.finite(p_value)]; padj <- p.adjust(sub$p_value,"BH")
    data.table(true = sum(padj<0.1 & sub$kappa>0), false = sum(padj<0.1 & sub$kappa==0)) }
  rbind(cbind(filter="current (n_nonzero_trt≥7)", f(dt[pass_old==TRUE])),
        cbind(filter="margin (m≥7)",              f(dt[pass_new==TRUE])))[, side := side]
}
db <- rbind(disc(both,"two-sided test"), disc(left,"left-sided test (CRISPRi)"))
dbm <- melt(db, id.vars = c("filter","side"), variable.name = "kind", value.name = "n")
dbm[, kind := factor(kind, levels=c("true","false"),
        labels=c("true discoveries","false discoveries"))]
p5 <- ggplot(dbm, aes(filter, n, fill = filter)) +
  geom_col(width = 0.65) + geom_text(aes(label = n), vjust = -0.3, size = 3.6) +
  facet_grid(kind ~ side, scales = "free_y") +
  scale_fill_manual(values = c(col_old, col_new), guide = "none") +
  labs(x = NULL, y = "number of pairs (BH-adjusted p < 0.1)",
       title = "Bottom line: the margin filter finds ~3× more true effects at the same controlled FDR",
       subtitle = "Realized FDR stays < 10% under both filters; the margin filter simply stops throwing away true positives.") +
  th + theme(axis.text.x = element_text(angle = 12, hjust = 1))
ggsave(file.path(fig, "fig5_discoveries.png"), p5, width = 9, height = 6, dpi = 140)

cat("Figures written to", fig, "\n")
cat(sprintf("Null-pair signed-z correlations: current n_nonzero_trt %+.3f | margin m %+.3f\n", c1, c2))
