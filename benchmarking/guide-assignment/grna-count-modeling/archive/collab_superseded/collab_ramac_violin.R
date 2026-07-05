## Violin plot of RAMAC target-gene expression across RAMAC gRNA-count bins, for the collaborator
## writeup. Same Replogle data + normalization as scripts/collab_dose_response_fig.R, but showing the
## FULL per-cell distribution of target expression at each count bin (not just the per-bin mean).
## Bins: 0, 1, 2, 3, 4, 5-10, 50+ (RAMAC's gap is 8-69, so 5-10 is the top of the low mode and 50+ is
## the integration/high mode). y = target expression relative to the count-0 (no-guide) baseline.
suppressMessages({library(Matrix); library(ondisc); library(sceptre); library(ggplot2)})
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Dm <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")

GUIDE <- "2828_RAMAC_P1P2_ENSG00000169612"
mc   <- as(readRDS(Dm), "CsparseMatrix")
resp <- initialize_odm_from_backing_file(file.path(Dp, "response.odm"))
so   <- readRDS(file.path(Dp, "sceptre_object.rds"))
lib  <- exp(so@covariate_matrix[, "log(response_n_umis)"])
tdf  <- so@grna_target_data_frame
target <- tdf$grna_target[match(GUIDE, tdf$grna_id)]
stopifnot(!is.na(target), target %in% rownames(resp))
cat(sprintf("RAMAC guide=%s  target=%s\n", GUIDE, target))

v  <- as.numeric(mc[GUIDE, ])                          # per-cell RAMAC gRNA UMI count
cp <- as.numeric(resp[target, ]) / lib * 1e4           # normalized target expression (per 10k)
base <- mean(cp[v == 0])                               # no-guide baseline
cat(sprintf("baseline (count-0) mean expr = %.3f per 10k\n", base))

## ---- bin the counts: 0, 1, 2, 3, 4, 5-10, 50+ -----------------------------------
bin_of <- function(x) {
  ifelse(x == 0, "0",
  ifelse(x == 1, "1",
  ifelse(x == 2, "2",
  ifelse(x == 3, "3",
  ifelse(x == 4, "4",
  ifelse(x >= 5 & x <= 10, "5-10",
  ifelse(x >= 50, "50+", NA)))))))
}
lev <- c("0", "1", "2", "3", "4", "5-10", "50+")
df <- data.frame(count = v, expr = cp, rel = cp / base)
df$bin <- bin_of(df$count)
df <- df[!is.na(df$bin), ]
df$bin <- factor(df$bin, levels = lev)

summ <- do.call(rbind, lapply(lev, function(b) {
  s <- df$rel[df$bin == b]
  data.frame(bin = b, n = length(s), mean = mean(s), median = median(s))
}))
cat("per-bin (relative to count-0 baseline):\n"); print(summ, row.names = FALSE)

## ---- violin ---------------------------------------------------------------------
p <- ggplot(df, aes(bin, rel)) +
  geom_violin(fill = "#0072B2", colour = "grey30", alpha = 0.55, scale = "width", linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", colour = "#D55E00", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_text(data = summ, aes(bin, y = -Inf, label = paste0("n=", n)),
            vjust = -0.6, size = 4, colour = "grey35", inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, quantile(df$rel, 0.995))) +
  labs(x = "RAMAC gRNA UMI count in the cell",
       y = "RAMAC target expression / no-guide baseline") +
  theme_bw(base_size = 16)
ggsave(file.path(OUT, "ramac_expression_violin.png"), p, width = 9, height = 5.2, dpi = 130)
cat("wrote", file.path(OUT, "ramac_expression_violin.png"), "\n")
