## Single-guide expression violin that PROVES the low-mode-is-partly-real point, using a Replogle
## guide with a well-expressed target and a densely-populated low mode: ENO1 (glycolytic enzyme,
## baseline ~16.6 per 10k, 3% dropout; low mode counts 1-17 with hundreds of cells; integration-mode
## knockdown 0.10x). Per-cell target expression / no-guide baseline, by gRNA-count bin.
## Bins: 0, 1, 2, 3, 4, 5-10, 50+ (ENO1's gap is 18-62, so 5-10 is mid-low-mode and 50+ is the
## integration mode). Same Replogle data + normalization as collab_dose_response_fig.R.
suppressMessages({library(Matrix); library(ondisc); library(sceptre); library(ggplot2)})
source("scripts/datasets.R")   # load_replogle_rd7_de()
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GUIDE <- "2640_ENO1_P1P2_ENSG00000074800"; GENE <- "ENO1"
rd <- load_replogle_rd7_de(); mc <- rd$mc; resp <- rd$resp; so <- rd$so
lib  <- exp(so@covariate_matrix[, "log(response_n_umis)"])
tdf  <- so@grna_target_data_frame
target <- tdf$grna_target[match(GUIDE, tdf$grna_id)]
stopifnot(!is.na(target), target %in% rownames(resp))

v  <- as.numeric(mc[GUIDE, ])
cp <- as.numeric(resp[target, ]) / lib * 1e4
base <- mean(cp[v == 0])
cat(sprintf("%s guide=%s target=%s  baseline=%.2f per10k\n", GENE, GUIDE, target, base))

bin_of <- function(x) ifelse(x <= 3, "0-3", ifelse(x == 4, "4",
  ifelse(x >= 5 & x <= 10, "5-10", ifelse(x >= 50, "50+", NA))))
lev <- c("0-3", "4", "5-10", "50+")
df <- data.frame(count = v, rel = cp / base); df$bin <- bin_of(df$count)
df <- df[!is.na(df$bin), ]; df$bin <- factor(df$bin, levels = lev)

summ <- do.call(rbind, lapply(lev, function(b) { s <- df$rel[df$bin == b]
  data.frame(bin = b, n = length(s), mean = round(mean(s), 3), median = round(median(s), 3)) }))
cat("per-bin (relative to count-0 baseline):\n"); print(summ, row.names = FALSE)

ytop <- as.numeric(quantile(df$rel, 0.99))
p <- ggplot(df, aes(bin, rel)) +
  geom_violin(fill = "#0072B2", colour = "grey30", alpha = 0.55, scale = "width", linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", colour = "#D55E00", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_text(data = summ, aes(bin, y = ytop, label = paste0("n=", n)),
            vjust = 1, size = 4, colour = "grey35", inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, ytop)) +
  labs(x = sprintf("%s gRNA UMI count in the cell", GENE),
       y = sprintf("%s target expression / no-guide baseline", GENE)) +
  theme_bw(base_size = 16)
ggsave(file.path(OUT, "eno1_expression_violin.png"), p, width = 9, height = 5.2, dpi = 130)
cat("wrote", file.path(OUT, "eno1_expression_violin.png"), "\n")
