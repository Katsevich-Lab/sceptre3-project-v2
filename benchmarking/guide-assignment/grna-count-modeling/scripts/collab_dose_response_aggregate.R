## A CAREFUL aggregate dose-response for the collaborator writeup (companion to the single-guide ENO1
## proof). Fixes the three flaws of the old dose_response.png:
##   (1) old strong-integration marker floated at an INVENTED x (max(xs)*30). Here the strong-
##       integration effect is a horizontal REFERENCE BAND (a y-level), tied to no x at all.
##   (2) old figure POOLED all cells by raw count -> pseudoreplication (high-cell guides dominate).
##       Here we average PER GUIDE first, then across guides (equal weight), with +/-SE across guides.
##   (3) old figure hid that high low-mode counts are only sampled by wide-gap guides. Here each
##       low-mode bin uses only a guide's OWN low-mode cells (count <= its gap_lo) and we annotate
##       how many guides contribute to each bin, so the selection is explicit.
## Same Replogle data + normalization + power-positive filter as collab_dose_response_fig.R.
suppressMessages({library(Matrix); library(ondisc); library(sceptre); library(ggplot2)})
OUT_SRC <- "results/global_ambient_poisson"; OUT <- "results/collaborator_writeup"
Dm <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")

mc  <- as(readRDS(Dm), "CsparseMatrix")
resp <- initialize_odm_from_backing_file(file.path(Dp, "response.odm"))
so  <- readRDS(file.path(Dp, "sceptre_object.rds"))
lib <- exp(so@covariate_matrix[, "log(response_n_umis)"]); tdf <- so@grna_target_data_frame; gids <- rownames(resp)
ntg <- tdf$grna_id[tdf$grna_target == "non-targeting"]
ntmask <- rep(FALSE, ncol(mc)); for (g in ntg) ntmask <- ntmask | (as.numeric(mc[g, ]) >= 30)

pg <- read.csv(file.path(OUT_SRC, "perguide_replogle_rd7.csv"), stringsAsFactors = FALSE)
pg <- pg[pg$gap_lo >= 2, ]; pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target != "non-targeting" & pg$target %in% gids, ]

# low-mode count bins, capped at ambient-scale counts. We deliberately STOP at 8: higher low-mode
# counts are sampled only by a handful of wide-gap guides, where a count of 15-30 is a moderate real
# integration (still below the gap, but not an ambiguous ambient-scale count), so those few-guide bins
# would do unearned work. Every retained count is far below any guide's integration onset (gap_hi>=56).
bins <- list(`1`=1, `2`=2, `3`=3, `4`=4, `5`=5, `6-8`=6:8)
MIN_CELLS_PER_GUIDE <- 3        # a guide contributes to a bin only with >= this many low-mode cells there
MIN_GUIDES_PER_BIN  <- 20       # a bin is shown only if >= this many guides contribute (keeps it broad)
MIN_SIG_CELLS       <- 3        # min integration cells for a guide to contribute to the band

# accumulate per-guide means: perbin[[bin]] = vector of per-guide mean(rel); sigvec = per-guide integration mean
perbin <- setNames(vector("list", length(bins)), names(bins)); sigvec <- c(); nguide <- 0
for (i in seq_len(nrow(pg))) {
  r <- match(pg$guide[i], rownames(mc)); v <- as.numeric(mc[r, ])
  cp <- as.numeric(resp[pg$target[i], ]) / lib * 1e4
  base <- mean(cp[ntmask & v == 0]); gh <- pg$gap_hi[i]; lo <- pg$gap_lo[i]; st <- which(v >= gh)
  if (base < 0.5 || !length(st) || mean(cp[st]) / base > 0.7) next          # power-positive guides only
  nguide <- nguide + 1
  rel <- cp / base
  if (length(st) >= MIN_SIG_CELLS) sigvec <- c(sigvec, mean(rel[st]))       # per-guide integration effect
  for (nm in names(bins)) {                                                  # LOW MODE only: v in bin AND v<=lo
    idx <- which(v %in% bins[[nm]] & v <= lo)
    if (length(idx) >= MIN_CELLS_PER_GUIDE) perbin[[nm]] <- c(perbin[[nm]], mean(rel[idx]))
  }
}
gm <- function(x) exp(mean(log(x)))                                          # geometric mean count for the bin x-pos
agg <- do.call(rbind, lapply(names(bins), function(nm) {
  m <- perbin[[nm]]; if (length(m) < MIN_GUIDES_PER_BIN) return(NULL)
  data.frame(bin = nm, x = min(bins[[nm]]), n_guides = length(m),
             mean = mean(m), se = sd(m) / sqrt(length(m)))
}))
sig_mean <- mean(sigvec); sig_se <- sd(sigvec) / sqrt(length(sigvec))
cat(sprintf("power-positive guides: %d ; integration-band from %d guides\n", nguide, length(sigvec)))
cat("per-bin (per-guide-averaged):\n"); print(agg, row.names = FALSE)
cat(sprintf("strong-integration level: %.3f +/- %.3f (across %d guides)\n", sig_mean, sig_se, length(sigvec)))

## ---- plot: low-mode dose-response (points +/-SE, n_guides annotated) + integration reference band --
xr <- range(agg$x)
p <- ggplot(agg, aes(x, mean)) +
  # strong-integration reference band (a LEVEL, spanning the whole width -- tied to no count)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = sig_mean - sig_se, ymax = sig_mean + sig_se,
           fill = "#D55E00", alpha = 0.15) +
  geom_hline(yintercept = sig_mean, colour = "#D55E00", linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = xr[2], y = sig_mean, hjust = 1, vjust = -0.6, size = 4.4, colour = "#B5490A",
           label = sprintf("strong-integration level (%.2f)", sig_mean)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  annotate("text", x = xr[2], y = 1, hjust = 1, vjust = -0.6, size = 4.2, colour = "grey40",
           label = "no-guide baseline") +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.14, colour = "#0072B2") +
  geom_point(colour = "#0072B2", size = 3) +
  geom_text(aes(label = paste0(n_guides, " guides"), y = mean + se), vjust = -0.9, size = 3.5, colour = "grey35") +
  scale_x_continuous(breaks = agg$x, labels = agg$bin) +
  coord_cartesian(xlim = c(xr[1] - 0.4, xr[2] + 0.4), ylim = c(min(sig_mean - sig_se, 0) , 1.12)) +
  labs(x = "gRNA UMI count in the cell (low mode, below each guide's gap)",
       y = "target expression / no-guide baseline\n(mean over guides ± SE)") +
  theme_bw(base_size = 15)
ggsave(file.path(OUT, "dose_response_aggregate.png"), p, width = 9.2, height = 5.4, dpi = 130)
cat("wrote", file.path(OUT, "dose_response_aggregate.png"), "\n")
