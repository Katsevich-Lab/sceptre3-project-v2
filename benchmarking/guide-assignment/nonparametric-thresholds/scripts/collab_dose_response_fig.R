## Clean Replogle dose-response figure for the collaborator writeup, with collaborator-facing
## labels (no "below-gap" / "NT baseline" jargon). Faithful to scripts/lowermode_dose_response.R
## block (B): per-count target knockdown for the low mode, vs the non-targeting baseline, with the
## integration-mode effect marked. Uses the real Replogle expression (ondisc).
suppressMessages({library(Matrix); library(ondisc); library(sceptre)})
OUT_SRC <- "results/global_ambient_poisson"; OUT <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Dm <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")

mc  <- as(readRDS(Dm), "CsparseMatrix")
resp <- initialize_odm_from_backing_file(file.path(Dp, "response.odm"))
so  <- readRDS(file.path(Dp, "sceptre_object.rds"))
lib <- exp(so@covariate_matrix[, "log(response_n_umis)"]); tdf <- so@grna_target_data_frame; gids <- rownames(resp)
ntg <- tdf$grna_id[tdf$grna_target == "non-targeting"]
ntmask <- rep(FALSE, ncol(mc)); for (g in ntg) ntmask <- ntmask | (as.numeric(mc[g, ]) >= 30)

pg <- read.csv(file.path(OUT_SRC, "perguide_replogle_rd7.csv"), stringsAsFactors = FALSE)
pg <- pg[pg$gap_lo >= 2, ]
pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target != "non-targeting" & pg$target %in% gids, ]

cnt <- c(); rat <- c(); sigr <- c()
for (i in seq_len(nrow(pg))) {
  r <- match(pg$guide[i], rownames(mc)); v <- as.numeric(mc[r, ]); cp <- as.numeric(resp[pg$target[i], ]) / lib * 1e4
  base <- mean(cp[ntmask & v == 0]); gh <- pg$gap_hi[i]; st <- which(v >= gh)
  if (base < 0.5 || !length(st) || mean(cp[st]) / base > 0.7) next          # power-positive guides only
  bg <- which(v >= 1 & v <= pg$gap_lo[i])
  cnt <- c(cnt, v[bg]); rat <- c(rat, cp[bg] / base); sigr <- c(sigr, cp[st] / base)
}
brk <- list(`1`=1, `2`=2, `3`=3, `4`=4, `5`=5, `6-7`=6:7, `8-11`=8:11, `12-19`=12:19, `20-40`=20:40)
xs <- c(); ys <- c()
for (nm in names(brk)) { s <- cnt %in% brk[[nm]]; if (sum(s) < 20) next; xs <- c(xs, mean(cnt[s])); ys <- c(ys, mean(rat[s])) }
sig <- mean(sigr)

png(file.path(OUT, "dose_response.png"), width = 1050, height = 680, res = 140)
par(mar = c(4.6, 4.9, 3.6, 1), mgp = c(2.9, 0.8, 0), cex.axis = 1.0, cex.lab = 1.1)
plot(xs, ys, log = "x", type = "b", pch = 16, col = "firebrick", cex = 1.4, ylim = c(0, 1.08),
     xlim = c(1, max(xs) * 60),
     xlab = "gRNA UMI count (the low mode, below the integration mode)",
     ylab = "target expression / non-targeting baseline",
     main = "Dose-response: knockdown deepens with the low-mode count\n(real low-dose integrations, approaching the integration effect)")
abline(h = 1, col = "grey45", lty = 2)
text(xs, ys, sprintf("%.2f", ys), pos = 3, cex = 0.72, col = "firebrick")
points(max(xs) * 30, sig, pch = 17, col = "steelblue4", cex = 1.8)
text(max(xs) * 30, sig, sprintf("integration\nmode %.2f", sig), pos = 3, cex = 0.82, col = "steelblue4")
dev.off()
cat("dose-response counts/knockdown:\n"); print(data.frame(count = round(xs,1), knockdown = round(ys,3)))
cat(sprintf("integration-mode effect: %.3f\n", sig))
cat("wrote", file.path(OUT, "dose_response.png"), "\n")
