## Search Replogle clean-gap targeting guides for one that makes a GOOD single-guide expression
## violin: (a) well-expressed target (high baseline -> per-cell expression not all zeros), (b) a
## populated low mode (cells across counts 2,3,4,5+ so per-count bins aren't tiny), (c) clear
## knockdown at the integration mode (a real hit). Prints a ranked candidate table.
suppressMessages({library(Matrix); library(ondisc); library(sceptre)})
source("scripts/datasets.R")   # load_replogle_rd7_de()
OUT_SRC <- "results/global_ambient_poisson"

rd <- load_replogle_rd7_de(); mc <- rd$mc; resp <- rd$resp; so <- rd$so
lib  <- exp(so@covariate_matrix[, "log(response_n_umis)"]); gids <- rownames(resp)
tdf  <- so@grna_target_data_frame

pg <- read.csv(file.path(OUT_SRC, "perguide_replogle_rd7.csv"), stringsAsFactors = FALSE)
pg$target <- tdf$grna_target[match(pg$guide, tdf$grna_id)]
pg <- pg[!is.na(pg$target) & pg$target != "non-targeting" & pg$target %in% gids, ]
# prefilter to reasonably wide low modes so per-count bins can be populated
pg <- pg[pg$gap_lo >= 3 & pg$guide %in% rownames(mc), ]
cat(sprintf("%d clean-gap targeting guides with gap_lo>=3 to evaluate\n", nrow(pg)))

rows <- list()
for (i in seq_len(nrow(pg))) {
  g <- pg$guide[i]; tgt <- pg$target[i]; lo <- pg$gap_lo[i]; hi <- pg$gap_hi[i]
  v  <- as.numeric(mc[g, ])
  cp <- as.numeric(resp[tgt, ]) / lib * 1e4
  base <- mean(cp[v == 0]); if (!is.finite(base) || base <= 0) next
  st <- which(v >= hi); kd <- if (length(st)) mean(cp[st]) / base else NA   # integration-mode knockdown ratio
  # low-mode population across counts 2..lo and the 5-10 bin
  n2 <- sum(v == 2); n3 <- sum(v == 3); n4 <- sum(v == 4); n5_10 <- sum(v >= 5 & v <= 10)
  n_low <- sum(v >= 2 & v <= lo)                       # weak-integration region
  frac0 <- mean(cp[v == 0] == 0)                        # dropout fraction at baseline (lower = better)
  rows[[length(rows)+1]] <- data.frame(
    guide = g, target = tgt, gap_lo = lo, gap_hi = hi,
    baseline = round(base, 2), dropout0 = round(frac0, 2), knockdown = round(kd, 2),
    n_sig = length(st), n2 = n2, n3 = n3, n4 = n4, n5_10 = n5_10, n_low = n_low)
}
df <- do.call(rbind, rows)

# "good violin" = well-expressed (high baseline, low dropout), populated low mode (cells at 3 and 4
# too, decent n_low), and a real hit (knockdown <= ~0.5). Rank by a simple composite.
cand <- df[df$baseline >= 3 & df$knockdown <= 0.6 & df$n3 >= 5 & df$n4 >= 2 & df$n_sig >= 30, ]
cand <- cand[order(-cand$baseline * (1 - cand$knockdown) * pmin(cand$n_low, 200)), ]
cat("\n===== TOP CANDIDATES (well-expressed target + populated low mode + real knockdown) =====\n")
print(head(cand, 12), row.names = FALSE)
cat("\n(for reference) highest-baseline clean-gap guides regardless of low-mode population:\n")
print(head(df[order(-df$baseline), ], 8), row.names = FALSE)
write.csv(df, file.path(OUT_SRC, "violin_guide_search.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT_SRC, "violin_guide_search.csv"), "\n")
