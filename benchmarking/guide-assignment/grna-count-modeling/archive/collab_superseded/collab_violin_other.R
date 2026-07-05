## Single-guide expression violin for the OTHER clean-gap cell types (EndoC, CD4) -- the companion of
## scripts/collab_eno1_violin.R + find_violin_guide.R (Replogle). We (1) SEARCH each dataset's clean-gap
## targeting guides for one that makes a good violin -- well-expressed target, a populated low mode across
## a few counts, and clear knockdown at the integration mode -- and (2) draw the violin only where such a
## guide exists (survey low modes are shallow, so a dataset may not qualify; we say so).
##
## Same survey data + normalization as collab_dose_response_aggregate_other.R / dose_survey():
##   grna_matrix_aligned.rds (guides x cells), response_matrix.rds (genes x cells), grna_target_data_frame.csv;
##   library size = colSums(response); normalized expr = response[target,]/lib*1e4.
##   EndoC baseline over NT-positive cells with count 0; CD4 (no NT) complement baseline (cells with count 0).
##
## Run from nonparametric-thresholds/.
suppressMessages({library(Matrix); library(ggplot2); library(patchwork)})
HERE    <- "/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/nonparametric-thresholds"
OUT_SRC <- file.path(HERE, "results/global_ambient_poisson")
OUT     <- file.path(HERE, "results/collaborator_writeup"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SV      <- path.expand("~/data/external/perturbseq-survey")

## returns list(td/g/r/lib + ntmask + candidate search table) for a survey dataset
load_survey <- function(dir, has_nt, pgfile) {
  g  <- as(readRDS(file.path(dir, "grna_matrix_aligned.rds")), "CsparseMatrix")
  r  <- as(readRDS(file.path(dir, "response_matrix.rds")),      "CsparseMatrix")
  td <- read.csv(file.path(dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
  lib  <- as.numeric(Matrix::colSums(r)); rn <- rownames(r); tmap <- setNames(td$grna_target, td$grna_id)
  nt <- td$grna_id[td$grna_target == "non-targeting"]; ntr <- match(nt, rownames(g)); ntr <- ntr[!is.na(ntr)]
  ntmask <- if (has_nt && length(ntr)) Matrix::colSums(g[ntr, , drop = FALSE] >= 30) >= 1 else rep(TRUE, ncol(g))
  pg <- read.csv(file.path(OUT_SRC, pgfile), stringsAsFactors = FALSE); pg <- pg[pg$gap_lo >= 2, ]
  pg$target <- tmap[pg$guide]; pg$symbol <- td$grna_target_symbol[match(pg$guide, td$grna_id)]
  pg <- pg[!is.na(pg$target) & pg$target != "non-targeting" & pg$target %in% rn & pg$guide %in% rownames(g), ]
  Gsub <- as.matrix(g[pg$guide, , drop = FALSE]); utg <- unique(pg$target); Est <- t(r[utg, , drop = FALSE]); colnames(Est) <- utg
  rows <- list()
  for (i in seq_len(nrow(pg))) {
    v <- as.numeric(Gsub[i, ]); cp <- as.numeric(Est[, pg$target[i]]) / lib * 1e4
    base <- mean(cp[ntmask & v == 0]); lo <- pg$gap_lo[i]; hi <- pg$gap_hi[i]; st <- which(v >= hi)
    kd <- if (length(st)) mean(cp[st]) / base else NA
    rows[[i]] <- data.frame(guide = pg$guide[i], symbol = pg$symbol[i], target = pg$target[i],
      gap_lo = lo, gap_hi = hi, baseline = round(base, 2), dropout0 = round(mean(cp[v == 0] == 0), 2),
      knockdown = round(kd, 2), n_sig = length(st),
      n1 = sum(v == 1 & v <= lo), n2 = sum(v == 2 & v <= lo), n3 = sum(v == 3 & v <= lo),
      n4p = sum(v >= 4 & v <= lo), n_low = sum(v >= 1 & v <= lo))
  }
  list(g = g, Gsub = Gsub, Est = Est, lib = lib, ntmask = ntmask, pg = pg, cand = do.call(rbind, rows))
}

## "good violin" criteria, RELAXED vs Replogle because survey low modes are shallow: well-expressed
## (baseline>=1, low dropout), a real hit (knockdown<=0.6), enough integration cells (n_sig>=5), and a
## low mode with cells at >=2 distinct counts and >=6 cells total so per-count bins aren't degenerate.
pick_guide <- function(cand) {
  ok <- cand$baseline >= 1 & cand$knockdown <= 0.6 & cand$n_sig >= 5 &
        cand$n_low >= 6 & cand$n2 >= 2 & (cand$n3 + cand$n4p) >= 1
  c2 <- cand[ok, ]
  if (!nrow(c2)) return(NULL)
  c2 <- c2[order(-c2$baseline * (1 - c2$knockdown) * pmin(c2$n_low, 60)), ]
  c2
}

## draw the violin for a chosen guide (per-cell target expr / baseline, by gRNA-count bin). Bins adapt
## to the guide's gap: low-mode singles up to gap_lo, then the integration mode (>=gap_hi).
make_violin <- function(D, row, title) {
  v <- as.numeric(D$Gsub[match(row$guide, D$pg$guide), ])
  cp <- as.numeric(D$Est[, row$target]) / D$lib * 1e4; base <- mean(cp[D$ntmask & v == 0])
  lo <- row$gap_lo; hi <- row$gap_hi
  ## Three clean, well-populated bins (survey low modes are shallow -- counts 2,3 alone are near-empty,
  ## so splitting them per-count would give degenerate violins): the no-guide mode (0), the whole
  ## below-gap LOW mode (1..gap_lo, the candidate weak-integration cells), and the integration mode
  ## (>=gap_hi, full knockdown). This directly shows the low mode sitting between baseline and knockdown.
  lo_lab <- if (lo == 1) "1" else paste0("1-", lo); hi_lab <- paste0(">=", hi)
  bin <- rep(NA_character_, length(v))
  bin[v == 0] <- "0"
  bin[v >= 1 & v <= lo] <- lo_lab
  bin[v >= hi] <- hi_lab
  df <- data.frame(count = v, rel = cp / base, bin = bin); df <- df[!is.na(df$bin), ]
  lev <- c("0", lo_lab, hi_lab); lev <- lev[lev %in% df$bin]; df$bin <- factor(df$bin, levels = lev)
  summ <- do.call(rbind, lapply(lev, function(b) { s <- df$rel[df$bin == b]
    data.frame(bin = b, n = length(s), mean = round(mean(s), 3), median = round(median(s), 3)) }))
  ytop <- as.numeric(quantile(df$rel, 0.99))
  cat(sprintf("[violin %s] guide=%s target=%s baseline=%.2f gap=%d-%d\n", title, row$guide, row$symbol, base, lo, hi))
  print(summ, row.names = FALSE)
  p <- ggplot(df, aes(bin, rel)) +
    geom_violin(fill = "#0072B2", colour = "grey30", alpha = 0.55, scale = "width", linewidth = 0.3) +
    stat_summary(fun = mean, geom = "point", colour = "#D55E00", size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_text(data = summ, aes(bin, y = ytop, label = paste0("n=", n)), vjust = 1, size = 3.6,
              colour = "grey35", inherit.aes = FALSE) +
    coord_cartesian(ylim = c(0, ytop)) +
    labs(title = title, x = sprintf("%s gRNA UMI count in the cell", row$symbol),
         y = sprintf("%s target expression / no-guide baseline", row$symbol)) +
    theme_bw(base_size = 16)
  list(p = p, summ = summ)
}

cat("################## searching EndoC ##################\n")
Den <- load_survey(file.path(SV, "endoc_t2d_perturbseq_GSE273677/sceptre"), TRUE, "perguide_endoc_t2d.csv")
en_c <- pick_guide(Den$cand)
cat("EndoC candidate table (top by composite; NULL if none qualify):\n")
if (is.null(en_c)) cat("  -- NO EndoC guide qualifies (dense-soup: shallow low modes / few real knockdowns)\n") else
  print(head(en_c[, c("guide","symbol","gap_lo","gap_hi","baseline","dropout0","knockdown","n_sig","n1","n2","n3","n4p","n_low")], 8), row.names = FALSE)

cat("\n################## searching CD4 (genome-scale) ##################\n")
Dtc <- load_survey(file.path(SV, "tcell_cd4_perturbseq_GSE314342/sceptre"), FALSE, "perguide_tcell_cd4.csv")
tc_c <- pick_guide(Dtc$cand)
cat("CD4 candidate table (top by composite; NULL if none qualify):\n")
if (is.null(tc_c)) cat("  -- NO CD4 guide qualifies\n") else
  print(head(tc_c[, c("guide","symbol","gap_lo","gap_hi","baseline","dropout0","knockdown","n_sig","n1","n2","n3","n4p","n_low")], 8), row.names = FALSE)

## write the full search tables for provenance
write.csv(Den$cand, file.path(OUT, "violin_guide_search_endoc.csv"), row.names = FALSE)
write.csv(Dtc$cand, file.path(OUT, "violin_guide_search_cd4.csv"),   row.names = FALSE)

## Build violins only for datasets with a qualifying guide. If both qualify -> two-panel; if one ->
## single-panel; if neither -> no figure (report it).
panels <- list(); chosen <- list()
if (!is.null(en_c)) { r <- en_c[1, ]; vp <- make_violin(Den, r, sprintf("EndoC — %s (gap %d–%d)", r$symbol, r$gap_lo, r$gap_hi))
  panels[["endoc"]] <- vp$p; chosen[["endoc"]] <- r }
if (!is.null(tc_c)) { r <- tc_c[1, ]; vp <- make_violin(Dtc, r, sprintf("CD4 — %s (gap %d–%d)", r$symbol, r$gap_lo, r$gap_hi))
  panels[["cd4"]] <- vp$p; chosen[["cd4"]] <- r }

if (length(panels) == 0) {
  cat("\nNo single-guide violin written: no qualifying guide in either dataset.\n")
} else if (length(panels) == 1) {
  ds <- names(panels)[1]; g <- chosen[[ds]]$guide
  ggsave(file.path(OUT, sprintf("expression_violin_%s.png", ds)), panels[[1]], width = 9, height = 5.2, dpi = 130)
  cat(sprintf("\nwrote %s  (single qualifying guide: %s = %s)\n",
              file.path(OUT, sprintf("expression_violin_%s.png", ds)), g, chosen[[ds]]$symbol))
} else {
  fig <- panels[["endoc"]] + panels[["cd4"]] +
    plot_annotation(title = "Single-guide knockdown across the low mode (other cell types)",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  ggsave(file.path(OUT, "expression_violin_other.png"), fig, width = 15, height = 5.4, dpi = 130)
  cat("\nwrote", file.path(OUT, "expression_violin_other.png"), "(EndoC + CD4)\n")
}
cat("done\n")
