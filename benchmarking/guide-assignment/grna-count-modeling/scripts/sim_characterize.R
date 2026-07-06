#!/usr/bin/env Rscript
# =============================================================================
# Systematic characterization of gRNA UMI distributions across ALL real datasets,
# to derive a data-driven grid of simulation regimes (one per dataset).
#
# Per guide we fit the bimodal structure on the raw counts (smoothed-valley) and
# save the parameters needed to RE-SIMULATE each guide (sim_regimes.R draws G
# real per-guide fits and builds Model B from them):
#   mu          : mean of above-valley cells   (NB mean for the sim signal)
#   theta       : NB dispersion (method-of-moments)
#   signal_frac : fraction of cells above valley (the per-guide perturbation rate)
#   left_mode/right_mode/separation : descriptive stats for plotting/ordering
#
# KNOWN BIASES (documented; partially compensated):
#  (a) mu/theta are fit on the left-TRUNCATED above-valley sample, so mu is
#      slightly upward-biased and theta upward-biased (variance underestimated).
#      Bias is small for well-separated guides; we set theta to a large cap (1e3)
#      when the sample is too small (< MIN_SIG_CELLS) to be informative.
#  (b) signal_frac counts ambient leakage above the valley as signal; sim_regimes.R
#      DEFLATES signal_frac by the dataset's ambient-above-valley mass before use.
#
# Dataset-level: MOI (ambient test over ALL guides), library, and the ambient
# "ambiguous middle" stats (mean nonzero, % >= 3, p99, max) on below-valley counts.
# Output: results/sim_framework/real_characterization.csv + real_perguide.rds + .png
# =============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
SURV <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
DATA <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
REPRO<- file.path(HERE, "external", "repro_work")
OUT  <- SIMFW()

MIN_NZ          <- 30   # min nonzero cells to attempt a fit (shared with sim_regimes.R)
MIN_SIG_CELLS   <- 10   # min above-valley cells for a meaningful theta MoM fit

# ---- assemble the dataset registry (loaders -> guides x cells matrix) --------
# Barnyard uses the authors' QC (purity 0.9), matching the cohort Model A
# consumes -- otherwise MOI/ambient stats would be apples-to-oranges.
loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) { f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) }) }
add_rds <- function(nm, f) if (file.exists(f)) loaders[[nm]] <<- local({ ff <- f; function() readRDS(ff) })
add_rds("gasperini", file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
add_rds("replogle",  file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { suppressPackageStartupMessages(library(fishash)); data(crispat_schraivogel); crispat_schraivogel }
for (s in BARN_SAMPLES) {
  loaders[[paste0("barnyard_", sub("mix", "", s))]] <- local({ ss <- s
    function() load_barnyard(ss, purity = 0.9)$counts })   # QC-matched to Model A
}

# ---- per-dataset characterization -------------------------------------------
characterize <- function(gm, max_cells = 40000, max_guides = 500, seed = 1) {
  gm <- as(gm, "CsparseMatrix"); storage.mode(gm@x) <- "double"
  set.seed(seed)
  if (ncol(gm) > max_cells) gm <- gm[, sort(sample(ncol(gm), max_cells))]
  gm <- gm[, Matrix::colSums(gm) > 0, drop = FALSE]
  G <- nrow(gm); C <- ncol(gm)
  gmr <- as(gm, "RsparseMatrix")
  gsel <- if (G > max_guides) sort(sample(G, max_guides)) else seq_len(G)
  valley <- rep(Inf, G)                        # per-guide signal threshold (for MOI)
  rows <- list(); below <- c()
  for (g in gsel) {
    cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < MIN_NZ) next
    v <- smoothed_valley_threshold(cv)
    if (isTRUE(v$ok)) {
      valley[g] <- v$t
      bl <- cv[cv > 0 & cv < v$t]; below <- c(below, bl)
      sig <- cv[cv >= v$t]; ns <- length(sig); ms <- mean(sig)
      vs <- if (ns > 1) var(sig) else ms
      # theta = NB dispersion (method of moments).  Cap when sample too small
      # to estimate; the cap is later subjected to a pmin in sim_regimes.R.
      theta_g <- if (ns >= MIN_SIG_CELLS && vs > ms && ms > 0) ms^2 / (vs - ms) else 1e3
      rows[[length(rows) + 1]] <- data.frame(
        left_mode = v$mode1, right_mode = v$mode2,
        separation = log1p(v$mode2) - log1p(v$mode1),
        signal_frac = mean(cv >= v$t),
        mu = ms, theta = theta_g)
    }
  }
  pg <- if (length(rows)) do.call(rbind, rows) else
    data.frame(left_mode=NA, right_mode=NA, separation=NA, signal_frac=NA, mu=NA, theta=NA)[0,]
  keep_g <- is.finite(pg$mu) & pg$mu > 0 & is.finite(pg$signal_frac) & is.finite(pg$theta)
  # MOI over ALL guides via the ambient test (the valley fit subsamples guides,
  # which would under-count MOI for many-guide datasets).
  A <- suppressMessages(ambient_test_assign(gm, q = 0.05, model = "hypergeometric",
                                            n_iter = 1)$assignment_matrix)
  moi <- Matrix::colSums(A)
  bv <- if (length(below)) c(mean = mean(below), p3 = mean(below >= 3),
                             p99 = as.numeric(quantile(below, 0.99)), mx = max(below))
        else c(mean = NA, p3 = NA, p99 = NA, mx = NA)
  out <- data.frame(n_guides = G, n_cells = C, n_bimodal = sum(keep_g),
             left_mode = median(pg$left_mode, na.rm = TRUE),
             right_mode = median(pg$right_mode, na.rm = TRUE),
             separation = median(pg$separation, na.rm = TRUE),
             signal_frac = median(pg$signal_frac, na.rm = TRUE),
             moi = median(moi), lib_med = median(Matrix::colSums(gm)),
             amb_mean = round(bv["mean"], 2), amb_pct3 = round(100 * bv["p3"], 1),
             amb_p99 = round(bv["p99"], 1), amb_max = bv["mx"])
  attr(out, "perguide") <- pg[keep_g, c("mu", "theta", "signal_frac", "right_mode")]
  out
}

rows <- list(); perguide <- list()
for (nm in names(loaders)) {
  cr <- tryCatch(characterize(loaders[[nm]]()),
                 error = function(e) { cat(sprintf("%-28s ERROR %s\n", nm, conditionMessage(e))); NULL })
  if (!is.null(cr)) {
    rows[[nm]] <- cbind(dataset = nm, cr); perguide[[nm]] <- attr(cr, "perguide")
    cat(sprintf("%-28s G=%-6d C=%-6d  Lmode=%-3g Rmode=%-5g sep=%-4.2f sigfrac=%-5.3f MOI=%-4.1f amb%%>=3=%-4.1f ambp99=%g ambmax=%g  (n_perguide=%d)\n",
        nm, cr$n_guides, cr$n_cells, cr$left_mode, cr$right_mode, cr$separation, cr$signal_frac,
        cr$moi, cr$amb_pct3, cr$amb_p99, cr$amb_max, nrow(perguide[[nm]])))
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(OUT, "real_characterization.csv"), row.names = FALSE)
saveRDS(perguide, file.path(OUT, "real_perguide.rds"))
cat(sprintf("\nwrote real_characterization.csv + real_perguide.rds (%d datasets)\n", nrow(tab)))

# ---- plot: where each dataset sits in (separation x signal level x MOI) ------
pd <- tab %>% mutate(dataset = reorder(dataset, separation))
p <- ggplot(pd, aes(separation, right_mode, size = moi, colour = amb_pct3, label = dataset)) +
  geom_point(alpha = 0.85) + geom_text(size = 2.6, vjust = -0.9, check_overlap = TRUE) +
  scale_y_log10() + scale_colour_viridis_c(option = "C", name = "ambient %>=3") +
  scale_size_continuous(name = "MOI") +
  labs(title = sprintf("Real gRNA parameter space: %d datasets span mode separation x signal level x MOI x ambient", nrow(tab)),
       x = "mode separation  log1p(right) - log1p(left)  (low = hard)",
       y = "right (signal) mode, raw UMI (log)") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "real_characterization.png"), p, width = 11, height = 6.5, dpi = 140)
cat("wrote real_characterization.png\n")
