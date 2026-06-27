#!/usr/bin/env Rscript
# =============================================================================
# Systematic characterization of gRNA UMI distributions across ALL real datasets,
# to derive a data-driven grid of simulation regimes (one per dataset).
#
# Per guide we fit the bimodal structure on the raw counts (smoothed-valley):
#   left_mode  : the lower (ambient/background) mode, raw count
#   right_mode : the upper (signal/perturbed) mode, raw count
#   separation : log1p(right_mode) - log1p(left_mode)  -- how separated the modes are
#   signal_frac: fraction of cells above the per-guide valley (~ how often perturbed)
#   ambient below-valley: mean nonzero, % >= 3, max  -- the "ambiguous middle"/contamination
# Per dataset we also estimate MOI (median guides/cell above their valley) + library.
# No ground truth needed: these are approximate, mode-based, method-light estimates.
# Output: results/sim_framework/real_characterization.csv + .png
# =============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
SURV <- path.expand("~/data/external/perturbseq-survey")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
REPRO<- file.path(HERE, "external", "repro_work")
OUT  <- SIMFW()

# ---- assemble the dataset registry (loaders -> guides x cells matrix) --------
loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) { f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) }) }
add_rds <- function(nm, f) if (file.exists(f)) loaders[[nm]] <<- local({ ff <- f; function() readRDS(ff) })
add_rds("gasperini", file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
add_rds("replogle",  file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { suppressPackageStartupMessages(library(fishash)); data(crispat_schraivogel); crispat_schraivogel }
for (s in BARN_SAMPLES) { mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (file.exists(mtx)) loaders[[paste0("barnyard_", sub("mix","",s))]] <- local({ mm <- mtx
    function() as(readMM(mm), "CsparseMatrix") }) }

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
    cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < 30) next
    v <- smoothed_valley_threshold(cv)
    if (isTRUE(v$ok)) {
      valley[g] <- v$t
      bl <- cv[cv > 0 & cv < v$t]; below <- c(below, bl)
      sig <- cv[cv >= v$t]; ms <- mean(sig); vs <- if (length(sig) > 1) var(sig) else ms
      rows[[length(rows) + 1]] <- data.frame(
        left_mode = v$mode1, right_mode = v$mode2,
        separation = log1p(v$mode2) - log1p(v$mode1),
        signal_frac = mean(cv >= v$t),
        mu = ms,                                              # signal mean (above valley)
        theta = if (vs > ms && ms > 0) ms^2 / (vs - ms) else 1e3)  # NB dispersion (method of moments)
    }
  }
  pg <- if (length(rows)) do.call(rbind, rows) else
    data.frame(left_mode=NA, right_mode=NA, separation=NA, signal_frac=NA, mu=NA, theta=NA)[0,]
  # MOI over ALL guides via the ambient test (the valley fit subsamples guides,
  # which would under-count MOI for many-guide datasets).
  A <- suppressMessages(ambient_test_assign(gm, q = 0.05, model = "hypergeometric", n_iter = 1)$assignment_matrix)
  moi <- Matrix::colSums(A)
  bv <- if (length(below)) c(mean = mean(below), p3 = mean(below >= 3), mx = max(below)) else c(mean=NA,p3=NA,mx=NA)
  out <- data.frame(n_guides = G, n_cells = C, n_bimodal = nrow(pg),
             left_mode = median(pg$left_mode, na.rm = TRUE),
             right_mode = median(pg$right_mode, na.rm = TRUE),
             separation = median(pg$separation, na.rm = TRUE),
             signal_frac = median(pg$signal_frac, na.rm = TRUE),
             moi = median(moi), lib_med = median(Matrix::colSums(gm)),
             amb_mean = round(bv["mean"], 2), amb_pct3 = round(100 * bv["p3"], 1), amb_max = bv["mx"])
  keep_g <- is.finite(pg$mu) & pg$mu > 0 & is.finite(pg$signal_frac) & is.finite(pg$theta)
  attr(out, "perguide") <- pg[keep_g, c("mu", "theta", "signal_frac", "right_mode")]  # per-guide fit
  out
}

rows <- list(); perguide <- list()
for (nm in names(loaders)) {
  cr <- tryCatch(characterize(loaders[[nm]]()),
                 error = function(e) { cat(sprintf("%-28s ERROR %s\n", nm, conditionMessage(e))); NULL })
  if (!is.null(cr)) {
    rows[[nm]] <- cbind(dataset = nm, cr); perguide[[nm]] <- attr(cr, "perguide")
    cat(sprintf("%-28s G=%-6d C=%-6d  Lmode=%-3g Rmode=%-5g sep=%-4.2f sigfrac=%-5.3f MOI=%-4.1f amb%%>=3=%-4.1f ambmax=%g  (n_perguide=%d)\n",
        nm, cr$n_guides, cr$n_cells, cr$left_mode, cr$right_mode, cr$separation, cr$signal_frac, cr$moi, cr$amb_pct3, cr$amb_max, nrow(perguide[[nm]])))
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(OUT, "real_characterization.csv"), row.names = FALSE)
saveRDS(perguide, file.path(OUT, "real_perguide.rds"))            # per-dataset per-guide signal modes
cat("\nwrote real_characterization.csv + real_perguide.rds (", nrow(tab), "datasets )\n")

# ---- plot: where each dataset sits in (separation x signal level x MOI) ------
pd <- tab %>% mutate(dataset = reorder(dataset, separation))
p <- ggplot(pd, aes(separation, right_mode, size = moi, colour = amb_pct3, label = dataset)) +
  geom_point(alpha = 0.85) + geom_text(size = 2.6, vjust = -0.9, check_overlap = TRUE) +
  scale_y_log10() + scale_colour_viridis_c(option = "C", name = "ambient %>=3") +
  scale_size_continuous(name = "MOI") +
  labs(title = "Real gRNA parameter space: 17 datasets span mode separation x signal level x MOI x ambient",
       x = "mode separation  log1p(right) - log1p(left)  (low = hard)", y = "right (signal) mode, raw UMI (log)") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "real_characterization.png"), p, width = 11, height = 6.5, dpi = 140)
cat("wrote real_characterization.png\n")
