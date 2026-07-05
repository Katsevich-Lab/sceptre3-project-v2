#!/usr/bin/env Rscript
# ============================================================================
# Comprehensive look at ambient DEPTH (per cell) and COMPOSITION (per guide)
# across EVERY real dataset in this effort:
#   - 10 perturbseq-survey datasets (2024-26)
#   - Gasperini + Replogle references
#   - Schraivogel (fishash package)
#   - 4 barnyard samples (Liu 2025: CROP-seq / direct-capture x 0/72hr)
#
# Ambient depth is defined THRESHOLD-FREE: run our ambient test, then
#   ambient_depth_c = library_size_c - (UMIs in assigned guides)_c
# and ambient composition gamma_g = ambient UMIs of guide g / all ambient UMIs.
# Reports, per dataset: how much library size overstates ambient depth, how
# uneven the ambient pool is, and how contaminated the naive composition is.
# ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(sparseMatrixStats)
})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
GA   <- normalizePath(file.path(HERE, ".."))
SURV <- path.expand("~/data/external/perturbseq-survey")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
REPRO<- file.path(HERE, "external", "repro_work")
OUT  <- file.path(HERE, "results", "ambient_intuition"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

# ---- assemble the dataset registry (loaders return a guides x cells matrix) --
loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) {
  f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) })
}
add_rds <- function(nm, f) if (file.exists(f)) loaders[[nm]] <<- local({ ff <- f; function() readRDS(ff) })
add_rds("gasperini_ref",  file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
add_rds("replogle_ref",   file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { library(fishash); data(crispat_schraivogel); crispat_schraivogel }
for (s in c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")) {
  mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (file.exists(mtx)) loaders[[paste0("barnyard_", s)]] <- local({ mm <- mtx
    function() as(readMM(mm), "CsparseMatrix") })
}

platform <- function(nm) {
  if (grepl("dctap|DirectCapture|tcell|cd8|cd4|hipsci|ipsc", nm, ignore.case = TRUE) &
      !grepl("Cropseq", nm)) "direct-capture-ish" else "cropseq/perturbseq"
}

MAXC <- 25000
metrics <- function(gm) {
  gm <- as(gm, "CsparseMatrix"); storage.mode(gm@x) <- "double"
  if (ncol(gm) > MAXC) { set.seed(1); gm <- gm[, sort(sample(ncol(gm), MAXC))] }
  gm <- gm[, colSums(gm) > 0, drop = FALSE]                 # drop empty cells
  G <- nrow(gm); C <- ncol(gm)
  lib  <- colSums(gm); nnz <- colSums(gm > 0)
  A <- suppressMessages(ambient_test_assign(gm, q = 0.05, model = "hypergeometric",
                                            n_iter = 1)$assignment_matrix)
  A <- as(A, "CsparseMatrix") * 1                          # logical -> numeric sparse
  sig_cell  <- colSums(gm * A)                              # signal UMIs per cell
  amb_depth <- pmax(lib - sig_cell, 0)
  top <- colMaxs(gm); top_frac <- top / pmax(lib, 1)
  gsub <- gm; gsub@x[gsub@x > 2] <- 0                       # CLEANSER sub-threshold (<=2)
  kappa_clean <- colSums(gsub)
  q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
  gini <- function(x) { x <- sort(x); n <- length(x); 2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n }
  # composition: restrict to WELL-OBSERVED guides (>=10 ambient UMIs) so the
  # ambient rate gamma_g is estimable and ratios don't blow up on empty guides.
  tot_g <- rowSums(gm); sig_g <- rowSums(gm * A); amb_g <- tot_g - sig_g  # ambient UMIs/guide
  well  <- amb_g >= 10
  rate  <- amb_g[well] / sum(amb_g[well])                  # ambient composition (well guides)
  s2a   <- tot_g[well] / amb_g[well]                       # total / ambient per guide (>=1)
  # depth ratio per cell (how many-fold library overstates ambient depth)
  dr <- lib / pmax(amb_depth, 0.5)
  assigned_cells <- sig_cell > 0
  data.frame(
    n_guides = G, n_cells = C, lib_med = median(lib), nnz_med = median(nnz),
    moi = median(colSums(A)), assigned_frac = mean(assigned_cells),
    signal_frac_med = median((sig_cell / lib)[assigned_cells]),
    depth_ratio_med = median(dr), depth_ratio_p90 = q(dr, .90),
    depth_cor  = suppressWarnings(cor(log1p(lib), log1p(amb_depth))),
    clean_cor  = suppressWarnings(cor(log1p(kappa_clean), log1p(amb_depth))),
    n_well = sum(well), comp_gini = gini(amb_g[well]),
    comp_fold_90_10 = q(rate, .90) / q(rate, .10),
    guide_contam_spread = q(s2a, .90) / q(s2a, .10))
}

rows <- list()
for (nm in names(loaders)) {
  r <- tryCatch({
    gm <- loaders[[nm]]()
    m <- metrics(gm); m$dataset <- nm; m$platform <- platform(nm)
    cat(sprintf("%-42s G=%-5d C=%-6d MOI=%-4.1f lib=%-5.0f  libsize/depth=%5.1fx (p90 %5.0fx)  sigfrac=%.2f  Gini=%.2f  comp10-90fold=%4.1fx  contam=%4.1fx\n",
        nm, m$n_guides, m$n_cells, m$moi, m$lib_med, m$depth_ratio_med, m$depth_ratio_p90,
        m$signal_frac_med, m$comp_gini, m$comp_fold_90_10, m$guide_contam_spread))
    m
  }, error = function(e) { cat(sprintf("%-42s ERROR: %s\n", nm, conditionMessage(e))); NULL })
  if (!is.null(r)) rows[[nm]] <- r
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
tab <- tab[, c("dataset","platform","n_guides","n_cells","moi","lib_med","nnz_med",
               "assigned_frac","signal_frac_med","depth_ratio_med","depth_ratio_p90",
               "depth_cor","clean_cor","n_well","comp_gini","comp_fold_90_10","guide_contam_spread")]
write.csv(tab, file.path(OUT, "survey_depth_composition.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "survey_depth_composition.csv"), "(", nrow(tab), "datasets )\n")
