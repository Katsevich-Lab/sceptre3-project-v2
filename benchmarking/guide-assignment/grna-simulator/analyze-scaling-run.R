#!/usr/bin/env Rscript
# Summarize pipeline runs on resampled scaling datasets: cost against size and
# nonzero structure, and the validation arms. PROTOTYPE.
#
#   Rscript analyze-scaling-run.R <out-dir-name> <run_id> [<run_id> ...]
#   e.g. Rscript analyze-scaling-run.R scaling-2026-09-30 rs_overnight rs_rays2
#
# Reads each <outputs>/<run_id>/{<run_id>_config.csv, trace.tsv, assignments_*}
# and each dataset's counts.rds, and writes to
# <LOCAL_BENCHMARKING_DIR>/guide_assignment/scaling-eda/<out-dir-name>/:
#   covariates.csv  per dataset: size, ray (c), nonzeros per cell and per guide
#   tasks.csv       per task: run, status, runtime (h), peak memory (GB), covariates
#   fits.csv        per source x method, over the grid datasets:
#                   log(y) = a + b log(N) + c log(G), y = runtime or memory
#   rays.csv        per source x method x ray with >= 3 grid points: the slope of
#                   log(y) against log(G) along the ray (N/G fixed), i.e. how cost
#                   grows as the whole screen grows
#   arms.csv        validation arms (same cells in A, B, C), where present: per
#                   method, mean assignments per cell, and A vs B agreement
#
# c is N/G relative to the real screen's N/G; a dataset is on ray c when c is
# within 1% of 0.5, 1 or 2. Grid datasets are the weighted ones whose G is not an
# arm's (G0 or 1000). Fits and ray slopes leave out G below FIT_MIN_G (Replogle
# G < 400, where ignoring ambient collisions overstates nonzeros per cell by
# roughly 8%); those tasks stay in tasks.csv.
# Assignment rules: crispat one row per assignment; pertpy "+"-joined guides
# per cell ("negative" = none); cleanser posterior >= 0.8 (as in
# assignment-analysis/assignment-helpers.R), 1-based indices.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: analyze-scaling-run.R <out-dir-name> <run_id> [<run_id> ...]")
runs <- args[-1]
bench    <- .get_config_path("LOCAL_BENCHMARKING_DIR")
in_root  <- file.path(bench, "guide_assignment", "input_data")
run_dirs <- setNames(file.path(bench, "guide_assignment", "outputs", runs), runs)
out_dir  <- file.path(bench, "guide_assignment", "scaling-eda", args[1])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

THRESHOLD <- c(gasperini = 5, `replogle-rd7` = 10)       # as RESAMPLE_SOURCES
ARM_G     <- c(gasperini = 13077, `replogle-rd7` = 2666)
REAL_N_OVER_G <- c(gasperini = 207324 / 13077, `replogle-rd7` = 616184 / 2666)
RAYS      <- c(0.5, 1, 2)
FIT_MIN_G <- c(gasperini = 0, `replogle-rd7` = 400)
CLEANSER_POSTERIOR <- 0.8

parse_id <- function(d) {
  m <- regmatches(d, regexec("^rs_(.+)_(uniform|weighted|none)_g([0-9]+)_c([0-9]+)$", d))[[1]]
  data.frame(dataset = d, source = m[2], placement = m[3],
             G = as.numeric(m[4]), N = as.numeric(m[5]))
}
hours <- function(s) vapply(s, function(x) {
  p <- regmatches(x, gregexpr("[0-9.]+(ms|s|m|h|d)", x))[[1]]
  if (!length(p)) return(NA_real_)
  v <- as.numeric(sub("[a-z]+$", "", p)); u <- sub("^[0-9.]+", "", p)
  sum(v * c(ms = 1e-3, s = 1, m = 60, h = 3600, d = 86400)[u]) / 3600
}, numeric(1), USE.NAMES = FALSE)
gb <- function(s) vapply(s, function(x) {
  m <- regmatches(x, regexec("^([0-9.]+) (KB|MB|GB|TB)$", x))[[1]]
  if (length(m) < 3) return(NA_real_)
  as.numeric(m[2]) * c(KB = 1e-6, MB = 1e-3, GB = 1, TB = 1e3)[[m[3]]]
}, numeric(1), USE.NAMES = FALSE)

cfg <- do.call(rbind, lapply(runs, function(r)
  cbind(run = r, read.csv(file.path(run_dirs[[r]], paste0(r, "_config.csv"))))))
trace <- do.call(rbind, lapply(runs, function(r) {
  t <- read.delim(file.path(run_dirs[[r]], "trace.tsv"), check.names = FALSE)
  data.frame(run = r, process = t$process, tag = t$tag, status = t$status,
             realtime = t$realtime, peak_rss = t$peak_rss, native_id = t$native_id)
}))
ds <- do.call(rbind, lapply(unique(cfg$dataset), parse_id))
ds$c <- ds$N / (ds$G * REAL_N_OVER_G[ds$source])
ds$ray <- vapply(ds$c, function(v) {
  hit <- RAYS[abs(v / RAYS - 1) < 0.01]
  if (length(hit)) hit else NA_real_ }, numeric(1))

# ---- covariates per dataset --------------------------------------------------
cov <- do.call(rbind, lapply(seq_len(nrow(ds)), function(r) {
  d <- ds[r, ]
  m <- readRDS(file.path(in_root, d$dataset, "counts.rds"))
  k <- diff(m@p)                                             # nonzeros per cell
  g <- tabulate(m@i + 1L, nbins = nrow(m))                   # nonzeros per guide
  q <- function(v, p) stats::quantile(v, p, names = FALSE)
  cat(sprintf("covariates: %s\n", d$dataset))
  data.frame(d, nnz = length(m@x),
             moi = sum(m@x >= THRESHOLD[[d$source]]) / ncol(m),
             cell_nnz_mean = mean(k), cell_nnz_q99 = q(k, .99), cell_nnz_max = max(k),
             guide_nnz_mean = mean(g), guide_nnz_cv = stats::sd(g) / mean(g),
             guide_nnz_q99 = q(g, .99), guide_nnz_max = max(g),
             guide_nnz_sumsq = sum(as.numeric(g)^2))
}))
write.csv(cov, file.path(out_dir, "covariates.csv"), row.names = FALSE)

# ---- tasks -----------------------------------------------------------------
tr <- data.frame(run = trace$run, method = tolower(sub("_ASSIGN$", "", trace$process)),
                 dataset = trace$tag, status = trace$status, hours = hours(trace$realtime),
                 peak_gb = gb(trace$peak_rss), native_id = trace$native_id)
tasks <- merge(cfg[, c("run", "method", "dataset", "cpus", "memory", "time")], tr,
               by = c("run", "method", "dataset"), all.x = TRUE)
tasks$status[is.na(tasks$status)] <- "NOT_FINISHED"
tasks <- merge(tasks, cov, by = "dataset")
write.csv(tasks, file.path(out_dir, "tasks.csv"), row.names = FALSE)

# ---- fits over the grid --------------------------------------------------------
grid <- subset(tasks, placement == "weighted" & G != 1000 & G != ARM_G[source] &
                 G >= FIT_MIN_G[source] & status == "COMPLETED")
fits <- do.call(rbind, lapply(split(grid, list(grid$source, grid$method), drop = TRUE), function(s) {
  do.call(rbind, lapply(c("hours", "peak_gb"), function(y) {
    ok <- is.finite(s[[y]]) & s[[y]] > 0
    if (sum(ok) < 4) return(NULL)
    f <- stats::lm(log(s[[y]][ok]) ~ log(s$N[ok]) + log(s$G[ok]))
    co <- stats::coef(f); se <- sqrt(diag(stats::vcov(f)))
    data.frame(source = s$source[1], method = s$method[1], response = y, n = sum(ok),
               b_N = co[2], se_N = se[2], b_G = co[3], se_G = se[3],
               r2 = summary(f)$r.squared, max_resid = max(abs(stats::resid(f))))
  }))
}))
write.csv(fits, file.path(out_dir, "fits.csv"), row.names = FALSE)

# ---- slopes along each ray ---------------------------------------------------
on_ray <- subset(grid, !is.na(ray))
rays <- do.call(rbind, lapply(split(on_ray, list(on_ray$source, on_ray$method, on_ray$ray), drop = TRUE),
                              function(s) {
  do.call(rbind, lapply(c("hours", "peak_gb"), function(y) {
    ok <- is.finite(s[[y]]) & s[[y]] > 0
    if (sum(ok) < 3) return(NULL)
    f <- stats::lm(log(s[[y]][ok]) ~ log(s$G[ok]))
    data.frame(source = s$source[1], method = s$method[1], ray = s$ray[1], response = y,
               n = sum(ok), G_min = min(s$G[ok]), G_max = max(s$G[ok]),
               slope = stats::coef(f)[2], se = sqrt(diag(stats::vcov(f)))[2],
               r2 = summary(f)$r.squared)
  }))
}))
write.csv(rays, file.path(out_dir, "rays.csv"), row.names = FALSE)

# ---- validation arms -------------------------------------------------------
per_cell <- function(method, dataset, n_cells) {
  r <- tasks$run[tasks$method == method & tasks$dataset == dataset][1]
  if (is.na(r)) return(NULL)
  fp <- file.path(run_dirs[[r]], sprintf("assignments_%s_%s.csv", method, dataset))
  if (!file.exists(fp)) return(NULL)
  a <- read.csv(fp)
  idx <- switch(method,
    crispat  = as.integer(sub("^CELL_", "", a$cell_id)),
    pertpy   = { keep <- a$grna_id != "negative"
                 rep(as.integer(sub("^CELL_", "", a$cell_id[keep])),
                     lengths(strsplit(a$grna_id[keep], "+", fixed = TRUE))) },
    cleanser = a$cell_id[a$posterior >= CLEANSER_POSTERIOR])
  tabulate(idx, nbins = n_cells)
}
arm_sources <- names(ARM_G)[vapply(names(ARM_G), function(src)
  any(ds$source == src & ds$placement == "none"), logical(1))]
arms <- do.call(rbind, lapply(arm_sources, function(src) {
  n <- unique(ds$N[ds$source == src & ds$placement == "none"])
  ids <- c(A = sprintf("rs_%s_none_g%d_c%d", src, ARM_G[[src]], n),
           B = sprintf("rs_%s_weighted_g%d_c%d", src, ARM_G[[src]], n),
           C = sprintf("rs_%s_weighted_g1000_c%d", src, n))
  do.call(rbind, lapply(c("crispat", "cleanser", "pertpy"), function(m) {
    pc <- lapply(ids, function(d) per_cell(m, d, n))
    rt <- vapply(ids, function(d) {
      h <- tasks$hours[tasks$method == m & tasks$dataset == d & tasks$status == "COMPLETED"]
      if (length(h)) h else NA_real_ }, numeric(1))
    mean_pc <- vapply(pc, function(v) if (is.null(v)) NA_real_ else mean(v), numeric(1))
    ab <- !is.null(pc$A) && !is.null(pc$B)
    data.frame(source = src, method = m, n_cells = n,
               per_cell_A = mean_pc[["A"]], per_cell_B = mean_pc[["B"]], per_cell_C = mean_pc[["C"]],
               hours_A = rt[["A"]], hours_B = rt[["B"]], hours_C = rt[["C"]],
               cor_AB = if (ab) stats::cor(pc$A, pc$B) else NA_real_,
               same_count_AB = if (ab) mean(pc$A == pc$B) else NA_real_)
  }))
}))
if (!is.null(arms)) write.csv(arms, file.path(out_dir, "arms.csv"), row.names = FALSE)

options(width = 200)
cat("\n=== fits (grid) ===\n");  print(fits, row.names = FALSE, digits = 3)
cat("\n=== ray slopes ===\n");   print(rays, row.names = FALSE, digits = 3)
if (!is.null(arms)) { cat("\n=== arms ===\n"); print(arms, row.names = FALSE, digits = 3) }
cat("\nwrote", out_dir, "\n")
