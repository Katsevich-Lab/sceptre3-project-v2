#!/usr/bin/env Rscript
# Time fishash and fishash+ on the full real matrices, on this machine.
#
#   Rscript time-fishash-variants.R                       # both datasets
#   Rscript time-fishash-variants.R gasperini             # one dataset
#   Rscript time-fishash-variants.R --skip-odm            # skip the .odm reads
#
# Needs the R_461 renv library, with fishash, ondisc and fishashplus installed
# into it:
#   LIB=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#   R_LIBS=$LIB R CMD INSTALL --preclean --no-docs --library=$LIB \
#     ~/katsevich-lab/code/fishash_plus_cpp
#
# --preclean matters: pkgload::load_all() compiles at -O0, and an install that
# finds those objects newer than the sources just relinks them, so the timing
# would be of an unoptimized build. This script refuses to run under load_all for
# the same reason, but it cannot detect stale objects in an installed package.
#
# Reports read time and assignment time separately, because they answer different
# questions. fishash+ is run on the matrix from each format in turn, so that
#   mtx vs odm      isolates the input format, the method being identical
#   fishash vs +    isolates the method, the input being identical
# The two formats are checked for equality before anything is timed: a speed
# comparison across formats is meaningless if they do not hold the same matrix.

suppressPackageStartupMessages({
  library(Matrix); library(fishash)
})
suppressPackageStartupMessages(library(fishashplus))
if (!is.null(pkgload::dev_meta("fishashplus")))
  stop("fishashplus is loaded via load_all (compiled -O0); install it instead")
source("~/.Rprofile")                                    # .get_config_path()

args     <- commandArgs(trailingOnly = TRUE)
skip_odm <- "--skip-odm" %in% args
picked   <- setdiff(args, "--skip-odm")

DATASETS <- c(gasperini = "gasperini", replogle = "replogle-rd7")
if (length(picked)) DATASETS <- DATASETS[picked]
if (!length(DATASETS)) stop("no dataset matched; use gasperini and/or replogle")

IN <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")

# Wall time and the peak R heap the call added, in GiB. gc() max-used is R's own
# accounting, not RSS -- it undercounts anything allocated outside R's heap.
timed <- function(expr) {
  invisible(gc(full = TRUE, reset = TRUE))
  t0 <- Sys.time()
  value <- force(expr)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  peak <- sum(gc(full = TRUE)[, "max used"] * c(56, 8)) / 2^30
  list(value = value, secs = secs, peak_gib = peak)
}

rows <- list()
record <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...)

for (nm in names(DATASETS)) {
  ds <- DATASETS[[nm]]
  mtx_fp <- file.path(IN, ds, "cleanser", "grna_matrix.mtx")
  odm_fp <- file.path(IN, ds, "sceptre-pipeline", "grna.odm")
  cat(sprintf("\n================ %s ================\n", nm))

  cat("reading .mtx ...\n")
  r_mtx <- timed(fishashplus::read_grna_matrix(mtx_fp, format = "mtx"))
  counts_mtx <- r_mtx$value
  cat(sprintf("  %d guides x %d cells, nnz = %d  [%.1f s, %.2f GiB]\n",
              nrow(counts_mtx), ncol(counts_mtx), nnzero(counts_mtx),
              r_mtx$secs, r_mtx$peak_gib))
  record(dataset = nm, step = "read", input = "mtx",
         seconds = r_mtx$secs, peak_gib = r_mtx$peak_gib)

  counts_odm <- NULL
  if (!skip_odm) {
    if (!file.exists(odm_fp)) {
      cat("  no .odm at", odm_fp, "-- skipping\n")
    } else {
      cat("reading .odm (one feature at a time; this is the slow one) ...\n")
      r_odm <- timed(fishashplus::read_grna_matrix(odm_fp, format = "odm"))
      counts_odm <- r_odm$value
      cat(sprintf("  [%.1f s, %.2f GiB]  %.1fx the .mtx read\n",
                  r_odm$secs, r_odm$peak_gib, r_odm$secs / r_mtx$secs))
      record(dataset = nm, step = "read", input = "odm",
             seconds = r_odm$secs, peak_gib = r_odm$peak_gib)

      # The comparison below is only meaningful if the two formats agree.
      a <- counts_mtx; b <- counts_odm
      dimnames(a) <- dimnames(b) <- NULL
      stopifnot(identical(dim(a), dim(b)), identical(a@i, b@i),
                identical(a@p, b@p), isTRUE(all.equal(a@x, b@x)))
      cat("  .mtx and .odm hold the same matrix\n")
    }
  }

  # fishash indexes rownames(counts) to build its assignment strings, so a matrix
  # without them errors inside tapply. fishash+ passes dimnames through untouched.
  named <- counts_mtx
  rownames(named) <- paste0("grna_", seq_len(nrow(named)))
  colnames(named) <- paste0("CELL_", seq_len(ncol(named)))

  cat("running fishash (refit = 10) ...\n")
  f <- timed(fishash(named, refit = 10, padj_cutoff = 0.05, exclude_empty = TRUE))
  n_f <- sum(SummarizedExperiment::assay(f$value, "assigned"))
  cat(sprintf("  [%.1f s, %.2f GiB]  %d assignments\n", f$secs, f$peak_gib, n_f))
  record(dataset = nm, step = "fishash", input = "mtx",
         seconds = f$secs, peak_gib = f$peak_gib, n_assigned = n_f)

  for (fmt in c("mtx", "odm")) {
    m <- if (fmt == "mtx") counts_mtx else counts_odm
    if (is.null(m)) next
    cat(sprintf("running fishash+ on the %s matrix ...\n", fmt))
    p <- timed(fishash_plus(m, q = 0.05, refit = 10, min_count = 2))
    n_p <- sum(p$value)
    cat(sprintf("  [%.1f s, %.2f GiB]  %d assignments\n", p$secs, p$peak_gib, n_p))
    record(dataset = nm, step = "fishash_plus", input = fmt,
           seconds = p$secs, peak_gib = p$peak_gib, n_assigned = n_p)
  }

  rm(counts_mtx, counts_odm, named); invisible(gc())
}

res <- do.call(rbind, lapply(rows, function(r) {
  if (is.null(r$n_assigned)) r$n_assigned <- NA_integer_
  r
}))
res$seconds  <- round(res$seconds, 2)
res$peak_gib <- round(res$peak_gib, 2)
cat("\n=== summary ===\n")
print(res, row.names = FALSE)

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
out <- file.path(dirname(normalizePath(sub("^--file=", "", .self))),
                 sprintf("fishash-variant-timings-%s.csv", format(Sys.Date(), "%Y%m%d")))
write.csv(res, out, row.names = FALSE)
cat("\nwrote", out, "\n")
