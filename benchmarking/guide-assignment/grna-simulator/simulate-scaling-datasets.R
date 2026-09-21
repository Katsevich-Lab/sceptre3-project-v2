#!/usr/bin/env Rscript
# Generate the computational-scaling datasets for one regime.
#
#   Rscript simulate-scaling-datasets.R gasperini                   # every rung
#   Rscript simulate-scaling-datasets.R replogle 1 2                # rungs 1 and 2 only
#   Rscript simulate-scaling-datasets.R gasperini --calibrate-only  # measure + calibrate
#
# Needs the R_461 renv library (fishash, zellkonverter, Matrix):
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# Measures the real matrix, calibrates snr and lambda against it, then generates
# each rung. Regimes, methods, fixed parameters and every step are defined in the
# "GuideBender scaling datasets" section of grna-sim-utils.R.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/input_data/:
#   <dataset>/                              one directory per rung
#   sim_scaling_calibration_<regime>.csv    calibrated lambda and snr
#   sim_scaling_manifest_<regime>.csv       real matrix and every rung, measured alike

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
dp    <- normalizePath(file.path(here, "..", "data-preprocessing"))
source(file.path(here, "grna-sim-utils.R"))
source(file.path(dp, "convert_odm_to_h5ad.R"))           # R_to_h5ad()
source(file.path(dp, "lib_make_guide_data.R"))           # write_h5ad_methods(), write_cleanser_method()

SEED_CAL <- 11L
SEED_GEN <- 1L

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% names(SCALING_REGIMES))
  stop("usage: simulate-scaling-datasets.R <", paste(names(SCALING_REGIMES), collapse = "|"),
       "> [rung numbers] [--calibrate-only]")
name <- args[1]
reg  <- SCALING_REGIMES[[name]]
calibrate_only <- "--calibrate-only" %in% args
pick <- suppressWarnings(as.integer(setdiff(args[-1], "--calibrate-only")))

out_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data")


# ---- 1. measure ------------------------------------------------------------
cat(sprintf("\n=== %s: measuring %s at threshold %d ===\n", name, reg$real_dataset, reg$threshold))
tg <- measure_real(reg)
cat(sprintf(paste0("  n_cells/n_guides %.3f | moi %.3f | nnz/cell %.3f | median UMIs/cell %.0f\n",
                   "  hurdle_prob %.4f | guide_infection_alpha %.3f | snr start %.2f\n"),
            tg$n_over_g, tg$moi, tg$nnz_per_cell, tg$count_per_cell,
            tg$hurdle_prob, tg$alpha, tg$snr_start))

rungs <- regime_rungs(reg, tg$n_over_g, name)
if (length(pick)) rungs <- rungs[pick, , drop = FALSE]
print(rungs, row.names = FALSE)


# ---- 2. calibrate ----------------------------------------------------------
t0  <- Sys.time()
cal <- calibrate(reg, tg, rungs, SEED_CAL)
cat(sprintf("  calibration took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cal_fp <- file.path(out_root, sprintf("sim_scaling_calibration_%s.csv", name))
write.csv(cal, cal_fp, row.names = FALSE)
cat("calibration ->", cal_fp, "\n")
if (calibrate_only) quit(save = "no")


# ---- 3. generate -----------------------------------------------------------
manifest <- list(tg$stats[, SCALING_MANIFEST_COLS])
for (i in seq_len(nrow(rungs))) {
  r <- rungs[i, ]
  k <- cal$role == "rung" & cal$dataset == r$dataset
  cat(sprintf("\n=== %s (%d guides x %d cells, lambda %.4f, snr %.3f) ===\n",
              r$dataset, r$n_guides, r$n_cells, cal$lambda[k], cal$snr[k]))
  manifest[[i + 1L]] <- generate_rung(name, reg, tg, r$dataset, r$n_guides, r$n_cells,
                                      cal$lambda[k], cal$snr[k], SEED_GEN, out_root)
  invisible(gc())
}

# The real matrix's row has no truth-based or generation columns.
extra <- setdiff(names(manifest[[length(manifest)]]), SCALING_MANIFEST_COLS)
manifest[[1]][extra] <- NA_real_
manifest <- do.call(rbind, manifest)
mf <- file.path(out_root, sprintf("sim_scaling_manifest_%s.csv", name))
write.csv(manifest, mf, row.names = FALSE)

cat("\n=== manifest (first column is the real matrix) ===\n")
show <- manifest
num  <- vapply(show, is.numeric, logical(1))
show[num] <- lapply(show[num], signif, 4)
rownames(show) <- show$label
print(t(show[, -1]), quote = FALSE)
cat("\nmanifest ->", mf, "\n")
