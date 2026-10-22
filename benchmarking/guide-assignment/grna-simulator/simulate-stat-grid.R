#!/usr/bin/env Rscript
# Generate the MOI x signal grid for statistical benchmarking.
#
#   sbatch -c 4 -t 4:00:00 run-in-image.sh simulate-stat-grid.R
#
# The design, the derivation of the simulator's inputs and the simulator call are
# in stat-grid-lib.R, shared with eda-stat-grid.R (which chose the background
# settings by looking at per-guide histograms).
#
# Two sets of datasets:
#   main      every (moi, s, rep), Phi_noise = 0. Seed = base_seed + row of the
#             expand.grid(rep, s, moi) grid, so adding a level at the END of a
#             STAT_GRID vector keeps the existing datasets' seeds.
#   phi_noise STAT_GRID$phi_noise_variant: its reps of every (moi, s) with ambient
#             overdispersed per entry, reusing the matching main dataset's seed
#             (paired: identical truth and signal, different ambient).
# A dataset whose sim_params.rds already exists is skipped, so levels can be added
# without regenerating, or overwriting the inputs of, datasets already in use.
#
# Dataset ids start with "stat_grid" so run_cleanser.py can pick its model (--dc).
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/input_data/:
#   stat_grid_moi<m>_s<s>_rep<NN>/, <phi prefix>_moi<m>_s<s>_rep<NN>/
#       inputs for STAT_GRID$methods, true_pert_matrix.rds, sim_params.rds
#   stat_grid_manifest.csv   design values and realized levels per dataset; rows
#                            for skipped datasets are carried over from the old file

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
dp    <- normalizePath(file.path(here, "..", "data-preprocessing"))
source(file.path(here, "grna-sim-utils.R"))             # write_method_inputs()
source(file.path(dp, "convert_odm_to_h5ad.R"))           # R_to_h5ad()
source(file.path(dp, "lib_make_guide_data.R"))           # write_h5ad_methods(), write_cleanser_method()
source(file.path(here, "stat-grid-lib.R"))               # STAT_GRID, simulate_grid(), realized(), name_dims()

out_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data")
B <- STAT_GRID$ambient_per_cell

main <- expand.grid(rep = seq_len(STAT_GRID$n_reps), s = STAT_GRID$signal_per_integration,
                    moi = STAT_GRID$moi)
main$seed      <- STAT_GRID$base_seed + seq_len(nrow(main))
main$phi_noise <- STAT_GRID_PARAMS$Phi_noise
main$dataset   <- sprintf("stat_grid_moi%g_s%g_rep%02d", main$moi, main$s, main$rep)

pv  <- STAT_GRID$phi_noise_variant
phi <- main[main$rep <= pv$reps, ]
phi$phi_noise <- pv$Phi_noise
phi$dataset   <- sprintf("%s_moi%g_s%g_rep%02d", pv$prefix, phi$moi, phi$s, phi$rep)

grid <- rbind(main, phi)
stopifnot(!anyDuplicated(grid$dataset))

mf <- file.path(out_root, "stat_grid_manifest.csv")
old_manifest <- if (file.exists(mf)) read.csv(mf) else NULL
if (!is.null(old_manifest) && !"phi_noise" %in% names(old_manifest))
  old_manifest$phi_noise <- 0                            # written before the Phi_noise variant

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  ds_dir <- file.path(out_root, g$dataset)
  if (file.exists(file.path(ds_dir, "sim_params.rds"))) {
    prev <- old_manifest[!is.null(old_manifest) & old_manifest$dataset == g$dataset, ]
    if (is.null(old_manifest) || nrow(prev) != 1)
      stop(g$dataset, " exists but has no row in ", mf)
    cat("skipping", g$dataset, "(exists)\n")
    rows[[length(rows) + 1L]] <- prev
    next
  }
  t0 <- Sys.time()
  r <- simulate_grid(g$moi, g$s, B, g$seed, Phi_noise = g$phi_noise)
  cat(sprintf("\n=== %s: Phi_noise %g, rho ~ Beta(%.3f, %.3f), mean %.4f, seed %d ===\n", g$dataset,
              g$phi_noise, r$inputs$rho_alpha, r$inputs$rho_beta, r$inputs$rho_mean, g$seed))
  counts <- name_dims(assay(r$sim, "counts"))
  truth  <- name_dims(assay(r$sim, "ground_truth"))

  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  write_method_inputs(counts, ds_dir, STAT_GRID$methods)
  saveRDS(truth, file.path(ds_dir, "true_pert_matrix.rds"))
  params <- STAT_GRID_PARAMS; params$Phi_noise <- g$phi_noise
  saveRDS(c(params, r$inputs,
            list(design = "stat_grid", simulator = "simulate_guidebender",
                 n_guides = STAT_GRID$n_guides, n_cells = STAT_GRID$n_cells,
                 hurdle_prob = STAT_GRID$hurdle_prob, moi = g$moi,
                 signal_per_integration = g$s, ambient_per_cell = B,
                 frac_noise_endo = STAT_GRID$frac_noise_endo,
                 rep = g$rep, seed = g$seed, methods = STAT_GRID$methods)),
          file.path(ds_dir, "sim_params.rds"))

  rows[[length(rows) + 1L]] <- cbind(
    data.frame(dataset = g$dataset, moi = g$moi, signal_per_integration = g$s,
               ambient_per_cell = B, rho_alpha = r$inputs$rho_alpha,
               rho_beta = r$inputs$rho_beta, rho_mean = r$inputs$rho_mean,
               expected_total_per_cell = r$inputs$expected_total_per_cell, seed = g$seed),
    realized(r$sim),
    data.frame(gen_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
               phi_noise = g$phi_noise))
  print(signif(Filter(is.numeric, rows[[length(rows)]]), 4), row.names = FALSE)
  rm(r, counts, truth); invisible(gc())
}

manifest <- do.call(rbind, rows)
write.csv(manifest, mf, row.names = FALSE)
cat("\nmanifest ->", mf, "(", nrow(manifest), "datasets )\n")

agg <- aggregate(manifest[, c("true_entries_per_cell", "signal_per_true_entry",
                              "ambient_per_cell_realized", "umis_per_cell_median",
                              "frac_true_count_lt2")],
                 by = manifest[, c("phi_noise", "moi", "signal_per_integration")], FUN = mean)
print(signif(agg, 4), row.names = FALSE)
