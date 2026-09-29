#!/usr/bin/env Rscript
# Generate resampled computational-scaling datasets. PROTOTYPE.
#
#   Rscript make-resampled-datasets.R <gasperini|replogle-rd7> \
#     --placement <uniform|weighted|none> --n-cells 20000,40000 --n-guides 500,1000 \
#     [--seed 1] [--method-inputs] [--overwrite] [--out-root DIR]
#
# Every combination of --n-cells and --n-guides is generated. --placement has no
# default, so every dataset's guide placement is a stated choice. On Betty, run it
# inside the image with run-in-image.sh. The generator is resample_counts() in
# resample-utils.R.
#
# Writes to <out-root>/<dataset>/, where out-root defaults to
# <LOCAL_BENCHMARKING_DIR>/guide_assignment/input_data:
#   counts.rds        guides x cells dgCMatrix, for analysis
#   sim_params.rds    source, placement, n_cells, n_guides, seed, the real column
#                     each simulated cell came from, and (weighted) the real guide
#                     each simulated guide took its weight from, and the weights
#   stats.csv         guide_matrix_stats() of the simulated matrix, and of the
#                     same real cells before their guides were redrawn
#   cleanser/, crispat/, pertpy/   each method's input, only with --method-inputs
# An existing dataset directory is an error unless --overwrite is given.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
dp    <- normalizePath(file.path(here, "..", "data-preprocessing"))
source(file.path(here, "grna-sim-utils.R"))              # guide_matrix_stats(), write_method_inputs()
source(file.path(here, "resample-utils.R"))
source(file.path(dp, "convert_odm_to_h5ad.R"))           # R_to_h5ad()
source(file.path(dp, "lib_make_guide_data.R"))           # write_h5ad_methods(), write_cleanser_method()


# ---- arguments ---------------------------------------------------------------
usage <- paste("usage: make-resampled-datasets.R <", paste(names(RESAMPLE_SOURCES), collapse = "|"),
               "> --placement <uniform|weighted|none> --n-cells N[,N...] --n-guides G[,G...]",
               "[--seed S] [--method-inputs] [--overwrite] [--out-root DIR]")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% names(RESAMPLE_SOURCES)) stop(usage)
source_name <- args[1]
opt <- function(flag, default = NULL) {
  k <- match(flag, args)
  if (is.na(k)) return(default)
  if (k == length(args)) stop(flag, " needs a value\n", usage)
  args[k + 1L]
}
int_list <- function(s) {
  v <- suppressWarnings(as.integer(strsplit(s, ",", fixed = TRUE)[[1]]))
  if (anyNA(v) || any(v <= 0)) stop("not a list of positive integers: ", s)
  v
}
if (is.null(opt("--n-cells")) || is.null(opt("--n-guides"))) stop(usage)
placement <- opt("--placement")
if (is.null(placement) || !placement %in% c("uniform", "weighted", "none"))
  stop("--placement must be uniform, weighted or none\n", usage)
n_cells_all  <- int_list(opt("--n-cells"))
n_guides_all <- int_list(opt("--n-guides"))
seed         <- as.integer(opt("--seed", "1"))
method_inputs <- "--method-inputs" %in% args
overwrite     <- "--overwrite" %in% args
out_root <- opt("--out-root", file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                                        "guide_assignment", "input_data"))
threshold <- RESAMPLE_SOURCES[[source_name]]$threshold

grid <- expand.grid(n_cells = n_cells_all, n_guides = n_guides_all)
grid$dataset <- resample_dataset_id(source_name, placement, grid$n_guides, grid$n_cells)
exists <- dir.exists(file.path(out_root, grid$dataset))
if (any(exists) && !overwrite)
  stop("already exist (pass --overwrite to replace): ",
       paste(grid$dataset[exists], collapse = ", "))


# ---- generate ----------------------------------------------------------------
real <- load_real_counts(source_name, "rds")
cat(sprintf("%s: %d guides x %d cells, %d nonzeros; busiest cell has %d nonzeros\n",
            source_name, nrow(real), ncol(real), length(real@x), max(diff(real@p))))

for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  cat(sprintf("\n=== %s (%s placement, seed %d) ===\n", g$dataset, placement, seed))
  t0 <- Sys.time()
  sim <- resample_counts(real, g$n_cells, g$n_guides, seed, placement)
  gen_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ds_dir <- file.path(out_root, g$dataset)
  if (dir.exists(ds_dir)) unlink(ds_dir, recursive = TRUE)   # only reached with --overwrite
  dir.create(ds_dir, recursive = TRUE)
  saveRDS(sim$counts, file.path(ds_dir, "counts.rds"))
  saveRDS(list(generator = "resample_counts", source = source_name, placement = placement,
               n_cells = g$n_cells, n_guides = g$n_guides, seed = seed,
               threshold_for_stats = threshold, source_cells = sim$source_cells,
               source_guides = sim$source_guides, guide_weights = sim$guide_weights),
          file.path(ds_dir, "sim_params.rds"))
  if (method_inputs) write_method_inputs(sim$counts, ds_dir)

  # The same real cells before their guides were redrawn: per-cell statistics
  # must match the simulated row exactly; per-guide ones show what redrawing did.
  st <- rbind(guide_matrix_stats(sim$counts, threshold, label = g$dataset),
              guide_matrix_stats(real[, sim$source_cells], threshold,
                                 label = paste0(g$dataset, "_source_cells")))
  st$gen_sec <- round(gen_sec, 2)
  write.csv(st, file.path(ds_dir, "stats.csv"), row.names = FALSE)
  cat(sprintf(paste0("  moi %.3f (source cells %.3f) | nnz/cell %.3f (%.3f) | ",
                     "perturbed cells/guide %.1f | %.1f s\n"),
              st$moi[1], st$moi[2], st$nnz_per_cell[1], st$nnz_per_cell[2],
              st$pert_per_guide[1], gen_sec))
  rm(sim); invisible(gc())
}
cat("\nwrote", nrow(grid), "dataset(s) under", out_root, "\n")
