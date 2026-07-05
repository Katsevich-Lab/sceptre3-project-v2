#!/usr/bin/env Rscript
# eda/cells_over_target_sweep.R
#
# Sizing EDA for the computational benchmark: how does the candidate cell count
# grow with #targets? For each num_targets we sample targets exactly as the
# generator does (select_targets_random, excluding NT/off-target/unknown), take
# their guides (+ all NT guides in low-MOI, per force_nt_inclusion), and count
# candidate cells = cells expressing >=1 of those guides in the assignment matrix.
#
# This is the FRONT of select_comp() with the response-gene + write steps removed
# -- i.e. the beginning of a generator. #genes never affects #cells, so it isn't
# swept here. Output: a per-(dataset, num_targets, rep) table + a summary, so you
# can see the cell curve (and its spread across sampled target sets) per dataset.

library(tidyverse)
library(Matrix)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "..", "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) source(file.path(lib, f))

# ---- config ---------------------------------------------------------------
DATASETS         <- list(dataset_gasperini, dataset_replogle)
NUM_TARGETS_GRID <- c(25, 50, 100, 200, 400, 800, 1600, 3200)
N_REPS           <- 5                 # sampled target set varies -> replicate to see spread
NT_NAME          <- "non-targeting"
IGNORE           <- c(NT_NAME, "nt_off_target", "unknown")
out_csv <- file.path(script_dir, "cells_over_target_sweep.csv")

# ---- one measurement (front of select_comp) -------------------------------
count_cells <- function(assn, grna_target_df, dataset, num_targets, seed) {
  targets <- select_targets_random(grna_target_df, num_targets, exclude = IGNORE, random_seed = seed)
  if (length(targets) < num_targets) return(NULL)               # ran out of targets

  target_guides <- grna_target_df$grna_id[grna_target_df$grna_target %in% targets]
  guides <- if (dataset$moi$force_nt_inclusion) {               # low-MOI folds in the NT control cells
    union(target_guides, grna_target_df$grna_id[grna_target_df$grna_target == NT_NAME])
  } else target_guides

  ncells <- sum(Matrix::colSums(assn[guides, , drop = FALSE]) > 0)
  data.frame(num_targets = num_targets, nguides = length(guides), ncells = ncells)
}

# ---- sweep ----------------------------------------------------------------
rows <- list()
for (dataset in DATASETS) {
  cat("\n================ ", dataset$name, " ================\n")
  src <- read_source_data(dataset$name)
  stopifnot(src$low_moi == dataset$moi$low_moi)
  al  <- load_assignment(dataset$name, dataset$assignment_method, src$cell_covariates,
                         excluded_required = dataset$excluded_required)

  n_usable    <- length(al$usable_cells)
  n_targeting <- length(setdiff(unique(src$grna_target_df$grna_target), IGNORE))
  n_nt_cells  <- if (dataset$moi$force_nt_inclusion) {
    sum(Matrix::colSums(al$assn[src$grna_target_df$grna_id[src$grna_target_df$grna_target == NT_NAME], ,
                                drop = FALSE]) > 0)
  } else 0L
  cat("usable cells:", n_usable, "| targeting targets:", n_targeting,
      "| NT-cell floor (low-MOI):", n_nt_cells, "\n")

  grid <- NUM_TARGETS_GRID[NUM_TARGETS_GRID <= n_targeting]
  for (nt in grid) for (rep in seq_len(N_REPS)) {
    r <- count_cells(al$assn, src$grna_target_df, dataset, nt, seed = 1000L * rep + nt)
    if (!is.null(r)) {
      r$dataset    <- dataset$name
      r$rep        <- rep
      r$n_usable   <- n_usable
      r$frac_cells <- r$ncells / n_usable
      rows[[length(rows) + 1]] <- r
    }
  }
}

results <- dplyr::bind_rows(rows)
readr::write_csv(results, out_csv)
cat("\nWrote", nrow(results), "rows to", out_csv, "\n\n")

# ---- summary --------------------------------------------------------------
results |>
  dplyr::group_by(dataset, num_targets) |>
  dplyr::summarize(
    nguides_med = median(nguides),
    ncells_med  = median(ncells),
    ncells_min  = min(ncells),
    ncells_max  = max(ncells),
    frac_med    = round(median(frac_cells), 3),
    .groups = "drop"
  ) |>
  print(n = 200)
