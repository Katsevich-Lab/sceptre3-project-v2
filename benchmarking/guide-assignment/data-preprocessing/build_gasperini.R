#!/usr/bin/env Rscript
# ============================================================================
# DORMANT (paused 2026-08-30). Part of the NESTED GUIDE-SUBSET benchmark:
# datasets of 100/200/400/800 guides selected by perturbed sample size, built to
# measure how method memory and runtime SCALE with guide count.
#
# That is a different question from the one currently being run, which is "which
# methods can handle the ENTIRE grna matrix at all" -- see build_full_h5ad.R.
# Nothing outside this directory references these files.
#
# Paused, not abandoned: this code was unit- and integration-tested and works.
# ============================================================================
# build_gasperini.R
# Thin driver: build the nested guide-assignment input datasets for Gasperini.
# crispat/pertpy/cleanser only (sceptre deferred). Replaces old make_gasperini.R.

rm(list = ls())
source("~/.Rprofile")                      # provides .get_config_path()

# Resolve own dir under both `Rscript` (--file=) and interactive source().
.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "convert_odm_to_h5ad.R"))     # R_to_h5ad
source(file.path(script_dir, "lib_make_guide_data.R"))

build_guide_assignment_datasets(
  source_name = "gasperini",
  out_prefix  = "gasperini",
  threshold   = 5,
  K           = 10,
  sizes       = c(100, 200, 400, 800),   # max(sizes) is the anchor (final dataset)
  methods     = c("crispat", "pertpy", "cleanser")
)
