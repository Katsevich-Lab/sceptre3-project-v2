# ============================================================================
# scripts/datasets.R — shared registry + loaders for the ambient-ceiling / method-comparison suite.
#
# Consolidates the ~17-dataset path list, the grna-matrix loader, the CLEANSER <=2 depth idiom, and
# the fishash / fishash+ contingency wrappers that were previously copy-pasted across many scripts.
#
# Source it AFTER scripts/contingency_method.R if you use the fishash*_assign() wrappers, e.g.
#     source("scripts/contingency_method.R"); source("scripts/datasets.R")
# Kept dependency-light: only Matrix is needed for the loaders / kappa_leq2().
#
# Data locations resolve through the lab config (~/.research_config via .get_config_path in ~/.Rprofile),
# NOT hardcoded ~/data paths -- matches the rest of the repo (e.g. association/data-preparation/lib/io.R).
# The config vars end in "/", so concatenate with paste0 (file.path would insert a double slash).
# ============================================================================

.GA_EXT <- .get_config_path("LOCAL_EXTERNAL_DATA_DIR")                       # ~/data/external/
.GA_IN  <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
.GA_SUR <- paste0(.GA_EXT, "perturbseq-survey")

# The 17 real gRNA-count datasets (guides x cells). Named list of grna_matrix.rds paths.
dataset_paths <- function() list(
  barnyard_lrb100_0hr  = file.path(.GA_IN, "barnyard_lrb100_0hr/sceptre/grna_matrix.rds"),
  barnyard_lrb100_72hr = file.path(.GA_IN, "barnyard_lrb100_72hr/sceptre/grna_matrix.rds"),
  barnyard_mch2_0hr    = file.path(.GA_IN, "barnyard_mch2_0hr/sceptre/grna_matrix.rds"),
  barnyard_mch2_72hr   = file.path(.GA_IN, "barnyard_mch2_72hr/sceptre/grna_matrix.rds"),
  gasperini            = file.path(.GA_IN, "gasperini/sceptre/grna_matrix.rds"),
  `replogle-rd7`       = file.path(.GA_IN, "replogle-rd7/sceptre/grna_matrix.rds"),
  mccutcheon           = paste0(.GA_EXT, "mccutcheon-2023-GSE218988/grna_matrix.rds"),
  a549                 = file.path(.GA_SUR, "a549_crispri_perturbseq_GSE236304/grna_matrix.rds"),
  cd8_tcell            = file.path(.GA_SUR, "cd8_tcell_perturbcite_GSE279498/grna_matrix.rds"),
  dctap_k562_highmoi   = file.path(.GA_SUR, "dctap_k562_highmoi_GSE303901/grna_matrix.rds"),
  dctap_k562_lowmoi    = file.path(.GA_SUR, "dctap_k562_lowmoi_GSE303901/grna_matrix.rds"),
  endoc_t2d            = file.path(.GA_SUR, "endoc_t2d_perturbseq_GSE273677/grna_matrix.rds"),
  gastric_organoid     = file.path(.GA_SUR, "gastric_organoid_cropseq_GSE280506/grna_matrix.rds"),
  invivo_cortex        = file.path(.GA_SUR, "invivo_cortex_perturbseq_GSE249416/grna_matrix.rds"),
  ipsc                 = file.path(.GA_SUR, "ipsc_crispri_hipsci_figshare27989294/grna_matrix.rds"),
  erythroid_multiome   = file.path(.GA_SUR, "perturb_multiome_erythroid_GSE274113/grna_matrix.rds"),
  cd4_tcell            = file.path(.GA_SUR, "tcell_cd4_perturbseq_GSE314342/grna_matrix.rds"))

# DC-TAP datasets also ship a sceptre-processed dir: aligned gRNA matrix + response (GEX) matrix +
# grna_target_data_frame.csv (used by the positive-control / per-guide analyses).
dctap_sceptre_dir <- function(which = c("highmoi", "lowmoi")) {
  which <- match.arg(which)
  sub <- if (which == "highmoi") "dctap_k562_highmoi_GSE303901" else "dctap_k562_lowmoi_GSE303901"
  file.path(.GA_SUR, sub, "sceptre")
}

# Load a guide x cell count matrix by dataset NAME (from dataset_paths()) or by explicit path,
# returned as a numeric CsparseMatrix (the form every downstream script expects).
# NOTE: several matrices are stored on disk as dgRMatrix (row-compressed), so the
# as(...,"CsparseMatrix") is a load-bearing class conversion, not cosmetic -- a bare readRDS()
# would return the wrong class. (storage.mode(@x)<-"double" IS a no-op; @x is already double.)
load_grna_matrix <- function(ds_or_path) {
  reg <- dataset_paths()
  p <- if (ds_or_path %in% names(reg)) reg[[ds_or_path]] else ds_or_path
  if (is.null(p) || !file.exists(p)) stop("cannot resolve gRNA matrix: ", ds_or_path)
  m <- as(readRDS(p), "CsparseMatrix"); storage.mode(m@x) <- "double"; m
}

# Directory holding the replogle-rd7 sceptre-pipeline DE artifacts (response.odm,
# grna.odm, sceptre_object.rds).  Single source of truth for the several scripts
# that read those files directly.
replogle_rd7_pipeline_dir <- function() file.path(.GA_IN, "replogle-rd7", "sceptre-pipeline")

# Replogle-rd7 with its downstream sceptre-DE artifacts: the gRNA count matrix, the
# response ODM, and the fitted sceptre_object from the sceptre-pipeline run.  The
# collaborator dose-response / ENO1 figures share this triple.  Requires ondisc +
# sceptre to be attached by the caller (for initialize_odm_from_backing_file).
load_replogle_rd7_de <- function() {
  Dp <- replogle_rd7_pipeline_dir()
  list(mc   = load_grna_matrix("replogle-rd7"),
       resp = initialize_odm_from_backing_file(file.path(Dp, "response.odm")),
       so   = readRDS(file.path(Dp, "sceptre_object.rds")))
}

# CLEANSER-style per-cell ambient depth L_i: sum of the cell's <=2 counts.
kappa_leq2 <- function(gm) { s <- gm; s@x[s@x > 2] <- 0; Matrix::colSums(s) }

# Thin wrappers over contingency_assign with the finalized config (q=0.05, refit=10, min_count=2,
# tail="hyper", fdr="GS"). fishash = observed cell margin; fishash+ = ambient (denoised) cell margin.
# Requires scripts/contingency_method.R to have been sourced.
fishash_assign     <- function(gm, ...) contingency_assign(gm, cell_margin = "observed", tail = "hyper", fdr = "GS", ...)
fishashplus_assign <- function(gm, ...) contingency_assign(gm, cell_margin = "ambient",  tail = "hyper", fdr = "GS", ...)

# CLEANSER MCMC assignments (produced by scripts/run_cleanser_batch.sh) live here — gitignored,
# regenerable. Read a <ds>.csv (columns guide,cell) into a logical assignment matrix aligned to gm.
CLEANSER_CALLS_DIR <- "results/ambient_ceiling/cleanser_calls"
load_cleanser_assignment <- function(gm, ds) {
  f <- file.path(CLEANSER_CALLS_DIR, paste0(ds, ".csv"))
  if (!file.exists(f)) stop("no CLEANSER calls for '", ds, "' at ", f, " (run scripts/run_cleanser_batch.sh)")
  calls <- utils::read.csv(f, stringsAsFactors = FALSE)
  gi <- match(calls$guide, rownames(gm)); ci <- match(calls$cell, colnames(gm)); ok <- !is.na(gi) & !is.na(ci)
  Matrix::sparseMatrix(i = gi[ok], j = ci[ok], x = TRUE, dims = dim(gm), dimnames = dimnames(gm))
}
