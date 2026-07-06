#!/usr/bin/env Rscript
# make_barnyard.R -- build the benchmarking barnyard guide matrices from the RAW import-liu-2025 output.
#
# QC applied here (the project repo, NOT the import repo): a species-PURITY gate. Each mix sample is a
# 293T (human) + NIH3T3 (mouse) barnyard; a cell's frac_homo is its human-gene UMI fraction (from the
# human/mouse GENE rows, staged as cell_species.csv by import-liu-2025). Keep only species-PURE cells
# (frac_homo > 0.9 OR < 0.1), dropping the 0.1-0.9 doublets/ambiguous. No mito/feature/depth gate.
#
# Reproduces input_data/barnyard_*/sceptre/grna_matrix.rds byte-identically (verified: 8265->7453,
# 8959->8793, 11512->8902, 8207->8057). Companion to make_gasperini.R / make_replogle-rd7.R.
suppressMessages(library(Matrix))

IN  <- paste0(.get_config_path("LOCAL_LIU_2025_DATA_DIR"), "processed/")            # raw import output
OUT <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/")
keys <- c("barnyard_lrb100_0hr", "barnyard_lrb100_72hr", "barnyard_mch2_0hr", "barnyard_mch2_72hr")

for (key in keys) {
  grna <- readRDS(paste0(IN, key, "/grna_matrix.rds"))                              # 200 x all cells (raw)
  sp   <- read.csv(paste0(IN, key, "/cell_species.csv"))
  keep <- !is.na(sp$species)                                                        # species-purity QC
  out  <- as(grna[, keep, drop = FALSE], "CsparseMatrix")
  od <- paste0(OUT, key, "/sceptre/"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, paste0(od, "grna_matrix.rds"))
  cat(sprintf("%-22s %d -> %d species-pure cells\n", key, ncol(grna), ncol(out)))
}
