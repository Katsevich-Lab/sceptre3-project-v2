#!/usr/bin/env Rscript

# Reproduce the Fishash Table 2 SCEPTRE-mixture assignment for one Liu et al.
# barnyard sample. This is a direct, file-oriented adaptation of
# jackkamm/fishash_analysis bin/run_sceptre_mixture.R. Unlike the simulation
# shortcut elsewhere in this repository, it uses the matched GEX count matrix.

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(sceptre))

parse_args <- function(x) {
  stopifnot(length(x) %% 2L == 0L)
  keys <- sub("^--", "", x[seq(1L, length(x), by = 2L)])
  vals <- x[seq(2L, length(x), by = 2L)]
  stats::setNames(as.list(vals), keys)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("raw-prefix", "out", "cpus")
missing <- setdiff(required, names(args))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

if (as.character(packageVersion("sceptre")) != "0.10.3") {
  stop("Table 2 requires sceptre 0.10.3; found ", packageVersion("sceptre"))
}

prefix <- args[["raw-prefix"]]
features <- read.delim(
  gzfile(paste0(prefix, "_features.tsv.gz")), header = FALSE,
  col.names = c("ID", "Symbol", "type"), check.names = FALSE
)
barcodes <- scan(
  gzfile(paste0(prefix, "_barcodes.tsv.gz")), what = character(), quiet = TRUE
)
counts <- as(readMM(gzfile(paste0(prefix, "_matrix.mtx.gz"))), "CsparseMatrix")
stopifnot(nrow(counts) == nrow(features), ncol(counts) == length(barcodes))

rownames(counts) <- make.unique(features$Symbol)
colnames(counts) <- barcodes
guide_counts <- counts[features$type == "CRISPR Guide Capture", , drop = FALSE]
gex_counts <- counts[features$type == "Gene Expression", , drop = FALSE]
stopifnot(nrow(guide_counts) == 200L)

# Match the authors' deterministic seed. Nextflow staged the processed RDS with
# the batch-derived basename, so run_sceptre_mixture.R saw this filename.
batch <- sub("^.*/GSE272457_", "", prefix)
seed_key <- paste0("run_sceptre_mixture.R_", batch, ".Rds")
set.seed(abs(digest::digest2int(seed_key)))

sceptre_object <- sceptre::import_data(
  response_matrix = gex_counts,
  grna_matrix = guide_counts,
  grna_target_data_frame = data.frame(
    grna_id = rownames(guide_counts), grna_target = "non-targeting"
  ),
  moi = "high"
)
sceptre_object <- sceptre::set_analysis_parameters(sceptre_object)

cpus <- as.integer(args$cpus)
if (cpus == 1L) {
  sceptre_object <- sceptre::assign_grnas(
    sceptre_object, method = "mixture", parallel = FALSE
  )
} else {
  sceptre_object <- sceptre::assign_grnas(
    sceptre_object, method = "mixture", parallel = TRUE,
    n_processors = cpus
  )
}

assignment <- sceptre::get_grna_assignments(sceptre_object)[rownames(guide_counts), ]
colnames(assignment) <- colnames(guide_counts)
writeMM(assignment, args$out)

cat(sprintf(
  "sceptre=%s sample=%s guides=%d cells=%d calls=%d seed_key=%s\n",
  packageVersion("sceptre"), batch, nrow(assignment), ncol(assignment),
  sum(assignment), seed_key
))
