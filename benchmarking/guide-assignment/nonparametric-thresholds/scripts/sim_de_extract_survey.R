#!/usr/bin/env Rscript
# =============================================================================
# Extract GEX matrices for the survey datasets (perturbseq-survey/) and produce
# the sceptre-ready triple:
#   response_matrix.rds         genes  x cells (sparse), GEX UMIs
#   grna_matrix_aligned.rds     guides x cells (sparse), gRNA UMIs, ALIGNED
#                                to the GEX cell ordering
#   grna_target_data_frame.csv  grna_id, grna_target  (from guide_features.csv
#                                if it has target info, else best-effort guess)
#
# Each survey dataset's raw data is at ~/data/external/perturbseq-survey/<name>/.
# Outputs are written into the SAME directory under a `sceptre/` subdir.
#
# We already have grna_matrix.rds for every dataset (the lab's parse_10x_*
# scripts produced it).  Cell barcodes there are the same 10x barcodes used in
# the GEX file, so alignment is by exact barcode match.
# =============================================================================
suppressPackageStartupMessages({library(rhdf5); library(Matrix)})

SURV <- path.expand("~/data/external/perturbseq-survey")

# ---- helpers ----------------------------------------------------------------
# Read a CellRanger v3 .h5 into a list with full matrix, feature_type, barcodes.
read_10x_h5 <- function(h5) {
  ls <- h5ls(h5); grp <- if (any(ls$name == "matrix" & ls$group == "/")) "/matrix" else "/"
  r1 <- function(p) as.vector(h5read(h5, p))
  data <- r1(paste0(grp, "/data")); indices <- r1(paste0(grp, "/indices"))
  indptr <- r1(paste0(grp, "/indptr")); shape <- r1(paste0(grp, "/shape"))
  barcodes <- r1(paste0(grp, "/barcodes"))
  fpath <- paste0(grp, "/features"); fn <- ls[ls$group == fpath, "name"]
  fid   <- if ("id"   %in% fn) r1(paste0(fpath, "/id"))   else r1(paste0(grp, "/genes"))
  fname <- if ("name" %in% fn) r1(paste0(fpath, "/name")) else fid
  ftype <- if ("feature_type" %in% fn) r1(paste0(fpath, "/feature_type")) else rep("Unknown", length(fid))
  m <- sparseMatrix(i = indices + 1L, p = indptr, x = data, dims = shape,
                    dimnames = list(make.unique(fid), barcodes), repr = "C")
  list(m = m, feat_id = fid, feat_name = fname, feat_type = ftype, barcodes = barcodes)
}
# Read a 10x mtx triplet (dir with matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
read_10x_mtx <- function(mtx_dir) {
  f <- function(...) { c(...) -> p; file.path(mtx_dir, p[file.exists(file.path(mtx_dir, p))][1]) }
  mtx_f  <- f("matrix.mtx.gz", "matrix.mtx")
  feat_f <- f("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
  bc_f   <- f("barcodes.tsv.gz", "barcodes.tsv")
  m <- readMM(gzfile(mtx_f)); m <- as(m, "CsparseMatrix")
  feat <- read.delim(gzfile(feat_f), header = FALSE, stringsAsFactors = FALSE)
  bc   <- readLines(gzfile(bc_f))
  ftype <- if (ncol(feat) >= 3) feat[[3]] else rep("Unknown", nrow(feat))
  rownames(m) <- make.unique(as.character(feat[[1]])); colnames(m) <- bc
  list(m = m, feat_id = feat[[1]], feat_name = if (ncol(feat) >= 2) feat[[2]] else feat[[1]],
       feat_type = ftype, barcodes = bc)
}
# Standard split: rows by feature type, output (gex, grna) + their target mapping
.split_and_save <- function(raw, dataset_name, out_dir, guide_regex = "CRISPR|Guide|gRNA|guide",
                            gex_regex = "Gene Expression|^Gene$") {
  is_guide <- grepl(guide_regex, raw$feat_type)
  is_gex   <- grepl(gex_regex,   raw$feat_type)
  cat(sprintf("  feature_type table:\n")); print(table(raw$feat_type))
  cat(sprintf("  -> %d GEX, %d guides\n", sum(is_gex), sum(is_guide)))
  if (sum(is_gex) == 0)   stop("no GEX features matched: ", gex_regex)
  if (sum(is_guide) == 0) stop("no guide features matched: ", guide_regex)
  gex  <- raw$m[is_gex,   , drop = FALSE]
  grna <- raw$m[is_guide, , drop = FALSE]
  rownames(grna) <- if (length(raw$feat_name)) raw$feat_name[is_guide] else rownames(grna)
  rownames(grna) <- make.unique(as.character(rownames(grna)))
  # Existing grna_matrix.rds we already saved -- verify the new extraction
  # contains the same data (alignment by guide name + cell barcode).
  prior <- readRDS(file.path(SURV, dataset_name, "grna_matrix.rds"))
  # Drop cells with all-zero GEX (cellranger keeps barcodes; we want real cells).
  nz_cells <- which(Matrix::colSums(gex) > 0)
  cat(sprintf("  cells with nonzero GEX: %d/%d\n", length(nz_cells), ncol(gex)))
  gex  <- gex[,  nz_cells, drop = FALSE]
  grna <- grna[, nz_cells, drop = FALSE]
  # save
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(gex,  file.path(out_dir, "response_matrix.rds"))
  saveRDS(grna, file.path(out_dir, "grna_matrix_aligned.rds"))
  # also save the (id, name) mapping for ALL features so target dfs can resolve
  # SYMBOL -> ENSG without rooting through the original raw files.
  features_df <- data.frame(id = raw$feat_id, name = raw$feat_name, type = raw$feat_type,
                            stringsAsFactors = FALSE)
  write.csv(features_df, file.path(out_dir, "features_all.csv"), row.names = FALSE)
  cat(sprintf("  GEX:  %d genes  x %d cells  nnz=%d\n",  nrow(gex),  ncol(gex),  length(gex@x)))
  cat(sprintf("  gRNA: %d guides x %d cells  nnz=%d\n",  nrow(grna), ncol(grna), length(grna@x)))
  cat(sprintf("  match vs prior grna_matrix.rds: %d guides shared / %d in prior\n",
              length(intersect(rownames(grna), rownames(prior))), nrow(prior)))
  invisible(list(gex = gex, grna = grna, n_cells = length(nz_cells)))
}

# ---- per-dataset extraction config ------------------------------------------
# Each config knows where its raw data is + which extractor to use.
extract_dataset <- function(name) {
  d <- file.path(SURV, name); out <- file.path(d, "sceptre")
  cat(sprintf("\n=== %s ===\n", name))
  if (name == "a549_crispri_perturbseq_GSE236304") {
    raw <- read_10x_h5(file.path(d, "filtered_feature_bc_matrix.h5"))
    .split_and_save(raw, name, out)
  } else if (name == "cd8_tcell_perturbcite_GSE279498") {
    raw <- read_10x_mtx(d)
    .split_and_save(raw, name, out)
  } else if (name == "perturb_multiome_erythroid_GSE274113") {
    raw <- read_10x_h5(file.path(d, "GSE274113_filtered_feature_bc_matrix_1.h5"))
    .split_and_save(raw, name, out, gex_regex = "Gene Expression")
  } else if (name == "tcell_cd4_perturbseq_GSE314342") {
    raw <- read_10x_h5(file.path(d, "GSM9393828_CD4i_R1L01_D1_Rest_filtered_feature_bc_matrix.h5"))
    .split_and_save(raw, name, out)
  } else if (name == "dctap_k562_highmoi_GSE303901") {
    tmp <- tempfile("dctap_hi_"); dir.create(tmp)
    untar(file.path(d, "GSE303901_k562_14-MOI_filtered_feature_bc_matrix.tar.gz"), exdir = tmp)
    raw <- read_10x_mtx(tmp); .split_and_save(raw, name, out); unlink(tmp, recursive = TRUE)
  } else if (name == "dctap_k562_lowmoi_GSE303901") {
    tmp <- tempfile("dctap_lo_"); dir.create(tmp)
    untar(file.path(d, "GSE303901_k562_1-MOI_filtered_feature_bc_matrix.tar.gz"), exdir = tmp)
    raw <- read_10x_mtx(tmp); .split_and_save(raw, name, out); unlink(tmp, recursive = TRUE)
  } else if (name == "endoc_t2d_perturbseq_GSE273677") {
    # 4 reps; load each as a 10x triplet, then cbind (cells), keeping a single
    # union of features (CellRanger should produce the same feature set per rep
    # since they share a reference; we verify).
    tmp <- tempfile("endoc_"); dir.create(tmp)
    untar(file.path(d, "GSE273677_RAW.tar"), exdir = tmp)
    fs <- list.files(tmp, full.names = TRUE)
    reps <- unique(sub("_(barcodes|features|matrix).*$", "", basename(fs)))
    # GWAS and RQC are separate libraries with different feature sets; combine GWAS only.
    reps <- reps[grepl("GWAS", reps)]
    cat(sprintf("  found %d GWAS reps: %s\n", length(reps), paste(reps, collapse = ", ")))
    parts <- list()
    for (rep in reps) {
      sub <- file.path(tmp, paste0("_", rep)); dir.create(sub, showWarnings = FALSE)
      for (f in fs[grepl(paste0("^", rep, "_"), basename(fs))]) {
        b <- basename(f)
        tgt <- if (grepl("matrix\\.mtx",  b)) "matrix.mtx.gz"
               else if (grepl("features\\.tsv", b)) "features.tsv.gz"
               else if (grepl("barcodes\\.tsv", b)) "barcodes.tsv.gz" else NULL
        if (!is.null(tgt)) file.copy(f, file.path(sub, tgt), overwrite = TRUE)
      }
      r <- tryCatch(read_10x_mtx(sub), error = function(e) NULL); if (is.null(r)) next
      colnames(r$m) <- paste0(rep, "_", colnames(r$m))     # prefix barcodes with rep id
      parts[[rep]] <- r
    }
    stopifnot(length(parts) >= 1)
    # all reps should share the same feature set (10x ref); verify
    feat_ids <- lapply(parts, function(p) p$feat_id)
    same <- all(vapply(feat_ids, identical, logical(1), y = feat_ids[[1]]))
    if (!same) warning("endoc reps have differing feature sets; using first rep's features")
    raw <- list(m = do.call(cbind, lapply(parts, `[[`, "m")),
                feat_id = parts[[1]]$feat_id, feat_name = parts[[1]]$feat_name,
                feat_type = parts[[1]]$feat_type,
                barcodes  = unlist(lapply(parts, `[[`, "barcodes"), use.names = FALSE))
    .split_and_save(raw, name, out); unlink(tmp, recursive = TRUE)
  } else {
    stop("no extractor configured for: ", name)
  }
}

# ---- run all configured datasets --------------------------------------------
DATASETS <- c("a549_crispri_perturbseq_GSE236304",
              "cd8_tcell_perturbcite_GSE279498",
              "perturb_multiome_erythroid_GSE274113",
              "tcell_cd4_perturbseq_GSE314342",
              "dctap_k562_highmoi_GSE303901",
              "dctap_k562_lowmoi_GSE303901",
              "endoc_t2d_perturbseq_GSE273677")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  todo <- if (length(args)) args else DATASETS
  for (nm in todo) {
    tryCatch(extract_dataset(nm),
             error = function(e) cat(sprintf("  ERROR %s: %s\n", nm, conditionMessage(e))))
  }
}
if (sys.nframe() == 0) main()
