#!/usr/bin/env Rscript
# Score the out-of-process methods (geomux, crispat) from their per-dataset
# *_calls.csv (guide,cell), and combine with scores_R.csv into scores_all.csv.
suppressPackageStartupMessages({ library(Matrix) })
source(file.path(getwd(), "scripts", "sim_lib.R"))
OUT <- SIMFW()
man <- read.csv(file.path(OUT, "manifest_all.csv"))

calls_to_matrix <- function(f, gnames, cnames) {
  if (!file.exists(f)) return(NULL)
  cl <- tryCatch(read.csv(f, colClasses = "character"), error = function(e) NULL)
  if (is.null(cl) || nrow(cl) == 0) return(Matrix::sparseMatrix(i = integer(0), j = integer(0),
      dims = c(length(gnames), length(cnames)), dimnames = list(gnames, cnames)))
  gi <- match(cl$guide, gnames); cj <- match(cl$cell, cnames)
  ok <- !is.na(gi) & !is.na(cj)
  Matrix::sparseMatrix(i = gi[ok], j = cj[ok], x = TRUE,
                       dims = c(length(gnames), length(cnames)), dimnames = list(gnames, cnames))
}

rows <- list()
for (k in seq_len(nrow(man))) {
  id <- man$id[k]; d <- load_dataset(id)
  counts <- as(d$counts, "CsparseMatrix")
  if (is.null(rownames(counts))) rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  tr <- truth_observed(d$Z, counts); dimnames(tr) <- dimnames(counts)
  for (mn in c("geomux", "crispat")) {
    A <- calls_to_matrix(file.path(d$dir, paste0(mn, "_calls.csv")), rownames(counts), colnames(counts))
    if (is.null(A)) next
    s <- score_assignment(.as_lgl(A), tr)$summary; s$method <- mn
    rows[[length(rows) + 1]] <- cbind(man[k, ], s[, c("method","jaccard","precision","recall",
                                                       "fdr_pooled","recall_pooled","n_pred","n_true")])
  }
}
ext <- if (length(rows)) do.call(rbind, rows) else NULL
scR <- read.csv(file.path(OUT, "scores_R.csv"))
all <- if (is.null(ext)) scR else rbind(scR, ext)
all$model_chem <- paste(all$model, all$chemistry)
write.csv(all, file.path(OUT, "scores_all.csv"), row.names = FALSE)
cat("wrote scores_all.csv:", nrow(all), "rows;",
    "methods:", paste(sort(unique(all$method)), collapse = ", "), "\n")
cat("external scored:", if (is.null(ext)) 0 else nrow(ext), "rows\n")
