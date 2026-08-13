#!/usr/bin/env Rscript

# Score freshly generated SCEPTRE and CLEANSER assignments with the Fishash
# authors' exact Table 2 cohort and binary species-assignment metric.

suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
samples <- c("mix0hr_Cropseq", "mix72hr_Cropseq", "mix0hr_DirectCapture", "mix72hr_DirectCapture")
`%||%` <- function(x, y) if (is.null(x)) y else x

published <- rbind(
  `cleanser cs (Table 2, 0.80)` = c(.9414, .9410, .8612, .9221),
  `cleanser dc (Table 2, 0.50)` = c(.9178, .8459, .8614, .9308),
  `sceptre mixture` = c(.9491, .9472, .8953, .9408)
)
colnames(published) <- samples
published_se <- rbind(
  `cleanser cs (Table 2, 0.80)` = c(.0030, .0027, .0039, .0028),
  `cleanser dc (Table 2, 0.50)` = c(.0034, .0041, .0039, .0027),
  `sceptre mixture` = c(.0028, .0026, .0036, .0025)
)
colnames(published_se) <- samples

binary_accuracy_metrics <- function(pred, truth) {
  tp <- sum(pred & truth, na.rm = TRUE); fp <- sum(pred & !truth, na.rm = TRUE)
  tn <- sum(!pred & !truth, na.rm = TRUE); fn <- sum(!pred & truth, na.rm = TRUE)
  na_true <- sum(is.na(pred) & truth); na_false <- sum(is.na(pred) & !truth)
  total <- tp + fp + tn + fn + na_true + na_false
  total_true <- tp + fn + na_true; total_false <- tn + fp + na_false
  acc_true <- tp / total_true; acc_false <- tn / total_false
  stderr <- sqrt((total_true * acc_true * (1 - acc_true) +
                  total_false * acc_false * (1 - acc_false)) / total^2)
  c(accuracy = (tp + tn) / total, stderr = stderr, n = total,
    both_or_neither = na_true + na_false)
}

read_mtx_assignment <- function(path) as(readMM(path) > 0, "CsparseMatrix")

# CLEANSER 1.2.x writes Matrix-Market-like coordinate output without a valid
# banner. The first non-comment row contains dimensions; later rows are
# 1-based guide, cell, P(native).
read_cleanser_assignment <- function(path, cutoff) {
  x <- read.table(path, comment.char = "%")
  dims <- as.integer(x[1L, 1:2])
  x <- x[-1L, , drop = FALSE]
  keep <- x[, 3] >= cutoff
  sparseMatrix(i = x[keep, 1], j = x[keep, 2], x = TRUE, dims = dims)
}

score_assignment <- function(A, guide_type, truth, selected) {
  stopifnot(nrow(A) == length(guide_type), ncol(A) == length(truth))
  has_human <- colSums(A[guide_type == "homo_guide", , drop = FALSE]) > 0
  has_mouse <- colSums(A[guide_type == "mus_guide", , drop = FALSE]) > 0
  assigned <- has_human
  assigned[has_human & has_mouse] <- NA
  assigned[!has_human & !has_mouse] <- NA
  binary_accuracy_metrics(assigned[selected], truth[selected])
}

specs <- list(
  `cleanser cs (Table 2, 0.80)` = list(kind = "cleanser", mode = "cs", cutoff = .80),
  `cleanser dc (Table 2, 0.50)` = list(kind = "cleanser", mode = "dc", cutoff = .50),
  `sceptre mixture` = list(kind = "mtx"),
  `cleanser cs (public-code 0.95)` = list(kind = "cleanser", mode = "cs", cutoff = .95),
  `cleanser dc (public-code 0.95)` = list(kind = "cleanser", mode = "dc", cutoff = .95)
)

out <- list()
for (sample in samples) {
  meta <- read.csv(file.path(HERE, paste0(sample, "_meta.csv")))
  guides <- read.csv(file.path(HERE, paste0(sample, "_guides.csv")))
  sum_gex <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_human <- meta$homo_sum_gex / sum_gex
  selected <- meta$mito_sum / sum_gex < .15 & meta$features_gex <= 6000 &
    sum_gex <= 20000 & meta$features_gex >= 1500 & sum_gex >= 3500 &
    (frac_human < .1 | frac_human > .9)
  truth <- frac_human >= .1 & sum_gex > 100

  for (method in names(specs)) {
    spec <- specs[[method]]
    if (spec$kind == "cleanser") {
      path <- file.path(HERE, paste0(sample, "_cleanser_", spec$mode, "_posterior.mtx"))
      log_path <- file.path(HERE, paste0(sample, "_cleanser_", spec$mode, ".log"))
      log_lines <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
      complete <- any(grepl("Random seed:", log_lines, fixed = TRUE))
      if (!file.exists(path) || file.info(path)$size == 0L || !complete) next
      A <- read_cleanser_assignment(path, spec$cutoff)
      nonfatal <- sum(grepl("Non-fatal error during sampling", log_lines, fixed = TRUE))
      convergence <- sum(grepl("Some chains may have failed to converge", log_lines, fixed = TRUE))
      divergent <- sum(grepl("divergent transitions", log_lines, fixed = TRUE))
    } else {
      path <- file.path(HERE, paste0(sample, "_sceptre_mixture.mtx"))
      if (!file.exists(path) || file.info(path)$size == 0L) next
      A <- read_mtx_assignment(path)
      nonfatal <- convergence <- divergent <- NA_integer_
    }
    result <- score_assignment(A, guides$guide_type, truth, selected)
    expected <- if (method %in% rownames(published)) published[method, sample] else NA_real_
    expected_se <- if (method %in% rownames(published_se)) published_se[method, sample] else NA_real_
    accuracy <- unname(result["accuracy"])
    stderr <- unname(result["stderr"])
    out[[length(out) + 1L]] <- data.frame(
      sample = sample, method = method, cutoff = spec$cutoff %||% NA_real_,
      accuracy = accuracy, stderr = stderr,
      n_qc = unname(result["n"]), both_or_neither = unname(result["both_or_neither"]),
      nonfatal_sampler_warnings = nonfatal, convergence_warnings = convergence,
      divergent_chain_warnings = divergent,
      published_accuracy = expected, published_stderr = expected_se,
      accuracy_delta = if (is.na(expected)) NA_real_ else accuracy - expected,
      stderr_delta = if (is.na(expected_se)) NA_real_ else stderr - expected_se,
      accuracy_match_4dp = if (is.na(expected)) NA else round(accuracy, 4L) == expected,
      stderr_match_4dp = if (is.na(expected_se)) NA else round(stderr, 4L) == expected_se,
      accuracy_within_0_001 = if (is.na(expected)) NA else abs(accuracy - expected) <= .001
    )
  }
}

result <- do.call(rbind, out)
write.csv(result, file.path(HERE, "extended_table2_reproduction.csv"), row.names = FALSE)
print(result, row.names = FALSE, digits = 5)
