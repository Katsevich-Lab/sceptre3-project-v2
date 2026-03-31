#!/usr/bin/env Rscript
# Usage:
#   Rscript run-all-datasets.R [datasets_config.csv] [batch_size]
#
# Examples:
#   Rscript run-all-datasets.R
#   Rscript run-all-datasets.R datasets_config.csv 3
#
# Run from the replogle-rd7/ directory.
# Submits dataset jobs to SGE in batches.
# Within each dataset:
#   run-set-analysis.sh -> run-sceptre-pipeline.sh
# Each new batch waits for all pipeline jobs from the previous batch.

args <- commandArgs(trailingOnly = TRUE)
csv_file   <- if (length(args) >= 1) args[1] else "datasets_config.csv"
batch_size <- if (length(args) >= 2) as.integer(args[2]) else 3L

if (is.na(batch_size) || batch_size < 1) {
  stop("batch_size must be a positive integer")
}

config <- read.csv(csv_file, stringsAsFactors = FALSE)

if (!"dataset" %in% names(config)) {
  stop("Expected a column named 'dataset' in ", csv_file)
}

submit_qsub <- function(args) {
  out <- system2("qsub", args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")

  if (!is.null(status) && status != 0) {
    stop("qsub failed:\n", paste(out, collapse = "\n"))
  }

  out_text <- paste(out, collapse = " ")
  job_id <- sub(".*Your job ([0-9]+).*", "\\1", out_text)

  if (!grepl("^[0-9]+$", job_id)) {
    stop("Could not parse job ID from qsub output:\n", paste(out, collapse = "\n"))
  }

  list(job_id = job_id, raw_output = out)
}

prev_batch_pipe_ids <- character(0)

for (batch_start in seq(1, nrow(config), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, nrow(config))
  current_batch_pipe_ids <- character(0)

  cat(sprintf(
    "\n=== Submitting batch %d-%d of %d ===\n",
    batch_start, batch_end, nrow(config)
  ))

  for (i in batch_start:batch_end) {
    assoc_dataset_name <- config$dataset[i]
    data_name <- strsplit(assoc_dataset_name, "_", fixed = TRUE)[[1]][1]

    cat(sprintf("  [%d/%d] %s\n", i, nrow(config), assoc_dataset_name))

    env_vars <- c(
      paste0("DATA_NAME=", data_name),
      paste0("ASSOC_DATASET_NAME=", assoc_dataset_name)
    )

    # Submit set-analysis job.
    # If this is not the first batch, hold on all pipeline jobs from the previous batch.
    set_args <- c()
    if (length(prev_batch_pipe_ids) > 0) {
      set_args <- c(set_args, "-hold_jid", paste(prev_batch_pipe_ids, collapse = ","))
    }
    set_args <- c(
      set_args,
      "-v", env_vars[1],
      "-v", env_vars[2],
      "run-set-analysis.sh"
    )

    set_res <- submit_qsub(set_args)
    cat("    qsub output:", paste(set_res$raw_output, collapse = " "), "\n")
    cat("    set-analysis job:", set_res$job_id, "\n")

    # Submit pipeline job, holding on its own set-analysis job
    pipe_args <- c(
      "-hold_jid", set_res$job_id,
      "-v", env_vars[1],
      "-v", env_vars[2],
      "run-sceptre-pipeline.sh"
    )

    pipe_res <- submit_qsub(pipe_args)
    cat("    qsub output:", paste(pipe_res$raw_output, collapse = " "), "\n")
    cat("    pipeline job:", pipe_res$job_id, "\n")

    current_batch_pipe_ids <- c(current_batch_pipe_ids, pipe_res$job_id)
  }

  prev_batch_pipe_ids <- current_batch_pipe_ids
}

cat("\nAll jobs submitted. Monitor with: qstat\n")
