#!/usr/bin/env Rscript
# Usage: Rscript run-all-datasets.R [datasets_config.csv]
# Run from the replogle-rd7/ directory.
# Submits all dataset jobs to SGE with dependency chaining.
# If a job fails, dependent jobs are held in Eqw state (not run).

args <- commandArgs(trailingOnly = TRUE)
csv_file <- if (length(args) >= 1) args[1] else "datasets_config.csv"

config <- read.csv(csv_file, stringsAsFactors = FALSE)

prev_job_id <- ""

for (i in seq_len(nrow(config))) {
  assoc_dataset_name <- config$dataset[i]
  data_name <- strsplit(assoc_dataset_name, "_")[[1]][1]

  cat(sprintf("=== [%d/%d] %s ===\n", i, nrow(config), assoc_dataset_name))

  env_vars <- c(
    paste0("DATA_NAME=", data_name),
    paste0("ASSOC_DATASET_NAME=", assoc_dataset_name)
  )

  # Submit set-analysis job (hold on previous pipeline job if any)
  set_args <- c(
    if (nchar(prev_job_id) > 0) c("-hold_jid", prev_job_id),
    "-v", env_vars[1], "-v", env_vars[2],
    "run-set-analysis.sh"
  )
  set_out <- system2("qsub", args = set_args, stdout = TRUE)
  cat("  qsub output:", set_out, "\n")
  set_job_id <- sub(".*Your job ([0-9]+).*", "\\1", paste(set_out, collapse = " "))
  cat("  set-analysis job:", set_job_id, "\n")

  # Submit pipeline job, holding on set-analysis
  pipe_args <- c(
    "-hold_jid", set_job_id,
    "-v", env_vars[1], "-v", env_vars[2],
    "run-sceptre-pipeline.sh"
  )
  pipe_out <- system2("qsub", args = pipe_args, stdout = TRUE)
  cat("  qsub output:", pipe_out, "\n")
  prev_job_id <- sub(".*Your job ([0-9]+).*", "\\1", paste(pipe_out, collapse = " "))
  cat("  pipeline job:", prev_job_id, "\n")
}

cat("All jobs submitted. Monitor with: qstat\n")
