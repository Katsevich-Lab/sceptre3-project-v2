#!/usr/bin/env Rscript
# Insert Otsu(log1p) and smoothed-valley into the "Metrics vs mu_pert" line plots
# from sceptre-nb/nb-bench_v3.Rmd. Uses that framework's OWN helpers
# (make_guide_metric_tbl / summarize_metric_tbl / make_avg_metrics_plot) and its
# load_run logic, so baseline numbers match the canonical plots; the two new
# methods are added to `all_assns` and scored identically.

suppressPackageStartupMessages({library(tidyverse); library(Matrix)})
HERE  <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA    <- normalizePath(file.path(HERE, ".."))
BENCH <- path.expand("~/data/projects/sceptre3/benchmarking")
source(file.path(GA, "sceptre-nb", "nb-bench_v3-helpers.R"))
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

# ---- load_run (faithful copy of nb-bench_v3.Rmd:35-126, base-R paths) -------
load_run <- function(run_name, dataset_name) {
  in_p  <- file.path(BENCH, "guide_assignment/input_data", dataset_name)
  out_p <- file.path(BENCH, "guide_assignment/outputs", run_name)
  true_perts  <- readRDS(file.path(in_p, "true_pert_matrix.rds"))
  grna_matrix <- readRDS(file.path(in_p, "sceptre/grna_matrix.rds"))
  guide_names <- rownames(grna_matrix)

  fnames <- dir(out_p)
  amf <- fnames[grepl("^assignment_matrix", fnames) & grepl(dataset_name, fnames, fixed = TRUE)]
  all_assns <- list()
  for (f in amf) {
    nm <- gsub("^assignment_matrix_", "", f) |>
      gsub(pattern = "^script_", replacement = "", fixed = TRUE) |>   # fixed=TRUE: "^script_" literal, so kept (matches Rmd)
      gsub(pattern = paste0("_", dataset_name, ".rds"), replacement = "", fixed = TRUE)
    mat <- readRDS(file.path(out_p, f))
    all_assns[[nm]] <- setNames(lapply(guide_names, function(g) as.integer(which(mat[g, ]))), guide_names)
  }
  # crispat (cell names -> indices)
  cf <- fnames[grepl("crispat", fnames) & grepl(dataset_name, fnames)]
  if (length(cf)) {
    raw <- read.csv(file.path(out_p, cf[1]))
    idx <- setNames(seq_len(ncol(grna_matrix)), colnames(grna_matrix))
    raw$cell_idx <- idx[as.character(raw$cell_id)]
    ul <- split(raw$cell_idx, raw$grna_id)
    cr <- setNames(vector("list", length(guide_names)), guide_names)
    cr[names(ul)] <- ul
    all_assns$crispat <- cr
  }
  true_assns <- setNames(lapply(guide_names, function(g)
    as.integer(which(true_perts[g, ] == 1 & grna_matrix[g, ] > 0))), guide_names)
  list(grna_matrix = grna_matrix, all_assns = all_assns, true_assns = true_assns)
}

# ---- add the two nonparametric threshold methods to all_assns ---------------
add_threshold_methods <- function(sim) {
  gm <- as(sim$grna_matrix, "RsparseMatrix")
  gn <- rownames(gm)
  thr <- function(fn) setNames(lapply(gn, function(g) {
    counts <- as.numeric(gm[g, ]); t <- fn(counts)$t
    if (is.finite(t)) as.integer(which(counts >= t)) else integer(0)
  }), gn)
  sim$all_assns[["otsu_log1p"]]      <- thr(otsu_threshold_log1p)
  sim$all_assns[["smoothed_valley"]] <- thr(smoothed_valley_threshold)
  sim
}

# ---- generic panel runner (one mu_pert-sweep plot, mirrors the Rmd chunks) ---
run_panel <- function(run_name, dataset_name, P_to_mean, x_breaks, keep,
                      title, out_png, rename = NULL, height = 4.5) {
  sim <- load_run(run_name, dataset_name) |> add_threshold_methods()
  guide_meta <- tibble(grna_id = rownames(sim$grna_matrix),
                       np = extract_part(grna_id, 2), p = extract_part(grna_id, 3)) |>
    left_join(P_to_mean, by = "p")
  gmt <- make_guide_metric_tbl(sim) |> left_join(guide_meta, by = "grna_id")
  avg <- summarize_metric_tbl(gmt, group_vars = c("method_name", "np", "p", "mu_pert"))
  plot_df <- avg |> mutate(ci = 2 * sd / sqrt(n_nonmissing),
                           ymin = pmax(0, mean - ci), ymax = pmin(1, mean + ci))
  if (!is.null(rename))                      # e.g. disp run has no sceptre matrix
    plot_df <- plot_df |> mutate(method_name = ifelse(method_name %in% names(rename),
                                                      rename[method_name], method_name))
  keep <- intersect(keep, unique(plot_df$method_name))
  p <- plot_df |> filter(method_name %in% keep) |>
    make_avg_metrics_plot(x_breaks = x_breaks) + ggtitle(title)
  ggsave(file.path(HERE, "results", out_png), p, width = 9, height = height, dpi = 120)
  cat("\n=====", title, "=====\n")
  print(plot_df |> filter(method_name %in% keep) |>
          select(method_name, np, p, mu_pert, metric, mean) |>
          pivot_wider(names_from = metric, values_from = mean) |>
          arrange(np, p, method_name), n = 200)
  cat("Wrote results/", out_png, "\n", sep = "")
}

# Panel 1: old simulations (NPpois background, one np group) ------------------
run_panel("pois0nb_sims_repeat_old", "sims_sum_repeat_old",
          tibble(p = c("Plarge", "Pmed", "Psmall"), mu_pert = c(970, 970/4, 970/8)),
          x_breaks = c(120, 240, 480, 970),
          keep = c("crispat", "sceptre", "script_threshglmpois1000_pois0nb",
                   "script_threshglmpois_pois0nb", "otsu_log1p", "smoothed_valley"),
          title = "Metrics for old simulations (+ nonparametric thresholds)",
          out_png = "framework_musweep_old_sims.png")

# Panel 2: sims_sum_2np_3p (two np groups, moderate effects) ------------------
run_panel("pois0nb_sims_sum_2np_3p", "sims_sum_2np_3p",
          tibble(p = c("P1", "P2", "P3"), mu_pert = c(500, 750, 1000)),
          x_breaks = c(500, 750, 1000),
          keep = c("crispat", "sceptre", "script_threshglmpois_pois0nb",
                   "otsu_log1p", "smoothed_valley"),
          title = "Metrics for sims_sum_2np_3p (+ nonparametric thresholds)",
          out_png = "framework_musweep_2np_3p.png", height = 6)

# Panel 3: sims_sum_1np_3p_disp (large, overdispersed; glmpois_poisbug=sceptre)
run_panel("pois0nb_sims_sum_1np_3p_disp", "sims_sum_1np_3p_disp",
          tibble(p = c("Phighdisp", "Pmeddisp", "Plowdisp"),
                 mu_pert = c(10001000, 1001000, 101000)),
          x_breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
          keep = c("crispat", "sceptre", "script_threshglmpois1000_pois0nb",
                   "otsu_log1p", "smoothed_valley"),
          rename = c(script_glmpois_poisbug = "sceptre"),
          title = "Metrics for sims_sum_1np_3p_disp (+ nonparametric thresholds)",
          out_png = "framework_musweep_1np_3p_disp.png")
