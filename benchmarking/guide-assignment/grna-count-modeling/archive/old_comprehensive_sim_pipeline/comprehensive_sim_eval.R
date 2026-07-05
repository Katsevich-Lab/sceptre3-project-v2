#!/usr/bin/env Rscript
# Decisive test: does the principled ambient-proportion test perform well on
# realistic MANY-guide simulations (Gasperini- and Replogle-based), alongside
# the established parametric baselines and the simple Otsu/valley thresholds?
# Per-guide precision/recall/Jaccard vs ground truth; overall + by mu_pert level.

suppressPackageStartupMessages({library(tidyverse); library(Matrix)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
BENCH<- path.expand("~/data/projects/sceptre3/benchmarking")
source(file.path(GA, "sceptre-nb", "nb-bench_v3-helpers.R"))           # set_metrics, make_guide_metric_tbl, summarize
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

load_run <- function(run_name, dataset_name) {
  in_p  <- file.path(BENCH, "guide_assignment/input_data", dataset_name)
  grna  <- readRDS(file.path(in_p, "sceptre/grna_matrix.rds")); gn <- rownames(grna)
  true_perts <- readRDS(file.path(in_p, "true_pert_matrix.rds"))
  all_assns <- list()
  if (!is.null(run_name)) {
    out_p <- file.path(BENCH, "guide_assignment/outputs", run_name); fn <- dir(out_p)
    for (f in fn[grepl("^assignment_matrix", fn) & grepl(dataset_name, fn, fixed = TRUE)]) {
      nm <- gsub("^assignment_matrix_", "", f) |>
        gsub(pattern = paste0("_", dataset_name, ".rds"), replacement = "", fixed = TRUE)
      mat <- readRDS(file.path(out_p, f))
      all_assns[[nm]] <- setNames(lapply(gn, function(g) as.integer(which(mat[g, ]))), gn)
    }
    cf <- fn[grepl("crispat", fn) & grepl(dataset_name, fn)]
    if (length(cf)) { raw <- read.csv(file.path(out_p, cf[1]))
      idx <- setNames(seq_len(ncol(grna)), colnames(grna)); raw$ci <- idx[as.character(raw$cell_id)]
      ul <- split(raw$ci, raw$grna_id); cr <- setNames(vector("list", length(gn)), gn)
      cr[names(ul)] <- ul; all_assns$crispat <- cr }
  }
  true_assns <- setNames(lapply(gn, function(g)
    as.integer(which(true_perts[g, ] == 1 & grna[g, ] > 0))), gn)
  list(grna_matrix = grna, all_assns = all_assns, true_assns = true_assns)
}

mat_to_list <- function(A, gn) setNames(lapply(gn, function(g) as.integer(which(A[g, ]))), gn)
add_methods <- function(sim) {
  gn <- rownames(sim$grna_matrix)
  sim$all_assns[["otsu"]]   <- mat_to_list(assign_by_threshold(sim$grna_matrix, otsu_threshold_log1p)$assignment_matrix, gn)
  sim$all_assns[["valley"]] <- mat_to_list(assign_by_threshold(sim$grna_matrix, smoothed_valley_threshold)$assignment_matrix, gn)
  sim$all_assns[["ambient.05"]] <- mat_to_list(ambient_test_assign(sim$grna_matrix, q = 0.05)$assignment_matrix, gn)
  sim$all_assns[["ambient.01"]] <- mat_to_list(ambient_test_assign(sim$grna_matrix, q = 0.01)$assignment_matrix, gn)
  sim
}

sim_cfg <- list(
  list(label = "2np_3p (Repl, clean)",   run = "pois0nb_sims_sum_2np_3p",      ds = "sims_sum_2np_3p"),
  list(label = "repeat_old (Repl, hard)",run = "pois0nb_sims_repeat_old",      ds = "sims_sum_repeat_old"),
  list(label = "1np_3p_disp (overdisp)", run = "pois0nb_sims_sum_1np_3p_disp", ds = "sims_sum_1np_3p_disp"),
  list(label = "gasperini_calibrated",   run = NULL,                           ds = "sims_gasperini_calibrated"))

summ <- list()
for (cc in sim_cfg) {
  sim <- add_methods(load_run(cc$run, cc$ds))
  gt <- sim$true_assns
  for (mn in names(sim$all_assns)) {
    preds <- sim$all_assns[[mn]]
    m <- do.call(rbind, lapply(names(gt), function(g) set_metrics(preds[[g]], gt[[g]])))
    summ[[length(summ)+1]] <- data.frame(sim = cc$label, method = mn,
      precision = round(mean(m$precision, na.rm=TRUE), 3),
      recall    = round(mean(m$recall,    na.rm=TRUE), 3),
      jaccard   = round(mean(m$jaccard,   na.rm=TRUE), 3))
  }
  cat("done:", cc$label, "\n")
}
df <- do.call(rbind, summ)
cat("\n===== per-guide metrics across simulations (sorted by Jaccard within sim) =====\n")
for (s in unique(df$sim)) { cat("\n--", s, "--\n"); d <- df[df$sim == s, ]; print(d[order(-d$jaccard), -1], row.names = FALSE) }
write.csv(df, file.path(HERE, "results", "comprehensive_sim_eval.csv"), row.names = FALSE)

# plot: Jaccard by method, faceted by sim; nonparametric/ambient highlighted
df$kind <- ifelse(grepl("ambient", df$method), "ambient test",
            ifelse(df$method %in% c("otsu","valley"), "nonparam threshold", "baseline"))
p <- ggplot(df, aes(reorder(method, jaccard), jaccard, fill = kind)) +
  geom_col() + coord_flip() + facet_wrap(~sim, scales = "free_y") +
  scale_fill_manual(values = c(`ambient test`="#2c7fb8", `nonparam threshold`="#7fbf7b", baseline="grey75")) +
  labs(title = "Guide-assignment Jaccard vs ground truth across realistic sims",
       x = NULL, y = "mean per-guide Jaccard", fill = NULL) + theme_bw(base_size = 9) +
  theme(legend.position = "top")
ggsave(file.path(HERE, "results", "comprehensive_sim_eval.png"), p, width = 11, height = 6.5, dpi = 120)
cat("\nWrote results/comprehensive_sim_eval.{csv,png}\n")
