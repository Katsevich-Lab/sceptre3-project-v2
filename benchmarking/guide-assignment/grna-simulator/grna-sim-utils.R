#!/usr/bin/env Rscript
# Shared helpers for gRNA simulation work.
#
# Sourced by the simulation scripts in this directory and by
# grna-count-modeling/scripts/{sim_validate_hist,sim_poc_paired}.R.
#
# Everything above the "GuideBender scaling datasets" section is
# simulator-agnostic. That section is specific to fishash::simulate_guidebender2
# and is used by simulate-scaling-datasets.R. Packages are called with `::`
# throughout, so sourcing this file attaches nothing.


#' Gini coefficient of a non-negative vector. 0 = every element equal.
gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x); s <- sum(x)
  if (n == 0L || s == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * s) - (n + 1) / n
}

cv <- function(x) {
  m <- mean(x)
  if (!is.finite(m) || m == 0) NA_real_ else stats::sd(x) / m
}


#' Statistics of a gRNA count matrix.
#'
#' Computed the same way for real and simulated matrices so the two can be put
#' side by side and compared row for row. Nothing here assumes a generative
#' model: every number is a property of the matrix it was given.
#'
#' Entries with count >= `threshold` are called PERTURBED and the remaining
#' nonzeros AMBIENT. On real data that split is a proxy, and it errs in both
#' directions: a genuine integration carrying few UMIs is counted as ambient,
#' and an ambient entry that happens to clear the threshold is counted as
#' signal. On simulated data the ground truth is known, so pass `is_pert` to use
#' the exact split -- running both on the same simulated matrix measures how far
#' the threshold proxy is off, which is the error the real-data numbers inherit.
#'
#' @param counts     guides x cells.
#' @param threshold  count at or above which an entry is called perturbed.
#' @param is_pert    optional logical matrix, same dims as `counts`; when given,
#'                   it replaces the threshold split.
#' @param label      row label in the returned data.frame.
#' @return a one-row data.frame.
guide_matrix_stats <- function(counts, threshold, is_pert = NULL,
                               label = NA_character_) {
  ct <- as(as(counts, "CsparseMatrix"), "TsparseMatrix")
  G <- nrow(ct); N <- ncol(ct)
  gi <- ct@i + 1L; cj <- ct@j + 1L; x <- as.numeric(ct@x)

  pert <- if (is.null(is_pert)) x >= threshold
          else as.logical(is_pert[cbind(gi, cj)])
  pert[is.na(pert)] <- FALSE
  amb <- !pert
  ones <- rep(1, length(x))

  agg <- function(sel, val, margin) {
    if (!any(sel)) return(numeric(if (margin == "col") N else G))
    m <- Matrix::sparseMatrix(i = gi[sel], j = cj[sel], x = val[sel], dims = c(G, N))
    if (margin == "col") Matrix::colSums(m) else Matrix::rowSums(m)
  }

  cell_total <- Matrix::colSums(ct)
  cell_sig   <- agg(pert, x,    "col")
  cell_amb   <- agg(amb,  x,    "col")
  cell_npert <- agg(pert, ones, "col")
  cell_namb  <- agg(amb,  ones, "col")
  guide_pert <- agg(pert, ones, "row")
  guide_amb  <- agg(amb,  x,    "row")

  both <- cell_sig > 0 & cell_amb > 0

  data.frame(
    label        = label,
    n_guides     = G,
    n_cells      = N,
    threshold    = threshold,

    # size and sparsity
    nnz          = length(x),
    zero_frac    = 1 - length(x) / (as.numeric(G) * N),
    nnz_per_cell = length(x) / N,

    # perturbation
    moi            = mean(cell_npert),
    pert_per_guide = mean(guide_pert),
    pert_rate      = mean(guide_pert) / N,
    frac_cells_unpert = mean(cell_npert == 0),

    # depth, split by the perturbed/ambient call
    umis_per_cell_median = stats::median(cell_total),
    umis_per_cell_mean   = mean(cell_total),
    sig_umis_per_cell    = mean(cell_sig),
    amb_umis_per_cell    = mean(cell_amb),
    snr_observed         = mean(cell_sig) / mean(cell_amb),

    # what an ambient nonzero looks like
    amb_nnz_per_cell      = mean(cell_namb),
    amb_count_mean        = if (any(amb)) mean(x[amb]) else NA_real_,
    amb_frac_eq1          = if (any(amb)) mean(x[amb] == 1) else NA_real_,
    amb_count_var_o_mean  = if (any(amb)) stats::var(x[amb]) / mean(x[amb]) else NA_real_,

    # what a perturbed nonzero looks like
    pert_count_mean       = if (any(pert)) mean(x[pert]) else NA_real_,
    pert_count_var_o_mean = if (any(pert)) stats::var(x[pert]) / mean(x[pert]) else NA_real_,

    # guide-level spread. Under a symmetric Dirichlet(1) over G guides each
    # weight is Beta(1, G-1), whose CV is sqrt((G-1)/(G+1)) -> 1, and whose
    # Gini -> 0.5 (the Exp(1) limit). Those are properties of that distribution,
    # to be compared against the measured values here.
    guide_amb_cv    = cv(guide_amb),
    guide_amb_gini  = gini(guide_amb),
    guide_pert_cv   = cv(guide_pert),
    guide_pert_gini = gini(guide_pert),

    # per-cell spread and the coupling between the two components
    cell_sig_cv      = cv(cell_sig),
    cell_amb_cv      = cv(cell_amb),
    cor_log_sig_amb  = if (sum(both) > 2)
      stats::cor(log(cell_sig[both]), log(cell_amb[both])) else NA_real_,

    stringsAsFactors = FALSE
  )
}


#' Side-by-side UMI histograms, one real guide against one simulated guide.
#'
#' Shared log1p bins across both panels; simulated bars are stacked
#' perturbed/non-perturbed. The y-axis is log10, so bars are drawn from 1 rather
#' than 0 -- bar tops are still the true bin counts.
#'
#' @param umis_real  UMI counts for a real guide, across cells.
#' @param umis_sim   UMI counts for a simulated guide, across cells.
#' @param is_pert    Logical, same length as `umis_sim`: ground-truth perturbation.
#' @param bins       Number of shared bins.
#' @param title      Plot title.
#' @param x_breaks   UMI values to label on the (log1p) x-axis.
plot_umi_histogram_real_vs_sim <- function(
    umis_real,
    umis_sim,
    is_pert,
    bins = 80,
    title = NULL,
    x_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                 2000, 5000, 10000, 20000, 50000, 100000)
) {
  stopifnot(all(umis_real >= 0, na.rm = TRUE))
  stopifnot(all(umis_sim >= 0, na.rm = TRUE))
  stopifnot(length(umis_sim) == length(is_pert))

  umis_real <- umis_real[is.finite(umis_real)]

  keep_sim <- is.finite(umis_sim) & !is.na(is_pert)
  umis_sim <- umis_sim[keep_sim]
  is_pert <- is_pert[keep_sim]

  x_real <- log1p(umis_real)
  x_sim <- log1p(umis_sim)
  x_all <- c(x_real, x_sim)

  # Shared x-axis bins across real and simulated data.
  h_all <- hist(x_all, breaks = bins, plot = FALSE)
  breaks <- h_all$breaks

  make_bin_df <- function(x, group, panel) {
    bin_id <- cut(
      x,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    )

    df <- data.frame(
      bin_id = bin_id,
      group = group
    )

    df <- df[!is.na(df$bin_id), , drop = FALSE]

    tab <- as.data.frame(table(df$bin_id, df$group))
    names(tab) <- c("bin_id", "group", "count")

    tab$bin_id <- as.integer(as.character(tab$bin_id))
    tab$count <- as.integer(tab$count)
    tab <- tab[tab$count > 0, , drop = FALSE]

    tab$xmin <- breaks[tab$bin_id]
    tab$xmax <- breaks[tab$bin_id + 1L]
    tab$panel <- panel

    tab
  }

  df_real <- make_bin_df(
    x = x_real,
    group = "Real",
    panel = "Real"
  )

  df_sim <- make_bin_df(
    x = x_sim,
    group = ifelse(is_pert, "Perturbed", "Non-perturbed"),
    panel = "Simulated"
  )

  # Stack simulated bars manually so total bar height is the true bin count.
  # Real panel has only one group.
  df <- rbind(df_real, df_sim)

  df$group <- factor(
    df$group,
    levels = c("Real", "Non-perturbed", "Perturbed")
  )

  df <- df[order(df$panel, df$bin_id, df$group), , drop = FALSE]

  df$ymin <- 0
  df$ymax <- 0

  split_ids <- split(seq_len(nrow(df)), paste(df$panel, df$bin_id, sep = "___"))

  for (idx in split_ids) {
    counts <- df$count[idx]
    cum_counts <- cumsum(counts)

    df$ymin[idx] <- c(0, head(cum_counts, -1))
    df$ymax[idx] <- cum_counts
  }

  # Log y-scale cannot display 0, so bars start visually at 1.
  # The top of each stacked bar is still the true bin count.
  df$ymin_plot <- pmax(df$ymin, 1)
  df$ymax_plot <- df$ymax

  df <- df[df$ymax_plot > df$ymin_plot, , drop = FALSE]

  # Keep only breaks inside the observed range.
  max_count <- max(c(umis_real, umis_sim), na.rm = TRUE)
  x_breaks <- x_breaks[x_breaks <= max_count]
  if (!0 %in% x_breaks) {
    x_breaks <- c(0, x_breaks)
  }

  facet_scales <- if (length(umis_real) == length(umis_sim)) {
    "fixed"
  } else {
    "free_y"
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin_plot,
        ymax = ymax_plot,
        fill = group
      )
    ) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = facet_scales) +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      limits = range(breaks),
      name = "UMI count"
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      name = "Number of cells"
    ) +
    ggplot2::labs(
      title = title,
      fill = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 7,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )
}


# ============================================================================
# GuideBender scaling datasets
# ============================================================================
#
# A regime is defined by four inputs: the real dataset it mimics, the count
# threshold at which an entry is called perturbed, the n_cells of each rung, and
# the number of cells used in calibration runs. Every target is measured from the
# real matrix at run time (measure_real), so changing the threshold changes all of
# them consistently.
#
# Two knobs are calibrated rather than set, because the simulator does not take
# the measured quantity directly (calibrate):
#   snr    - held equal across rungs, so ambient UMIs per cell are equal across
#            rungs. Solved once, by simulating at the real dataset's number of
#            guides and matching its nonzeros per cell.
#   lambda - simulate_guidebender2's `moi`: the rate of the Poisson draw of
#            infections per cell, before zero truncation, the hurdle, and repeated
#            draws of the same guide. Solved at each rung so the thresholded MOI
#            of the simulated counts equals the real one.
#
# Cells are drawn independently given the guide-level latents, and chunk-level
# quantities depend only on the cells in the chunk, so the distribution of a
# per-cell statistic does not depend on n_cells. Calibration runs therefore use
# n_cal cells whatever the rung, and every evaluation reuses one seed.

# Dataset ids keep the regime name because run_cleanser.py chooses its CROP-seq
# (--cs) or direct-capture (--dc) model by matching "gasperini" / "replogle" in
# the id. The threshold is higher for replogle because its perturbed entries carry
# far more UMIs: measure-real-targets.R reports a mean of 1,405 at a perturbed
# entry for replogle-rd7 against 21 for gasperini.
SCALING_REGIMES <- list(
  gasperini = list(real_dataset = "gasperini",    threshold = 5L,
                   n_cells = c(3000L, 6000L, 12000L, 24000L, 48000L),
                   n_cal   = 4000L),
  replogle  = list(real_dataset = "replogle-rd7", threshold = 10L,
                   n_cells = c(12000L, 24000L, 48000L, 96000L, 192000L),
                   n_cal   = 10000L)
)

# The methods the datasets are built for, and the input each reads:
#   "h5ad" - <dataset>/<method>/grna_matrix.h5ad, cells x guides
#   "mtx"  - <dataset>/cleanser/grna_matrix.mtx, guides x cells
# fishash has no input directory of its own: main.nf maps it to cleanser/, since
# fishash() takes the same Matrix Market count matrix. A new method is added here,
# and write_method_inputs() stops on any input type it has no writer for.
SCALING_METHODS <- list(
  crispat  = list(input = "h5ad"),
  pertpy   = list(input = "h5ad"),
  cleanser = list(input = "mtx"),
  fishash  = list(input = "mtx")
)

# Held fixed in both regimes: the fishash package defaults, except
# frac_noise_endo (no default; set to the paper's high-expression regime,
# Section 3.2) and chunk_cells (default NULL).
SCALING_PARAMS_SHARED <- list(
  Phi_cell        = 1,
  Phi_noise       = 0,
  frac_noise_endo = 0.75,   # 25% of the noise exogenous
  rho_sum         = 10,
  eps_alpha       = 50,
  d_sigma_guide   = 0.5,    # standard deviations of the log size factors; the
  d_sigma_cell    = 0.5,    # paper states a variance of 1/2, the code takes an
  d_sigma_drop    = 0.5,    # SD and defaults it to 0.5
  endo_shape_sum  = 1,
  endo_shape_flat = 0,
  use_median      = TRUE,   # count_per_cell is a median, as it is measured
  chunk_cells     = 1000    # the exogenous profile and the median rescaling are
                            # computed per chunk, so this is held fixed for every
                            # dataset, calibration runs included
)

SCALING_CAL_TOL    <- 0.005   # calibration stops within 0.5% of each target
SCALING_CAL_PASSES <- 6L


#' guide_infection_alpha from the measured spread of perturbed cells per guide.
#'
#' Under a symmetric Dirichlet(alpha) over G guides each weight has
#' CV^2 = (G-1)/(G*alpha + 1). A count drawn around each weight adds sampling
#' variance of about 1/mean in CV^2 units, which is removed first.
alpha_from_cv <- function(cv, G, mean_count) {
  cv2 <- cv^2 - 1 / mean_count
  ((G - 1) / cv2 - 1) / G
}

#' Measure a real matrix and derive everything the simulation is matched to.
#' Needs .get_config_path() (from ~/.Rprofile).
measure_real <- function(reg) {
  fp <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data",
                  reg$real_dataset, "cleanser", "grna_matrix.mtx")
  cat("reading", fp, "\n")
  s <- guide_matrix_stats(as(Matrix::readMM(fp), "CsparseMatrix"), reg$threshold,
                          label = paste0(reg$real_dataset, "_real"))
  list(
    stats          = s,
    n_guides       = s$n_guides,
    n_over_g       = s$n_cells / s$n_guides,
    moi            = s$moi,
    nnz_per_cell   = s$nnz_per_cell,
    count_per_cell = s$umis_per_cell_median,
    hurdle_prob    = s$frac_cells_unpert,
    alpha          = alpha_from_cv(s$guide_pert_cv, s$n_guides, s$pert_per_guide),
    snr_start      = s$snr_observed
  )
}

#' The ray. n_cells is the input and n_guides is derived, because chunk_cells
#' must divide n_cells while n_guides has no such constraint.
regime_rungs <- function(reg, n_over_g, name) {
  stopifnot(reg$n_cells %% SCALING_PARAMS_SHARED$chunk_cells == 0,
            reg$n_cal   %% SCALING_PARAMS_SHARED$chunk_cells == 0)
  rungs <- data.frame(n_cells  = reg$n_cells,
                      n_guides = as.integer(round(reg$n_cells / n_over_g)))
  rungs$dataset <- sprintf("sim_%s_g%d_c%d", name, rungs$n_guides, rungs$n_cells)
  rungs
}

#' The single call into the simulator, used by calibration and generation alike.
#' `tg` is the output of measure_real().
simulate_regime <- function(tg, n_guides, n_cells, lambda, snr, seed) {
  set.seed(seed)
  do.call(fishash::simulate_guidebender2, c(
    list(n_guides = n_guides, n_cells = n_cells, moi = lambda, snr = snr,
         count_per_cell        = tg$count_per_cell,
         hurdle_prob           = tg$hurdle_prob,
         guide_infection_alpha = tg$alpha,
         return_sparse_only    = TRUE),
    SCALING_PARAMS_SHARED
  ))
}

#' Thresholded MOI and nonzeros per cell, defined as in guide_matrix_stats().
cal_stats <- function(counts, threshold) {
  counts <- as(counts, "CsparseMatrix")
  x <- counts@x[counts@x > 0]
  c(moi = sum(x >= threshold) / ncol(counts), nnz_per_cell = length(x) / ncol(counts))
}

#' Root of f(x) = target, searched on a log scale from x0. `dir` is "upX" when f
#' increases in x and "downX" when it decreases, so uniroot extends the bracket
#' the right way.
solve_log <- function(f, target, x0, dir) {
  exp(stats::uniroot(function(lx) f(exp(lx)) - target,
                     lower = log(x0) - 0.2, upper = log(x0) + 0.2,
                     extendInt = dir, tol = 1e-3)$root)
}

#' Calibrate snr (once, at the real number of guides) and lambda (per rung).
#' Returns one row per calibration: the real-G solve, then each rung.
calibrate <- function(reg, tg, rungs, seed) {
  evaluate <- function(G, lambda, snr) {
    sim <- simulate_regime(tg, G, reg$n_cal, lambda, snr, seed)
    cal_stats(SummarizedExperiment::assay(sim, "counts"), reg$threshold)
  }
  target <- c(moi = tg$moi, nnz_per_cell = tg$nnz_per_cell)
  within <- function(st, which = names(target))
    all(abs(st[which] / target[which] - 1) < SCALING_CAL_TOL)
  row <- function(role, dataset, G, lambda, snr, st)
    data.frame(role = role, dataset = dataset, n_guides = G, lambda = lambda, snr = snr,
               cal_moi = st[["moi"]], cal_nnz_per_cell = st[["nnz_per_cell"]])

  # snr and lambda together, at the real number of guides
  G0 <- tg$n_guides
  snr <- tg$snr_start
  lambda <- tg$moi / (1 - tg$hurdle_prob)
  cat(sprintf("\n[calibrate] snr at the real G = %d (%d cells per run)\n", G0, reg$n_cal))
  for (pass in seq_len(SCALING_CAL_PASSES)) {
    lambda <- solve_log(function(l) evaluate(G0, l, snr)[["moi"]],
                        target[["moi"]], lambda, "upX")
    snr    <- solve_log(function(s) evaluate(G0, lambda, s)[["nnz_per_cell"]],
                        target[["nnz_per_cell"]], snr, "downX")
    st <- evaluate(G0, lambda, snr)
    cat(sprintf("  pass %d: lambda %.4f  snr %.3f  ->  moi %.3f (target %.3f)  nnz/cell %.3f (target %.3f)\n",
                pass, lambda, snr, st[["moi"]], target[["moi"]],
                st[["nnz_per_cell"]], target[["nnz_per_cell"]]))
    if (within(st)) break
  }
  if (!within(st))
    stop("snr calibration did not reach ", SCALING_CAL_TOL * 100, "% in ",
         SCALING_CAL_PASSES, " passes")
  cal <- row("real_G", NA_character_, G0, lambda, snr, st)

  # lambda at each rung, snr held at the value above
  cat(sprintf("\n[calibrate] lambda per rung (snr fixed at %.3f)\n", snr))
  lam <- lambda
  for (i in seq_len(nrow(rungs))) {
    G <- rungs$n_guides[i]
    lam <- solve_log(function(l) evaluate(G, l, snr)[["moi"]], target[["moi"]], lam, "upX")
    st <- evaluate(G, lam, snr)
    if (!within(st, "moi")) stop("lambda calibration missed the MOI target at G = ", G)
    cat(sprintf("  G = %5d: lambda %.4f  ->  moi %.3f  nnz/cell %.3f\n",
                G, lam, st[["moi"]], st[["nnz_per_cell"]]))
    cal <- rbind(cal, row("rung", rungs$dataset[i], G, lam, snr, st))
  }
  cal
}

#' Write each method's input for one dataset, according to SCALING_METHODS.
#' Needs write_cleanser_method() and write_h5ad_methods() from
#' data-preprocessing/lib_make_guide_data.R.
write_method_inputs <- function(counts, ds_dir, methods = names(SCALING_METHODS)) {
  unknown <- setdiff(methods, names(SCALING_METHODS))
  if (length(unknown)) stop("not in SCALING_METHODS: ", paste(unknown, collapse = ", "))
  inputs <- vapply(SCALING_METHODS[methods], `[[`, character(1), "input")
  writers <- list(
    h5ad = function(m) write_h5ad_methods(counts, ds_dir, m),
    mtx  = function(m) write_cleanser_method(counts, ds_dir)
  )
  missing <- setdiff(unique(inputs), names(writers))
  if (length(missing)) stop("no writer for input type: ", paste(missing, collapse = ", "))
  for (type in unique(inputs)) writers[[type]](names(inputs)[inputs == type])
}

SCALING_MANIFEST_COLS <- c(
  "label", "n_guides", "n_cells", "nnz", "zero_frac", "nnz_per_cell", "moi",
  "pert_per_guide", "frac_cells_unpert", "umis_per_cell_median", "umis_per_cell_mean",
  "snr_observed", "amb_count_mean", "amb_frac_eq1", "guide_pert_cv", "guide_amb_cv"
)

#' Simulate one rung, write every method's input for it, and return its manifest
#' row.
generate_rung <- function(name, reg, tg, dataset, n_guides, n_cells, lambda, snr,
                          seed, out_root, methods = names(SCALING_METHODS)) {
  t0 <- Sys.time()
  sim <- simulate_regime(tg, n_guides, n_cells, lambda, snr, seed)
  gen_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  counts <- SummarizedExperiment::assay(sim, "counts")
  truth  <- SummarizedExperiment::assay(sim, "ground_truth")
  colnames(counts) <- colnames(truth) <- paste0("CELL_", seq_len(n_cells))   # pipeline convention
  rownames(counts) <- rownames(truth) <- paste0("grna_", seq_len(n_guides))

  ds_dir <- file.path(out_root, dataset)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  write_method_inputs(counts, ds_dir, methods)
  saveRDS(truth, file.path(ds_dir, "true_pert_matrix.rds"))
  saveRDS(c(SCALING_PARAMS_SHARED,
            list(regime = name, real_dataset = reg$real_dataset, threshold = reg$threshold,
                 n_guides = n_guides, n_cells = n_cells, moi = lambda, snr = snr,
                 count_per_cell = tg$count_per_cell, hurdle_prob = tg$hurdle_prob,
                 guide_infection_alpha = tg$alpha, seed = seed, methods = methods,
                 targets = tg[setdiff(names(tg), "stats")])),
          file.path(ds_dir, "sim_params.rds"))

  s_thr   <- guide_matrix_stats(counts, reg$threshold, label = dataset)
  s_truth <- guide_matrix_stats(counts, reg$threshold, is_pert = truth)
  out <- s_thr[, SCALING_MANIFEST_COLS]
  out$moi_truth            <- s_truth$moi
  out$pert_per_guide_truth <- s_truth$pert_per_guide
  out$lambda  <- lambda
  out$snr     <- snr
  out$gen_min <- round(gen_min, 2)
  out
}
