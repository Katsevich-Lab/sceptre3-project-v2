library(dplyr)
library(tidyr)
library(purrr)



## this block is for visualizing UMI counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_umi_histogram <- function(
    counts,
    bins = 80,
    title = NULL,
    x_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
) {
  stopifnot(is.numeric(counts))
  stopifnot(all(counts >= 0, na.rm = TRUE))
  
  counts <- counts[is.finite(counts)]
  x <- log1p(counts)
  
  # Keep only breaks inside the observed range.
  x_breaks <- x_breaks[x_breaks <= max(counts, na.rm = TRUE)]
  if (!0 %in% x_breaks) {
    x_breaks <- c(0, x_breaks)
  }
  
  # Compute histogram manually so that log-y bars start at 1, not 0.
  h <- hist(x, breaks = bins, plot = FALSE)
  
  df <- data.frame(
    xmin = head(h$breaks, -1),
    xmax = tail(h$breaks, -1),
    count = h$counts
  )
  
  df <- df[df$count > 0, , drop = FALSE]
  
  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 1,
        ymax = count
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      name = "UMI count"
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      name = "Number of cells"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_classic()
}



plot_umi_histogram_real_vs_sim <- function(
    umis_real,
    umis_sim,
    is_pert,
    bins = 80,
    title = NULL,
    x_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                 2000, 5000, 10000, 20000, 50000, 100000)
    # x_breaks = c(0,1,2,3,4,5,6,7,8,9,10 * 2^(0:14))
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


plot_umi_histogram_real_vs_sim_list <- function(
    umis_real,
    sim_list,
    bins = 40,
    title = NULL,
    x_breaks = c(
      0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
      2000, 5000, 10000, 20000, 50000, 100000
    )
) {
  stopifnot(all(umis_real >= 0, na.rm = TRUE))
  stopifnot(is.list(sim_list))
  
  if (is.null(names(sim_list)) || any(names(sim_list) == "")) {
    names(sim_list) <- paste0("Simulated ", seq_along(sim_list))
  }
  
  umis_real <- umis_real[is.finite(umis_real)]
  
  sim_list <- lapply(sim_list, function(x) {
    stopifnot(is.list(x))
    stopifnot(all(c("umis_sim", "is_pert") %in% names(x)))
    stopifnot(length(x$umis_sim) == length(x$is_pert))
    stopifnot(all(x$umis_sim >= 0, na.rm = TRUE))
    
    keep <- is.finite(x$umis_sim) & !is.na(x$is_pert)
    
    list(
      umis_sim = x$umis_sim[keep],
      is_pert = x$is_pert[keep]
    )
  })
  
  x_real <- log1p(umis_real)
  
  x_sim_all <- unlist(
    lapply(sim_list, function(x) log1p(x$umis_sim)),
    use.names = FALSE
  )
  
  x_all <- c(x_real, x_sim_all)
  
  # Shared x-axis bins across real and all simulated datasets.
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
  
  df_sims <- lapply(names(sim_list), function(sim_name) {
    curr <- sim_list[[sim_name]]
    
    make_bin_df(
      x = log1p(curr$umis_sim),
      group = ifelse(curr$is_pert, "Perturbed", "Non-perturbed"),
      panel = sim_name
    )
  })
  
  df <- do.call(rbind, c(list(df_real), df_sims))
  
  df$group <- factor(
    df$group,
    levels = c("Real", "Non-perturbed", "Perturbed")
  )
  
  df$panel <- factor(
    df$panel,
    levels = c("Real", names(sim_list))
  )
  
  df <- df[order(df$panel, df$bin_id, df$group), , drop = FALSE]
  
  # Manually stack bars within each panel/bin.
  df$ymin <- 0
  df$ymax <- 0
  
  split_ids <- split(
    seq_len(nrow(df)),
    paste(df$panel, df$bin_id, sep = "___")
  )
  
  for (idx in split_ids) {
    counts <- df$count[idx]
    cum_counts <- cumsum(counts)
    
    df$ymin[idx] <- c(0, head(cum_counts, -1))
    df$ymax[idx] <- cum_counts
  }
  
  # Log y-scale cannot display 0.
  df$ymin_plot <- pmax(df$ymin, 1)
  df$ymax_plot <- df$ymax
  
  df <- df[df$ymax_plot > df$ymin_plot, , drop = FALSE]
  
  # Keep only breaks inside the observed raw UMI range.
  max_count <- max(
    c(
      umis_real,
      unlist(lapply(sim_list, function(x) x$umis_sim), use.names = FALSE)
    ),
    na.rm = TRUE
  )
  
  x_breaks <- x_breaks[x_breaks <= max_count]
  
  if (!0 %in% x_breaks) {
    x_breaks <- c(0, x_breaks)
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
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      limits = range(breaks),
      name = "UMI count"
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      name = "Number of cells",
      expand = ggplot2::expansion(mult = c(0, 0.05))
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













## this block is for making my plots of average metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_metrics <- function(pred, truth) {
  pred <- unique(pred)
  truth <- unique(truth)
  
  n_inter <- length(intersect(pred, truth))
  n_pred <- length(pred)
  n_truth <- length(truth)
  n_union <- length(union(pred, truth))
  
  tibble(
    precision = ifelse(n_pred == 0, NA_real_, n_inter / n_pred),
    recall    = ifelse(n_truth == 0, NA_real_, n_inter / n_truth),
    jaccard   = ifelse(n_union == 0, NA_real_, n_inter / n_union),
    n_pred = n_pred,
    n_truth = n_truth,
    n_inter = n_inter,
    n_union = n_union
  )
}

extract_part <- function(x, k) {
  vapply(
    strsplit(x, "_", fixed = TRUE),
    function(z) if (length(z) >= k) z[[k]] else NA_character_,
    character(1)
  )
}


make_guide_metric_tbl <- function(sim_results) {
  guide_names <- rownames(sim_results$grna_matrix)
  
  imap_dfr(
    sim_results$all_assns,
    function(method_assns, method_name) {
      map_dfr(
        guide_names,
        function(grna_id) {
          set_metrics(
            pred = method_assns[[grna_id]],
            truth = sim_results$true_assns[[grna_id]]
          ) %>%
            mutate(
              method_name = method_name,
              grna_id = grna_id,
              .before = 1
            )
        }
      )
    }
  )
}


# Per-guide "tail mean": mean of the counts that are >= thresh. This is the
# grouping covariate -- a robust proxy for a guide's perturbed-cell expression
# level that needs no ground truth, so it works on real data. NA when no count
# reaches thresh.
guide_tail_mean <- function(counts, thresh = 10) {
  tail <- counts[counts >= thresh]
  if (length(tail) == 0L) NA_real_ else mean(tail)
}

# Cut a numeric vector into k groups. method = "quantile" (default) gives
# roughly equal-count groups (sensible for a grouping variable used to average
# metrics across guides); "interval" gives k equal-width bins (base cut(x, k)).
# Quantile breaks are de-duplicated, so heavily tied data may yield < k groups.
cut_into_groups <- function(x, k, method = c("quantile", "interval")) {
  method <- match.arg(method)
  if (method == "interval") {
    return(cut(x, breaks = k))
  }
  qs <- stats::quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  qs <- unique(qs)
  if (length(qs) < 2L) {
    return(factor(rep(NA_character_, length(x))))
  }
  cut(x, breaks = qs, include.lowest = TRUE)
}

# Per-guide precision / recall / jaccard with ONE run treated as ground truth,
# so it can be run on real data (no true_assns needed). For every method other
# than `ref_method`, metrics are computed against `ref_method`'s assignments.
# The grouping covariate is each guide's tail mean (guide_tail_mean over the
# shared grna_matrix), cut into k groups; the numeric `tail_mean` is kept too.
#
#   all_assns    : named list of methods; each element a named-by-grna_id list
#                  of assigned cell ids (same shape as sim_results$all_assns).
#   grna_matrix  : shared gRNA count matrix (rownames = grna_ids), the tail-mean
#                  source.
#   ref_method   : name in all_assns to treat as ground truth.
#   k            : number of tail-mean groups.
#   tail_thresh  : threshold for guide_tail_mean (default 10).
#   group_method : "quantile" (equal-count, default) or "interval" (equal-width).
#   include_ref  : keep the reference method in the output (trivially metric 1);
#                  default FALSE.
#
# Feed the result to summarize_metric_tbl(., group_vars =
# c("method_name", "tail_mean_group")) for grouped averages.
make_guide_metric_tbl_vs_ref <- function(all_assns, grna_matrix, ref_method,
                                         k = 4L, tail_thresh = 10,
                                         group_method = c("quantile", "interval"),
                                         include_ref = FALSE) {
  group_method <- match.arg(group_method)
  stopifnot(ref_method %in% names(all_assns))
  guide_names <- rownames(grna_matrix)
  ref_assns   <- all_assns[[ref_method]]

  # per-guide tail mean (the grouping covariate) + its k groups
  tail_mean <- vapply(
    guide_names,
    function(grna_id) guide_tail_mean(as.numeric(grna_matrix[grna_id, ]), tail_thresh),
    numeric(1)
  )
  names(tail_mean) <- guide_names
  tail_mean_group <- cut_into_groups(tail_mean, k = k, method = group_method)
  names(tail_mean_group) <- guide_names

  methods <- names(all_assns)
  if (!include_ref) methods <- setdiff(methods, ref_method)

  imap_dfr(
    all_assns[methods],
    function(method_assns, method_name) {
      map_dfr(
        guide_names,
        function(grna_id) {
          set_metrics(
            pred  = method_assns[[grna_id]],
            truth = ref_assns[[grna_id]]
          ) %>%
            mutate(
              method_name     = method_name,
              grna_id         = grna_id,
              tail_mean       = tail_mean[[grna_id]],
              tail_mean_group = tail_mean_group[[grna_id]],
              .before = 1
            )
        }
      )
    }
  )
}


mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
}

summarize_metric_tbl <- function(guide_metric_tbl, group_vars = "method_name") {
  guide_metric_tbl %>%
    pivot_longer(
      cols = c(precision, recall, jaccard),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(across(all_of(c(group_vars, "metric")))) %>%
    summarize(
      mean = mean_or_na(value),
      sd = sd_or_na(value),
      n_guides = n(),
      n_nonmissing = sum(!is.na(value)),
      .groups = "drop"
    )
}

make_avg_metrics_plot <- function(df, x_breaks=NULL) {
  if(is.null(x_breaks)) x_breaks <- c(200,400,800,1600)
  ggplot(
    df,
    aes(
      x = mu_pert,
      y = mean,
      group = method_name,
      color = method_name,
      shape = method_name
      # linetype = method_name
    )
  ) +
    geom_line() +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.04,
      alpha = 0.7
    ) +
    geom_point(size = 2) +
    facet_grid(np ~ metric) +
    scale_x_continuous(
      trans = "log2",
      breaks = x_breaks
    ) +
    scale_shape_manual(
      values = c(16, 17, 15, 3, 7, 8, 4, 5, 6, 9,10,11)
    ) +
    labs(
      x = expression(mu[pert]),
      y = "Mean",
      color = "Method",
      shape = "Method"
      # linetype = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom"
    )
}


# Average number of cells each method assigns per guide, within each level of
# the tail-mean grouping variable (guide_tail_mean over the shared grna_matrix,
# cut into k groups). No reference assignments -- this just counts assignments,
# so it needs no ground truth and runs on real data. Returns a ggplot.
#
#   all_assns    : named list of methods; each element a named-by-grna_id list
#                  of assigned cell ids (same shape as elsewhere).
#   grna_matrix  : shared gRNA count matrix (rownames = grna_ids), tail-mean source.
#   k            : number of tail-mean groups.
#   tail_thresh  : threshold for guide_tail_mean (default 10).
#   group_method : "quantile" (equal-count, default) or "interval" (equal-width).
plot_avg_n_assignments_by_group <- function(all_assns, grna_matrix,
                                            k = 4L, tail_thresh = 10,
                                            group_method = c("quantile", "interval")) {
  group_method <- match.arg(group_method)
  guide_names <- rownames(grna_matrix)

  tail_mean <- vapply(
    guide_names,
    function(grna_id) guide_tail_mean(as.numeric(grna_matrix[grna_id, ]), tail_thresh),
    numeric(1)
  )
  names(tail_mean) <- guide_names
  tail_mean_group <- cut_into_groups(tail_mean, k = k, method = group_method)
  names(tail_mean_group) <- guide_names

  tbl <- imap_dfr(
    all_assns,
    function(method_assns, method_name) {
      map_dfr(
        guide_names,
        function(grna_id) {
          tibble(
            method_name     = method_name,
            grna_id         = grna_id,
            tail_mean_group = tail_mean_group[[grna_id]],
            n_assigned      = length(method_assns[[grna_id]])
          )
        }
      )
    }
  )

  summ <- tbl %>%
    filter(!is.na(tail_mean_group)) %>%
    group_by(method_name, tail_mean_group) %>%
    summarize(
      mean_n_assigned = mean_or_na(n_assigned),
      sd              = sd_or_na(n_assigned),
      n_guides        = n(),
      .groups = "drop"
    )

  ggplot(
    summ,
    aes(
      x = tail_mean_group,
      y = mean_n_assigned,
      group = method_name,
      color = method_name,
      shape = method_name
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    scale_shape_manual(
      values = c(16, 17, 15, 3, 7, 8, 4, 5, 6, 9, 10, 11)
    ) +
    labs(
      x = "guide tail-mean group",
      y = "Mean # assignments per guide",
      color = "Method",
      shape = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


## this block is for seeing how a guide's UMIs vary with covariates ~~~~~~~~~~~~

# For ONE guide, show how its UMI counts vary across the levels of each covariate.
# Each ROW is a covariate (a column of `covariate_matrix`); within a row, the
# covariate is cut into `k` bins and each COLUMN is one bin, holding a histogram
# of the guide's log1p(UMI) over the cells in that bin. So you read down a row to
# see whether the guide's count distribution shifts as that covariate increases.
#
#   covariate_matrix : cells x covariates (e.g. the design matrix, or a raw
#                      covariate matrix with columns grna_n_nonzero, grna_n_umis,
#                      ...). Rows must align with `guide_umis`.
#   guide_umis       : the guide's per-cell UMI count vector (length = n_cells =
#                      nrow(covariate_matrix)).
#   label            : optional guide name for the plot title.
#   k                : number of covariate bins per row (default 4).
#   covariates       : which covariate columns to use (names); default = all
#                      non-constant columns (so a model.matrix intercept is
#                      dropped automatically).
#   group_method     : "quantile" (equal-count bins, default) or "interval"
#                      (equal-width bins, i.e. base cut(x, k)). Quantile breaks
#                      are de-duplicated, so heavily-tied covariates may yield
#                      fewer than k bins.
#   bins             : histogram bins (default 40).
#   x_breaks         : raw-UMI tick positions (placed on the log1p scale).
#   facet_scales     : facet_grid scales (default "free_y"; y frees per row).
#   log_y            : log10 the count axis (default FALSE).
#
# Each panel is annotated with that bin's covariate range and its cell count n.
plot_umi_by_covariate_bins <- function(covariate_matrix, guide_umis,
                                       label        = NULL,
                                       k            = 4L,
                                       covariates   = NULL,
                                       group_method = c("quantile", "interval"),
                                       bins         = 40,
                                       x_breaks     = c(0, 1, 2, 5, 10, 20, 50,
                                                        100, 200, 500, 1000, 2000,
                                                        5000, 10000, 20000, 50000,
                                                        100000),
                                       facet_scales = "free_y",
                                       log_y        = FALSE) {
  group_method <- match.arg(group_method)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_umi_by_covariate_bins requires the ggplot2 package.")
  }
  # Accept a data.frame or a matrix. We pull columns individually rather than
  # as.matrix()-ing the whole thing, so mixed-type data frames (e.g. a batch
  # factor alongside numeric covariates) don't get coerced to a character matrix.
  n_cells <- nrow(covariate_matrix)
  guide_umis <- as.numeric(guide_umis)
  if (length(guide_umis) != n_cells) {
    stop("length(guide_umis) (", length(guide_umis), ") != nrow(covariate_matrix) (",
         n_cells, "); cells must align.")
  }

  is_df    <- is.data.frame(covariate_matrix)
  all_cols <- colnames(covariate_matrix)
  getcol   <- function(cv) if (is_df) covariate_matrix[[cv]] else covariate_matrix[, cv]

  # default: every numeric, non-constant column (drops a model.matrix intercept
  # and any non-numeric columns such as a batch factor).
  if (is.null(covariates)) {
    keep <- vapply(all_cols, function(cv) {
      col <- getcol(cv)
      is.numeric(col) && length(unique(col[is.finite(col)])) > 1L
    }, logical(1))
    covariates <- all_cols[keep]
    dropped <- setdiff(all_cols, covariates)
    if (length(dropped)) {
      message("Using numeric, non-constant covariates; dropping: ",
              paste(dropped, collapse = ", "))
    }
  } else {
    missing_cols <- setdiff(covariates, all_cols)
    if (length(missing_cols)) {
      stop("Covariate column(s) not found: ", paste(missing_cols, collapse = ", "))
    }
    nonnum <- covariates[!vapply(covariates, function(cv) is.numeric(getcol(cv)),
                                 logical(1))]
    if (length(nonnum)) {
      stop("Covariate column(s) not numeric (cannot bin): ",
           paste(nonnum, collapse = ", "))
    }
  }
  if (length(covariates) == 0L) stop("No numeric, non-constant covariate columns to bin.")

  y      <- guide_umis
  logumi <- log1p(y)

  make_groups <- function(x) {
    if (group_method == "interval") {
      cut(x, breaks = k, include.lowest = TRUE)
    } else {
      qs <- stats::quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
      qs <- unique(qs)
      if (length(qs) < 2L) return(factor(rep(NA_character_, length(x))))
      cut(x, breaks = qs, include.lowest = TRUE)
    }
  }

  parts <- lapply(covariates, function(cv) {
    f  <- make_groups(as.numeric(getcol(cv)))
    data.frame(
      covariate   = cv,
      group_idx   = as.integer(f),
      range_label = as.character(f),
      logumi      = logumi,
      row.names   = NULL,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, parts)
  df <- df[!is.na(df$group_idx), , drop = FALSE]

  df$covariate <- factor(df$covariate, levels = covariates)
  lvls <- sort(unique(df$group_idx))
  df$group_idx <- factor(df$group_idx, levels = lvls,
                         labels = paste0("bin ", lvls))

  # per-panel annotation: that bin's covariate range + cell count.
  ann <- df %>%
    dplyr::group_by(covariate, group_idx, range_label) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop")

  # keep only UMI breaks inside the observed range.
  max_count <- max(y, na.rm = TRUE)
  x_breaks  <- x_breaks[x_breaks <= max_count]
  if (!0 %in% x_breaks) x_breaks <- c(0, x_breaks)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = logumi)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = -Inf, y = Inf, label = paste0(range_label, "\nn = ", n)),
      hjust = -0.04, vjust = 1.15, size = 2.6, lineheight = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(covariate ~ group_idx, scales = facet_scales) +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      name   = "UMI count"
    ) +
    ggplot2::labs(
      y     = "number of cells",
      title = sprintf("UMI distribution of %s by covariate bin (%s, k = %d)",
                      label %||% "guide", group_method, k)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      axis.text.x  = ggplot2::element_text(size = 7, angle = 90,
                                           vjust = 0.5, hjust = 1)
    )

  if (log_y) {
    p <- p + ggplot2::scale_y_log10(name = "number of cells")
  }
  p
}


# For a "middling" UMI range, see what the covariates look like. Each ROW is a
# guide (from `grna_ids`). The FIRST panel of a row is that guide's full UMI
# histogram, exactly as plot_umi_histogram() draws it (no filtering -- the whole
# guide). Each SUBSEQUENT panel is a log-log histogram of one provided covariate,
# restricted to the cells whose UMI count y_i falls in `y_interval`. So you scan
# left (the whole count distribution, to pick/confirm the middling band) then
# right (what those middling cells' covariates look like). Returns a cowplot.
#
#   covariate_matrix : cells x covariates (data.frame or matrix); rows align with
#                      the columns of grna_matrix.
#   grna_matrix      : gRNAs x cells count matrix (rownames = grna_ids).
#   grna_ids         : guides to show (one row each).
#   covariates       : covariate column names (one covariate panel each).
#   y_interval       : length-2 c(lo, hi); covariate panels use cells with
#                      lo <= y_i <= hi (inclusive).
#   bins             : UMI-histogram bins for the first panel (default 80).
#   cov_bins         : covariate-histogram bins (default 40).
#   x_breaks         : raw-UMI tick positions for the first panel (log1p scale).
#   rel_widths       : optional cowplot rel_widths for the panels in a row
#                      (length 1 + length(covariates)); default equal widths.
plot_covariates_for_middling_y <- function(covariate_matrix, grna_matrix, grna_ids,
                                           covariates, y_interval,
                                           bins       = 80,
                                           cov_bins   = 40,
                                           x_breaks   = c(0, 1, 2, 5, 10, 20, 50,
                                                          100, 200, 500, 1000, 2000,
                                                          5000, 10000, 20000, 50000,
                                                          100000),
                                           rel_widths = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_covariates_for_middling_y requires the ggplot2 package.")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("plot_covariates_for_middling_y requires the cowplot package.")
  }
  stopifnot(length(y_interval) == 2L, y_interval[1] <= y_interval[2])
  lo <- y_interval[1]; hi <- y_interval[2]

  is_df    <- is.data.frame(covariate_matrix)
  n_cells  <- nrow(covariate_matrix)
  all_cols <- colnames(covariate_matrix)
  getcol   <- function(cv) if (is_df) covariate_matrix[[cv]] else covariate_matrix[, cv]

  if (ncol(grna_matrix) != n_cells) {
    stop("ncol(grna_matrix) (", ncol(grna_matrix), ") != nrow(covariate_matrix) (",
         n_cells, "); cells must align.")
  }
  miss_g <- setdiff(grna_ids, rownames(grna_matrix))
  if (length(miss_g)) stop("grna_ids not in grna_matrix: ", paste(miss_g, collapse = ", "))
  miss_c <- setdiff(covariates, all_cols)
  if (length(miss_c)) stop("covariate column(s) not found: ", paste(miss_c, collapse = ", "))
  nonnum <- covariates[!vapply(covariates, function(cv) is.numeric(getcol(cv)), logical(1))]
  if (length(nonnum)) stop("covariate column(s) not numeric: ", paste(nonnum, collapse = ", "))

  if (is.null(rel_widths)) rel_widths <- rep(1, 1L + length(covariates))

  # log-log histogram of a positive numeric vector (manual bars from 1, like
  # plot_umi_histogram, so the log10 y-axis has no log(0) trouble).
  loglog_hist <- function(v, title, n_sel) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) == 0L) {
      # text-only placeholder: an empty ggplot would train scales on no data
      # (non-finite viewport / min-max warnings), so draw a label instead.
      return(cowplot::ggdraw() +
               cowplot::draw_label(sprintf("%s\nn = %d (no cells > 0)", title, n_sel),
                                   size = 10))
    }
    h <- hist(log10(v), breaks = cov_bins, plot = FALSE)
    d <- data.frame(xmin  = 10^head(h$breaks, -1),
                    xmax  = 10^tail(h$breaks, -1),
                    count = h$counts)
    d <- d[d$count > 0, , drop = FALSE]
    ggplot2::ggplot(d) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = count)) +
      ggplot2::scale_x_log10(name = title) +
      ggplot2::scale_y_log10(labels = scales::label_number(), name = "number of cells") +
      ggplot2::labs(subtitle = sprintf("n = %d cells", n_sel)) +
      ggplot2::theme_classic()
  }

  grna_rownames <- rownames(grna_matrix)
  rows <- lapply(grna_ids, function(id) {
    # Extract by INTEGER row index, one whole row at a time: ondisc ODM objects
    # only support single whole-row access, and character-rowname indexing can
    # silently return nothing -- match() to an integer index sidesteps that.
    ri <- match(id, grna_rownames)
    y  <- as.numeric(grna_matrix[ri, ])
    if (length(y) != n_cells) {
      stop("Guide '", id, "' (row ", ri, ") returned ", length(y),
           " counts but there are ", n_cells, " cells; row extraction failed.")
    }
    umi   <- plot_umi_histogram(counts = y, bins = bins, title = id, x_breaks = x_breaks)
    sel   <- which(y >= lo & y <= hi)
    n_sel <- length(sel)
    cov_panels <- lapply(covariates, function(cv) {
      loglog_hist(as.numeric(getcol(cv))[sel], title = cv, n_sel = n_sel)
    })
    cowplot::plot_grid(plotlist = c(list(umi), cov_panels), nrow = 1,
                       rel_widths = rel_widths)
  })

  cowplot::plot_grid(plotlist = rows, ncol = 1,
                     rel_heights = rep(1, length(rows)))
}


# Like plot_covariates_for_middling_y, but each covariate panel OVERLAYS two
# histograms split by a UMI threshold: the covariate values for cells with
# y_min <= y_i <= y_thresh (the "middling" band) vs. for cells with y_i >
# y_thresh (the "high" cells). Lets you see whether a covariate distinguishes
# middling-count cells from clearly-perturbed ones. The first panel of each row
# is still the guide's full UMI histogram (plot_umi_histogram). Returns a cowplot.
#
#   covariate_matrix : cells x covariates (data.frame or matrix); rows align with
#                      the columns of grna_matrix.
#   grna_matrix      : gRNAs x cells count matrix (rownames = grna_ids).
#   grna_ids         : guides to show (one row each).
#   covariates       : covariate column names (one overlay panel each).
#   y_min, y_thresh  : the band is y_min <= y_i <= y_thresh; "high" is y_i >
#                      y_thresh (defaults y_min = 1, y_thresh = 50).
#   bins             : UMI-histogram bins for the first panel (default 80).
#   cov_bins         : covariate-histogram bins (default 40).
#   x_breaks         : raw-UMI tick positions for the first panel (log1p scale).
#   rel_widths       : optional cowplot rel_widths for the panels in a row.
plot_covariates_split_by_y <- function(covariate_matrix, grna_matrix, grna_ids,
                                       covariates,
                                       y_min      = 1,
                                       y_thresh   = 50,
                                       bins       = 80,
                                       cov_bins   = 40,
                                       x_breaks   = c(0, 1, 2, 5, 10, 20, 50,
                                                      100, 200, 500, 1000, 2000,
                                                      5000, 10000, 20000, 50000,
                                                      100000),
                                       rel_widths = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_covariates_split_by_y requires the ggplot2 package.")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("plot_covariates_split_by_y requires the cowplot package.")
  }
  stopifnot(length(y_min) == 1L, length(y_thresh) == 1L, y_min <= y_thresh)

  is_df    <- is.data.frame(covariate_matrix)
  n_cells  <- nrow(covariate_matrix)
  all_cols <- colnames(covariate_matrix)
  getcol   <- function(cv) if (is_df) covariate_matrix[[cv]] else covariate_matrix[, cv]

  if (ncol(grna_matrix) != n_cells) {
    stop("ncol(grna_matrix) (", ncol(grna_matrix), ") != nrow(covariate_matrix) (",
         n_cells, "); cells must align.")
  }
  miss_g <- setdiff(grna_ids, rownames(grna_matrix))
  if (length(miss_g)) stop("grna_ids not in grna_matrix: ", paste(miss_g, collapse = ", "))
  miss_c <- setdiff(covariates, all_cols)
  if (length(miss_c)) stop("covariate column(s) not found: ", paste(miss_c, collapse = ", "))
  nonnum <- covariates[!vapply(covariates, function(cv) is.numeric(getcol(cv)), logical(1))]
  if (length(nonnum)) stop("covariate column(s) not numeric: ", paste(nonnum, collapse = ", "))

  if (is.null(rel_widths)) rel_widths <- rep(1, 1L + length(covariates))

  lab_lo <- sprintf("%g <= y <= %g", y_min, y_thresh)
  lab_hi <- sprintf("y > %g", y_thresh)
  fills  <- stats::setNames(c("#4C72B0", "#DD8452"), c(lab_lo, lab_hi))

  # Two overlaid log-log histograms (manual bars from 1, shared breaks so the two
  # groups align). `v_lo` / `v_hi` are the covariate values for the band / high
  # cells; n_lo / n_hi are the group sizes (pre >0 filter) for the subtitle.
  overlap_hist <- function(v_lo, v_hi, title, n_lo, n_hi, show_legend) {
    d <- rbind(
      if (length(v_lo)) data.frame(v = v_lo, grp = lab_lo) else NULL,
      if (length(v_hi)) data.frame(v = v_hi, grp = lab_hi) else NULL
    )
    if (is.null(d)) d <- d[0, , drop = FALSE]
    d <- d[is.finite(d$v) & d$v > 0, , drop = FALSE]
    if (nrow(d) == 0L) {
      return(cowplot::ggdraw() +
               cowplot::draw_label(sprintf("%s\nno cells > 0", title), size = 10))
    }
    shared <- hist(log10(d$v), breaks = cov_bins, plot = FALSE)$breaks
    d$bin  <- cut(log10(d$v), breaks = shared, include.lowest = TRUE, labels = FALSE)
    rects  <- d %>%
      dplyr::group_by(grp, bin) %>%
      dplyr::summarize(count = dplyr::n(), .groups = "drop")
    rects$xmin <- 10^shared[rects$bin]
    rects$xmax <- 10^shared[rects$bin + 1L]
    rects$grp  <- factor(rects$grp, levels = c(lab_lo, lab_hi))

    p <- ggplot2::ggplot(rects) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = count, fill = grp),
        alpha = 0.5) +
      ggplot2::scale_x_log10(name = title) +
      ggplot2::scale_y_log10(labels = scales::label_number(), name = "number of cells") +
      ggplot2::scale_fill_manual(values = fills, drop = FALSE, name = NULL) +
      ggplot2::labs(subtitle = sprintf("n: %d / %d", n_lo, n_hi)) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = if (show_legend) "bottom" else "none")
    p
  }

  grna_rownames <- rownames(grna_matrix)
  rows <- lapply(grna_ids, function(id) {
    ri <- match(id, grna_rownames)
    y  <- as.numeric(grna_matrix[ri, ])
    if (length(y) != n_cells) {
      stop("Guide '", id, "' (row ", ri, ") returned ", length(y),
           " counts but there are ", n_cells, " cells; row extraction failed.")
    }
    umi    <- plot_umi_histogram(counts = y, bins = bins, title = id, x_breaks = x_breaks)
    sel_lo <- which(y >= y_min & y <= y_thresh)
    sel_hi <- which(y > y_thresh)
    cov_panels <- lapply(seq_along(covariates), function(j) {
      cv <- covariates[[j]]
      x  <- as.numeric(getcol(cv))
      overlap_hist(x[sel_lo], x[sel_hi], title = cv,
                   n_lo = length(sel_lo), n_hi = length(sel_hi),
                   show_legend = (j == 1L))
    })
    cowplot::plot_grid(plotlist = c(list(umi), cov_panels), nrow = 1,
                       rel_widths = rel_widths)
  })

  cowplot::plot_grid(plotlist = rows, ncol = 1,
                     rel_heights = rep(1, length(rows)))
}
