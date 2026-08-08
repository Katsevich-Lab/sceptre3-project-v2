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
