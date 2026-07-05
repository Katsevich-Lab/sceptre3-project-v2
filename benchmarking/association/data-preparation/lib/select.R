# lib/select.R
# Per-type selectors: turn source data + a loaded assignment into a DatasetSpec
# describing WHAT goes into a dataset, without extraction or writing (that's
# assemble.R). A DatasetSpec is a plain list:
#
#   genes             character  candidate response genes (pre-filter)
#   candidate_cells   integer    cell indices into the original ODM columns
#   guides            character  guides to place in the grna indicator matrix
#   grna_target_df    df         grna_id -> grna_target for those guides
#   discovery_targets character  target names to test; NULL => no discovery_pairs
#
# All selectors draw cells only from the assignment matrix passed in, which the
# caller has already had load_assignment() zero on excluded cells -- so excluded
# cells can never be selected and exclusion never appears here.
#
# Currently implements the NEGATIVE-CONTROL selector + helpers; computational
# and positive-control selectors are added alongside later.


# ---- create_pseudo_targets ------------------------------------------------
# Each NT guide becomes its own singleton pseudo-target (non-targeting1, ...),
# so the negative control tests each NT guide as if it were a distinct target.
create_pseudo_targets <- function(guide_ids) {
  df <- data.frame(
    grna_id = guide_ids,
    grna_target = paste0("non-targeting", seq_along(guide_ids)),
    stringsAsFactors = FALSE
  )
  cat("Created", nrow(df), "pseudo-targets.\n")
  stopifnot(!any(duplicated(df$grna_target)))
  df
}


# ---- identify_nt_only_cells (low-MOI neg-control cells) --------------------
# Cells expressing NO targeting guide AND at least one NT guide. Under "maximum"
# assignment every usable cell has exactly one guide, so this is precisely the
# cells whose single guide is non-targeting. Ported from legacy.
identify_nt_only_cells <- function(assn, grna_target_df,
                                   nt_target_name = "non-targeting",
                                   targets_to_ignore = character(0)) {
  nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == nt_target_name]
  cat("Found", length(nt_guides), "NT guides total.\n")
  if (length(nt_guides) == 0) {
    stop("No non-targeting guides found with target name '", nt_target_name, "'")
  }

  non_targeting_targets <- c(nt_target_name, targets_to_ignore)
  targeting_guides <- grna_target_df$grna_id[!grna_target_df$grna_target %in% non_targeting_targets]
  stopifnot(targeting_guides %in% rownames(assn))
  cat("Found", length(targeting_guides), "targeting guides (excluding",
      paste(non_targeting_targets, collapse = ", "), ")\n")

  cell_idx <- intersect(
    which(Matrix::colSums(assn[targeting_guides, , drop = FALSE]) == 0),
    which(Matrix::colSums(assn[nt_guides, , drop = FALSE]) > 0)
  )
  cat("Found", length(cell_idx), "cells expressing no targeting guides and at least one NT guide.\n")
  if (length(cell_idx) == 0) stop("No cells found expressing only NT guides")
  cell_idx
}


# ---- cells_expressing_these_guides_and_not_gene_targets (high-MOI) ---------
# Cells expressing at least one of `guides_to_use` (the NT guides) but NOT any
# guide whose target is among the sampled `genes` -- so no retained cell is
# perturbed for any gene we will actually test. Ported from legacy.
cells_expressing_these_guides_and_not_gene_targets <- function(guides_to_use, genes, assn, grna_target_df) {
  guides_for_sampled_genes <- grna_target_df$grna_id[grna_target_df$grna_target %in% genes]
  cells_targeting_genes <- which(Matrix::colSums(assn[guides_for_sampled_genes, , drop = FALSE]) > 0)
  cells_for_guides_to_use <- which(Matrix::colSums(assn[guides_to_use, , drop = FALSE]) > 0)
  setdiff(cells_for_guides_to_use, cells_targeting_genes)
}


# ---- select_neg -----------------------------------------------------------
#' Negative-control DatasetSpec. Seeds and samples genes FIRST (matching the
#' legacy RNG order), then selects NT cells by the dataset's MOI regime.
#'
#' @param assn              zeroed assignment matrix from load_assignment()$assn
#' @param grna_target_df    the full source grna_target_df
#' @param genes_passing_qc  gene pool from compute_genes_passing_qc()
#' @param num_genes         number of genes to sample
#' @param dataset           a dataset config (dataset_gasperini / dataset_replogle)
#' @param nt_name           non-targeting target name
#' @param random_seed       seed set before the gene sample (legacy parity)
select_neg <- function(assn, grna_target_df, genes_passing_qc, num_genes,
                       dataset, nt_name = "non-targeting", random_seed = 243535) {
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  # Gene selection FIRST -- the only RNG draw here, mirroring legacy order.
  genes <- sample(genes_passing_qc, num_genes)
  cat("Selected", length(genes), "genes.\n")

  nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == nt_name]
  cat("There are", length(nt_guides), "NT guides.\n")

  # Cell selection -- MOI-specific, both correct (see selector notes above).
  if (dataset$moi$low_moi) {
    candidate_cells <- identify_nt_only_cells(assn, grna_target_df, nt_target_name = nt_name)
  } else {
    candidate_cells <- cells_expressing_these_guides_and_not_gene_targets(
      guides_to_use = nt_guides, genes = genes, assn = assn, grna_target_df = grna_target_df
    )
    cat("There are", length(candidate_cells),
        "cells expressing any NT guide and none for these genes.\n")
  }

  pseudo_df <- create_pseudo_targets(nt_guides)

  list(
    genes             = genes,
    candidate_cells   = candidate_cells,
    guides            = nt_guides,
    grna_target_df    = pseudo_df,
    discovery_targets = pseudo_df$grna_target
  )
}


# ---- select_targets_random ------------------------------------------------
# Ported from legacy. Reseeds internally (matches legacy RNG order). Excludes
# the given targets (we always pass the three-way NT/off-target/unknown set).
select_targets_random <- function(grna_target_df, num_targets, exclude = NULL, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)

  all_targets <- unique(grna_target_df$grna_target)
  if (!is.null(exclude)) {
    exclude <- intersect(exclude, grna_target_df$grna_target)   # absent names are a no-op safety net
    all_targets <- setdiff(all_targets, exclude)
  }
  cat("Found", length(all_targets), "candidate targets.\n")

  if (length(all_targets) < num_targets) {
    warning("Requested ", num_targets, " targets but only ", length(all_targets),
            " available. Using all available.")
    return(all_targets)
  }
  targets <- sample(all_targets, num_targets)
  cat("Selected", num_targets, "targets.\n")
  targets
}


# ---- sample_genes_with_expression -----------------------------------------
# Ported from legacy: shuffle all genes, keep those with any UMI in the given
# cells until num_genes reached. Only used when a driver passes a gene COUNT
# instead of an explicit gene vector.
sample_genes_with_expression <- function(response_odm, num_genes, cell_idx, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)
  genes <- character(0)
  all_genes <- sample(rownames(response_odm))
  for (curr_gene in all_genes) {
    if (sum(response_odm[curr_gene, ][cell_idx]) > 0) {
      genes <- c(genes, curr_gene)
      if (length(genes) >= num_genes) break
    }
  }
  genes
}


# ---- select_comp ----------------------------------------------------------
#' Computational-benchmarking DatasetSpec: random real targets, their guides
#' (+ NT if force_nt), all cells expressing them downsampled to max_num_cells,
#' random genes. NT keeps its real "non-targeting" label (not pseudo-targets).
#'
#' RNG order matches legacy: select_targets_random reseeds+samples, then the
#' cell downsample draws. `genes` may be an explicit character vector (already
#' sampled by the driver) or a single integer count (then sampled here).
#'
#' @param force_nt add all NT guides to the guide set (default = regime flag)
select_comp <- function(assn, response_odm, grna_target_df, genes, num_targets, max_num_cells,
                        dataset, nt_name = "non-targeting", random_seed = 243535,
                        force_nt = dataset$moi$force_nt_inclusion) {
  targets <- select_targets_random(
    grna_target_df, num_targets,
    exclude = c(nt_name, "nt_off_target", "unknown"), random_seed = random_seed
  )

  target_guides <- grna_target_df$grna_id[grna_target_df$grna_target %in% targets]
  guides <- if (force_nt) {
    union(target_guides, grna_target_df$grna_id[grna_target_df$grna_target == nt_name])
  } else target_guides

  # unname: which() carries cell barcodes as names (Gasperini); cell_idx is a
  # bare position vector (barcodes live in cell_info$cell_name).
  candidate_cells <- unname(which(Matrix::colSums(assn[guides, , drop = FALSE]) > 0))
  cat("Initially", length(candidate_cells), "cells found.\n")
  if (is.finite(max_num_cells)) {
    if (length(candidate_cells) < max_num_cells) stop("Not enough cells for these parameters.")
    cat("Downsampling to", max_num_cells, "cells.\n")
    candidate_cells <- sample(candidate_cells, max_num_cells)
  } else {
    cat("max_num_cells = Inf: keeping all cells for these targets.\n")
  }

  if (is.numeric(genes) && length(genes) == 1) {
    genes <- sample_genes_with_expression(response_odm, genes, candidate_cells)
  } else if (!is.character(genes)) {
    stop("`genes` must be a character vector or a single integer count.")
  }
  cat("Using", length(genes), "genes.\n")

  list(
    genes             = genes,
    candidate_cells   = candidate_cells,
    guides            = guides,
    grna_target_df    = grna_target_df[grna_target_df$grna_id %in% guides, , drop = FALSE],
    discovery_targets = targets                       # real sampled targets (not NT)
  )
}


# ---- select_pos -----------------------------------------------------------
#' Positive-control DatasetSpec: on-targets (targets that are also QC-passing
#' genes) with enough cells; genes == targets (each perturbation scored against
#' its OWN target gene); guides = target guides (+ NT if force_nt); cells
#' downsampled to max_num_cells. No discovery_pairs -- sceptre builds
#' positive_control_pairs itself, so discovery_targets = NULL.
select_pos <- function(assn, grna_target_df, genes_passing_qc, num_targets, max_num_cells,
                       only_consider_targets_with_this_many_cells, dataset,
                       nt_name = "non-targeting", random_seed = 34534,
                       force_nt = dataset$moi$force_nt_inclusion) {
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  all_on_targets <- intersect(grna_target_df$grna_target, genes_passing_qc)
  all_on_target_guides <- grna_target_df$grna_id[grna_target_df$grna_target %in% all_on_targets]

  # only on-targets with enough cells survive the later downsample
  enough <- data.frame(
    grna_id = all_on_target_guides,
    cells_per_guide = Matrix::rowSums(assn[all_on_target_guides, , drop = FALSE]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(grna_target_df[, c("grna_id", "grna_target")], by = "grna_id") |>
    dplyr::group_by(grna_target) |>
    dplyr::summarize(cells_per_target = sum(cells_per_guide), .groups = "drop") |>
    dplyr::filter(cells_per_target >= only_consider_targets_with_this_many_cells)
  if (nrow(enough) < num_targets) stop("Not enough on-targets have enough cells per target.")

  targets <- sample(enough$grna_target, num_targets)
  genes <- targets                                   # on-target: gene == target
  cat("There are", length(all_on_targets), "on-targets;",
      nrow(enough), "have >=", only_consider_targets_with_this_many_cells, "cells;",
      length(targets), "sampled.\n")

  target_guides <- grna_target_df$grna_id[grna_target_df$grna_target %in% targets]
  guides <- if (force_nt) {
    union(target_guides, grna_target_df$grna_id[grna_target_df$grna_target == nt_name])
  } else target_guides

  candidate_cells <- unname(which(Matrix::colSums(assn[guides, , drop = FALSE]) > 0))
  cat("Initially", length(candidate_cells), "cells found.\n")
  if (is.finite(max_num_cells)) {
    if (length(candidate_cells) < max_num_cells) stop("Not enough cells for these parameters.")
    cat("Downsampling to", max_num_cells, "cells.\n")
    candidate_cells <- sample(candidate_cells, max_num_cells)
  } else {
    cat("max_num_cells = Inf: keeping all", length(candidate_cells), "cells for these targets.\n")
  }

  list(
    genes             = genes,
    candidate_cells   = candidate_cells,
    guides            = guides,
    grna_target_df    = grna_target_df[grna_target_df$grna_id %in% guides, , drop = FALSE],
    discovery_targets = NULL                          # sceptre builds PC pairs itself
  )
}
