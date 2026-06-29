#!/usr/bin/env Rscript
# =============================================================================
# sim_de_lib.R  --  foundation for the downstream SCEPTRE DE pushthrough.
#
# Goal: take each of our 11 assignment methods' (guide, cell) calls on a real
# dataset (Gasperini, Replogle, ...) and run sceptre's positive-control (power)
# + negative-control (calibration) DE tests on the resulting binary indicator
# matrix.  Compare methods on the metric that matters: NT-pair type-I error +
# on-target gene-effect recovery.
#
# The lab already has a complete DE infrastructure:
#   benchmarking/association/data-preparation/  -- builds pos/neg ctrl datasets
#   benchmarking/association/pos-control-pipeline/bin/run_sceptre.R  -- pos run
#   benchmarking/association/neg-control-pipeline/bin/run_sceptre.R  -- neg run
# Their pipelines do `assign_grnas(method="thresholding", threshold=1)` on an
# already-BINARY grna_matrix, so we override that grna_matrix with each method's
# logical assignment matrix and re-run.
#
# Inputs live in lab-owned ondisc-backed sceptre objects:
#   ~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/<ds>/
#     sceptre-pipeline/sceptre_object_initial.rds  (analysis status = imported)
#     sceptre-pipeline/response.odm                (real GEX, ondisc)
#     sceptre-pipeline/grna.odm                    (real guide counts, ondisc)
#
# This library provides:
#   load_real_dataset(name)         load ondisc-backed sceptre object + the raw
#                                   grna count matrix (for OUR methods to consume)
#   run_methods_on_real(counts)     run the 11-method panel; returns a list of
#                                   logical assignment matrices (guides x cells)
#   build_de_inputs(<assignment>,   build a self-contained pos-control +
#                   <real_data>,     neg-control input directory (see lab's
#                   <method_name>)   data-preparation/ for the exact layout)
#   run_de(<input_dir>, ...)        invoke the lab's run_sceptre.R wrappers
#   score_de(<output_dir>)          read pos/neg results -> calibration + power
#
# This file is the FOUNDATION; concrete runs live in scripts/sim_de_<dataset>.R.
# Pipeline status is upstream-only assignment scoring until this is exercised.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(sceptre)
})

.here <- function() {
  if (file.exists(file.path(getwd(), "scripts", "sim_de_lib.R"))) return(normalizePath(getwd()))
  if (file.exists(file.path(getwd(), "..", "scripts", "sim_de_lib.R"))) return(normalizePath(file.path(getwd(), "..")))
  getwd()
}
HERE <- .here()
GA   <- normalizePath(file.path(HERE, ".."))
LAB_DATA  <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
LAB_ASSOC <- path.expand("~/data/projects/sceptre3/benchmarking/association")
DE_OUT    <- file.path(HERE, "results", "sim_framework", "de")
dir.create(DE_OUT, recursive = TRUE, showWarnings = FALSE)

# also borrow the lab's data-prep helpers (odm_to_sparse_matrix, etc.)
.lab_prep <- file.path(GA, "..", "association", "data-preparation")
if (file.exists(file.path(.lab_prep, "utils_data_prep.R"))) {
  source(file.path(.lab_prep, "utils_data_prep.R"))
  source(file.path(.lab_prep, "neg_control_functions.R"))
}

# Dataset registry.  Two storage modes:
#  - "ondisc": lab's pre-built sceptre_object_initial.rds + response.odm +
#    grna.odm (the lab's reference datasets Gasperini and Replogle).
#  - "survey": in-memory triple in <dir>/sceptre/ -- response_matrix.rds +
#    grna_matrix_aligned.rds, produced by sim_de_extract_survey.R from the raw
#    perturbseq-survey/ files (smaller in-memory-friendly datasets).
# See DATASETS_STATUS.md for the per-dataset audit.
SURV_BASE <- path.expand("~/data/external/perturbseq-survey")
REAL_DATASETS <- list(
  gasperini = list(mode = "ondisc", moi = "high",
    dir = file.path(LAB_DATA, "gasperini", "sceptre-pipeline")),
  replogle = list(mode = "ondisc", moi = "low",
    dir = file.path(LAB_DATA, "replogle-rd7", "sceptre-pipeline")),
  a549                = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "a549_crispri_perturbseq_GSE236304", "sceptre")),
  cd8_tcell           = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "cd8_tcell_perturbcite_GSE279498", "sceptre")),
  dctap_highmoi       = list(mode = "survey", moi = "high",
    dir = file.path(SURV_BASE, "dctap_k562_highmoi_GSE303901", "sceptre")),
  dctap_lowmoi        = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "dctap_k562_lowmoi_GSE303901", "sceptre")),
  endoc               = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "endoc_t2d_perturbseq_GSE273677", "sceptre")),
  multiome_erythroid  = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "perturb_multiome_erythroid_GSE274113", "sceptre")),
  tcell_cd4           = list(mode = "survey", moi = "low",
    dir = file.path(SURV_BASE, "tcell_cd4_perturbseq_GSE314342", "sceptre"))
)

# ---- (1) load a real dataset --------------------------------------------------
# Loads the lab's ondisc-backed sceptre_object_initial.rds + the real grna count
# matrix (as sparse for our methods).  Optional row/col subsetting at load time
# for smoke tests and downsampling on a memory-constrained host (full Replogle
# grna_mat is ~10 GB densified; with default subsetting it's <100 MB sparse).
# Returns:
#   $obj         the ondisc-backed sceptre object (analysis_status = imported)
#   $grna_mat    guides x cells sparse: what OUR methods consume
#   $cells_idx   the cell indices in the FULL dataset (for later GEX pulls)
#   $guides_idx  the guide indices in the FULL dataset
#   $meta        MOI + formula + extra covariates
load_real_dataset <- function(name, guides_idx = NULL, cells_idx = NULL,
                              max_cells = NULL, max_guides = NULL, seed = 1) {
  spec <- REAL_DATASETS[[name]]; if (is.null(spec)) stop("unknown dataset: ", name)
  mode <- spec$mode %||% "ondisc"
  if (mode == "survey") return(.load_survey(name, spec, guides_idx, cells_idx, max_cells, max_guides, seed))
  # ondisc path (lab's Gasperini / Replogle)
  obj  <- sceptre::read_ondisc_backed_sceptre_object(
            sceptre_object_fp     = file.path(spec$dir, "sceptre_object_initial.rds"),
            response_odm_file_fp  = file.path(spec$dir, "response.odm"),
            grna_odm_file_fp      = file.path(spec$dir, "grna.odm"))
  grna_odm <- obj@grna_matrix[[1]]
  G_full <- nrow(grna_odm); C_full <- ncol(grna_odm)
  if (is.null(guides_idx)) guides_idx <- seq_len(G_full)
  if (is.null(cells_idx))  cells_idx  <- seq_len(C_full)
  if (!is.null(max_guides) && length(guides_idx) > max_guides) {
    set.seed(seed); guides_idx <- sort(sample(guides_idx, max_guides))
  }
  if (!is.null(max_cells) && length(cells_idx) > max_cells) {
    set.seed(seed); cells_idx <- sort(sample(cells_idx, max_cells))
  }
  G <- length(guides_idx); C <- length(cells_idx)
  cat(sprintf("[%s] loading gRNA submatrix: %d guides x %d cells\n", name, G, C))
  rn <- rownames(grna_odm); cn <- colnames(grna_odm)
  guide_names <- if (!is.null(rn)) rn[guides_idx] else as.character(guides_idx)
  ilist <- jlist <- xlist <- vector("list", G)
  for (k in seq_len(G)) {
    v_full <- as.numeric(grna_odm[guide_names[k], ])
    v <- v_full[cells_idx]; nz <- which(v > 0)
    if (length(nz)) { ilist[[k]] <- rep.int(k, length(nz)); jlist[[k]] <- nz; xlist[[k]] <- v[nz] }
  }
  grna_mat <- Matrix::sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
    dims = c(G, C),
    dimnames = list(guide_names, if (!is.null(cn)) cn[cells_idx] else as.character(cells_idx)))
  list(obj = obj, grna_mat = grna_mat, cells_idx = cells_idx, guides_idx = guides_idx,
       mode = "ondisc", meta = spec, name = name)
}

# Survey-mode loader: in-memory response + aligned guide matrix already on disk.
# Returns a structure compatible with the ondisc one: a placeholder sceptre_object
# (built by import_data on the in-memory matrices) + the guide matrix + subset indices.
.load_survey <- function(name, spec, guides_idx, cells_idx, max_cells, max_guides, seed) {
  resp <- readRDS(file.path(spec$dir, "response_matrix.rds"))
  grna <- readRDS(file.path(spec$dir, "grna_matrix_aligned.rds"))
  stopifnot(ncol(resp) == ncol(grna))   # aligned by construction
  G_full <- nrow(grna); C_full <- ncol(grna)
  if (is.null(guides_idx)) guides_idx <- seq_len(G_full)
  if (is.null(cells_idx))  cells_idx  <- seq_len(C_full)
  if (!is.null(max_guides) && length(guides_idx) > max_guides) {
    set.seed(seed); guides_idx <- sort(sample(guides_idx, max_guides))
  }
  if (!is.null(max_cells) && length(cells_idx) > max_cells) {
    set.seed(seed); cells_idx <- sort(sample(cells_idx, max_cells))
  }
  grna_mat <- grna[guides_idx, cells_idx, drop = FALSE]
  resp_sub <- resp[, cells_idx, drop = FALSE]
  cat(sprintf("[%s/survey] %d genes x %d cells, %d guides\n",
              name, nrow(resp_sub), ncol(resp_sub), nrow(grna_mat)))
  # Stash response in a virtual obj struct so build_pos/neg_control_input has a
  # uniform interface.  We do NOT build a full sceptre_object here -- the
  # downstream builders read $resp_mem directly when mode=="survey".
  list(grna_mat = grna_mat, cells_idx = cells_idx, guides_idx = guides_idx,
       mode = "survey", resp_mem = resp_sub, meta = spec, name = name,
       # placeholders so other code that reads $obj@... doesn't break (it should
       # use mode to branch).
       obj = NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- (2) run our methods on the real grna_mat -------------------------------
# Reuses sim_lib.R::run_methods.  The dataset's grna_target_df + the ondisc
# response are not needed at this stage -- we only need the count matrix.
run_methods_on_real <- function(real, methods = c("thresh3","thresh10","ambient",
                                "otsu","valley","fishash","depth_fix","sceptre"),
                                external = list(), verbose = TRUE) {
  source(file.path(HERE, "scripts", "sim_lib.R"))   # provides run_methods()
  sc_moi <- if (real$meta$moi == "low") "low" else "high"
  run_methods(real$grna_mat, methods = methods, external = external,
              sceptre_moi = sc_moi, verbose = verbose)
}

# ---- (3) save assignment matrices for later DE inputs -----------------------
# Lives at results/sim_framework/de/<dataset>/assignments/<method>.rds
save_assignments <- function(name, assignments) {
  d <- file.path(DE_OUT, name, "assignments"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (m in names(assignments)) {
    saveRDS(assignments[[m]], file.path(d, paste0(m, ".rds")))
    cat(sprintf("  saved %s/%s.rds  (%d calls)\n", name, m, sum(assignments[[m]])))
  }
  invisible(d)
}

# ---- (4) build minimal pos/neg control input directories --------------------
# These produce the exact files the lab's run_sceptre.R runners expect, but
# minimally (no Mixscale/FRPerturb outputs needed for our DE pushthrough).
# PER METHOD: we substitute the lab's `assign_grnas(threshold=1)` step by
# providing a grna_matrix that IS the method's binary assignment.

# Pull a small dense response submatrix from the ondisc-backed gene matrix.
.pull_response <- function(response_odm, gene_ids, cell_idx) {
  ilist <- jlist <- xlist <- vector("list", length(gene_ids))
  for (k in seq_along(gene_ids)) {
    v <- as.numeric(response_odm[gene_ids[k], ])[cell_idx]; nz <- which(v > 0)
    if (length(nz)) { ilist[[k]] <- rep.int(k, length(nz)); jlist[[k]] <- nz; xlist[[k]] <- v[nz] }
  }
  Matrix::sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
                       dims = c(length(gene_ids), length(cell_idx)),
                       dimnames = list(gene_ids, NULL))
}

# build_pos_control_input(): for the on-target power test.
#   - Sample up to n_targets on-target guide groups that the assignment hits in
#     >= min_cells_per_target cells (so power isn't bounded by sample size).
#   - Take every gRNA in those targets + all non-targeting gRNAs.
#   - Take all cells assigned to any sampled guide (downsample to max_cells).
#   - Pull the corresponding gene from the ondisc response matrix.
# Writes sceptre/{response_matrix.rds, grna_matrix.rds, grna_target_data_frame.csv,
#                 cell_covariates.csv} into <out_dir>.
build_pos_control_input <- function(real, assignment, out_dir,
                                    n_targets = 50, max_cells = 20000,
                                    min_cells_per_target = 30, nt_name = "non-targeting",
                                    grna_target_df = NULL, seed = 1) {
  set.seed(seed)
  # branch by storage mode -- ondisc has the rich sceptre_object; survey has
  # in-memory response + a user-supplied or default target df.
  if (real$mode == "ondisc") {
    scep_obj <- real$obj
    grna_tdf <- scep_obj@grna_target_data_frame
    resp_odm <- scep_obj@response_matrix[[1]]
    cov_df   <- scep_obj@covariate_data_frame
    cells_idx <- real$cells_idx
    pull_resp <- function(gene_ids, cells) .pull_response(resp_odm, gene_ids, cells_idx[cells])
    all_resp_rownames <- rownames(resp_odm)
    have_cov <- TRUE
  } else {
    # survey mode: caller must pass a grna_target_df (we don't have one canonical).
    if (is.null(grna_target_df))
      stop("build_pos_control_input(): survey-mode dataset requires `grna_target_df = ...`; ",
           "construct guide->gene-target mapping from the dataset's guide_features.csv.")
    grna_tdf <- grna_target_df
    resp_mat <- real$resp_mem
    pull_resp <- function(gene_ids, cells) resp_mat[gene_ids, cells, drop = FALSE]
    all_resp_rownames <- rownames(resp_mat)
    # compute basic covariates from the in-memory matrices
    cov_df <- data.frame(
      response_n_umis    = Matrix::colSums(real$resp_mem),
      response_n_nonzero = Matrix::colSums(real$resp_mem > 0),
      grna_n_umis        = Matrix::colSums(real$grna_mat),
      grna_n_nonzero     = Matrix::colSums(real$grna_mat > 0)
    )
    cells_idx <- seq_len(ncol(real$resp_mem))
    have_cov <- TRUE
  }
  # 1. Restrict targets to those whose guides are in our assignment AND that the
  #    method has called in enough cells; sample n_targets of them.
  guide_names_in_A <- rownames(assignment)
  tdf_sub <- grna_tdf[grna_tdf$grna_id %in% guide_names_in_A & grna_tdf$grna_target != nt_name, ]
  on_targets_pool <- intersect(unique(tdf_sub$grna_target), all_resp_rownames)
  guide_cells <- Matrix::rowSums(assignment)                   # cells per guide
  per_target <- tapply(guide_cells[tdf_sub$grna_id], tdf_sub$grna_target, sum)
  ok <- intersect(names(per_target)[per_target >= min_cells_per_target], on_targets_pool)
  if (length(ok) < 1) stop("no on-targets with >= min_cells_per_target cells in this assignment")
  targets <- sample(ok, min(n_targets, length(ok)))
  cat(sprintf("  pos: %d targets selected (from %d eligible)\n", length(targets), length(ok)))
  # 2. Guides for those targets + NTs WITH >=min_cells_per_target cells
  # assigned (very-rare NTs trip sceptre's NT-cell update step in run_qc).
  use_targets <- c(targets, nt_name)
  cand_grna <- grna_tdf$grna_id[grna_tdf$grna_target %in% use_targets]
  cand_grna <- intersect(cand_grna, rownames(assignment))
  # filter NTs to those with enough cells assigned
  nt_in <- cand_grna[grna_tdf$grna_target[match(cand_grna, grna_tdf$grna_id)] == nt_name]
  nt_keep <- nt_in[Matrix::rowSums(assignment[nt_in, , drop = FALSE]) >= min_cells_per_target]
  use_guides <- c(setdiff(cand_grna, nt_in), nt_keep)
  cat(sprintf("  pos: keeping %d NTs with >= %d cells (of %d total NTs in panel)\n",
              length(nt_keep), min_cells_per_target, length(nt_in)))
  # 3. Cells assigned to ANY of those guides (max_cells).
  A_sub <- assignment[use_guides, , drop = FALSE]
  cand_cells <- which(Matrix::colSums(A_sub) > 0)
  if (length(cand_cells) > max_cells) cand_cells <- sort(sample(cand_cells, max_cells))
  cat(sprintf("  pos: %d guides x %d cells\n", length(use_guides), length(cand_cells)))
  # 4. Pull the gene matrix for the selected targets (genes), all chosen cells.
  resp <- pull_resp(targets, cand_cells)
  nz_cells <- which(Matrix::colSums(resp) > 0)
  if (length(nz_cells) < length(cand_cells)) {
    cat(sprintf("  pos: dropping %d cells with all-zero gene expression\n", length(cand_cells) - length(nz_cells)))
    cand_cells <- cand_cells[nz_cells]; resp <- resp[, nz_cells, drop = FALSE]
  }
  # 5. Final pieces.
  grna_bin <- as(A_sub[, cand_cells, drop = FALSE], "lgCMatrix")
  storage.mode(grna_bin@x) <- "logical"
  # cov_df rows correspond to full-dataset cells in ondisc; to in-memory cells
  # in survey -- normalize via cells_idx.
  cov_sub  <- .rename_full(cov_df[cells_idx[cand_cells], , drop = FALSE])
  tdf_out  <- grna_tdf[match(use_guides, grna_tdf$grna_id), ]
  dir.create(file.path(out_dir, "sceptre"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(as(resp, "dgCMatrix"),   file.path(out_dir, "sceptre", "response_matrix.rds"))
  saveRDS(as(grna_bin + 0, "dgCMatrix"), file.path(out_dir, "sceptre", "grna_matrix.rds"))  # numeric 0/1
  write.csv(tdf_out, file.path(out_dir, "sceptre", "grna_target_data_frame.csv"), row.names = FALSE)
  write.csv(cov_sub, file.path(out_dir, "sceptre", "cell_covariates.csv"),       row.names = FALSE)
  invisible(list(n_targets = length(targets), n_guides = length(use_guides), n_cells = length(cand_cells)))
}

# The lab's pos-control runner expects covariate names with the "_full" suffix
# (their data-prep renames them).  We rename here so the lab's runner works as-is.
.rename_full <- function(df) {
  for (b in c("response_n_nonzero", "response_n_umis", "response_p_mito",
              "grna_n_nonzero", "grna_n_umis")) {
    if (b %in% names(df)) names(df)[names(df) == b] <- paste0(b, "_full")
  }
  df
}

# build_neg_control_input(): for the calibration test.
#   - Take all NT guides assigned by the method.
#   - Pick n_genes random "untargeted" genes (not in any guide target).
#   - All (NT_guide, random_gene) pairs are the discovery_pairs.rds; under the
#     null (no real perturbation) the p-values should be ~Uniform[0,1].
build_neg_control_input <- function(real, assignment, out_dir, n_genes = 100, max_cells = 30000,
                                    nt_name = "non-targeting", grna_target_df = NULL, seed = 2) {
  set.seed(seed)
  if (real$mode == "ondisc") {
    scep_obj <- real$obj
    grna_tdf <- scep_obj@grna_target_data_frame
    resp_odm <- scep_obj@response_matrix[[1]]
    cov_df   <- scep_obj@covariate_data_frame
    cells_idx <- real$cells_idx
    pull_resp <- function(gene_ids, cells) .pull_response(resp_odm, gene_ids, cells_idx[cells])
    all_genes <- rownames(resp_odm)
  } else {
    if (is.null(grna_target_df))
      stop("build_neg_control_input(): survey-mode requires `grna_target_df = ...`")
    grna_tdf <- grna_target_df
    resp_mat <- real$resp_mem
    pull_resp <- function(gene_ids, cells) resp_mat[gene_ids, cells, drop = FALSE]
    cov_df <- data.frame(
      response_n_umis    = Matrix::colSums(real$resp_mem),
      response_n_nonzero = Matrix::colSums(real$resp_mem > 0),
      grna_n_umis        = Matrix::colSums(real$grna_mat),
      grna_n_nonzero     = Matrix::colSums(real$grna_mat > 0)
    )
    cells_idx <- seq_len(ncol(real$resp_mem))
    all_genes <- rownames(real$resp_mem)
  }
  nt_guides <- intersect(grna_tdf$grna_id[grna_tdf$grna_target == nt_name], rownames(assignment))
  if (length(nt_guides) < 2) stop("need >= 2 NT guides for negative control")
  # Untargeted genes: those NOT named as any grna_target.  `all_genes` was set
  # to the relevant gene-id source per mode above.
  targeted  <- unique(grna_tdf$grna_target); targeted <- targeted[targeted != nt_name]
  untar     <- setdiff(all_genes, targeted)
  pick_genes <- sample(untar, min(n_genes, length(untar)))
  # Cells assigned to any NT or to NOT-assigned-anywhere
  A_nt <- assignment[nt_guides, , drop = FALSE]
  cand_cells <- which(Matrix::colSums(A_nt) > 0)
  if (length(cand_cells) > max_cells) cand_cells <- sort(sample(cand_cells, max_cells))
  cat(sprintf("  neg: %d NT guides, %d untargeted genes, %d cells\n",
              length(nt_guides), length(pick_genes), length(cand_cells)))
  resp <- pull_resp(pick_genes, cand_cells)
  # DO NOT filter cells based on the picked-gene subset library size -- the
  # picked genes are a random sample of untargeted genes, so a cell can have
  # 0 UMIs across this subset while being a perfectly good cell overall.
  # (pos filters globally on its target-gene subset because all-zero target
  # expression is a real signal floor; here we're testing whether p-values
  # are calibrated, which doesn't require any minimum gene expression.)
  grna_bin <- as(A_nt[, cand_cells, drop = FALSE], "lgCMatrix")
  cov_sub  <- .rename_full(cov_df[cells_idx[cand_cells], , drop = FALSE])  # avoid sceptre's reserved names
  # Relabel each NT guide as its own pseudo-target so sceptre's run_discovery_analysis
  # accepts them (NTs are forbidden in discovery_pairs; this is the lab's workaround
  # from neg_control_functions.R::create_pseudo_targets).  Under the null all p-values
  # should still be ~Uniform[0,1].
  pseudo_targets <- paste0("nt_pseudo_", nt_guides)
  tdf_out <- data.frame(grna_id = nt_guides, grna_target = pseudo_targets, stringsAsFactors = FALSE)
  disc_pairs <- expand.grid(grna_target = pseudo_targets, response_id = pick_genes, stringsAsFactors = FALSE)
  # Formula uses the renamed (_full) covariates.
  fmla_str <- paste0("~ log(response_n_umis_full + 1) + log(response_n_nonzero_full + 1) + ",
                     "log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1)")
  dir.create(file.path(out_dir, "sceptre"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(as(resp, "dgCMatrix"),       file.path(out_dir, "sceptre", "response_matrix.rds"))
  saveRDS(as(grna_bin + 0, "dgCMatrix"), file.path(out_dir, "sceptre", "grna_matrix.rds"))
  write.csv(tdf_out,    file.path(out_dir, "sceptre", "grna_target_data_frame.csv"), row.names = FALSE)
  write.csv(cov_sub,    file.path(out_dir, "sceptre", "cell_covariates.csv"),       row.names = FALSE)
  saveRDS(disc_pairs,   file.path(out_dir, "sceptre", "discovery_pairs.rds"))
  saveRDS(fmla_str,     file.path(out_dir, "sceptre", "formula_object.rds"))
  invisible(list(n_guides = length(nt_guides), n_genes = length(pick_genes), n_cells = length(cand_cells)))
}

# ---- (5) run lab's sceptre DE runners on a built input ----------------------
# Returns the path to the results CSV.  Set NCPUS env to enable parallel
# discovery analysis (neg control); pos is fast enough serial.
run_de_sceptre <- function(input_dir, kind = c("pos", "neg"), dataset_id) {
  kind <- match.arg(kind)
  # Two runners: the lab's dataset_id-aware one (Gasperini/Replogle) and our
  # generic survey runner.  Pick by whether dataset_id is one of the canonical
  # lab names.
  is_lab <- dataset_id %in% c("gasperini", "replogle")
  bin <- if (is_lab) {
    if (kind == "pos") file.path(GA, "..", "association", "pos-control-pipeline", "bin", "run_sceptre.R")
    else               file.path(GA, "..", "association", "neg-control-pipeline", "bin", "run_sceptre.R")
  } else {
    file.path(HERE, "scripts", "run_sceptre_survey.R")
  }
  stopifnot(file.exists(bin))
  # Write meta.json with kind + moi so the survey runner has them.
  meta <- list(kind = kind, dataset_id = dataset_id,
               moi = if (!is.null(REAL_DATASETS[[dataset_id]])) REAL_DATASETS[[dataset_id]]$moi else "low")
  writeLines(sprintf('{"kind":"%s","moi":"%s","dataset_id":"%s"}',
                     meta$kind, meta$moi, meta$dataset_id),
             file.path(input_dir, "sceptre", "meta.json"))
  out_csv <- if (kind == "pos") "association_on_target_sceptre.csv"
             else                "association_neg_control_sceptre.csv"
  wd0 <- getwd(); on.exit(setwd(wd0), add = TRUE)
  setwd(input_dir)
  cmd_args <- if (is_lab) c("--vanilla", bin, file.path(input_dir, "sceptre"), dataset_id)
              else        c("--vanilla", bin, file.path(input_dir, "sceptre"), kind)
  ec <- system2("Rscript", cmd_args, stdout = "_run.log", stderr = "_run.err")
  if (ec != 0) stop("run_sceptre.R (", kind, ") failed; see ", file.path(input_dir, "_run.err"))
  file.path(input_dir, out_csv)
}

# ---- (6) score the DE outputs -----------------------------------------------
# For pos: rejection rate at q = 0.1 (rough power; lower bound for any method).
# For neg: realized type-I error at nominal q = 0.05 and overall p-value
# distribution skew (KS distance from Uniform).
score_de <- function(csv_path, kind = c("pos", "neg"), q = 0.1) {
  kind <- match.arg(kind)
  r <- read.csv(csv_path)
  if (!"p_value" %in% names(r)) stop("expected 'p_value' column in ", csv_path)
  if (kind == "pos") {
    fdr <- stats::p.adjust(r$p_value, "BH")
    list(kind = "pos", n_pairs = nrow(r), power_q = mean(fdr < q),
         median_p = median(r$p_value, na.rm = TRUE))
  } else {
    list(kind = "neg", n_pairs = nrow(r),
         realized_t1e_005 = mean(r$p_value < 0.05, na.rm = TRUE),
         realized_t1e_001 = mean(r$p_value < 0.01, na.rm = TRUE),
         ks_uniform = ks.test(r$p_value[!is.na(r$p_value)], "punif")$statistic)
  }
}

cat("sim_de_lib.R loaded: load_real_dataset, run_methods_on_real, save_assignments,\n",
    "                     build_pos_control_input, build_neg_control_input,\n",
    "                     run_de_sceptre, score_de.\n", sep = "")
