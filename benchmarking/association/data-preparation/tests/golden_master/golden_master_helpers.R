# golden_master_helpers.R
# Shared plumbing for the per-pipeline golden-master harnesses
# (golden_master_neg.R / _comp.R / _pos.R). Each of those is a thin, pipeline-
# specific entry point; the fake-data builder, legacy/new sourcing, FR-Perturb
# stub, and file comparison live here.
#
# The caller must set `.gm_dir` (this file's directory) before sourcing.

stopifnot(exists(".gm_dir"))
PREP_DIR <- normalizePath(file.path(.gm_dir, "..", ".."))            # .../data-preparation
FAKE_SRC <- normalizePath(file.path(PREP_DIR, "..", "..",
                          "guide-assignment", "data-preprocessing",
                          "create_fake_sceptre_data.R"))

# fresh temp output root for one harness run
make_out <- function(name) {
  o <- file.path(tempdir(), name)
  unlink(o, recursive = TRUE); dir.create(o, recursive = TRUE, showWarnings = FALSE)
  o
}

# ---- deterministic fake source bundle + assignment ------------------------
# for_pos=TRUE adds a response row per target_* name so that on-targets exist
# (positive control requires targets that are also response genes).
make_fake <- function(regime, for_pos = FALSE) {
  source(FAKE_SRC, local = TRUE); d <- make_data()   # internal set.seed(42)
  resp <- d$response_matrix; grna <- d$grna_matrix
  gtdf <- d$grna_target_data_frame; cov <- d$covariate_data_frame
  gtdf$grna_group <- gtdf$grna_target   # real Replogle gtdf has this; legacy pos-replogle does select(-grna_group)

  if (for_pos) {
    tnames <- setdiff(unique(gtdf$grna_target), "non-targeting")
    add <- matrix(rpois(length(tnames) * ncol(resp), 3), nrow = length(tnames),
                  dimnames = list(tnames, colnames(resp)))
    resp <- rbind(resp, add)
  }

  # nt_off_target & unknown targets exist in REPLOGLE only (legacy Replogle
  # select_targets_random() hard-errors without them; legacy Gasperini excludes
  # only non-targeting). Add them to the LOW regime alone -- adding them to
  # Gasperini would inflate its legacy target pool and make selection diverge
  # from the (correct) new three-way exclusion. Tiny counts => never assigned
  # (0 cells), never sampled (excluded); they only need to EXIST.
  if (regime == "low") {
    for (tt in c("nt_off_target", "unknown")) for (k in 1:2) {
      gid <- paste0(tt, "_", k)
      grna <- rbind(grna, matrix(1L, 1, ncol(grna), dimnames = list(gid, colnames(grna))))
      gtdf <- rbind(gtdf, data.frame(grna_id = gid, grna_target = tt, grna_group = tt,
                                     stringsAsFactors = FALSE))
    }
  }

  cov$response_n_nonzero <- as.integer(colSums(resp > 0))
  cov$response_n_umis    <- as.integer(colSums(resp))
  cov$grna_n_nonzero     <- as.integer(colSums(grna > 0))
  cov$grna_n_umis        <- as.integer(colSums(grna))
  cov$prep_batch         <- cov$batch
  rownames(cov) <- colnames(resp)

  if (regime == "low") {
    amax <- apply(grna, 2, which.max)                # primary guide is the unique max
    assn <- Matrix::sparseMatrix(i = amax, j = seq_len(ncol(grna)), x = 1, dims = dim(grna))
  } else {
    assn <- as(grna >= 5, "CsparseMatrix")
  }
  dimnames(assn) <- dimnames(grna)
  assn <- as(as(assn, "lMatrix"), "CsparseMatrix")   # lgCMatrix, like load_assignment

  list(resp = resp, grna = grna, gtdf = gtdf, cov = cov, assn = assn, genes = rownames(resp))
}

# stub for BOTH runs -- writes the frperturb covariate df + cell names, no h5ad/conda
frperturb_stub <- function(response_matrix, cell_names, cell_covariates_frpert, output_path, ...) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  saveRDS(as.data.frame(cell_covariates_frpert), file.path(output_path, "frperturb_covariates.rds"))
  saveRDS(cell_names, file.path(output_path, "frperturb_cell_names.rds"))
}

# extract & eval ONE `name <- function(...)` def from a script that also runs
# top-level code (the legacy pos builders live inside their driver scripts).
eval_fn_def <- function(file, fn_name) {
  for (e in parse(file)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        length(e) >= 2 && identical(as.character(e[[2]]), fn_name)) {
      eval(e, envir = globalenv()); return(invisible(TRUE))
    }
  }
  stop("definition of ", fn_name, " not found in ", file)
}

# source the legacy function set + mock config path + stub frperturb.
# `out` = harness output root; legacy writes under <out>/legacy.
source_legacy <- function(out) {
  source(file.path(PREP_DIR, "utils_data_prep.R"), local = FALSE)
  source(file.path(PREP_DIR, "neg_control_functions.R"), local = FALSE)
  source(file.path(PREP_DIR, "computational_benchmarking_functions.R"), local = FALSE)
  eval_fn_def(file.path(PREP_DIR, "make_pos_control_gasperini.R"), "make_pos_control_gasperini")
  eval_fn_def(file.path(PREP_DIR, "make_pos-control_replogle-rd7.R"), "make_pos_control_replogle_rd7")
  assign(".get_config_path", function(key) file.path(out, "legacy"), envir = globalenv())
  assign("write_frperturb_output", frperturb_stub, envir = globalenv())
}

source_new <- function() {
  for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R"))
    source(file.path(PREP_DIR, "lib", f), local = FALSE)
  assign("write_frperturb_output", frperturb_stub, envir = globalenv())
}

# ---- comparison -----------------------------------------------------------
load_any <- function(p) if (grepl("\\.rds$", p)) readRDS(p) else read.csv(p, stringsAsFactors = FALSE)
.nrows   <- function(x) if (is.data.frame(x) || is.matrix(x) || isS4(x)) nrow(x) else length(x)

.describe <- function(a, b) {
  na <- .nrows(a); nb <- .nrows(b)
  if (!is.null(na) && !is.null(nb) && na != nb) {
    tag <- if (nb < na) "new SUBSET (filter)" else "new LARGER (!)"
    cat(sprintf("        rows legacy %s / new %s  <- %s\n", na, nb, tag))
  }
  if (is.data.frame(a) && is.data.frame(b) && !identical(names(a), names(b))) {
    cat("        cols legacy:", paste(names(a), collapse = ","), "\n")
    cat("        cols new   :", paste(names(b), collapse = ","), "\n")
  }
}

# diff every shared output file between two dataset dirs
compare_dirs <- function(legacy_dir, new_dir) {
  files <- c("cell_info.csv", "cell_covariates.csv",
             "sceptre/response_matrix.rds", "sceptre/grna_matrix.rds",
             "sceptre/grna_target_data_frame.csv", "sceptre/discovery_pairs.rds",
             "sceptre/formula_object.rds", "sceptre/cell_covariates.csv",
             "mixscale/response_matrix.rds", "mixscale/assignments.rds",
             "mixscale/mixscale_nt_targets.rds",
             "frperturb/frperturb_covariates.rds", "frperturb/frperturb_cell_names.rds")
  for (f in files) {
    lp <- file.path(legacy_dir, f); np <- file.path(new_dir, f)
    le <- file.exists(lp); ne <- file.exists(np)
    if (!le && !ne) next
    if (le != ne) { cat(sprintf("  [ONLY %s] %s\n", if (le) "LEGACY" else "NEW", f)); next }
    lo <- load_any(lp); no <- load_any(np)
    if (isTRUE(all.equal(lo, no))) cat(sprintf("  [same]   %s\n", f))
    else { cat(sprintf("  [DIFFER] %s\n", f)); .describe(lo, no) }
  }
}
