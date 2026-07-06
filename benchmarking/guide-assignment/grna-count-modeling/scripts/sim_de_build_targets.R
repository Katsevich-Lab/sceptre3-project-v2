#!/usr/bin/env Rscript
# =============================================================================
# Build grna_target_data_frame.csv for each survey dataset that has usable GEX.
# Saves next to the existing sceptre/{response_matrix,grna_matrix_aligned}.rds.
#
# Sceptre needs (grna_id, grna_target) per guide, with "non-targeting" used as
# the canonical NT label.  For positive-control DE pushthrough, the target must
# match a gene id/name in the dataset's response_matrix.
#
# Per dataset, source of target info:
#   a549              h5 has /matrix/features/target_gene_name (canonical)
#   cd8_tcell         regex on guide name:    ^g(.+)_\d+$       (mouse)
#   dctap_*           prefix mapping (safe/nontarget/NC -> NT;
#                                     MYC-TSS / e-GATA1 -> MYC / GATA1 ...)
#   endoc             NT = ^Non-Targeting-*; targeting = ^([^-]+)-\d+$
#   invivo_cortex     NT = explicit list; targeting = ^([^-]+)-\d+$ (mouse)
#   tcell_cd4         targeting = ^([^-]+)-\d+$; NO NTs (genome-wide screen)
#   multiome_erythroid  guide names are pure sequences -> SKIPPED here
#
# After building, verify each target gene is present in the dataset's GEX rows
# (by SYMBOL or ENSG depending on what the response_matrix rownames are).
# =============================================================================
suppressPackageStartupMessages({library(Matrix); library(rhdf5)})

SURV <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
NT <- "non-targeting"

# Build a SYMBOL -> ENSG mapping for a dataset using its raw 10x features list.
# This is the lookup used to resolve guide-target symbols to GEX rownames.
# Returns named char vector (names = SYMBOL, values = ENSG / feature id).
build_symbol_to_id <- function(name_dir) {
  d <- file.path(SURV, name_dir)
  # First choice: the features_all.csv saved by sim_de_extract_survey.R
  feat_csv <- file.path(d, "sceptre", "features_all.csv")
  if (file.exists(feat_csv)) {
    f <- read.csv(feat_csv, stringsAsFactors = FALSE)
    gex <- f$type == "Gene Expression"
    map <- f$id[gex]; names(map) <- f$name[gex]
    return(map[!duplicated(names(map))])
  }
  # Second: h5
  h5 <- list.files(d, pattern = "\\.h5$", full.names = TRUE)[1]
  if (!is.na(h5)) {
    suppressMessages(library(rhdf5))
    fid   <- as.vector(h5read(h5, "/matrix/features/id"))
    fname <- as.vector(h5read(h5, "/matrix/features/name"))
    ftype <- as.vector(h5read(h5, "/matrix/features/feature_type"))
    gex <- ftype == "Gene Expression"
    map <- fid[gex]; names(map) <- fname[gex]
    return(map[!duplicated(names(map))])
  }
  warning("could not build symbol->id map for ", name_dir); return(character(0))
}

# Case-flexible lookup: try exact, then UPPER(symbol) -> UPPER(map names),
# then title-case (mouse).  Returns named char vec same length as symbols,
# values = the looked-up id or NA.
case_flex_lookup <- function(symbols, sym2id) {
  ids <- sym2id[symbols]                                      # exact match
  miss <- is.na(ids)
  if (any(miss) && length(sym2id)) {
    up_map <- sym2id; names(up_map) <- toupper(names(sym2id))
    up_map <- up_map[!duplicated(names(up_map))]
    ids2 <- up_map[toupper(symbols)]; ids[miss] <- ids2[miss]
    miss <- is.na(ids)
  }
  if (any(miss) && length(sym2id)) {
    # title-case: first letter upper, rest lower (mouse Mgi convention)
    tc <- function(x) ifelse(nchar(x) == 0, x,
                              paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2))))
    tc_map <- sym2id; names(tc_map) <- tc(names(sym2id))
    tc_map <- tc_map[!duplicated(names(tc_map))]
    ids3 <- tc_map[tc(symbols)]; ids[miss] <- ids3[miss]
  }
  ids
}

# Resolve symbols to GEX rownames; entries that don't resolve become NA.
# Pass-through for the NT label.  Uses case-flexible lookup so mouse/human
# capitalization differences don't break the resolution.
resolve_to_gex_rownames <- function(symbols, sym2id, resp_rownames) {
  out <- ifelse(symbols == NT, NT, case_flex_lookup(symbols, sym2id))
  # if a symbol happens to already match resp_rownames (resp uses symbols, rare),
  # accept that too
  fall_through <- is.na(out) & !is.na(symbols) & symbols %in% resp_rownames
  out[fall_through] <- symbols[fall_through]
  ok <- (out == NT) | (out %in% resp_rownames)
  out[!ok] <- NA_character_
  out
}

# ---- per-dataset extractors -------------------------------------------------
build_targets_a549 <- function() {
  dir <- file.path(SURV, "a549_crispri_perturbseq_GSE236304")
  resp <- readRDS(file.path(dir, "sceptre", "response_matrix.rds"))
  h5 <- file.path(dir, "filtered_feature_bc_matrix.h5")
  fid   <- as.vector(h5read(h5, "/matrix/features/id"))
  fname <- as.vector(h5read(h5, "/matrix/features/name"))
  ftype <- as.vector(h5read(h5, "/matrix/features/feature_type"))
  tname <- as.vector(h5read(h5, "/matrix/features/target_gene_name"))
  tid   <- as.vector(h5read(h5, "/matrix/features/target_gene_id"))
  g <- which(ftype == "CRISPR Guide Capture")
  # The lab's response_matrix uses GENE NAMES (the 10x features.tsv column 2),
  # so we resolve to gene symbol via target_gene_name.  Map "Non-Targeting" -> NT.
  target <- ifelse(tname[g] == "Non-Targeting", NT, tname[g])
  data.frame(grna_id = fid[g], grna_target = target,
             grna_group = NA_character_, stringsAsFactors = FALSE)
}

# Regex extractor: returns the captured group, or NT if name is in nt_set.
.parse_target_by_regex <- function(grna_id, regex, nt_set, toupper = TRUE) {
  target <- ifelse(grna_id %in% nt_set, NT, sub(regex, "\\1", grna_id))
  # unmatched names (no regex match) remain unchanged -> mark NA later
  unmatched <- !(target %in% c(NT)) & target == grna_id
  if (toupper) target <- ifelse(target == NT, NT, toupper(target))
  target[unmatched] <- NA_character_
  target
}

build_targets_simple <- function(name_dir, regex, nt_patterns, upper = TRUE) {
  # Use the SAME ids that are in the saved grna_matrix_aligned.rds -- the
  # lab's `guide_features.csv` is sometimes a supplemental list that disagrees
  # with the 10x feature list the matrix was built from (e.g. cd8_tcell where
  # gScramble3_1/5_2 are NTs in the matrix but NEG_CTRL-* in the lab csv).
  grna <- readRDS(file.path(SURV, name_dir, "sceptre", "grna_matrix_aligned.rds"))
  ids  <- rownames(grna)
  is_nt <- Reduce(`|`, lapply(nt_patterns, function(p) grepl(p, ids)))
  target <- sub(regex, "\\1", ids)
  if (upper) target <- toupper(target)
  target <- ifelse(is_nt, NT, target)
  matched <- grepl(regex, ids) | is_nt
  target[!matched] <- NA_character_
  data.frame(grna_id = ids, grna_target = target,
             grna_group = NA_character_, stringsAsFactors = FALSE)
}

build_targets_dctap <- function(subdir) {
  dir <- file.path(SURV, subdir)
  gf <- read.csv(file.path(dir, "guide_features.csv"), stringsAsFactors = FALSE)
  # NT controls
  is_nt <- grepl("^(safe|nontarget|NC)[_-]", gf$id, ignore.case = TRUE)
  # Targeting guides: extract gene from prefix.  Patterns: MYC-TSS-*, MYC-e2-*,
  # GATA1-TSS-*, HDAC6-TSS-*, e-GATA1-*, e-HDAC6-*.
  target <- gf$id
  target[grepl("^MYC[-_]",   target)] <- "MYC"
  target[grepl("^GATA1[-_]", target)] <- "GATA1"
  target[grepl("^HDAC6[-_]", target)] <- "HDAC6"
  target[grepl("^e-GATA1",   target)] <- "GATA1"
  target[grepl("^e-HDAC6",   target)] <- "HDAC6"
  target[is_nt] <- NT
  unmatched <- !is_nt & !(target %in% c("MYC", "GATA1", "HDAC6"))
  if (any(unmatched)) { cat("  dctap unmatched (set NA):", paste(gf$id[unmatched], collapse=", "), "\n")
    target[unmatched] <- NA_character_ }
  data.frame(grna_id = gf$id, grna_target = target,
             grna_group = NA_character_, stringsAsFactors = FALSE)
}

# ---- run + verify -----------------------------------------------------------
verify_and_save <- function(tdf, name_dir) {
  dir <- file.path(SURV, name_dir)
  resp <- readRDS(file.path(dir, "sceptre", "response_matrix.rds"))
  grna <- readRDS(file.path(dir, "sceptre", "grna_matrix_aligned.rds"))
  # Make sure every grna in tdf is also in grna_matrix
  in_grna <- tdf$grna_id %in% rownames(grna)
  if (!all(in_grna)) {
    cat(sprintf("  WARN: %d/%d guides in target_df not in grna_matrix (dropping)\n",
                sum(!in_grna), nrow(tdf)))
    tdf <- tdf[in_grna, ]
  }
  # Resolve SYMBOL -> GEX rowname (ENSG).  We keep BOTH columns: grna_target
  # is sceptre's contract (any string), grna_target_symbol is human-readable.
  sym2id <- build_symbol_to_id(name_dir)
  cat(sprintf("  symbol->id map: %d entries (first few: %s)\n",
              length(sym2id), paste(head(names(sym2id), 3), collapse=", ")))
  tdf$grna_target_symbol <- tdf$grna_target
  tdf$grna_target <- resolve_to_gex_rownames(tdf$grna_target, sym2id, rownames(resp))
  # diagnostics
  nt_n <- sum(tdf$grna_target == NT, na.rm = TRUE)
  targeting_n <- sum(tdf$grna_target != NT & !is.na(tdf$grna_target))
  na_n <- sum(is.na(tdf$grna_target))
  ut <- unique(tdf$grna_target[tdf$grna_target != NT & !is.na(tdf$grna_target)])
  cat(sprintf("  %d guides | %d NT | %d targeting -> %d distinct ENSG resolved | %d unresolved\n",
              nrow(tdf), nt_n, targeting_n, length(ut), na_n))
  if (na_n > 0) {
    nas <- unique(tdf$grna_target_symbol[is.na(tdf$grna_target)])
    cat(sprintf("    unresolved symbols (sample): %s\n",
                paste(head(nas, 10), collapse = ", ")))
  }
  # Save with both columns (lab convention uses grna_target as the canonical)
  write.csv(tdf, file.path(dir, "sceptre", "grna_target_data_frame.csv"),
            row.names = FALSE, quote = TRUE)
  cat(sprintf("  saved %s/sceptre/grna_target_data_frame.csv\n", name_dir))
  invisible(tdf)
}

# ---- run all configured -----------------------------------------------------
main <- function() {
  configs <- list(
    a549 = list(dir = "a549_crispri_perturbseq_GSE236304", fn = build_targets_a549),
    cd8_tcell = list(dir = "cd8_tcell_perturbcite_GSE279498",
      fn = function() build_targets_simple("cd8_tcell_perturbcite_GSE279498",
              regex = "^g(.+)_\\d+$",
              nt_patterns = c("(?i)nontarget", "(?i)^gScramble"), upper = FALSE)),
    dctap_highmoi = list(dir = "dctap_k562_highmoi_GSE303901",
      fn = function() build_targets_dctap("dctap_k562_highmoi_GSE303901")),
    dctap_lowmoi  = list(dir = "dctap_k562_lowmoi_GSE303901",
      fn = function() build_targets_dctap("dctap_k562_lowmoi_GSE303901")),
    endoc = list(dir = "endoc_t2d_perturbseq_GSE273677",
      fn = function() build_targets_simple("endoc_t2d_perturbseq_GSE273677",
              regex = "^([^-]+)-\\d+$", nt_patterns = c("^Non-?Targeting"), upper = TRUE)),
    # invivo_cortex: mouse target symbols (title case); SKIP until GEX is extracted
    # (R 4.4 + qs v1 needed; see DE_DATASETS_STATUS.md).
    # invivo_cortex = list(dir = "invivo_cortex_perturbseq_GSE249416",
    #   fn = function() build_targets_simple("invivo_cortex_perturbseq_GSE249416",
    #           regex = "^([^-]+)-\\d+$",
    #           nt_patterns = c("^Neg-backbone$", "^NonTarget[12]$", "^SafeTarget$", "^GFP$"),
    #           upper = FALSE)),
    tcell_cd4 = list(dir = "tcell_cd4_perturbseq_GSE314342",
      fn = function() build_targets_simple("tcell_cd4_perturbseq_GSE314342",
              regex = "^([^-]+)-\\d+$", nt_patterns = c("^___NEVER___$"), upper = TRUE))
  )
  for (nm in names(configs)) {
    cat(sprintf("\n=== %s ===\n", nm))
    tdf <- tryCatch(configs[[nm]]$fn(), error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(tdf)) verify_and_save(tdf, configs[[nm]]$dir)
  }
  cat("\n=== multiome_erythroid: SKIPPED -- guide names are 20-mer sequences;\n")
  cat("    target mapping requires the Liu et al 2024 paper's supplementary table.\n")
}
if (sys.nframe() == 0) main()
