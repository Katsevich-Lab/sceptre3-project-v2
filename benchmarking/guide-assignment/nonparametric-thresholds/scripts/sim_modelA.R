#!/usr/bin/env Rscript
# =============================================================================
# MODEL A  --  semi-synthetic: REAL barnyard ambient + INJECTED synthetic signal
#
# The most faithful object we can build: the ambient is literally real (the
# wrong-species submatrix of a QC'd barnyard sample, which is provably ambient),
# so its shape -- Poisson in CROP-seq, doublet-heavy in direct-capture -- comes
# along for free, with no parametric ambient model to mis-specify.  Only the
# signal is synthetic, and therefore the ground truth is exact.
#
# Construction, per (sample, purity, mu_pert, MOI):
#   1. Load barnyard at the chosen purity (0.90 = authors' QC; 0.99 = strict).
#   2. Focal cells = the larger pure-species group.  Their WRONG-species guide
#      counts form the real ambient backbone B (100 guides x N cells).
#   3. Integrations Z: n_c ~ Poisson(MOI) guides per cell, drawn (without
#      replacement) from a mild-Dirichlet guide-popularity vector.
#   4. Signal: at each integrated (g,c), add NB(mu = mu_pert * r_c, size = theta),
#      r_c = cell capture size factor (real library / median).  Some draws are 0
#      => realistic "integrated but not observed" dropout.
#   5. Observed Y = B + signal.  Non-integrated entries are pure real ambient.
#
# mu_pert grids are chemistry-specific so the moderate level lands near the real
# separation; the realism gate (sim_realism_gate.R) confirms this.
# Outputs: results/sim_framework/datasets/A__*  +  manifest_modelA.csv
# =============================================================================
suppressPackageStartupMessages({ library(Matrix) })
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))

build_modelA <- function(sample, purity, mu_pert, moi, theta_sig = 1,
                         max_cells = 8000, seed = 1) {
  b <- load_barnyard(sample, purity = purity)
  # Focal species = the one whose wrong-species ambient carries the HEAVIER tail
  # (so the real doublet tail at low purity is preserved, not averaged away),
  # unless that group is too small, in which case use the larger group.
  hum_c <- which(b$cell_homo); mus_c <- which(!b$cell_homo)
  max_h <- if (length(mus_c)) max(b$counts[which(!b$guide_homo), hum_c, drop = FALSE]) else 0  # human cells, mouse-guide ambient
  max_m <- if (length(hum_c)) max(b$counts[which(b$guide_homo),  mus_c, drop = FALSE]) else 0  # mouse cells, human-guide ambient
  focal_human <- if (min(length(hum_c), length(mus_c)) < 1000) length(hum_c) >= length(mus_c) else max_h >= max_m
  cell_sel   <- if (focal_human) hum_c else mus_c
  amb_guides <- if (focal_human) which(!b$guide_homo) else which(b$guide_homo)  # wrong-species => ambient
  B    <- b$counts[amb_guides, cell_sel, drop = FALSE]                # real ambient backbone
  capt <- Matrix::colSums(b$counts[, cell_sel, drop = FALSE])          # capture proxy (full real library)
  set.seed(seed)
  if (ncol(B) > max_cells) { keep <- sort(sample(ncol(B), max_cells)); B <- B[, keep]; capt <- capt[keep] }
  G <- nrow(B); C <- ncol(B)
  rownames(B) <- paste0("g", seq_len(G)); colnames(B) <- paste0("c", seq_len(C))
  r_c <- capt / median(capt[capt > 0]); r_c[!is.finite(r_c) | r_c <= 0] <- min(r_c[r_c > 0])
  p_g <- rdirichlet1(rep(3, G))                                        # mild guide-popularity unevenness
  n_c <- rpois(C, moi)                                                 # integrations per cell
  zi <- integer(0); zj <- integer(0)
  for (cc in which(n_c > 0)) {
    k <- min(n_c[cc], G)
    g_hit <- sample.int(G, size = k, replace = FALSE, prob = p_g)
    zi <- c(zi, g_hit); zj <- c(zj, rep.int(cc, k))
  }
  Z <- sparseMatrix(i = zi, j = zj, x = TRUE, dims = c(G, C), dimnames = dimnames(B))
  sig <- rnbinom(length(zi), mu = pmax(mu_pert * r_c[zj], 1e-6), size = theta_sig)
  S <- sparseMatrix(i = zi, j = zj, x = as.numeric(sig), dims = c(G, C), dimnames = dimnames(B))
  Y <- drop0(B + S)
  list(counts = as(Y, "CsparseMatrix"), Z = as(Z, "lgCMatrix"),
       n_guides = G, n_cells = C, focal = if (focal_human) "human" else "mouse",
       ambient_nnz = sum(B > 0), signal_nnz = length(zi))
}

# ---- grid -------------------------------------------------------------------
mu_grids <- list("CROP-seq"       = c(weak = 5,  mod = 15,  strong = 50),
                 "direct-capture" = c(weak = 30, mod = 100, strong = 400))
moi_grid <- c(low = 1, high = 5, vhigh = 15)   # vhigh = genuinely high MOI (Gasperini-like)
purities <- c(0.90, 0.99)
THETA    <- 1

dir.create(DSDIR(), recursive = TRUE, showWarnings = FALSE)
# Two passes so the original low/high datasets keep their seeds/ids byte-for-byte
# (pass 1 == the original loop order) and vhigh is appended (pass 2).
rows <- list(); k <- 0
for (pass_moi in list(c("low","high"), c("vhigh"))) {
for (s in BARN_SAMPLES) {
  chem <- if (grepl("Cropseq", s)) "CROP-seq" else "direct-capture"
  mus  <- mu_grids[[chem]]
  for (pur in purities) for (mn in names(mus)) for (moin in pass_moi) {
    k <- k + 1
    mu  <- mus[[mn]]; moi <- moi_grid[[moin]]
    id  <- sprintf("A__%s__pur%02d__mu%s__moi%s", s, round(pur * 100), mn, moin)
    seed <- 1000 + k
    d <- build_modelA(s, pur, mu, moi, theta_sig = THETA, seed = seed)
    meta <- list(model = "A", id = id, sample = s, chemistry = chem, purity = pur,
                 mu_level = mn, mu_pert = mu, moi_level = moin, moi = moi, theta_sig = THETA,
                 n_guides = d$n_guides, n_cells = d$n_cells, focal = d$focal, seed = seed)
    save_dataset(id, d$counts, d$Z, meta)
    rows[[id]] <- as.data.frame(meta[c("model","id","sample","chemistry","purity","mu_level","mu_pert",
                                       "moi_level","moi","theta_sig","n_guides","n_cells","focal")])
    cat(sprintf("[%2d] %-46s  G=%d C=%d  amb_nnz=%d sig_nnz=%d\n",
                k, id, d$n_guides, d$n_cells, d$ambient_nnz, d$signal_nnz))
  }
}
}
man <- do.call(rbind, rows); rownames(man) <- NULL
write.csv(man, file.path(SIMFW(), "manifest_modelA.csv"), row.names = FALSE)
cat("\nwrote", nrow(man), "Model A datasets +", file.path(SIMFW(), "manifest_modelA.csv"), "\n")
