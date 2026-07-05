#!/usr/bin/env Rscript
# Build a COMPREHENSIVE simulation spanning the real gRNA parameter ranges found
# across diverse 2024-26 datasets (characterize_grna.R):
#   mu_pert in {15,50,150,500}  (signal level; real span ~15-600)
#   theta_pert in {0.5,1.5}     (NB dispersion; real span ~0.4-2)
#   2 MOI regimes via REAL library sizes: low-MOI (Replogle-like) and high-MOI
#   (Gasperini-like). prob_pert ~0.02 so each guide has enough perturbed cells.
# Uses the paper's sum-process simulator. Saves grna_matrix.rds + true_pert_matrix.rds.

suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
source(file.path(GA, "grna-simulator", "sims-for-paper.R"))   # grna_simulator_iid_sum_process

NCELL <- 50000
get_lib <- function(ds) {
  gm <- readRDS(file.path(DATA, ds, "sceptre", "grna_matrix.rds"))
  set.seed(1); cs <- sort(sample(ncol(gm), min(NCELL, ncol(gm))))
  nnz <- Matrix::colSums(gm[, cs] > 0); nnz[nnz == 0] <- 1
  nnz / mean(nnz)
}
libs <- list(lowMOI = get_lib("replogle-rd7_medium"), highMOI = get_lib("gasperini"))
cat("library sizes -- lowMOI mean nnz:", round(mean(Matrix::colSums(readRDS(file.path(DATA,"replogle-rd7_medium/sceptre/grna_matrix.rds"))>0)),1),
    " highMOI:", round(mean(Matrix::colSums(readRDS(file.path(DATA,"gasperini/sceptre/grna_matrix.rds"))>0)),1), "\n")

bg <- list(prob_endo = 0.02, mu_endo = 3, theta_endo = 2, mu_exo = 0.07)  # realistic small ambient
mu_grid    <- c(5, 15, 50, 150, 500)   # 5 = CD4/CD8-like extreme low signal (real)
theta_grid <- c(0.15, 0.5, 1.5)        # 0.15 = iPSC-like extreme overdispersion (real theta~0.07-1.6)
NG <- 60   # guides per (moi, mu, theta) cell

pieces_g <- list(); pieces_t <- list(); seed <- 100
for (moi in names(libs)) for (mu in mu_grid) for (th in theta_grid) {
  seed <- seed + 1
  params <- c(bg, list(prob_pert = 0.02, mu_pert = mu, theta_pert = th))
  s <- grna_simulator_iid_sum_process(NG, libs[[moi]], params, seed = seed)
  tag <- sprintf("%s.mu%d.th%s", moi, mu, sub("\\.", "", as.character(th)))
  gn  <- paste0("sim_", moi, "_mu", mu, "_th", gsub("\\.","",th), "_", seq_len(NG))
  pieces_g[[tag]] <- s$grna_matrix_sim
  pieces_t[[tag]] <- s$true_perts
  attr(pieces_g[[tag]], "gn") <- gn
}
gns <- unlist(lapply(pieces_g, attr, "gn"))
cell_names <- paste0("cell_idx_", seq_len(NCELL))
grna_matrix      <- as(do.call(rbind, pieces_g), "CsparseMatrix")
true_pert_matrix <- as(do.call(rbind, pieces_t), "CsparseMatrix")
dimnames(grna_matrix)      <- list(gns, cell_names)
dimnames(true_pert_matrix) <- list(gns, cell_names)

out <- file.path(DATA, "sims_comprehensive"); dir.create(file.path(out, "sceptre"), recursive = TRUE, showWarnings = FALSE)
saveRDS(grna_matrix,      file.path(out, "sceptre", "grna_matrix.rds"))
saveRDS(true_pert_matrix, file.path(out, "true_pert_matrix.rds"))
saveRDS(list(bg = bg, mu_grid = mu_grid, theta_grid = theta_grid, moi = names(libs),
             prob_pert = 0.02, ng_per = NG, num_cells = NCELL),
        file.path(out, "sim_params.rds"))
write.csv(data.frame(grna_id = gns, grna_target = c("non-targeting", paste0("t_", seq_len(length(gns)-1)))),
          file.path(out, "sceptre", "grna_target_data_frame.csv"), row.names = FALSE)
cat("Saved sims_comprehensive:", nrow(grna_matrix), "guides x", ncol(grna_matrix), "cells (",
    length(mu_grid)*length(theta_grid)*length(libs), "param groups )\n")
