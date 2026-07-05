#!/usr/bin/env Rscript
# How much do Otsu and the smoothed valley actually disagree? Their disagreement
# on a gRNA is exactly the cells whose count lies between the two thresholds. We
# tally this across ALL gRNAs and report agreement as a fraction of the cells
# either method assigns (1 - Jaccard), since raw cell counts make it look ~100%.
#
# Usage: Rscript agreement.R [dataset_id]

suppressPackageStartupMessages(library(Matrix))
HERE     <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE))))
DATA_DIR <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(HERE, "..", "guide-assignment-pipeline", "bin", "script", "lib",
                 "threshold_methods.R"))

dataset_id <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "replogle-rd7_medium"
max_grna   <- if (length(commandArgs(TRUE)) >= 2) as.integer(commandArgs(TRUE)[2]) else Inf
gm <- as(readRDS(file.path(DATA_DIR, dataset_id, "sceptre", "grna_matrix.rds")),
         "RsparseMatrix")
n_cell <- ncol(gm)
if (nrow(gm) > max_grna) {            # subsample gRNAs (keeps memory bounded)
  set.seed(1); gm <- gm[sort(sample(nrow(gm), max_grna)), , drop = FALSE]
  cat(sprintf("(subsampled to %d gRNAs)\n", max_grna))
}
n_grna <- nrow(gm)
cat(sprintf("Dataset %s: %d gRNAs x %d cells\n", dataset_id, n_grna, n_cell))

to <- tv <- n_o <- n_v <- n_inter <- n_union <- n_disagree <- numeric(n_grna)
for (g in seq_len(n_grna)) {
  counts <- as.numeric(gm[g, ])
  t_o <- otsu_threshold_log1p(counts)$t
  t_v <- smoothed_valley_threshold(counts)$t
  a_o <- if (is.finite(t_o)) counts >= t_o else logical(length(counts))
  a_v <- if (is.finite(t_v)) counts >= t_v else logical(length(counts))
  to[g] <- t_o; tv[g] <- t_v
  n_o[g] <- sum(a_o); n_v[g] <- sum(a_v)
  n_inter[g]   <- sum(a_o & a_v)
  n_union[g]   <- sum(a_o | a_v)
  n_disagree[g]<- sum(xor(a_o, a_v))
}

jac    <- ifelse(n_union > 0, n_inter / n_union, NA)   # per-gRNA agreement
both0  <- n_union == 0
cat("\n--- threshold values ---\n")
cat(sprintf("  median Otsu t = %.0f, median valley t = %.0f\n",
            median(to[is.finite(to)]), median(tv[is.finite(tv)])))
cat(sprintf("  gRNAs where t_otsu == t_valley: %d / %d (%.1f%%)\n",
            sum(to == tv & is.finite(to)), n_grna,
            100 * mean(to == tv & is.finite(to))))
cat(sprintf("  Otsu abstains: %d | valley abstains: %d\n",
            sum(!is.finite(to)), sum(!is.finite(tv))))

cat("\n--- agreement on assigned cells (over gRNAs with any assignment) ---\n")
cat(sprintf("  total cells assigned: Otsu = %d, valley = %d\n", sum(n_o), sum(n_v)))
cat(sprintf("  pooled Jaccard (sum inter / sum union) = %.4f\n",
            sum(n_inter) / sum(n_union)))
cat(sprintf("  pooled disagreement = %d / %d union assignments = %.2f%%\n",
            sum(n_disagree), sum(n_union), 100 * sum(n_disagree) / sum(n_union)))
cat(sprintf("  mean per-gRNA Jaccard   = %.4f\n", mean(jac[!both0])))
cat(sprintf("  median per-gRNA Jaccard = %.4f\n", median(jac[!both0])))
cat(sprintf("  gRNAs with perfect agreement (Jaccard=1): %d / %d (%.1f%%)\n",
            sum(jac[!both0] == 1), sum(!both0), 100 * mean(jac[!both0] == 1)))

cat("\n--- disagreement relative to ALL cell-gRNA pairs ---\n")
cat(sprintf("  %d disagreeing pairs / %.3g total pairs = %.5f%%\n",
            sum(n_disagree), n_grna * n_cell,
            100 * sum(n_disagree) / (n_grna * n_cell)))

# distribution of per-gRNA disagreement counts
cat("\n--- per-gRNA disagreement (cells) quantiles ---\n")
print(quantile(n_disagree, c(0, .5, .9, .99, 1)))
saveRDS(data.frame(grna = rownames(gm), t_otsu = to, t_valley = tv,
                   n_otsu = n_o, n_valley = n_v, n_disagree = n_disagree,
                   jaccard = jac),
        file.path(HERE, "results", paste0("agreement_", dataset_id, ".rds")))
