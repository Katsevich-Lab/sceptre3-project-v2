#!/usr/bin/env Rscript
# =============================================================================
# REALISM GATE  --  do the simulators look like real data?
# Computes the method-agnostic realism summaries for:
#   - REAL barnyard (4 samples x 2 purities), ambient mask = wrong-species
#   - Model A and Model B (all datasets), ambient mask = non-integrated (!Z)
#   - the lab's OLD comprehensive sim, ambient mask = non-perturbed (!true_pert)
# and writes a comparison table + figure.  The headline: A/B land in the real
# bands; the old comprehensive sim's ambient is far lumpier than any real data.
# Outputs: results/sim_framework/realism_gate.{csv,png}
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(dplyr); library(tidyr) })
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
OUT <- SIMFW(); amask <- function(Z) as(!as.matrix(as(Z, "lgCMatrix")), "lgCMatrix")

rows <- list()

# ---- REAL barnyard, both purities ------------------------------------------
for (s in BARN_SAMPLES) for (pur in c(0.90, 0.99)) {
  b <- load_barnyard(s, purity = pur)
  r <- realism_stats(b$counts, ambient_mask = b$ambient_mask,
                     label = sprintf("%s pur%02d", s, round(pur * 100)))
  r$source <- "real barnyard"; r$chemistry <- b$chemistry; r$purity <- pur
  rows[[length(rows) + 1]] <- r
  cat("real:", s, pur, "\n")
}

# ---- Model A and Model B ----------------------------------------------------
for (mf in c("manifest_modelA.csv", "manifest_modelB.csv")) {
  man <- read.csv(file.path(OUT, mf))
  for (i in seq_len(nrow(man))) {
    d <- load_dataset(man$id[i])
    r <- realism_stats(d$counts, ambient_mask = amask(d$Z), label = man$id[i])
    r$source <- if (man$model[i] == "A") "Model A (semi-synthetic)" else "Model B (mechanistic)"
    r$chemistry <- man$chemistry[i]; r$purity <- man$purity[i]
    rows[[length(rows) + 1]] <- r
  }
  cat("done", mf, "\n")
}

# ---- lab's OLD comprehensive sim (the cautionary baseline) ------------------
csim <- file.path(path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data"),
                  "sims_comprehensive")
if (file.exists(file.path(csim, "sceptre", "grna_matrix.rds"))) {
  gm <- readRDS(file.path(csim, "sceptre", "grna_matrix.rds"))
  tp <- readRDS(file.path(csim, "true_pert_matrix.rds"))
  set.seed(1); gs <- sort(sample(nrow(gm), 100)); cspl <- sort(sample(ncol(gm), 8000))
  gm <- gm[gs, cspl]; tp <- tp[gs, cspl]
  r <- realism_stats(gm, ambient_mask = as(!as.matrix(as(tp, "lgCMatrix")), "lgCMatrix"),
                     label = "lab comprehensive sim")
  r$source <- "lab comprehensive sim"; r$chemistry <- "(sum-process)"; r$purity <- NA
  rows[[length(rows) + 1]] <- r
  cat("done lab comprehensive sim\n")
}

tab <- bind_rows(rows)
keep <- c("source","chemistry","purity","label","n_cells","moi_med","depth_ratio_med","depth_cor",
          "comp_gini","mode_gap_med","ambient_fano","ambient_vmr_nz","ambient_max")
write.csv(tab[, keep], file.path(OUT, "realism_gate.csv"), row.names = FALSE)

# ---- summary by source x chemistry -----------------------------------------
summ <- tab %>% group_by(source, chemistry) %>%
  summarise(n = n(),
            depth_ratio = round(median(depth_ratio_med, na.rm = TRUE), 1),
            gini = round(median(comp_gini, na.rm = TRUE), 2),
            mode_gap = round(median(mode_gap_med, na.rm = TRUE), 2),
            ambient_vmr_nz = round(median(ambient_vmr_nz, na.rm = TRUE), 1),
            ambient_max = round(median(ambient_max, na.rm = TRUE)), .groups = "drop")
cat("\n=== realism summary (median over datasets within source x chemistry) ===\n")
print(as.data.frame(summ), row.names = FALSE)

# ---- figure: ambient lumpiness (the headline) ------------------------------
pd <- tab %>% filter(!is.na(ambient_vmr_nz)) %>%
  mutate(grp = ifelse(source == "lab comprehensive sim", "lab comprehensive sim", source))
p <- ggplot(pd, aes(x = source, y = ambient_vmr_nz, colour = chemistry)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  scale_y_log10() + coord_flip() +
  labs(title = "Ambient lumpiness: var/mean of nonzero ambient counts (dashed = Poisson)",
       subtitle = "Model A inherits real shape; Model B is clean Poisson; the lab sum-process sim is far lumpier",
       y = "var / mean of nonzero ambient counts (log scale)", x = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "realism_gate.png"), p, width = 10, height = 5.5, dpi = 140)
cat("\nwrote realism_gate.{csv,png}\n")
