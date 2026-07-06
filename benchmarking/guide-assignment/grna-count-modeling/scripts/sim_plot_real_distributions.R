#!/usr/bin/env Rscript
# =============================================================================
# Gallery of sample gRNA UMI count distributions across cells, for all 17 real
# datasets (the parameter space made tangible).  For each dataset we pick a
# representative bimodal guide (separation closest to the dataset median) and
# plot its count distribution across cells, coloured ambient (below valley) vs
# signal (above valley).  Panels ordered hard -> easy by mode separation.
# Output: results/sim_framework/real_distributions.png
# =============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
SURV <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
DATA <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
REPRO<- file.path(HERE, "external", "repro_work")
OUT  <- SIMFW()
ch <- read.csv(file.path(OUT, "real_characterization.csv"))

loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) { f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) }) }
loaders[["gasperini"]] <- function() readRDS(file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
loaders[["replogle"]]  <- function() readRDS(file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { suppressPackageStartupMessages(library(fishash)); data(crispat_schraivogel); crispat_schraivogel }
for (s in BARN_SAMPLES) { mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (file.exists(mtx)) loaders[[paste0("barnyard_", sub("mix","",s))]] <- local({ mm <- mtx
    function() as(readMM(mm), "CsparseMatrix") }) }

# representative guide = bimodal guide whose separation is closest to the dataset median
pick_repr <- function(gm, max_cells = 30000, seed = 1) {
  set.seed(seed)                                      # seed BOTH samples (cells + guides)
  gm <- as(gm, "CsparseMatrix"); if (ncol(gm) > max_cells) gm <- gm[, sort(sample(ncol(gm), max_cells))]
  gmr <- as(gm, "RsparseMatrix"); gi <- if (nrow(gmr) > 400) sort(sample(nrow(gmr), 400)) else seq_len(nrow(gmr))
  cand <- list(); seps <- c()
  for (g in gi) { cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < 50) next
    v <- smoothed_valley_threshold(cv); if (isTRUE(v$ok) && sum(cv >= v$t) < 0.3 * length(cv)) {
      cand[[length(cand)+1]] <- list(cv = cv, t = v$t); seps <- c(seps, log1p(v$mode2) - log1p(v$mode1)) } }
  if (!length(cand)) { cv <- as.numeric(gmr[which.max(Matrix::rowSums(gmr > 0)), ]); return(list(cv = cv, t = Inf)) }
  cand[[which.min(abs(seps - median(seps)))]]
}

df <- list()
for (nm in names(loaders)) {
  r <- tryCatch(pick_repr(loaders[[nm]]()), error = function(e) NULL); if (is.null(r)) next
  sep <- ch$separation[ch$dataset == nm]; sep <- if (length(sep)) sep[1] else NA
  df[[nm]] <- data.frame(dataset = nm, sep = sep, count = r$cv,
                         mode = ifelse(r$cv >= r$t, "signal", "ambient"))
}
D <- bind_rows(df)
# order facets hard -> easy by separation; label with separation
ord <- ch %>% arrange(separation) %>% pull(dataset)
D$dataset <- factor(D$dataset, levels = ord)
labs <- setNames(sprintf("%s\n(sep %.1f)", ord, ch$separation[match(ord, ch$dataset)]), ord)
D$z <- log1p(D$count)

brks <- c(0,1,2,5,10,20,50,100,200,500,1000,2000)
p <- ggplot(D, aes(z, fill = mode)) +
  geom_histogram(bins = 38, position = "stack") +
  facet_wrap(~dataset, scales = "free", ncol = 4, labeller = as_labeller(labs)) +
  scale_y_log10() +
  scale_x_continuous(breaks = log1p(brks), labels = brks) +
  scale_fill_manual(values = c(ambient = "#C9ADA7", signal = "#1D9E75"), name = NULL) +
  labs(title = "Sample gRNA UMI count distributions across cells — all 17 real datasets",
       subtitle = "one representative guide per dataset; coloured ambient (below valley) vs signal (above); ordered hard -> easy",
       x = "UMI count (log1p axis)", y = "number of cells (log)") +
  theme_bw(base_size = 8) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        strip.text = element_text(size = 7))
ggsave(file.path(OUT, "real_distributions.png"), p, width = 12, height = 11, dpi = 140)
cat("wrote real_distributions.png (", length(unique(D$dataset)), "datasets )\n")
