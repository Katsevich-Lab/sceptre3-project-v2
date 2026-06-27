#!/usr/bin/env Rscript
# Show that each regime captures the WITHIN-dataset variety of guide types: the
# distribution of per-guide ABOVE-VALLEY MEANS (mu), real vs simulated, for every
# regime.  mu is what the per-guide fit consumes and the sim simulates from, so
# pairing real-mu to sim-mu is the apples-to-apples check; the sim spread should
# overlap the real spread (the regime then contains the real mix of strong/weak
# guides).
# Output: results/sim_framework/within_regime_variation.png
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
OUT <- SIMFW()
pg <- readRDS(file.path(OUT, "real_perguide.rds"))           # per-guide df: mu, theta, signal_frac, right_mode
ch <- read.csv(file.path(OUT, "real_characterization.csv"))

# per-guide above-valley MEAN of a simulated matrix (the SAME quantity sim_characterize
# saves as `mu` for the real data; consistent comparison).
sim_mus <- function(counts) {
  smr <- as(counts, "RsparseMatrix"); m <- c()
  for (g in seq_len(nrow(smr))) { cv <- as.numeric(smr[g, ]); if (sum(cv > 0) < 30) next
    v <- smoothed_valley_threshold(cv)
    if (isTRUE(v$ok)) { sig <- cv[cv >= v$t]; if (length(sig)) m <- c(m, mean(sig)) } }
  m
}

df <- list()
for (ds in names(pg)) {
  real_mu <- pg[[ds]]$mu; if (length(real_mu) < 3) next
  sim <- tryCatch(load_dataset(paste0("B__regime__", ds)), error = function(e) NULL); if (is.null(sim)) next
  sim_mu <- sim_mus(sim$counts)
  df[[ds]] <- rbind(data.frame(dataset = ds, source = "real",      mu = real_mu),
                    data.frame(dataset = ds, source = "simulated", mu = sim_mu))
}
D <- bind_rows(df) %>% filter(mu > 0)
ord <- ch %>% arrange(separation) %>% pull(dataset)
ord <- ord[ord %in% unique(D$dataset)]
labs <- setNames(sprintf("%s (sep %.1f)", ord, ch$separation[match(ord, ch$dataset)]), ord)
D$dataset <- factor(D$dataset, levels = ord)

p <- ggplot(D, aes(mu, fill = source, colour = source)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~dataset, scales = "free", ncol = 4, labeller = as_labeller(labs)) +
  scale_x_log10() +
  scale_fill_manual(values = c(real = "#C9897B", simulated = "#1D9E75"), name = NULL) +
  scale_colour_manual(values = c(real = "#C9897B", simulated = "#1D9E75"), name = NULL) +
  labs(title = sprintf("Within-regime guide variety: per-guide above-valley MEAN (mu), real vs simulated (%d regimes)",
                       length(unique(D$dataset))),
       subtitle = "sim draws (mu, theta, signal_frac) per guide from the dataset's empirical fits -- the sim spread (green) should overlap the real spread (orange)",
       x = "per-guide above-valley mean (UMI, log scale)", y = "density") +
  theme_bw(base_size = 8) +
  theme(legend.position = "top", strip.text = element_text(size = 7),
        axis.text.x = element_text(size = 6))
ggsave(file.path(OUT, "within_regime_variation.png"), p, width = 12, height = 10, dpi = 140)
cat(sprintf("wrote within_regime_variation.png (%d datasets)\n", length(unique(D$dataset))))
