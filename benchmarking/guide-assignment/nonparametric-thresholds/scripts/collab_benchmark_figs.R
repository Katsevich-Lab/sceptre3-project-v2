## Benchmark figures for the collaborator writeup:
##  (1) barnyard ground-truth Table-2 accuracy: fishash+ vs fishash vs published methods
##      (fishash+ from scripts/barnyard_benchmark_fishashplus.R; published from Final_Barnyard_ranking.csv,
##       whose fishash number = our exact reproduction, so the scales are identical).
##  (2) high-MOI simulation: recall vs MOI for fishash vs fishash+ (from benchmark_update/canonical_bench.csv).
suppressPackageStartupMessages({library(ggplot2)})
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- (1) barnyard bar chart -----------------------------------------------------
bp  <- read.csv(file.path(OUT, "barnyard_benchmark_fishashplus.csv"), check.names = FALSE)
fap <- bp$`fishash+`[bp$sample == "MEAN"]
fa  <- bp$fishash[bp$sample == "MEAN"]
rk  <- read.csv("results/Final_Barnyard_ranking.csv")   # published methods, same Table-2 metric
pub <- data.frame(
  method = c("fishash+", "fishash", "SCEPTRE mixture", "CLEANSER", "crispat (Poisson)", "geomux"),
  acc    = c(fap, fa,
             rk$mean_acc[rk$method == "sceptre mixture"],
             rk$mean_acc[rk$method == "cleanser cs"],
             rk$mean_acc[rk$method == "crispat poisson"],
             rk$mean_acc[rk$method == "geomux"]),
  ours   = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
pub$method <- factor(pub$method, levels = pub$method[order(pub$acc)])
pB <- ggplot(pub, aes(method, acc, fill = ours)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.3f", acc)), hjust = -0.15, size = 3.5) +
  coord_flip(ylim = c(0.75, 0.985)) +
  scale_fill_manual(values = c(`FALSE` = "grey72", `TRUE` = "#D55E00"), guide = "none") +
  labs(title = "Barnyard ground-truth accuracy (per-cell correct-species, 4 samples)",
       subtitle = "fishash+ matches fishash (no regression); both lead published methods",
       x = NULL, y = "assignment accuracy") +
  theme_bw(base_size = 12) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "barnyard_benchmark.png"), pB, width = 8.5, height = 3.8, dpi = 130)

## ---- (2) high-MOI simulation ----------------------------------------------------
cb <- read.csv("results/benchmark_update/canonical_bench.csv")
cb <- cb[cb$regime %in% c("MOI 0.3","MOI 1","MOI 3","MOI 5") & cb$method %in% c("fishash","+ depth"), ]
cb$moi <- as.numeric(sub("MOI ", "", cb$regime))
ag <- aggregate(cbind(recall, fdr) ~ moi + method, cb, mean)
ag$method <- factor(ifelse(ag$method == "+ depth", "fishash+", "fishash"), levels = c("fishash","fishash+"))

lab <- ag[ag$moi == 5, ]
pS <- ggplot(ag, aes(moi, recall, colour = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  geom_text(data = lab, aes(label = sprintf("%.2f", recall)), hjust = -0.3, size = 3.6, show.legend = FALSE) +
  scale_colour_manual(values = c("fishash" = "#0072B2", "fishash+" = "#D55E00")) +
  scale_x_continuous(breaks = c(0.3, 1, 3, 5), limits = c(0.2, 5.6)) +
  labs(title = "Recall of true integrations vs multiplicity of infection (simulation)",
       subtitle = "As MOI rises, fishash+ recovers the recall fishash loses to self-masking (both control FDR < 0.02)",
       x = "MOI (mean guides per cell)", y = "recall of true (guide, cell) integrations", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = c(0.14, 0.22),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))
ggsave(file.path(OUT, "sim_highmoi_recall.png"), pS, width = 8.5, height = 4.6, dpi = 130)

cat("barnyard: fishash+ =", fap, " fishash =", fa, "\n")
cat("sim recall @ MOI5: fishash+ =", ag$recall[ag$moi==5 & ag$method=="fishash+"],
    " fishash =", ag$recall[ag$moi==5 & ag$method=="fishash"], "\n")
cat("wrote barnyard_benchmark.png + sim_highmoi_recall.png\n")
