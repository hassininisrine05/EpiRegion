# Reproduce the benchmark figures: Rscript scripts/run_benchmark.R
library(hifiepi); library(ggplot2)
res <- run_benchmark(coverages = c(10, 20, 30), deltas = c(0.05, 0.10, 0.15),
                     n_rep = 10, seed = 2026)
dir.create("docs", showWarnings = FALSE)
write.csv(res, "docs/benchmark_results.csv", row.names = FALSE)
for (m in c("tpr", "fdr", "auc"))
  ggsave(sprintf("docs/benchmark_%s.png", m), plot_benchmark(res, m),
         width = 8, height = 3.6, dpi = 150)
print(aggregate(value ~ coverage + delta + method + metric, res, mean))
