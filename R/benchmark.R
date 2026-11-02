#' Benchmark CpG-level vs region-level strategies
#'
#' For each coverage x effect-size setting, simulates `n_rep` datasets and
#' records (i) region-level sensitivity / FDR of a CpG-level strategy (region
#' called if any CpG is significant at BH 5%) vs a region-level strategy, and
#' (ii) held-out AUC of an episignature classifier trained on CpG features vs
#' region features (50/50 stratified split).
#'
#' @param coverages,deltas Grids of mean coverage and effect sizes.
#' @param n_rep Replicates per setting.
#' @param n_per_group Samples per group (before the train/test split).
#' @param ... Passed to [simulate_methylation()].
#' @param seed RNG seed.
#' @param progress Print progress?
#' @return A long data.frame (coverage, delta, rep, method, metric, value).
#' @export
run_benchmark <- function(coverages = c(10, 20, 30), deltas = c(0.05, 0.10, 0.15),
                          n_rep = 5, n_per_group = 60, ..., seed = 1,
                          progress = TRUE) {
  set.seed(seed)
  out <- list()
  for (cv in coverages) for (d in deltas) for (r in seq_len(n_rep)) {
    if (progress) message(sprintf("coverage=%g delta=%g rep=%d", cv, d, r))
    sim <- simulate_methylation(n_cases = n_per_group, n_controls = n_per_group,
                                coverage = cv, delta = d, ...)
    truth <- sim$truth
    tc <- test_cpgs(sim)
    hit <- unique(tc$region[!is.na(tc$padj) & tc$padj < 0.05])
    m_cpg <- detection_metrics(seq_along(truth) %in% hit, truth)
    tr <- test_regions(sim)
    m_reg <- detection_metrics(!is.na(tr$padj) & tr$padj < 0.05, truth)
    y <- sim$pheno$group
    tr_idx <- unlist(lapply(levels(y), function(l) {
      i <- which(y == l); sample(i, floor(length(i) / 2))
    }))
    te_idx <- setdiff(seq_along(y), tr_idx)
    Xc <- t(sim$meth / pmax(sim$cov, 1))
    Xr <- t(aggregate_regions(sim$meth, sim$cov, sim$cpg$region)$beta)
    auc <- function(X) {
      f <- fit_episignature(X[tr_idx, ], y[tr_idx])
      as.numeric(pROC::auc(pROC::roc(y[te_idx],
                  predict_episignature(f, X[te_idx, ]), quiet = TRUE,
                  levels = levels(y), direction = "<")))
    }
    cell <- function(method, m, a) data.frame(
      coverage = cv, delta = d, rep = r, method = method,
      metric = c(names(m), "auc"), value = c(unname(m), a),
      stringsAsFactors = FALSE)
    out[[length(out) + 1]] <- rbind(cell("CpG-level", m_cpg, auc(Xc)),
                                    cell("Region-level", m_reg, auc(Xr)))
  }
  do.call(rbind, out)
}

#' Plot benchmark results
#' @param res Output of [run_benchmark()].
#' @param metric One of "tpr", "fdr", "auc".
#' @export
plot_benchmark <- function(res, metric = c("tpr", "fdr", "auc")) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  d <- res[res$metric == metric, ]
  lab <- c(tpr = "Sensitivity (region level)", fdr = "False discovery rate",
           auc = "Held-out AUC")[metric]
  ggplot2::ggplot(d, ggplot2::aes(factor(coverage), value, fill = method)) +
    ggplot2::geom_boxplot(outlier.size = 0.8) +
    ggplot2::facet_wrap(~ paste0("Delta = ", delta * 100, "%")) +
    ggplot2::labs(x = "Mean coverage (x)", y = lab, fill = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
