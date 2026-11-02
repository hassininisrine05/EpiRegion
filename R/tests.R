#' Quasi-binomial differential methylation test
#'
#' For each row (CpG or region) fits a quasi-binomial GLM of methylated vs
#' unmethylated counts on group. The quasi-binomial dispersion absorbs
#' inter-individual variability beyond binomial sampling noise, which is what
#' keeps type I error under control when pooling reads across CpGs.
#'
#' @param meth,cov Feature x sample count matrices.
#' @param group Two-level factor (control, case).
#' @return data.frame: id, delta_beta (case - control, pooled), stat, p, padj.
#' @export
test_differential <- function(meth, cov, group) {
  group <- factor(group)
  stopifnot(nlevels(group) == 2, length(group) == ncol(meth))
  X <- cbind(1, as.numeric(group) - 1)
  g1 <- group == levels(group)[2]
  res <- t(vapply(seq_len(nrow(meth)), function(i) {
    m <- meth[i, ]; n <- cov[i, ]
    out <- c(NA_real_, NA_real_)
    ok <- n > 0
    if (sum(ok) < 4 || length(unique(g1[ok])) < 2) return(out)
    tryCatch(suppressWarnings({
      fit <- stats::glm.fit(X[ok, , drop = FALSE], cbind(m[ok], n[ok] - m[ok]),
                            family = stats::quasibinomial())
      df <- fit$df.residual
      disp <- max(sum(fit$weights * fit$residuals^2) / df, 1e-8)
      R <- fit$qr$qr[1:2, 1:2, drop = FALSE]
      v <- chol2inv(R)[2, 2] * disp
      z <- unname(fit$coefficients[2]) / sqrt(v)
      c(z, 2 * stats::pt(-abs(z), df))
    }), error = function(e) out)
  }, numeric(2)))
  b1 <- rowSums(meth[, g1, drop = FALSE]) / pmax(rowSums(cov[, g1, drop = FALSE]), 1)
  b0 <- rowSums(meth[, !g1, drop = FALSE]) / pmax(rowSums(cov[, !g1, drop = FALSE]), 1)
  data.frame(id = rownames(meth), delta_beta = b1 - b0, stat = res[, 1],
             p = res[, 2], padj = stats::p.adjust(res[, 2], "BH"),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Single-CpG differential test on a `hifi_sim` object
#' @param sim Object from [simulate_methylation()].
#' @export
test_cpgs <- function(sim) {
  r <- test_differential(sim$meth, sim$cov, sim$pheno$group)
  r$region <- sim$cpg$region
  r
}

#' Region-level differential test on a `hifi_sim` object
#' @param sim Object from [simulate_methylation()].
#' @export
test_regions <- function(sim) {
  a <- aggregate_regions(sim$meth, sim$cov, sim$cpg$region)
  test_differential(a$meth, a$cov, sim$pheno$group)
}

#' Detection metrics at region level
#'
#' @param called Logical vector, one entry per region (TRUE = called).
#' @param truth Logical vector, true signal regions.
#' @return Named vector: tpr (sensitivity), fdr, n_called.
#' @export
detection_metrics <- function(called, truth) {
  tp <- sum(called & truth); fp <- sum(called & !truth)
  c(tpr = if (sum(truth)) tp / sum(truth) else NA_real_,
    fdr = if (tp + fp) fp / (tp + fp) else 0,
    n_called = tp + fp)
}
