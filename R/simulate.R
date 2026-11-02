#' Simulate long-read (HiFi-like) CpG methylation counts
#'
#' Generates per-CpG methylated / total read counts for a case-control design.
#' A subset of regions carries a coherent methylation shift (`delta`) in cases,
#' mimicking an episignature. Read-level counts are beta-binomial: they combine
#' binomial sampling at a given coverage with extra-binomial variability, and
#' each individual has its own random shift per region (inter-individual
#' variability).
#'
#' @param n_cases,n_controls Number of cases / controls.
#' @param n_regions Number of genomic regions.
#' @param n_signal Number of regions carrying a true shift in cases.
#' @param cpgs_per_region Number of CpGs per region.
#' @param coverage Mean per-CpG read coverage (Poisson).
#' @param delta Absolute methylation difference (0-1 scale) in signal regions.
#' @param sd_ind SD (logit scale) of the individual-by-region random effect.
#' @param rho Beta-binomial overdispersion (0 = pure binomial).
#' @param seed Optional RNG seed.
#' @return A list of class `hifi_sim` with matrices `meth` and `cov`
#'   (CpG x sample), `cpg` (data.frame), `pheno` (data.frame) and `truth`
#'   (logical per region).
#' @export
simulate_methylation <- function(n_cases = 40, n_controls = 40, n_regions = 200,
                                 n_signal = 20, cpgs_per_region = 10,
                                 coverage = 30, delta = 0.10, sd_ind = 0.15,
                                 rho = 0.02, seed = NULL) {
  stopifnot(n_signal <= n_regions, delta >= 0, delta < 0.5, coverage > 0)
  if (!is.null(seed)) set.seed(seed)
  n <- n_cases + n_controls
  group <- factor(rep(c("control", "case"), c(n_controls, n_cases)),
                  levels = c("control", "case"))
  region_id <- rep(seq_len(n_regions), each = cpgs_per_region)
  n_cpg <- length(region_id)

  truth <- seq_len(n_regions) %in% sample.int(n_regions, n_signal)
  sign_r <- ifelse(stats::runif(n_regions) < 0.5, -1, 1)
  base_r <- stats::rbeta(n_regions, 0.8, 0.8)
  base_r[truth] <- stats::runif(sum(truth), 0.2, 0.8)  # room for the shift
  base_r <- pmin(pmax(base_r, 0.03), 0.97)
  base_c <- stats::plogis(stats::qlogis(base_r[region_id]) +
                            stats::rnorm(n_cpg, 0, 0.3))

  u <- matrix(stats::rnorm(n_regions * n, 0, sd_ind), n_regions, n)
  p <- stats::plogis(stats::qlogis(base_c) + u[region_id, , drop = FALSE])
  shift <- outer((truth * sign_r)[region_id], as.numeric(group == "case")) * delta
  p <- pmin(pmax(p + shift, 0.005), 0.995)

  cov <- matrix(pmax(stats::rpois(n_cpg * n, coverage), 1L), n_cpg, n)
  if (rho > 0) {
    k <- 1 / rho - 1
    p <- matrix(stats::rbeta(n_cpg * n, p * k, (1 - p) * k), n_cpg, n)
  }
  meth <- matrix(stats::rbinom(n_cpg * n, cov, p), n_cpg, n)

  ids <- paste0("S", seq_len(n))
  dimnames(meth) <- dimnames(cov) <- list(paste0("cpg", seq_len(n_cpg)), ids)
  structure(list(
    meth = meth, cov = cov,
    cpg = data.frame(cpg = rownames(meth), region = region_id,
                     stringsAsFactors = FALSE),
    pheno = data.frame(sample = ids, group = group, stringsAsFactors = FALSE),
    truth = stats::setNames(truth, seq_len(n_regions))
  ), class = "hifi_sim")
}
