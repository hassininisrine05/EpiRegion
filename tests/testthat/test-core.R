test_that("simulation has the expected structure", {
  s <- simulate_methylation(10, 10, n_regions = 30, n_signal = 5,
                            cpgs_per_region = 6, seed = 1)
  expect_equal(dim(s$meth), c(180, 20))
  expect_true(all(s$meth <= s$cov))
  expect_equal(sum(s$truth), 5)
})

test_that("aggregation sums counts", {
  s <- simulate_methylation(5, 5, n_regions = 4, n_signal = 1,
                            cpgs_per_region = 3, seed = 2)
  a <- aggregate_regions(s$meth, s$cov, s$cpg$region)
  expect_equal(unname(a$meth[1, ]), unname(colSums(s$meth[1:3, ])))
  expect_true(all(a$beta >= 0 & a$beta <= 1))
})

test_that("region test finds strong signal and controls FDR", {
  s <- simulate_methylation(40, 40, n_regions = 100, n_signal = 20,
                            delta = 0.15, seed = 3)
  m <- detection_metrics(test_regions(s)$padj < 0.05, s$truth)
  expect_gt(m[["tpr"]], 0.6)
  expect_lt(m[["fdr"]], 0.2)
})

test_that("null data yields few calls", {
  s <- simulate_methylation(30, 30, n_regions = 100, n_signal = 0, seed = 4)
  expect_lte(sum(test_regions(s)$padj < 0.05, na.rm = TRUE), 3)
})

test_that("classifier separates groups", {
  s <- simulate_methylation(40, 40, n_regions = 60, n_signal = 15,
                            delta = 0.15, seed = 5)
  X <- t(aggregate_regions(s$meth, s$cov, s$cpg$region)$beta); y <- s$pheno$group
  idx <- c(1:20, 41:60)
  f <- fit_episignature(X[idx, ], y[idx])
  p <- predict_episignature(f, X[-idx, ])
  au <- as.numeric(pROC::auc(pROC::roc(y[-idx], p, quiet = TRUE,
                                       levels = levels(y), direction = "<")))
  expect_gt(au, 0.8)
})

test_that("bedMethyl reader filters by mod code", {
  f <- tempfile()
  rows <- rbind(
    c("chr1", 10, 11, "m", 5, "+", 10, 11, "255,0,0", 30, 50, 15, 15, 0, 0, 0, 0, 0),
    c("chr1", 10, 11, "a", 5, "+", 10, 11, "255,0,0", 20, 10, 2, 18, 0, 0, 0, 0, 0))
  write.table(rows, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  d <- read_bedmethyl(f, "m")
  expect_equal(nrow(d), 1); expect_equal(d$meth, 15)
  expect_equal(read_bedmethyl(f, "a")$cov, 20)
})

test_that("tile_regions bins positions", {
  expect_equal(tile_regions(c("chr1", "chr1"), c(5, 1500), 1000),
               c("chr1:0", "chr1:1000"))
})
