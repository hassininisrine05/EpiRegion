#' Aggregate CpG counts into regions
#'
#' Sums methylated and total read counts over the CpGs of each region. Pooling
#' counts (rather than averaging per-CpG proportions) weights CpGs by coverage
#' and gains precision, which is the key idea for modest-depth long-read data.
#'
#' @param meth,cov CpG x sample matrices of methylated / total counts.
#' @param region Vector (length nrow(meth)) giving the region of each CpG.
#' @return List with region x sample matrices `meth`, `cov`, `beta`.
#' @export
aggregate_regions <- function(meth, cov, region) {
  stopifnot(identical(dim(meth), dim(cov)), length(region) == nrow(meth))
  m <- rowsum(meth, region, reorder = FALSE)
  c <- rowsum(cov, region, reorder = FALSE)
  list(meth = m, cov = c, beta = m / pmax(c, 1))
}

#' Tile CpG positions into fixed-size windows
#' @param chrom,pos Chromosome and position vectors.
#' @param width Window width in bp.
#' @return Character vector of region identifiers.
#' @export
tile_regions <- function(chrom, pos, width = 1000) {
  paste0(chrom, ":", floor(pos / width) * width)
}

#' Read a modkit bedMethyl file
#'
#' Reads the standard `modkit pileup` bedMethyl output
#' (column 4 = modification code, 10 = valid coverage, 12 = modified count).
#' Use `mod_code = "m"` for 5mC and `"a"` for 6mA. Column layout should be
#' checked against the modkit version in use.
#'
#' @param file Path to a (possibly gzipped) bedMethyl file.
#' @param mod_code Modification code to keep.
#' @param min_cov Minimum valid coverage.
#' @return data.frame with chrom, start, end, cov, meth.
#' @export
read_bedmethyl <- function(file, mod_code = "m", min_cov = 1) {
  d <- utils::read.table(file, sep = "\t", header = FALSE, quote = "",
                         comment.char = "", stringsAsFactors = FALSE,
                         fill = TRUE)
  stopifnot(ncol(d) >= 12)
  d <- d[d[[4]] == mod_code & d[[10]] >= min_cov, , drop = FALSE]
  data.frame(chrom = d[[1]], start = d[[2]], end = d[[3]],
             cov = d[[10]], meth = d[[12]], stringsAsFactors = FALSE)
}
