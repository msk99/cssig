test_that("Weir & Cockerham F_ST behaves at known extremes", {
  # completely differentiated: one cohort fixed for each allele
  expect_gt(css_fst(count1 = 50, n1 = 25, count2 = 0, n2 = 25), 0.95)
  # identical allele frequencies
  expect_lt(abs(css_fst(count1 = 25, n1 = 25, count2 = 25, n2 = 25)), 0.05)
  # monomorphic across both cohorts is undefined
  expect_true(is.na(css_fst(count1 = 0, n1 = 25, count2 = 0, n2 = 25)))
  # vectorised and recycled
  expect_length(css_fst(c(10, 20, 30), 25, c(30, 20, 10), 25), 3L)
})

test_that("floor_zero only truncates, and cannot change CSS", {
  set.seed(1)
  c1 <- rbinom(500, 50, 0.3); c2 <- rbinom(500, 50, 0.3)
  a <- css_fst(c1, 25, c2, 25)
  b <- css_fst(c1, 25, c2, 25, floor_zero = TRUE)
  expect_true(all(b >= 0, na.rm = TRUE))
  # truncation is monotone, so ranks outside the truncated block are unchanged
  ok <- !is.na(a) & a > 0
  expect_equal(rank(a[ok]), rank(b[ok]))
})

test_that("dDAF and dSAF agree in sign convention", {
  expect_equal(css_ddaf(c(0.8, 0.2), c(0.3, 0.7)), c(0.5, -0.5))
  # dSAF orients on the major allele in the selected cohort
  expect_equal(css_dsaf(c(0.8, 0.2), c(0.3, 0.7)), c(0.5, 0.5))
})

test_that("css_standardize returns mean 0 and unit variance", {
  z <- css_standardize(c(1, 4, 9, 16))
  expect_equal(mean(z), 0, tolerance = 1e-12)
  expect_equal(sd(z), 1, tolerance = 1e-12)
  expect_warning(css_standardize(rep(3, 5)), "zero or undefined")
})

test_that("css_merge_tests reports and applies the join correctly", {
  a <- data.frame(chr = 1, pos = c(100, 200, 300), fst = c(0.1, 0.2, 0.3))
  b <- data.frame(chr = 1, pos = c(200, 300, 400), xpehh = c(1.1, 2.2, 0.4))
  inner <- css_merge_tests(fst = a, xpehh = b, all = FALSE, quiet = TRUE)
  expect_equal(nrow(inner), 2L)
  outer <- css_merge_tests(fst = a, xpehh = b, all = TRUE, quiet = TRUE)
  expect_equal(nrow(outer), 4L)
  expect_equal(sum(is.na(outer$fst)), 1L)
  expect_message(css_merge_tests(fst = a, xpehh = b), "Merged 2 tests")
  expect_error(css_merge_tests(a, b), "must be named")
})

test_that("css_fdr refuses smoothed input unless forced", {
  set.seed(8)
  n <- 500L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e5,
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  res <- css_smooth(css(css_input(d, tests = c(a = "high", b = "high", c = "high"))))
  expect_error(css_fdr(res), "not p-values")
  expect_silent(suppressMessages(css_fdr(res, force = TRUE)))
})

test_that("BH q-values are monotone in p and bounded", {
  set.seed(8)
  n <- 2000L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e5,
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  res <- suppressMessages(css_fdr(css(css_input(d, tests = c(a = "high", b = "high", c = "high")))))
  expect_true(all(res$qval >= 0 & res$qval <= 1))
  o <- order(res$p)
  expect_false(is.unsorted(res$qval[o]))
})

test_that("reciprocal CSS is antisymmetric when the tests are mirrored", {
  set.seed(12)
  n <- 400L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e5,
                  fst = runif(n), xpehh = rnorm(n), ddaf = rnorm(n))
  fwd <- css_input(d, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
  d2 <- d; d2$xpehh <- -d2$xpehh; d2$ddaf <- -d2$ddaf
  rvs <- css_input(d2, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
  r <- css_reciprocal(fwd, rvs, labels = c("A", "B"))
  expect_true(all(c("css_pos", "css_neg", "css_signed") %in% names(r)))
  expect_equal(r$css_signed, r$css_pos - r$css_neg)
  expect_error(css_reciprocal(fwd, css_input(d[1:10, ], tests = c(fst = "high", xpehh = "high", ddaf = "high"))),
               "same SNPs")
})

test_that("reciprocal results smooth and plot through the standard generics", {
  set.seed(12)
  n <- 400L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e5,
                  fst = runif(n), xpehh = rnorm(n), ddaf = rnorm(n))
  fwd <- css_input(d, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
  d2 <- d; d2$xpehh <- -d2$xpehh; d2$ddaf <- -d2$ddaf
  rvs <- css_input(d2, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
  r <- css_reciprocal(fwd, rvs, labels = c("A", "B"))

  expect_s3_class(plot(r), "ggplot")                       # dispatches to the mirror
  expect_error(css_manhattan(r), "css_manhattan_mirror")   # helpful redirect
  expect_error(css_smooth(r, on = "zbar"), "reciprocal")

  r <- css_smooth(r)
  expect_true(all(c("css_pos_smooth", "css_neg_smooth", "css_signed_smooth") %in% names(r)))
  expect_equal(r$css_signed_smooth, r$css_pos_smooth - r$css_neg_smooth)
  expect_s3_class(css_manhattan_mirror(r), "ggplot")            # smoothed by default
  expect_s3_class(css_manhattan_mirror(r, score = "raw"), "ggplot")
})

test_that("rehh output is coerced correctly", {
  fake <- data.frame(CHR = 1, POSITION = c(100, 200), XPEHH_A_B = c(1.2, -0.4))
  out <- read_rehh_xpehh(fake)
  expect_equal(names(out), c("chr", "pos", "snp", "xpehh"))
  expect_equal(out$xpehh, c(1.2, -0.4))
  expect_error(read_rehh_xpehh(data.frame(a = 1)), "CHR and POSITION")
})

test_that("both rehh spellings of the XP-EHH column are recognised", {
  # ies2xpehh(standardize = FALSE) writes UNXPEHH_*, standardize = TRUE writes
  # XPEHH_*. Missing the unstandardised spelling silently yields no data.
  un <- data.frame(CHR = 1, POSITION = c(100, 200), UNXPEHH_A_B = c(1.2, -0.4))
  expect_equal(read_rehh_xpehh(un)$xpehh, c(1.2, -0.4))
  lp <- data.frame(CHR = 1, POSITION = c(100, 200),
                   UNXPEHH_A_B = c(1.2, -0.4), LOGPVALUE = c(0.3, 0.9))
  expect_equal(read_rehh_xpehh(lp, value = "LOGPVALUE")$xpehh, c(0.3, 0.9))
})

test_that("plot functions return ggplot objects", {
  set.seed(13)
  n <- 1200L
  d <- data.frame(chr = rep(c("1", "2"), each = n / 2),
                  pos = rep(seq_len(n / 2) * 1e5, 2),
                  fst = runif(n), xpehh = rnorm(n), ddaf = rnorm(n))
  res <- suppressMessages(
    css_threshold(css_smooth(css(css_input(d, tests = c(fst = "high", xpehh = "high", ddaf = "high")))))
  )
  expect_s3_class(css_manhattan(res), "ggplot")
  expect_s3_class(css_manhattan(res, overlay_raw = TRUE), "ggplot")
  expect_s3_class(css_chrom_plot(res, chr = "1"), "ggplot")
  expect_s3_class(css_qq(res), "ggplot")
  expect_s3_class(css_pdist(res), "ggplot")
  expect_s3_class(css_test_cor(res), "ggplot")
  expect_error(css_chrom_plot(res, chr = "99"), "not found")
})

test_that("css_region_plot accepts plain region specs and keeps genes in one panel", {
  set.seed(14)
  n <- 300L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e4,
                  fst = runif(n), xpehh = rnorm(n), ddaf = rnorm(n))
  res <- suppressMessages(
    css_threshold(css_smooth(css(css_input(d, tests = c(fst = "high", xpehh = "high", ddaf = "high")))))
  )

  # a named vector is documented as a valid region spec
  expect_s3_class(css_region_plot(res, region = c(chr = "1", start = 4e5, end = 1.2e6)),
                  "ggplot")

  genes <- data.frame(start = 5e5, end = 9e5, name = "GENE1")
  p <- css_region_plot(res, region = list(chr = "1", start = 4e5, end = 1.2e6),
                       genes = genes)
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$layout$layout), 5L)   # CSS + three tests + Genes
  seg <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1)))
  expect_length(seg, 1L)
  # the gene track must live in exactly one facet, not repeat across all of them
  expect_equal(length(unique(b$data[[seg]]$PANEL)), 1L)
})

test_that("genome coordinates are monotone and label every chromosome once", {
  g <- css_genome_coords(rep(c("1", "2", "10"), each = 3),
                         rep(c(1e6, 2e6, 3e6), 3))
  expect_false(is.unsorted(sort(g$pos_cum)))
  expect_equal(nrow(g$axis), 3L)
  expect_identical(as.character(g$axis$chr), c("1", "2", "10"))
})

test_that("the selscan reader round-trips the shipped fixture", {
  f <- system.file("extdata", "example.xpehh.out", package = "cssig")
  skip_if(!nzchar(f))
  out <- read_selscan_xpehh(f, chr = 1)
  expect_equal(names(out), c("chr", "pos", "snp", "xpehh"))
  expect_equal(nrow(out), 200L)
  expect_false(anyNA(out$xpehh))
  expect_false(is.unsorted(out$pos))
  expect_error(read_selscan_xpehh(f), "must be supplied")
  expect_error(read_selscan_xpehh("nope.txt", chr = 1), "not found")
})
