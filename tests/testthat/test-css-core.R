micro <- function() {
  data.frame(
    chr = 1, pos = (1:5) * 1e5,
    t1 = c(0.1, 0.5, 0.9, 0.3, 0.7),
    t2 = c(5, 1, 9, 3, 7),
    t3 = c(2, 4, 10, 1, 6)
  )
}
tests3 <- c(t1 = "high", t2 = "high", t3 = "high")

test_that("CSS matches the published formula computed by hand", {
  d <- micro()
  r <- css(css_input(d, tests = tests3))

  n <- 5L; m <- 3L
  z <- sapply(c("t1", "t2", "t3"), function(k) qnorm(rank(d[[k]]) / (n + 1)))
  zbar <- rowMeans(z)
  p <- 1 - pnorm(sqrt(m) * zbar)

  expect_equal(r$zbar, zbar, tolerance = 1e-12)
  expect_equal(r$p, p, tolerance = 1e-12)
  expect_equal(r$css, -log10(p), tolerance = 1e-12)
})

test_that("CSS is invariant to strictly monotone transforms of any test", {
  d <- micro()
  base <- css(css_input(d, tests = tests3))$css

  d2 <- d
  d2$t1 <- exp(3 * d2$t1)
  d2$t2 <- log(d2$t2)
  d2$t3 <- d2$t3^3
  expect_equal(css(css_input(d2, tests = tests3))$css, base, tolerance = 1e-12)
})

test_that("standardising a constituent test does not change CSS", {
  d <- micro()
  base <- css(css_input(d, tests = tests3))$css
  d$t2 <- css_standardize(d$t2)
  expect_equal(css(css_input(d, tests = tests3))$css, base, tolerance = 1e-12)
})

test_that("direction = 'low' on a negated test equals 'high' on the original", {
  d <- micro()
  base <- css(css_input(d, tests = tests3))$css
  d$t1 <- -d$t1
  flipped <- css(css_input(d, tests = c(t1 = "low", t2 = "high", t3 = "high")))$css
  expect_equal(flipped, base, tolerance = 1e-12)
})

test_that("the SNP ranked first in every test attains the maximum CSS", {
  d <- micro()
  r <- css(css_input(d, tests = tests3))
  expect_equal(which.max(r$css), 3L)      # SNP 3 is top of all three tests
  expect_true(all(is.finite(r$css)))
})

test_that("CSS p-values are uniform under the null", {
  skip_on_cran()
  set.seed(11)
  n <- 50000L
  d <- data.frame(chr = 1, pos = seq_len(n),
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  r <- css(css_input(d, tests = c(a = "high", b = "high", c = "high")))
  expect_gt(suppressWarnings(ks.test(r$p, "punif")$p.value), 0.01)
  # E[-log10(U)] = 1/ln(10)
  expect_equal(mean(r$css), 1 / log(10), tolerance = 0.02)
})

test_that("extreme p-values are floored rather than becoming infinite", {
  set.seed(3)
  n <- 200000L
  # perfectly concordant tests give the most extreme attainable CSS
  v <- rnorm(n)
  d <- data.frame(chr = 1, pos = seq_len(n), a = v, b = v, c = v)
  r <- suppressMessages(css(css_input(d, tests = c(a = "high", b = "high", c = "high"))))
  expect_true(all(is.finite(r$css)))
  expect_true(all(r$p > 0))
})

test_that("m = 1 warns that CSS is degenerate", {
  d <- micro()
  expect_warning(css_input(d, tests = c(t1 = "high")), "one constituent test")
})

test_that("a large single-value block is reported, ordinary discreteness is not", {
  d <- data.frame(chr = 1, pos = 1:100,
                  a = c(rep(0, 60), rnorm(40)), b = rnorm(100), c = rnorm(100))
  # Discrete-but-spread values must NOT warn: statistics built from allele
  # counts are inherently discrete, and a naive tied-fraction check would fire
  # on every real dataset. Modelled on a frequency difference between two
  # cohorts of 400 haplotypes.
  set.seed(21)
  freq_diff <- function(n) (sample.int(401, n, TRUE) - sample.int(401, n, TRUE)) / 400
  spread <- data.frame(chr = 1, pos = 1:2000,
                       a = freq_diff(2000), b = freq_diff(2000), c = freq_diff(2000))
  expect_no_warning(css(css_input(spread, tests = c(a = "high", b = "high", c = "high"))))
  expect_warning(css(css_input(d, tests = c(a = "high", b = "high", c = "high"))),
                 "share a single value")
})

test_that("frank reproduces base rank for every supported ties method", {
  set.seed(7)
  z <- round(rnorm(500), 1)
  for (tm in c("average", "first", "dense")) {
    expect_equal(as.numeric(data.table::frank(z, ties.method = tm)),
                 as.numeric(rank(z, ties.method = if (tm == "dense") "min" else tm)),
                 tolerance = if (tm == "dense") Inf else 1e-12,
                 info = tm)
  }
  expect_equal(as.numeric(data.table::frank(z, ties.method = "average")),
               as.numeric(rank(z)))
})

test_that("na_action = 'pairwise' uses each SNP's own m", {
  d <- micro()
  d$t2[2] <- NA
  r <- css(css_input(d, tests = tests3))
  expect_true("n_tests" %in% names(r))
  expect_equal(r$n_tests, c(3L, 2L, 3L, 3L, 3L))
  expect_false(is.na(r$css[2]))
})

test_that("na_action = 'omit' drops incomplete SNPs", {
  d <- micro()
  d$t2[2] <- NA
  r <- suppressMessages(css(css_input(d, tests = tests3), na_action = "omit"))
  expect_equal(nrow(r), 4L)
})

test_that("weights are rejected when malformed and change the result when used", {
  d <- micro()
  x <- css_input(d, tests = tests3)
  expect_error(css(x, weights = c(1, 2)), "one entry per test")
  expect_error(css(x, weights = c(1, -1, 1)), "strictly positive")
  a <- css(x, .copy = TRUE)$css
  b <- css(x, weights = c(3, 1, 1), .copy = TRUE)$css
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("equal weights reduce exactly to the published unweighted statistic", {
  # The weighted branch computes sum(w*Z)/sqrt(sum(w^2)); at w = 1 this is
  # sum(Z)/sqrt(m) = sqrt(m)*mean(Z), the published statistic. If these ever
  # diverge, `weights` has stopped being a strict generalisation.
  d <- micro()
  x <- css_input(d, tests = tests3)
  plain <- css(x, .copy = TRUE)
  eq    <- css(x, weights = c(1, 1, 1), .copy = TRUE)
  expect_equal(plain$css, eq$css, tolerance = 1e-12)
  expect_equal(plain$p,   eq$p,   tolerance = 1e-12)
})

test_that("weights are scale-invariant", {
  # Only the ratios matter; multiplying every weight by a constant must not
  # change the statistic.
  d <- micro()
  x <- css_input(d, tests = tests3)
  a <- css(x, weights = c(1, 2, 3), .copy = TRUE)$css
  b <- css(x, weights = c(10, 20, 30), .copy = TRUE)$css
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("a down-weighted test moves CSS towards the others", {
  d <- micro()
  d$t3 <- rev(d$t3)   # make t3 disagree with t1/t2
  x <- css_input(d, tests = tests3)
  full <- css(x, .copy = TRUE)$css
  muted <- css(x, weights = c(1, 1, 0.01), .copy = TRUE)$css
  two_only <- css(css_input(d, tests = c(t1 = "high", t2 = "high")), .copy = TRUE)$css
  # nearly zero weight on t3 should approach the two-test answer more closely
  expect_lt(max(abs(muted - two_only)), max(abs(full - two_only)))
})
