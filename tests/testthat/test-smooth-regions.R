test_that("the window smoother matches a naive reference implementation", {
  set.seed(42)
  n <- 1500L
  d <- data.table::data.table(
    chr = rep(c("1", "2"), each = n / 2),
    pos = c(cumsum(sample(1:80000, n / 2, TRUE)),
            cumsum(sample(1:80000, n / 2, TRUE))),
    v = rnorm(n))
  d[sample(n, 40), v := NA]
  data.table::setkey(d, chr, pos)
  hw <- 5e5

  naive_m <- naive_n <- numeric(n)
  for (i in seq_len(n)) {
    sel <- d$chr == d$chr[i] & d$pos >= d$pos[i] - hw &
           d$pos <= d$pos[i] + hw & !is.na(d$v)
    naive_m[i] <- mean(d$v[sel])
    naive_n[i] <- sum(sel)
  }
  res <- d[, cssig:::.window_mean(pos, v, hw), by = chr]
  expect_equal(res$mean, naive_m)
  expect_identical(as.integer(res$n), as.integer(naive_n))
})

test_that("windows never span a chromosome boundary", {
  b <- data.table::data.table(chr = c("1", "1", "2", "2"),
                              pos = c(1e6, 1.1e6, 1.05e6, 1.2e6),
                              v = c(1, 1, 100, 100))
  data.table::setkey(b, chr, pos)
  r <- b[, cssig:::.window_mean(pos, v, 5e5), by = chr]
  expect_equal(r$mean, c(1, 1, 100, 100))
})

test_that("duplicate positions and all-NA windows are handled", {
  b <- data.table::data.table(chr = "1", pos = c(1, 1, 2, 3) * 1e5,
                              v = c(1, 3, NA, 5))
  data.table::setkey(b, chr, pos)
  r <- b[, cssig:::.window_mean(pos, v, 5e5), by = chr]
  expect_equal(r$n, rep(3L, 4))
  expect_equal(r$mean, rep(3, 4))

  allna <- data.table::data.table(chr = "1", pos = c(1, 2) * 1e5,
                                  v = c(NA_real_, NA_real_))
  r2 <- allna[, cssig:::.window_mean(pos, v, 5e5), by = chr]
  expect_true(all(is.na(r2$mean)))
  expect_equal(r2$n, c(0L, 0L))
})

test_that("min_snps masks sparse windows", {
  set.seed(2)
  d <- data.frame(chr = "1",
                  pos = c(seq(1e6, 2e6, length.out = 40), 50e6, 90e6),
                  a = rnorm(42), b = rnorm(42), c = rnorm(42))
  x <- css(css_input(d, tests = c(a = "high", b = "high", c = "high")))
  r <- suppressMessages(css_smooth(x, min_snps = 5))
  # the two isolated SNPs sit alone in their windows
  expect_true(all(is.na(tail(r$css_smooth, 2))))
  expect_false(any(is.na(head(r$css_smooth, 40))))
})

test_that("on = 'zbar' refuses weighted results", {
  set.seed(3)
  d <- data.frame(chr = "1", pos = seq_len(50) * 1e5,
                  a = rnorm(50), b = rnorm(50), c = rnorm(50))
  r <- css(css_input(d, tests = c(a = "high", b = "high", c = "high")),
           weights = c(1, 1, 2))
  expect_error(css_smooth(r, on = "zbar"), "weighted")
})

test_that("region calling recovers an implanted signal and respects the merge gap", {
  set.seed(9)
  n <- 4000L
  d <- data.frame(chr = "1", pos = seq_len(n) * 5e4,
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  hit <- 2000:2020                        # a 1 Mb block of concordant signal
  d$a[hit] <- d$a[hit] + 6
  d$b[hit] <- d$b[hit] + 6
  d$c[hit] <- d$c[hit] + 6

  res <- suppressMessages(
    css_threshold(css_smooth(css(css_input(d, tests = c(a = "high", b = "high", c = "high")))))
  )
  reg <- css_regions(res, method = "cluster")
  expect_gte(nrow(reg), 1L)
  true_pos <- d$pos[2010]
  expect_true(any(reg$start <= true_pos & reg$end >= true_pos))

  reg2 <- css_regions(res, method = "flank")
  expect_true(any(reg2$start <= true_pos & reg2$end >= true_pos))
})

# Two signal blocks 8.5 Mb apart, used by the merge and padding tests below.
two_block_result <- function(seed = 4) {
  set.seed(seed)
  n <- 4000L
  d <- data.frame(chr = "1", pos = seq_len(n) * 5e4,
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  for (blk in list(1000:1030, 1170:1200)) {   # pos 50.0-51.5 Mb and 58.5-60.0 Mb
    d$a[blk] <- d$a[blk] + 6; d$b[blk] <- d$b[blk] + 6; d$c[blk] <- d$c[blk] + 6
  }
  suppressMessages(
    css_threshold(css_smooth(css(css_input(d, tests = c(a = "high", b = "high", c = "high")))),
                  top = 0.02, top2 = 0.05)
  )
}

test_that("merge_gap controls whether nearby clusters become one region", {
  res <- two_block_result()
  merged <- css_regions(res, merge_gap = 20e6)   # wider than the 8.5 Mb gap
  split  <- css_regions(res, merge_gap = 1e6)    # narrower than it
  expect_equal(nrow(merged), 1L)
  expect_equal(nrow(split), 2L)
  # merging spans both blocks
  expect_gt(merged$end - merged$start, 8e6)
})

test_that("padded boundaries are reported separately and never inflate start/end", {
  res <- two_block_result()
  reg <- css_regions(res, merge_gap = 1e6, flank = 5e5)
  expect_gt(nrow(reg), 0L)
  expect_true(all(reg$start_padded <= reg$start))
  expect_true(all(reg$end_padded >= reg$end))
  expect_equal(reg$end_padded - reg$end, rep(5e5, nrow(reg)))
  expect_true(all(reg$start_padded >= 0))
})

test_that("thresholds keep approximately the requested tail fraction", {
  set.seed(6)
  n <- 20000L
  d <- data.frame(chr = "1", pos = seq_len(n) * 1e4,
                  a = rnorm(n), b = rnorm(n), c = rnorm(n))
  res <- suppressMessages(
    css_threshold(css_smooth(css(css_input(d, tests = c(a = "high", b = "high", c = "high")))),
                  top = 0.01))
  frac <- mean(res$significant, na.rm = TRUE)
  expect_gt(frac, 0.005)
  expect_lt(frac, 0.02)
})
