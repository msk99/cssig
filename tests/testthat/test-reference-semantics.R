# The single most important regression suite in the package: data.table's `:=`
# modifies in place, so a careless implementation would make columns appear in
# the user's own object as a side effect of calling css(). See PLAN.md 3.7.

make_df <- function(n = 200) {
  set.seed(5)
  data.frame(chr = rep(c("1", "2"), each = n / 2),
             pos = c(seq_len(n / 2), seq_len(n / 2)) * 1e4,
             fst = runif(n), xpehh = rnorm(n), ddaf = rnorm(n))
}
tests3 <- c(fst = "high", xpehh = "high", ddaf = "high")

test_that("css_input() does not modify a data.frame the user passes in", {
  d <- make_df()
  before <- d
  css_input(d, tests = tests3)
  expect_identical(d, before)
})

test_that("css_input() does not modify a data.table the user passes in", {
  dt <- data.table::as.data.table(make_df())
  before <- data.table::copy(dt)
  addr <- data.table::address(dt)

  x <- css_input(dt, tests = tests3)
  css(x)
  css_smooth(x)
  suppressMessages(css_threshold(x))

  expect_identical(dt, before)
  expect_identical(names(dt), names(before))
  expect_identical(data.table::address(dt), addr)
})

test_that("pipeline stages add columns by reference to their own object", {
  x <- css_input(make_df(), tests = tests3)
  res <- css(x)
  # `res` and `x` are the same object: the pipeline avoids copying
  expect_true("css" %in% names(x))
  expect_identical(data.table::address(res), data.table::address(x))
})

test_that(".copy = TRUE leaves the input stage untouched", {
  x <- css_input(make_df(), tests = tests3)
  res <- css(x, .copy = TRUE)
  expect_false("css" %in% names(x))
  expect_true("css" %in% names(res))
})

test_that("all input dialects give identical results", {
  d <- make_df()
  a <- css(css_input(d, tests = tests3))
  b <- css(css_input(data.table::as.data.table(d), tests = tests3))
  expect_equal(a$css, b$css)

  skip_if_not_installed("tibble")
  c3 <- css(css_input(tibble::as_tibble(d), tests = tests3))
  expect_equal(a$css, c3$css)
})

test_that("results do not depend on input row order", {
  d <- make_df()
  a <- css(css_input(d, tests = tests3))
  shuffled <- d[sample(nrow(d)), ]
  b <- css(css_input(shuffled, tests = tests3))
  expect_equal(a$css, b$css)
  expect_equal(as.character(a$chr), as.character(b$chr))
  expect_equal(a$pos, b$pos)
})

test_that("the (chr, pos) key survives every stage", {
  x <- css_input(make_df(), tests = tests3)
  for (stage in list(css, css_smooth, function(z) suppressMessages(css_threshold(z)))) {
    x <- stage(x)
    expect_true(data.table::haskey(x))
    expect_identical(data.table::key(x), c("chr", "pos"))
  }
})

test_that("results are identical regardless of data.table thread count", {
  old <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old))
  d <- make_df(2000)
  data.table::setDTthreads(1)
  a <- css_smooth(css(css_input(d, tests = tests3)))$css_smooth
  data.table::setDTthreads(4)
  b <- css_smooth(css(css_input(d, tests = tests3)))$css_smooth
  expect_identical(a, b)
})

test_that("column names that collide with CSS output are rejected", {
  d <- make_df()
  names(d)[names(d) == "fst"] <- "css"
  expect_error(css_input(d, tests = c(css = "high", xpehh = "high", ddaf = "high")),
               "reserved")
})

test_that("duplicate SNP positions are rejected", {
  d <- make_df(10)
  d$pos[2] <- d$pos[1]
  expect_error(css_input(d, tests = tests3), "duplicated")
})

test_that("chromosomes sort numerically, not alphabetically", {
  d <- data.frame(chr = c("10", "2", "1", "X"), pos = 1:4 * 1e5,
                  a = 1:4, b = 4:1, c = c(2, 3, 1, 4))
  x <- css_input(d, tests = c(a = "high", b = "high", c = "high"))
  expect_identical(levels(x$chr), c("1", "2", "10", "X"))
})

test_that("an existing `snp` column is used rather than discarded", {
  d <- make_df(10)
  d$snp <- paste0("rs", seq_len(nrow(d)))
  x <- css_input(d, tests = tests3)
  expect_identical(x$snp, d$snp[order(d$chr, d$pos)])

  # and chr:pos ids are still generated when there is no snp column
  d2 <- make_df(10)
  x2 <- css_input(d2, tests = tests3)
  expect_true(all(grepl("^[12]:", x2$snp)))
})

test_that("chromosome factor levels give genomic, not alphabetical, order", {
  # With chr as character a keyed table sorts 1, 10, 2 -- which would scramble
  # the Manhattan axis and the order of called regions.
  d <- data.frame(chr = c("10", "2", "1", "X"), pos = 1:4 * 1e5,
                  a = 1:4, b = 4:1, c = c(2, 3, 1, 4))
  x <- css_input(d, tests = c(a = "high", b = "high", c = "high"))
  expect_identical(as.character(x$chr), c("1", "2", "10", "X"))
  # filtering by chromosome still behaves naturally despite the factor
  expect_equal(nrow(x[as.character(chr) == "2"]), 1L)
})
