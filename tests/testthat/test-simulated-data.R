# Tests against the shipped simulated data, where the answers are known.
#
# Deliberately absent: an assertion that the founder-effect trap escapes
# detection. A founder-effect haplotype block has genuinely elevated F_ST *and*
# genuinely extended haplotype homozygosity, which are the same two signatures
# a sweep leaves. No combination of these constituent tests separates them.
# Asserting otherwise would encode a claim the method cannot support.

tests3 <- c(fst = "high", xpehh = "high", ddaf = "high")

analyse <- function(dat, top = 0.001) {
  suppressWarnings(suppressMessages(
    css_threshold(css_smooth(css(css_input(dat, tests = tests3))), top = top)
  ))
}

test_that("the shipped datasets are well formed", {
  data(css_sim); data(css_sim_small); data(css_sim_truth)
  expect_s3_class(css_sim, "data.table")
  expect_true(all(c("chr", "pos", "fst", "xpehh", "ddaf", "dsaf") %in% names(css_sim)))
  expect_equal(data.table::uniqueN(css_sim$chr), 29L)
  expect_false(anyNA(css_sim$fst))
  expect_false(anyNA(css_sim$xpehh))
  expect_true(all(css_sim$maf >= 0.01 & css_sim$maf <= 0.5))
  expect_equal(nrow(css_sim_truth), 7L)
  expect_true(all(as.character(css_sim_small$chr) %in% c("1", "2")))
})

test_that("positions are unique and sorted within chromosome", {
  data(css_sim)
  d <- data.table::as.data.table(css_sim)
  expect_false(any(duplicated(d, by = c("chr", "pos"))))
  expect_true(all(d[, all(diff(pos) > 0), by = chr]$V1))
})

test_that("causal variants are excluded from the panel, as on a real chip", {
  data(css_sim); data(css_sim_truth)
  d <- data.table::as.data.table(css_sim)
  tr <- data.table::as.data.table(css_sim_truth)[selected == TRUE]
  for (i in seq_len(nrow(tr))) {
    expect_equal(nrow(d[as.character(chr) == tr$chr[i] & pos == tr$pos[i]]), 0L)
  }
})

# Percentile of the 1 Mb bin containing each locus, ranked by its maximum
# smoothed CSS against every other 1 Mb bin in the genome. Taking the maximum
# over a window and comparing it to the *per-SNP* distribution would not be
# calibrated: the maximum of ~18 SNPs sits near the 95th percentile of that
# distribution by chance alone. Comparing bin maxima against bin maxima is.
bin_percentile <- function(res, chr_i, pos_i, col = "css_smooth") {
  d <- data.table::as.data.table(res)[!is.na(get(col))]
  b <- d[, .(m = max(get(col))), by = .(chr, bin = floor(pos / 1e6))]
  b[, pct := data.table::frank(m) / .N]
  v <- b[as.character(chr) == chr_i & bin == floor(pos_i / 1e6), pct]
  if (!length(v)) NA_real_ else v
}

test_that("the strongest sweep is recovered and well localised", {
  skip_on_cran()
  data(css_sim); data(css_sim_truth)
  res <- analyse(css_sim)
  reg <- css_regions(res, method = "cluster")
  expect_gt(nrow(reg), 0L)

  strong <- data.table::as.data.table(css_sim_truth)[sweep_id == "1"]
  # Padded boundaries, because the causal variant is excluded from the panel
  # exactly as it would be from a real SNP chip, so the nearest significant SNP
  # generally falls just outside the region. Both source papers add 0.5 Mb on
  # each side for precisely this reason when mining regions for genes.
  hit <- reg[as.character(chr) == strong$chr &
             start_padded <= strong$pos & end_padded >= strong$pos]
  expect_equal(nrow(hit), 1L)
  expect_lt(min(abs(hit$peak_pos - strong$pos)), 1e6)

  # and it is the top-ranked megabase in the genome
  expect_gt(bin_percentile(res, strong$chr, strong$pos), 0.999)
})

test_that("sweep strength orders the CSS signal", {
  skip_on_cran()
  data(css_sim); data(css_sim_truth)
  res <- analyse(css_sim)
  tr <- data.table::as.data.table(css_sim_truth)[cohort == "selected"]

  pctl <- vapply(seq_len(nrow(tr)),
                 function(i) bin_percentile(res, tr$chr[i], tr$pos[i]), numeric(1))

  # strongest sweep outranks weakest
  s_order <- order(tr$s)
  expect_gt(pctl[s_order[length(s_order)]], pctl[s_order[1]])
  # the weakest (s = 0.10) is genuinely not detectable, by design
  expect_lt(min(pctl, na.rm = TRUE), 0.6)
  # correlation between selection strength and signal
  expect_gt(stats::cor(tr$s, pctl, method = "spearman"), 0.5)
})

test_that("XP-EHH points at the reference cohort for the reference-cohort sweep", {
  skip_on_cran()
  data(css_sim); data(css_sim_truth)
  d <- data.table::as.data.table(css_sim)
  tr <- data.table::as.data.table(css_sim_truth)[sweep_id == "6"]
  w <- d[as.character(chr) == tr$chr & abs(pos - tr$pos) < 5e5]

  # XP-EHH is computed selected-against-reference, so a sweep in the reference
  # cohort must drive it negative. This is the constituent that carries
  # reliable direction.
  expect_lt(mean(w$xpehh), -0.5)
  expect_gt(mean(w$xpehh < 0), 0.7)

  # dDAF does NOT carry reliable direction away from the causal site: at a
  # hitchhiking SNP the allele riding the swept haplotype is ancestral or
  # derived roughly at random, so the sign of dDAF is scrambled. Assert only
  # that its magnitude is elevated, not its sign.
  expect_gt(mean(abs(w$ddaf)), mean(abs(d$ddaf)))
})

test_that("reciprocal CSS produces both directions over the same SNPs", {
  skip_on_cran()
  data(css_sim)
  rev_dat <- data.table::copy(data.table::as.data.table(css_sim))
  rev_dat[, xpehh := -xpehh][, ddaf := -ddaf]
  recip <- suppressWarnings(suppressMessages(css_reciprocal(
    css_input(css_sim, tests = tests3), css_input(rev_dat, tests = tests3))))
  expect_equal(nrow(recip), nrow(css_sim))
  expect_equal(recip$css_signed, recip$css_pos - recip$css_neg)
  # neither direction is systematically inflated
  expect_equal(mean(recip$css_pos), mean(recip$css_neg), tolerance = 0.05)
})

test_that("F_ST is not systematically inflated across sweep chromosomes", {
  # The uniform-forward-treatment check: if only swept breeds had been run
  # forward, drift would raise F_ST across the whole chromosome, not just at
  # the sweep. Compare background F_ST away from sweeps on sweep-carrying
  # chromosomes against chromosomes with no sweep at all.
  data(css_sim); data(css_sim_truth)
  d <- data.table::as.data.table(css_sim)
  tr <- data.table::as.data.table(css_sim_truth)

  d[, near_signal := FALSE]
  for (i in seq_len(nrow(tr))) {
    d[as.character(chr) == tr$chr[i] & abs(pos - tr$pos[i]) < 5e6,
      near_signal := TRUE]
  }
  sweep_chr <- unique(tr$chr)
  bg_on_sweep_chr <- d[as.character(chr) %in% sweep_chr & !near_signal, mean(fst, na.rm = TRUE)]
  bg_on_clean_chr <- d[!as.character(chr) %in% sweep_chr, mean(fst, na.rm = TRUE)]

  expect_equal(bg_on_sweep_chr, bg_on_clean_chr, tolerance = 0.25)
})

test_that("the null is calibrated on permuted data", {
  skip_on_cran()
  data(css_sim)
  set.seed(101)
  # Permuting each test independently destroys any shared signal while keeping
  # each test's marginal distribution intact.
  d <- data.table::copy(data.table::as.data.table(css_sim))
  d[, fst := sample(fst)][, xpehh := sample(xpehh)][, ddaf := sample(ddaf)]
  r <- suppressWarnings(suppressMessages(css(css_input(d, tests = tests3))))
  expect_equal(mean(r$css), 1 / log(10), tolerance = 0.03)
  expect_lt(max(r$css), 8)
})

test_that("shipped positions are whole base pairs", {
  # scrm reports positions on [0, 1]; scaling to bp without rounding produced
  # nonsense like 31501.53 bp and SNP ids to match. Regression guard.
  data(css_sim); data(css_sim_small); data(css_sim_truth)
  expect_true(all(css_sim$pos == floor(css_sim$pos)))
  expect_true(all(css_sim_small$pos == floor(css_sim_small$pos)))
  expect_true(all(css_sim_truth$pos == floor(css_sim_truth$pos)))
  expect_false(any(grepl("\\.", css_sim$snp)))
})

test_that("shipped statistics stay in their valid ranges", {
  data(css_sim)
  expect_true(all(css_sim$maf >= 0 & css_sim$maf <= 0.5))
  expect_true(all(css_sim$daf_selected  >= 0 & css_sim$daf_selected  <= 1))
  expect_true(all(css_sim$daf_reference >= 0 & css_sim$daf_reference <= 1))
  expect_true(all(css_sim$fst <= 1, na.rm = TRUE))
  expect_true(all(is.finite(css_sim$xpehh)))
  expect_true(all(is.finite(css_sim$ddaf)))
  expect_true(all(is.finite(css_sim$dsaf)))
})
