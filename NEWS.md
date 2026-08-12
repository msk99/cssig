# cssig 1.2.0

The published method remains the default throughout: every addition below is
opt-in, and `print()` flags runs that deviate from the papers.

## New

* `css(calibrate = TRUE)` corrects the CSS p-value for cross-test correlation:
  a Stouffer-style variance correction (Brown 1975; Kost & McDermott 2002)
  with the correlation matrix estimated on the z scale, computed per
  missingness pattern under pairwise NA handling. For complete data it is a
  monotone rescaling, so the CSS ranking and empirical thresholds are
  unchanged; only `p`, `css` and downstream q-values move. Off by default.
* `css_fdr()` gains `method = "empirical-null"` — Efron-style q-values from
  fdrtool's normal mode fitted to the z statistic, which repairs the
  documented failure where the p-value route estimates a null proportion of 1
  and returns q = 1 everywhere — and `method = "BY"`, valid under arbitrary
  dependence. The default remains `"BH"`.
* `css_input()` gains `keep_cols` to carry annotation columns (MAF, gene
  names) through the pipeline.
* `?css_fdr` now has a section on what a per-SNP q-value can and cannot say
  about sweeps, and why permutation cannot provide region-level error rates.

## Fixed

* `read_selscan_xpehh(normalised = TRUE)` errors when the file has no
  `normxpehh` column instead of silently using the unnormalised score.
* `n_tests`, written by `css()` under pairwise NA handling, is now a reserved
  column name in `css_input()`.
* Duplicate names in `tests` are rejected with a clear message instead of an
  internal `setnames()` error.
* `css_manhattan()` no longer builds a zero-length subtitle when the stored
  threshold belongs to the other score; it says so instead.
* `css_smooth(on = "zbar")` refuses weighted results, whose statistic is not
  `sqrt(m) * zbar`, and documents the pairwise-m approximation.
* `css_fst()` rejects allele counts outside `[0, ploidy * n]`, catching
  swapped count/sample-size arguments.
* `graphics` declared in Imports; stale `globalVariables()` entries pruned;
  dead attribute check removed from the chromosome-factor helper.
* data-raw: the XP-EHH column regex now truly admits both rehh spellings;
  the recombination rate is passed from the parameter list instead of being
  hardcoded; a stale comment about a dropped dataset removed. Shipped data
  unchanged.
* README: dataset SNP counts corrected to the shipped data (46,976 and
  5,502).

## Infrastructure

* pkgdown workflow added (deploys to msk99.github.io/cssig; enable GitHub
  Pages on the `gh-pages` branch once); site URL added to DESCRIPTION.
* `workflow_dispatch` trigger added to R-CMD-check.
* Tests added for `css_fdr_density()` and `css_circos()`; the unused vdiffr
  dropped from Suggests.
* `?css_regions` documents that the flank rule groups top-1% SNPs by
  distance, slightly more permissive than a literal reading of the 2015 rule.

# cssig 1.1.0

## Fixed

* Help pages now render markdown properly (code formatting and working
  cross-references). Documentation was previously generated without roxygen2's
  markdown mode, so every page displayed literal backticks and `[function()]`
  link syntax.
* `css_input()` no longer corrupts factor columns. `as.numeric()` on a factor
  returns the internal level codes, so a factor `pos` or test column was
  silently replaced with garbage; factors are now coerced through their labels.
* `css_region_plot()` draws the gene track in a dedicated bottom panel. It was
  drawn into every facet, dragging each panel's free y-scale down to the
  track's baseline. The `name` column, documented but previously ignored, now
  labels the genes.
* `css_region_plot()` accepts a named vector for `region`, as documented; it
  previously failed with "$ operator is invalid for atomic vectors".
* `plot()` on a `css_reciprocal()` result dispatches to
  `css_manhattan_mirror()` instead of failing with a misleading "Run `css()`
  first" error, and `css_manhattan()` redirects reciprocal input there.
* `?css_sim` now states the correct number of animals sampled per breed
  (50, not 20; the README's figure of 800 animals was already correct).
* Removed a duplicated "Sweep 6" section from `vignette("simulation")`.

## New

* `css_smooth()` supports reciprocal results, adding `css_pos_smooth`,
  `css_neg_smooth` and `css_signed_smooth`.
* `css_manhattan_mirror()` gained a `score` argument and, like
  `css_manhattan()`, plots the smoothed scores by default once they exist.

# cssig 1.0.0

* First release.
