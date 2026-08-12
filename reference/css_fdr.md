# Estimate false discovery rates for CSS p-values

Estimate false discovery rates for CSS p-values

## Usage

``` r
css_fdr(
  x,
  method = c("BH", "BY", "fdrtool", "empirical-null", "isotonic"),
  force = FALSE,
  .copy = FALSE
)
```

## Arguments

- x:

  A `css_result` from
  [`css()`](https://msk99.github.io/cssig/reference/css.md).

- method:

  One of `"BH"`, `"BY"`, `"fdrtool"`, `"empirical-null"`, `"isotonic"`;
  see Details.

- force:

  Allow FDR estimation even if the object has been smoothed. Default
  `FALSE`.

- .copy:

  If `TRUE`, work on a copy and leave `x` untouched.

## Value

`x` with a `qval` column added (and `p_adj` for `"BH"` and `"BY"`).

## Details

FDR is computed from the **unsmoothed** CSS p-values. Smoothing averages
\\-\log\_{10}(p)\\ over neighbouring SNPs, and the result is no longer a
p-value, so feeding smoothed scores into an FDR procedure produces
numbers that look like q-values but are not. `css_fdr()` therefore
refuses smoothed input unless `force = TRUE`.

Five methods are available:

- `"BH"`:

  Benjamini-Hochberg via
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). No extra
  dependency; the default. Its FDR guarantee assumes positive-regression
  dependence, which genome-wide rank statistics with mixed-sign
  cross-test correlation do not obviously satisfy.

- `"BY"`:

  Benjamini-Yekutieli (2001) via
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Valid
  under arbitrary dependence, at a known cost in power; the conservative
  bound.

- `"fdrtool"`:

  Tail-area q-values from the fdrtool package with
  `statistic = "pvalue"`. This is the route used by Randhawa et al.
  (2014). Note the failure mode below: on correlated constituents it can
  estimate the null proportion as 1 and return q = 1 everywhere.

- `"empirical-null"`:

  Efron-style empirical null (Efron 2004): the *z statistic*
  \\\sqrt{m}\bar{Z}\\ is handed to `fdrtool(statistic = "normal")`,
  which estimates the null sd from the central bulk of the data instead
  of assuming the independence value 1. This directly repairs the
  `"fdrtool"` failure mode. It is two-sided, so the low tail – selection
  toward the reference cohort – is scored as well, which the upper-tail
  CSS p-value ignores. The fitted null sd and null proportion are stored
  in the `css_fdr` attribute.

- `"isotonic"`:

  Recalibrates the empirical p-value distribution with an isotonic
  regression of observed on expected quantiles before passing the result
  to fdrtool. Randhawa et al. (2014) use ConReg-R for this step. **This
  is a reimplementation in the same spirit, not ConReg-R**, and will not
  reproduce it exactly. The recalibration map is returned in the
  `css_fdr` attribute so it can be inspected.

## The CSS p-value is not calibrated on real data

\\p = 1 - \Phi(m^{1/2}\bar{Z})\\ assumes the constituent tests are
independent, so that \\\bar{Z}\\ has variance \\1/m\\. Real constituent
tests are correlated, and the correlation need not be positive: on the
shipped [css_sim](https://msk99.github.io/cssig/reference/css_sim.md)
data XP-EHH and \\\Delta\\DAF have Spearman correlation about -0.32,
which *deflates* the variance of \\\bar{Z}\\. The observed standard
deviation of \\m^{1/2}\bar{Z}\\ there is 0.90 rather than 1; permuting
each test independently restores it to 1.00.

The practical consequences:

- CSS p-values, and therefore CSS scores, are useful as a *ranking* but
  should not be read as calibrated tail probabilities. This is why both
  source papers threshold on the empirical top 0.1% of the genome-wide
  distribution rather than on a p-value, and why
  [`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md)
  does the same.

- FDR estimated from these p-values inherits the miscalibration. On
  [css_sim](https://msk99.github.io/cssig/reference/css_sim.md),
  `fdrtool` estimates the null proportion as 1 and returns q = 1
  everywhere, while `"BH"` returns a usable ordering. Two repairs are
  available: `css(calibrate = TRUE)` corrects the p-values themselves,
  and `method = "empirical-null"` fits the null scale when computing
  q-values. Check
  [`css_qq()`](https://msk99.github.io/cssig/reference/css_qq.md) before
  relying on any of them.

## What a per-SNP q-value can and cannot say

Even a perfectly calibrated per-SNP q-value answers "is this SNP's
cross-test rank concordance beyond global-null chance?", not "is this a
sweep?". A neutral genome violates the independent-ranks null *locally*:
shared genealogy produces real, selection-free concordance in
drift-driven blocks, of which the founder-effect trap in
[css_sim_truth](https://msk99.github.io/cssig/reference/css_sim_truth.md)
is the extreme case. On
[css_sim](https://msk99.github.io/cssig/reference/css_sim.md), the
majority of q \<= 0.05 SNPs under every method fall outside the
implanted sweep regions for exactly this reason. Sweep evidence remains
what the source papers say it is: clustered, smoothed signal called into
regions against an empirical threshold
([`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md),
[`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md),
[`css_regions()`](https://msk99.github.io/cssig/reference/css_regions.md)).
Treat q-values as a calibrated *descriptive* summary of the per-SNP
evidence, and note that permuting SNP labels cannot rescue region-level
error rates – permutation destroys the neutral autocorrelation that
generates false clusters, so it is anti-conservative for regions; an
honest region-level null needs neutral simulation or replication across
independent cohorts.

## References

Benjamini Y, Yekutieli D (2001). The control of the false discovery rate
in multiple testing under dependency. *Annals of Statistics*
29:1165-1188.

Efron B (2004). Large-scale simultaneous hypothesis testing: the choice
of a null hypothesis. *JASA* 99:96-104.

Strimmer K (2008). fdrtool: a versatile R package for estimating local
and tail area-based false discovery rates. *Bioinformatics*
24:1461-1462.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
res <- css_fdr(res, method = "BH")
#> FDR (BH): 30 SNPs with q <= 0.05, 13 with q <= 0.01.
sum(res$qval <= 0.05)
#> [1] 30
```
