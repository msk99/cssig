# Smooth CSS scores in sliding genomic windows

Averages a score over a window of fixed *physical* width centred on each
SNP, as in both source papers: Randhawa et al. (2014) describe 1 Mb
windows centred at each SNP, and Randhawa et al. (2015) describe
averaging 0.5 Mb on each side, which is the same window. Smoothing
suppresses single-SNP noise and exploits the fact that hitchhiking makes
true signals cluster.

## Usage

``` r
css_smooth(
  x,
  half_width = 5e+05,
  min_snps = 5L,
  on = c("css", "zbar"),
  cols = NULL,
  .copy = FALSE
)
```

## Arguments

- x:

  A `css_result` from
  [`css()`](https://msk99.github.io/cssig/reference/css.md) or
  [`css_reciprocal()`](https://msk99.github.io/cssig/reference/css_reciprocal.md),
  or any keyed table with `chr`, `pos` and the columns named in `cols`.

- half_width:

  Half the window width in base pairs. Default `5e5`, giving the 1 Mb
  window of the papers.

- min_snps:

  Minimum number of SNPs in a window for it to be retained. Default `5`.

- on:

  What to smooth. `"css"` (default) averages the \\-\log\_{10}(p)\\ CSS
  score, matching the published figures. `"zbar"` averages the mean z
  instead and re-derives a p-value from the smoothed value, which is not
  what the papers plot but is available for users who prefer it.
  `"zbar"` is refused for weighted and reciprocal results, whose
  statistic is not `sqrt(m) * zbar`; under `na_action = "pairwise"` its
  re-derived p-value uses the global `m`, an approximation for SNPs
  scored by fewer tests.

- cols:

  Optional character vector of additional columns to smooth, for example
  the constituent tests, for side-by-side comparison. Smoothed columns
  are suffixed `_smooth`.

- .copy:

  If `TRUE`, work on a copy and leave `x` untouched.

## Value

`x` with `css_smooth` and `n_window` added, plus one `<col>_smooth`
column for each entry of `cols`. For a reciprocal result the smoothed
columns are `css_pos_smooth`, `css_neg_smooth` and `css_signed_smooth`.

## Details

Windows never span a chromosome boundary. Windows containing fewer than
`min_snps` SNPs are returned as `NA` and are excluded from thresholding
and region calling; Randhawa et al. (2014) prune such low-density
windows for the same reason.

The implementation is vectorised with
[`findInterval()`](https://rdrr.io/r/base/findInterval.html) and
[`cumsum()`](https://rdrr.io/r/base/cumsum.html) per chromosome, so cost
does not grow with window width.
[`data.table::frollmean()`](https://rdrr.io/pkg/data.table/man/froll.html)
is deliberately not used: it windows by *row count*, whereas CSS windows
are defined by base-pair distance over irregularly spaced SNPs.

A reciprocal result from
[`css_reciprocal()`](https://msk99.github.io/cssig/reference/css_reciprocal.md)
carries two directed scores rather than one `css` column, so both are
smoothed: `css_pos_smooth`, `css_neg_smooth` and their difference
`css_signed_smooth` are added, and
[`css_manhattan_mirror()`](https://msk99.github.io/cssig/reference/css_manhattan_mirror.md)
then plots the smoothed scores by default. `on = "zbar"` is not
available in that case.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
res <- css_smooth(res)
summary(res$css_smooth)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.09917 0.25985 0.35899 0.42048 0.48566 2.24335 
```
