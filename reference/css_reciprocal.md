# Reciprocal (signed) composite selection signals

Randhawa et al. (2015) compare two cohorts in both directions: each
cohort in turn plays the selected population, with the other as
reference. The two runs are plotted mirrored about zero (their Figure
2), so that a signal's sign says which cohort carries the selection
signature.

## Usage

``` r
css_reciprocal(forward, reverse, labels = c("forward", "reverse"), ...)
```

## Arguments

- forward, reverse:

  `css_input` objects for the two directions.

- labels:

  Length-2 character vector naming the two cohorts, used by
  [`css_manhattan_mirror()`](https://msk99.github.io/cssig/reference/css_manhattan_mirror.md).
  Default `c("forward", "reverse")`.

- ...:

  Passed to [`css()`](https://msk99.github.io/cssig/reference/css.md).

## Value

A `css_result` keyed on `(chr, pos)` with columns `css_pos` (forward
direction), `css_neg` (reverse) and `css_signed`, the signed score
`css_pos - css_neg` used for mirrored plotting.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
draws
[`css_manhattan_mirror()`](https://msk99.github.io/cssig/reference/css_manhattan_mirror.md),
and
[`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md)
smooths both directions (`css_pos_smooth`, `css_neg_smooth`,
`css_signed_smooth`).

## Details

The two inputs must describe the same SNPs. `forward` should hold the
statistics computed with cohort A as the selected population and B as
reference; `reverse` the same tests with the roles swapped. For
directional tests such as XP-EHH and \\\Delta\\DAF this means
recomputing them; for symmetric tests such as \\F\_{ST}\\ the same
column is reused, which is correct — \\F\_{ST}\\ carries no direction,
and it is the directional tests that break the tie.

CSS is computed independently in each direction, so the ranks in one
direction do not depend on the other. Thresholds are likewise per
direction.

## Examples

``` r
data(css_sim_small)
fwd <- css_input(css_sim_small,
                 tests = c(fst = "high", xpehh = "high", ddaf = "high"))
# Reverse direction: negate the directional tests, keep F_ST as is.
rev_dat <- data.table::copy(css_sim_small)
rev_dat$xpehh <- -rev_dat$xpehh
rev_dat$ddaf  <- -rev_dat$ddaf
rvs <- css_input(rev_dat,
                 tests = c(fst = "high", xpehh = "high", ddaf = "high"))
res <- css_reciprocal(fwd, rvs, labels = c("selected", "reference"))
head(res)
#> <css_result> 6 SNPs, m = 3 constituent tests
#>   ties = "average", na_action = "pairwise"
#> 
#> Key: <chr, pos>
#>    chr    pos         snp    css_pos   css_neg  css_signed     p_pos      p_neg
#> 1:   1  31502  chr1_31502 0.10945766 1.3301010 -1.22064332 0.7772171 0.04676264
#> 2:   1 185282 chr1_185282 0.04336324 0.3451362 -0.30177293 0.9049754 0.45171429
#> 3:   1 265069 chr1_265069 0.07677106 1.3031801 -1.22640902 0.8379709 0.04975307
#> 4:   1 391046 chr1_391046 0.33743248 0.2817206  0.05571188 0.4597985 0.52273238
#> 5:   1 442469 chr1_442469 0.25330478 0.2742968 -0.02099199 0.5580784 0.53174478
#> 6:   1 577946 chr1_577946 0.77410004 0.1331530  0.64094700 0.1682287 0.73594772
```
