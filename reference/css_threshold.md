# Apply an empirical significance threshold to CSS scores

Both source papers declare SNPs significant when their score falls in
the top 0.1% of the genome-wide empirical distribution, rather than
referring to any theoretical null. This computes that threshold and
flags the SNPs above it.

## Usage

``` r
css_threshold(
  x,
  top = 0.001,
  top2 = 0.01,
  on = c("smoothed", "raw"),
  .copy = FALSE
)
```

## Arguments

- x:

  A `css_result`, normally after
  [`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md).

- top:

  Upper tail fraction to declare significant. Default `0.001` (top
  0.1%), as used in both papers.

- top2:

  A second, more permissive fraction used by the `"flank"` region caller
  of
  [`css_regions()`](https://msk99.github.io/cssig/reference/css_regions.md).
  Default `0.01` (top 1%).

- on:

  Which score to threshold: `"smoothed"` (default) or `"raw"`.

- .copy:

  If `TRUE`, work on a copy and leave `x` untouched.

## Value

`x` with logical columns `significant` and `significant2` added. The
numeric cut-offs are stored in the `css_threshold` attribute.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
res <- css_threshold(css_smooth(res))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
attr(res, "css_threshold")
#> $top
#> [1] 0.001
#> 
#> $top2
#> [1] 0.01
#> 
#> $on
#> [1] "smoothed"
#> 
#> $score_col
#> [1] "css_smooth"
#> 
#> $cut
#> [1] 2.131988
#> 
#> $cut2
#> [1] 1.55651
#> 
#> $n_scored
#> [1] 5502
#> 
#> $n_significant
#> [1] 6
#> 
#> $n_significant2
#> [1] 56
#> 
```
