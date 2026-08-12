# Mirrored Manhattan plot for reciprocal cohort contrasts

Reproduces Figure 2 of Randhawa et al. (2015): the two directions of a
reciprocal comparison drawn above and below zero, with a separate
empirical threshold for each.

## Usage

``` r
css_manhattan_mirror(
  x,
  score = c("smoothed", "raw"),
  top = 0.001,
  gap = 2e+07,
  point_size = 0.7,
  point_alpha = 0.85,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A `css_result` from
  [`css_reciprocal()`](https://msk99.github.io/cssig/reference/css_reciprocal.md).

- score:

  Which scores to plot: `"smoothed"` (the default once
  [`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md)
  has been run on the reciprocal result) or `"raw"`.

- top:

  Upper tail fraction used for the two threshold lines, applied to each
  direction's own distribution of the plotted score. Default `0.001`.

- gap, point_size, point_alpha, title, subtitle:

  As for
  [`css_manhattan()`](https://msk99.github.io/cssig/reference/css_manhattan.md).

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
data(css_sim_small)
fwd <- css_input(css_sim_small,
                 tests = c(fst = "high", xpehh = "high", ddaf = "high"))
rd <- data.table::copy(css_sim_small)
rd$xpehh <- -rd$xpehh; rd$ddaf <- -rd$ddaf
rvs <- css_input(rd, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
recip <- css_reciprocal(fwd, rvs, labels = c("large", "small"))
css_manhattan_mirror(recip)              # raw scores

css_manhattan_mirror(css_smooth(recip))  # smoothed, as the 2015 figure plots

```
