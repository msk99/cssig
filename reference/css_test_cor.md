# Correlation among constituent tests and CSS

Randhawa et al. (2014) report low to moderate correlation between pairs
of constituent tests but high correlation between CSS and each of them,
which is the evidence that CSS captures information from across the
tests rather than tracking any single one.

## Usage

``` r
css_test_cor(
  x,
  method = c("spearman", "pearson", "kendall"),
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A `css_result` still carrying its constituent test columns.

- method:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). Default
  `"spearman"`, which is the appropriate choice here because CSS is
  itself a rank statistic.

- title, subtitle:

  Plot labels.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The underlying correlation matrix is attached as the `"cor"`
attribute.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
css_test_cor(res)

```
