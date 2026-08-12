# Histogram of CSS p-values

Histogram of CSS p-values

## Usage

``` r
css_pdist(x, bins = 50L, title = NULL, subtitle = NULL)
```

## Arguments

- x:

  A `css_result`.

- bins:

  Number of histogram bins. Default `50`.

- title, subtitle:

  Plot labels.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
css_pdist(res)
```
