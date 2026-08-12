# QQ plot of CSS p-values

Under the null the CSS p-values should be uniform, so the points follow
the diagonal. Departure at the tail is the signal; departure along the
whole line indicates the null is miscalibrated, usually because the
constituent tests are strongly correlated.

## Usage

``` r
css_qq(x, max_points = 50000L, title = NULL, subtitle = NULL)
```

## Arguments

- x:

  A `css_result`.

- max_points:

  Thin to at most this many points for rendering. The smallest p-values
  are always kept. Default `50000`.

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
css_qq(res)

```
