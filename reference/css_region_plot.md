# Regional plot with constituent test tracks

Zooms into one region and shows each constituent selection test
alongside CSS, which is how the source papers diagnose whether a
composite signal is driven by one test or supported by several.

## Usage

``` r
css_region_plot(
  x,
  region,
  pad = 5e+05,
  tests = NULL,
  smooth = TRUE,
  genes = NULL,
  units = c("Mb", "kb", "bp"),
  title = NULL
)
```

## Arguments

- x:

  A `css_result` that still carries its constituent test columns.

- region:

  Either a single row of a `css_regions` object, or a list/vector with
  `chr`, `start` and `end`.

- pad:

  Extra base pairs shown on each side. Default `5e5`.

- tests:

  Optional character vector selecting which constituent tests to show.
  Defaults to all of them.

- smooth:

  Draw the smoothed track for each test where available. Default `TRUE`.

- genes:

  Optional `data.frame` with `start`, `end` and `name` columns, drawn as
  a gene track in its own bottom panel, with `name` as the label.

- units:

  Axis units: `"Mb"` (default), `"kb"` or `"bp"`.

- title:

  Plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object, faceted by test with CSS on top.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
res <- css_threshold(css_smooth(res))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
reg <- css_regions(res)
if (nrow(reg)) css_region_plot(res, reg[1])

```
