# Genome-wide Manhattan plot of composite selection signals

Reproduces the layout of Figure 3 of Randhawa et al. (2014): CSS against
genomic position, chromosomes in alternating colours, with the empirical
significance threshold drawn as a dashed line.

## Usage

``` r
css_manhattan(
  x,
  score = c("smoothed", "raw"),
  overlay_raw = FALSE,
  regions = NULL,
  label = NULL,
  thin = FALSE,
  thin_below = 1,
  thin_frac = 0.1,
  gap = 2e+07,
  point_size = 0.7,
  point_alpha = 0.85,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A `css_result`. If
  [`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md)
  has been run the threshold line is drawn automatically.

- score:

  Which score to plot: `"smoothed"` (default when available) or `"raw"`.

- overlay_raw:

  Draw the raw per-SNP scores behind the smoothed line, as in Figure 1
  of the 2014 paper. Only meaningful when `score = "smoothed"`.

- regions:

  Optional `css_regions` object; significant regions are highlighted
  and, if `label` is set, annotated.

- label:

  Optional column name in `regions` to use as a text label, for example
  a gene name column the user has added.

- thin:

  Drop a random subset of points below `thin_below` to keep rendering
  fast on very dense data. Off by default so that nothing is hidden
  without the user asking.

- thin_below:

  CSS value under which thinning applies. Default `1`.

- thin_frac:

  Fraction of low points to keep when thinning. Default `0.1`.

- gap:

  Gap between chromosomes in base pairs, passed to
  [`css_genome_coords()`](https://msk99.github.io/cssig/reference/css_genome_coords.md).

- point_size, point_alpha:

  Point aesthetics.

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
res <- css_threshold(css_smooth(res))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
css_manhattan(res)

```
