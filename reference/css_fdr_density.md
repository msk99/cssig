# Density of q-values inside called regions versus the rest of the genome

Reproduces Figure 2 of Randhawa et al. (2014). A clear separation
between the two densities is the evidence that called regions carry a
real excess of selection signal rather than the tail of the genome-wide
distribution.

## Usage

``` r
css_fdr_density(x, regions, fdr = 0.05, title = NULL, subtitle = NULL)
```

## Arguments

- x:

  A `css_result` that has been through
  [`css_fdr()`](https://msk99.github.io/cssig/reference/css_fdr.md).

- regions:

  A `css_regions` object.

- fdr:

  FDR level to mark with a vertical line. Default `0.05`.

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
res <- css_fdr(res, method = "BH")
#> FDR (BH): 30 SNPs with q <= 0.05, 13 with q <= 0.01.
res <- css_threshold(css_smooth(res))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
reg <- css_regions(res)
if (nrow(reg)) css_fdr_density(res, reg)

```
