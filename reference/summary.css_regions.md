# Summarise called regions

Prints the table in the shape of Table 2 of Randhawa et al. (2014): one
row per region with its span, peak and the number of significant SNPs
contributed by CSS and by each constituent test.

## Usage

``` r
# S3 method for class 'css_regions'
summary(object, ...)
```

## Arguments

- object:

  A `css_regions` object.

- ...:

  Unused.

## Value

Invisibly, a `data.table` of the printed summary.

## Examples

``` r
data(css_sim_small)
res <- css_threshold(css_smooth(css(css_input(css_sim_small,
         tests = c(fst = "high", xpehh = "high", ddaf = "high")))))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
summary(css_regions(res))
#> Key: <chr>
#>    region chr    position span_kb peak_Mb peak_css n_snps n_sig n_fst n_xpehh
#> 1:      1   2 60.08-60.46     375  60.217     2.24      6     6     0       0
#>    n_ddaf
#> 1:      0
#> 
#> Positions in Mb. n_* columns count SNPs in each test's own top fraction.
```
