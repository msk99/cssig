# Summarise a CSS result

Summarise a CSS result

## Usage

``` r
# S3 method for class 'css_result'
summary(object, ...)
```

## Arguments

- object:

  A `css_result`.

- ...:

  Unused.

## Value

Invisibly, a list of summary components.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
summary(res)
#> Composite selection signals: 5,502 SNPs on 2 chromosomes, m = 3 tests
#> 
#> CSS distribution:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.002005 0.133424 0.272591 0.420046 0.510309 5.871482 
#> 
#> Spearman correlation (constituents and CSS):
#>         fst  xpehh   ddaf   css
#> fst   1.000  0.014  0.019 0.647
#> xpehh 0.014  1.000 -0.308 0.457
#> ddaf  0.019 -0.308  1.000 0.395
#> css   0.647  0.457  0.395 1.000
```
