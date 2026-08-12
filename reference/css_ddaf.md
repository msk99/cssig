# Change in derived allele frequency between cohorts

\\\Delta\\DAF, the difference in derived allele frequency between the
putatively selected cohort and the reference cohort, as used by Grossman
et al. (2010) and as a CSS constituent in Randhawa et al. (2014).

## Usage

``` r
css_ddaf(daf_selected, daf_reference, standardize = FALSE)
```

## Arguments

- daf_selected, daf_reference:

  Derived allele frequency in each cohort.

- standardize:

  Rescale the result to mean 0 and unit variance, as both source papers
  do. Note this has no effect on CSS, which is rank-based; it matters
  only when plotting constituent tests on a shared axis.

## Value

Numeric vector of \\\Delta\\DAF values.

## Examples

``` r
css_ddaf(c(0.8, 0.2, 0.5), c(0.3, 0.25, 0.5))
#> [1]  0.50 -0.05  0.00
```
