# Change in selected allele frequency between cohorts

\\\Delta\\SAF, the ancestral-allele-free alternative to \\\Delta\\DAF
introduced by Randhawa et al. (2014) for datasets where the ancestral
allele cannot be inferred, such as their sheep panels. It is the
frequency difference for the allele that is the *major* allele in the
selected cohort.

## Usage

``` r
css_dsaf(freq_selected, freq_reference, standardize = FALSE)
```

## Arguments

- freq_selected, freq_reference:

  Frequency of the same reference allele in each cohort.

- standardize:

  Rescale to mean 0 and unit variance. See
  [`css_ddaf()`](https://msk99.github.io/cssig/reference/css_ddaf.md).

## Value

Numeric vector of \\\Delta\\SAF values.

## Examples

``` r
css_dsaf(c(0.8, 0.2, 0.5), c(0.3, 0.25, 0.5))
#> [1] 0.50 0.05 0.00
```
