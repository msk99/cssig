# Standardise a statistic to mean zero and unit variance

Both source papers standardise XP-EHH and \\\Delta\\DAF this way. It
leaves CSS unchanged, because CSS depends only on ranks; it is done so
that constituent tests share an axis when plotted.

## Usage

``` r
css_standardize(x, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector.

- na.rm:

  Ignore missing values when computing the mean and standard deviation.
  Default `TRUE`.

## Value

Numeric vector of z-scores.

## Examples

``` r
css_standardize(c(1, 4, 9, 16))
#> [1] -0.9912407 -0.5337450  0.2287479  1.2962378
```
