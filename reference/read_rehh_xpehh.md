# Coerce rehh XP-EHH output for use with CSS

Takes the data frame returned by
[`rehh::ies2xpehh()`](https://rdrr.io/pkg/rehh/man/ies2xpehh.html) and
renames its columns to the shape
[`css_input()`](https://msk99.github.io/cssig/reference/css_input.md)
expects.

## Usage

``` r
read_rehh_xpehh(x, value = c("XPEHH", "LOGPVALUE"))
```

## Arguments

- x:

  A data frame from
  [`rehh::ies2xpehh()`](https://rdrr.io/pkg/rehh/man/ies2xpehh.html).

- value:

  Which column to use: `"XPEHH"` (default) picks the log-ratio column,
  whose exact name depends on the population labels and on whether
  `ies2xpehh(standardize = )` was `TRUE` (`XPEHH_a_b`) or `FALSE`
  (`UNXPEHH_a_b`); both spellings are recognised. `"LOGPVALUE"` uses
  rehh's -log10(p) column instead.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with `chr`, `pos`, `snp` and `xpehh`.

## Examples

``` r
fake <- data.frame(CHR = 1, POSITION = c(100, 200),
                   XPEHH_A_B = c(1.2, -0.4))
read_rehh_xpehh(fake)
#> Key: <chr, pos>
#>       chr   pos    snp xpehh
#>    <char> <num> <char> <num>
#> 1:      1   100   <NA>   1.2
#> 2:      1   200   <NA>  -0.4
```
