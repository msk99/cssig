# Join per-test statistic tables into a single CSS input

Constituent statistics usually come from different tools, on SNP sets
that do not fully agree. This joins them on `(chr, pos)` and reports
what each test contributed, so that a silent loss of SNPs at the merge
step becomes visible.

## Usage

``` r
css_merge_tests(
  ...,
  chr = "chr",
  pos = "pos",
  value = NULL,
  all = TRUE,
  quiet = FALSE
)
```

## Arguments

- ...:

  Two or more `data.frame`s or `data.table`s, ideally named
  (`fst = fst_table, xpehh = xpehh_table`). Each must contain the
  chromosome and position columns plus exactly one statistic column,
  unless `value` is given.

- chr, pos:

  Column names for chromosome and position, recycled across all inputs.

- value:

  Optional character vector naming the statistic column in each input.
  If `NULL`, each input must have exactly one column beyond `chr`, `pos`
  and `snp`.

- all:

  If `TRUE` (default) keep the union of SNPs, leaving `NA` where a test
  is missing. If `FALSE` keep only SNPs scored by every test, which is
  the behaviour assumed in the source papers.

- quiet:

  Suppress the merge report.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
suitable for
[`css_input()`](https://msk99.github.io/cssig/reference/css_input.md).

## Examples

``` r
a <- data.frame(chr = 1, pos = c(100, 200, 300), fst   = c(0.1, 0.2, 0.3))
b <- data.frame(chr = 1, pos = c(200, 300, 400), xpehh = c(1.1, 2.2, 0.4))
css_merge_tests(fst = a, xpehh = b, all = FALSE, quiet = TRUE)
#> Key: <chr, pos>
#>       chr   pos   fst xpehh
#>    <char> <num> <num> <num>
#> 1:      1   200   0.2   1.1
#> 2:      1   300   0.3   2.2
```
