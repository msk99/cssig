# Read XP-EHH output from selscan

Reads a selscan `.xpehh.out` file, or the normalised `.xpehh.out.norm`
produced by `norm`, into the shape
[`css_input()`](https://msk99.github.io/cssig/reference/css_input.md)
expects.

## Usage

``` r
read_selscan_xpehh(file, chr, normalised = NULL)
```

## Arguments

- file:

  Path to the selscan output file.

- chr:

  Chromosome label to attach. selscan output covers one chromosome and
  does not record which, so this must be supplied.

- normalised:

  Set `TRUE` when reading a `.norm` file, which carries the extra
  `normxpehh` and `crit` columns. If `NULL` (default) the format is
  detected from the header.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with `chr`, `pos`, `snp` and `xpehh`.

## Examples

``` r
f <- system.file("extdata", "example.xpehh.out", package = "cssig")
if (nzchar(f)) read_selscan_xpehh(f, chr = 1)
#> Key: <chr, pos>
#>         chr     pos    snp    xpehh
#>      <char>   <num> <char>    <num>
#>   1:      1  157717    rs1  1.13271
#>   2:      1  231071    rs2  0.33827
#>   3:      1  254249    rs3 -1.69323
#>   4:      1  262906    rs4  1.06877
#>   5:      1  328548    rs5  1.31290
#>  ---                               
#> 196:      1 4921130  rs196  0.54597
#> 197:      1 4931082  rs197 -1.45258
#> 198:      1 4948462  rs198  2.00715
#> 199:      1 4958766  rs199  0.55957
#> 200:      1 4980387  rs200  0.72937
```
