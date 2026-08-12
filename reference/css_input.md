# Assemble constituent selection statistics for a CSS analysis

Validates a table of pre-computed constituent selection-test statistics
and prepares it for
[`css()`](https://msk99.github.io/cssig/reference/css.md). One row per
SNP; one column per constituent test.

## Usage

``` r
css_input(
  data,
  chr = "chr",
  pos = "pos",
  snp = NULL,
  tests,
  keep_cols = NULL,
  drop_na_pos = TRUE
)
```

## Arguments

- data:

  A `data.frame`, `tibble`, `data.table` or matrix with one row per SNP.

- chr, pos:

  Column names holding the chromosome and the base-pair position.

- snp:

  Optional column name holding a SNP identifier. If `NULL`, a column
  named `snp` is used when present; failing that, an identifier of the
  form `"chr:pos"` is generated.

- tests:

  A named character vector mapping column names to directions, for
  example `c(fst = "high", xpehh = "high", ddaf = "high")`. Names are
  the columns of `data`; values are one of `"high"`, `"low"` or `"abs"`.

- keep_cols:

  Optional character vector of additional columns of `data` to carry
  through unchanged, for example MAF or gene annotation. They ride along
  untouched through every pipeline stage, and can be smoothed with the
  `cols` argument of
  [`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md).

- drop_na_pos:

  Drop rows with a missing chromosome or position. Default `TRUE`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of class `css_input`, keyed on `(chr, pos)`, with `chr` stored as an
ordered factor.

## Reference semantics

`css_input()` takes a defensive
[`data.table::copy()`](https://rdrr.io/pkg/data.table/man/copy.html) of
`data`, so the object you pass in is never modified. From that point on
the pipeline stages
([`css()`](https://msk99.github.io/cssig/reference/css.md),
[`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md),
[`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md),
[`css_fdr()`](https://msk99.github.io/cssig/reference/css_fdr.md)) add
columns *by reference* to the object `css_input()` returned, which
avoids copying large tables between stages. Each stage also returns the
object, so both the piped and the in-place idiom work. Pass
`.copy = TRUE` to a stage to get a copy instead. See
[`vignette("cssig")`](https://msk99.github.io/cssig/articles/cssig.md).

## Direction

CSS ranks every constituent test so that *larger means more evidence of
selection*. `tests` records, per test, how to get there:

- `"high"`:

  Use as-is. Correct for \\F\_{ST}\\, XP-EHH computed with the selected
  population as the target, and \\\Delta\\DAF computed as selected minus
  reference.

- `"low"`:

  Negate before ranking. For tests where small values indicate
  selection, such as Tajima's D or a nucleotide-diversity ratio.

- `"abs"`:

  Use the absolute value. For a signed test used without regard to
  direction.

## See also

[`css()`](https://msk99.github.io/cssig/reference/css.md),
[`css_merge_tests()`](https://msk99.github.io/cssig/reference/css_merge_tests.md)

## Examples

``` r
data(css_sim_small)
x <- css_input(css_sim_small,
               tests = c(fst = "high", xpehh = "high", ddaf = "high"))
x
#> <css_input> 5,502 SNPs on 2 chromosomes, 3 constituent tests
#>   fst          direction = high   missing = 0
#>   xpehh        direction = high   missing = 0
#>   ddaf         direction = high   missing = 0
#> 
#> Key: <chr, pos>
#>       chr       pos            snp          fst     xpehh      ddaf
#>    1:   1     31502     chr1_31502  0.033731900 -0.400611 -1.601780
#>    2:   1    185282    chr1_185282 -0.000442191 -0.789998 -0.258591
#>    3:   1    265069    chr1_265069  0.024843400 -0.905406 -1.168970
#>    4:   1    391046    chr1_391046  0.009940270 -0.602531  0.830885
#>    5:   1    442469    chr1_442469  0.005942830  0.781350 -0.631699
#>   ---                                                              
#> 5498:   2 136780586 chr2_136780586  0.003893630  0.461470 -0.601850
#> 5499:   2 136862246 chr2_136862246  0.002223530 -0.152700  0.218988
#> 5500:   2 136895797 chr2_136895797  0.073797600 -0.845844 -2.303220
#> 5501:   2 136937931 chr2_136937931  0.004824170  0.489281  0.607020
#> 5502:   2 137022913 chr2_137022913 -0.000698308  0.804152  0.144366
```
