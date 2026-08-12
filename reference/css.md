# Compute composite selection signals

The core CSS statistic of Randhawa et al. (2014, 2015). For each
constituent test the SNPs are ranked genome-wide, the ranks are rescaled
to fractional ranks, those are mapped to standard normal quantiles, and
the resulting z-scores are averaged across tests at each SNP. The mean z
is converted to an upper-tail p-value and reported as \\CSS =
-\log\_{10}(p)\\.

## Usage

``` r
css(
  x,
  ties = c("average", "first", "random", "dense"),
  na_action = c("pairwise", "omit"),
  weights = NULL,
  calibrate = FALSE,
  .copy = FALSE
)
```

## Arguments

- x:

  A `css_input` object from
  [`css_input()`](https://msk99.github.io/cssig/reference/css_input.md).

- ties:

  Tie-handling for the rank step, passed to
  [`data.table::frank()`](https://rdrr.io/pkg/data.table/man/frank.html):
  one of `"average"`, `"first"`, `"random"`, `"dense"`.

- na_action:

  `"pairwise"` (default) or `"omit"`; see Details.

- weights:

  Optional numeric vector of per-test weights, one per test. `NULL`
  (default) gives the equally weighted mean of the papers.

- calibrate:

  If `TRUE`, correct the p-value for the estimated cross-test
  correlation; see the Calibration section. Default `FALSE`, the
  published method.

- .copy:

  If `TRUE`, work on a copy and leave `x` untouched. Default `FALSE`,
  which adds the result columns to `x` by reference.

## Value

`x`, with columns `zbar`, `p` and `css` added, classed `css_result`.

## Details

With \\m\\ tests and \\n\\ SNPs, for test \\i\\ at SNP \\j\\:
\$\$R\_{ij} = \mathrm{rank}(T\_{ij}), \quad R'\_{ij} = R\_{ij}/(n+1),
\quad Z\_{ij} = \Phi^{-1}(R'\_{ij})\$\$ \$\$\bar{Z}\_j = m^{-1}\sum_i
Z\_{ij}, \quad p_j = 1 - \Phi\\\left(m^{1/2}\bar{Z}\_j\right), \quad
CSS_j = -\log\_{10}(p_j)\$\$ since \\\bar{Z}\_j \sim N(0, m^{-1})\\
under the null.

## CSS is rank-based

Only the *order* of each constituent test matters. Standardising XP-EHH
or \\\Delta\\DAF to \\N(0,1)\\, as both source papers do, leaves CSS
exactly unchanged; it matters only for plotting the constituent tests on
a common axis. Equally, any strictly increasing transformation of a
constituent test gives identical CSS.

## Deviations from the published method

Two situations the papers do not specify are handled explicitly, and
both defaults are flagged by the
[`print()`](https://rdrr.io/r/base/print.html) method when they are in
use:

- Ties:

  `ties = "average"` by default. Statistics computed from allele counts
  are discrete, so ties are normal and harmless. What is not harmless is
  a single value shared by a large block of SNPs, such as \\F\_{ST} =
  0\\ wherever a cohort is monomorphic: that block collapses to one
  averaged rank. A warning fires when any one value covers more than 5%
  of a test.

- Missing values:

  `na_action = "pairwise"` averages over the tests available at each SNP
  and uses that SNP's own \\m\\ in the p-value. This is an extension: a
  SNP scored by fewer tests is not strictly comparable to one scored by
  all of them. Use `"omit"` for the strict complete-case behaviour of
  the papers.

Non-`NULL` `weights` give a Stouffer-style weighted mean z, which is
*not* the published method; the default `NULL` reproduces the papers
exactly.

## Calibration

`p = 1 - \Phi(\sqrt{m}\bar{Z})` assumes the constituent tests are
independent. They are correlated in practice, which moves the true null
variance of \\\sqrt{m}\bar{Z}\\ away from 1 (on
[css_sim](https://msk99.github.io/cssig/reference/css_sim.md) the
observed sd is 0.90, explained to three decimals by the cross-test
correlations). `calibrate = TRUE` estimates the correlation matrix of
the per-test z columns genome-wide and divides the statistic by the
implied null sd, \\\sigma^2 = (w' R w) / \sum w^2\\, the Stouffer
analogue of Brown's correction for combining correlated tests (Brown
1975; Kost & McDermott 2002). Under `na_action = "pairwise"` the sd is
computed per missingness pattern. For complete data the correction is a
monotone rescaling, so the CSS *ranking* and the empirical thresholds of
[`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md)
are unchanged; only `p`, `css` and downstream q-values move. It is off
by default because it is not the published method, and
[`print()`](https://rdrr.io/r/base/print.html) flags a calibrated run.

## References

Randhawa IAS, Khatkar MS, Thomson PC, Raadsma HW (2014). Composite
selection signals can localize the trait specific genomic regions in
multi-breed populations of cattle and sheep. *BMC Genetics* 15:34.
[doi:10.1186/1471-2156-15-34](https://doi.org/10.1186/1471-2156-15-34)

Randhawa IAS, Khatkar MS, Thomson PC, Raadsma HW (2015). Composite
selection signals for complex traits exemplified through bovine stature
using multibreed cohorts of European and African *Bos taurus*. *G3*
5:1391-1401.
[doi:10.1534/g3.115.017772](https://doi.org/10.1534/g3.115.017772)

For `calibrate = TRUE`: Brown MB (1975). A method for combining
non-independent, one-sided tests of significance. *Biometrics*
31:987-992. Kost JT, McDermott MP (2002). Combining dependent p-values.
*Statistics & Probability Letters* 60:183-190.

## See also

[`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md),
[`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md),
[`css_regions()`](https://msk99.github.io/cssig/reference/css_regions.md)

## Examples

``` r
data(css_sim_small)
x <- css_input(css_sim_small,
               tests = c(fst = "high", xpehh = "high", ddaf = "high"))
res <- css(x)
head(res)
#> <css_result> 6 SNPs, m = 3 constituent tests
#>   ties = "average", na_action = "pairwise"
#> 
#> Key: <chr, pos>
#>    chr    pos         snp          fst     xpehh      ddaf        zbar
#> 1:   1  31502  chr1_31502  0.033731900 -0.400611 -1.601780 -0.44041910
#> 2:   1 185282 chr1_185282 -0.000442191 -0.789998 -0.258591 -0.75657906
#> 3:   1 265069 chr1_265069  0.024843400 -0.905406 -1.168970 -0.56935553
#> 4:   1 391046 chr1_391046  0.009940270 -0.602531  0.830885  0.05827858
#> 5:   1 442469 chr1_442469  0.005942830  0.781350 -0.631699 -0.08435034
#> 6:   1 577946 chr1_577946  0.015554700  0.400111  1.054750  0.55494255
#>            p        css
#> 1: 0.7772171 0.10945766
#> 2: 0.9049754 0.04336324
#> 3: 0.8379709 0.07677106
#> 4: 0.4597985 0.33743248
#> 5: 0.5580784 0.25330478
#> 6: 0.1682287 0.77410004
```
