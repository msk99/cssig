# Weir and Cockerham's F_ST between two cohorts

Per-SNP \\\theta\\ of Weir & Cockerham (1984) for two populations, the
estimator used for the \\F\_{ST}\\ constituent of CSS.

## Usage

``` r
css_fst(count1, n1, count2, n2, ploidy = 2L, floor_zero = FALSE)
```

## Arguments

- count1, count2:

  Integer vector of copies of the reference allele in cohort 1 and
  cohort 2.

- n1, n2:

  Number of *individuals* genotyped in each cohort, per SNP (scalars are
  recycled).

- ploidy:

  Ploidy, default `2`.

- floor_zero:

  Truncate negative estimates at zero. Default `FALSE`.

## Value

Numeric vector of \\F\_{ST}\\ estimates, `NA` where a SNP is monomorphic
across both cohorts.

## Details

Inputs are allele counts, not frequencies, because the estimator needs
the sample sizes: it corrects for the finite number of individuals
sampled, which a frequency-only calculation cannot do. Supply the count
of one allele and the number of genotyped individuals in each cohort.

Negative values are a normal outcome of the estimator when true
differentiation is near zero. They are returned as-is by default; set
`floor_zero = TRUE` to truncate at zero. Truncation biases the
genome-wide distribution upward, but because CSS is rank-based and
truncation is monotone it has no effect on CSS itself.

## References

Weir BS, Cockerham CC (1984). Estimating F-statistics for the analysis
of population structure. *Evolution* 38:1358-1370.

## Examples

``` r
css_fst(count1 = c(10, 40, 25), n1 = 25,
        count2 = c(40, 8, 26),  n2 = 25)
#> [1]  0.51960784  0.57310606 -0.02001667
```
