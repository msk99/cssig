# Call genomic regions under selection

Groups significant SNPs into contiguous regions. Two rules are provided,
one from each source paper; they differ enough to change results, so
both are available rather than silently blended.

## Usage

``` r
css_regions(
  x,
  method = c("cluster", "flank"),
  min_snps = 3L,
  min_flank = 5L,
  merge_gap = 1e+06,
  flank = 5e+05
)
```

## Arguments

- x:

  A `css_result` that has been through
  [`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md).

- method:

  `"cluster"` (default) or `"flank"`; see Details.

- min_snps:

  Minimum significant SNPs per cluster for `"cluster"`. Default `3`.

- min_flank:

  Minimum adjoining top-`top2` SNPs for `"flank"`. Default `5`.

- merge_gap:

  Maximum distance in base pairs between clusters that are merged into
  one region. Default `1e6`.

- flank:

  Padding in base pairs reported for gene mining. Default `5e5`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of class `css_regions`, one row per region, with the peak SNP, span, SNP
counts and, when the constituent tests are present, the number of SNPs
in each test's own top `top` fraction that fall inside the region.

## Details

- `method = "cluster"`:

  Randhawa et al. (2014). A cluster of at least `min_snps` significant
  SNPs (top 0.1% of smoothed CSS) spanning a window of `merge_gap`
  around the core SNP. Region boundaries are the first and last
  significant SNP. Clusters closer than `merge_gap` are merged.

- `method = "flank"`:

  Randhawa et al. (2015). At least one SNP in the top 0.1%, flanked by
  at least `min_flank` adjoining SNPs in the top 1%. Boundaries are the
  first and last top-1% SNP of the cluster. Clusters closer than
  `merge_gap` are merged.

"Adjoining" is implemented by distance: top-`top2` SNPs within
`merge_gap` of each other form one cluster even when SNPs below the
top-`top2` cut lie between them, which is slightly more permissive than
a literal reading of the 2015 rule. Both papers additionally extend
regions by 0.5 Mb on each side *for gene mining only*. That padding is
reported in `start_padded` / `end_padded` and never applied to `start` /
`end`, so region boundaries are not silently inflated.

## Examples

``` r
data(css_sim_small)
res <- css(css_input(css_sim_small,
                     tests = c(fst = "high", xpehh = "high", ddaf = "high")))
res <- css_threshold(css_smooth(res))
#> Threshold on smoothed CSS: top 0.1% at 2.1320 (6 SNPs), top 1.0% at 1.5565 (56 SNPs).
css_regions(res)
#> <css_regions> 1 region, method = "cluster"
#>   span: median 0.37 Mb, range 0.37-0.37 Mb; 6 significant SNPs total
#> 
#> Key: <chr>
#>    region_id chr    start      end start_padded end_padded n_snps n_sig
#> 1:         1   2 60083571 60458346     59583571   60958346      6     6
#>    peak_pos      peak_snp peak_css n_fst n_xpehh n_ddaf
#> 1: 60217074 chr2_60217074 2.243345     0       0      0
```
