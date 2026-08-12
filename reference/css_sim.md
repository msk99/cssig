# Simulated multi-breed cattle panel with known selection signatures

**These data are simulated. They describe no real animal, breed or
genome.**

## Usage

``` r
css_sim
```

## Format

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with one row per SNP:

- chr:

  Chromosome, an ordered factor with the 29 bovine autosomes.

- pos:

  Position in base pairs on a UMD3.1-shaped coordinate system.

- snp:

  SNP identifier.

- ancestral, derived:

  Allele coding. The ancestral allele is known exactly, because the
  simulation records it.

- maf:

  Minor allele frequency across the whole panel.

- daf_selected, daf_reference:

  Derived allele frequency in each cohort.

- fst:

  Weir & Cockerham \\F\_{ST}\\ between the cohorts, from
  [`css_fst()`](https://msk99.github.io/cssig/reference/css_fst.md).

- xpehh:

  Cross-population extended haplotype homozygosity, selected against
  reference, computed with rehh and standardised.

- ddaf:

  Standardised \\\Delta\\DAF, from
  [`css_ddaf()`](https://msk99.github.io/cssig/reference/css_ddaf.md).

- dsaf:

  Standardised \\\Delta\\SAF, from
  [`css_dsaf()`](https://msk99.github.io/cssig/reference/css_dsaf.md).
  Provided so that the ancestral-allele-free route of Randhawa et
  al. (2014) can be compared against \\\Delta\\DAF on the same data.

## Details

Constituent selection statistics for a simulated multi-breed cattle
panel, computed from simulated genotypes rather than drawn directly from
a distribution. That distinction matters: \\F\_{ST}\\, XP-EHH and
\\\Delta\\DAF are correlated with each other because they read the same
underlying genealogy, and CSS exists precisely to exploit that shared
structure. Statistics simulated independently would carry a correlation
structure chosen by the package author, and any demonstration built on
them would be circular.

## How it was generated

Sixteen breeds descend from a common ancestral population that passes
through a domestication bottleneck, simulated with the
sequentially-Markovian coalescent (scrm). Breeds are then formed by a
forward Wright-Fisher simulation with recombination, run for 40
generations at an effective size of 150, during which the sweeps in
[css_sim_truth](https://msk99.github.io/cssig/reference/css_sim_truth.md)
act. Eight breeds form the selected cohort and eight the reference
cohort, with 50 animals sampled per breed. SNPs are ascertained
chip-style, by discovery in three breeds followed by thinning to roughly
uniform spacing, which reproduces the shifted frequency spectrum of the
BovineSNP50 data the source papers analyse.

Every breed is run forward for the same number of generations at the
same effective size, whether or not it carries a sweep. Without that,
the swept breeds would accumulate extra drift and \\F\_{ST}\\ would be
raised across the whole region rather than at the sweep.

The causal variants themselves are excluded from the panel, as they
would be on a real SNP chip; the signal comes from hitchhiking
neighbours.

Full parameters and the generating scripts are in `data-raw/` in the
package source.

## See also

[css_sim_truth](https://msk99.github.io/cssig/reference/css_sim_truth.md)
for the known answers,
[css_sim_small](https://msk99.github.io/cssig/reference/css_sim_small.md)
for a subset sized for quick examples.

## Examples

``` r
data(css_sim)
str(css_sim)
#> Classes ‘data.table’ and 'data.frame':   46976 obs. of  12 variables:
#>  $ chr          : Factor w/ 29 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ pos          : num  31502 185282 265069 391046 442469 ...
#>  $ snp          : chr  "chr1_31502" "chr1_185282" "chr1_265069" "chr1_391046" ...
#>  $ ancestral    : chr  "A" "A" "A" "A" ...
#>  $ derived      : chr  "G" "G" "G" "G" ...
#>  $ maf          : num  0.444 0.444 0.227 0.344 0.239 ...
#>  $ daf_selected : num  0.378 0.434 0.725 0.38 0.214 ...
#>  $ daf_reference: num  0.51 0.454 0.821 0.309 0.265 ...
#>  $ fst          : num  0.033732 -0.000442 0.024843 0.00994 0.005943 ...
#>  $ xpehh        : num  -0.401 -0.79 -0.905 -0.603 0.781 ...
#>  $ ddaf         : num  -1.602 -0.259 -1.169 0.831 -0.632 ...
#>  $ dsaf         : num  1.517 0.171 -1.22 -0.921 0.545 ...
#>  - attr(*, ".internal.selfref")=<pointer: (nil)> 
```
