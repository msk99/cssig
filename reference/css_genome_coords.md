# Build cumulative genome coordinates for a genome-wide plot

Converts per-chromosome positions to a single running coordinate, and
returns the axis tick positions. Exported because users assembling
custom plots need the same mapping the package's own Manhattan plots
use.

## Usage

``` r
css_genome_coords(chr, pos, gap = 2e+07)
```

## Arguments

- chr:

  Chromosome labels.

- pos:

  Positions within chromosome.

- gap:

  Gap in base pairs inserted between chromosomes. Default `2e7`.

## Value

A list with `pos_cum` (numeric vector aligned to the inputs), `axis` (a
`data.table` of chromosome label and tick midpoint) and `bounds`
(chromosome start offsets).

## Examples

``` r
g <- css_genome_coords(rep(1:2, each = 3), c(1e6, 2e6, 3e6, 1e6, 2e6, 3e6))
g$axis
#>       chr     mid      lo      hi
#>    <fctr>   <num>   <num>   <num>
#> 1:      1 2.0e+06 1.0e+06 3.0e+06
#> 2:      2 2.5e+07 2.4e+07 2.6e+07
```
