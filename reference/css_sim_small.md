# Two-chromosome subset of the simulated panel

Chromosomes 1 and 2 of
[css_sim](https://msk99.github.io/cssig/reference/css_sim.md), carrying
the incomplete and complete hard sweeps respectively. Sized so that
examples and tests run quickly.

## Usage

``` r
css_sim_small
```

## Format

As [css_sim](https://msk99.github.io/cssig/reference/css_sim.md),
restricted to two chromosomes.

## Examples

``` r
data(css_sim_small)
table(css_sim_small$chr)
#> 
#>    1    2 
#> 2962 2540 
```
