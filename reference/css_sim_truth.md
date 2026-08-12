# Known selection signatures in the simulated panel

The ground truth behind
[css_sim](https://msk99.github.io/cssig/reference/css_sim.md): where the
sweeps are, how strong they were, and which cohort carried them. Because
the answers are known, this supports genuine power and false-positive
comparisons of CSS against its constituent tests.

## Usage

``` r
css_sim_truth
```

## Format

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with one row per implanted signature:

- sweep_id:

  Identifier; `"trap"` marks the deliberate false positive.

- chr, pos:

  Location of the causal site.

- scenario:

  What kind of signature this is.

- s:

  Selection coefficient, additive with fitnesses \\1, 1+s/2, 1+s\\.

- p0:

  Starting frequency; `NA` means a new mutation at \\1/(2N_e)\\.

- target_freq:

  Frequency the sweep was designed to reach.

- cohort:

  Which breeds were under selection.

- selected:

  `FALSE` for the trap, `TRUE` for genuine sweeps.

- t_generations, ne_breed:

  Duration and effective size of the forward phase.

## Why the selection coefficients are large

The forward phase runs for only 40 generations, roughly 200 years of
breed formation. Under additive fitness an allele needs \\t =
(2/s)\log\[(p_1/(1-p_1))/(p_0/(1-p_0))\]\\ generations to travel from
\\p_0\\ to \\p_1\\, so a complete sweep from a single copy inside 40
generations requires \\s\\ near 0.5. These are values appropriate to
intense artificial selection, not to slow natural selection.

## The trap

Row `"trap"` carries no selection at all. In four selected-cohort breeds
the whole 40-45 Mb region of chromosome 5 descends from a single founder
haplotype, which is what happens when a breed is founded from few
animals. It produces elevated \\F\_{ST}\\ and long-range haplotype
homozygosity through drift alone. It is included so that users can see
what the method does when differentiation is real but selection is
absent, not because it is guaranteed to escape detection.

## Examples

``` r
data(css_sim_truth)
css_sim_truth
#>    sweep_id    chr      pos                                      scenario     s
#>      <char> <char>    <num>                                        <char> <num>
#> 1:        1      2 60000000                           complete hard sweep  0.50
#> 2:        2      1 95000000                         incomplete hard sweep  0.30
#> 3:        3      6 40000000            soft sweep from standing variation  0.20
#> 4:        4     14 30000000              weak sweep on standing variation  0.10
#> 5:        5     10 55000000                   breed-restricted hard sweep  0.45
#> 6:        6     13 45000000                 sweep in the reference cohort  0.45
#> 7:     trap      5 42500000 founder-effect haplotype block (no selection)  0.00
#>       p0 target_freq                   cohort selected t_generations ne_breed
#>    <num>       <num>                   <char>   <lgcl>         <int>    <int>
#> 1:    NA        0.99                 selected     TRUE            40      150
#> 2:    NA        0.70                 selected     TRUE            40      150
#> 3:  0.10        0.85                 selected     TRUE            40      150
#> 4:  0.15        0.60                 selected     TRUE            40      150
#> 5:    NA        0.98 selected (4 of 8 breeds)     TRUE            40      150
#> 6:    NA        0.98                reference     TRUE            40      150
#> 7:    NA          NA    selected (breeds 1-4)    FALSE            40      150
```
