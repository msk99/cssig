# Package index

## Input

Getting constituent statistics into shape

- [`css_input()`](https://msk99.github.io/cssig/reference/css_input.md)
  : Assemble constituent selection statistics for a CSS analysis
- [`css_merge_tests()`](https://msk99.github.io/cssig/reference/css_merge_tests.md)
  : Join per-test statistic tables into a single CSS input
- [`read_selscan_xpehh()`](https://msk99.github.io/cssig/reference/read_selscan_xpehh.md)
  : Read XP-EHH output from selscan
- [`read_rehh_xpehh()`](https://msk99.github.io/cssig/reference/read_rehh_xpehh.md)
  : Coerce rehh XP-EHH output for use with CSS

## The CSS pipeline

Score, smooth, threshold, call regions, estimate FDR

- [`css()`](https://msk99.github.io/cssig/reference/css.md) : Compute
  composite selection signals
- [`css_smooth()`](https://msk99.github.io/cssig/reference/css_smooth.md)
  : Smooth CSS scores in sliding genomic windows
- [`css_threshold()`](https://msk99.github.io/cssig/reference/css_threshold.md)
  : Apply an empirical significance threshold to CSS scores
- [`css_regions()`](https://msk99.github.io/cssig/reference/css_regions.md)
  : Call genomic regions under selection
- [`css_fdr()`](https://msk99.github.io/cssig/reference/css_fdr.md) :
  Estimate false discovery rates for CSS p-values
- [`css_reciprocal()`](https://msk99.github.io/cssig/reference/css_reciprocal.md)
  : Reciprocal (signed) composite selection signals

## Constituent statistics

Helpers for the constituents computable from allele frequencies

- [`css_fst()`](https://msk99.github.io/cssig/reference/css_fst.md) :
  Weir and Cockerham's F_ST between two cohorts
- [`css_ddaf()`](https://msk99.github.io/cssig/reference/css_ddaf.md) :
  Change in derived allele frequency between cohorts
- [`css_dsaf()`](https://msk99.github.io/cssig/reference/css_dsaf.md) :
  Change in selected allele frequency between cohorts
- [`css_standardize()`](https://msk99.github.io/cssig/reference/css_standardize.md)
  : Standardise a statistic to mean zero and unit variance

## Visualisation

- [`css_manhattan()`](https://msk99.github.io/cssig/reference/css_manhattan.md)
  : Genome-wide Manhattan plot of composite selection signals
- [`css_manhattan_mirror()`](https://msk99.github.io/cssig/reference/css_manhattan_mirror.md)
  : Mirrored Manhattan plot for reciprocal cohort contrasts
- [`css_chrom_plot()`](https://msk99.github.io/cssig/reference/css_chrom_plot.md)
  : Single-chromosome CSS plot
- [`css_region_plot()`](https://msk99.github.io/cssig/reference/css_region_plot.md)
  : Regional plot with constituent test tracks
- [`css_qq()`](https://msk99.github.io/cssig/reference/css_qq.md) : QQ
  plot of CSS p-values
- [`css_pdist()`](https://msk99.github.io/cssig/reference/css_pdist.md)
  : Histogram of CSS p-values
- [`css_fdr_density()`](https://msk99.github.io/cssig/reference/css_fdr_density.md)
  : Density of q-values inside called regions versus the rest of the
  genome
- [`css_test_cor()`](https://msk99.github.io/cssig/reference/css_test_cor.md)
  : Correlation among constituent tests and CSS
- [`css_circos()`](https://msk99.github.io/cssig/reference/css_circos.md)
  : Circular genome plot of CSS and its constituent tests
- [`css_theme()`](https://msk99.github.io/cssig/reference/css_theme.md)
  : A minimal theme for CSS plots
- [`css_genome_coords()`](https://msk99.github.io/cssig/reference/css_genome_coords.md)
  : Build cumulative genome coordinates for a genome-wide plot

## Data

- [`css_sim`](https://msk99.github.io/cssig/reference/css_sim.md) :
  Simulated multi-breed cattle panel with known selection signatures
- [`css_sim_small`](https://msk99.github.io/cssig/reference/css_sim_small.md)
  : Two-chromosome subset of the simulated panel
- [`css_sim_truth`](https://msk99.github.io/cssig/reference/css_sim_truth.md)
  : Known selection signatures in the simulated panel

## Methods

- [`summary(`*`<css_result>`*`)`](https://msk99.github.io/cssig/reference/summary.css_result.md)
  : Summarise a CSS result
- [`summary(`*`<css_regions>`*`)`](https://msk99.github.io/cssig/reference/summary.css_regions.md)
  : Summarise called regions
