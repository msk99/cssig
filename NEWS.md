# cssig 1.1.0

## Fixed

* Help pages now render markdown properly (code formatting and working
  cross-references). Documentation was previously generated without roxygen2's
  markdown mode, so every page displayed literal backticks and `[function()]`
  link syntax.
* `css_input()` no longer corrupts factor columns. `as.numeric()` on a factor
  returns the internal level codes, so a factor `pos` or test column was
  silently replaced with garbage; factors are now coerced through their labels.
* `css_region_plot()` draws the gene track in a dedicated bottom panel. It was
  drawn into every facet, dragging each panel's free y-scale down to the
  track's baseline. The `name` column, documented but previously ignored, now
  labels the genes.
* `css_region_plot()` accepts a named vector for `region`, as documented; it
  previously failed with "$ operator is invalid for atomic vectors".
* `plot()` on a `css_reciprocal()` result dispatches to
  `css_manhattan_mirror()` instead of failing with a misleading "Run `css()`
  first" error, and `css_manhattan()` redirects reciprocal input there.
* `?css_sim` now states the correct number of animals sampled per breed
  (50, not 20; the README's figure of 800 animals was already correct).
* Removed a duplicated "Sweep 6" section from `vignette("simulation")`.

## New

* `css_smooth()` supports reciprocal results, adding `css_pos_smooth`,
  `css_neg_smooth` and `css_signed_smooth`.
* `css_manhattan_mirror()` gained a `score` argument and, like
  `css_manhattan()`, plots the smoothed scores by default once they exist.

# cssig 1.0.0

* First release.
