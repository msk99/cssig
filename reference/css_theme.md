# A minimal theme for CSS plots

A minimal theme for CSS plots

## Usage

``` r
css_theme(base_size = 11, base_family = "")
```

## Arguments

- base_size:

  Base font size in points.

- base_family:

  Base font family.

## Value

A [ggplot2::theme](https://ggplot2.tidyverse.org/reference/theme.html)
object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + css_theme()
```
