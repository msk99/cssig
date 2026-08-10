#' A minimal theme for CSS plots
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @return A [ggplot2::theme] object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + css_theme()
#' @export
css_theme <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(fill = NA, colour = "grey70"),
      axis.ticks         = ggplot2::element_line(colour = "grey70"),
      strip.text         = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      plot.title          = ggplot2::element_text(face = "bold"),
      plot.subtitle       = ggplot2::element_text(colour = "grey30"),
      plot.caption        = ggplot2::element_text(colour = "grey45", size = ggplot2::rel(0.8)),
      legend.position     = "bottom",
      legend.key.height   = ggplot2::unit(0.8, "lines")
    )
}

# Alternating chromosome colours. Colour-blind safe by construction: the two
# tones differ in lightness as well as hue, so the plot still reads in
# greyscale and under the common forms of colour vision deficiency.
.css_chrom_colours <- c("#2C3E6B", "#E08A1E")
.css_sig_colour    <- "#C0392B"

#' Build cumulative genome coordinates for a genome-wide plot
#'
#' Converts per-chromosome positions to a single running coordinate, and returns
#' the axis tick positions. Exported because users assembling custom plots need
#' the same mapping the package's own Manhattan plots use.
#'
#' @param chr Chromosome labels.
#' @param pos Positions within chromosome.
#' @param gap Gap in base pairs inserted between chromosomes. Default `2e7`.
#'
#' @return A list with `pos_cum` (numeric vector aligned to the inputs),
#'   `axis` (a `data.table` of chromosome label and tick midpoint) and
#'   `bounds` (chromosome start offsets).
#'
#' @examples
#' g <- css_genome_coords(rep(1:2, each = 3), c(1e6, 2e6, 3e6, 1e6, 2e6, 3e6))
#' g$axis
#'
#' @export
css_genome_coords <- function(chr, pos, gap = 2e7) {
  f <- .as_chrom_factor(chr)
  d <- data.table::data.table(chr = f, pos = as.numeric(pos),
                              .row = seq_along(pos))
  spans <- d[, .(len = max(pos, na.rm = TRUE)), by = chr]
  data.table::setorderv(spans, "chr")
  spans[, offset := c(0, cumsum(utils::head(len, -1L) + gap))]

  d <- merge(d, spans[, .(chr, offset)], by = "chr", sort = FALSE)
  data.table::setorderv(d, ".row")
  pos_cum <- d$pos + d$offset

  axis <- merge(spans, d[, .(lo = min(pos + offset), hi = max(pos + offset)), by = chr],
                by = "chr", sort = FALSE)
  axis[, mid := (lo + hi) / 2]

  list(pos_cum = pos_cum,
       axis = axis[, .(chr, mid, lo, hi)],
       bounds = spans[, .(chr, len, offset)])
}
