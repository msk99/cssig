#' Single-chromosome CSS plot
#'
#' Reproduces Figure 1 of Randhawa et al. (2014): raw CSS scores as points with
#' the smoothed score overlaid as a line, the empirical threshold as a dashed
#' line, and candidate gene positions marked with vertical rules.
#'
#' @param x A `css_result`.
#' @param chr Chromosome to plot. A single value matching a level of `x$chr`.
#' @param highlight Optional named numeric vector of positions to mark, for
#'   example `c(MSTN = 6213566)`. Names are used as labels.
#' @param xlim Optional length-2 numeric giving a position range in base pairs.
#' @param units Axis units for position: `"Mb"` (default), `"kb"` or `"bp"`.
#' @param show_raw Draw the unsmoothed per-SNP scores. Default `TRUE`.
#' @param title,subtitle Plot labels.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_threshold(css_smooth(res))
#' css_chrom_plot(res, chr = 2)
#'
#' @export
css_chrom_plot <- function(x, chr, highlight = NULL, xlim = NULL,
                           units = c("Mb", "kb", "bp"),
                           show_raw = TRUE, title = NULL, subtitle = NULL) {
  units <- match.arg(units)
  .require_col(x, c("chr", "pos", "css"), "css_chrom_plot")
  if (missing(chr) || length(chr) != 1L) .stopf("`chr` must be a single chromosome.")

  div <- switch(units, Mb = 1e6, kb = 1e3, bp = 1)
  chr_target <- as.character(chr)
  d <- data.table::as.data.table(x)[as.character(chr) == chr_target]
  if (!nrow(d)) {
    .stopf("Chromosome `%s` not found. Available: %s.",
           as.character(chr), paste(levels(x$chr), collapse = ", "))
  }
  if (!is.null(xlim)) d <- d[pos >= xlim[1] & pos <= xlim[2]]
  d[, x := pos / div]

  thr <- attr(x, "css_threshold")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = x))

  if (!is.null(highlight)) {
    hl <- data.table::data.table(x = unname(highlight) / div,
                                 label = names(highlight) %||%
                                   rep("", length(highlight)))
    p <- p + ggplot2::geom_vline(data = hl, ggplot2::aes(xintercept = x),
                                 colour = "#1B7837", linewidth = 0.6, alpha = 0.9)
  }

  if (show_raw) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = css), colour = .css_chrom_colours[1],
                                 size = 0.7, alpha = 0.55)
  }
  if ("css_smooth" %in% names(d)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = css_smooth),
                                colour = .css_chrom_colours[2], linewidth = 0.7,
                                na.rm = TRUE)
  }
  if (!is.null(thr)) {
    p <- p + ggplot2::geom_hline(yintercept = thr$cut, linetype = "dashed",
                                 colour = .css_sig_colour, linewidth = 0.4)
  }
  if (!is.null(highlight) && !is.null(names(highlight))) {
    hl <- data.table::data.table(x = unname(highlight) / div, label = names(highlight))
    ytop <- max(c(d$css, d$css_smooth, thr$cut), na.rm = TRUE)
    p <- p + ggplot2::geom_text(data = hl, ggplot2::aes(x = x, label = label),
                                y = ytop, inherit.aes = FALSE,
                                vjust = 1.1, hjust = -0.1,
                                size = 3, fontface = "italic", colour = "#1B7837")
  }

  p +
    ggplot2::labs(
      x = sprintf("Chromosome %s position (%s)", as.character(chr), units),
      y = expression(CSS ~ (-log[10] * italic(p))),
      title = title,
      subtitle = subtitle %||% if ("css_smooth" %in% names(d))
        "Points: per-SNP CSS. Line: sliding-window mean." else NULL
    ) +
    css_theme()
}

#' Regional plot with constituent test tracks
#'
#' Zooms into one region and shows each constituent selection test alongside
#' CSS, which is how the source papers diagnose whether a composite signal is
#' driven by one test or supported by several.
#'
#' @param x A `css_result` that still carries its constituent test columns.
#' @param region Either a single row of a `css_regions` object, or a list/vector
#'   with `chr`, `start` and `end`.
#' @param pad Extra base pairs shown on each side. Default `5e5`.
#' @param tests Optional character vector selecting which constituent tests to
#'   show. Defaults to all of them.
#' @param smooth Draw the smoothed track for each test where available.
#'   Default `TRUE`.
#' @param genes Optional `data.frame` with `start`, `end` and `name` columns,
#'   drawn as a gene track.
#' @param units Axis units: `"Mb"` (default), `"kb"` or `"bp"`.
#' @param title Plot title.
#'
#' @return A [ggplot2::ggplot] object, faceted by test with CSS on top.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_threshold(css_smooth(res))
#' reg <- css_regions(res)
#' if (nrow(reg)) css_region_plot(res, reg[1])
#'
#' @export
css_region_plot <- function(x, region, pad = 5e5, tests = NULL,
                            smooth = TRUE, genes = NULL,
                            units = c("Mb", "kb", "bp"), title = NULL) {
  units <- match.arg(units)
  div <- switch(units, Mb = 1e6, kb = 1e3, bp = 1)

  rchr <- as.character(region$chr[1])
  rstart <- as.numeric(region$start[1])
  rend <- as.numeric(region$end[1])
  if (is.na(rchr) || is.na(rstart) || is.na(rend)) {
    .stopf("`region` must supply `chr`, `start` and `end`.")
  }

  all_tests <- names(attr(x, "css_call")$tests)
  if (is.null(tests)) tests <- all_tests
  tests <- intersect(tests, names(x))

  d <- data.table::as.data.table(x)[
    as.character(chr) == rchr & pos >= (rstart - pad) & pos <= (rend + pad)]
  if (!nrow(d)) .stopf("No SNPs in the requested region.")

  panels <- c("CSS", tests)
  long <- data.table::rbindlist(lapply(panels, function(tc) {
    val <- if (tc == "CSS") d$css else d[[tc]]
    sm_col <- if (tc == "CSS") "css_smooth" else paste0(tc, "_smooth")
    data.table::data.table(
      pos = d$pos, test = tc, value = val,
      smooth_value = if (smooth && sm_col %in% names(d)) d[[sm_col]] else NA_real_
    )
  }))
  long[, test := factor(test, levels = panels)]
  long[, x := pos / div]

  # Each facet has its own y scale, so the shaded band needs per-facet bounds
  # rather than a single infinite rectangle.
  bands <- long[, .(ymin = min(value, na.rm = TRUE),
                    ymax = max(value, na.rm = TRUE)), by = test]
  bands[, `:=`(xmin = rstart / div, xmax = rend / div)]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = value)) +
    ggplot2::geom_rect(
      data = bands, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = .css_sig_colour, alpha = 0.06
    ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7, colour = .css_chrom_colours[1]) +
    ggplot2::facet_wrap(~ test, ncol = 1, scales = "free_y", strip.position = "left")

  if (smooth && any(!is.na(long$smooth_value))) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = smooth_value),
                                colour = .css_chrom_colours[2],
                                linewidth = 0.6, na.rm = TRUE)
  }

  if (!is.null(genes) && nrow(genes)) {
    gd <- data.table::as.data.table(genes)
    gd <- gd[start <= (rend + pad) & end >= (rstart - pad)]
    if (nrow(gd)) {
      gd[, `:=`(ybase = min(bands$ymin))]
      p <- p + ggplot2::geom_segment(
        data = gd, inherit.aes = FALSE,
        ggplot2::aes(x = start / div, xend = end / div, y = ybase, yend = ybase),
        colour = "#1B7837", linewidth = 2
      )
    }
  }

  p +
    ggplot2::labs(
      x = sprintf("Chromosome %s position (%s)", rchr, units),
      y = NULL,
      title = title %||% sprintf("Chromosome %s: %.2f-%.2f %s",
                                 rchr, rstart / div, rend / div, units),
      subtitle = "Shaded band: called region. Line: sliding-window mean."
    ) +
    css_theme() +
    ggplot2::theme(strip.placement = "outside")
}
