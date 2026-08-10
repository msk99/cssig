#' Circular genome plot of CSS and its constituent tests
#'
#' Reproduces the layout of Figure 4 of Randhawa et al. (2014): concentric
#' tracks showing the smoothed CSS and each constituent test around the genome,
#' with called regions highlighted.
#'
#' @details
#' Drawn with \pkg{circlize}, which must be installed. The source papers used
#' RCircos; the arrangement here is equivalent, not identical.
#'
#' Unlike the other plotting functions this one draws to the active graphics
#' device rather than returning a \pkg{ggplot2} object, because \pkg{circlize}
#' is built on base graphics.
#'
#' @param x A `css_result`, normally after [css_smooth()].
#' @param regions Optional `css_regions` object; highlighted on the CSS track.
#' @param tests Constituent tests to draw as inner tracks. Defaults to all.
#' @param smoothed Use the smoothed columns where available. Default `TRUE`.
#' @param track_height Height of each track as a fraction of the radius.
#' @param labels Optional `data.frame` with `chr`, `pos` and `name` for outer
#'   text labels, for example candidate genes.
#' @param title Optional title drawn in the centre.
#'
#' @return Invisibly `NULL`, called for its side effect of drawing a plot.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("circlize", quietly = TRUE)) {
#'   data(css_sim)
#'   res <- css(css_input(css_sim,
#'                        tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#'   res <- css_threshold(css_smooth(res))
#'   css_circos(res, css_regions(res))
#' }
#' }
#'
#' @export
css_circos <- function(x, regions = NULL, tests = NULL, smoothed = TRUE,
                       track_height = 0.13, labels = NULL, title = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    .stopf(paste0("css_circos() needs the circlize package.\n",
                  "Install it with install.packages(\"circlize\")."))
  }
  .require_col(x, c("chr", "pos", "css"), "css_circos")

  all_tests <- names(attr(x, "css_call")$tests)
  if (is.null(tests)) tests <- all_tests
  tests <- intersect(tests, names(x))

  pick <- function(base) {
    sm <- if (base == "CSS") "css_smooth" else paste0(base, "_smooth")
    raw <- if (base == "CSS") "css" else base
    if (smoothed && sm %in% names(x)) x[[sm]] else x[[raw]]
  }

  d <- data.table::data.table(chr = as.character(x$chr), pos = as.numeric(x$pos))
  cyto <- d[, .(start = 0, end = max(pos, na.rm = TRUE)), by = chr]
  cyto <- as.data.frame(cyto[, .(chr, start, end)])

  op <- circlize::circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1.2,
                             start.degree = 90, track.margin = c(0.005, 0.005))
  on.exit({circlize::circos.clear(); suppressWarnings(circlize::circos.par(op))},
          add = TRUE)

  circlize::circos.initializeWithIdeogram(cyto, plotType = c("labels"))

  panels <- c("CSS", tests)
  for (tc in panels) {
    v <- pick(tc)
    dd <- data.frame(chr = d$chr, start = d$pos, end = d$pos, value = v)
    dd <- dd[!is.na(dd$value), , drop = FALSE]
    circlize::circos.genomicTrack(
      dd, numeric.column = "value", track.height = track_height,
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicLines(region, value, type = "l",
                                      col = "grey20", lwd = 0.5, ...)
      })
    circlize::circos.text(0, 0.5, tc, sector.index = circlize::get.all.sector.index()[1],
                          track.index = circlize::get.current.track.index(),
                          facing = "downward", adj = c(1.2, 0.5), cex = 0.6)
  }

  if (!is.null(regions) && nrow(regions)) {
    rd <- data.frame(chr = as.character(regions$chr),
                     start = regions$start, end = regions$end)
    circlize::circos.genomicTrackPlotRegion(
      rd, ylim = c(0, 1), track.height = 0.04,
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                                     col = .css_sig_colour, border = NA, ...)
      })
  }

  if (!is.null(labels) && nrow(labels)) {
    ld <- data.frame(chr = as.character(labels$chr),
                     start = labels$pos, end = labels$pos)
    circlize::circos.genomicLabels(ld, labels = labels$name, side = "outside",
                                   cex = 0.55)
  }

  if (!is.null(title)) text(0, 0, title, cex = 0.9, font = 2)
  invisible(NULL)
}
