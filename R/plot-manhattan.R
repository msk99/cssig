#' Genome-wide Manhattan plot of composite selection signals
#'
#' Reproduces the layout of Figure 3 of Randhawa et al. (2014): CSS against
#' genomic position, chromosomes in alternating colours, with the empirical
#' significance threshold drawn as a dashed line.
#'
#' @param x A `css_result`. If [css_threshold()] has been run the threshold line
#'   is drawn automatically.
#' @param score Which score to plot: `"smoothed"` (default when available) or
#'   `"raw"`.
#' @param overlay_raw Draw the raw per-SNP scores behind the smoothed line, as
#'   in Figure 1 of the 2014 paper. Only meaningful when `score = "smoothed"`.
#' @param regions Optional `css_regions` object; significant regions are
#'   highlighted and, if `label` is set, annotated.
#' @param label Optional column name in `regions` to use as a text label, for
#'   example a gene name column the user has added.
#' @param thin Drop a random subset of points below `thin_below` to keep
#'   rendering fast on very dense data. Off by default so that nothing is hidden
#'   without the user asking.
#' @param thin_below CSS value under which thinning applies. Default `1`.
#' @param thin_frac Fraction of low points to keep when thinning. Default `0.1`.
#' @param gap Gap between chromosomes in base pairs, passed to
#'   [css_genome_coords()].
#' @param point_size,point_alpha Point aesthetics.
#' @param title,subtitle Plot labels.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_threshold(css_smooth(res))
#' css_manhattan(res)
#'
#' @export
css_manhattan <- function(x,
                          score = c("smoothed", "raw"),
                          overlay_raw = FALSE,
                          regions = NULL,
                          label = NULL,
                          thin = FALSE,
                          thin_below = 1,
                          thin_frac = 0.1,
                          gap = 2e7,
                          point_size = 0.7,
                          point_alpha = 0.85,
                          title = NULL,
                          subtitle = NULL) {
  score <- if (missing(score) && !"css_smooth" %in% names(x)) "raw" else match.arg(score)
  score_col <- if (score == "smoothed") "css_smooth" else "css"
  if (!"css" %in% names(x) && !is.null(attr(x, "css_reciprocal"))) {
    .stopf(paste0("`x` is a reciprocal result with two directed scores, not one `css`.\n",
                  "Plot it with css_manhattan_mirror()."))
  }
  .require_col(x, c("chr", "pos", score_col), "css_manhattan")

  d <- data.table::data.table(
    chr = x$chr, pos = x$pos, snp = x$snp,
    value = x[[score_col]],
    css_raw = if ("css" %in% names(x)) x$css else NA_real_
  )
  d <- d[!is.na(value)]
  if (!nrow(d)) .stopf("No non-missing values in `%s` to plot.", score_col)

  g <- css_genome_coords(d$chr, d$pos, gap = gap)
  d[, pos_cum := g$pos_cum]
  d[, shade := factor(as.integer(chr) %% 2L)]

  if (thin) d <- .thin_points(d, thin_below, thin_frac)

  thr <- attr(x, "css_threshold")
  cut <- if (!is.null(thr) && thr$score_col == score_col) thr$cut else NULL

  # Full-height bands and edge-anchored labels are conventionally drawn with
  # ymin = -Inf / y = Inf. Infinite extents are avoided here in favour of
  # finite bounds taken from the data: they render identically, and they keep
  # grid from having to resolve infinite coordinates on interactive devices.
  yr <- range(d$value, na.rm = TRUE)
  if (overlay_raw) yr <- range(c(yr, d$css_raw), na.rm = TRUE)
  if (!is.null(cut)) yr <- range(c(yr, cut))
  pad <- if (diff(yr) > 0) diff(yr) * 0.04 else 0.5
  ylo <- yr[1] - pad
  yhi <- yr[2] + pad

  p <- ggplot2::ggplot(d, ggplot2::aes(x = pos_cum, y = value))

  if (!is.null(regions) && nrow(regions)) {
    rb <- .region_bands(regions, g)
    p <- p + ggplot2::geom_rect(
      data = rb, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax), ymin = ylo, ymax = yhi,
      fill = .css_sig_colour, alpha = 0.12
    )
  }

  if (overlay_raw && score == "smoothed" && !all(is.na(d$css_raw))) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(y = css_raw, colour = shade),
      size = point_size, alpha = 0.25, show.legend = FALSE
    )
    p <- p + ggplot2::geom_line(ggplot2::aes(group = chr), colour = "grey15", linewidth = 0.4)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = shade),
                                 size = point_size, alpha = point_alpha,
                                 show.legend = FALSE)
  }

  if (!is.null(cut)) {
    p <- p + ggplot2::geom_hline(yintercept = cut, linetype = "dashed",
                                 colour = .css_sig_colour, linewidth = 0.4)
  }

  if (!is.null(regions) && !is.null(label) && nrow(regions) && label %in% names(regions)) {
    rb <- .region_bands(regions, g)
    rb[, label := as.character(regions[[label]])]
    p <- p + ggplot2::geom_text(
      data = rb, inherit.aes = FALSE,
      ggplot2::aes(x = (xmin + xmax) / 2, label = label), y = yhi,
      vjust = 1.1, size = 2.8, colour = "grey20"
    )
  }

  if (is.null(subtitle) && !is.null(thr)) {
    subtitle <- sprintf("Dashed line: top %s of the genome-wide %s distribution (CSS = %.2f)",
                        .pct(thr$top), score, cut)
  }

  p +
    ggplot2::scale_colour_manual(values = .css_chrom_colours) +
    ggplot2::scale_x_continuous(breaks = g$axis$mid, labels = as.character(g$axis$chr),
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.08))) +
    ggplot2::labs(x = "Chromosome",
                  y = expression(CSS ~ (-log[10] * italic(p))),
                  title = title, subtitle = subtitle) +
    css_theme()
}

#' Mirrored Manhattan plot for reciprocal cohort contrasts
#'
#' Reproduces Figure 2 of Randhawa et al. (2015): the two directions of a
#' reciprocal comparison drawn above and below zero, with a separate empirical
#' threshold for each.
#'
#' @param x A `css_result` from [css_reciprocal()].
#' @param score Which scores to plot: `"smoothed"` (the default once
#'   [css_smooth()] has been run on the reciprocal result) or `"raw"`.
#' @param top Upper tail fraction used for the two threshold lines, applied to
#'   each direction's own distribution of the plotted score. Default `0.001`.
#' @param gap,point_size,point_alpha,title,subtitle As for [css_manhattan()].
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' data(css_sim_small)
#' fwd <- css_input(css_sim_small,
#'                  tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' rd <- data.table::copy(css_sim_small)
#' rd$xpehh <- -rd$xpehh; rd$ddaf <- -rd$ddaf
#' rvs <- css_input(rd, tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' recip <- css_reciprocal(fwd, rvs, labels = c("large", "small"))
#' css_manhattan_mirror(recip)              # raw scores
#' css_manhattan_mirror(css_smooth(recip))  # smoothed, as the 2015 figure plots
#'
#' @export
css_manhattan_mirror <- function(x, score = c("smoothed", "raw"),
                                 top = 0.001, gap = 2e7,
                                 point_size = 0.7, point_alpha = 0.85,
                                 title = NULL, subtitle = NULL) {
  score <- if (missing(score) && !"css_pos_smooth" %in% names(x)) "raw" else match.arg(score)
  up_col   <- if (score == "smoothed") "css_pos_smooth" else "css_pos"
  down_col <- if (score == "smoothed") "css_neg_smooth" else "css_neg"
  .require_col(x, c("chr", "pos", up_col, down_col), "css_manhattan_mirror")
  labels <- attr(x, "css_reciprocal")$labels
  if (is.null(labels)) labels <- c("forward", "reverse")

  d <- data.table::data.table(chr = x$chr, pos = x$pos,
                              up = x[[up_col]], down = x[[down_col]])
  d <- d[!is.na(up) & !is.na(down)]
  g <- css_genome_coords(d$chr, d$pos, gap = gap)
  d[, pos_cum := g$pos_cum]
  d[, shade := factor(as.integer(chr) %% 2L)]

  cut_up   <- stats::quantile(d$up,   1 - top, na.rm = TRUE, names = FALSE)
  cut_down <- stats::quantile(d$down, 1 - top, na.rm = TRUE, names = FALSE)

  long <- data.table::rbindlist(list(
    d[, .(pos_cum, chr, shade, value = up)],
    d[, .(pos_cum, chr, shade, value = -down)]
  ))

  ggplot2::ggplot(long, ggplot2::aes(x = pos_cum, y = value, colour = shade)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = c(cut_up, -cut_down), linetype = "dashed",
                        colour = .css_sig_colour, linewidth = 0.4) +
    ggplot2::annotate("text", x = min(long$pos_cum), y = max(long$value),
                      hjust = -0.1, vjust = 1.2,
                      label = labels[1], size = 3.2, colour = "grey20") +
    ggplot2::annotate("text", x = min(long$pos_cum), y = min(long$value),
                      hjust = -0.1, vjust = -0.5,
                      label = labels[2], size = 3.2, colour = "grey20") +
    ggplot2::scale_colour_manual(values = .css_chrom_colours) +
    ggplot2::scale_x_continuous(breaks = g$axis$mid, labels = as.character(g$axis$chr),
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(labels = function(v) format(abs(v))) +
    ggplot2::labs(x = "Chromosome",
                  y = expression(CSS ~ (-log[10] * italic(p))),
                  title = title,
                  subtitle = subtitle %||%
                    sprintf("Each cohort as the selected population in turn; dashed lines at the top %s of each direction's %s scores",
                            .pct(top), score)) +
    css_theme()
}

.thin_points <- function(d, below, frac) {
  low <- which(d$value < below)
  if (!length(low)) return(d)
  drop <- low[stats::runif(length(low)) > frac]
  if (!length(drop)) return(d)
  .msgf("Thinned %d of %d points below CSS = %g for rendering.",
        length(drop), nrow(d), below)
  d[-drop]
}

.region_bands <- function(regions, g) {
  b <- data.table::as.data.table(regions)[, .(chr, start, end)]
  off <- g$bounds
  b <- merge(b, off[, .(chr, offset)], by = "chr", sort = FALSE)
  b[, `:=`(xmin = start + offset, xmax = end + offset)]
  b[]
}
