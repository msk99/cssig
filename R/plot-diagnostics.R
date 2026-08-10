#' QQ plot of CSS p-values
#'
#' Under the null the CSS p-values should be uniform, so the points follow the
#' diagonal. Departure at the tail is the signal; departure along the whole line
#' indicates the null is miscalibrated, usually because the constituent tests
#' are strongly correlated.
#'
#' @param x A `css_result`.
#' @param max_points Thin to at most this many points for rendering. The
#'   smallest p-values are always kept. Default `50000`.
#' @param title,subtitle Plot labels.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' css_qq(res)
#'
#' @export
css_qq <- function(x, max_points = 50000L, title = NULL, subtitle = NULL) {
  .require_col(x, "p", "css_qq")
  pv <- sort(x$p[!is.na(x$p)])
  n <- length(pv)
  if (!n) .stopf("No non-missing p-values to plot.")
  expected <- -log10(seq_len(n) / (n + 1))
  observed <- -log10(pv)

  keep <- if (n > max_points) {
    unique(c(seq_len(min(1000L, n)),
             round(seq(1001, n, length.out = max_points - 1000L))))
  } else seq_len(n)

  d <- data.table::data.table(x = expected[keep], y = observed[keep])

  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = .css_sig_colour,
                         linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(size = 0.8, alpha = 0.8, colour = .css_chrom_colours[1]) +
    ggplot2::labs(
      x = expression(Expected ~ -log[10] * italic(p)),
      y = expression(Observed ~ -log[10] * italic(p)),
      title = title,
      subtitle = subtitle %||% sprintf("%s SNPs", format(n, big.mark = ","))
    ) +
    css_theme()
}

#' Histogram of CSS p-values
#'
#' @param x A `css_result`.
#' @param bins Number of histogram bins. Default `50`.
#' @param title,subtitle Plot labels.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' css_pdist(res)
#' @export
css_pdist <- function(x, bins = 50L, title = NULL, subtitle = NULL) {
  .require_col(x, "p", "css_pdist")
  d <- data.table::data.table(p = x$p[!is.na(x$p)])
  ggplot2::ggplot(d, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(bins = bins, fill = .css_chrom_colours[1],
                            colour = "white", linewidth = 0.15) +
    ggplot2::geom_hline(yintercept = nrow(d) / bins, linetype = "dashed",
                        colour = .css_sig_colour, linewidth = 0.4) +
    ggplot2::labs(x = expression(italic(p)), y = "SNPs", title = title,
                  subtitle = subtitle %||%
                    "Dashed line: the uniform expectation under the null") +
    css_theme()
}

#' Density of q-values inside called regions versus the rest of the genome
#'
#' Reproduces Figure 2 of Randhawa et al. (2014). A clear separation between the
#' two densities is the evidence that called regions carry a real excess of
#' selection signal rather than the tail of the genome-wide distribution.
#'
#' @param x A `css_result` that has been through [css_fdr()].
#' @param regions A `css_regions` object.
#' @param fdr FDR level to mark with a vertical line. Default `0.05`.
#' @param title,subtitle Plot labels.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_fdr(res, method = "BH")
#' res <- css_threshold(css_smooth(res))
#' reg <- css_regions(res)
#' if (nrow(reg)) css_fdr_density(res, reg)
#'
#' @export
css_fdr_density <- function(x, regions, fdr = 0.05,
                            title = NULL, subtitle = NULL) {
  .require_col(x, "qval", "css_fdr_density")
  if (!nrow(regions)) .stopf("`regions` is empty; nothing to compare against.")

  d <- data.table::data.table(chr = x$chr, pos = x$pos, qval = x$qval)[!is.na(qval)]
  d[, in_region := .in_regions(chr, pos, regions)]
  d[, cohort := factor(ifelse(in_region, "In called regions", "Rest of genome"),
                       levels = c("In called regions", "Rest of genome"))]

  ggplot2::ggplot(d, ggplot2::aes(x = qval, fill = cohort)) +
    ggplot2::geom_density(alpha = 0.65, colour = NA, adjust = 1) +
    ggplot2::geom_vline(xintercept = fdr, linetype = "dashed",
                        colour = "grey25", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("#E08A1E", "#9E9E9E"), name = NULL) +
    ggplot2::labs(x = expression(italic(q) * "-value"), y = "Density",
                  title = title,
                  subtitle = subtitle %||%
                    sprintf("Dashed line at FDR = %.2f. %d SNPs in %d regions.",
                            fdr, sum(d$in_region), nrow(regions))) +
    css_theme()
}

#' Correlation among constituent tests and CSS
#'
#' Randhawa et al. (2014) report low to moderate correlation between pairs of
#' constituent tests but high correlation between CSS and each of them, which is
#' the evidence that CSS captures information from across the tests rather than
#' tracking any single one.
#'
#' @param x A `css_result` still carrying its constituent test columns.
#' @param method Correlation method passed to [stats::cor()]. Default
#'   `"spearman"`, which is the appropriate choice here because CSS is itself a
#'   rank statistic.
#' @param title,subtitle Plot labels.
#'
#' @return A [ggplot2::ggplot] object. The underlying correlation matrix is
#'   attached as the `"cor"` attribute.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' css_test_cor(res)
#'
#' @export
css_test_cor <- function(x, method = c("spearman", "pearson", "kendall"),
                         title = NULL, subtitle = NULL) {
  method <- match.arg(method)
  tests <- names(attr(x, "css_call")$tests)
  cols <- c(intersect(tests, names(x)), "css")
  if (length(cols) < 2L) .stopf("Need at least two columns to correlate.")

  m <- stats::cor(as.matrix(as.data.frame(x)[, cols, drop = FALSE]),
                  use = "pairwise.complete.obs", method = method)

  d <- data.table::as.data.table(as.table(m))
  data.table::setnames(d, c("test1", "test2", "r"))
  d[, test1 := factor(test1, levels = cols)]
  d[, test2 := factor(test2, levels = rev(cols))]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = test1, y = test2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3.1,
                       colour = "grey15") +
    ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                                  midpoint = 0, limits = c(-1, 1),
                                  name = sprintf("%s r", method)) +
    ggplot2::labs(x = NULL, y = NULL, title = title,
                  subtitle = subtitle %||%
                    "CSS should correlate with every constituent; constituents need not correlate with each other") +
    css_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  attr(p, "cor") <- m
  p
}
