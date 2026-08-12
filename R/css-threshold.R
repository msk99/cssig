#' Apply an empirical significance threshold to CSS scores
#'
#' Both source papers declare SNPs significant when their score falls in the top
#' 0.1% of the genome-wide empirical distribution, rather than referring to any
#' theoretical null. This computes that threshold and flags the SNPs above it.
#'
#' @param x A `css_result`, normally after [css_smooth()].
#' @param top Upper tail fraction to declare significant. Default `0.001`
#'   (top 0.1%), as used in both papers.
#' @param top2 A second, more permissive fraction used by the `"flank"` region
#'   caller of [css_regions()]. Default `0.01` (top 1%).
#' @param on Which score to threshold: `"smoothed"` (default) or `"raw"`.
#' @param .copy If `TRUE`, work on a copy and leave `x` untouched.
#'
#' @return `x` with logical columns `significant` and `significant2` added. The
#'   numeric cut-offs are stored in the `css_threshold` attribute.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_threshold(css_smooth(res))
#' attr(res, "css_threshold")
#'
#' @export
css_threshold <- function(x, top = 0.001, top2 = 0.01,
                          on = c("smoothed", "raw"), .copy = FALSE) {
  on <- match.arg(on)
  if (top <= 0 || top >= 1) .stopf("`top` must be strictly between 0 and 1.")
  if (top2 <= 0 || top2 >= 1) .stopf("`top2` must be strictly between 0 and 1.")
  if (top2 < top) .stopf("`top2` (%g) must be at least as permissive as `top` (%g).", top2, top)

  score_col <- if (on == "smoothed") "css_smooth" else "css"
  .require_col(x, score_col, "css_threshold")

  x <- .maybe_copy(x, .copy)

  v <- x[[score_col]]
  n_ok <- sum(!is.na(v))
  if (n_ok == 0L) .stopf("All values of `%s` are missing.", score_col)

  cut1 <- stats::quantile(v, probs = 1 - top,  na.rm = TRUE, names = FALSE, type = 7)
  cut2 <- stats::quantile(v, probs = 1 - top2, na.rm = TRUE, names = FALSE, type = 7)

  data.table::set(x, j = "significant",  value = !is.na(v) & v >= cut1)
  data.table::set(x, j = "significant2", value = !is.na(v) & v >= cut2)

  data.table::setattr(x, "css_threshold", list(
    top = top, top2 = top2, on = on, score_col = score_col,
    cut = cut1, cut2 = cut2, n_scored = n_ok,
    n_significant = sum(x$significant), n_significant2 = sum(x$significant2)
  ))

  .msgf("Threshold on %s CSS: top %s at %.4f (%d SNPs), top %s at %.4f (%d SNPs).",
        on, .pct(top), cut1, sum(x$significant), .pct(top2), cut2, sum(x$significant2))
  x[]
}
