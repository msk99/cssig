#' Reciprocal (signed) composite selection signals
#'
#' Randhawa et al. (2015) compare two cohorts in both directions: each cohort in
#' turn plays the selected population, with the other as reference. The two runs
#' are plotted mirrored about zero (their Figure 2), so that a signal's sign
#' says which cohort carries the selection signature.
#'
#' @details
#' The two inputs must describe the same SNPs. `forward` should hold the
#' statistics computed with cohort A as the selected population and B as
#' reference; `reverse` the same tests with the roles swapped. For directional
#' tests such as XP-EHH and \eqn{\Delta}DAF this means recomputing them; for
#' symmetric tests such as \eqn{F_{ST}} the same column is reused, which is
#' correct — \eqn{F_{ST}} carries no direction, and it is the directional tests
#' that break the tie.
#'
#' CSS is computed independently in each direction, so the ranks in one
#' direction do not depend on the other. Thresholds are likewise per direction.
#'
#' @param forward,reverse `css_input` objects for the two directions.
#' @param labels Length-2 character vector naming the two cohorts, used by
#'   [css_manhattan_mirror()]. Default `c("forward", "reverse")`.
#' @param ... Passed to [css()].
#'
#' @return A `css_result` keyed on `(chr, pos)` with columns `css_pos` (forward
#'   direction), `css_neg` (reverse) and `css_signed`, the signed score
#'   `css_pos - css_neg` used for mirrored plotting. `plot()` on the result
#'   draws [css_manhattan_mirror()], and [css_smooth()] smooths both directions
#'   (`css_pos_smooth`, `css_neg_smooth`, `css_signed_smooth`).
#'
#' @examples
#' data(css_sim_small)
#' fwd <- css_input(css_sim_small,
#'                  tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' # Reverse direction: negate the directional tests, keep F_ST as is.
#' rev_dat <- data.table::copy(css_sim_small)
#' rev_dat$xpehh <- -rev_dat$xpehh
#' rev_dat$ddaf  <- -rev_dat$ddaf
#' rvs <- css_input(rev_dat,
#'                  tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' res <- css_reciprocal(fwd, rvs, labels = c("selected", "reference"))
#' head(res)
#'
#' @export
css_reciprocal <- function(forward, reverse,
                           labels = c("forward", "reverse"), ...) {
  if (!inherits(forward, "css_input") || !inherits(reverse, "css_input")) {
    .stopf("`forward` and `reverse` must both come from `css_input()`.")
  }
  if (length(labels) != 2L) .stopf("`labels` must have exactly two entries.")

  if (nrow(forward) != nrow(reverse) ||
      !identical(as.character(forward$chr), as.character(reverse$chr)) ||
      !isTRUE(all.equal(forward$pos, reverse$pos))) {
    .stopf(paste0("`forward` and `reverse` must cover the same SNPs in the same order.\n",
                  "Use css_merge_tests() to align them first."))
  }

  f <- css(forward, .copy = TRUE, ...)
  r <- css(reverse, .copy = TRUE, ...)

  out <- f[, .(chr, pos, snp)]
  out[, `:=`(
    css_pos    = f$css,
    css_neg    = r$css,
    css_signed = f$css - r$css,
    p_pos      = f$p,
    p_neg      = r$p
  )]
  data.table::setkeyv(out, c("chr", "pos"))

  data.table::setattr(out, "class", c("css_result", class(out)))
  data.table::setattr(out, "css_call", attr(f, "css_call"))
  data.table::setattr(out, "css_reciprocal", list(labels = labels))
  out[]
}
