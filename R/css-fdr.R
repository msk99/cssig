#' Estimate false discovery rates for CSS p-values
#'
#' @details
#' FDR is computed from the **unsmoothed** CSS p-values. Smoothing averages
#' \eqn{-\log_{10}(p)} over neighbouring SNPs, and the result is no longer a
#' p-value, so feeding smoothed scores into an FDR procedure produces numbers
#' that look like q-values but are not. `css_fdr()` therefore refuses smoothed
#' input unless `force = TRUE`.
#'
#' Three methods are available:
#' \describe{
#'   \item{`"BH"`}{Benjamini-Hochberg via [stats::p.adjust()]. No extra
#'     dependency; the default.}
#'   \item{`"fdrtool"`}{Tail-area q-values from the \pkg{fdrtool} package with
#'     `statistic = "pvalue"`. This is the route used by Randhawa et al.
#'     (2014).}
#'   \item{`"isotonic"`}{Recalibrates the empirical p-value distribution with an
#'     isotonic regression of observed on expected quantiles before passing the
#'     result to \pkg{fdrtool}. Randhawa et al. (2014) use ConReg-R for this
#'     step. **This is a reimplementation in the same spirit, not ConReg-R**,
#'     and will not reproduce it exactly. The recalibration map is returned in
#'     the `css_fdr` attribute so it can be inspected.}
#' }
#'
#' @section The CSS p-value is not calibrated on real data:
#' \eqn{p = 1 - \Phi(m^{1/2}\bar{Z})} assumes the constituent tests are
#' independent, so that \eqn{\bar{Z}} has variance \eqn{1/m}. Real constituent
#' tests are correlated, and the correlation need not be positive: on the
#' shipped [css_sim] data XP-EHH and \eqn{\Delta}DAF have Spearman correlation
#' about -0.32, which *deflates* the variance of \eqn{\bar{Z}}. The observed
#' standard deviation of \eqn{m^{1/2}\bar{Z}} there is 0.90 rather than 1;
#' permuting each test independently restores it to 1.00.
#'
#' The practical consequences:
#' \itemize{
#'   \item CSS p-values, and therefore CSS scores, are useful as a *ranking*
#'     but should not be read as calibrated tail probabilities. This is why both
#'     source papers threshold on the empirical top 0.1% of the genome-wide
#'     distribution rather than on a p-value, and why [css_threshold()] does the
#'     same.
#'   \item FDR estimated from these p-values inherits the miscalibration. On
#'     [css_sim], `fdrtool` estimates the null proportion as 1 and returns
#'     q = 1 everywhere, while `"BH"` returns a usable ordering. Treat q-values
#'     as a descriptive summary of the p-value distribution, not as a
#'     guaranteed error rate, and check [css_qq()] before relying on them.
#' }
#'
#' @param x A `css_result` from [css()].
#' @param method One of `"BH"`, `"fdrtool"`, `"isotonic"`.
#' @param force Allow FDR estimation even if the object has been smoothed.
#'   Default `FALSE`.
#' @param .copy If `TRUE`, work on a copy and leave `x` untouched.
#'
#' @return `x` with a `qval` column added (and `p_adj` for `method = "BH"`).
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_fdr(res, method = "BH")
#' sum(res$qval <= 0.05)
#'
#' @export
css_fdr <- function(x, method = c("BH", "fdrtool", "isotonic"),
                    force = FALSE, .copy = FALSE) {
  method <- match.arg(method)
  .require_col(x, "p", "css_fdr")

  if ("css_smooth" %in% names(x) && !force) {
    .stopf(paste0(
      "`x` has been smoothed. Smoothed CSS scores are averages of -log10(p) and\n",
      "are not p-values, so FDR estimated from them is not interpretable.\n",
      "Run css_fdr() before css_smooth(), or pass force = TRUE if you are sure."))
  }

  x <- .maybe_copy(x, .copy)
  pv <- x$p
  ok <- !is.na(pv)
  q <- rep(NA_real_, length(pv))
  info <- list(method = method, n = sum(ok))

  if (method == "BH") {
    q[ok] <- stats::p.adjust(pv[ok], method = "BH")
    data.table::set(x, j = "p_adj", value = q)
  } else {
    if (!requireNamespace("fdrtool", quietly = TRUE)) {
      .stopf(paste0("method = \"%s\" needs the fdrtool package.\n",
                    "Install it with install.packages(\"fdrtool\"), or use method = \"BH\"."),
             method)
    }
    p_use <- pv[ok]
    if (method == "isotonic") {
      cal <- .recalibrate_p(p_use)
      p_use <- cal$calibrated
      info$calibration <- cal$map
    }
    ft <- fdrtool::fdrtool(p_use, statistic = "pvalue", plot = FALSE, verbose = FALSE)
    q[ok] <- ft$qval
    info$eta0 <- ft$param[1, "eta0"]
  }

  data.table::set(x, j = "qval", value = q)
  data.table::setattr(x, "css_fdr", info)
  .msgf("FDR (%s): %d SNPs with q <= 0.05, %d with q <= 0.01.",
        method, sum(q <= 0.05, na.rm = TRUE), sum(q <= 0.01, na.rm = TRUE))
  x[]
}

# Monotone recalibration of an empirical p-value distribution.
#
# Under a true null the ordered p-values should lie on the uniform quantiles
# k/(n+1). Observed CSS p-values depart from this because the constituent tests
# are correlated along the genome. We fit an isotonic regression of the observed
# ordered p-values on those expected quantiles and use it as a monotone
# calibration map. This is the same intent as the ConReg-R step in Randhawa
# et al. (2014) but is not that algorithm.
.recalibrate_p <- function(p) {
  n <- length(p)
  o <- order(p)
  expected <- seq_len(n) / (n + 1)
  fit <- stats::isoreg(expected, p[o])
  calibrated <- numeric(n)
  # Invert: map each observed p onto the expected quantile scale, holding the
  # fitted curve monotone and bounded in (0, 1).
  yf <- pmin(pmax(fit$yf, .Machine$double.eps), 1 - .Machine$double.eps)
  calibrated[o] <- stats::approx(yf, expected, xout = p[o],
                                 rule = 2, ties = "ordered")$y
  calibrated <- pmin(pmax(calibrated, .Machine$double.xmin), 1)
  list(calibrated = calibrated,
       map = data.table::data.table(observed = p[o], expected = expected, fitted = yf))
}
