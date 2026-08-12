#' Smooth CSS scores in sliding genomic windows
#'
#' Averages a score over a window of fixed *physical* width centred on each SNP,
#' as in both source papers: Randhawa et al. (2014) describe 1 Mb windows
#' centred at each SNP, and Randhawa et al. (2015) describe averaging 0.5 Mb on
#' each side, which is the same window. Smoothing suppresses single-SNP noise
#' and exploits the fact that hitchhiking makes true signals cluster.
#'
#' @details
#' Windows never span a chromosome boundary. Windows containing fewer than
#' `min_snps` SNPs are returned as `NA` and are excluded from thresholding and
#' region calling; Randhawa et al. (2014) prune such low-density windows for the
#' same reason.
#'
#' The implementation is vectorised with [findInterval()] and [cumsum()] per
#' chromosome, so cost does not grow with window width.
#' [data.table::frollmean()] is deliberately not used: it windows by *row
#' count*, whereas CSS windows are defined by base-pair distance over
#' irregularly spaced SNPs.
#'
#' A reciprocal result from [css_reciprocal()] carries two directed scores
#' rather than one `css` column, so both are smoothed: `css_pos_smooth`,
#' `css_neg_smooth` and their difference `css_signed_smooth` are added, and
#' [css_manhattan_mirror()] then plots the smoothed scores by default.
#' `on = "zbar"` is not available in that case.
#'
#' @param x A `css_result` from [css()] or [css_reciprocal()], or any keyed
#'   table with `chr`, `pos` and the columns named in `cols`.
#' @param half_width Half the window width in base pairs. Default `5e5`, giving
#'   the 1 Mb window of the papers.
#' @param min_snps Minimum number of SNPs in a window for it to be retained.
#'   Default `5`.
#' @param on What to smooth. `"css"` (default) averages the
#'   \eqn{-\log_{10}(p)} CSS score, matching the published figures. `"zbar"`
#'   averages the mean z instead and re-derives a p-value from the smoothed
#'   value, which is not what the papers plot but is available for users who
#'   prefer it.
#' @param cols Optional character vector of additional columns to smooth, for
#'   example the constituent tests, for side-by-side comparison. Smoothed
#'   columns are suffixed `_smooth`.
#' @param .copy If `TRUE`, work on a copy and leave `x` untouched.
#'
#' @return `x` with `css_smooth` and `n_window` added, plus one `<col>_smooth`
#'   column for each entry of `cols`. For a reciprocal result the smoothed
#'   columns are `css_pos_smooth`, `css_neg_smooth` and `css_signed_smooth`.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_smooth(res)
#' summary(res$css_smooth)
#'
#' @export
css_smooth <- function(x,
                       half_width = 5e5,
                       min_snps = 5L,
                       on = c("css", "zbar"),
                       cols = NULL,
                       .copy = FALSE) {
  on <- match.arg(on)
  if (!is.numeric(half_width) || length(half_width) != 1L || half_width <= 0) {
    .stopf("`half_width` must be a single positive number of base pairs.")
  }
  # A reciprocal result has no single score column; smooth both directions.
  recip <- !is.null(attr(x, "css_reciprocal"))
  if (recip && on == "zbar") {
    .stopf("`on = \"zbar\"` is not available for a reciprocal result; its directed scores are already -log10(p).")
  }
  primary <- if (recip) "css_pos" else on
  .require_col(x, c("chr", "pos", if (recip) c("css_pos", "css_neg") else on),
               "css_smooth")

  x <- .maybe_copy(x, .copy)
  if (!data.table::haskey(x) || !identical(data.table::key(x), c("chr", "pos"))) {
    data.table::setkeyv(x, c("chr", "pos"))
  }

  target <- c(if (recip) c("css_pos", "css_neg") else on, cols)
  missing_cols <- setdiff(target, names(x))
  if (length(missing_cols)) {
    .stopf("Column%s to smooth not found: %s.",
           if (length(missing_cols) > 1) "s" else "",
           paste0("`", missing_cols, "`", collapse = ", "))
  }

  out_name_for <- function(cl) {
    if (!recip && identical(cl, on)) "css_smooth" else paste0(cl, "_smooth")
  }
  for (cl in target) {
    res <- x[, .window_mean(pos, get(cl), half_width), by = chr]
    data.table::set(x, j = out_name_for(cl), value = res$mean)
    if (identical(cl, primary)) data.table::set(x, j = "n_window", value = res$n)
  }

  # Prune sparse windows. Done once, on the window count from the primary
  # column, so every smoothed column is masked consistently.
  if (min_snps > 1L) {
    sparse <- which(x$n_window < min_snps)
    if (length(sparse)) {
      for (cl in target) {
        data.table::set(x, i = sparse, j = out_name_for(cl), value = NA_real_)
      }
      .msgf("Masked %d of %d windows (%s) containing fewer than %d SNPs.",
            length(sparse), nrow(x), .pct(length(sparse) / nrow(x)), min_snps)
    }
  }

  if (recip) {
    x[, css_signed_smooth := css_pos_smooth - css_neg_smooth]
  } else if (on == "zbar") {
    m <- attr(x, "css_call")$m
    if (!is.null(m)) {
      x[, css_smooth := -log10(stats::pnorm(sqrt(m) * css_smooth, lower.tail = FALSE))]
    }
  }

  data.table::setattr(x, "css_smooth_call",
                      list(half_width = half_width, min_snps = min_snps,
                           on = on, cols = cols))
  x[]
}

# Mean of `value` over all SNPs within +/- half_width of each position.
# `pos` is assumed sorted ascending (guaranteed by the (chr, pos) key).
# NA values are excluded from both the sum and the count.
.window_mean <- function(pos, value, half_width) {
  n <- length(pos)
  if (n == 0L) return(list(mean = numeric(0), n = integer(0)))

  ok <- !is.na(value)
  v <- ifelse(ok, value, 0)

  # Half-open interval search: lo = number of SNPs strictly before the window,
  # hi = number of SNPs up to and including the window's right edge.
  lo <- findInterval(pos - half_width, pos, left.open = TRUE)
  hi <- findInterval(pos + half_width, pos, left.open = FALSE)

  # cs[k + 1] is the sum of the first k elements, so the window spanning
  # elements (lo + 1) .. hi sums to cs[hi + 1] - cs[lo + 1]. Indexing with `lo`
  # directly would be an off-by-one, and silently so: R drops zero subscripts
  # rather than returning zero, which shortens the result instead of erroring.
  cs_v <- c(0, cumsum(v))
  cs_n <- c(0L, cumsum(as.integer(ok)))

  tot <- cs_v[hi + 1L] - cs_v[lo + 1L]
  cnt <- cs_n[hi + 1L] - cs_n[lo + 1L]

  list(mean = ifelse(cnt > 0L, tot / cnt, NA_real_), n = cnt)
}
