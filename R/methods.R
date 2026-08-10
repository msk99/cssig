#' @export
print.css_input <- function(x, ...) {
  tests <- attr(x, "css_tests")
  na <- attr(x, "css_na")
  cat(sprintf("<css_input> %s SNPs on %d chromosome%s, %d constituent test%s\n",
              format(nrow(x), big.mark = ","),
              data.table::uniqueN(x$chr),
              if (data.table::uniqueN(x$chr) > 1) "s" else "",
              length(tests), if (length(tests) > 1) "s" else ""))
  for (tc in names(tests)) {
    cat(sprintf("  %-12s direction = %-5s  missing = %d\n",
                tc, tests[[tc]], na[[tc]]))
  }
  cat("\n")
  print(data.table::as.data.table(x), class = FALSE, ...)
  invisible(x)
}

#' @export
print.css_result <- function(x, ...) {
  cl <- attr(x, "css_call")
  thr <- attr(x, "css_threshold")
  sm <- attr(x, "css_smooth_call")

  cat(sprintf("<css_result> %s SNPs, m = %d constituent tests\n",
              format(nrow(x), big.mark = ","), cl$m %||% NA))
  if (!is.null(cl)) {
    cat(sprintf("  ties = \"%s\", na_action = \"%s\"%s\n",
                cl$ties, cl$na_action,
                if (!is.null(cl$weights)) "  [WEIGHTED: not the published method]" else ""))
  }
  if (!is.null(sm)) {
    cat(sprintf("  smoothed: %.0f kb window (half-width %.0f kb), min %d SNPs, on \"%s\"\n",
                2 * sm$half_width / 1000, sm$half_width / 1000, sm$min_snps, sm$on))
  }
  if (!is.null(thr)) {
    cat(sprintf("  threshold: top %s of %s CSS at %.3f -> %d significant SNPs\n",
                .pct(thr$top), thr$on, thr$cut, thr$n_significant))
  }
  if ("qval" %in% names(x)) {
    fd <- attr(x, "css_fdr")
    cat(sprintf("  FDR (%s): %d SNPs with q <= 0.05\n",
                fd$method %||% "?", sum(x$qval <= 0.05, na.rm = TRUE)))
  }
  cat("\n")
  print(data.table::as.data.table(x), class = FALSE, ...)
  invisible(x)
}

#' @export
print.css_regions <- function(x, ...) {
  cl <- attr(x, "css_regions_call")
  cat(sprintf("<css_regions> %d region%s, method = \"%s\"\n",
              nrow(x), if (nrow(x) == 1) "" else "s", cl$method %||% "?"))
  if (nrow(x)) {
    span <- x$end - x$start
    cat(sprintf("  span: median %.2f Mb, range %.2f-%.2f Mb; %d significant SNPs total\n",
                stats::median(span) / 1e6, min(span) / 1e6, max(span) / 1e6,
                sum(x$n_sig)))
  }
  cat("\n")
  print(data.table::as.data.table(x), class = FALSE, ...)
  invisible(x)
}

#' Summarise a CSS result
#'
#' @param object A `css_result`.
#' @param ... Unused.
#' @return Invisibly, a list of summary components.
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' summary(res)
#' @export
summary.css_result <- function(object, ...) {
  cl <- attr(object, "css_call")
  tests <- names(cl$tests)
  present <- intersect(tests, names(object))

  out <- list(
    n_snps = nrow(object),
    n_chr = data.table::uniqueN(object$chr),
    m = cl$m,
    css = summary(object$css),
    per_chr = data.table::as.data.table(object)[
      , .(n = .N, max_css = max(css, na.rm = TRUE)), by = chr]
  )
  if (length(present) > 1) {
    out$cor <- stats::cor(as.matrix(as.data.frame(object)[, c(present, "css")]),
                          use = "pairwise.complete.obs", method = "spearman")
  }

  cat(sprintf("Composite selection signals: %s SNPs on %d chromosomes, m = %d tests\n\n",
              format(out$n_snps, big.mark = ","), out$n_chr, out$m))
  cat("CSS distribution:\n"); print(out$css)
  if (!is.null(out$cor)) {
    cat("\nSpearman correlation (constituents and CSS):\n")
    print(round(out$cor, 3))
  }
  invisible(out)
}

#' Summarise called regions
#'
#' Prints the table in the shape of Table 2 of Randhawa et al. (2014): one row
#' per region with its span, peak and the number of significant SNPs
#' contributed by CSS and by each constituent test.
#'
#' @param object A `css_regions` object.
#' @param ... Unused.
#' @return Invisibly, a `data.table` of the printed summary.
#' @examples
#' data(css_sim_small)
#' res <- css_threshold(css_smooth(css(css_input(css_sim_small,
#'          tests = c(fst = "high", xpehh = "high", ddaf = "high")))))
#' summary(css_regions(res))
#' @export
summary.css_regions <- function(object, ...) {
  if (!nrow(object)) {
    cat("No regions called.\n")
    return(invisible(object))
  }
  d <- data.table::as.data.table(object)
  n_cols <- grep("^n_", names(d), value = TRUE)
  out <- d[, c(list(
    region = region_id, chr = chr,
    position = sprintf("%.2f-%.2f", start / 1e6, end / 1e6),
    span_kb = round((end - start) / 1e3),
    peak_Mb = round(peak_pos / 1e6, 3),
    peak_css = round(peak_css, 2)
  ), lapply(.SD, identity)), .SDcols = n_cols]
  print(out, class = FALSE)
  cat(sprintf("\nPositions in Mb. n_* columns count SNPs in each test's own top fraction.\n"))
  invisible(out)
}

#' @export
as.data.frame.css_result <- function(x, ...) {
  as.data.frame(data.table::as.data.table(x), ...)
}

#' @export
plot.css_result <- function(x, ...) css_manhattan(x, ...)

#' @export
plot.css_regions <- function(x, y, ...) {
  .stopf("Plot regions by passing them to css_manhattan(res, regions = ...) or css_region_plot().")
}
