#' Internal helpers
#' @keywords internal
#' @noRd
NULL

# Reserved output column names. css_input() refuses input that would collide
# with these, so that pipeline stages can never silently overwrite user data.
.css_reserved <- c(
  "zbar", "p", "css", "n_tests", "css_smooth", "n_window", "significant",
  "significant2", "qval", "p_adj", "pos_cum", "css_pos", "css_neg",
  "css_signed", "css_pos_smooth", "css_neg_smooth", "css_signed_smooth"
)

.stopf <- function(...) stop(sprintf(...), call. = FALSE)
.warnf <- function(...) warning(sprintf(...), call. = FALSE)
.msgf  <- function(...) message(sprintf(...))

# Order chromosome labels the way a genomicist expects: numeric chromosomes
# first in numeric order, then non-numeric ones (X, Y, MT, scaffolds)
# alphabetically. Plain alphabetical sorting puts chromosome 10 before 2,
# which silently scrambles every genome-wide plot.
.order_chrom_levels <- function(x) {
  u <- unique(as.character(x))
  u <- u[!is.na(u)]
  num <- suppressWarnings(as.numeric(sub("^(chr|Chr|CHR|BTA|OAR)", "", u)))
  is_num <- !is.na(num)
  c(u[is_num][order(num[is_num])], sort(u[!is_num]))
}

.as_chrom_factor <- function(x) {
  factor(as.character(x), levels = .order_chrom_levels(x))
}

# Weighted-or-plain check that an object is a css_result carrying `col`.
.require_col <- function(x, col, fun) {
  missing <- setdiff(col, names(x))
  if (length(missing)) {
    .stopf(
      "`%s()` needs the column%s %s, which %s not present.\nRun %s first.",
      fun,
      if (length(missing) > 1) "s" else "",
      paste0("`", missing, "`", collapse = ", "),
      if (length(missing) > 1) "are" else "is",
      switch(missing[1],
        css        = "`css()`",
        css_smooth = "`css_smooth()`",
        significant = "`css_threshold()`",
        css_pos    = "`css_reciprocal()`",
        css_pos_smooth = "`css_smooth()`",
        "the preceding pipeline stage"
      )
    )
  }
  invisible(TRUE)
}

# data.table's `:=` modifies in place. Every exported stage funnels through
# this so the `.copy` contract is implemented in exactly one place.
.maybe_copy <- function(x, .copy) {
  if (isTRUE(.copy)) data.table::copy(x) else x
}

.pct <- function(x) sprintf("%.1f%%", 100 * x)

# Defined locally rather than relying on base R, which only gained `%||%` in
# 4.4.0 while this package supports 4.1.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Logical vector: is each (chr, pos) inside any region? Uses an interval join
# rather than a per-region scan, so cost is n log n in the number of SNPs.
.in_regions <- function(chr, pos, regions) {
  n <- length(pos)
  if (!nrow(regions)) return(rep(FALSE, n))
  pts <- data.table::data.table(
    chr = as.character(chr), start = as.numeric(pos), end = as.numeric(pos),
    .idx = seq_len(n)
  )
  iv <- data.table::as.data.table(regions)[, .(chr = as.character(chr),
                                               start = as.numeric(start),
                                               end = as.numeric(end))]
  data.table::setkeyv(iv, c("chr", "start", "end"))
  ov <- data.table::foverlaps(pts, iv, type = "within", nomatch = NULL)
  out <- rep(FALSE, n)
  if (nrow(ov)) out[ov$.idx] <- TRUE
  out
}
