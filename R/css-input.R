#' Assemble constituent selection statistics for a CSS analysis
#'
#' Validates a table of pre-computed constituent selection-test statistics and
#' prepares it for [css()]. One row per SNP; one column per constituent test.
#'
#' @section Reference semantics:
#' `css_input()` takes a defensive [data.table::copy()] of `data`, so the object
#' you pass in is never modified. From that point on the pipeline stages
#' ([css()], [css_smooth()], [css_threshold()], [css_fdr()]) add columns *by
#' reference* to the object `css_input()` returned, which avoids copying large
#' tables between stages. Each stage also returns the object, so both the piped
#' and the in-place idiom work. Pass `.copy = TRUE` to a stage to get a copy
#' instead. See `vignette("cssig")`.
#'
#' @section Direction:
#' CSS ranks every constituent test so that *larger means more evidence of
#' selection*. `tests` records, per test, how to get there:
#' \describe{
#'   \item{`"high"`}{Use as-is. Correct for \eqn{F_{ST}}, XP-EHH computed with
#'     the selected population as the target, and \eqn{\Delta}DAF computed as
#'     selected minus reference.}
#'   \item{`"low"`}{Negate before ranking. For tests where small values
#'     indicate selection, such as Tajima's D or a nucleotide-diversity ratio.}
#'   \item{`"abs"`}{Use the absolute value. For a signed test used without
#'     regard to direction.}
#' }
#'
#' @param data A `data.frame`, `tibble`, `data.table` or matrix with one row
#'   per SNP.
#' @param chr,pos Column names holding the chromosome and the base-pair
#'   position.
#' @param snp Optional column name holding a SNP identifier. If `NULL`, a
#'   column named `snp` is used when present; failing that, an identifier of the
#'   form `"chr:pos"` is generated.
#' @param tests A named character vector mapping column names to directions,
#'   for example `c(fst = "high", xpehh = "high", ddaf = "high")`. Names are the
#'   columns of `data`; values are one of `"high"`, `"low"` or `"abs"`.
#' @param drop_na_pos Drop rows with a missing chromosome or position.
#'   Default `TRUE`.
#'
#' @return A [data.table::data.table] of class `css_input`, keyed on
#'   `(chr, pos)`, with `chr` stored as an ordered factor.
#'
#' @examples
#' data(css_sim_small)
#' x <- css_input(css_sim_small,
#'                tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' x
#'
#' @seealso [css()], [css_merge_tests()]
#' @export
css_input <- function(data,
                      chr = "chr",
                      pos = "pos",
                      snp = NULL,
                      tests,
                      drop_na_pos = TRUE) {
  if (missing(tests) || !length(tests)) {
    .stopf("`tests` must name at least two constituent test columns.")
  }
  if (is.null(names(tests))) {
    .stopf("`tests` must be a *named* vector, e.g. c(fst = \"high\", xpehh = \"high\").")
  }
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) {
    .stopf("`data` must be a data.frame, tibble, data.table or matrix; got <%s>.",
           paste(class(data), collapse = "/"))
  }

  # A column literally called "snp" is what the user means by SNP identifiers,
  # so use it rather than silently discarding it and generating "chr:pos" ids.
  if (is.null(snp) && "snp" %in% names(data)) snp <- "snp"

  test_cols <- names(tests)
  dirs <- unname(as.character(tests))
  bad_dir <- setdiff(dirs, c("high", "low", "abs"))
  if (length(bad_dir)) {
    .stopf("Unknown direction%s %s. Must be one of \"high\", \"low\", \"abs\".",
           if (length(bad_dir) > 1) "s" else "",
           paste0("\"", unique(bad_dir), "\"", collapse = ", "))
  }

  needed <- c(chr, pos, if (!is.null(snp)) snp, test_cols)
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    .stopf("Column%s not found in `data`: %s.\nAvailable columns: %s.",
           if (length(absent) > 1) "s" else "",
           paste0("`", absent, "`", collapse = ", "),
           paste0("`", names(data), "`", collapse = ", "))
  }

  # Refuse names that pipeline stages would later overwrite (PLAN.md 3.7 #5).
  clash <- intersect(test_cols, .css_reserved)
  if (length(clash)) {
    .stopf(paste0("These test columns use names reserved for CSS output: %s.\n",
                  "Rename them before calling `css_input()`."),
           paste0("`", clash, "`", collapse = ", "))
  }

  if (length(test_cols) < 2L) {
    .warnf(paste0("Only one constituent test supplied. CSS with m = 1 is just a ",
                  "rank transform of that test and carries no composite information."))
  }

  # --- defensive copy: the caller's object is never touched ------------------
  dt <- data.table::as.data.table(data)   # copies a data.frame; may alias a data.table
  cols <- c(chr, pos, if (!is.null(snp)) snp, test_cols)
  dt <- dt[, cols, with = FALSE]          # subset always allocates, so we now own `dt`
  new_names <- c("chr", "pos", if (!is.null(snp)) "snp", test_cols)
  data.table::setnames(dt, cols, new_names)
  if (is.null(snp)) dt[, snp := NA_character_]

  # Coerce factors through their labels: as.numeric() on a factor returns the
  # internal level codes, silently replacing every value.
  as_num <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    suppressWarnings(as.numeric(v))
  }
  if (!is.numeric(dt$pos)) {
    data.table::set(dt, j = "pos", value = as_num(dt$pos))
    if (all(is.na(dt$pos))) .stopf("`%s` could not be coerced to numeric positions.", pos)
  }
  for (tc in test_cols) {
    if (!is.numeric(dt[[tc]])) {
      data.table::set(dt, j = tc, value = as_num(dt[[tc]]))
    }
  }

  n0 <- nrow(dt)
  if (drop_na_pos) dt <- dt[!is.na(chr) & !is.na(pos)]
  n_dropped <- n0 - nrow(dt)

  dt[, chr := .as_chrom_factor(chr)]
  if (is.null(snp) || anyNA(dt$snp)) {
    dt[is.na(snp), snp := paste0(as.character(chr), ":", format(pos, scientific = FALSE, trim = TRUE))]
  }

  dup <- duplicated(dt, by = c("chr", "pos"))
  if (any(dup)) {
    .stopf("%d duplicated chr:pos combination%s found (first: %s:%s). Each SNP must appear once.",
           sum(dup), if (sum(dup) > 1) "s" else "",
           as.character(dt$chr[which(dup)[1]]), format(dt$pos[which(dup)[1]], scientific = FALSE))
  }

  data.table::setkeyv(dt, c("chr", "pos"))
  data.table::setcolorder(dt, c("chr", "pos", "snp", test_cols))

  na_counts <- vapply(test_cols, function(tc) sum(is.na(dt[[tc]])), integer(1))

  data.table::setattr(dt, "class", c("css_input", class(dt)))
  data.table::setattr(dt, "css_tests", stats::setNames(dirs, test_cols))
  data.table::setattr(dt, "css_na", na_counts)
  data.table::setattr(dt, "css_n_dropped", n_dropped)

  if (n_dropped > 0L) {
    .msgf("Dropped %d row%s with a missing chromosome or position.",
          n_dropped, if (n_dropped > 1) "s" else "")
  }
  dt
}

#' Join per-test statistic tables into a single CSS input
#'
#' Constituent statistics usually come from different tools, on SNP sets that do
#' not fully agree. This joins them on `(chr, pos)` and reports what each test
#' contributed, so that a silent loss of SNPs at the merge step becomes visible.
#'
#' @param ... Two or more `data.frame`s or `data.table`s, ideally named
#'   (`fst = fst_table, xpehh = xpehh_table`). Each must contain the chromosome
#'   and position columns plus exactly one statistic column, unless `value` is
#'   given.
#' @param chr,pos Column names for chromosome and position, recycled across all
#'   inputs.
#' @param value Optional character vector naming the statistic column in each
#'   input. If `NULL`, each input must have exactly one column beyond `chr`,
#'   `pos` and `snp`.
#' @param all If `TRUE` (default) keep the union of SNPs, leaving `NA` where a
#'   test is missing. If `FALSE` keep only SNPs scored by every test, which is
#'   the behaviour assumed in the source papers.
#' @param quiet Suppress the merge report.
#'
#' @return A [data.table::data.table] suitable for [css_input()].
#'
#' @examples
#' a <- data.frame(chr = 1, pos = c(100, 200, 300), fst   = c(0.1, 0.2, 0.3))
#' b <- data.frame(chr = 1, pos = c(200, 300, 400), xpehh = c(1.1, 2.2, 0.4))
#' css_merge_tests(fst = a, xpehh = b, all = FALSE, quiet = TRUE)
#'
#' @export
css_merge_tests <- function(..., chr = "chr", pos = "pos", value = NULL,
                            all = TRUE, quiet = FALSE) {
  inputs <- list(...)
  if (length(inputs) < 2L) .stopf("Supply at least two tables to merge.")
  nms <- names(inputs)
  if (is.null(nms) || any(!nzchar(nms))) {
    .stopf("All inputs must be named, e.g. css_merge_tests(fst = a, xpehh = b).")
  }
  if (!is.null(value) && length(value) != length(inputs)) {
    .stopf("`value` must have one entry per input table (%d supplied, %d needed).",
           length(value), length(inputs))
  }

  prepped <- vector("list", length(inputs))
  for (i in seq_along(inputs)) {
    d <- data.table::as.data.table(inputs[[i]])
    if (!all(c(chr, pos) %in% names(d))) {
      .stopf("Input `%s` is missing `%s` or `%s`.", nms[i], chr, pos)
    }
    vcol <- if (!is.null(value)) value[i] else setdiff(names(d), c(chr, pos, "snp"))
    if (length(vcol) != 1L) {
      .stopf(paste0("Input `%s` has %d candidate statistic columns (%s). ",
                    "Pass `value` to say which to use."),
             nms[i], length(vcol), paste(vcol, collapse = ", "))
    }
    d <- d[, c(chr, pos, vcol), with = FALSE]
    data.table::setnames(d, c("chr", "pos", nms[i]))
    d[, chr := as.character(chr)]
    prepped[[i]] <- unique(d, by = c("chr", "pos"))
  }

  out <- Reduce(function(a, b) merge(a, b, by = c("chr", "pos"), all = all), prepped)
  data.table::setkeyv(out, c("chr", "pos"))

  if (!quiet) {
    n_in <- vapply(prepped, nrow, integer(1))
    kept <- vapply(nms, function(n) sum(!is.na(out[[n]])), integer(1))
    .msgf("Merged %d tests on (chr, pos), keeping the %s: %d SNPs.",
          length(inputs), if (all) "union" else "intersection", nrow(out))
    for (i in seq_along(nms)) {
      .msgf("  %-12s %7d in -> %7d retained (%s)",
            nms[i], n_in[i], kept[i], .pct(kept[i] / max(1L, n_in[i])))
    }
  }
  out[]
}
