#' Call genomic regions under selection
#'
#' Groups significant SNPs into contiguous regions. Two rules are provided,
#' one from each source paper; they differ enough to change results, so both
#' are available rather than silently blended.
#'
#' @details
#' \describe{
#'   \item{`method = "cluster"`}{Randhawa et al. (2014). A cluster of at least
#'     `min_snps` significant SNPs (top 0.1% of smoothed CSS) spanning a window
#'     of `merge_gap` around the core SNP. Region boundaries are the first and
#'     last significant SNP. Clusters closer than `merge_gap` are merged.}
#'   \item{`method = "flank"`}{Randhawa et al. (2015). At least one SNP in the
#'     top 0.1%, flanked by at least `min_flank` adjoining SNPs in the top 1%.
#'     Boundaries are the first and last top-1% SNP of the cluster. Clusters
#'     closer than `merge_gap` are merged.}
#' }
#' "Adjoining" is implemented by distance: top-`top2` SNPs within `merge_gap`
#' of each other form one cluster even when SNPs below the top-`top2` cut lie
#' between them, which is slightly more permissive than a literal reading of
#' the 2015 rule.
#' Both papers additionally extend regions by 0.5 Mb on each side *for gene
#' mining only*. That padding is reported in `start_padded` / `end_padded` and
#' never applied to `start` / `end`, so region boundaries are not silently
#' inflated.
#'
#' @param x A `css_result` that has been through [css_threshold()].
#' @param method `"cluster"` (default) or `"flank"`; see Details.
#' @param min_snps Minimum significant SNPs per cluster for `"cluster"`.
#'   Default `3`.
#' @param min_flank Minimum adjoining top-`top2` SNPs for `"flank"`.
#'   Default `5`.
#' @param merge_gap Maximum distance in base pairs between clusters that are
#'   merged into one region. Default `1e6`.
#' @param flank Padding in base pairs reported for gene mining.
#'   Default `5e5`.
#'
#' @return A [data.table::data.table] of class `css_regions`, one row per
#'   region, with the peak SNP, span, SNP counts and, when the constituent tests
#'   are present, the number of SNPs in each test's own top `top` fraction that
#'   fall inside the region.
#'
#' @examples
#' data(css_sim_small)
#' res <- css(css_input(css_sim_small,
#'                      tests = c(fst = "high", xpehh = "high", ddaf = "high")))
#' res <- css_threshold(css_smooth(res))
#' css_regions(res)
#'
#' @export
css_regions <- function(x,
                        method = c("cluster", "flank"),
                        min_snps = 3L,
                        min_flank = 5L,
                        merge_gap = 1e6,
                        flank = 5e5) {
  method <- match.arg(method)
  .require_col(x, c("chr", "pos", "significant"), "css_regions")
  thr <- attr(x, "css_threshold")
  score_col <- if (!is.null(thr$score_col)) thr$score_col else "css"

  d <- data.table::data.table(
    chr   = x$chr,
    pos   = x$pos,
    snp   = x$snp,
    score = x[[score_col]],
    sig   = x$significant,
    sig2  = if ("significant2" %in% names(x)) x$significant2 else x$significant
  )
  data.table::setkeyv(d, c("chr", "pos"))

  seed_col <- if (method == "cluster") "sig" else "sig2"
  d[, keep := get(seed_col) & !is.na(score)]

  if (!any(d$keep)) {
    return(.empty_regions(method, min_snps, min_flank, merge_gap, flank))
  }

  # Contiguous runs of flagged SNPs, split wherever the gap to the next flagged
  # SNP exceeds merge_gap. rleid over the flag handles the runs; the cumsum
  # over the gap indicator does the splitting and the merging in one step.
  s <- d[keep == TRUE]
  s[, gap := c(0, diff(pos)), by = chr]
  s[, run := cumsum(gap > merge_gap | gap < 0), by = chr]
  s[, region_id := .GRP, by = .(chr, run)]

  regs <- s[, .(
    start   = min(pos),
    end     = max(pos),
    n_snps  = .N,
    n_sig   = sum(sig, na.rm = TRUE)
  ), by = .(chr, region_id)]

  # Peak taken over the top-`top` SNPs where any exist, otherwise over the run.
  peaks <- s[, {
    idx <- if (any(sig)) which(sig) else seq_len(.N)
    k <- idx[which.max(score[idx])]
    .(peak_pos = pos[k], peak_snp = snp[k], peak_css = score[k])
  }, by = .(chr, region_id)]
  regs <- merge(regs, peaks, by = c("chr", "region_id"))

  regs <- switch(
    method,
    cluster = regs[n_sig >= min_snps],
    flank   = regs[n_sig >= 1L & n_snps >= (min_flank + 1L)]
  )

  if (!nrow(regs)) {
    return(.empty_regions(method, min_snps, min_flank, merge_gap, flank))
  }

  regs[, `:=`(start_padded = pmax(0, start - flank), end_padded = end + flank)]

  # How many SNPs in each constituent test's own top fraction fall in each
  # region? This is the comparison reported in Table 2 of the 2014 paper.
  tests <- attr(x, "css_call")$tests
  if (!is.null(tests)) {
    top <- if (!is.null(thr$top)) thr$top else 0.001
    iv <- regs[, .(chr, start, end, region_id)]
    data.table::setkeyv(iv, c("chr", "start", "end"))
    for (tc in names(tests)) {
      if (!tc %in% names(x)) next
      v <- switch(tests[[tc]], high = x[[tc]], low = -x[[tc]], abs = abs(x[[tc]]))
      cutoff <- stats::quantile(v, 1 - top, na.rm = TRUE, names = FALSE)
      pts <- data.table::data.table(chr = x$chr, start = x$pos, end = x$pos)[
        !is.na(v) & v >= cutoff]
      cnt <- if (nrow(pts)) {
        ov <- data.table::foverlaps(pts, iv, type = "within", nomatch = NULL)
        ov[, .N, by = region_id]
      } else {
        data.table::data.table(region_id = integer(0), N = integer(0))
      }
      regs[, (paste0("n_", tc)) := 0L]
      if (nrow(cnt)) {
        regs[cnt, on = "region_id", (paste0("n_", tc)) := i.N]
      }
    }
  }

  data.table::setorderv(regs, c("chr", "start"))
  regs[, region_id := seq_len(.N)]
  data.table::setcolorder(regs, c("region_id", "chr", "start", "end",
                                  "start_padded", "end_padded",
                                  "n_snps", "n_sig",
                                  "peak_pos", "peak_snp", "peak_css"))

  data.table::setattr(regs, "class", c("css_regions", class(regs)))
  data.table::setattr(regs, "css_regions_call", list(
    method = method, min_snps = min_snps, min_flank = min_flank,
    merge_gap = merge_gap, flank = flank, score_col = score_col,
    threshold = thr
  ))
  regs[]
}

.empty_regions <- function(method, min_snps, min_flank, merge_gap, flank) {
  .warnf("No regions met the `%s` criteria. Consider a more permissive `top` in css_threshold().",
         method)
  out <- data.table::data.table(
    region_id = integer(0), chr = factor(), start = numeric(0), end = numeric(0),
    start_padded = numeric(0), end_padded = numeric(0),
    n_snps = integer(0), n_sig = integer(0),
    peak_pos = numeric(0), peak_snp = character(0), peak_css = numeric(0)
  )
  data.table::setattr(out, "class", c("css_regions", class(out)))
  data.table::setattr(out, "css_regions_call", list(
    method = method, min_snps = min_snps, min_flank = min_flank,
    merge_gap = merge_gap, flank = flank
  ))
  out
}
