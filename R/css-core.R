#' Compute composite selection signals
#'
#' The core CSS statistic of Randhawa et al. (2014, 2015). For each constituent
#' test the SNPs are ranked genome-wide, the ranks are rescaled to fractional
#' ranks, those are mapped to standard normal quantiles, and the resulting
#' z-scores are averaged across tests at each SNP. The mean z is converted to an
#' upper-tail p-value and reported as \eqn{CSS = -\log_{10}(p)}.
#'
#' @details
#' With \eqn{m} tests and \eqn{n} SNPs, for test \eqn{i} at SNP \eqn{j}:
#' \deqn{R_{ij} = \mathrm{rank}(T_{ij}), \quad R'_{ij} = R_{ij}/(n+1), \quad
#'       Z_{ij} = \Phi^{-1}(R'_{ij})}
#' \deqn{\bar{Z}_j = m^{-1}\sum_i Z_{ij}, \quad
#'       p_j = 1 - \Phi\!\left(m^{1/2}\bar{Z}_j\right), \quad
#'       CSS_j = -\log_{10}(p_j)}
#' since \eqn{\bar{Z}_j \sim N(0, m^{-1})} under the null.
#'
#' @section CSS is rank-based:
#' Only the *order* of each constituent test matters. Standardising XP-EHH or
#' \eqn{\Delta}DAF to \eqn{N(0,1)}, as both source papers do, leaves CSS exactly
#' unchanged; it matters only for plotting the constituent tests on a common
#' axis. Equally, any strictly increasing transformation of a constituent test
#' gives identical CSS.
#'
#' @section Deviations from the published method:
#' Two situations the papers do not specify are handled explicitly, and both
#' defaults are flagged by the `print()` method when they are in use:
#' \describe{
#'   \item{Ties}{`ties = "average"` by default. Statistics computed from allele
#'     counts are discrete, so ties are normal and harmless. What is not
#'     harmless is a single value shared by a large block of SNPs, such as
#'     \eqn{F_{ST} = 0} wherever a cohort is monomorphic: that block collapses
#'     to one averaged rank. A warning fires when any one value covers more
#'     than 5% of a test.}
#'   \item{Missing values}{`na_action = "pairwise"` averages over the tests
#'     available at each SNP and uses that SNP's own \eqn{m} in the p-value.
#'     This is an extension: a SNP scored by fewer tests is not strictly
#'     comparable to one scored by all of them. Use `"omit"` for the strict
#'     complete-case behaviour of the papers.}
#' }
#' Non-`NULL` `weights` give a Stouffer-style weighted mean z, which is *not*
#' the published method; the default `NULL` reproduces the papers exactly.
#'
#' @section Calibration:
#' `p = 1 - \Phi(\sqrt{m}\bar{Z})` assumes the constituent tests are
#' independent. They are correlated in practice, which moves the true null
#' variance of \eqn{\sqrt{m}\bar{Z}} away from 1 (on [css_sim] the observed sd
#' is 0.90, explained to three decimals by the cross-test correlations).
#' `calibrate = TRUE` estimates the correlation matrix of the per-test z
#' columns genome-wide and divides the statistic by the implied null sd,
#' \eqn{\sigma^2 = (w' R w) / \sum w^2}, the Stouffer analogue of Brown's
#' correction for combining correlated tests (Brown 1975; Kost & McDermott
#' 2002). Under `na_action = "pairwise"` the sd is computed per missingness
#' pattern. For complete data the correction is a monotone rescaling, so the
#' CSS *ranking* and the empirical thresholds of [css_threshold()] are
#' unchanged; only `p`, `css` and downstream q-values move. It is off by
#' default because it is not the published method, and `print()` flags a
#' calibrated run.
#'
#' @param x A `css_input` object from [css_input()].
#' @param ties Tie-handling for the rank step, passed to
#'   [data.table::frank()]: one of `"average"`, `"first"`, `"random"`,
#'   `"dense"`.
#' @param na_action `"pairwise"` (default) or `"omit"`; see Details.
#' @param weights Optional numeric vector of per-test weights, one per test.
#'   `NULL` (default) gives the equally weighted mean of the papers.
#' @param calibrate If `TRUE`, correct the p-value for the estimated
#'   cross-test correlation; see the Calibration section. Default `FALSE`,
#'   the published method.
#' @param .copy If `TRUE`, work on a copy and leave `x` untouched. Default
#'   `FALSE`, which adds the result columns to `x` by reference.
#'
#' @return `x`, with columns `zbar`, `p` and `css` added, classed
#'   `css_result`.
#'
#' @references
#' Randhawa IAS, Khatkar MS, Thomson PC, Raadsma HW (2014).
#' Composite selection signals can localize the trait specific genomic regions
#' in multi-breed populations of cattle and sheep. *BMC Genetics* 15:34.
#' \doi{10.1186/1471-2156-15-34}
#'
#' Randhawa IAS, Khatkar MS, Thomson PC, Raadsma HW (2015).
#' Composite selection signals for complex traits exemplified through bovine
#' stature using multibreed cohorts of European and African *Bos taurus*.
#' *G3* 5:1391-1401. \doi{10.1534/g3.115.017772}
#'
#' For `calibrate = TRUE`:
#' Brown MB (1975). A method for combining non-independent, one-sided tests of
#' significance. *Biometrics* 31:987-992.
#' Kost JT, McDermott MP (2002). Combining dependent p-values.
#' *Statistics & Probability Letters* 60:183-190.
#'
#' @examples
#' data(css_sim_small)
#' x <- css_input(css_sim_small,
#'                tests = c(fst = "high", xpehh = "high", ddaf = "high"))
#' res <- css(x)
#' head(res)
#'
#' @seealso [css_smooth()], [css_threshold()], [css_regions()]
#' @export
css <- function(x,
                ties = c("average", "first", "random", "dense"),
                na_action = c("pairwise", "omit"),
                weights = NULL,
                calibrate = FALSE,
                .copy = FALSE) {
  if (!inherits(x, "css_input")) {
    .stopf("`x` must come from `css_input()`; got <%s>.",
           paste(class(x), collapse = "/"))
  }
  ties <- match.arg(ties)
  na_action <- match.arg(na_action)

  tests <- attr(x, "css_tests")
  test_cols <- names(tests)
  m_all <- length(test_cols)

  if (!is.null(weights)) {
    if (length(weights) != m_all) {
      .stopf("`weights` must have one entry per test (%d needed, %d supplied).",
             m_all, length(weights))
    }
    if (anyNA(weights) || any(weights <= 0)) {
      .stopf("`weights` must be finite and strictly positive.")
    }
  }

  x <- .maybe_copy(x, .copy)

  if (na_action == "omit") {
    keep <- stats::complete.cases(as.data.frame(x)[, test_cols, drop = FALSE])
    if (!all(keep)) {
      .msgf("na_action = \"omit\": dropping %d of %d SNPs with a missing constituent test.",
            sum(!keep), length(keep))
      # Subsetting allocates a new object, so the by-reference contract for the
      # remaining rows still holds on the returned object.
      x <- x[keep]
      data.table::setkeyv(x, c("chr", "pos"))
    }
  }

  n <- nrow(x)
  if (n < 2L) .stopf("Need at least two SNPs to rank; got %d.", n)

  # --- steps 1-3: rank -> fractional rank -> normal quantile ----------------
  zmat <- matrix(NA_real_, nrow = n, ncol = m_all,
                 dimnames = list(NULL, test_cols))
  for (i in seq_along(test_cols)) {
    tc <- test_cols[i]
    v <- x[[tc]]
    v <- switch(tests[[tc]],
                high = v,
                low  = -v,
                abs  = abs(v))

    ok <- !is.na(v)
    n_i <- sum(ok)
    if (n_i < 2L) .stopf("Test `%s` has fewer than two non-missing values.", tc)

    # What distorts ranks is not ties in general -- statistics built from
    # allele counts are discrete, so with tens of thousands of SNPs almost
    # every value is shared with some other SNP -- but a single value held by a
    # large block of SNPs, such as F_ST = 0 wherever a cohort is monomorphic.
    # That block collapses to one averaged rank and flattens real signal. So
    # the check is on the largest tie group, not the overall tied fraction.
    biggest_tie <- max(tabulate(match(v[ok], unique(v[ok]))))
    if (biggest_tie > 1L && biggest_tie / n_i > 0.05) {
      .warnf(paste0("Test `%s`: %s of non-missing values share a single value. ",
                    "That block collapses to one rank, so CSS there is driven by ",
                    "the other tests and by `ties = \"%s\"`."),
             tc, .pct(biggest_tie / n_i), ties)
    }

    r <- rep(NA_real_, n)
    r[ok] <- data.table::frank(v[ok], ties.method = ties)
    zmat[, i] <- stats::qnorm(r / (n_i + 1))
  }

  # --- steps 4-6: mean z -> p -> CSS ----------------------------------------
  if (is.null(weights)) {
    m_j <- rowSums(!is.na(zmat))
    zbar <- rowMeans(zmat, na.rm = TRUE)
    zbar[m_j == 0L] <- NA_real_
    stat <- sqrt(m_j) * zbar
  } else {
    w <- matrix(weights, nrow = n, ncol = m_all, byrow = TRUE)
    w[is.na(zmat)] <- NA_real_
    zsum <- rowSums(w * zmat, na.rm = TRUE)
    wsum <- rowSums(w, na.rm = TRUE)
    w2   <- rowSums(w^2, na.rm = TRUE)
    m_j  <- rowSums(!is.na(zmat))
    zbar <- ifelse(wsum > 0, zsum / wsum, NA_real_)
    # Var(sum w Z / sum w) = sum(w^2) / (sum w)^2  under independence
    stat <- ifelse(wsum > 0, zsum / sqrt(w2), NA_real_)
    stat[m_j == 0L] <- NA_real_
  }

  # --- optional calibration for cross-test correlation ----------------------
  cal <- NULL
  if (isTRUE(calibrate) && m_all >= 2L) {
    # The z columns are marginally N(0,1) by construction, so their Pearson
    # correlation is a rank-based (van der Waerden) correlation, robust to the
    # tests' marginal distributions.
    R <- stats::cor(zmat, use = "pairwise.complete.obs")
    if (anyNA(R)) {
      .warnf("Some test pairs share no SNPs; their correlation is treated as zero for calibration.")
      R[is.na(R)] <- 0
      diag(R) <- 1
    }
    w_cal <- if (is.null(weights)) rep(1, m_all) else weights
    # Var(sum_S w Z / sqrt(sum_S w^2)) = (w_S' R_S w_S) / sum_S w^2, computed
    # once per missingness pattern S rather than per SNP.
    pat <- !is.na(zmat)
    pat_id <- as.integer(pat %*% 2^(seq_len(m_all) - 1L))
    first <- which(!duplicated(pat_id))
    sigma_u <- vapply(first, function(r) {
      S <- which(pat[r, ])
      if (length(S) < 2L) return(1)
      wS <- w_cal[S]
      sqrt(max(drop(t(wS) %*% R[S, S, drop = FALSE] %*% wS) / sum(wS^2), 1e-8))
    }, numeric(1))
    stat <- stat / sigma_u[match(pat_id, pat_id[first])]
    cal <- list(
      R = R,
      sigma = stats::setNames(sigma_u, vapply(first, function(r) {
        paste(test_cols[pat[r, ]], collapse = "+")
      }, character(1)))
    )
  } else if (isTRUE(calibrate)) {
    .warnf("`calibrate = TRUE` needs at least two tests; ignored.")
    calibrate <- FALSE
  }

  pv <- stats::pnorm(stat, lower.tail = FALSE)
  # Guard the extreme tail: pnorm() underflows to 0 beyond about z = 38, which
  # would make CSS infinite. Floor at the smallest representable double.
  n_floor <- sum(pv == 0, na.rm = TRUE)
  if (n_floor > 0L) {
    pv[!is.na(pv) & pv == 0] <- .Machine$double.xmin
    .msgf("%d p-value%s underflowed to zero and %s floored at %.3g.",
          n_floor, if (n_floor > 1) "s" else "",
          if (n_floor > 1) "were" else "was", .Machine$double.xmin)
  }

  x[, `:=`(zbar = zbar, p = pv, css = -log10(pv))]
  if (na_action == "pairwise" && any(m_j < m_all)) {
    x[, n_tests := m_j]
  }

  if (!inherits(x, "css_result")) {
    data.table::setattr(x, "class", c("css_result", class(x)))
  }
  data.table::setattr(x, "css_call", list(
    ties = ties, na_action = na_action,
    weights = weights, calibrate = isTRUE(calibrate), calibration = cal,
    m = m_all, n = n, tests = tests
  ))
  x[]
}
