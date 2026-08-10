#' Weir and Cockerham's F_ST between two cohorts
#'
#' Per-SNP \eqn{\theta} of Weir & Cockerham (1984) for two populations, the
#' estimator used for the \eqn{F_{ST}} constituent of CSS.
#'
#' @details
#' Inputs are allele counts, not frequencies, because the estimator needs the
#' sample sizes: it corrects for the finite number of individuals sampled, which
#' a frequency-only calculation cannot do. Supply the count of one allele and
#' the number of genotyped individuals in each cohort.
#'
#' Negative values are a normal outcome of the estimator when true
#' differentiation is near zero. They are returned as-is by default; set
#' `floor_zero = TRUE` to truncate at zero. Truncation biases the genome-wide
#' distribution upward, but because CSS is rank-based and truncation is monotone
#' it has no effect on CSS itself.
#'
#' @param count1,count2 Integer vector of copies of the reference allele in
#'   cohort 1 and cohort 2.
#' @param n1,n2 Number of *individuals* genotyped in each cohort, per SNP
#'   (scalars are recycled).
#' @param ploidy Ploidy, default `2`.
#' @param floor_zero Truncate negative estimates at zero. Default `FALSE`.
#'
#' @return Numeric vector of \eqn{F_{ST}} estimates, `NA` where a SNP is
#'   monomorphic across both cohorts.
#'
#' @references Weir BS, Cockerham CC (1984). Estimating F-statistics for the
#'   analysis of population structure. *Evolution* 38:1358-1370.
#'
#' @examples
#' css_fst(count1 = c(10, 40, 25), n1 = 25,
#'         count2 = c(40, 8, 26),  n2 = 25)
#'
#' @export
css_fst <- function(count1, n1, count2, n2, ploidy = 2L, floor_zero = FALSE) {
  n <- max(length(count1), length(count2))
  count1 <- rep_len(count1, n); count2 <- rep_len(count2, n)
  n1 <- rep_len(n1, n);         n2 <- rep_len(n2, n)

  # Sample sizes in individuals; allele counts out of ploidy * n
  p1 <- count1 / (ploidy * n1)
  p2 <- count2 / (ploidy * n2)

  r <- 2
  nbar <- (n1 + n2) / r
  nc <- (r * nbar - (n1^2 + n2^2) / (r * nbar)) / (r - 1)
  pbar <- (n1 * p1 + n2 * p2) / (r * nbar)
  s2 <- (n1 * (p1 - pbar)^2 + n2 * (p2 - pbar)^2) / ((r - 1) * nbar)

  # Heterozygosity term. Without individual genotypes we cannot observe hbar,
  # so assume Hardy-Weinberg within cohorts, which is the standard substitution
  # when working from allele counts.
  hbar <- (n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / (r * nbar)

  a <- (nbar / nc) * (s2 - (1 / (nbar - 1)) *
        (pbar * (1 - pbar) - ((r - 1) / r) * s2 - hbar / 4))
  b <- (nbar / (nbar - 1)) *
        (pbar * (1 - pbar) - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  cc <- hbar / 2

  denom <- a + b + cc
  fst <- ifelse(denom == 0 | !is.finite(denom), NA_real_, a / denom)
  if (floor_zero) fst <- pmax(fst, 0)
  fst
}

#' Change in derived allele frequency between cohorts
#'
#' \eqn{\Delta}DAF, the difference in derived allele frequency between the
#' putatively selected cohort and the reference cohort, as used by Grossman
#' et al. (2010) and as a CSS constituent in Randhawa et al. (2014).
#'
#' @param daf_selected,daf_reference Derived allele frequency in each cohort.
#' @param standardize Rescale the result to mean 0 and unit variance, as both
#'   source papers do. Note this has no effect on CSS, which is rank-based; it
#'   matters only when plotting constituent tests on a shared axis.
#'
#' @return Numeric vector of \eqn{\Delta}DAF values.
#'
#' @examples
#' css_ddaf(c(0.8, 0.2, 0.5), c(0.3, 0.25, 0.5))
#'
#' @export
css_ddaf <- function(daf_selected, daf_reference, standardize = FALSE) {
  d <- daf_selected - daf_reference
  if (standardize) d <- css_standardize(d)
  d
}

#' Change in selected allele frequency between cohorts
#'
#' \eqn{\Delta}SAF, the ancestral-allele-free alternative to \eqn{\Delta}DAF
#' introduced by Randhawa et al. (2014) for datasets where the ancestral allele
#' cannot be inferred, such as their sheep panels. It is the frequency
#' difference for the allele that is the *major* allele in the selected cohort.
#'
#' @param freq_selected,freq_reference Frequency of the same reference allele in
#'   each cohort.
#' @param standardize Rescale to mean 0 and unit variance. See [css_ddaf()].
#'
#' @return Numeric vector of \eqn{\Delta}SAF values.
#'
#' @examples
#' css_dsaf(c(0.8, 0.2, 0.5), c(0.3, 0.25, 0.5))
#'
#' @export
css_dsaf <- function(freq_selected, freq_reference, standardize = FALSE) {
  # Orient each SNP on the allele that is major in the selected cohort.
  flip <- !is.na(freq_selected) & freq_selected < 0.5
  fs <- ifelse(flip, 1 - freq_selected, freq_selected)
  fr <- ifelse(flip, 1 - freq_reference, freq_reference)
  d <- fs - fr
  if (standardize) d <- css_standardize(d)
  d
}

#' Standardise a statistic to mean zero and unit variance
#'
#' Both source papers standardise XP-EHH and \eqn{\Delta}DAF this way. It leaves
#' CSS unchanged, because CSS depends only on ranks; it is done so that
#' constituent tests share an axis when plotted.
#'
#' @param x Numeric vector.
#' @param na.rm Ignore missing values when computing the mean and standard
#'   deviation. Default `TRUE`.
#'
#' @return Numeric vector of z-scores.
#'
#' @examples
#' css_standardize(c(1, 4, 9, 16))
#'
#' @export
css_standardize <- function(x, na.rm = TRUE) {
  s <- stats::sd(x, na.rm = na.rm)
  if (is.na(s) || s == 0) {
    .warnf("Standard deviation is zero or undefined; returning centred values only.")
    return(x - mean(x, na.rm = na.rm))
  }
  (x - mean(x, na.rm = na.rm)) / s
}
