#' Read XP-EHH output from selscan
#'
#' Reads a selscan `.xpehh.out` file, or the normalised `.xpehh.out.norm`
#' produced by `norm`, into the shape [css_input()] expects.
#'
#' @param file Path to the selscan output file.
#' @param chr Chromosome label to attach. selscan output covers one chromosome
#'   and does not record which, so this must be supplied.
#' @param normalised Set `TRUE` when reading a `.norm` file, which carries the
#'   extra `normxpehh` and `crit` columns. If `NULL` (default) the format is
#'   detected from the header.
#'
#' @return A [data.table::data.table] with `chr`, `pos`, `snp` and `xpehh`.
#'
#' @examples
#' f <- system.file("extdata", "example.xpehh.out", package = "cssig")
#' if (nzchar(f)) read_selscan_xpehh(f, chr = 1)
#'
#' @export
read_selscan_xpehh <- function(file, chr, normalised = NULL) {
  if (!file.exists(file)) .stopf("File not found: %s", file)
  if (missing(chr)) .stopf("`chr` must be supplied; selscan output does not record the chromosome.")

  d <- data.table::fread(file)
  nm <- tolower(names(d))

  # selscan writes: id pos gpos p1 ihh1 p2 ihh2 xpehh  (+ normxpehh crit if normalised)
  has_norm <- any(grepl("^normxpehh", nm))
  if (is.null(normalised)) normalised <- has_norm
  if (normalised && !has_norm) {
    .stopf(paste0("`normalised = TRUE` but no `normxpehh` column found; this does not ",
                  "look like a selscan `norm` output file. Columns: %s."),
           paste(names(d), collapse = ", "))
  }
  val_col <- if (normalised) {
    names(d)[grep("^normxpehh", nm)[1]]
  } else if (any(nm == "xpehh")) {
    names(d)[which(nm == "xpehh")[1]]
  } else {
    .stopf("Could not find an `xpehh` or `normxpehh` column. Columns: %s.",
           paste(names(d), collapse = ", "))
  }

  id_col  <- names(d)[which(nm %in% c("id", "locus", "snp"))[1]]
  pos_col <- names(d)[which(nm %in% c("pos", "physpos", "position"))[1]]
  if (is.na(pos_col)) .stopf("Could not find a position column. Columns: %s.",
                             paste(names(d), collapse = ", "))

  out <- data.table::data.table(
    chr = as.character(chr),
    pos = as.numeric(d[[pos_col]]),
    snp = if (!is.na(id_col)) as.character(d[[id_col]]) else NA_character_,
    xpehh = as.numeric(d[[val_col]])
  )
  data.table::setkeyv(out, c("chr", "pos"))
  out[]
}

#' Coerce rehh XP-EHH output for use with CSS
#'
#' Takes the data frame returned by [rehh::ies2xpehh()] and renames its columns
#' to the shape [css_input()] expects.
#'
#' @param x A data frame from `rehh::ies2xpehh()`.
#' @param value Which column to use: `"XPEHH"` (default) picks the log-ratio
#'   column, whose exact name depends on the population labels and on whether
#'   `ies2xpehh(standardize = )` was `TRUE` (`XPEHH_a_b`) or `FALSE`
#'   (`UNXPEHH_a_b`); both spellings are recognised. `"LOGPVALUE"` uses rehh's
#'   -log10(p) column instead.
#'
#' @return A [data.table::data.table] with `chr`, `pos`, `snp` and `xpehh`.
#'
#' @examples
#' fake <- data.frame(CHR = 1, POSITION = c(100, 200),
#'                    XPEHH_A_B = c(1.2, -0.4))
#' read_rehh_xpehh(fake)
#'
#' @export
read_rehh_xpehh <- function(x, value = c("XPEHH", "LOGPVALUE")) {
  value <- match.arg(value)
  d <- data.table::as.data.table(x)
  nm <- toupper(names(d))

  chr_col <- names(d)[which(nm %in% c("CHR", "CHROMOSOME"))[1]]
  pos_col <- names(d)[which(nm %in% c("POSITION", "POS"))[1]]
  if (is.na(chr_col) || is.na(pos_col)) {
    .stopf("Could not find CHR and POSITION columns. Columns: %s.",
           paste(names(d), collapse = ", "))
  }

  # rehh names the statistic column XPEHH_<pop1>_<pop2> when ies2xpehh() was
  # called with standardize = TRUE and UNXPEHH_<pop1>_<pop2> when it was not,
  # so the pattern has to admit both spellings.
  val_col <- if (value == "LOGPVALUE") {
    names(d)[grep("^LOGPVALUE", nm)[1]]
  } else {
    names(d)[grep("^(UN)?XPEHH", nm)[1]]
  }
  if (is.na(val_col)) {
    .stopf("Could not find a %s column. Columns: %s.", value,
           paste(names(d), collapse = ", "))
  }

  snp_col <- if (!is.null(rownames(x)) && !identical(rownames(x), as.character(seq_len(nrow(x))))) {
    rownames(x)
  } else NA_character_

  out <- data.table::data.table(
    chr = as.character(d[[chr_col]]),
    pos = as.numeric(d[[pos_col]]),
    snp = if (length(snp_col) > 1) snp_col else NA_character_,
    xpehh = as.numeric(d[[val_col]])
  )
  data.table::setkeyv(out, c("chr", "pos"))
  out[]
}
