# ---------------------------------------------------------------------------
# Assemble the shipped datasets and write them to data/.
#
# Size budget: the installed package must stay under CRAN's 5 MB. Statistics
# are rounded to six significant figures and everything is xz-compressed.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(usethis)
})

round_stats <- function(d, digits = 6) {
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  num <- setdiff(num, "pos")     # positions stay exact
  for (cl in num) set(d, j = cl, value = signif(d[[cl]], digits))
  d
}

build_truth <- function(P) {
  sw <- as.data.table(P$sweeps)
  tr <- P$drift_trap
  truth <- rbind(
    sw[, .(
      sweep_id = as.character(sweep_id),
      chr = chr, pos = pos, scenario = scenario,
      s = s, p0 = p0, target_freq = target_p, cohort = cohort,
      selected = TRUE
    )],
    data.table(
      sweep_id = "trap", chr = tr$chr, pos = (tr$start + tr$end) / 2,
      scenario = "founder-effect haplotype block (no selection)",
      s = 0, p0 = NA_real_, target_freq = NA_real_,
      cohort = sprintf("selected (breeds %s)", paste(range(tr$breeds), collapse = "-")),
      selected = FALSE
    )
  )
  truth[, t_generations := P$t_sel]
  truth[, ne_breed := P$ne_breed]
  truth[]
}

write_datasets <- function(panel, P) {
  css_sim <- copy(panel)
  round_stats(css_sim)
  setDF(css_sim); setDT(css_sim)

  css_sim_truth <- build_truth(P)

  # A small subset for examples and tests that must run in well under 5 s.
  # BTA-1 and BTA-2 carry sweeps 2 and 1 respectively.
  css_sim_small <- css_sim[chr %in% c("1", "2")]
  css_sim_small[, chr := droplevels(chr)]

  usethis::use_data(css_sim, overwrite = TRUE, compress = "xz")
  usethis::use_data(css_sim_small, overwrite = TRUE, compress = "xz")
  usethis::use_data(css_sim_truth, overwrite = TRUE, compress = "xz")

  sizes <- file.info(list.files("data", full.names = TRUE))["size"]
  message("data/ contents:")
  print(round(sizes / 1024))
  message(sprintf("total data/ = %.2f MB", sum(sizes) / 1024^2))
  invisible(list(css_sim = css_sim, css_sim_small = css_sim_small,
                 css_sim_truth = css_sim_truth))
}
