# ---------------------------------------------------------------------------
# Rebuild the cssig example datasets from scratch.
#
#   Rscript data-raw/run-all.R
#
# Run from the package root. Stage 1 caches one .rds per chromosome in
# data-raw/cache/, so an interrupted run resumes where it stopped; delete the
# cache to force a full rebuild.
# ---------------------------------------------------------------------------

stopifnot(file.exists("DESCRIPTION"))
t_start <- Sys.time()

source("data-raw/00-parameters.R")
source("data-raw/01-simulate-genome.R")
source("data-raw/02-constituent-stats.R")
source("data-raw/03-package-data.R")

pkgload::load_all(".", quiet = TRUE)   # for css_fst / css_ddaf / css_dsaf

message("== stage 1: simulate genomes ==")
timings <- simulate_all(PARAMS)

message("== stage 2: ascertain and compute constituent statistics ==")
panel <- build_panel(PARAMS)

message("== stage 3: write package data ==")
out <- write_datasets(panel, PARAMS)

# --- provenance ------------------------------------------------------------
restarts <- vapply(names(PARAMS$chrom_len), function(ch) {
  f <- file.path(PARAMS$cache_dir, sprintf("chr%s.rds", ch))
  if (file.exists(f)) as.integer(readRDS(f)$restarts %||% 0L) else NA_integer_
}, integer(1))

sink("data-raw/provenance.txt")
cat("cssig example data provenance\n")
cat("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("runtime  :", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), "\n\n")
cat("RNGkind:", paste(RNGkind(), collapse = ", "), "\n")
cat("seeds  : background", PARAMS$seed_background,
    "forward", PARAMS$seed_forward, "ascertain", PARAMS$seed_ascertain, "\n\n")
cat("SNPs in css_sim      :", nrow(out$css_sim), "\n")
cat("SNPs in css_sim_small:", nrow(out$css_sim_small), "\n\n")
cat("sweep restarts per chromosome (conditioning on non-loss):\n")
print(restarts[restarts > 0])
cat("\nper-chromosome timings:\n"); print(timings)
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

message(sprintf("done in %.1f min; see data-raw/provenance.txt",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
