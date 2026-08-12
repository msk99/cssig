# ---------------------------------------------------------------------------
# SNP ascertainment and computation of the three constituent statistics.
#
# Ascertainment mimics a SNP chip: variants are discovered in a small panel of
# breeds and then thinned to roughly uniform spacing. This matters because the
# source papers analyse BovineSNP50 genotypes, whose frequency spectrum is
# strongly shifted towards common variants; simulating sequence-like data would
# give the method an easier problem than it actually faces.
#
# F_ST, dDAF and dSAF are computed with the package's own exported functions,
# so building the example data also exercises them end to end. XP-EHH comes
# from rehh, on the perfectly phased simulated haplotypes.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rehh)
})

#' Chip-style ascertainment for one chromosome.
ascertain <- function(fw, P) {
  n_site <- length(fw$pos)
  hap <- fw$hap
  breed <- fw$breed

  # 1. discovery panel: a few breeds, a few animals each
  disc_hap <- unlist(lapply(P$discovery_breeds, function(b) {
    idx <- which(breed == b)
    idx[seq_len(min(length(idx), P$discovery_n * 2L))]
  }))
  af_disc <- rowMeans(hap[, disc_hap, drop = FALSE])
  maf_disc <- pmin(af_disc, 1 - af_disc)
  discovered <- which(maf_disc >= P$discovery_maf)

  # 2. thin the discovered sites to roughly uniform spacing
  target <- max(1L, round(P$panel_target * P$chrom_len[[fw$chr]] / P$genome_bp))
  if (length(discovered) > target) {
    bins <- cut(fw$pos[discovered], breaks = target, labels = FALSE)
    discovered <- sort(vapply(split(discovered, bins),
                              function(i) if (length(i) == 1L) i else sample(i, 1L),
                              integer(1)))
  }

  # 3. final QC on the full panel, as in the source papers
  af_all <- rowMeans(hap[, , drop = FALSE])
  maf_all <- pmin(af_all, 1 - af_all)
  keep <- discovered[maf_all[discovered] >= P$panel_maf]

  # 4. drop the causal site itself. Real SNP chips almost never carry the
  #    functional variant, and both source papers detect regions flanking the
  #    known genes rather than the mutations themselves.
  if (!is.na(fw$focal)) keep <- setdiff(keep, fw$focal)

  sort(keep)
}

#' Constituent statistics for one chromosome.
chromosome_stats <- function(fw, P, chr) {
  keep <- ascertain(fw, P)
  if (length(keep) < 10L) return(NULL)

  hap <- fw$hap[keep, , drop = FALSE]
  pos <- fw$pos[keep]
  sel <- fw$breed %in% P$selected_breeds
  ref <- !sel

  n_sel_ind <- sum(sel) / 2
  n_ref_ind <- sum(ref) / 2

  cnt_sel <- rowSums(hap[, sel, drop = FALSE])   # copies of the derived allele
  cnt_ref <- rowSums(hap[, ref, drop = FALSE])
  daf_sel <- cnt_sel / sum(sel)
  daf_ref <- cnt_ref / sum(ref)

  fst <- cssig::css_fst(cnt_sel, n_sel_ind, cnt_ref, n_ref_ind)
  ddaf <- cssig::css_ddaf(daf_sel, daf_ref)
  dsaf <- cssig::css_dsaf(daf_sel, daf_ref)

  # --- XP-EHH via rehh -----------------------------------------------------
  # rehh wants haplotypes in rows and markers in columns, an integer matrix
  # coded 0 = ancestral / 1 = derived, and positions as doubles.
  mk <- paste0("chr", chr, "_", format(pos, scientific = FALSE, trim = TRUE))
  build <- function(cols) {
    m <- t(hap[, cols, drop = FALSE])
    storage.mode(m) <- "integer"
    colnames(m) <- mk
    new("haplohh", haplo = m, positions = as.numeric(pos), chr.name = as.character(chr))
  }
  s_sel <- rehh::scan_hh(build(which(sel)), polarized = TRUE,
                         discard_integration_at_border = FALSE)
  s_ref <- rehh::scan_hh(build(which(ref)), polarized = TRUE,
                         discard_integration_at_border = FALSE)
  xp <- rehh::ies2xpehh(s_sel, s_ref, popname1 = "selected", popname2 = "reference",
                        standardize = FALSE, verbose = FALSE)
  # rehh names the column UNXPEHH_* when standardize = FALSE and XPEHH_* when
  # TRUE, so the pattern must admit both. Standardising happens once,
  # genome-wide, in build_panel() rather than per chromosome.
  xcol <- grep("^(UN)?XPEHH", names(xp), value = TRUE)[1]
  stopifnot(!is.na(xcol), nrow(xp) == length(pos))

  out <- data.table(
    chr = as.character(chr),
    pos = as.numeric(pos),
    snp = mk,
    ancestral = "A", derived = "G",
    maf = pmin(rowMeans(hap), 1 - rowMeans(hap)),
    daf_selected = daf_sel,
    daf_reference = daf_ref,
    fst = fst,
    xpehh = as.numeric(xp[[xcol]]),
    ddaf = ddaf,
    dsaf = dsaf
  )
  stopifnot(nrow(out) == length(keep))
  out[]
}

#' Assemble all chromosomes and standardise the directional statistics.
build_panel <- function(P, chroms = names(P$chrom_len), verbose = TRUE) {
  set.seed(P$seed_ascertain)
  parts <- vector("list", length(chroms))
  for (i in seq_along(chroms)) {
    chr <- chroms[i]
    f <- file.path(P$cache_dir, sprintf("chr%s.rds", chr))
    if (!file.exists(f)) { warning("missing cache for chr ", chr); next }
    fw <- readRDS(f)
    fw$chr <- chr
    parts[[i]] <- chromosome_stats(fw, P, chr)
    if (verbose) message(sprintf("chr %-3s -> %5d SNPs", chr, nrow(parts[[i]])))
    rm(fw); gc(verbose = FALSE)
  }
  d <- rbindlist(parts)

  # Standardised genome-wide, as both source papers do. This has no effect on
  # CSS, which depends only on ranks; it makes the constituent tests
  # comparable when plotted on a shared axis.
  d[, xpehh := cssig::css_standardize(xpehh)]
  d[, ddaf  := cssig::css_standardize(ddaf)]
  d[, dsaf  := cssig::css_standardize(dsaf)]

  d[, chr := factor(chr, levels = names(P$chrom_len))]
  setkeyv(d, c("chr", "pos"))
  d[]
}
