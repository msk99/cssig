# ---------------------------------------------------------------------------
# Stage A: neutral coalescent background (scrm)
# Stage B: forward Wright-Fisher breed formation, with sweeps
#
# Produces one cached .rds per chromosome holding the final phased haplotypes
# for all 16 breeds, plus the site positions and ancestral/derived coding.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(scrm)
  library(data.table)
})

# --- stage A ---------------------------------------------------------------

#' Simulate neutral founder haplotypes for one chromosome.
#'
#' The chromosome is simulated in `chunk_bp` pieces to keep the segregating-site
#' matrix within memory. Chunks are independent, so linkage disequilibrium does
#' not carry across a chunk boundary. With 10 Mb chunks and cattle LD decaying
#' within roughly 1-2 Mb this affects only the windows straddling a boundary.
simulate_background <- function(chr, P) {
  L <- P$chrom_len[[chr]]
  nsam <- P$n_breed * P$founder_hap
  n_chunk <- ceiling(L / P$chunk_bp)

  # Demography in scrm's coalescent time units of 4 * ne_anc generations.
  t_bott_start <- P$bott_start / (4 * P$ne_anc)
  t_bott_end   <- (P$bott_start + P$bott_len) / (4 * P$ne_anc)
  demog <- sprintf("-eN %.6f %.6f -eN %.6f 1.0",
                   t_bott_start, P$bott_ne / P$ne_anc, t_bott_end)

  pos_all <- numeric(0)
  hap_all <- vector("list", n_chunk)

  for (k in seq_len(n_chunk)) {
    chunk_start <- (k - 1) * P$chunk_bp
    chunk_len <- min(P$chunk_bp, L - chunk_start)
    if (chunk_len <= 0) next

    theta <- 4 * P$ne_anc * P$mu  * chunk_len
    rho   <- 4 * P$ne_anc * P$rec * chunk_len

    out <- scrm::scrm(sprintf("%d 1 -t %.4f -r %.4f %d %s",
                              nsam, theta, rho, round(chunk_len), demog))
    ss <- out$seg_sites[[1]]
    if (is.null(ss) || ncol(ss) == 0L) next

    # scrm reports positions on [0, 1] within the chunk; 0 = ancestral allele.
    # Scale to base pairs and round: a physical position is an integer number
    # of bases, and leaving it continuous produces nonsense like 31501.53 bp
    # and SNP identifiers to match.
    p <- round(as.numeric(colnames(ss)) * chunk_len + chunk_start)

    # Two segregating sites can round onto the same base pair; keep one, since
    # downstream every SNP must occupy a distinct (chr, pos).
    dup <- duplicated(p)
    if (any(dup)) {
      ss <- ss[, !dup, drop = FALSE]
      p <- p[!dup]
    }
    storage.mode(ss) <- "integer"

    # Thin inside the chunk immediately, so peak memory is one chunk.
    keep <- prethin_sites(ss, p, P, chunk_len)
    hap_all[[k]] <- ss[, keep, drop = FALSE]
    pos_all <- c(pos_all, p[keep])
    rm(ss); gc(verbose = FALSE)
  }

  hap <- do.call(cbind, hap_all[!vapply(hap_all, is.null, logical(1))])
  list(hap = hap, pos = pos_all)   # hap: nsam x nsites, 0 = ancestral
}

# Drop very rare sites and thin to roughly uniform spacing, so the forward
# simulation carries a manageable number of sites.
prethin_sites <- function(ss, p, P, chunk_len) {
  af <- colMeans(ss)
  maf <- pmin(af, 1 - af)
  ok <- which(maf >= P$prethin_maf)
  if (!length(ok)) return(integer(0))

  target <- max(1L, round(P$prethin_target * chunk_len / P$genome_bp))
  if (length(ok) <= target) return(ok)

  # Even spatial coverage: one site per equal-width bin, chosen at random.
  # Deliberately *not* the commonest site in the bin -- biasing the frequency
  # spectrum here would pre-empt the discovery-panel ascertainment in
  # 02-constituent-stats.R, which is where chip-like bias belongs.
  bins <- cut(p[ok], breaks = target, labels = FALSE)
  idx <- tapply(seq_along(ok), bins, function(i) ok[if (length(i) == 1L) i else sample(i, 1L)])
  sort(unlist(idx, use.names = FALSE))
}

# --- stage B ---------------------------------------------------------------

#' One Wright-Fisher generation with recombination and optional selection.
#'
#' @param H integer matrix, nsites x nhap, 0/1. Haplotypes 2i-1 and 2i form
#'   individual i.
#' @param pos site positions in bp, ascending.
#' @param n_out number of haplotypes to produce (must be even).
#' @param s selection coefficient; fitness is additive, (1, 1+s/2, 1+s).
#' @param focal row index of the site under selection, or NA.
#' @param rec recombination rate in Morgans per bp (P$rec).
wf_generation <- function(H, pos, n_out, s = 0, focal = NA_integer_, rec = 1e-8) {
  S <- nrow(H); N <- ncol(H); nind <- N %/% 2L

  if (!is.na(focal) && s != 0) {
    g <- H[focal, seq.int(1L, N, by = 2L)] + H[focal, seq.int(2L, N, by = 2L)]
    w <- 1 + s * g / 2
  } else {
    w <- NULL
  }

  par <- if (is.null(w)) sample.int(nind, n_out, replace = TRUE)
         else sample.int(nind, n_out, replace = TRUE, prob = w)
  pA <- 2L * par - 1L
  pB <- 2L * par

  # Crossovers: Poisson in number, uniform in position. Far cheaper than
  # drawing a Bernoulli per site, and identical in distribution.
  span_morgan <- (pos[S] - pos[1]) * rec
  nxo <- stats::rpois(n_out, span_morgan)
  tot <- sum(nxo)

  if (tot > 0L) {
    col <- rep.int(seq_len(n_out), nxo)
    xpos <- stats::runif(tot, pos[1], pos[S])
    row <- pmin(pmax(findInterval(xpos, pos) + 1L, 1L), S)
    cnt <- tabulate((col - 1L) * S + row, nbins = S * n_out)
    cs <- cumsum(cnt)
    cs <- matrix(cs, S, n_out)
    if (n_out > 1L) cs <- cs - rep(c(0, cs[S, -n_out]), each = S)
    phase <- cs %% 2L                      # 0 = from parent A, 1 = from parent B
    src <- ifelse(phase == 0L, rep(pA, each = S), rep(pB, each = S))
  } else {
    src <- rep(pA, each = S)
  }

  lin <- (as.integer(src) - 1L) * S + rep.int(seq_len(S), n_out)
  matrix(H[lin], S, n_out)
}

#' Forward-simulate all breeds on one chromosome.
#'
#' Every breed runs for the same number of generations at the same effective
#' size, whether or not it carries a sweep. Without that, breeds carrying a
#' sweep would accumulate extra drift and F_ST would be inflated across the
#' whole region rather than at the sweep.
forward_chromosome <- function(bg, chr, P, sweeps, trap) {
  S <- length(bg$pos)
  pos <- bg$pos
  n_hap <- P$ne_breed * 2L

  sw <- sweeps[sweeps$chr == chr, , drop = FALSE]
  focal <- NA_integer_
  if (nrow(sw)) {
    focal <- which.min(abs(pos - sw$pos[1]))
  }

  breeds <- vector("list", P$n_breed)
  restarts <- 0L

  for (b in seq_len(P$n_breed)) {
    idx <- ((b - 1L) * P$founder_hap + 1L):(b * P$founder_hap)
    H0 <- t(bg$hap[idx, , drop = FALSE])       # nsites x founder_hap

    # Founder-effect haplotype block: the deliberate neutral false positive.
    if (!is.null(trap) && trap$chr == chr && b %in% trap$breeds) {
      inblock <- which(pos >= trap$start & pos <= trap$end)
      if (length(inblock)) {
        donor <- sample.int(ncol(H0), 1L)
        H0[inblock, ] <- H0[inblock, donor]
      }
    }

    s_b <- 0
    carries <- FALSE
    if (nrow(sw) && !is.na(focal)) {
      carries <- switch(
        sw$cohort[1],
        "selected"                  = b %in% P$selected_breeds,
        "selected (4 of 8 breeds)"  = b %in% P$selected_breeds[1:4],
        "reference"                 = b %in% P$reference_breeds,
        FALSE
      )
      if (carries) s_b <- sw$s[1] else H0[focal, ] <- 0L
    }

    # A beneficial mutation entering at a single copy is lost with probability
    # roughly 1 - s even when selection is strong, so a run that loses it says
    # nothing about what a sweep looks like. Condition on non-loss by retrying,
    # and record how many retries that took: conditioning is itself a modelling
    # choice and should be visible in the provenance, not hidden.
    max_try <- 200L
    for (attempt in seq_len(max_try)) {
      H <- H0
      if (carries) {
        H[focal, ] <- 0L
        n_carry <- if (is.na(sw$p0[1])) 1L else max(1L, round(sw$p0[1] * ncol(H)))
        H[focal, sample.int(ncol(H), n_carry)] <- 1L
      }
      H <- wf_generation(H, pos, n_hap, s_b, focal, rec = P$rec)
      for (g in seq_len(P$t_sel - 1L)) {
        H <- wf_generation(H, pos, n_hap, s_b, focal, rec = P$rec)
      }
      if (!carries || mean(H[focal, ]) > 0) break
      restarts <- restarts + 1L
      if (attempt == max_try) {
        warning(sprintf("chr %s breed %d: sweep lost in %d attempts; keeping last run.",
                        chr, b, max_try))
      }
    }

    take <- sort(sample.int(ncol(H), P$n_sample * 2L))
    breeds[[b]] <- H[, take, drop = FALSE]
    rm(H, H0)
  }

  list(hap = do.call(cbind, breeds),   # nsites x (n_breed * 2 * n_sample)
       pos = pos,
       focal = focal,
       restarts = restarts,
       breed = rep(seq_len(P$n_breed), each = P$n_sample * 2L))
}

# --- driver ----------------------------------------------------------------

simulate_all <- function(P, chroms = names(P$chrom_len), verbose = TRUE) {
  dir.create(P$cache_dir, showWarnings = FALSE, recursive = TRUE)
  timings <- data.table()

  for (chr in chroms) {
    f <- file.path(P$cache_dir, sprintf("chr%s.rds", chr))
    if (file.exists(f)) {
      if (verbose) message(sprintf("chr %-3s cached, skipping", chr))
      next
    }
    t0 <- Sys.time()
    set.seed(P$seed_background + as.integer(chr))
    bg <- simulate_background(chr, P)

    set.seed(P$seed_forward + as.integer(chr))
    fw <- forward_chromosome(bg, chr, P, P$sweeps, P$drift_trap)
    rm(bg); gc(verbose = FALSE)

    saveRDS(fw, f, compress = "xz")
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    timings <- rbind(timings, data.table(chr = chr, sites = length(fw$pos), secs = dt))
    if (verbose) {
      message(sprintf("chr %-3s %6d sites  %5.1f s", chr, length(fw$pos), dt))
    }
    rm(fw); gc(verbose = FALSE)
  }
  timings
}
