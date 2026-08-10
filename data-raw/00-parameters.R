# ---------------------------------------------------------------------------
# Every parameter of the cssig example simulation lives here.
# Nothing downstream hard-codes a constant; change it here and re-run run-all.R.
#
# See PLAN-simulation.R in the project root for the reasoning behind each value.
# ---------------------------------------------------------------------------

PARAMS <- list(

  ## --- genome -------------------------------------------------------------
  # UMD3.1 autosome lengths (bp). Total 2512.08 Mb, matching the assembly
  # figure reported by Randhawa et al. (2015).
  chrom_len = c(
    `1` = 158337067, `2` = 137060424, `3` = 121430405, `4` = 120829699,
    `5` = 121191424, `6` = 119458736, `7` = 112638659, `8` = 113384836,
    `9` = 105708250, `10` = 104305016, `11` = 107310763, `12` =  91163125,
    `13` =  84240350, `14` =  84648390, `15` =  85296676, `16` =  81724687,
    `17` =  75158596, `18` =  66004023, `19` =  64057457, `20` =  72042655,
    `21` =  71599096, `22` =  61435874, `23` =  52530062, `24` =  62714930,
    `25` =  42904170, `26` =  51681464, `27` =  45407902, `28` =  46312546,
    `29` =  51505224
  ),

  mu = 1.2e-8,          # mutation rate per bp per generation
  rec = 1e-8,           # recombination rate per bp per generation (1 cM/Mb)

  ## --- coalescent background (stage A, scrm) ------------------------------
  ne_anc = 2000,        # ancestral effective population size
  bott_start = 2000,    # generations ago the domestication bottleneck begins
  bott_len = 40,        # its duration in generations
  bott_ne = 400,        # effective size during the bottleneck
  chunk_bp = 1e7,       # simulate the coalescent in 10 Mb chunks (see README)

  ## --- breeds and cohorts -------------------------------------------------
  n_breed = 16L,
  n_selected = 8L,      # breeds 1-8 are the selected cohort, 9-16 the reference
  founder_hap = 100L,   # haplotypes drawn per breed to found it
  ne_breed = 150L,      # diploid effective size of each breed during selection
  t_sel = 40L,          # generations of forward simulation (~200 years)
  # Animals sampled per breed. 50 x 16 = 800 animals, 400 per cohort, which
  # sits inside the range of the source papers (P1 dataset A: 212 animals over
  # 14 breeds; P2: 1106 over 55 breeds). An earlier draft used 20/breed, which
  # left the analysis limited by sampling noise in F_ST rather than by the
  # strength of the sweeps. Sampling more of the already-simulated population
  # costs nothing and sharpens every constituent test equally, so it does not
  # tilt the CSS-versus-constituents comparison in the vignette.
  n_sample = 50L,

  ## --- SNP thinning and ascertainment -------------------------------------
  prethin_target = 150000L,  # sites carried into the forward simulation
  prethin_maf = 0.02,        # minimum founder MAF to be carried forward
  discovery_breeds = c(1L, 5L, 12L),  # artificial SNP-discovery panel
  discovery_n = 8L,          # animals per discovery breed
  discovery_maf = 0.05,      # MAF required in the discovery panel
  panel_target = 50000L,     # final SNP count
  panel_maf = 0.01,          # final QC MAF, as in the source papers

  ## --- sweeps -------------------------------------------------------------
  # Selection coefficients are large because the window is short: 40
  # generations of intense artificial selection, not slow natural selection.
  # With additive fitness (1, 1+s/2, 1+s) the time to travel from p0 to p1 is
  # (2/s) * log[(p1/(1-p1)) / (p0/(1-p0))], so a complete sweep from a single
  # copy inside 40 generations needs s of roughly 0.5. Values below were solved
  # from that relation for the target frequencies.
  sweeps = data.frame(
    sweep_id  = 1:6,
    chr       = c("2", "1", "6", "14", "10", "13"),
    pos       = c(60e6, 95e6, 40e6, 30e6, 55e6, 45e6),
    scenario  = c("complete hard sweep", "incomplete hard sweep",
                  "soft sweep from standing variation", "weak sweep on standing variation",
                  "breed-restricted hard sweep", "sweep in the reference cohort"),
    s         = c(0.50, 0.30, 0.20, 0.10, 0.45, 0.45),
    p0        = c(NA, NA, 0.10, 0.15, NA, NA),   # NA = new mutation at 1/(2Ne)
    target_p  = c(0.99, 0.70, 0.85, 0.60, 0.98, 0.98),
    cohort    = c("selected", "selected", "selected", "selected",
                  "selected (4 of 8 breeds)", "reference"),
    stringsAsFactors = FALSE
  ),

  # Deliberate false positive: a founder-effect haplotype block. In four
  # selected-cohort breeds the whole region is descended from a single founder
  # haplotype, which is what happens when a breed is founded from few animals.
  # It produces high F_ST and long-range homozygosity with no selection at all.
  drift_trap = list(chr = "5", start = 40e6, end = 45e6, breeds = 1:4),

  ## --- reproducibility ----------------------------------------------------
  seed_background = 20240809L,
  seed_forward    = 20240810L,
  seed_ascertain  = 20240811L,

  cache_dir = file.path("data-raw", "cache")
)

PARAMS$selected_breeds  <- seq_len(PARAMS$n_selected)
PARAMS$reference_breeds <- seq.int(PARAMS$n_selected + 1L, PARAMS$n_breed)
PARAMS$genome_bp        <- sum(PARAMS$chrom_len)
