#' Simulated multi-breed cattle panel with known selection signatures
#'
#' **These data are simulated. They describe no real animal, breed or genome.**
#'
#' Constituent selection statistics for a simulated multi-breed cattle panel,
#' computed from simulated genotypes rather than drawn directly from a
#' distribution. That distinction matters: \eqn{F_{ST}}, XP-EHH and
#' \eqn{\Delta}DAF are correlated with each other because they read the same
#' underlying genealogy, and CSS exists precisely to exploit that shared
#' structure. Statistics simulated independently would carry a correlation
#' structure chosen by the package author, and any demonstration built on them
#' would be circular.
#'
#' @format A [data.table::data.table] with one row per SNP:
#' \describe{
#'   \item{chr}{Chromosome, an ordered factor with the 29 bovine autosomes.}
#'   \item{pos}{Position in base pairs on a UMD3.1-shaped coordinate system.}
#'   \item{snp}{SNP identifier.}
#'   \item{ancestral, derived}{Allele coding. The ancestral allele is known
#'     exactly, because the simulation records it.}
#'   \item{maf}{Minor allele frequency across the whole panel.}
#'   \item{daf_selected, daf_reference}{Derived allele frequency in each cohort.}
#'   \item{fst}{Weir & Cockerham \eqn{F_{ST}} between the cohorts, from
#'     [css_fst()].}
#'   \item{xpehh}{Cross-population extended haplotype homozygosity, selected
#'     against reference, computed with \pkg{rehh} and standardised.}
#'   \item{ddaf}{Standardised \eqn{\Delta}DAF, from [css_ddaf()].}
#'   \item{dsaf}{Standardised \eqn{\Delta}SAF, from [css_dsaf()]. Provided so
#'     that the ancestral-allele-free route of Randhawa et al. (2014) can be
#'     compared against \eqn{\Delta}DAF on the same data.}
#' }
#'
#' @section How it was generated:
#' Sixteen breeds descend from a common ancestral population that passes
#' through a domestication bottleneck, simulated with the sequentially-Markovian
#' coalescent (\pkg{scrm}). Breeds are then formed by a forward Wright-Fisher
#' simulation with recombination, run for 40 generations at an effective size of
#' 150, during which the sweeps in [css_sim_truth] act. Eight breeds form the
#' selected cohort and eight the reference cohort, with 20 animals sampled per
#' breed. SNPs are ascertained chip-style, by discovery in three breeds followed
#' by thinning to roughly uniform spacing, which reproduces the shifted
#' frequency spectrum of the BovineSNP50 data the source papers analyse.
#'
#' Every breed is run forward for the same number of generations at the same
#' effective size, whether or not it carries a sweep. Without that, the swept
#' breeds would accumulate extra drift and \eqn{F_{ST}} would be raised across
#' the whole region rather than at the sweep.
#'
#' The causal variants themselves are excluded from the panel, as they would be
#' on a real SNP chip; the signal comes from hitchhiking neighbours.
#'
#' Full parameters and the generating scripts are in `data-raw/` in the package
#' source.
#'
#' @seealso [css_sim_truth] for the known answers, [css_sim_small] for a subset
#'   sized for quick examples.
#' @examples
#' data(css_sim)
#' str(css_sim)
"css_sim"

#' Known selection signatures in the simulated panel
#'
#' The ground truth behind [css_sim]: where the sweeps are, how strong they
#' were, and which cohort carried them. Because the answers are known, this
#' supports genuine power and false-positive comparisons of CSS against its
#' constituent tests.
#'
#' @format A [data.table::data.table] with one row per implanted signature:
#' \describe{
#'   \item{sweep_id}{Identifier; `"trap"` marks the deliberate false positive.}
#'   \item{chr, pos}{Location of the causal site.}
#'   \item{scenario}{What kind of signature this is.}
#'   \item{s}{Selection coefficient, additive with fitnesses
#'     \eqn{1, 1+s/2, 1+s}.}
#'   \item{p0}{Starting frequency; `NA` means a new mutation at \eqn{1/(2N_e)}.}
#'   \item{target_freq}{Frequency the sweep was designed to reach.}
#'   \item{cohort}{Which breeds were under selection.}
#'   \item{selected}{`FALSE` for the trap, `TRUE` for genuine sweeps.}
#'   \item{t_generations, ne_breed}{Duration and effective size of the forward
#'     phase.}
#' }
#'
#' @section Why the selection coefficients are large:
#' The forward phase runs for only 40 generations, roughly 200 years of breed
#' formation. Under additive fitness an allele needs
#' \eqn{t = (2/s)\log[(p_1/(1-p_1))/(p_0/(1-p_0))]} generations to travel from
#' \eqn{p_0} to \eqn{p_1}, so a complete sweep from a single copy inside 40
#' generations requires \eqn{s} near 0.5. These are values appropriate to
#' intense artificial selection, not to slow natural selection.
#'
#' @section The trap:
#' Row `"trap"` carries no selection at all. In four selected-cohort breeds the
#' whole 40-45 Mb region of chromosome 5 descends from a single founder
#' haplotype, which is what happens when a breed is founded from few animals. It
#' produces elevated \eqn{F_{ST}} and long-range haplotype homozygosity through
#' drift alone. It is included so that users can see what the method does when
#' differentiation is real but selection is absent, not because it is guaranteed
#' to escape detection.
#'
#' @examples
#' data(css_sim_truth)
#' css_sim_truth
"css_sim_truth"

#' Two-chromosome subset of the simulated panel
#'
#' Chromosomes 1 and 2 of [css_sim], carrying the incomplete and complete hard
#' sweeps respectively. Sized so that examples and tests run quickly.
#'
#' @format As [css_sim], restricted to two chromosomes.
#' @examples
#' data(css_sim_small)
#' table(css_sim_small$chr)
"css_sim_small"
