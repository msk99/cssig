#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @importFrom stats qnorm pnorm quantile p.adjust cor complete.cases setNames
#'   isoreg approx sd median
#' @importFrom utils head tail globalVariables packageVersion
#' @importFrom graphics text
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_vline
#'   geom_rect geom_segment geom_text geom_density geom_histogram geom_tile
#'   geom_abline scale_x_continuous scale_y_continuous scale_colour_manual
#'   scale_fill_manual scale_fill_gradient2 labs facet_wrap facet_grid theme
#'   theme_minimal element_blank element_text element_line element_rect
#'   expansion coord_cartesian guides guide_legend margin unit waiver
NULL

# Internal column names created by the package. Declared so that R CMD check
# does not report them as undefined globals when used in data.table's
# non-standard evaluation. See PLAN.md section 3.7.
utils::globalVariables(c(
  ".", ".N", ".I", ".SD", ".GRP",
  "chr", "pos", "snp", "css", "css_smooth", "zbar", "p", "p_adj", "qval",
  "n_window", "significant", "significant2", "run", "grp", "gap",
  "chr_cum", "pos_cum", "peak_css", "peak_pos", "peak_snp", "n_snps",
  "start", "end", "start_padded", "end_padded", "region_id", "direction",
  "css_pos", "css_neg", "css_signed", "cohort", "test", "value", "rank_i",
  "smooth_value", "in_region", "x", "y", "xend", "yend", "label", "mid",
  "chr_index", "shade", "n_sig", "prop", "test1", "test2", "r",
  "css_raw", "threshold", "is_top", "is_top2", "keep", "..cols",
  "offset", "len", "lo", "hi", "xmin", "xmax", "up", "down",
  "score", "sig", "sig2", "i.N", "n_tests", ".idx", ".row",
  "daf_selected", "daf_reference", "ancestral", "selected", "near_signal",
  "ymin", "ymax", "ybase"
))
