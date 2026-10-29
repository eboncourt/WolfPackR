#' @title Which territory does each sample fall in?
#' @description Every sample goes to the territory whose UD is the densest where the sample lies, or to
#'   `"outside"` when no territory beats a uniform use of the study area (the grid of `ud`).
#'   A territory's UD is the mix of its members' UDs, weighted by their numbers of samples. With
#'   `leave_one_out = TRUE`, a member's own samples are judged against the territory without that member,
#'   so that an individual doesn't make its own territory look right.
#' @param ud A `wolfpack_ud` object from [space_use()].
#' @param territories Either a named list of member IDs (one element per territory), or a data.frame with
#'   `Individual` and a second column giving the territory (`NA` = no territory).
#' @param leave_one_out See description (default: TRUE).
#' @param outside A sample is "outside" when its best territory density is below `outside` times the uniform
#'   density (default: 1).
#' @return The samples of `ud` (data.frame with their columns, projected `x`, `y`, `row` in the original data)
#'   plus `Territory`, the territory each sample is attributed to.
#' @export
#' @examples
#' data("samples_Carpathians")
#' ud <- space_use(samples_Carpathians)
#' terr <- list(c = c("CC00TC", "CC00X0", "CC04K4"), d = c("CC025E", "CC024X", "CC03E3"))
#' att <- attribute_samples(ud, terr)
#' table(att$Territory)
attribute_samples <- function(ud, territories, leave_one_out = TRUE, outside = 1) {
  if (!inherits(ud, "wolfpack_ud")) stop("'ud' must come from space_use().")
  if (is.data.frame(territories)) {
    lab <- as.character(territories[[2]])
    territories <- split(as.character(territories$Individual)[!is.na(lab)], lab[!is.na(lab)])
  }
  territories <- lapply(territories, function(m) intersect(as.character(m), ud$ids))
  territories <- territories[lengths(territories) > 0]
  s <- ud$samples
  if (!length(territories)) { s$Territory <- "outside"; return(s) }

  # Each territory as a weighted sum of its members' UDs ("mass"), then densities at each sample's cell
  W <- vapply(territories, function(m) ifelse(ud$ids %in% m, ud$n, 0), numeric(length(ud$ids)))
  W <- matrix(W, nrow = length(ud$ids), dimnames = list(ud$ids, names(territories)))
  M <- ud$U %*% W                                                     # cells x territories
  tot <- colSums(W)
  D <- M[s$cell, , drop = FALSE]
  if (leave_one_out) {                                               # remove the sample's own individual from its territories
    own <- W[s$Individual, , drop = FALSE] > 0
    ni <- ud$n[s$Individual]
    ui <- ud$U[cbind(s$cell, match(s$Individual, ud$ids))]
    D <- D - own * (ni * ui)
    den <- matrix(tot, nrow(D), ncol(D), byrow = TRUE) - own * ni
    D <- ifelse(den > 0, D / den, 0)
  } else {
    D <- sweep(D, 2, tot, "/")
  }
  D <- matrix(D, nrow = nrow(s), dimnames = list(NULL, names(territories)))

  best <- max.col(D, ties.method = "first")
  s$Territory <- ifelse(D[cbind(seq_len(nrow(D)), best)] < outside / nrow(ud$grid), "outside", colnames(D)[best])
  s
}
