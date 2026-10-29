#' @title Overlap between the space use of individuals
#' @description Pairwise overlap between the utilisation distributions (UD) of a [space_use()] object.
#'   All indices come from Fieberg & Kochanny (2005), the same ones as `kerneloverlap()` in adehabitatHR:
#'   * `"BA"` (default): Bhattacharyya's affinity, sum of sqrt(UD_i x UD_j). Symmetric, from 0 (no overlap) to 1 (same UD);
#'   * `"VI"`: volume of intersection, sum of min(UD_i, UD_j). Symmetric, from 0 to 1;
#'   * `"PHR"`: probability of finding j in the home range of i (isopleth `level` of i). Not symmetric;
#'   * `"HR"`: share of the home range of i covered by the home range of j. Not symmetric.
#' @param ud A `wolfpack_ud` object from [space_use()].
#' @param method One of `"BA"`, `"VI"`, `"PHR"`, `"HR"`.
#' @param individuals Individuals to compare (default: all).
#' @param level Isopleth defining the home range for `"PHR"` and `"HR"` (default: 0.95).
#' @return A square matrix (individuals x individuals). For PHR and HR, row i is the home range of i.
#' @references Fieberg, J. & Kochanny, C. O. (2005). Quantifying home-range overlap: the importance of the
#'   utilization distribution. Journal of Wildlife Management, 69, 1346-1359.
#' @export
#' @examples
#' data("samples_Carpathians")
#' ud <- space_use(samples_Carpathians)
#' ba <- ud_overlap(ud)
#' round(ba[1:4, 1:4], 2)
ud_overlap <- function(ud, method = c("BA", "VI", "PHR", "HR"), individuals = NULL, level = 0.95) {
  method <- match.arg(method)
  if (!inherits(ud, "wolfpack_ud")) stop("'ud' must come from space_use().")
  ids <- if (is.null(individuals)) ud$ids else intersect(individuals, ud$ids)
  U <- ud$U[, ids, drop = FALSE]

  # Home range of each individual = cells of its isopleth (for PHR and HR)
  in_hr <- function() vapply(ids, function(i) .hr_cells(U[, i], level), logical(nrow(U)))

  out <- switch(method,
    BA = crossprod(sqrt(U)),
    VI = vapply(ids, function(j) colSums(pmin(U, U[, j])), numeric(length(ids))),
    PHR = { H <- in_hr(); crossprod(H * 1, U) },                          # row i: UD of each j inside the HR of i
    HR = { H <- in_hr() * 1; crossprod(H) / colSums(H) }                   # row i: share of the HR of i shared with j
  )
  out <- matrix(out, length(ids), length(ids), dimnames = list(ids, ids))
  out
}

# Cells inside the isopleth `level` of a UD vector: the densest cells holding `level` of the volume
.hr_cells <- function(u, level) {
  o <- order(u, decreasing = TRUE)
  keep <- o[seq_len(min(length(o), sum(cumsum(u[o]) < level) + 1))]
  inside <- logical(length(u)); inside[keep] <- TRUE
  inside
}
