#' @title Territory outline from kernel densities
#' @description Polygon(s) holding a given share of the utilisation distribution (UD) of one individual or of a
#'   group of individuals (their UDs mixed, weighted by their numbers of samples): the 95% isopleth is the usual
#'   home range or territory, the 50% isopleth its core area.
#' @param ud A `wolfpack_ud` object from [space_use()].
#' @param individuals The individual(s) whose UDs are combined (e.g. the members of a pack).
#' @param level One or several levels between 0 and 1 (default: `c(0.95, 0.5)`).
#' @return An `sf` object (in the projected CRS of `ud`) with one (multi)polygon per level and its `area_km2`.
#' @export
#' @examples
#' data("samples_Carpathians")
#' ud <- space_use(samples_Carpathians)
#' terr <- isopleth(ud, c("CC00TC", "CC00X0", "CC04K4"))
#' terr
isopleth <- function(ud, individuals, level = c(0.95, 0.5)) {
  if (!inherits(ud, "wolfpack_ud")) stop("'ud' must come from space_use().")
  if (!any(individuals %in% ud$ids)) stop("None of 'individuals' is in 'ud'.")
  u <- .mix_ud(ud, individuals)

  # For each level: the cells of the isopleth, drawn as squares and merged into a polygon
  polys <- lapply(level, function(l) {
    cells <- which(.hr_cells(u, l))
    pts <- sf::st_as_sf(as.data.frame(ud$grid[cells, , drop = FALSE]), coords = c("x", "y"), crs = ud$crs)
    sf::st_union(sf::st_buffer(pts, ud$res / 2, endCapStyle = "SQUARE"))
  })
  out <- sf::st_sf(level = level, geometry = do.call(c, polys))
  out$area_km2 <- as.numeric(sf::st_area(out)) / 1e6
  out
}
