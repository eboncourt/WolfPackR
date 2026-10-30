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

  # For each level: the cells of the isopleth, drawn as squares and merged into a polygon.
  # Squares are drawn in grid units (whole numbers) so that neighbours share exactly the same edges, then put back
  # in metres: drawn in metres, rounding left hairline gaps between columns and the outline came in strips
  polys <- lapply(level, function(l) {
    cells <- which(.hr_cells(u, l))
    ij <- round(cbind((ud$grid[cells, 1] - ud$x0) / ud$res, (ud$grid[cells, 2] - ud$y0) / ud$res))
    sq <- sf::st_buffer(sf::st_as_sf(as.data.frame(ij), coords = 1:2), 0.5, endCapStyle = "SQUARE")
    sf::st_union(sq) * ud$res + c(ud$x0, ud$y0)
  })
  geom <- do.call(c, polys)
  sf::st_crs(geom) <- ud$crs
  out <- sf::st_sf(level = level, geometry = geom)
  out$area_km2 <- as.numeric(sf::st_area(out)) / 1e6
  out
}
