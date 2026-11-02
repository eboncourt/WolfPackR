#' @title Territories of several packs, with or without overlap
#' @description Kernel isopleths of several packs at once, each one from the UDs of its members mixed (weighted by
#'   their numbers of samples, as in [isopleth()]).
#'
#'   The isopleths of neighbouring packs usually overlap, even when their samples don't mix: the kernel spreads every
#'   sample over a distance of about h, so each territory spills over its edge. With `exclusive = TRUE`, a cell claimed
#'   by several packs (inside their largest isopleth) goes to the pack whose UD is the highest there: the boundary runs
#'   where both packs use the space as much, the rule [attribute_samples()] uses for the samples. The smaller isopleths
#'   (core areas) are cut the same way, so a core stays inside its own territory.
#'
#'   The shares are made on the UDs: they tell where each pack is more likely to be, not where wolves defend a border.
#'   Contested cells are best reported too, e.g. with `shared_pct` of [summarise_packs()], measured on the plain isopleths.
#' @param ud A `wolfpack_ud` object from [space_use()].
#' @param packs Named list of the members of each pack (names = pack labels), e.g. `split(ind$Individual, ind$Pack)`.
#' @param levels Isopleth levels (default: 95% territory and 50% core area).
#' @param exclusive Share out the cells claimed by several packs (default: `FALSE`, the plain isopleths of [isopleth()]).
#' @return An `sf` object (in the projected CRS of `ud`) with one row per pack and level: `Pack`, `level` and
#'   `area_km2`. With `exclusive = TRUE`, a pack that loses every cell of a level has no row for that level.
#' @seealso [isopleth()] for a single set of individuals.
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25, period = "Campaign")
#' members <- split(res$individuals$Individual, res$individuals$Pack)
#' plain <- pack_territories(res$ud, members)
#' excl <- pack_territories(res$ud, members, exclusive = TRUE)
#' merge(sf::st_drop_geometry(plain), sf::st_drop_geometry(excl), by = c("Pack", "level"),
#'       suffixes = c("_plain", "_exclusive"))
pack_territories <- function(ud, packs, levels = c(0.95, 0.5), exclusive = FALSE) {
  if (!inherits(ud, "wolfpack_ud")) stop("'ud' must come from space_use().")
  if (!is.list(packs) || !length(packs) || is.null(names(packs)) || any(is.na(names(packs)) | names(packs) == ""))
    stop("'packs' must be a named list: the members of each pack, named after the pack.")
  packs <- lapply(packs, function(m) intersect(as.character(m), ud$ids))
  if (any(lengths(packs) == 0)) stop("No member of ", paste(names(packs)[lengths(packs) == 0], collapse = ", "), " is in 'ud'.")

  # UD of each pack (one column per pack) and, if exclusive, the owner of each cell: the pack with the highest UD
  # among the packs whose largest isopleth holds the cell
  U <- vapply(packs, function(m) .mix_ud(ud, m), numeric(nrow(ud$grid)))
  owner <- NULL
  if (exclusive && length(packs) > 1) {
    claim <- vapply(seq_along(packs), function(k) .hr_cells(U[, k], max(levels)), logical(nrow(U)))
    D <- U; D[!claim] <- -1
    owner <- max.col(D, ties.method = "first")
  }

  # One polygon per pack and level
  rows <- list()
  for (k in seq_along(packs)) for (l in levels) {
    cells <- .hr_cells(U[, k], l)
    if (!is.null(owner)) cells <- cells & owner == k
    if (any(cells)) rows[[length(rows) + 1]] <- sf::st_sf(Pack = names(packs)[k], level = l,
                                                          geometry = .cells_polygon(ud, which(cells)))
  }
  out <- do.call(rbind, rows)
  out$area_km2 <- as.numeric(sf::st_area(out)) / 1e6
  rownames(out) <- NULL
  out
}

# Grid cells drawn as squares and merged into a (multi)polygon, in the CRS of ud.
# Squares are drawn in grid units (whole numbers) so that neighbours share exactly the same edges, then put back
# in metres: drawn in metres, rounding left hairline gaps between columns and the outline came in strips
.cells_polygon <- function(ud, cells) {
  ij <- round(cbind((ud$grid[cells, 1] - ud$x0) / ud$res, (ud$grid[cells, 2] - ud$y0) / ud$res))
  sq <- sf::st_buffer(sf::st_as_sf(as.data.frame(ij), coords = 1:2), 0.5, endCapStyle = "SQUARE")
  g <- sf::st_union(sq) * ud$res + c(ud$x0, ud$y0)
  sf::st_crs(g) <- ud$crs
  g
}
