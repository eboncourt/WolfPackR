#' @title Calculate Minimum Convex Polygon (MCP) with Percentile Filtering
#' @description Computes the MCP for a set of spatial points using `sf`, keeping only the points within a given percentile of distance from the centroid. Largely inspired by juoe (https://rdrr.io/github/juoe/sdmflow/).
#' @param data An `sf` object or a data.frame with columns `X` and `Y`.
#' @param percentile The percentage of points (the closest to the centroid) kept to calculate the MCP. Default is 100 (all points).
#' @param buffer_radius The buffer distance used around the point(s) when there are fewer than 3 points left, or when all points are aligned. It is in metres for geographic (lon/lat) data and in map units otherwise (default: 0.01).
#' @return An `sfc` polygon with the same CRS as `data`: the MCP, or a buffer around the points when no MCP can be computed.
#' @importFrom sf st_as_sf st_centroid st_distance st_union st_convex_hull st_buffer st_geometry st_cast st_combine st_geometry_type
#' @importFrom stats quantile
#' @export
#' @examples
#' obs <- data.frame(X = c(1,2,2,5,5,6,7,8,9), Y = c(1,2,9,4,5,6,1,8,9))
#' obs_sf <- sf::st_as_sf(obs, coords = c("X", "Y"))
#' mcp_sf(obs_sf, percentile=100, buffer_radius=0.01)
mcp_sf <- function(data, percentile=100, buffer_radius=0.01){

  if (!inherits(data, "sf")) {
    if (!all(c("X", "Y") %in% colnames(data))) {
      stop("If 'data' is a data.frame, it must contain columns 'X' and 'Y'.")
    }
    data <- st_as_sf(data, coords = c("X", "Y"))  # a plain data.frame has no CRS, so none is set
  }

  if (nrow(data) == 0) {
    stop("'data' contains no points.")
  }

  # Percentile filter: keep the points closest to the centroid (only meaningful with 3+ points)
  if (nrow(data) >= 3) {
    centroid <- st_centroid(st_union(data))
    dist <- as.numeric(st_distance(data, centroid))
    data <- data[dist <= quantile(dist, percentile / 100), ]
  }
  n <- nrow(data)  # nrow, not length: length() of an sf object is its number of columns

  # Fewer than 3 points: no MCP possible, so we use a buffer around the point or the segment.
  # Everything stays in the CRS of the data, so the buffer has the same units for 1 and 2 points
  if (n < 3) {
    warning(paste("Insufficient points (", n, ") to calculate MCP. Use of a buffer."))
    geom <- st_geometry(data)
    if (n == 2) geom <- st_cast(st_combine(geom), "LINESTRING")  # segment between the 2 points
    return(st_buffer(geom, buffer_radius))
  }

  # Regular case: convex hull of the remaining points
  hull <- st_convex_hull(st_union(data))

  # Aligned or duplicated points give a line or a point with zero area: buffer it too
  if (!all(st_geometry_type(hull) %in% c("POLYGON", "MULTIPOLYGON"))) {
    warning("Points are aligned or duplicated: MCP has no area. Use of a buffer.")
    hull <- st_buffer(hull, buffer_radius)
  }

  hull
}
