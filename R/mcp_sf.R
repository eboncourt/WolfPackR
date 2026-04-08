#' @title Calculate Minimum Convex Polygon (MCP) with Percentile Filtering
#' @description Computes the MCP for a set of spatial points using `sf`, with an option to exclude points beyond a given percentile of distance from the centroid. Largely inspired by juoe (https://rdrr.io/github/juoe/sdmflow/).
#' @param data An `sf` object or a data.frame with columns `X` and `Y`.
#' @param percentile The percentile of points (distance from centroid) excluded before calculating the MCP. Default is 100.
#' @param buffer_radius The buffer distance around the coordinates if there is insufficient data to calculate a MCP (less than 3 data points) (default: 0.01).
#' @return An `sf` object representing the MCP, or `NULL` if insufficient points.
#' @importFrom sf st_as_sf st_centroid st_distance st_union st_convex_hull
#' @export
#' @examples
#' obs <- data.frame(X = c(1,2,2,5,5,6,7,8,9), Y = c(1,2,9,4,5,6,1,8,9))
#' obs_sf <- st_as_sf(obs, coords = c("X", "Y"))
#' mcp_sf(obs_sf, percentile=100, buffer_radius=0.01)
mcp_sf <- function(data, percentile=100, buffer_radius=0.01){

  if (!inherits(data, "sf")) {
    if (!all(c("X", "Y") %in% colnames(data))) {
      stop("If 'data' is a data.frame, it must contain columns 'X' and 'Y'.")
    }
    data <- st_as_sf(data, coords = c("X", "Y"), crs = st_crs(data))
  }

  if (nrow(data) == 2) {
    warning(paste("Insufficient points (", nrow(data), ") to calculate MCP. Use of a buffer."))
    return(st_geometry(st_buffer(st_linestring(st_coordinates(data)), buffer_radius)))
  }

  if (nrow(data) == 1) {
    warning(paste("Insufficient points (", nrow(data), ") to calculate MCP. Use of a buffer."))
    return(st_geometry(st_buffer(data, buffer_radius)))
  }

  centroid <- st_centroid(st_union(data))
  dist <- as.numeric(st_distance(data, centroid))
  within_percentile_range <- dist <= quantile(dist, percentile/100)

  if (nrow(data[within_percentile_range,]) < 3) {
    warning(paste("Insufficient points (", length(data[within_percentile_range,]), ") after filtering by percentile to calculate MCP. Minimum required: 3"))
      if (length(data[within_percentile_range,]) == 2) {
      return(st_geometry(st_buffer(st_linestring(st_coordinates(data[within_percentile_range,])), buffer_radius)))
    }
  if (nrow(data[within_percentile_range,]) == 1) {
      return(st_geometry(st_buffer(data[within_percentile_range,], buffer_radius)))
    }
  }

  data_filter <- st_union(data[within_percentile_range,])

  st_convex_hull(data_filter)
}
