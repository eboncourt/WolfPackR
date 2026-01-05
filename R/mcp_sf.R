#' @title Calculate Minimum Convex Polygon (MCP) with Percentile Filtering
#' @description Computes the MCP for a set of spatial points using `sf`, with an option to exclude points beyond a given percentile of distance from the centroid. Largely inspired by juoe (https://rdrr.io/github/juoe/sdmflow/).
#' @param data An `sf` object or a data.frame with columns `X` and `Y`.
#' @param percentile The percentile of points (distance from centroid) excluded before calculating the MCP. Default is 100.
#' @return An `sf` object representing the MCP, or `NULL` if insufficient points.
#' @examples
#' obs_i <- subset(obs, group == 1)
#' mcp_i <- mcp_sf(obs_i, percentile = 95)
mcp_sf <- function(data, percentile=100){
  
  if (!inherits(data, "sf")) {
    if (!all(c("X", "Y") %in% colnames(data))) {
      stop("If 'data' is a data.frame, it must contain columns 'X' and 'Y'.")
    }
    data <- st_as_sf(data, coords = c("X", "Y"), crs = st_crs(data))
  }
  
  if (nrow(data) < 3) {
    warning(paste("Insufficient points (", nrow(data), ") to calculate MCP. Minimum required: 3"))
    return(NULL)
  }
  
  centroid <- st_centroid(st_union(data))
  dist <- as.numeric(st_distance(data, centroid))
  within_percentile_range <- dist <= quantile(dist, percentile/100)
  
  if (length(data[within_percentile_range,]) < 3) {
    warning(paste("Insufficient points (", length(data[within_percentile_range,]), ") after filtering by percentile to calculate MCP. Minimum required: 3"))
    return(NULL)
  }
  
  data_filter <- st_union(data[within_percentile_range,])
  
  st_convex_hull(data_filter)
}