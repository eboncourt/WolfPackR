#' Calculate Overlap Between Two Polygons
#'
#' This function calculates the overlap ratio between two polygons.
#' It handles edge cases such as NULL inputs, empty intersections, and non-polygon geometries.
#'
#' @param mcp1 The first polygon (sf object).
#' @param mcp2 The second polygon (sf object).
#'
#' @return A numeric value between 0 and 1 representing the overlap ratio.
#'
#' @examples
#' # Example usage:
#' overlap_ratio <- calculate_overlap(mcp1 = polygon1, mcp2 = polygon2)
#'
#' @export
calculate_overlap <- function(mcp1, mcp2) {
  # Check inputs
  if (is.null(mcp1) || is.null(mcp2)) {
    return(0)
  }

  # Calculate intersection
  intersection <- tryCatch({
    st_intersection(mcp1, mcp2)
  }, error = function(e) {
    NULL
  })

  # Check intersection
  if (is.null(intersection)) {
    return(0)
  }

  # Check if intersection is empty
  if (inherits(intersection, "sf") && nrow(intersection) == 0) {
    return(0)
  }

  # Check if intersection is a GEOMETRYCOLLECTION or MULTILINESTRING
  if (inherits(intersection, "sf") && st_geometry_type(intersection) %in% c("GEOMETRYCOLLECTION", "MULTILINESTRING")) {
    intersection <- st_buffer(intersection, dist = 0.001)
  }

  # Calculate areas
  intersection_area <- tryCatch({
    area_value <- st_area(intersection)
    if (is.null(area_value) || length(area_value) == 0) {
      NA_real_
    } else {
      as.numeric(area_value)
    }
  }, error = function(e) {
    NA_real_
  })

  mcp2_area <- tryCatch({
    area_value <- st_area(mcp2)
    if (is.null(area_value) || length(area_value) == 0) {
      NA_real_
    } else {
      as.numeric(area_value)
    }
  }, error = function(e) {
    NA_real_
  })

  # Handle NA values and vectors
  if (length(intersection_area) > 1) {
    intersection_area <- sum(intersection_area, na.rm = TRUE)
  }
  if (length(mcp2_area) > 1) {
    mcp2_area <- sum(mcp2_area, na.rm = TRUE)
  }

  if (is.na(intersection_area) || is.na(mcp2_area)) {
    return(0)
  }

  # Check area validity
  if (mcp2_area <= 0) {
    return(0)
  }

  # Calculate overlap ratio
  overlap_ratio <- intersection_area / mcp2_area

  # Final check
  if (is.na(overlap_ratio) || !is.finite(overlap_ratio)) {
    return(0)
  }

  # Ensure ratio is between 0 and 1
  if (overlap_ratio < 0) {
    return(0)
  } else if (overlap_ratio > 1) {
    return(1)
  } else {
    return(overlap_ratio)
  }
}
