# load data
samples <- pkg_data("samples")

test_that("mcp_sf returns a polygon for a regular set of points", {
  result <- mcp_sf(samples)
  expect_s3_class(result, "sfc_POLYGON")
  expect_equal(sf::st_crs(result), sf::st_crs(samples))
})

test_that("mcp_sf returns a polygon even with 1 point", {
  expect_warning(result <- mcp_sf(samples[1, ]), "Insufficient points")
  expect_s3_class(result, "sfc_POLYGON")
})

test_that("mcp_sf returns a polygon even with 2 points", {
  expect_warning(result <- mcp_sf(samples[1:2, ]), "Insufficient points")
  expect_s3_class(result, "sfc_POLYGON")
})

test_that("mcp_sf accepts a data.frame with X and Y", {
  df <- data.frame(X = c(1, 2, 2, 5, 5, 6, 7, 8, 9), Y = c(1, 2, 9, 4, 5, 6, 1, 8, 9))
  expect_s3_class(mcp_sf(df), "sfc_POLYGON")
  expect_error(mcp_sf(data.frame(a = 1)), "X' and 'Y'")
})

test_that("mcp_sf refuses an empty input", {
  expect_error(mcp_sf(samples[0, ]), "no points")
})

# Regression test: with 1 or 2 points left after the percentile filter, the old code used
# length() (= number of columns) instead of nrow() and returned a LINESTRING with no area
test_that("mcp_sf returns a polygon when the percentile filter leaves fewer than 3 points", {
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(id = 1:10, a = 1, b = 2, X = stats::runif(10), Y = stats::runif(10)),
                      coords = c("X", "Y"))
  expect_warning(result <- mcp_sf(pts, percentile = 15), "Insufficient points \\( 2 \\)")
  expect_s3_class(result, "sfc_POLYGON")
  expect_gt(as.numeric(sf::st_area(result)), 0)
})

# Regression test: the 2-point buffer used to lose the CRS and be computed in degrees,
# while the 1-point buffer was in metres
test_that("mcp_sf keeps the CRS and uses metres for 1 and 2 points in lon/lat", {
  two <- samples[1:2, ]
  r <- 100                                                      # metres
  m1 <- suppressWarnings(mcp_sf(two[1, ], buffer_radius = r))
  m2 <- suppressWarnings(mcp_sf(two, buffer_radius = r))
  expect_equal(sf::st_crs(m2), sf::st_crs(samples))

  # 1 point: a disc of radius r; 2 points: a "sausage" of length L and width 2r.
  # s2 buffers are approximations, hence the loose tolerance (the old bug was off by ~10^7)
  L <- as.numeric(sf::st_distance(two[1, ], two[2, ]))
  expect_equal(as.numeric(sf::st_area(m1)), pi * r^2, tolerance = 0.25)
  expect_equal(as.numeric(sf::st_area(m2)), L * 2 * r + pi * r^2, tolerance = 0.25)
})

test_that("mcp_sf buffers aligned or duplicated points instead of returning a line", {
  aligned <- sf::st_as_sf(data.frame(X = c(0, 1, 2, 3), Y = c(0, 1, 2, 3)), coords = c("X", "Y"))
  expect_warning(result <- mcp_sf(aligned, buffer_radius = 0.1), "aligned or duplicated")
  expect_s3_class(result, "sfc_POLYGON")
  expect_gt(as.numeric(sf::st_area(result)), 0)
})

# The 1- and 2-point MCPs must be usable together (they used to have different CRS,
# so calculate_overlap() failed silently and returned 0)
test_that("MCPs of individuals with 1 and 2 points can be compared", {
  two <- samples[1:2, ]
  m1 <- suppressWarnings(mcp_sf(two[1, ], buffer_radius = 100))
  m2 <- suppressWarnings(mcp_sf(two, buffer_radius = 100))
  expect_gt(calculate_overlap(m2, m1), 0.9)  # the disc sits at the end of the sausage
})
