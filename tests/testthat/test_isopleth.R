test_that("isopleth gives one polygon per level, the core inside the territory", {
  ud <- space_use(sim_samples(), h = 2000)
  iso <- isopleth(ud, c("A1", "A2", "A3"), level = c(0.95, 0.5))
  expect_s3_class(iso, "sf")
  expect_equal(iso$level, c(0.95, 0.5))
  expect_gt(iso$area_km2[1], iso$area_km2[2])
  expect_true(sf::st_covered_by(iso[2, ], sf::st_buffer(iso[1, ], 1), sparse = FALSE)[1, 1])
  expect_error(isopleth(ud, "nobody"), "None")
})

test_that("isopleth outlines are in one piece whatever the bandwidth (no gaps between grid columns)", {
  ud <- space_use(sim_samples(), h = 1783.29)                     # a bandwidth that left hairline gaps in 1.1.0
  iso <- isopleth(ud, c("A1", "A2", "A3"), level = 0.95)
  n_parts <- length(sf::st_cast(sf::st_geometry(iso), "POLYGON"))
  n_closed <- length(sf::st_cast(sf::st_union(sf::st_buffer(sf::st_geometry(iso), 1)), "POLYGON"))
  expect_equal(n_parts, n_closed)                                 # a 1 m buffer doesn't merge anything
  expect_equal(sf::st_crs(iso), ud$crs)
})
