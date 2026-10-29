test_that("isopleth gives one polygon per level, the core inside the territory", {
  ud <- space_use(sim_samples(), h = 2000)
  iso <- isopleth(ud, c("A1", "A2", "A3"), level = c(0.95, 0.5))
  expect_s3_class(iso, "sf")
  expect_equal(iso$level, c(0.95, 0.5))
  expect_gt(iso$area_km2[1], iso$area_km2[2])
  expect_true(sf::st_covered_by(iso[2, ], sf::st_buffer(iso[1, ], 1), sparse = FALSE)[1, 1])
  expect_error(isopleth(ud, "nobody"), "None")
})
