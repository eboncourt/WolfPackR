square <- function(xmin, ymin, xmax, ymax, crs = 4326) {
  sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = crs))
}

test_that("calculate_overlap returns an overlapping value between 0 and 1", {
  result <- calculate_overlap(square(0, 0, 5, 5), square(4, 4, 9, 9))
  expect_type(result, "double")
  expect_gt(result, 0)
  expect_lt(result, 1)
})

test_that("calculate_overlap is relative to the area of mcp2", {
  big <- square(0, 0, 4, 4, crs = NA)
  small <- square(1, 1, 2, 2, crs = NA)
  expect_equal(calculate_overlap(big, small), 1)       # small entirely inside big
  expect_equal(calculate_overlap(small, big), 1 / 16)  # 1 unit2 out of 16
})

test_that("calculate_overlap returns 0 if no overlapping", {
  expect_equal(calculate_overlap(square(0, 0, 2, 2), square(3, 3, 5, 5)), 0)
})

test_that("calculate_overlap returns 0 if one polygon is null", {
  expect_equal(calculate_overlap(square(0, 0, 2, 2), NULL), 0)
})

# Regression test: different CRS used to be swallowed by a tryCatch and return 0
test_that("calculate_overlap stops when the CRS differ", {
  expect_error(calculate_overlap(square(0, 0, 2, 2), square(1, 1, 3, 3, crs = NA)), "same CRS")
})
