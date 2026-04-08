
test_that("calculate_overlap returns an overlapping value between 0 and 1", {
  mcp1 <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 5, ymax = 5), crs = 4326))
  mcp2 <- st_as_sfc(st_bbox(c(xmin = 4, ymin = 4, xmax = 9, ymax = 9), crs = 4326))
  result <- calculate_overlap(mcp1, mcp2)

  expect_is(result, "numeric")
  expect_true(result <= 1)
  expect_true(result >= 0)

})

test_that("calculate_overlap returns 0 if no overlapping", {
  mcp1 <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2), crs = 4326))
  mcp2 <- st_as_sfc(st_bbox(c(xmin = 3, ymin = 3, xmax = 5, ymax = 5), crs = 4326))

  result <- calculate_overlap(mcp1, mcp2)
  expect_true(result == 0)
})

