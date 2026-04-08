# load data
data("samples")
data("relate")

# Prepare data
gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
obs <- merge(samples, gengroups, by = "Individual")

# Test for sf object
test_that("mcp_sf returns an sf object with polygons", {

  result <- mcp_sf(samples)
  expect_is(result, "sfc_POLYGON")

})

# Test for small datasets
samples1 <- samples[c(1),]
test_that("mcp_sf returns an sf object with polygons even with 1 data", {

  result1 <- mcp_sf(samples1)
  expect_is(result, "sfc_POLYGON")

})

samples2 <- samples[c(1),]
test_that("mcp_sf returns an sf object with polygons even with 2 data", {

  result2 <- mcp_sf(samples2)
  expect_is(result, "sfc_POLYGON")

})

# Test for data.frame
samples3 <- data.frame(X = c(1,2,3,4,5,6,7,8,9), Y = c(1,2,3,4,5,6,7,8,9))
data <- st_as_sf(samples3, coords = c("X", "Y"), crs = st_crs(samples3))
test_that("mcp_sf returns an sf object with polygons even with 2 data", {

  result3 <- mcp_sf(samples3)
  expect_is(result3, "sfc_POLYGON")

})

