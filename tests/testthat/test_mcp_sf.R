library(testthat)
library(WolfPackR)
library(sf)

# load data
data("samples")

# Prepare data
gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
obs <- merge(samples, gengroups, by = "Individual")

# Tests
test_that("mcp_sf retourne un objet sf avec des polygones MCP", {

  result <- mcp_sf(samples)
  expect_is(result, "sfc_POLYGON")

})

