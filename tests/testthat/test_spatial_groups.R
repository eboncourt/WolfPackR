# Load data
data("samples")

# base test
test_that("spatial_groups groups individuals spatially", {


  # Execute the function
  result <- spatial_groups(
    samples, percentile = 100,
    buffer_radius = 0.00001, max_iterations = 20, min_mcp_overlap = 0.2
  )

  # Tests with testthat
  expect_is(result, "data.frame")
  expect_true("Subgroup" %in% names(result))
  expect_true(all(!is.na(result$Subgroup)))
})

# test with groups
samples2 <- samples %>% mutate(group= 1:47)

result2 <- spatial_groups(
  samples2, percentile = 100,
  group = "group",
  buffer_radius = 0.00001, max_iterations = 20, min_mcp_overlap = 0.2
)

test_that("spatial_groups can deal with given groups", {
  expect_is(result2, "data.frame")
  expect_true("Subgroup" %in% names(result2))
  expect_true(all(!is.na(result2$Subgroup)))
})
