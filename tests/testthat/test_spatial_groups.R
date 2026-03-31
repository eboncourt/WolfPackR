testthat::test_that("spatial_groups groups individuals spatially", {

  # Load data
  data("samples")

  # Execute the function
  result <- WolfPackR::spatial_groups(
    samples, percentile = 100,
    buffer_radius = 0.00001, max_iterations = 20, min_mcp_overlap = 0.2
  )

  # Tests
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("Subgroup" %in% names(result))
  testthat::expect_true(all(!is.na(result$Subgroup)))
})
