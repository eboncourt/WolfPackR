test_that("ugly_ducklings misclassified individuals", {

  # Load data
  data("relate")
  data("samples")

  # Prepare data
  gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
  obs <- merge(samples, gengroups, by = "Individual")

  # Execute the function
  result <- ugly_ducklings(obs, group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)

  # Tests with testthat
  expect_is(result, "data.frame")
  expect_true("new_group" %in% names(result))
  expect_equal(nrow(result), length(unique(samples$Individual)))
})
