  # Load data
  data("relate")
  data("samples")

  # Prepare data
  gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
  obs <- merge(samples, gengroups, by = "Individual")

  # Execute the function
  result <- ugly_ducklings(obs, group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)

  # Tests with testthat
test_that("ugly_ducklings misclassified individuals", {
  expect_is(result, "data.frame")
  expect_true("new_group" %in% names(result))
  expect_equal(nrow(result), length(unique(samples$Individual)))
})

# Test with 1 data
obs2 <- obs[-c(2:5),]
result2 <- ugly_ducklings(obs2, group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)

test_that("ugly_ducklings can deal with individuals with only 1 point", {
  expect_is(result2, "data.frame")
  expect_true("new_group" %in% names(result2))
  expect_equal(nrow(result2), length(unique(samples$Individual)))
})

# Test with 2 data
obs3 <- obs[-c(3:5),]
result3 <- ugly_ducklings(obs3, group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)

test_that("ugly_ducklings can deal with individuals with only 2 points", {
  expect_is(result3, "data.frame")
  expect_true("new_group" %in% names(result3))
  expect_equal(nrow(result3), length(unique(samples$Individual)))
})
