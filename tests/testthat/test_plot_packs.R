  # Load data
  data("relate")
  data("samples")

  # Prepare data
  gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
  obs <- merge(samples, gengroups, by = "Individual")

  # Execute the function

  plot_result <- plot_packs(obs, group)

# Test for output format
test_that("plot_packs returns a leaflet object", {

    # Tests with testthat
  expect_is(plot_result, "leaflet")  # The result is a leaflet object
})

# test for empty obs
obs2 <- obs[c(),]

test_that("plot_packs returns an error if the input is empty", {
  expect_error(plot_packs(obs2, group))  # expect an error message
})

# test for input of different class
obs3 <- data.frame()

test_that("plot_packs returns an error if the input is not an sf object", {
  expect_error(plot_packs(obs3, group))  # expect an error message
})

# test for lone individuals
obs
obs[obs$group == 1,]$group <- "Lone Individual"
plot_result <- plot_packs(obs, group)

test_that("plot_packs can deal with lone individuals", {

  # Tests with testthat
  expect_is(plot_result, "leaflet")  # The result is a leaflet object
})
