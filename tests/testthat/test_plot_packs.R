test_that("plot_packs retourne un objet ggplot", {

  # Load data
  data("relate")
  data("samples")

  # Prepare data
  gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
  obs <- merge(samples, gengroups, by = "Individual")

  # Execute the function

  plot_result <- plot_packs(obs, group)

  # Tests with testthat
  expect_is(plot_result, "leaflet")  # The result is a leaflet object
})
