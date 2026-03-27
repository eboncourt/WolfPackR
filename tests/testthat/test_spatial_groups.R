library(sf)

testthat::test_that("spatial_groups regroupe les individus spatialement", {
  # Créer un jeu de données fictif avec des coordonnées et des groupes
  samples <- data.frame(
    Individual = c("W1", "W2", "W3", "W4"),
    X = c(0, 1, 10, 11),
    Y = c(0, 1, 10, 11),
    new_group = c(1, 1, 2, 2)
  )

  # Exécuter la fonction
  result <- WolfPackR::spatial_groups(
    samples, group = "new_group", percentile = 100,
    buffer_radius = 1, max_iterations = 20, min_mcp_overlap = 0.5
  )

  # Vérifications
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("Subgroup" %in% names(result))
  testthat::expect_true(all(!is.na(result$Subgroup)))
})
