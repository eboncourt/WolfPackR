testthat::test_that("ugly_ducklings identifie les individus mal assignés", {
  # Créer un jeu de données fictif avec des coordonnées et des groupes
  samples <- data.frame(
    Individual = c("W1", "W2", "W3", "W4"),
    X = c(0, 1, 10, 11),
    Y = c(0, 1, 10, 11),
    group = c(1, 1, 2, 2)
  )

  # Exécuter la fonction
  result <- WolfPackR::ugly_ducklings(samples, group = "group", min_overlap = 0.5, buffer = 0, mcp.percent = 100)

  # Vérifications
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("new_group" %in% names(result))
  testthat::expect_equal(nrow(result), nrow(samples))
})
