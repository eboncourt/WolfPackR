testthat::test_that("plot_packs retourne un objet ggplot", {
  # Créer un jeu de données fictif
  samples <- data.frame(
    Individual = c("W1", "W2", "W3", "W4"),
    X = c(0, 1, 10, 11),
    Y = c(0, 1, 10, 11),
    Subgroup = c(1, 1, 2, 2)
  )

  # Exécuter la fonction (et capturer le résultat)
  plot_result <- WolfPackR::plot_packs(samples, Subgroup)

  # Vérifications
  testthat::expect_is(plot_result, "ggplot")  # Le résultat est un objet ggplot
})
