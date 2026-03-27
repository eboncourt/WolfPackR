# Charger les packages nécessaires
library(testthat)
library(WolfPackR)

  # Créer un jeu de données fictif de relatedness
create_fake_relate <- function() { data.frame(
    ind1 = c("W1", "W1", "W1", "W2", "W2", "W3"),
    ind2 = c("W1", "W2", "W3", "W2", "W3", "W3"),
    indicator = c(1.0, 0.5, 0.2, 1.0, 0.3, 1.0)
  )
}

  # Test 1 : Cas nominal (succès)
  test_that("genetic_groups retourne un data.frame avec des groupes", {
    relate <- create_fake_relate()
    result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
    expect_is(result, "data.frame")
    expect_true("Individual" %in% names(result))
    expect_true("group" %in% names(result))
    expect_false(any(is.na(result$group)))
  })

  # Test 2 : Ajout des diagonales
  test_that("genetic_groups ajoute les diagonales à la matrice", {
    relate <- create_fake_relate()
    result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
    expect_equal(nrow(result), length(unique(c(relate$ind1, relate$ind2))))
  })

  # Test 3 : Individus sans groupe
  test_that("genetic_groups gère les individus sans groupe", {
    relate <- create_fake_relate()
    result <- genetic_groups(relate, estimator = "indicator", threshold = 0.5)
    w5_group <- result$group[result$Individual == "W5"]
    expect_true(length(unique(w5_group)) == 1)
    expect_true(w5_group != result$group[result$Individual != "W5"])
  })

  # Test 4 : Application du seuil de relatedness
  test_that("genetic_groups applique le seuil de relatedness", {
    relate <- create_fake_relate()
    result_low <- genetic_groups(relate, estimator = "indicator", threshold = 0.1)
    expect_true(length(unique(result_low$group)) < length(unique(result_low$Individual)))

    result_high <- genetic_groups(relate, estimator = "indicator", threshold = 0.8)
    expect_true(length(unique(result_high$group)) > 1)
  })
