# Create a sample dataset of relatedness
create_fake_relate <- function() { data.frame(
    ind1 = c("W1", "W1", "W1", "W2", "W2", "W3"),
    ind2 = c("W1", "W2", "W3", "W2", "W3", "W3"),
    indicator = c(1.0, 0.5, 0.2, 1.0, 0.3, 1.0)
  )
}

# Test 1: Nominal case (success)
test_that("genetic_groups returns a data.frame with groups", {
  relate <- create_fake_relate()
  result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Individual", "group") %in% names(result)))
  expect_false(any(is.na(result$group)))
})

# Test 2: one row per individual
test_that("genetic_groups returns one row per individual", {
  relate <- create_fake_relate()
  result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
  expect_equal(nrow(result), length(unique(c(relate$ind1, relate$ind2))))
})

# Test 3: Individuals without a group
test_that("genetic_groups gives its own group to an unrelated individual", {
  relate <- create_fake_relate()
  result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
  w3_group <- result$group[result$Individual == "W3"]
  expect_length(w3_group, 1)
  expect_false(w3_group %in% result$group[result$Individual != "W3"])
  expect_equal(result$group[result$Individual == "W1"], result$group[result$Individual == "W2"])
})

# Test 4: Applying the relatedness threshold
test_that("genetic_groups applies the relatedness threshold", {
  relate <- create_fake_relate()
  result_low <- genetic_groups(relate, estimator = "indicator", threshold = 0.1)
  expect_equal(length(unique(result_low$group)), 1)   # everybody linked

  result_high <- genetic_groups(relate, estimator = "indicator", threshold = 0.8)
  expect_equal(length(unique(result_high$group)), 3)  # nobody linked
})

# Test 5: regression test, pairs used to be dropped when ind1 sorted after ind2
test_that("genetic_groups does not depend on the order of ind1 and ind2", {
  rel_ab <- data.frame(ind1 = c("A", "C"), ind2 = c("B", "D"), r = c(0.9, 0.9))
  rel_ba <- data.frame(ind1 = c("B", "D"), ind2 = c("A", "C"), r = c(0.9, 0.9))
  res_ab <- genetic_groups(rel_ab, "r", 0.4)
  res_ba <- genetic_groups(rel_ba, "r", 0.4)

  expect_equal(length(unique(res_ab$group)), 2)
  expect_equal(length(unique(res_ba$group)), 2)
  expect_equal(res_ba$group[res_ba$Individual == "A"], res_ba$group[res_ba$Individual == "B"])
  expect_equal(res_ba$group[res_ba$Individual == "C"], res_ba$group[res_ba$Individual == "D"])
})

# Test 6: one row per pair in random order (like coancestry output) gives the same groups as the full matrix
test_that("genetic_groups gives the same groups with one row per pair as with the full matrix", {
  relate <- pkg_data("relate")
  full <- genetic_groups(relate, "wang", 0.4)

  half <- relate[as.character(relate$ind1) < as.character(relate$ind2), ]  # one row per pair
  set.seed(1)
  flip <- stats::runif(nrow(half)) < 0.5                                  # then shuffle the orientation
  half[flip, c("ind1", "ind2")] <- half[flip, c("ind2", "ind1")]
  from_half <- genetic_groups(half, "wang", 0.4)

  # same partition (group numbers may differ, so we compare the co-membership)
  same_group <- function(res) outer(res$group, res$group, "==")
  full <- full[order(full$Individual), ]
  from_half <- from_half[order(from_half$Individual), ]
  expect_equal(same_group(full), same_group(from_half))
})
