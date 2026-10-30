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
  half <- pkg_data("relate_example")                                      # one row per pair, as relatedness() gives
  swap <- half
  swap[, c("ind1", "ind2")] <- half[, c("ind2", "ind1")]
  full <- genetic_groups(rbind(half, swap), "QG", 0.25)                   # both orientations of every pair

  set.seed(1)
  flip <- stats::runif(nrow(half)) < 0.5                                  # shuffle the orientation
  half[flip, c("ind1", "ind2")] <- half[flip, c("ind2", "ind1")]
  from_half <- genetic_groups(half, "QG", 0.25)

  # same partition (group numbers may differ, so we compare the co-membership)
  same_group <- function(res) outer(res$group, res$group, "==")
  full <- full[order(full$Individual), ]
  from_half <- from_half[order(from_half$Individual), ]
  expect_equal(same_group(full), same_group(from_half))
})

# Two families of 4 joined by a single weak link: the classic case where connected components merge them
two_families <- function() {
  fam <- function(ids) { p <- t(utils::combn(ids, 2)); data.frame(ind1 = p[, 1], ind2 = p[, 2], r = 0.5) }
  rbind(fam(c("A1", "A2", "A3", "A4")), fam(c("B1", "B2", "B3", "B4")),
        data.frame(ind1 = "A1", ind2 = "B1", r = 0.3),          # one spurious link between the families
        data.frame(ind1 = "Z", ind2 = "A2", r = 0.05))          # and someone unrelated to everybody
}

test_that("community detection splits two families that components merge", {
  rel <- two_families()
  comp <- genetic_groups(rel, "r", 0.25, method = "components")
  expect_equal(length(unique(comp$group[comp$Individual != "Z"])), 1)
  for (m in c("leiden", "louvain", "walktrap")) {
    res <- genetic_groups(rel, "r", 0.25, method = m)
    expect_equal(nrow(res), 9)
    g <- stats::setNames(res$group, res$Individual)
    expect_equal(length(unique(g[c("A1", "A2", "A3", "A4")])), 1)
    expect_equal(length(unique(g[c("B1", "B2", "B3", "B4")])), 1)
    expect_false(g["A1"] == g["B1"])
    expect_false(g["Z"] %in% g[names(g) != "Z"])           # the unrelated one stays alone
  }
})

test_that("genetic_groups is reproducible and leaves the user's random numbers alone", {
  rel <- two_families()
  set.seed(123); before <- stats::runif(1)
  set.seed(123)
  a <- genetic_groups(rel, "r", 0.25, method = "louvain")
  after <- stats::runif(1)
  expect_equal(before, after)
  expect_equal(a, genetic_groups(rel, "r", 0.25, method = "louvain"))
})

test_that("genetic_groups copes with no pair above the threshold, and can draw the graph", {
  rel <- two_families()
  res <- genetic_groups(rel, "r", 0.9)
  expect_equal(length(unique(res$group)), 9)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(genetic_groups(rel, "r", 0.25, plot = TRUE))
})

test_that("genetic_groups refuses non-positive links, whatever the method", {
  rel <- rbind(two_families(), data.frame(ind1 = "Z", ind2 = "B2", r = -0.2))   # QG or LR can be negative
  expect_error(genetic_groups(rel, "r", -1, method = "leiden"), "positive weights")
  expect_error(genetic_groups(rel, "r", -1), "positive weights")
})

test_that("leiden is the default method since 1.0.0", {
  rel <- two_families()
  expect_equal(genetic_groups(rel, "r", 0.25), genetic_groups(rel, "r", 0.25, method = "leiden"))
})

test_that("the null model gives a modularity well above random graphs for clear families", {
  rel <- pkg_data("relate_Carpathians")
  gg <- genetic_groups(rel, "QG", 0.25, n_null = 49)
  mod <- attr(gg, "modularity")
  expect_named(mod, c("Q", "null_mean", "null_sd", "z", "p"))
  expect_gt(mod[["z"]], 3)
  expect_lte(mod[["p"]], 0.05)
  expect_null(attr(genetic_groups(rel, "QG", 0.25), "modularity"))   # no test unless asked
})
