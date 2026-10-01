# Prepare data
samples <- pkg_data("samples")
obs <- toy_obs()
n_ind <- length(unique(samples$Individual))

# mcp_sf() warns every time it falls back on a buffer: expected here, so we keep the output quiet
quiet_ud <- function(...) suppressWarnings(ugly_ducklings(...))

test_that("ugly_ducklings returns one row per individual", {
  result <- quiet_ud(obs, group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Individual", "Genetic_Group", "Pack", "new_group") %in% names(result)))
  expect_equal(nrow(result), n_ind)
  expect_false(anyDuplicated(result$Individual) > 0)
})

test_that("ugly_ducklings can deal with individuals with only 1 point", {
  result <- quiet_ud(obs[-c(2:5), ], group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)
  expect_equal(nrow(result), n_ind)
  expect_false(anyDuplicated(result$Individual) > 0)
})

test_that("ugly_ducklings can deal with individuals with only 2 points", {
  result <- quiet_ud(obs[-c(3:5), ], group = group, min_overlap = 0.7, buffer = 0, mcp.percent = 100)
  expect_equal(nrow(result), n_ind)
  expect_false(anyDuplicated(result$Individual) > 0)
})

# Regression test: an isolated individual with 2 samples used to get two ids (Lone_1 and Lone_2)
test_that("each lone individual gets exactly one Lone_ id", {
  far_away <- fake_ind(obs, "W10", group = 1, n = 2, x = 7.5, y = 45.45)  # 2 samples, far from everyone
  result <- quiet_ud(rbind(obs, far_away), group = group, min_overlap = 0.5)

  w10 <- result[result$Individual == "W10", ]
  expect_equal(nrow(w10), 1)
  expect_equal(w10$Pack, "Lone Individual")
  expect_match(w10$new_group, "^Lone_[0-9]+$")

  lone_ids <- result$new_group[grepl("^Lone_", result$new_group)]
  expect_false(anyDuplicated(lone_ids) > 0)                      # one id = one individual
  expect_equal(nrow(merge(rbind(obs, far_away), result, by = "Individual")),
               nrow(obs) + 2)                                      # merging back doesn't duplicate samples
})

# Regression test: an ugly duckling with no overlap with any pack used to inherit
# the pack of the individual processed just before it (best_pack was never reset)
test_that("an ugly duckling with no overlap at all goes back to its genetic group", {
  in_pack_2 <- fake_ind(obs, "W12", group = 1, n = 5, x = 6.10, y = 45.45)  # group 1, inside pack 2's area
  nowhere   <- fake_ind(obs, "W11", group = 3, n = 3, x = 7.5, y = 45.45)  # group 3, overlaps nobody
  result <- quiet_ud(rbind(obs, in_pack_2, nowhere), group = group, min_overlap = 0.5)

  expect_equal(result$Pack[result$Individual == "W12"], "2")  # reassigned to the pack it overlaps
  expect_equal(result$Pack[result$Individual == "W11"], "3")  # and not to "2", W12's pack
})
