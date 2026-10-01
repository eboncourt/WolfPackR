# Load data
samples <- pkg_data("samples")
n_ind <- length(unique(samples$Individual))

# mcp_sf() warns when it falls back on a buffer: expected here
quiet_sg <- function(...) suppressWarnings(spatial_groups(...))

test_that("spatial_groups groups individuals spatially", {
  result <- quiet_sg(samples, percentile = 100, buffer_radius = 10,
                     max_iterations = 20, min_mcp_overlap = 0.2)
  expect_s3_class(result, "data.frame")
  expect_true("Subgroup" %in% names(result))
  expect_false(any(is.na(result$Subgroup)))
})

test_that("spatial_groups can deal with given groups", {
  samples2 <- samples
  samples2$group <- as.integer(factor(samples2$Individual))  # one group per individual
  result <- quiet_sg(samples2, group = "group", percentile = 100, buffer_radius = 10,
                     max_iterations = 20, min_mcp_overlap = 0.2)
  expect_s3_class(result, "data.frame")
  expect_false(any(is.na(result$Subgroup)))
  expect_true(all(result$Subgroup == "Lone Individual"))      # alone in its group = lone
})

# Regression test: the output used to have one row per sample, so merging it back
# with the samples gave a cartesian product (277 rows instead of 47 with the toy data)
test_that("spatial_groups returns one row per individual", {
  result <- quiet_sg(samples, percentile = 100, buffer_radius = 10,
                     max_iterations = 20, min_mcp_overlap = 0.2)
  expect_equal(nrow(result), n_ind)
  expect_false(anyDuplicated(result$Individual) > 0)
  expect_equal(nrow(merge(samples, result, by = "Individual")), nrow(samples))
})
