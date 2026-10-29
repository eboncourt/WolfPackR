test_that("pack_stability runs the settings and the bootstrap, and scores each individual", {
  st <- pack_stability(sim_samples(), sim_relate(), "r", 0.25, vary = list(h = c(1500, 2500)), n_boot = 2,
                       h = 2000, method = "components")
  expect_equal(nrow(st$runs), 1 + 2 + 2)                               # base + 2 settings + 2 replicates
  expect_equal(st$runs$type, c("base", "sensitivity", "sensitivity", "bootstrap", "bootstrap"))
  expect_true(all(st$individuals$stability >= 0 & st$individuals$stability <= 1))
  expect_equal(dim(st$coassignment), c(11, 11))
  expect_true(isSymmetric(st$coassignment))
  expect_error(pack_stability(sim_samples(), sim_relate(), "r", 0.25, h = 2000), "Nothing")
})

test_that("pack_stability subsamples the data to check the sampling effort", {
  st <- pack_stability(sim_samples(), sim_relate(), "r", 0.25, fraction = c(0.8, 0.5), n_sub = 3,
                       h = 2000, method = "components")
  expect_equal(nrow(st$runs), 1 + 2 * 3)
  expect_equal(unique(st$runs$type), c("base", "subsampling"))
  n <- nrow(sim_samples())
  expect_equal(unique(st$runs$n_samples[st$runs$type == "subsampling"]), c(round(0.8 * n), round(0.5 * n)))
  expect_equal(st$subsampling$fraction, c(0.8, 0.5))
  expect_true(all(st$subsampling$ARI_base_low <= st$subsampling$ARI_base & st$subsampling$ARI_base <= st$subsampling$ARI_base_high))
  expect_true(all(st$individuals$n_runs <= 6))                         # an individual can miss a subsample
  expect_true(all(st$individuals$stability >= 0 & st$individuals$stability <= 1, na.rm = TRUE))
  expect_error(pack_stability(sim_samples(), sim_relate(), "r", 0.25, fraction = 1.5, h = 2000), "between 0 and 1")
})

test_that("the same seed gives the same subsamples", {
  a <- pack_stability(sim_samples(), sim_relate(), "r", 0.25, fraction = 0.6, n_sub = 2, h = 2000, method = "components")
  b <- pack_stability(sim_samples(), sim_relate(), "r", 0.25, fraction = 0.6, n_sub = 2, h = 2000, method = "components")
  expect_identical(a$runs, b$runs)
  expect_null(pack_stability(sim_samples(), sim_relate(), "r", 0.25, n_boot = 1, h = 2000, method = "components")$subsampling)
})
