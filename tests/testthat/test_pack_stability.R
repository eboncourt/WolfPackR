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
