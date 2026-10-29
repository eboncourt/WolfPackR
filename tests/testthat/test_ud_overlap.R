test_that("BA and VI are symmetric, 1 on the diagonal, high within a family and low between families", {
  ud <- space_use(sim_samples(), h = 2000)
  for (m in c("BA", "VI")) {
    o <- ud_overlap(ud, m)
    expect_true(isSymmetric(o))
    expect_equal(unname(diag(o)), rep(1, 11), tolerance = 1e-8)
    expect_gt(o["A1", "A2"], 0.5)
    expect_lt(o["A1", "B1"], 0.05)
  }
})

test_that("PHR and HR are read row-wise and PHR's diagonal is the isopleth level", {
  ud <- space_use(sim_samples(), h = 2000)
  phr <- ud_overlap(ud, "PHR", level = 0.9)
  expect_equal(unname(diag(phr)), rep(0.9, 11), tolerance = 0.02)
  hr <- ud_overlap(ud, "HR", individuals = c("A1", "A2", "B1"))
  expect_equal(dim(hr), c(3, 3))
  expect_equal(unname(diag(hr)), rep(1, 3))
})
