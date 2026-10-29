test_that("samples go to the territory they fall in, or outside", {
  ud <- space_use(sim_samples(), h = 2000)
  att <- attribute_samples(ud, list(A = c("A1", "A2", "A3"), B = c("B1", "B2", "B3")))
  expect_true(all(att$Territory[att$Individual %in% c("A1", "A2", "A3", "U")] == "A"))
  expect_true(all(att$Territory[att$Individual %in% c("B1", "B2", "B3", "D")] == "B"))
  expect_gt(mean(att$Territory[att$Individual == "T"] == "outside"), 0.5)
  expect_equal(nrow(att), nrow(ud$samples))
})

test_that("an individual is judged against its territory without itself", {
  ud <- space_use(sim_samples(), h = 2000)
  alone <- attribute_samples(ud, list(A = "A1"))                        # A1 alone: nothing left to compare with
  expect_true(all(alone$Territory[alone$Individual == "A1"] == "outside"))
  kept <- attribute_samples(ud, list(A = "A1"), leave_one_out = FALSE)
  expect_true(all(kept$Territory[kept$Individual == "A1"] == "A"))
})

test_that("territories can be given as a data.frame", {
  ud <- space_use(sim_samples(), h = 2000)
  df <- data.frame(Individual = c("A1", "A2", "B1", "B2", "T"), terr = c("A", "A", "B", "B", NA))
  att <- attribute_samples(ud, df)
  expect_setequal(unique(att$Territory), c("A", "B", "outside"))
})
