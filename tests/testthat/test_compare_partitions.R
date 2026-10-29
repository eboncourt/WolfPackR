test_that("identical partitions agree fully, whatever the labels", {
  ref <- c(a = "1", b = "1", c = "2", d = "2")
  res <- compare_partitions(c(a = "x", b = "x", c = "y", d = "y"), ref)$summary
  expect_equal(res$ARI, 1); expect_equal(res$NMI, 1); expect_equal(res$n_correct, 4)
})

test_that("unassigned individuals count as alone, and data.frames work", {
  ref <- data.frame(Individual = c("a", "b", "c", "d"), Pack = c("1", "1", "2", "2"))
  x <- data.frame(Individual = c("a", "b", "c", "d", "e"), Pack = c("x", NA, "y", "y", "y"))
  res <- compare_partitions(x, ref)
  expect_equal(res$summary$n_reference, 4)
  expect_equal(res$summary$n_assigned, 3)
  expect_equal(res$summary$n_correct, 3)
  expect_lt(res$summary$ARI, 1)
  expect_equal(res$table["1", "none"], 1)
})
