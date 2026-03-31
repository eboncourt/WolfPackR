# Load data
data("relate")
data("samples")

# Prepare data
gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
obs <- merge(samples, gengroups, by = "Individual")

# Execute the function
result <- summarise_packs(obs, group)

test_that("summarise_packs returns a list with statistics by group", {

  expect_is(result, "list")
  expect_all_true(unique(obs$group) %in% names(result))

})


