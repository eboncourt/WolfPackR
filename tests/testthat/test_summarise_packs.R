# Load data
data("relate")
data("samples")

# Prepare data
gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
obs <- merge(samples, gengroups, by = "Individual")

# Execute the function
result <- summarise_packs(obs, group, sex_column = "Sex", male_pattern = "Male", female_pattern = "Female")

test_that("summarise_packs returns a list with statistics by group", {

  expect_is(result, "list")
  expect_all_true(unique(obs$group) %in% names(result))

})


# test with no sex column
obs2 <- obs[,c(1,3,4)]
obs2[obs2$Individual == "W1",]$Individual <- "M1"
obs2[obs2$Individual == "W2",]$Individual <- "M2"
obs2[obs2$Individual == "W3",]$Individual <- "M3"
obs2[obs2$Individual == "W4",]$Individual <- "M4"
obs2[obs2$Individual == "W5",]$Individual <- "F1"
obs2[obs2$Individual == "W6",]$Individual <- "F2"
obs2[obs2$Individual == "W7",]$Individual <- "F3"
obs2[obs2$Individual == "W8",]$Individual <- "F4"
obs2[obs2$Individual == "W9",]$Individual <- "F5"

result2 <- summarise_packs(obs2, group, sex_column = "Sex", male_pattern = "M", female_pattern = "F")
test_that("summarise_packs returns a list with statistics by group", {

  expect_is(result2, "list")
  expect_all_true(unique(obs$group) %in% names(result2))
  expect_all_true(!is.na(result2$`1`$PutativeDominantMale))
})
