# Prepare data
obs <- toy_obs()

test_that("summarise_packs returns a list with statistics by group", {
  result <- summarise_packs(obs, group, sex_column = "Sex", male_pattern = "^Male$", female_pattern = "^Female$")
  expect_type(result, "list")
  expect_setequal(names(result), as.character(unique(obs$group)))
  expect_equal(sum(sapply(result, `[[`, "TotalPoints")), nrow(obs))
})

# Regression test: with a sex column, patterns were compared with ==, so the default
# "M" / "F" never matched "Male" / "Female" and every dominant came out NA
test_that("summarise_packs uses the patterns as regex on the sex column", {
  result <- summarise_packs(obs, group)  # default patterns "M" and "F"
  sex_of <- function(ids) unique(obs$Sex[obs$Individual %in% strsplit(ids, ", ")[[1]]])

  for (p in result) {
    pack_sex <- unique(obs$Sex[obs$group == p$Pack])
    if ("Male" %in% pack_sex) {
      expect_false(is.na(p$PutativeDominantMale))
      expect_equal(sex_of(p$PutativeDominantMale), "Male")
    }
    if ("Female" %in% pack_sex) {
      expect_false(is.na(p$PutativeDominantFemale))
      expect_equal(sex_of(p$PutativeDominantFemale), "Female")
    }
  }
})

# Without a sex column, the sex is read from the IDs
test_that("summarise_packs falls back on the IDs without sex column", {
  obs2 <- obs[, c("Individual", "group")]
  ids <- c(W1 = "M1", W2 = "M2", W3 = "M3", W4 = "M4", W5 = "F1",
           W6 = "F2", W7 = "F3", W8 = "F4", W9 = "F5")
  obs2$Individual <- unname(ids[obs2$Individual])

  result <- summarise_packs(obs2, group, male_pattern = "^M", female_pattern = "^F")
  expect_setequal(names(result), as.character(unique(obs$group)))
  expect_match(result$`1`$PutativeDominantMale, "^M")
})
