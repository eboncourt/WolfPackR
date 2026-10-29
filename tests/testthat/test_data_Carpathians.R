# The Carpathian datasets: wolves only (no dogs any more), dated, and consistent with each other
test_that("samples_Carpathians holds the 263 wolf samples, dated", {
  s <- pkg_data("samples_Carpathians")
  expect_s3_class(s, "sf")
  expect_equal(nrow(s), 263)
  expect_equal(length(unique(s$Individual)), 47)
  expect_equal(names(s)[1], "Individual")                 # some functions want the ID first
  expect_s3_class(s$Date, "Date")
  expect_setequal(unique(s$Campaign), c("2017-18", "2018-19", "2019-20"))
  dogs <- c("CC066T", "CC065T", "CC062C", "CC0607", "CC05TU", "CC05T8", "CC05JA", "CC04JM", "CC04FE", "CC04AP", "CC01CK")
  expect_false(any(s$Individual %in% dogs))
})

test_that("genotypes, relatedness and Bruvo tables cover the same 47 individuals", {
  ids <- sort(unique(pkg_data("samples_Carpathians")$Individual))
  expect_equal(sort(pkg_data("genotypes_Carpathians")$Individual), ids)
  rel <- pkg_data("relate_Carpathians")
  expect_equal(nrow(rel), choose(47, 2))
  expect_setequal(unique(c(rel$ind1, rel$ind2)), ids)
  bru <- pkg_data("Bruvo_Carpathians")
  expect_setequal(unique(c(bru$ind1, bru$ind2)), ids)
})

test_that("effort_Carpathians holds every located sample, and its genotyped ones are samples_Carpathians", {
  e <- pkg_data("effort_Carpathians")
  s <- pkg_data("samples_Carpathians")
  expect_equal(nrow(e), 502)
  expect_setequal(e$Sample[e$Verdict == "Genotyped"], s$Sample)
  expect_true(all(is.na(e$Individual[e$Verdict != "Genotyped"])))
  expect_false("CC01CK" %in% e$Individual)               # the dog labelled "Wolf" in the source file
})

test_that("every sample sits in the study area (no 0, 0 coordinates)", {
  for (nm in c("samples_Carpathians", "effort_Carpathians")) {
    xy <- sf::st_coordinates(pkg_data(nm))
    expect_true(all(xy[, "X"] > 24.5 & xy[, "X"] < 25.5), label = paste(nm, "longitudes"))
    expect_true(all(xy[, "Y"] > 45 & xy[, "Y"] < 46), label = paste(nm, "latitudes"))
  }
})
