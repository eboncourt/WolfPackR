test_that("wolf_packs runs the whole workflow on the Carpathians and matches Iosif et al. reasonably", {
  res <- suppressMessages(wolf_packs(pkg_data("samples_Carpathians"), pkg_data("relate_Carpathians"), "QG", 0.25,
                                     period = "Campaign"))
  expect_s3_class(res, "wolfpack")
  expect_equal(nrow(res$individuals), 47)
  expect_true(all(c("Pack", "Status", "Genetic_Group", "Group") %in% names(res$samples)))
  cmp <- compare_partitions(res$individuals, pkg_data("reference_Carpathians"))$summary
  expect_gte(cmp$ARI, 0.7)                                             # 0.77 when written
  expect_gte(cmp$reference_groups_found, 5)
  expect_output(print(res), "packs")
})

test_that("wolf_packs works on the toy data", {
  res <- suppressMessages(wolf_packs(pkg_data("samples"), pkg_data("relate"), "wang", 0.4))
  expect_equal(nrow(res$individuals), 9)
  expect_true(any(!is.na(res$individuals$Pack)))
})
